"""Negative marker analysis: biology-grounded segmentation quality metric.

For each segmentation method, scores cells by the rate of biologically
impossible co-expression patterns. A cell co-expressing CD3E (T cell) and
GATA3 (luminal epithelial) is almost certainly a segmentation artifact where
two adjacent cells' transcripts were merged into one.

Unlike ARI-based metrics, this requires no single method as ground truth -
the biological exclusivity of lineage markers is the reference.

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/tables/negative_marker_violations.csv
        results/tables/negative_marker_summary.csv
        results/figures/negative_marker.png

Usage::

    conda run -n segbench python scripts/run_negative_marker.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns

from segbench.constants import (
    METHOD_COLORS,
    METHOD_LABELS,
    NEGATIVE_PAIRS_TIER1,
    NEGATIVE_PAIRS_TIER2,
    NUCLEAR_ONLY,
)
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES = Path("results/tables")
FIGURES = Path("results/figures")

MIN_COUNTS = 2

_ALL_METHODS = [
    ("10x_native",        "adata_10x.h5ad"),
    ("cellpose",          "adata_cellpose.h5ad"),
    ("stardist",          "adata_stardist.h5ad"),
    ("mesmer",            "adata_mesmer.h5ad"),
    ("voronoi",           "adata_voronoi.h5ad"),
    ("voronoi_stardist",  "adata_voronoi_stardist.h5ad"),
    ("voronoi_mesmer",    "adata_voronoi_mesmer.h5ad"),
    ("baysor",            "adata_baysor.h5ad"),
    ("baysor_prior_c08",  "adata_baysor_prior_c08.h5ad"),
    ("baysor_prior_c10",  "adata_baysor_prior_c10.h5ad"),
    ("bidcell",           "adata_bidcell.h5ad"),
    ("segger",            "adata_segger.h5ad"),
]


def _get_counts(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs)
    idx = list(adata.var_names).index(gene)
    col = adata.X[:, idx]
    if sp.issparse(col):
        return np.asarray(col.todense()).ravel()
    return np.asarray(col).ravel()


def score_method(
    adata: ad.AnnData,
    pairs: list[tuple[str, str, str, str]],
    min_counts: int = MIN_COUNTS,
) -> pd.DataFrame:
    rows = []
    for gene_a, gene_b, type_a, type_b in pairs:
        counts_a = _get_counts(adata, gene_a)
        counts_b = _get_counts(adata, gene_b)
        violations = (counts_a >= min_counts) & (counts_b >= min_counts)
        rows.append({
            "gene_a": gene_a,
            "gene_b": gene_b,
            "type_a": type_a,
            "type_b": type_b,
            "pair": f"{gene_a}+{gene_b}",
            "n_violations": int(violations.sum()),
            "violation_rate": violations.mean(),
        })
    return pd.DataFrame(rows)


def contamination_score(
    adata: ad.AnnData,
    pairs: list[tuple[str, str, str, str]],
) -> np.ndarray:
    total = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    total = np.maximum(total, 1.0)
    score = np.zeros(adata.n_obs)
    for gene_a, gene_b, _, _ in pairs:
        ca = _get_counts(adata, gene_a)
        cb = _get_counts(adata, gene_b)
        score += np.minimum(ca, cb)
    return score / total


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    all_pairs = NEGATIVE_PAIRS_TIER1 + NEGATIVE_PAIRS_TIER2

    summary_rows = []
    detail_rows = []

    for key, fname in _ALL_METHODS:
        path = ROI_DIR / fname
        if not path.exists():
            print(f"  {METHOD_LABELS[key]}: skipped (file not found)")
            continue

        label = METHOD_LABELS[key]
        adata = ad.read_h5ad(path)
        print(f"\n=== {label} ({adata.n_obs} cells) ===")

        pair_df = score_method(adata, all_pairs)
        pair_df["method"] = label
        detail_rows.append(pair_df)

        total_tx = float(np.asarray(adata.X.sum(axis=1)).ravel().sum())
        median_tx = float(np.median(np.asarray(adata.X.sum(axis=1)).ravel()))

        tier1_df = score_method(adata, NEGATIVE_PAIRS_TIER1)
        counts_a_all = np.zeros(adata.n_obs, dtype=bool)
        for gene_a, gene_b, _, _ in NEGATIVE_PAIRS_TIER1:
            ca = _get_counts(adata, gene_a)
            cb = _get_counts(adata, gene_b)
            counts_a_all |= (ca >= MIN_COUNTS) & (cb >= MIN_COUNTS)
        tier1_violation_rate = counts_a_all.mean()

        cont = contamination_score(adata, all_pairs)

        summary_rows.append({
            "method": label,
            "method_key": key,
            "n_cells": adata.n_obs,
            "median_tx_per_cell": median_tx,
            "tier1_violation_rate": tier1_violation_rate,
            "tier1_violations": int(counts_a_all.sum()),
            "mean_contamination_score": cont.mean(),
            "median_contamination_score": float(np.median(cont)),
            "violations_per_1000tx": tier1_violation_rate * 1000 / max(median_tx, 1),
        })
        print(f"  Tier 1 violation rate: {tier1_violation_rate:.4f} "
              f"({counts_a_all.sum()} / {adata.n_obs} cells)")
        print(f"  Median tx/cell: {median_tx:.0f}")
        print(f"  Violations per 1000 tx: {summary_rows[-1]['violations_per_1000tx']:.2f}")

        for _, row in tier1_df.iterrows():
            if row["n_violations"] > 0:
                print(f"    {row['pair']}: {row['n_violations']} violations "
                      f"({row['violation_rate']:.4f})")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(TABLES / "negative_marker_summary.csv", index=False)

    detail_df = pd.concat(detail_rows, ignore_index=True)
    detail_df.to_csv(TABLES / "negative_marker_violations.csv", index=False)

    print("\n=== Summary ===")
    print(summary_df[["method", "n_cells", "median_tx_per_cell",
                       "tier1_violation_rate", "violations_per_1000tx"]]
          .to_string(index=False))

    # ---------------------------------------------------------------- figure
    plot_summary = summary_df[~summary_df["method_key"].isin(NUCLEAR_ONLY)].copy()
    nuclear_labels = {METHOD_LABELS[k] for k in NUCLEAR_ONLY}
    plot_detail = detail_df[~detail_df["method"].isin(nuclear_labels)].copy()

    apply_style()
    fig = plt.figure(figsize=(30, 10))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1.5, 1], wspace=0.4)

    ax_a = fig.add_subplot(gs[0])
    order = plot_summary.sort_values("tier1_violation_rate", ascending=True)
    colors = [METHOD_COLORS.get(k, "#999999") for k in order["method_key"]]
    ax_a.barh(range(len(order)), order["tier1_violation_rate"] * 100, color=colors)
    ax_a.set_yticks(range(len(order)))
    ax_a.set_yticklabels(order["method"])
    ax_a.set_xlabel("Cells with Tier 1 violation (%)")
    ax_a.set_title("Violation rate", fontweight="bold")

    ax_b = fig.add_subplot(gs[1])
    tier1_pairs = [f"{a}+{b}" for a, b, _, _ in NEGATIVE_PAIRS_TIER1]
    heatmap_data = plot_detail[plot_detail["pair"].isin(tier1_pairs)].pivot_table(
        index="pair", columns="method", values="violation_rate", fill_value=0
    )
    method_order = [m for m in plot_summary.sort_values(
        "tier1_violation_rate", ascending=False
    )["method"] if m in heatmap_data.columns]
    heatmap_data = heatmap_data.reindex(columns=method_order)
    sns.heatmap(
        heatmap_data * 100, annot=True, fmt=".2f", cmap="YlOrRd",
        ax=ax_b, cbar_kws={"label": "Violation rate (%)"},
        linewidths=0.5, linecolor="white",
    )
    ax_b.set_title("Per-pair violation rate (%)", fontweight="bold")
    ax_b.set_xlabel("")
    ax_b.set_ylabel("")

    ax_c = fig.add_subplot(gs[2])
    for _, row in plot_summary.iterrows():
        ax_c.scatter(
            row["median_tx_per_cell"],
            row["tier1_violation_rate"] * 100,
            color=METHOD_COLORS.get(row["method_key"], "#999999"),
            s=160, zorder=5, edgecolors="white", linewidth=1.5,
        )
        ax_c.annotate(
            row["method"], (row["median_tx_per_cell"], row["tier1_violation_rate"] * 100),
            textcoords="offset points", xytext=(10, 0),
        )
    ax_c.set_xlabel("Median transcripts per cell")
    ax_c.set_ylabel("Cells with Tier 1 violation (%)")
    ax_c.set_title("Violation vs. transcript capture", fontweight="bold")

    fig.suptitle(
        "Negative marker analysis: biology-grounded segmentation quality",
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIGURES / "negative_marker.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {FIGURES / 'negative_marker.png'}")


if __name__ == "__main__":
    main()
