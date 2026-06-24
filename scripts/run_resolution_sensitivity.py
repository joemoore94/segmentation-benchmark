"""Leiden resolution sensitivity: does the ARI ordering survive resolution changes?

Runs the normalize -> PCA -> neighbors -> Leiden pipeline at multiple resolutions
for each method and the 10x-native reference, computes ARI vs. 10x native, and
checks whether the method ordering is stable.  Disagreement and Moran's I are
computed under both Hungarian (one-to-one) and argmax (many-to-one) cluster
matching to show how the choice of alignment algorithm affects the results.

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/tables/resolution_sensitivity.csv
        results/figures/resolution_sensitivity_hungarian.png
        results/figures/resolution_sensitivity_argmax.png

Usage::

    conda run -n segbench python scripts/run_resolution_sensitivity.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from segbench.compare import (
    cell_type_agreement,
    cluster_cell_types,
    match_cluster_labels,
    match_cluster_labels_argmax,
)
from segbench.constants import METHOD_COLORS as _MC, METHOD_LABELS, NUCLEAR_ONLY
from segbench.spatial import disagreement_table, morans_i
from segbench.style import apply_style

ROI_DIR   = Path("data/processed/roi")
TABLES    = Path("results/tables")
FIGURES   = Path("results/figures")

RESOLUTIONS = [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]
_TICK_LABELS = {r: str(r) for r in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]}

_ALL_METHODS = [
    "cellpose", "stardist", "mesmer",
    "voronoi", "voronoi_stardist", "voronoi_mesmer",
    "baysor", "baysor_prior_c08", "baysor_prior_c10", "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10", "bidcell", "segger",
]
METHODS = [(k, METHOD_LABELS[k]) for k in _ALL_METHODS]

METHOD_COLORS = {METHOD_LABELS[k]: _MC[k] for k in _ALL_METHODS}

PLOT_METHODS = [m for m in METHODS if m[0] not in NUCLEAR_ONLY]

_MATCHERS = {
    "hungarian": match_cluster_labels,
    "argmax": match_cluster_labels_argmax,
}


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style()

    print("Loading AnnData files...")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    available = [(m, l) for m, l in METHODS if (ROI_DIR / f"adata_{m}.h5ad").exists()]
    skipped = [(m, l) for m, l in METHODS if not (ROI_DIR / f"adata_{m}.h5ad").exists()]
    for m, l in skipped:
        print(f"  {l}: skipped (file not found)")
    adatas = {m: ad.read_h5ad(ROI_DIR / f"adata_{m}.h5ad") for m, _ in available}

    print("Loading match tables...")
    matches = {
        m: pd.read_csv(TABLES / f"disagreement_table_10x_{m}.csv")
        for m, _ in available
    }

    rows = []
    for res in RESOLUTIONS:
        print(f"\nResolution {res}:")
        labels_10x = cluster_cell_types(adata_10x, resolution=res)
        labels_10x.index = labels_10x.index.astype(str)
        n_10x = labels_10x.nunique()
        for method, label in available:
            labels_comp = cluster_cell_types(adatas[method], resolution=res)
            labels_comp.index = labels_comp.index.astype(str)
            n_comp = labels_comp.nunique()
            m = matches[method].copy()
            m["id_a"] = m["id_a"].astype(str)
            m["id_b"] = m["id_b"].astype(str)
            result = cell_type_agreement(labels_10x, labels_comp, m)
            ari = result["ari"]

            row_base = {
                "resolution": res,
                "method": label,
                "ari": round(ari, 4),
                "n_clusters_10x": n_10x,
                "n_clusters_comp": n_comp,
            }

            for matcher_name, matcher_fn in _MATCHERS.items():
                labels_aligned = matcher_fn(labels_10x, labels_comp, m)
                dt = disagreement_table(m, labels_10x, labels_aligned, adata_10x)
                coords = dt[["centroid_x", "centroid_y"]].to_numpy()
                mi = morans_i(coords, dt["disagree"].to_numpy(dtype=float))
                disagree_pct = dt["disagree"].mean() * 100
                row = {
                    **row_base,
                    "matcher": matcher_name,
                    "disagree_pct": round(disagree_pct, 2),
                    "morans_i": round(mi, 4),
                }
                rows.append(row)

            print(f"  {label:18s}: ARI={ari:.4f}  clusters: 10x={n_10x}, {label}={n_comp}")

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "resolution_sensitivity.csv", index=False)
    print("\nSaved resolution_sensitivity.csv")

    # ---------------------------------------------------------------- figures
    for matcher_name in _MATCHERS:
        _plot_resolution_sensitivity(df, matcher_name)

    print("Saved resolution_sensitivity_hungarian.png")
    print("Saved resolution_sensitivity_argmax.png")


def _plot_resolution_sensitivity(df: pd.DataFrame, matcher: str) -> None:
    matcher_label = "Hungarian (one-to-one)" if matcher == "hungarian" else "Argmax (many-to-one)"
    sub_df = df[df["matcher"] == matcher]

    fig, axes = plt.subplots(3, 1, figsize=(14, 24))

    # Panel 1: ARI vs resolution
    ax = axes[0]
    for _, label in PLOT_METHODS:
        s = sub_df[sub_df["method"] == label]
        if s.empty:
            continue
        ax.plot(s["resolution"], s["ari"], "o-", color=METHOD_COLORS[label],
                label=label, linewidth=2.5, markersize=8)
    ax.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.4, label="default (1.0)")
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("ARI vs. 10x native")
    ax.set_title("ARI across Leiden resolutions", fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_xticks(RESOLUTIONS)
    ax.set_xticklabels([_TICK_LABELS.get(r, "") for r in RESOLUTIONS])

    # Panel 2: Disagreement % vs resolution
    ax2 = axes[1]
    for _, label in PLOT_METHODS:
        s = sub_df[sub_df["method"] == label]
        if s.empty:
            continue
        ax2.plot(s["resolution"], s["disagree_pct"], "o-", color=METHOD_COLORS[label],
                 label=label, linewidth=2.5, markersize=8)
    ax2.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.4)
    ax2.set_xlabel("Leiden resolution")
    ax2.set_ylabel("Disagreement (%)")
    ax2.set_title("Cell-type disagreement across resolutions", fontweight="bold")
    ax2.set_xticks(RESOLUTIONS)
    ax2.set_xticklabels([_TICK_LABELS.get(r, "") for r in RESOLUTIONS])
    ax2.legend(fontsize=10)

    # Panel 3: Moran's I vs resolution
    ax3 = axes[2]
    for _, label in PLOT_METHODS:
        s = sub_df[sub_df["method"] == label]
        if s.empty:
            continue
        ax3.plot(s["resolution"], s["morans_i"], "o-", color=METHOD_COLORS[label],
                 label=label, linewidth=2.5, markersize=8)
    ax3.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.4)
    ax3.set_xlabel("Leiden resolution")
    ax3.set_ylabel("Global Moran's I of disagreement")
    ax3.set_title("Spatial structure of disagreement across resolutions", fontweight="bold")
    ax3.set_xticks(RESOLUTIONS)
    ax3.set_xticklabels([_TICK_LABELS.get(r, "") for r in RESOLUTIONS])
    ax3.legend(fontsize=10)

    fig.suptitle(
        f"Leiden resolution sensitivity — {matcher_label} cluster alignment",
        fontsize=18, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / f"resolution_sensitivity_{matcher}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
