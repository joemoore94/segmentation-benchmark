"""Pairwise inter-method ARI: does Voronoi raise agreement between methods?

The 10x-vs-X comparisons show how each method relates to the platform reference.
This script asks a complementary question: do the Voronoi variants also increase
agreement *among themselves* relative to their nuclear-only counterparts?

All 21 pairs of 7 methods are matched by nearest centroid (max 10µm), clustered
independently at Leiden resolution 1.0, and evaluated with ARI. The result is a
symmetric 7×7 consensus matrix ordered by method family.

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/figures/pairwise_consensus.png
        results/tables/pairwise_ari.csv
        results/tables/pairwise_ari_matrix.csv

Usage::

    conda run -n segbench python scripts/run_pairwise_consensus.py
"""

from __future__ import annotations

import matplotlib.patches as mpatches
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from segbench.compare import cell_type_agreement, cluster_cell_types, match_cells_by_centroid
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES  = Path("results/tables")
FIGURES = Path("results/figures")

MAX_MATCH_DIST = 10.0

# Ordered by family so the heatmap groups naturally.
METHODS = [
    ("10x",               "10x native",   "adata_10x.h5ad"),
    ("cellpose",          "CellPose",     "adata_cellpose.h5ad"),
    ("stardist",          "StarDist",     "adata_stardist.h5ad"),
    ("mesmer",            "Mesmer",       "adata_mesmer.h5ad"),
    ("voronoi",           "Voronoi (CP)", "adata_voronoi.h5ad"),
    ("voronoi_stardist",  "Voronoi (SD)", "adata_voronoi_stardist.h5ad"),
    ("voronoi_mesmer",    "Voronoi (M)",  "adata_voronoi_mesmer.h5ad"),
    ("baysor",            "Baysor",       "adata_baysor.h5ad"),
]

FAMILY = {
    "10x native":   "Reference",
    "CellPose":     "Nuclear",
    "StarDist":     "Nuclear",
    "Mesmer":       "Nuclear",
    "Voronoi (CP)": "Voronoi",
    "Voronoi (SD)": "Voronoi",
    "Voronoi (M)":  "Voronoi",
    "Baysor":       "Transcript-density",
}

FAMILY_COLORS = {
    "Reference":          "#55A868",
    "Nuclear":            "#4C72B0",
    "Voronoi":            "#17BECF",
    "Transcript-density": "#DD8452",
}


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style(scatter=True)

    print("Loading AnnData files...")
    adatas: dict[str, ad.AnnData] = {}
    for _, label, fname in METHODS:
        adatas[label] = ad.read_h5ad(ROI_DIR / fname)
        print(f"  {label}: {adatas[label].n_obs} cells")

    print("\nClustering each method (Leiden resolution 1.0)...")
    labels: dict[str, pd.Series] = {}
    for _, label, _ in METHODS:
        labels[label] = cluster_cell_types(adatas[label], resolution=1.0, seed=0)
        labels[label].index = labels[label].index.astype(str)
        print(f"  {label}: {labels[label].nunique()} clusters")

    method_labels = [m[1] for m in METHODS]
    n = len(method_labels)
    ari_matrix = np.full((n, n), np.nan)
    np.fill_diagonal(ari_matrix, 1.0)

    rows = []
    for i, j in combinations(range(n), 2):
        la, lb = method_labels[i], method_labels[j]
        print(f"\nMatching {la} vs. {lb}...")
        matches = match_cells_by_centroid(adatas[la], adatas[lb], max_dist=MAX_MATCH_DIST)
        matches["id_a"] = matches["id_a"].astype(str)
        matches["id_b"] = matches["id_b"].astype(str)
        result = cell_type_agreement(labels[la], labels[lb], matches)
        ari = result["ari"]
        ari_matrix[i, j] = ari
        ari_matrix[j, i] = ari
        print(f"  ARI={ari:.4f}, n_matched={result['n_matched']}")
        rows.append({"method_a": la, "method_b": lb, "ari": round(ari, 4),
                     "n_matched": result["n_matched"]})

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "pairwise_ari.csv", index=False)
    ari_df = pd.DataFrame(ari_matrix, index=method_labels, columns=method_labels)
    ari_df.to_csv(TABLES / "pairwise_ari_matrix.csv")
    print("\nPairwise ARI matrix:")
    print(ari_df.round(3).to_string())

    # ---------------------------------------------------------------- figure
    apply_style()
    fig, ax = plt.subplots(figsize=(15, 13))

    diag_mask = np.eye(n, dtype=bool)
    off_diag = ari_matrix[~diag_mask]
    vmin = np.nanmin(off_diag)

    sns.heatmap(
        ari_df,
        annot=True, fmt=".3f",
        cmap="YlOrRd",
        vmin=vmin, vmax=1.0,
        mask=diag_mask,
        linewidths=0.5, linecolor="white",
        ax=ax,
        cbar_kws={"label": "ARI", "shrink": 0.8},
        annot_kws={"size": 13, "weight": "bold"},
    )

    # Diagonal cells: grey fill + dash
    for k in range(n):
        ax.add_patch(plt.Rectangle((k, k), 1, 1, fill=True,
                                   color="#dddddd", lw=0, zorder=3))
        ax.text(k + 0.5, k + 0.5, "—", ha="center", va="center",
                fontsize=14, color="#888888", zorder=4)

    # Colour tick labels by family (no patches that cover the labels)
    ax.set_xticklabels(method_labels, rotation=35, ha="right", fontsize=13)
    ax.set_yticklabels(method_labels, rotation=0, fontsize=13)
    for tick in ax.get_xticklabels():
        label_text = tick.get_text()
        if label_text in FAMILY:
            tick.set_color(FAMILY_COLORS[FAMILY[label_text]])
            tick.set_fontweight("bold")
    for tick in ax.get_yticklabels():
        label_text = tick.get_text()
        if label_text in FAMILY:
            tick.set_color(FAMILY_COLORS[FAMILY[label_text]])
            tick.set_fontweight("bold")

    ax.set_title(
        "Pairwise ARI between all segmentation methods\n"
        "(Leiden resolution 1.0, nearest-centroid matching ≤ 10 µm)",
        fontweight="bold", fontsize=14,
    )

    handles = [mpatches.Patch(color=c, label=f) for f, c in FAMILY_COLORS.items()]
    ax.legend(handles=handles, loc="upper right",
              bbox_to_anchor=(1.35, 1.02), fontsize=12, title="Family",
              title_fontsize=12)

    fig.tight_layout()
    fig.savefig(FIGURES / "pairwise_consensus.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved pairwise_consensus.png")


if __name__ == "__main__":
    main()
