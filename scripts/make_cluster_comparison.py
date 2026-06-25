"""Compare Leiden clustering structure across segmentation methods.

Two-panel figure: number of clusters per method (left) and cells-per-cluster
distribution as a box plot (right).  Both panels use resolution 1.0.

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/figures/cluster_comparison.png

Usage::

    conda run -n segbench python scripts/make_cluster_comparison.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from segbench.compare import cluster_cell_types
from segbench.constants import (
    MAIN_METHODS,
    METHOD_COLORS,
    METHOD_FAMILIES,
    METHOD_LABELS,
)
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
FIGURES = Path("results/figures")

METHODS = [m for m in MAIN_METHODS if m != "10x_native"]

FILE_OVERRIDES = {
    "10x_native": "adata_10x.h5ad",
}

FAMILY_COLORS = {
    "Reference": "#333333",
    "Voronoi": "#1B9E77",
    "Transcript-density": "#D95F02",
}


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style()

    rows: list[dict] = []

    # Load 10x native first (different filename convention)
    h5_10x = ROI_DIR / FILE_OVERRIDES["10x_native"]
    adata_10x = ad.read_h5ad(h5_10x)
    labels_10x = cluster_cell_types(adata_10x, resolution=1.0)
    counts_10x = labels_10x.value_counts()
    for cluster_id, n in counts_10x.items():
        rows.append({
            "method": METHOD_LABELS["10x_native"],
            "method_key": "10x_native",
            "family": METHOD_FAMILIES["10x_native"],
            "cluster": cluster_id,
            "cells": n,
        })
    print(f"  {METHOD_LABELS['10x_native']}: {counts_10x.nunique()} clusters, {len(labels_10x)} cells")

    for method in METHODS:
        h5 = ROI_DIR / f"adata_{method}.h5ad"
        if not h5.exists():
            print(f"  {METHOD_LABELS.get(method, method)}: skipped (not found)")
            continue
        adata = ad.read_h5ad(h5)
        labels = cluster_cell_types(adata, resolution=1.0)
        counts = labels.value_counts()
        label = METHOD_LABELS[method]
        family = METHOD_FAMILIES[method]
        for cluster_id, n in counts.items():
            rows.append({
                "method": label,
                "method_key": method,
                "family": family,
                "cluster": cluster_id,
                "cells": n,
            })
        print(f"  {label}: {counts.nunique()} clusters, {len(labels)} cells")

    df = pd.DataFrame(rows)

    all_keys = ["10x_native"] + METHODS
    method_order = [
        METHOD_LABELS[m] for m in all_keys
        if METHOD_LABELS[m] in df["method"].values
    ]

    cluster_counts = (
        df.groupby("method")["cluster"]
        .nunique()
        .reindex(method_order)
    )
    families = (
        df.drop_duplicates("method")
        .set_index("method")["family"]
        .reindex(method_order)
    )
    colors = [FAMILY_COLORS[f] for f in families]

    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=(28, 10),
        gridspec_kw={"width_ratios": [1, 2.2]},
    )

    # --- Left panel: cluster count ---
    y_pos = np.arange(len(method_order))
    ax_left.barh(y_pos, cluster_counts.values, color=colors, edgecolor="white",
                 linewidth=0.5, height=0.7)
    ax_left.set_yticks(y_pos)
    ax_left.set_yticklabels(method_order)
    ax_left.set_xlabel("Leiden clusters")
    ax_left.set_title("Number of clusters", fontweight="bold")
    ax_left.invert_yaxis()
    for i, v in enumerate(cluster_counts.values):
        ax_left.text(v + 0.3, i, str(int(v)), va="center", fontsize=16)
    ax_left.set_xlim(0, cluster_counts.max() + 4)

    # --- Right panel: cells per cluster (box plot) ---
    box_data = [
        df[df["method"] == m]["cells"].values for m in method_order
    ]
    bp = ax_right.boxplot(
        box_data, vert=False, patch_artist=True,
        widths=0.6,
        medianprops=dict(color="black", linewidth=2),
        flierprops=dict(marker="o", markersize=5, alpha=0.6),
    )
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.8)

    ax_right.set_yticklabels(method_order)
    ax_right.set_xlabel("Cells per cluster")
    ax_right.set_title("Cluster size distribution", fontweight="bold")
    ax_right.invert_yaxis()

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=FAMILY_COLORS[f], label=f, alpha=0.7)
        for f in ["Reference", "Voronoi", "Transcript-density"]
    ]
    ax_right.legend(handles=legend_handles, loc="lower right", fontsize=18)

    fig.suptitle(
        "Leiden clustering comparison (resolution 1.0)",
        fontsize=28, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "cluster_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {FIGURES / 'cluster_comparison.png'}")


if __name__ == "__main__":
    main()
