"""Cross-lineage bleed spatial analysis.

For each segmentation method, identifies cells of one lineage (e.g. T cells)
that aberrantly express markers of a different lineage (e.g. CAF markers),
then maps those cells spatially alongside the "source" cell type to show that
high-bleed cells cluster at lineage boundaries — a signature of segmentation
mis-assignment rather than genuine co-expression.

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/figures/cross_lineage_bleed.png
        results/figures/cross_lineage_bleed_detail.png

Usage::

    conda run -n segbench python scripts/run_cross_lineage_bleed.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import cKDTree

from segbench.constants import (
    CELLTYPE_COLORS,
    CLUSTER_ANNOTATIONS,
    METHOD_LABELS,
)
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
FIGURES = Path("results/figures")
TABLES = Path("results/tables")

CAF_MARKERS = ["LUM", "SFRP4", "FBLN1", "CCDC80", "THBS2", "MMP2"]
T_CELL_MARKERS = ["CD3E", "CD3G", "TRAC", "TRBC1", "CD96", "IL7R", "CCL5"]

METHODS = [
    ("10x_native",              "adata_10x.h5ad"),
    ("voronoi",                 "adata_voronoi.h5ad"),
    ("voronoi_mesmer",          "adata_voronoi_mesmer.h5ad"),
    ("baysor_prior_c10",        "adata_baysor_prior_c10.h5ad"),
    ("watershed_10x",           "adata_watershed_10x.h5ad"),
    ("cellpose_cyto3_density",  "adata_cellpose_cyto3_density.h5ad"),
    ("mesmer_wholecell_density", "adata_mesmer_wholecell_density.h5ad"),
]


def _cluster_and_annotate(adata: ad.AnnData) -> ad.AnnData:
    sc.settings.verbosity = 0
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph")
    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(CLUSTER_ANNOTATIONS).astype("category")
    )
    return adata


def _bleed_score(adata: ad.AnnData, markers: list[str]) -> np.ndarray:
    """Mean log-normalized expression of foreign markers per cell."""
    available = [g for g in markers if g in adata.var_names]
    if not available:
        return np.zeros(adata.n_obs)
    X = adata[:, available].X
    if sp.issparse(X):
        X = X.toarray()
    return np.mean(X, axis=1)


def _nearest_source_dist(
    target_coords: np.ndarray,
    source_coords: np.ndarray,
) -> np.ndarray:
    """Distance from each target cell to its nearest source cell (microns)."""
    if len(source_coords) == 0:
        return np.full(len(target_coords), np.inf)
    tree = cKDTree(source_coords)
    dists, _ = tree.query(target_coords, k=1)
    return dists


def fig_bleed_spatial(
    adata: ad.AnnData,
    method_label: str,
    ax: plt.Axes,
    target_type: str = "T cells",
    source_type: str = "CAFs",
    bleed_markers: list[str] = CAF_MARKERS,
) -> dict:
    """Plot one spatial panel: target cells colored by bleed, source cells shown."""
    ct = adata.obs["cell_type"]
    cx = adata.obs["centroid_x"].values
    cy = adata.obs["centroid_y"].values

    target_mask = ct == target_type
    source_mask = ct == source_type
    other_mask = ~target_mask & ~source_mask

    bleed = _bleed_score(adata, bleed_markers)
    target_bleed = bleed[target_mask]

    bleed_thresh = float(np.median(target_bleed))
    high_bleed_mask = target_mask & (bleed > bleed_thresh)

    ax.scatter(
        cx[other_mask & ~source_mask], cy[other_mask & ~source_mask],
        c="#E8E8E8", s=1, alpha=0.1, rasterized=True,
    )
    ax.scatter(
        cx[source_mask], cy[source_mask],
        c=CELLTYPE_COLORS[source_type], s=8, alpha=0.45,
        label=source_type, rasterized=True,
    )

    if high_bleed_mask.any():
        sc = ax.scatter(
            cx[high_bleed_mask], cy[high_bleed_mask],
            c=bleed[high_bleed_mask],
            cmap="YlOrRd", s=18, alpha=0.9, zorder=5,
            vmin=bleed_thresh,
            vmax=np.percentile(bleed[high_bleed_mask], 98),
            edgecolors="black", linewidth=0.4,
            label=f"{target_type} (above-median {source_type} bleed)",
            rasterized=True,
        )
        cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label(f"Mean {source_type} marker expr", fontsize=12)

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (μm)")
    ax.set_ylabel("y (μm)")
    ax.set_title(method_label, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, markerscale=2)

    # Compute stats
    target_coords = np.column_stack([cx[target_mask], cy[target_mask]])
    source_coords = np.column_stack([cx[source_mask], cy[source_mask]])
    high_coords = np.column_stack([cx[high_bleed_mask], cy[high_bleed_mask]]) if high_bleed_mask.any() else np.empty((0, 2))

    dist_all_target = _nearest_source_dist(target_coords, source_coords)
    dist_high = _nearest_source_dist(high_coords, source_coords) if len(high_coords) > 0 else np.array([])

    return {
        "method": method_label,
        "n_target": int(target_mask.sum()),
        "n_source": int(source_mask.sum()),
        "n_high_bleed": int(high_bleed_mask.sum()),
        "pct_high_bleed": high_bleed_mask.sum() / max(target_mask.sum(), 1) * 100,
        "median_dist_all_target_to_source": float(np.median(dist_all_target)),
        "median_dist_high_bleed_to_source": float(np.median(dist_high)) if len(dist_high) > 0 else float("nan"),
    }


def fig_bleed_detail(
    adata: ad.AnnData,
    method_label: str,
    target_type: str = "T cells",
    source_type: str = "CAFs",
    bleed_markers: list[str] = CAF_MARKERS,
) -> plt.Figure:
    """Detailed single-method figure: spatial map + scatter + histogram."""
    ct = adata.obs["cell_type"]
    cx = adata.obs["centroid_x"].values
    cy = adata.obs["centroid_y"].values

    target_mask = ct == target_type
    source_mask = ct == source_type

    bleed = _bleed_score(adata, bleed_markers)
    target_bleed = bleed[target_mask]

    bleed_thresh = float(np.median(target_bleed))
    high = target_bleed > bleed_thresh
    low = ~high

    target_coords = np.column_stack([cx[target_mask], cy[target_mask]])
    source_coords = np.column_stack([cx[source_mask], cy[source_mask]])
    dist_to_source = _nearest_source_dist(target_coords, source_coords)

    apply_style(scatter=True)
    fig, axes = plt.subplots(1, 3, figsize=(30, 9))

    # Left: spatial map — above-median T cells + CAFs
    other_mask = ~target_mask & ~source_mask
    axes[0].scatter(cx[other_mask], cy[other_mask], c="#E8E8E8", s=1, alpha=0.1, rasterized=True)
    axes[0].scatter(
        cx[source_mask], cy[source_mask],
        c=CELLTYPE_COLORS[source_type], s=8, alpha=0.45,
        label=source_type, rasterized=True,
    )
    high_mask_full = target_mask & (bleed > bleed_thresh)
    sc = axes[0].scatter(
        cx[high_mask_full], cy[high_mask_full],
        c=bleed[high_mask_full], cmap="YlOrRd", s=18, alpha=0.9, zorder=5,
        vmin=bleed_thresh,
        vmax=np.percentile(bleed[high_mask_full], 98) if high_mask_full.any() else 1,
        edgecolors="black", linewidth=0.4, rasterized=True,
    )
    cbar = plt.colorbar(sc, ax=axes[0], shrink=0.6, pad=0.02)
    cbar.set_label(f"Mean {source_type} marker expr", fontsize=14)
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    axes[0].set_xlabel("x (μm)")
    axes[0].set_ylabel("y (μm)")
    axes[0].set_title(f"Above-median {target_type} + {source_type}", fontweight="bold")
    axes[0].legend(loc="upper right", fontsize=11, markerscale=2)

    # Center: bleed score vs. distance to nearest source cell
    axes[1].scatter(
        dist_to_source[low], target_bleed[low],
        c=CELLTYPE_COLORS[target_type], s=6, alpha=0.2,
        label="Below median", rasterized=True,
    )
    axes[1].scatter(
        dist_to_source[high], target_bleed[high],
        c="firebrick", s=12, alpha=0.5,
        label="Above median", rasterized=True,
    )
    axes[1].axhline(bleed_thresh, color="gray", ls="--", lw=1, label=f"Median ({bleed_thresh:.2f})")
    axes[1].set_xlabel(f"Distance to nearest {source_type} (μm)")
    axes[1].set_ylabel(f"Mean {source_type} marker expression")
    axes[1].set_title("Bleed vs. proximity", fontweight="bold")
    axes[1].legend(fontsize=11)

    # Right: distance histograms split by above/below median bleed
    cap = min(200, dist_to_source.max())
    bins = np.linspace(0, cap, 40)
    axes[2].hist(
        dist_to_source[low], bins=bins, alpha=0.5, density=True,
        color=CELLTYPE_COLORS[target_type], label="Below median bleed",
    )
    axes[2].hist(
        dist_to_source[high], bins=bins, alpha=0.5, density=True,
        color="firebrick", label="Above median bleed",
    )
    axes[2].set_xlabel(f"Distance to nearest {source_type} (μm)")
    axes[2].set_ylabel("Density")
    axes[2].set_title(f"{target_type} proximity to {source_type}", fontweight="bold")
    axes[2].legend(fontsize=11)

    fig.suptitle(
        f"Cross-lineage bleed: {target_type} × {source_type} markers — {method_label}",
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    # Multi-method comparison panel
    apply_style(scatter=True)
    _n = len(METHODS)
    _ncols = min(_n, 4)
    _nrows = (_n + _ncols - 1) // _ncols
    fig, axes = plt.subplots(_nrows, _ncols, figsize=(6 * _ncols, 5.5 * _nrows))
    axes = axes.ravel()

    stats_rows = []

    for i, (key, fname) in enumerate(METHODS):
        path = ROI_DIR / fname
        if not path.exists():
            print(f"  {key}: skipped (not found)")
            continue

        label = METHOD_LABELS[key]
        print(f"Processing {label}...")
        adata = ad.read_h5ad(path)
        adata = _cluster_and_annotate(adata)

        stats = fig_bleed_spatial(adata, label, axes[i])
        stats_rows.append(stats)

        print(f"  {stats['n_target']} T cells, {stats['n_high_bleed']} with high CAF bleed "
              f"({stats['pct_high_bleed']:.1f}%)")
        print(f"  Median dist to CAF: all T cells={stats['median_dist_all_target_to_source']:.1f} μm, "
              f"high-bleed={stats['median_dist_high_bleed_to_source']:.1f} μm")

    fig.suptitle(
        "T cells expressing CAF markers: spatial distribution across methods",
        fontweight="bold", fontsize=22,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES / "cross_lineage_bleed.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {FIGURES / 'cross_lineage_bleed.png'}")

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(TABLES / "cross_lineage_bleed_stats.csv", index=False)
    print(f"Saved {TABLES / 'cross_lineage_bleed_stats.csv'}")

    print("\n=== Summary ===")
    print(stats_df.to_string(index=False))

    # Detailed single-method figure (Voronoi CP — typically highest bleed)
    detail_key, detail_fname = "voronoi", "adata_voronoi.h5ad"
    detail_path = ROI_DIR / detail_fname
    if detail_path.exists():
        print(f"\nGenerating detail figure for {METHOD_LABELS[detail_key]}...")
        adata = ad.read_h5ad(detail_path)
        adata = _cluster_and_annotate(adata)
        fig_detail = fig_bleed_detail(adata, METHOD_LABELS[detail_key])
        fig_detail.savefig(
            FIGURES / "cross_lineage_bleed_detail.png",
            dpi=200, bbox_inches="tight",
        )
        plt.close(fig_detail)
        print(f"Saved {FIGURES / 'cross_lineage_bleed_detail.png'}")


if __name__ == "__main__":
    main()
