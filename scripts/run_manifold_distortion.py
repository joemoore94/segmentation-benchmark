"""Phenotypic manifold distortion induced by segmentation method.

Fits a PCA model on the 10x-native reference, projects all other methods into
that shared latent space, then runs a joint UMAP so every method's cells share
the same coordinate system. Compares phenotypic landscape coverage by computing
2D kernel density estimates per method and showing enrichment / depletion
relative to 10x native.

Key questions answered:
- Does segmentation alter the shape of the phenotypic landscape?
- Which cell states are gained or lost by each method?
- Are manifold distortions concentrated in specific biological populations?

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/figures/manifold_distortion.png
        results/figures/manifold_shared_umap.png

Usage::

    conda run -n segbench python scripts/run_manifold_distortion.py
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
import umap as umap_lib
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
FIGURES = Path("results/figures")

from segbench.constants import METHOD_COLORS as _MC, METHOD_LABELS, NUCLEAR_ONLY

_ALL_METHODS = [
    "cellpose", "stardist", "mesmer",
    "voronoi", "voronoi_stardist", "voronoi_mesmer",
    "baysor", "baysor_prior_c08", "bidcell", "segger",
]
METHODS = [(k, METHOD_LABELS[k]) for k in _ALL_METHODS]

METHOD_COLORS = {METHOD_LABELS[k]: _MC[k] for k in _ALL_METHODS}
METHOD_COLORS["10x native"] = _MC["10x_native"]

PLOT_METHODS = [m for m in METHODS if m[0] not in NUCLEAR_ONLY]

N_PCS     = 30
N_SAMPLE  = 5_000   # cells per method for UMAP (full set used for density)
GRID_SIZE = 100     # KDE grid resolution
RANDOM_STATE = 0


def normalize_log(X) -> np.ndarray:
    """Normalize-total then log1p, returns dense float32 array."""
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    X = X / totals * np.median(totals)
    return np.log1p(X)


def density_grid(xy: np.ndarray, grid_x: np.ndarray, grid_y: np.ndarray) -> np.ndarray:
    """Normalised 2D KDE on a precomputed meshgrid."""
    kde = gaussian_kde(xy.T, bw_method=0.15)
    positions = np.vstack([grid_x.ravel(), grid_y.ravel()])
    z = kde(positions).reshape(grid_x.shape)
    return z / z.sum()


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style(scatter=True)
    rng = np.random.default_rng(RANDOM_STATE)

    print("Loading AnnData files...")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    available_methods = [(m, l) for m, l in METHODS if (ROI_DIR / f"adata_{m}.h5ad").exists()]
    for m, l in METHODS:
        if not (ROI_DIR / f"adata_{m}.h5ad").exists():
            print(f"  {l}: skipped (file not found)")
    adatas = {label: ad.read_h5ad(ROI_DIR / f"adata_{m}.h5ad") for m, label in available_methods}

    # ---------------------------------------------------------------- shared gene set
    all_genes = set(adata_10x.var_names)
    for label, adata in adatas.items():
        all_genes &= set(adata.var_names)
    shared_genes = sorted(all_genes)
    print(f"Shared genes across all methods: {len(shared_genes)}")

    # ---------------------------------------------------------------- normalize + PCA
    print("Normalising and fitting shared PCA on 10x native...")
    X_10x = normalize_log(adata_10x[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    Z_10x = pca.fit_transform(X_10x)
    print(f"  Variance explained by {N_PCS} PCs: {pca.explained_variance_ratio_.sum():.1%}")

    print("Projecting comparison methods into shared PCA space...")
    projections: dict[str, np.ndarray] = {"10x native": Z_10x}
    for m, label in available_methods:
        X = normalize_log(adatas[label][:, shared_genes].X)
        projections[label] = pca.transform(X)
        print(f"  {label}: {X.shape[0]} cells projected")

    # ---------------------------------------------------------------- shared UMAP
    # Subsample for fitting UMAP; project full sets afterwards.
    print("\nFitting shared UMAP on 10x native subsample...")
    idx_10x = rng.choice(len(Z_10x), size=min(N_SAMPLE, len(Z_10x)), replace=False)
    reducer = umap_lib.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                             random_state=RANDOM_STATE)
    reducer.fit(Z_10x[idx_10x])

    print("Projecting all methods into shared UMAP...")
    umap_coords: dict[str, np.ndarray] = {}
    for method_label, Z in projections.items():
        umap_coords[method_label] = reducer.transform(Z)
        print(f"  {method_label}: done")

    # ---------------------------------------------------------------- figure 1: shared UMAP coloured by method
    print("\nRendering shared UMAP figure...")
    fig, ax = plt.subplots(figsize=(14, 11))
    for method_label in ["10x native"] + [lab for _, lab in available_methods]:
        xy = umap_coords[method_label]
        idx = rng.choice(len(xy), size=min(2000, len(xy)), replace=False)
        ax.scatter(xy[idx, 0], xy[idx, 1],
                   c=METHOD_COLORS[method_label], s=3, alpha=0.35,
                   label=method_label, rasterized=True)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("All methods in shared 10x-native PCA/UMAP space", fontweight="bold")
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=METHOD_COLORS[lab], markersize=10,
                          label=lab)
               for lab in ["10x native"] + [lab for _, lab in available_methods]]
    ax.legend(handles=handles, fontsize=10, loc="best")
    fig.suptitle("Phenotypic landscape: shared reference space", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES / "manifold_shared_umap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved manifold_shared_umap.png")

    # ---------------------------------------------------------------- figure 2: density distortion
    # Build a common grid from 10x UMAP extent.
    all_u = np.concatenate([umap_coords["10x native"]] +
                           [umap_coords[lab] for _, lab in available_methods])
    x_min, x_max = np.percentile(all_u[:, 0], [1, 99])
    y_min, y_max = np.percentile(all_u[:, 1], [1, 99])
    gx, gy = np.meshgrid(np.linspace(x_min, x_max, GRID_SIZE),
                         np.linspace(y_min, y_max, GRID_SIZE))

    print("\nComputing 2D KDE for each method...")
    dens_10x = density_grid(umap_coords["10x native"], gx, gy)
    dens: dict[str, np.ndarray] = {}
    for _, label in available_methods:
        dens[label] = density_grid(umap_coords[label], gx, gy)
        print(f"  {label}: done")

    # log2 density ratio: positive = enriched vs 10x, negative = depleted
    eps = dens_10x.max() * 1e-4
    log2_ratio: dict[str, np.ndarray] = {
        label: np.log2((d + eps) / (dens_10x + eps))
        for label, d in dens.items()
    }
    vmax = max(np.abs(v).max() for v in log2_ratio.values())
    vmax = min(vmax, 3.0)   # cap at ±3 log2 for colour saturation

    fig, axes = plt.subplots(2, 2, figsize=(20, 18))
    cmap = "RdBu_r"
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    for ax, (_, label) in zip(axes.flatten(), PLOT_METHODS):
        im = ax.pcolormesh(gx, gy, log2_ratio[label], cmap=cmap, norm=norm,
                           shading="auto", rasterized=True)
        # Overlay 10x native UMAP as faint grey dots for orientation
        idx = rng.choice(len(umap_coords["10x native"]),
                         size=min(1500, len(umap_coords["10x native"])), replace=False)
        ax.scatter(umap_coords["10x native"][idx, 0],
                   umap_coords["10x native"][idx, 1],
                   c="white", s=1, alpha=0.15, rasterized=True)
        plt.colorbar(im, ax=ax, label="log₂(method / 10x density)")
        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    fig.suptitle(
        "Phenotypic landscape distortion vs. 10x native\n"
        "Red = cell states enriched in comparison method  ·  Blue = depleted",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "manifold_distortion.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved manifold_distortion.png")


if __name__ == "__main__":
    main()
