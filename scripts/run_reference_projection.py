"""Reference transcriptome projection.

Fits PCA on companion scRNA-seq (GSE243275 3' sample) subsetted to the
380-gene Xenium panel, projects each segmentation method's cells into that
shared reference space, then clusters in that space and compares results.

Analysis levels:
  1. Cell-by-cell projection + density comparison vs reference
  2. Leiden clustering in reference PCA space (shared embedding across methods)
  3. ARI between reference-space clusterings across methods
  4. Spatial validation: reference-space cluster labels mapped back to tissue

Reads:  data/reference/scrna_3p_filtered_feature_bc_matrix.h5
        data/processed/roi/adata_*.h5ad
Writes: results/figures/ref_projection_umap.png
        results/figures/ref_projection_density.png
        results/figures/ref_projection_ref_clustering.png
        results/figures/ref_projection_spatial.png
        results/figures/ref_projection_ari.png

Usage::

    conda run -n segbench python scripts/run_reference_projection.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score

from segbench.constants import (
    CLUSTER_ANNOTATIONS,
    CELLTYPE_COLORS,
    METHOD_COLORS as _MC,
    NUCLEAR_ONLY,
)
from segbench.compare import cluster_cell_types
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
REF_PATH = Path("data/reference/scrna_3p_filtered_feature_bc_matrix.h5")
FIGURES = Path("results/figures")

# key → (file stem, display label)
# 10x native is adata_10x.h5ad (not adata_10x_native.h5ad)
EXPANSION_METHODS = [
    ("10x",                         "10x native"),
    ("voronoi",                     "Voronoi (CP)"),
    ("voronoi_stardist",            "Voronoi (SD)"),
    ("voronoi_mesmer",              "Voronoi (M)"),
    ("voronoi_10x_ranger",          "Voronoi (10x)"),
    ("baysor",                      "Baysor"),
    ("baysor_prior_c10",            "Baysor (CP prior 1.0)"),
    ("baysor_stardist_prior_c10",   "Baysor (SD prior 1.0)"),
    ("baysor_mesmer_prior_c10",     "Baysor (M prior 1.0)"),
    ("baysor_10x_ranger_prior_c10", "Baysor (10x prior 1.0)"),
]

METHOD_COLORS = {}
for key, label in EXPANSION_METHODS:
    mc_key = key if key != "10x" else "10x_native"
    if mc_key in _MC:
        METHOD_COLORS[label] = _MC[mc_key]

N_PCS = 30
LEIDEN_RES = 1.0
N_SAMPLE_UMAP = 5_000
GRID_SIZE = 100
RANDOM_STATE = 0


def normalize_log(X) -> np.ndarray:
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    X = X / totals * np.median(totals)
    return np.log1p(X)


def density_grid(xy: np.ndarray, grid_x: np.ndarray, grid_y: np.ndarray) -> np.ndarray:
    kde = gaussian_kde(xy.T, bw_method=0.15)
    positions = np.vstack([grid_x.ravel(), grid_y.ravel()])
    z = kde(positions).reshape(grid_x.shape)
    return z / z.sum()


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style(scatter=True)
    rng = np.random.default_rng(RANDOM_STATE)

    # ---------------------------------------------------------------- load reference
    print("Loading scRNA-seq reference...")
    ref = sc.read_10x_h5(str(REF_PATH))
    ref.var_names_make_unique()
    print(f"  Raw: {ref.n_obs} cells, {ref.n_vars} genes")

    sc.pp.filter_cells(ref, min_genes=200)
    sc.pp.filter_genes(ref, min_cells=3)
    ref = ref[ref.obs["n_genes"] < 6000].copy()
    print(f"  After QC: {ref.n_obs} cells, {ref.n_vars} genes")

    # ---------------------------------------------------------------- load Xenium
    print("Loading Xenium AnnData files...")
    adata_10x_raw = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    xenium_genes = sorted(set(adata_10x_raw.var_names))

    available: list[tuple[str, str]] = []
    adatas_raw: dict[str, ad.AnnData] = {}
    for key, label in EXPANSION_METHODS:
        path = ROI_DIR / f"adata_{key}.h5ad"
        if path.exists():
            adatas_raw[label] = ad.read_h5ad(path)
            available.append((key, label))
            print(f"  {label}: {adatas_raw[label].n_obs} cells")
        else:
            print(f"  {label}: skipped (file not found)")

    # ---------------------------------------------------------------- shared gene set
    ref_genes = set(ref.var_names)
    shared_genes = sorted(set(xenium_genes) & ref_genes)
    print(f"\nXenium panel genes in reference: {len(shared_genes)} / {len(xenium_genes)}")

    # ---------------------------------------------------------------- fit reference PCA
    print("Fitting PCA on scRNA-seq reference ({}-gene subset)...".format(len(shared_genes)))
    X_ref = normalize_log(ref[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    Z_ref = pca.fit_transform(X_ref)
    var_explained = pca.explained_variance_ratio_.sum()
    print(f"  {N_PCS} PCs explain {var_explained:.1%} of variance in reference")

    # ---------------------------------------------------------------- 10x native cell type labels (own-space)
    print("\nClustering 10x native in its own space for cell type labels...")
    adata_10x_own = adata_10x_raw.copy()
    sc.pp.normalize_total(adata_10x_own)
    sc.pp.log1p(adata_10x_own)
    sc.pp.pca(adata_10x_own, n_comps=30, random_state=0)
    sc.pp.neighbors(adata_10x_own, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata_10x_own, resolution=1.0, random_state=0, flavor="igraph")
    adata_10x_own.obs["cell_type"] = adata_10x_own.obs["leiden"].map(CLUSTER_ANNOTATIONS)
    ct_labels_10x = adata_10x_own.obs["cell_type"].values

    # ---------------------------------------------------------------- project all methods into ref PCA
    print("\nProjecting all methods into reference PCA space...")
    projections: dict[str, np.ndarray] = {}
    for key, label in available:
        adata = adatas_raw[label]
        X_norm = normalize_log(adata[:, shared_genes].X)
        projections[label] = pca.transform(X_norm)
        print(f"  {label}: {adata.n_obs} cells projected")

    # ---------------------------------------------------------------- cluster each method in ref PCA space
    print(f"\nClustering each method in reference PCA space (Leiden res={LEIDEN_RES})...")
    ref_leiden: dict[str, np.ndarray] = {}
    for label in projections:
        Z = projections[label]
        tmp = ad.AnnData(
            X=Z,
            obs=pd.DataFrame(index=np.arange(Z.shape[0]).astype(str)),
        )
        sc.pp.neighbors(tmp, use_rep="X", n_neighbors=15, random_state=0)
        sc.tl.leiden(tmp, resolution=LEIDEN_RES, random_state=0, flavor="igraph")
        ref_leiden[label] = tmp.obs["leiden"].values
        n_cl = len(set(ref_leiden[label]))
        print(f"  {label}: {n_cl} clusters")

    # ---------------------------------------------------------------- fit UMAP on reference
    print("\nFitting UMAP on reference subsample...")
    import umap as umap_lib
    idx_ref_sample = rng.choice(len(Z_ref), size=min(N_SAMPLE_UMAP, len(Z_ref)), replace=False)
    reducer = umap_lib.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                            random_state=RANDOM_STATE)
    reducer.fit(Z_ref[idx_ref_sample])

    print("Projecting all into shared UMAP...")
    umap_coords: dict[str, np.ndarray] = {"scRNA-seq reference": reducer.transform(Z_ref)}
    for label, Z in projections.items():
        umap_coords[label] = reducer.transform(Z)
        print(f"  {label}: done")

    ref_xy = umap_coords["scRNA-seq reference"]
    plot_methods = [(k, l) for k, l in available if k not in NUCLEAR_ONLY]
    n_methods = len(plot_methods)
    ncols = min(3, n_methods)
    nrows = (n_methods + ncols - 1) // ncols

    # ================================================================ FIGURE 1: cell-by-cell UMAP
    print("\nRendering cell-by-cell UMAP figure...")
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, (key, label) in enumerate(plot_methods):
        ax = axes[i]
        idx = rng.choice(len(ref_xy), size=min(3000, len(ref_xy)), replace=False)
        ax.scatter(ref_xy[idx, 0], ref_xy[idx, 1],
                   c="#DDDDDD", s=2, alpha=0.3, rasterized=True)
        method_xy = umap_coords[label]
        idx_m = rng.choice(len(method_xy), size=min(3000, len(method_xy)), replace=False)
        ax.scatter(method_xy[idx_m, 0], method_xy[idx_m, 1],
                   c=METHOD_COLORS.get(label, "#333333"), s=3, alpha=0.4,
                   rasterized=True)
        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Xenium cells projected into scRNA-seq reference PCA/UMAP space\n"
        f"(PCA fit on {ref.n_obs:,} scRNA-seq cells, {len(shared_genes)} genes, "
        f"{var_explained:.0%} variance explained)",
        fontsize=18, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_umap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_umap.png")

    # ================================================================ FIGURE 2: density distortion vs reference
    print("\nComputing density ratios vs scRNA-seq reference...")
    all_u = np.concatenate([ref_xy] + [umap_coords[l] for _, l in plot_methods])
    x_min, x_max = np.percentile(all_u[:, 0], [1, 99])
    y_min, y_max = np.percentile(all_u[:, 1], [1, 99])
    gx, gy = np.meshgrid(np.linspace(x_min, x_max, GRID_SIZE),
                         np.linspace(y_min, y_max, GRID_SIZE))
    dens_ref = density_grid(ref_xy, gx, gy)
    eps = dens_ref.max() * 1e-4

    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
    axes = np.atleast_1d(axes).flatten()
    vmax = 3.0

    for i, (key, label) in enumerate(plot_methods):
        ax = axes[i]
        d = density_grid(umap_coords[label], gx, gy)
        log2r = np.log2((d + eps) / (dens_ref + eps))
        norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
        im = ax.pcolormesh(gx, gy, log2r, cmap="RdBu_r", norm=norm,
                           shading="auto", rasterized=True)
        idx = rng.choice(len(ref_xy), size=min(1500, len(ref_xy)), replace=False)
        ax.scatter(ref_xy[idx, 0], ref_xy[idx, 1],
                   c="white", s=1, alpha=0.15, rasterized=True)
        plt.colorbar(im, ax=ax, label="log₂(Xenium / scRNA-seq density)")
        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Density distortion vs. scRNA-seq reference\n"
        "Red = Xenium enriched  |  Blue = Xenium depleted",
        fontsize=18, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_density.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_density.png")

    # ================================================================ FIGURE 3: reference-space clustering UMAP
    print("\nRendering reference-space clustering UMAPs...")
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, (key, label) in enumerate(plot_methods):
        ax = axes[i]
        method_xy = umap_coords[label]
        clusters = ref_leiden[label]
        n_cl = len(set(clusters))
        cmap = plt.colormaps.get_cmap("tab20").resampled(n_cl)
        cluster_ids = sorted(set(clusters), key=lambda c: int(c))
        color_map = {c: cmap(j) for j, c in enumerate(cluster_ids)}
        colors = [color_map[c] for c in clusters]

        idx_m = rng.choice(len(method_xy), size=min(5000, len(method_xy)), replace=False)
        ax.scatter(method_xy[idx_m, 0], method_xy[idx_m, 1],
                   c=[colors[k] for k in idx_m], s=3, alpha=0.5, rasterized=True)
        ax.set_title(f"{label} ({n_cl} clusters)", fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f"Leiden clustering in reference PCA space (resolution {LEIDEN_RES})\n"
        "Each method's cells clustered in the same scRNA-seq-derived coordinate system",
        fontsize=18, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_ref_clustering.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_ref_clustering.png")

    # ================================================================ FIGURE 4: pairwise ARI in reference space
    print("\nComputing pairwise ARI in reference PCA space...")
    from scipy.spatial import cKDTree

    method_labels = [l for _, l in plot_methods]
    n = len(method_labels)
    ari_matrix = np.zeros((n, n))

    for i_m, label_a in enumerate(method_labels):
        adata_a = adatas_raw[label_a]
        xy_a = np.column_stack([adata_a.obs["centroid_x"], adata_a.obs["centroid_y"]])
        tree_a = cKDTree(xy_a)
        clusters_a = ref_leiden[label_a]

        for j_m, label_b in enumerate(method_labels):
            if i_m == j_m:
                ari_matrix[i_m, j_m] = 1.0
                continue
            if j_m < i_m:
                ari_matrix[i_m, j_m] = ari_matrix[j_m, i_m]
                continue

            adata_b = adatas_raw[label_b]
            xy_b = np.column_stack([adata_b.obs["centroid_x"], adata_b.obs["centroid_y"]])
            clusters_b = ref_leiden[label_b]

            dists, idx_match = tree_a.query(xy_b)
            matched = dists < 15
            if matched.sum() < 100:
                ari_matrix[i_m, j_m] = np.nan
                continue

            ari_val = adjusted_rand_score(
                clusters_a[idx_match[matched]],
                clusters_b[matched],
            )
            ari_matrix[i_m, j_m] = ari_val
            ari_matrix[j_m, i_m] = ari_val
            print(f"  {label_a} vs {label_b}: ARI = {ari_val:.3f}")

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(ari_matrix, cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(n))
    ax.set_xticklabels(method_labels, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(method_labels)

    for i_m in range(n):
        for j_m in range(n):
            val = ari_matrix[i_m, j_m]
            if not np.isnan(val):
                ax.text(j_m, i_m, f"{val:.2f}", ha="center", va="center",
                        fontsize=12, fontweight="bold",
                        color="white" if val > 0.6 else "black")

    plt.colorbar(im, ax=ax, label="Adjusted Rand Index", shrink=0.8)
    ax.set_title(
        "Pairwise ARI: Leiden clustering in scRNA-seq reference PCA space\n"
        "(matched cells by nearest centroid, <15 µm)",
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_ari.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_ari.png")

    # ================================================================ FIGURE 5: spatial maps of ref-space clusters
    print("\nRendering spatial maps of reference-space clusters...")

    spatial_methods = [
        (k, l) for k, l in plot_methods
        if l in ("10x native", "Voronoi (M)", "Baysor", "Baysor (M prior 1.0)")
    ]
    if not spatial_methods:
        spatial_methods = plot_methods[:4]

    n_sp = len(spatial_methods)
    fig, axes = plt.subplots(1, n_sp, figsize=(7 * n_sp, 6))
    axes = np.atleast_1d(axes).flatten()

    for i, (key, label) in enumerate(spatial_methods):
        ax = axes[i]
        adata = adatas_raw[label]
        clusters = ref_leiden[label]
        n_cl = len(set(clusters))
        cmap = plt.colormaps.get_cmap("tab20").resampled(n_cl)
        cluster_ids = sorted(set(clusters), key=lambda c: int(c))
        color_map = {c: cmap(j) for j, c in enumerate(cluster_ids)}
        colors = [color_map[c] for c in clusters]

        ax.scatter(adata.obs["centroid_x"], adata.obs["centroid_y"],
                   c=colors, s=1, alpha=0.6, rasterized=True)
        ax.set_title(f"{label}\n({n_cl} ref-space clusters)", fontweight="bold")
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_aspect("equal")
        ax.invert_yaxis()

    fig.suptitle(
        "Reference-space Leiden clusters mapped back to tissue coordinates",
        fontsize=18, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_spatial.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_spatial.png")

    # ================================================================ FIGURE 6: matched-cell displacement
    # For each comparison method, match cells to 10x native by nearest centroid,
    # then measure how far the same cell moves in reference PCA space.
    print("\nComputing matched-cell displacement in reference PCA space...")
    from scipy.spatial import cKDTree

    Z_10x = projections["10x native"]
    xy_10x = np.column_stack([
        adata_10x_own.obs["centroid_x"], adata_10x_own.obs["centroid_y"]
    ])
    celltypes = sorted(set(ct_labels_10x))

    compare_methods = [(k, l) for k, l in plot_methods if l != "10x native"]

    disp_by_method: dict[str, np.ndarray] = {}
    disp_by_ct: dict[str, dict[str, np.ndarray]] = {}

    for key, label in compare_methods:
        adata_m = adatas_raw[label]
        Z_m = projections[label]
        xy_m = np.column_stack([adata_m.obs["centroid_x"], adata_m.obs["centroid_y"]])
        tree_m = cKDTree(xy_m)
        match_dists, match_idx = tree_m.query(xy_10x)
        matched = match_dists < 15

        displacement = np.sqrt(
            ((Z_10x[matched] - Z_m[match_idx[matched]]) ** 2).sum(axis=1)
        )
        disp_by_method[label] = displacement
        print(f"  {label}: {matched.sum()} matched, "
              f"median displacement = {np.median(displacement):.2f}")

        ct_disp = {}
        for ct in celltypes:
            ct_mask = (ct_labels_10x == ct) & matched
            if ct_mask.sum() < 20:
                continue
            ct_disp[ct] = np.sqrt(
                ((Z_10x[ct_mask] - Z_m[match_idx[ct_mask]]) ** 2).sum(axis=1)
            )
        disp_by_ct[label] = ct_disp

    # --- violin: displacement by method ---
    fig, ax = plt.subplots(figsize=(16, 8))
    labels_ordered = [l for l in [l for _, l in compare_methods]]
    data = [disp_by_method[l] for l in labels_ordered]
    parts = ax.violinplot(data, positions=range(len(labels_ordered)),
                          showmedians=True, showextrema=False)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(METHOD_COLORS.get(labels_ordered[i], "#999999"))
        pc.set_alpha(0.7)
    parts["cmedians"].set_color("black")
    ax.set_xticks(range(len(labels_ordered)))
    ax.set_xticklabels(labels_ordered, rotation=45, ha="right")
    ax.set_ylabel("Displacement from 10x native position\n(Euclidean distance in 30-PC reference space)")
    ax.set_title(
        "How far does the same cell move in reference space depending on segmentation method?\n"
        "(matched cells by nearest centroid, <15 µm)",
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_displacement.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_displacement.png")

    # --- grouped bar: displacement by cell type ---
    n_ct = len(celltypes)
    n_meth = len(disp_by_ct)
    fig, ax = plt.subplots(figsize=(18, 9))
    x_pos = np.arange(n_ct)
    width = 0.8 / max(n_meth, 1)

    for i, (label, ct_disp) in enumerate(disp_by_ct.items()):
        medians = [float(np.median(ct_disp[ct])) if ct in ct_disp else 0
                   for ct in celltypes]
        ax.bar(x_pos + i * width, medians, width,
               color=METHOD_COLORS.get(label, "#333333"),
               label=label, alpha=0.85)

    ax.set_xticks(x_pos + width * n_meth / 2)
    ax.set_xticklabels(celltypes, rotation=45, ha="right")
    ax.set_ylabel("Median displacement from 10x native\n(Euclidean in 30-PC reference space)")
    ax.set_title(
        "Per-cell-type displacement in scRNA-seq reference space\n"
        "(how far segmentation method choice moves matched cells from their 10x native position)",
        fontweight="bold",
    )
    ax.legend(fontsize=10, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURES / "ref_projection_displacement_by_ct.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved ref_projection_displacement_by_ct.png")

    # ================================================================ summary table
    # Internal-space clustering for before/after comparison
    print("\nClustering each method in its own PCA space...")
    internal_leiden: dict[str, pd.Series] = {}
    for key, label in available:
        internal_leiden[label] = cluster_cell_types(
            adatas_raw[label], resolution=LEIDEN_RES, seed=RANDOM_STATE,
        )
        internal_leiden[label].index = internal_leiden[label].index.astype(str)

    # ARI vs 10x native's internal clustering (fixed anchor).
    # "Before": method X in own PCA vs 10x native in own PCA.
    # "After":  method X in ref PCA vs 10x native in own PCA.
    # Delta isolates the effect of projecting method X into the reference space.
    xy_10x_sum = np.column_stack([
        adatas_raw["10x native"].obs["centroid_x"],
        adatas_raw["10x native"].obs["centroid_y"],
    ])
    tree_10x_sum = cKDTree(xy_10x_sum)
    tenx_internal = internal_leiden["10x native"]

    ari_before: dict[str, float] = {}
    ari_after: dict[str, float] = {}
    for key, label in available:
        if label == "10x native":
            ari_before[label] = 1.0
            ari_after[label] = adjusted_rand_score(
                tenx_internal.values,
                ref_leiden[label],
            )
            continue
        xy_m = np.column_stack([
            adatas_raw[label].obs["centroid_x"],
            adatas_raw[label].obs["centroid_y"],
        ])
        dists, idx_match = tree_10x_sum.query(xy_m)
        matched = dists < 15
        tenx_matched = tenx_internal.values[idx_match[matched]]
        ari_before[label] = adjusted_rand_score(
            tenx_matched,
            internal_leiden[label].values[matched],
        )
        ari_after[label] = adjusted_rand_score(
            tenx_matched,
            ref_leiden[label][matched],
        )

    print("\n" + "=" * 95)
    print("SUMMARY: ARI vs 10x native (own-PCA clustering) — before and after projection")
    print("=" * 95)
    print(f"{'Method':<30} {'Ref clusters':>12}  {'ARI (before)':>12}  {'ARI (after)':>12}  {'ΔARI':>7}  {'Med. disp.':>11}")
    print("-" * 90)
    for _, label in available:
        n_cl = len(set(ref_leiden[label]))
        med_d = (f"{np.median(disp_by_method[label]):.2f}"
                 if label in disp_by_method else "—")
        ab = ari_before[label]
        aa = ari_after[label]
        d_ari = aa - ab
        print(f"{label:<30} {n_cl:>12}  {ab:>12.3f}  {aa:>12.3f}  {d_ari:>+7.3f}  {med_d:>11}")

    print("\nDone.")


if __name__ == "__main__":
    main()
