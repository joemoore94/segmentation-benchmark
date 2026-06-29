"""Cluster centroid distance: all segmentation methods vs scRNA-seq reference.

For each segmentation method and Leiden resolution, projects cells into a
shared PCA space (fit on scRNA-seq reference, 374 shared genes, 30 PCs),
clusters in that space, computes cluster centroids, and measures the Euclidean
distance to the nearest scRNA-seq reference centroid. This answers: which
method and resolution produce clusters closest to the reference cell states?

Also computes per-cell-type centroid distances using marker-based annotation.

Reads:  data/reference/scrna_3p_filtered_feature_bc_matrix.h5
        data/processed/roi/adata_*.h5ad
Writes: results/tables/centroid_distance.csv
        results/tables/centroid_distance_summary.csv
        results/tables/celltype_centroid_distance.csv

Usage::

    conda run -n segbench python scripts/run_tangent_deviation.py
"""

from __future__ import annotations

from itertools import product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import PCA

from segbench.constants import METHOD_LABELS, NUCLEAR_ONLY

ROI_DIR  = Path("data/processed/roi")
REF_PATH = Path("data/reference/scrna_3p_filtered_feature_bc_matrix.h5")
TABLES   = Path("results/tables")

RESOLUTIONS = [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]
N_PCS = 30
RANDOM_STATE = 0

METHODS = [
    ("10x_native",                  "adata_10x.h5ad"),
    ("voronoi",                     "adata_voronoi.h5ad"),
    ("voronoi_stardist",            "adata_voronoi_stardist.h5ad"),
    ("voronoi_mesmer",              "adata_voronoi_mesmer.h5ad"),
    ("voronoi_10x_ranger",          "adata_voronoi_10x_ranger.h5ad"),
    ("baysor",                      "adata_baysor.h5ad"),
    ("baysor_prior_c10",            "adata_baysor_prior_c10.h5ad"),
    ("baysor_stardist_prior_c10",   "adata_baysor_stardist_prior_c10.h5ad"),
    ("baysor_mesmer_prior_c10",     "adata_baysor_mesmer_prior_c10.h5ad"),
    ("baysor_10x_ranger_prior_c10", "adata_baysor_10x_ranger_prior_c10.h5ad"),
]

CANONICAL_MARKERS: dict[str, list[str]] = {
    "Luminal epithelial": [
        "EPCAM", "KRT7", "KRT8", "GATA3", "ESR1", "PGR", "FOXA1",
        "MUC1", "ANKRD30A", "TACSTD2", "CCND1",
    ],
    "Myoepithelial": ["KRT14", "KRT5", "ACTA2", "MYLK", "DST"],
    "T cells": ["CD3E", "CD3G", "TRAC", "TRBC1", "CD96", "IL7R", "CCL5"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19"],
    "Plasma cells": ["MZB1", "SLAMF7", "TENT5C", "TNFRSF17"],
    "Macrophages": ["CD14", "CD68", "CD163", "AIF1", "LYZ", "FCER1G"],
    "CAFs": ["LUM", "SFRP4", "FBLN1", "CCDC80", "THBS2", "MMP2", "PDGFRA"],
    "Smooth muscle": ["MYH11", "ACTA2", "MYLK", "RGS5", "CAV1"],
    "Endothelial": ["PECAM1", "VWF", "AQP1", "CD93", "CLEC14A"],
    "Adipocytes": ["ADIPOQ", "PLIN1", "PPARG", "LPL", "G0S2"],
}


def normalize_log(X) -> np.ndarray:
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    X = X / totals * np.median(totals)
    return np.log1p(X)


def cluster_centroids(Z: np.ndarray, labels: np.ndarray) -> dict[str, np.ndarray]:
    centroids = {}
    for cl in sorted(set(labels), key=lambda c: int(c)):
        centroids[cl] = Z[labels == cl].mean(axis=0)
    return centroids


def score_and_assign(adata_raw, gene_universe):
    a = adata_raw.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    scores = {}
    for ct, markers in CANONICAL_MARKERS.items():
        valid = [g for g in markers if g in gene_universe]
        if len(valid) < 2:
            continue
        sc.tl.score_genes(a, valid, score_name=f"score_{ct}")
        scores[ct] = a.obs[f"score_{ct}"].values
    return pd.DataFrame(scores, index=a.obs_names).idxmax(axis=1)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0

    # ---- load and QC reference
    print("Loading scRNA-seq reference...")
    ref = sc.read_10x_h5(str(REF_PATH))
    ref.var_names_make_unique()
    sc.pp.filter_cells(ref, min_genes=200)
    sc.pp.filter_genes(ref, min_cells=3)
    ref = ref[ref.obs["n_genes"] < 6000].copy()
    print(f"  {ref.n_obs} cells after QC")

    # ---- load all segmentation methods
    print("\nLoading segmentation methods...")
    adatas_raw: dict[str, ad.AnnData] = {}
    available: list[tuple[str, str]] = []
    for key, fname in METHODS:
        path = ROI_DIR / fname
        if not path.exists():
            print(f"  {METHOD_LABELS[key]}: skipped")
            continue
        label = METHOD_LABELS[key]
        adatas_raw[label] = ad.read_h5ad(path)
        available.append((key, label))
        print(f"  {label}: {adatas_raw[label].n_obs} cells")

    # ---- shared gene set + PCA (fit on reference)
    first_label = available[0][1]
    shared_genes = sorted(set(adatas_raw[first_label].var_names) & set(ref.var_names))
    print(f"\nShared genes: {len(shared_genes)}")

    X_ref = normalize_log(ref[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    Z_ref = pca.fit_transform(X_ref)
    print(f"PCA: {N_PCS} PCs, {pca.explained_variance_ratio_.sum():.1%} variance")

    # ---- project all methods into reference PCA
    projections: dict[str, np.ndarray] = {}
    for key, label in available:
        X_m = normalize_log(adatas_raw[label][:, shared_genes].X)
        projections[label] = pca.transform(X_m)

    # ---- cluster reference at each resolution
    print("\nClustering reference...")
    ref_ad = ad.AnnData(X=Z_ref, obs=pd.DataFrame(index=ref.obs_names))
    sc.pp.neighbors(ref_ad, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)
    ref_clusterings: dict[float, np.ndarray] = {}
    for res in RESOLUTIONS:
        sc.tl.leiden(ref_ad, resolution=res, random_state=RANDOM_STATE, flavor="igraph")
        ref_clusterings[res] = ref_ad.obs["leiden"].values.copy()
    print(f"  {len(RESOLUTIONS)} resolutions, "
          f"{len(set(ref_clusterings[0.3]))}–{len(set(ref_clusterings[2.0]))} clusters")

    # ---- cluster each method at each resolution
    print("Clustering methods...")
    method_clusterings: dict[str, dict[float, np.ndarray]] = {}
    for key, label in available:
        Z_m = projections[label]
        m_ad = ad.AnnData(X=Z_m, obs=pd.DataFrame(index=adatas_raw[label].obs_names))
        sc.pp.neighbors(m_ad, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)
        method_clusterings[label] = {}
        for res in RESOLUTIONS:
            sc.tl.leiden(m_ad, resolution=res, random_state=RANDOM_STATE, flavor="igraph")
            method_clusterings[label][res] = m_ad.obs["leiden"].values.copy()
        n_lo = len(set(method_clusterings[label][0.3]))
        n_hi = len(set(method_clusterings[label][2.0]))
        print(f"  {label}: {n_lo}–{n_hi} clusters")

    # ---- compute centroid distances for all (method, ref_res, method_res) combos
    print("\nComputing centroid distances...")
    detail_rows = []
    summary_rows = []

    for key, label in available:
        Z_m = projections[label]
        for ref_res, m_res in product(RESOLUTIONS, RESOLUTIONS):
            ref_centroids = cluster_centroids(Z_ref, ref_clusterings[ref_res])
            m_centroids = cluster_centroids(Z_m, method_clusterings[label][m_res])
            m_labels = method_clusterings[label][m_res]

            ref_centers = np.array(list(ref_centroids.values()))
            ref_keys = list(ref_centroids.keys())

            distances = []
            for cl, centroid in m_centroids.items():
                dists = np.linalg.norm(ref_centers - centroid, axis=1)
                nearest_idx = int(np.argmin(dists))
                nearest_dist = float(dists[nearest_idx])
                n_cells = int((m_labels == cl).sum())

                detail_rows.append({
                    "method": label,
                    "ref_res": ref_res,
                    "method_res": m_res,
                    "cluster": cl,
                    "n_cells": n_cells,
                    "nearest_ref_cluster": ref_keys[nearest_idx],
                    "distance": round(nearest_dist, 4),
                })
                distances.append((nearest_dist, n_cells))

            dists_arr = np.array([d for d, _ in distances])
            sizes = np.array([n for _, n in distances])
            weights = sizes / sizes.sum()

            summary_rows.append({
                "method": label,
                "ref_res": ref_res,
                "method_res": m_res,
                "n_ref_clusters": len(ref_centroids),
                "n_method_clusters": len(m_centroids),
                "weighted_mean_dist": round(float((dists_arr * weights).sum()), 4),
                "unweighted_mean_dist": round(float(dists_arr.mean()), 4),
                "median_dist": round(float(np.median(dists_arr)), 4),
                "min_dist": round(float(dists_arr.min()), 4),
                "max_dist": round(float(dists_arr.max()), 4),
            })

    detail_df = pd.DataFrame(detail_rows)
    detail_df.to_csv(TABLES / "centroid_distance.csv", index=False)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(TABLES / "centroid_distance_summary.csv", index=False)

    # ---- best resolution per method
    print("\n" + "=" * 95)
    print("Best (ref_res, method_res) per method (minimizing weighted mean distance)")
    print("=" * 95)
    print(f"{'Method':<30} {'Ref res':>7} {'Meth res':>8} {'Ref cl':>6} {'Meth cl':>7} "
          f"{'Wt mean':>8} {'Mean':>8} {'Median':>8}")
    print("-" * 90)
    for key, label in available:
        sub = summary_df[summary_df["method"] == label]
        best_idx = sub["weighted_mean_dist"].idxmin()
        b = sub.loc[best_idx]
        print(f"{label:<30} {b['ref_res']:>7.1f} {b['method_res']:>8.1f} "
              f"{b['n_ref_clusters']:>6.0f} {b['n_method_clusters']:>7.0f} "
              f"{b['weighted_mean_dist']:>8.4f} {b['unweighted_mean_dist']:>8.4f} "
              f"{b['median_dist']:>8.4f}")

    # ---- at fixed ref_res=1.0, best method_res per method
    print("\n" + "=" * 95)
    print("Best method_res per method at ref_res=1.0")
    print("=" * 95)
    print(f"{'Method':<30} {'Meth res':>8} {'Meth cl':>7} "
          f"{'Wt mean':>8} {'Mean':>8}")
    print("-" * 65)
    for key, label in available:
        sub = summary_df[(summary_df["method"] == label) &
                         (summary_df["ref_res"] == 1.0)]
        best_idx = sub["weighted_mean_dist"].idxmin()
        b = sub.loc[best_idx]
        print(f"{label:<30} {b['method_res']:>8.1f} {b['n_method_clusters']:>7.0f} "
              f"{b['weighted_mean_dist']:>8.4f} {b['unweighted_mean_dist']:>8.4f}")

    print("\nSaved centroid_distance.csv and centroid_distance_summary.csv")

    # ================================================================
    # Cell type centroid comparison — all methods
    # ================================================================
    print("\n" + "=" * 80)
    print("Cell type centroid distance (marker-based annotation, all methods)")
    print("=" * 80)

    print("\nAnnotating reference cells...")
    ref_ct = score_and_assign(ref, set(ref.var_names))

    ct_rows = []
    for key, label in available:
        print(f"\n  {label}...")
        m_ct = score_and_assign(adatas_raw[label], set(adatas_raw[label].var_names))
        Z_m = projections[label]

        shared_types = sorted(set(ref_ct.unique()) & set(m_ct.unique()))
        for ct in shared_types:
            ref_mask = (ref_ct == ct).values
            m_mask = (m_ct == ct).values
            n_ref = int(ref_mask.sum())
            n_m = int(m_mask.sum())
            if n_ref < 5 or n_m < 5:
                continue
            ref_centroid = Z_ref[ref_mask].mean(axis=0)
            m_centroid = Z_m[m_mask].mean(axis=0)
            dist = float(np.linalg.norm(m_centroid - ref_centroid))
            ct_rows.append({
                "method": label,
                "cell_type": ct,
                "n_ref": n_ref,
                "n_method": n_m,
                "distance": round(dist, 4),
            })

    ct_df = pd.DataFrame(ct_rows)
    ct_df.to_csv(TABLES / "celltype_centroid_distance.csv", index=False)

    # Summary: mean cell type centroid distance per method
    print("\n" + "=" * 60)
    print("Mean cell type centroid distance per method")
    print("=" * 60)
    print(f"{'Method':<30} {'Mean dist':>10} {'Median':>8} {'Max':>8}")
    print("-" * 60)
    for key, label in available:
        sub = ct_df[ct_df["method"] == label]
        if sub.empty:
            continue
        print(f"{label:<30} {sub['distance'].mean():>10.4f} "
              f"{sub['distance'].median():>8.4f} {sub['distance'].max():>8.4f}")

    # Per cell type across methods
    print("\n" + "=" * 80)
    cell_types = sorted(ct_df["cell_type"].unique())
    method_labels = [l for _, l in available]
    header = f"{'Cell type':<25}" + "".join(f"{l:>12}" for l in method_labels)
    print(header)
    print("-" * len(header))
    for ct in cell_types:
        row = f"{ct:<25}"
        for label in method_labels:
            sub = ct_df[(ct_df["method"] == label) & (ct_df["cell_type"] == ct)]
            if sub.empty:
                row += f"{'—':>12}"
            else:
                row += f"{sub['distance'].values[0]:>12.2f}"
        print(row)

    print("\nSaved celltype_centroid_distance.csv")


if __name__ == "__main__":
    main()
