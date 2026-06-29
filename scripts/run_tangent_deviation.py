"""Cluster centroid distance: 10x native vs scRNA-seq reference.

For each Leiden resolution, clusters both the scRNA-seq reference and 10x
native in a shared PCA space (fit on the reference, 374 shared genes, 30 PCs),
computes cluster centroids, and measures the Euclidean distance from each 10x
native centroid to the nearest scRNA-seq centroid. This answers: does each
Xenium cluster correspond to a cell state present in the reference?

Reads:  data/reference/scrna_3p_filtered_feature_bc_matrix.h5
        data/processed/roi/adata_10x.h5ad
Writes: results/tables/centroid_distance.csv
        results/tables/centroid_distance_summary.csv

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

ROI_DIR  = Path("data/processed/roi")
REF_PATH = Path("data/reference/scrna_3p_filtered_feature_bc_matrix.h5")
TABLES   = Path("results/tables")

RESOLUTIONS = [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]
N_PCS = 30
RANDOM_STATE = 0


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

    # ---- load 10x native
    print("Loading 10x native...")
    xenium = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    print(f"  {xenium.n_obs} cells")

    # ---- shared gene set + PCA
    shared_genes = sorted(set(xenium.var_names) & set(ref.var_names))
    print(f"  Shared genes: {len(shared_genes)}")

    X_ref = normalize_log(ref[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    Z_ref = pca.fit_transform(X_ref)
    print(f"  PCA: {N_PCS} PCs, {pca.explained_variance_ratio_.sum():.1%} variance")

    X_xen = normalize_log(xenium[:, shared_genes].X)
    Z_xen = pca.transform(X_xen)

    # ---- build neighbor graphs (once per dataset, reuse for all resolutions)
    ref_ad = ad.AnnData(X=Z_ref, obs=pd.DataFrame(index=ref.obs_names))
    sc.pp.neighbors(ref_ad, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)

    xen_ad = ad.AnnData(X=Z_xen, obs=pd.DataFrame(index=xenium.obs_names))
    sc.pp.neighbors(xen_ad, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)

    # ---- sweep resolutions
    print("\nClustering at each resolution...")
    ref_clusterings: dict[float, np.ndarray] = {}
    xen_clusterings: dict[float, np.ndarray] = {}
    for res in RESOLUTIONS:
        sc.tl.leiden(ref_ad, resolution=res, random_state=RANDOM_STATE, flavor="igraph")
        ref_clusterings[res] = ref_ad.obs["leiden"].values.copy()

        sc.tl.leiden(xen_ad, resolution=res, random_state=RANDOM_STATE, flavor="igraph")
        xen_clusterings[res] = xen_ad.obs["leiden"].values.copy()

        n_ref = len(set(ref_clusterings[res]))
        n_xen = len(set(xen_clusterings[res]))
        print(f"  res {res}: ref={n_ref} clusters, xenium={n_xen} clusters")

    # ---- compute centroid distances for all (ref_res, xen_res) pairs
    print("\nComputing centroid distances...")
    detail_rows = []
    summary_rows = []

    for ref_res, xen_res in product(RESOLUTIONS, RESOLUTIONS):
        ref_centroids = cluster_centroids(Z_ref, ref_clusterings[ref_res])
        xen_centroids = cluster_centroids(Z_xen, xen_clusterings[xen_res])
        xen_labels = xen_clusterings[xen_res]

        ref_centers = np.array(list(ref_centroids.values()))
        ref_keys = list(ref_centroids.keys())

        distances = []
        for cl, centroid in xen_centroids.items():
            dists = np.linalg.norm(ref_centers - centroid, axis=1)
            nearest_idx = int(np.argmin(dists))
            nearest_dist = float(dists[nearest_idx])
            n_cells = int((xen_labels == cl).sum())

            detail_rows.append({
                "ref_res": ref_res,
                "xen_res": xen_res,
                "xen_cluster": cl,
                "n_cells": n_cells,
                "nearest_ref_cluster": ref_keys[nearest_idx],
                "distance": round(nearest_dist, 4),
            })
            distances.append((nearest_dist, n_cells))

        dists_arr = np.array([d for d, _ in distances])
        sizes = np.array([n for _, n in distances])
        weights = sizes / sizes.sum()

        summary_rows.append({
            "ref_res": ref_res,
            "xen_res": xen_res,
            "n_ref_clusters": len(ref_centroids),
            "n_xen_clusters": len(xen_centroids),
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

    # ---- print summary
    print("\n" + "=" * 95)
    print("Centroid distance: 10x native clusters → nearest scRNA-seq cluster (30-PC Euclidean)")
    print("=" * 95)
    print(f"{'Ref res':>7}  {'Xen res':>7}  {'Ref cl':>6}  {'Xen cl':>6}  "
          f"{'Wt mean':>8}  {'Mean':>8}  {'Median':>8}  {'Min':>8}  {'Max':>8}")
    print("-" * 85)
    for _, r in summary_df.iterrows():
        print(f"{r['ref_res']:>7.1f}  {r['xen_res']:>7.1f}  "
              f"{r['n_ref_clusters']:>6.0f}  {r['n_xen_clusters']:>6.0f}  "
              f"{r['weighted_mean_dist']:>8.4f}  {r['unweighted_mean_dist']:>8.4f}  "
              f"{r['median_dist']:>8.4f}  {r['min_dist']:>8.4f}  {r['max_dist']:>8.4f}")

    # Best combination
    best = summary_df.loc[summary_df["weighted_mean_dist"].idxmin()]
    print(f"\nBest: ref_res={best['ref_res']}, xen_res={best['xen_res']}, "
          f"weighted_mean_dist={best['weighted_mean_dist']:.4f}")

    print("\nSaved centroid_distance.csv and centroid_distance_summary.csv")

    # ================================================================
    # Cell type centroid comparison
    # ================================================================
    # Score every cell in both datasets with the same canonical marker
    # panels, assign cell type by highest score, then compare per-type
    # centroids in the shared PCA space.

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

    print("\n" + "=" * 70)
    print("Cell type centroid comparison (marker-based annotation)")
    print("=" * 70)

    def score_and_assign(adata_raw, gene_universe, label):
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
        score_df = pd.DataFrame(scores, index=a.obs_names)
        assignments = score_df.idxmax(axis=1)
        print(f"  {label}: annotated {len(assignments)} cells across "
              f"{assignments.nunique()} types")
        for ct in sorted(assignments.unique()):
            print(f"    {ct}: {(assignments == ct).sum()}")
        return assignments

    ref_genes = set(ref.var_names)
    xen_genes = set(xenium.var_names)
    gene_universe_ref = ref_genes
    gene_universe_xen = xen_genes

    print("\nAnnotating reference cells...")
    ref_ct = score_and_assign(ref, gene_universe_ref, "scRNA-seq")
    print("\nAnnotating 10x native cells...")
    xen_ct = score_and_assign(xenium, gene_universe_xen, "10x native")

    # Compute per-cell-type centroids in the shared PCA space
    shared_types = sorted(set(ref_ct.unique()) & set(xen_ct.unique()))

    ct_rows = []
    print(f"\n{'Cell type':<25} {'Ref cells':>9} {'Xen cells':>9} {'Distance':>10}")
    print("-" * 58)
    for ct in shared_types:
        ref_mask = (ref_ct == ct).values
        xen_mask = (xen_ct == ct).values
        n_ref = int(ref_mask.sum())
        n_xen = int(xen_mask.sum())
        if n_ref < 5 or n_xen < 5:
            continue
        ref_centroid = Z_ref[ref_mask].mean(axis=0)
        xen_centroid = Z_xen[xen_mask].mean(axis=0)
        dist = float(np.linalg.norm(xen_centroid - ref_centroid))
        print(f"{ct:<25} {n_ref:>9} {n_xen:>9} {dist:>10.4f}")
        ct_rows.append({
            "cell_type": ct, "n_ref": n_ref, "n_xen": n_xen,
            "distance": round(dist, 4),
        })

    ct_df = pd.DataFrame(ct_rows)
    ct_df.to_csv(TABLES / "celltype_centroid_distance.csv", index=False)
    if len(ct_df):
        print(f"\nMean cell type centroid distance: {ct_df['distance'].mean():.4f}")
    print("Saved celltype_centroid_distance.csv")


if __name__ == "__main__":
    main()
