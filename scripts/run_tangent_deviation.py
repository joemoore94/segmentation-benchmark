"""Cluster centroid distance: segmentation methods vs scRNA-seq and 10x native.

At each Leiden resolution, clusters all methods in a shared PCA space (fit on
scRNA-seq reference), Hungarian-matches cluster centroids one-to-one, and
reports matched centroid distances. Also computes cell-level ARI vs 10x native
and pairwise between methods.

Reads:  data/reference/scrna_3p_filtered_feature_bc_matrix.h5
        data/processed/roi/adata_*.h5ad
Writes: results/tables/centroid_distance.csv
        results/tables/centroid_distance_pairwise.csv
        results/tables/celltype_centroid_distance.csv

Usage::

    conda run -n segbench python scripts/run_tangent_deviation.py
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score

from segbench.constants import METHOD_LABELS

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


def cluster_centroids(Z: np.ndarray, labels: np.ndarray) -> tuple[list[str], np.ndarray]:
    ids = sorted(set(labels), key=lambda c: int(c))
    centers = np.array([Z[labels == cl].mean(axis=0) for cl in ids])
    return ids, centers


def hungarian_centroid_dist(centers_a: np.ndarray, centers_b: np.ndarray):
    """Hungarian-match centroids, return matched distances and assignment."""
    n_a, n_b = len(centers_a), len(centers_b)
    cost = np.zeros((n_a, n_b))
    for i in range(n_a):
        cost[i] = np.linalg.norm(centers_b - centers_a[i], axis=1)
    row_ind, col_ind = linear_sum_assignment(cost)
    matched_dists = cost[row_ind, col_ind]
    return matched_dists, row_ind, col_ind


def spatial_match(adata_a, adata_b, max_dist=15.0):
    """Mutual nearest-neighbor spatial matching. Returns matched index arrays."""
    xy_a = np.column_stack([adata_a.obs["centroid_x"], adata_a.obs["centroid_y"]])
    xy_b = np.column_stack([adata_b.obs["centroid_x"], adata_b.obs["centroid_y"]])
    tree_a = cKDTree(xy_a)
    tree_b = cKDTree(xy_b)
    dist_ab, idx_b = tree_b.query(xy_a)
    dist_ba, idx_a = tree_a.query(xy_b)
    matched_a, matched_b = [], []
    for i_a, (i_b, d) in enumerate(zip(idx_b, dist_ab)):
        if d <= max_dist and idx_a[i_b] == i_a:
            matched_a.append(i_a)
            matched_b.append(i_b)
    return np.array(matched_a), np.array(matched_b)


def match_and_ari(adata_a, adata_b, labels_a, labels_b, max_dist=15.0):
    """Match cells by spatial centroid, compute ARI."""
    matched_a, matched_b = spatial_match(adata_a, adata_b, max_dist)
    if len(matched_a) < 100:
        return float("nan"), 0
    return adjusted_rand_score(labels_a[matched_a], labels_b[matched_b]), len(matched_a)


def argmax_centroid_dist(
    Z_a, labels_a, Z_b, labels_b, matched_a, matched_b,
):
    """Argmax-match clusters via cell overlap, return per-cluster centroid distances.

    Each cluster in B is mapped to whichever A cluster contains the plurality
    of its matched cells (many-to-one). Returns distances for each B cluster.
    """
    ids_a, centers_a = cluster_centroids(Z_a, labels_a)
    ids_b, centers_b = cluster_centroids(Z_b, labels_b)
    center_map_a = {cl: c for cl, c in zip(ids_a, centers_a)}

    la = labels_a[matched_a]
    lb = labels_b[matched_b]
    confusion = pd.crosstab(pd.Series(la, name="a"), pd.Series(lb, name="b"))

    dists = []
    for i, cl_b in enumerate(ids_b):
        if cl_b in confusion.columns:
            best_a = confusion[cl_b].idxmax()
            d = float(np.linalg.norm(centers_b[i] - center_map_a[best_a]))
        else:
            d = float("nan")
        dists.append(d)
    return np.array(dists)


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

    # ---- load reference
    print("Loading scRNA-seq reference...")
    ref = sc.read_10x_h5(str(REF_PATH))
    ref.var_names_make_unique()
    sc.pp.filter_cells(ref, min_genes=200)
    sc.pp.filter_genes(ref, min_cells=3)
    ref = ref[ref.obs["n_genes"] < 6000].copy()
    print(f"  {ref.n_obs} cells")

    # ---- load methods
    print("Loading segmentation methods...")
    adatas: dict[str, ad.AnnData] = {}
    available: list[tuple[str, str]] = []
    for key, fname in METHODS:
        path = ROI_DIR / fname
        if not path.exists():
            continue
        label = METHOD_LABELS[key]
        adatas[label] = ad.read_h5ad(path)
        available.append((key, label))
        print(f"  {label}: {adatas[label].n_obs} cells")

    # ---- shared PCA
    shared_genes = sorted(set(adatas[available[0][1]].var_names) & set(ref.var_names))
    print(f"\nShared genes: {len(shared_genes)}")
    X_ref = normalize_log(ref[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    Z_ref = pca.fit_transform(X_ref)
    print(f"PCA: {N_PCS} PCs, {pca.explained_variance_ratio_.sum():.1%} variance")

    projections: dict[str, np.ndarray] = {}
    for key, label in available:
        projections[label] = pca.transform(
            normalize_log(adatas[label][:, shared_genes].X)
        )

    # ---- build neighbor graphs + cluster at all resolutions
    print("\nClustering...")
    ref_ad = ad.AnnData(X=Z_ref, obs=pd.DataFrame(index=ref.obs_names))
    sc.pp.neighbors(ref_ad, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)

    method_ads: dict[str, ad.AnnData] = {}
    for key, label in available:
        method_ads[label] = ad.AnnData(
            X=projections[label],
            obs=pd.DataFrame(index=adatas[label].obs_names),
        )
        sc.pp.neighbors(method_ads[label], use_rep="X", n_neighbors=15,
                        random_state=RANDOM_STATE)

    ref_cl: dict[float, np.ndarray] = {}
    method_cl: dict[str, dict[float, np.ndarray]] = {l: {} for _, l in available}
    for res in RESOLUTIONS:
        sc.tl.leiden(ref_ad, resolution=res, random_state=RANDOM_STATE, flavor="igraph")
        ref_cl[res] = ref_ad.obs["leiden"].values.copy()
        for key, label in available:
            sc.tl.leiden(method_ads[label], resolution=res,
                         random_state=RANDOM_STATE, flavor="igraph")
            method_cl[label][res] = method_ads[label].obs["leiden"].values.copy()

    method_labels = [l for _, l in available]

    # ==================================================================
    # 1. Method vs scRNA-seq reference — centroid distances
    # ==================================================================
    print("\n" + "=" * 95)
    print("1. Centroid distance to scRNA-seq reference (Hungarian matching, same resolution)")
    print("=" * 95)

    vs_ref_rows = []
    header = f"{'Method':<30}" + "".join(f"{r:>6}" for r in RESOLUTIONS)
    print(header)
    print("-" * len(header))
    for key, label in available:
        vals = []
        for res in RESOLUTIONS:
            _, ref_centers = cluster_centroids(Z_ref, ref_cl[res])
            _, m_centers = cluster_centroids(projections[label], method_cl[label][res])
            dists, _, _ = hungarian_centroid_dist(m_centers, ref_centers)
            n_matched = len(dists)
            sizes = np.array([int((method_cl[label][res] == cl).sum())
                              for cl in sorted(set(method_cl[label][res]), key=int)])
            weights = sizes / sizes.sum()
            wt_mean = float((dists * weights[:n_matched]).sum()) if n_matched == len(sizes) else float(dists.mean())
            mean_d = float(dists.mean())
            vals.append(mean_d)
            vs_ref_rows.append({
                "method": label, "resolution": res,
                "n_ref_cl": len(ref_centers), "n_method_cl": len(m_centers),
                "n_matched": n_matched,
                "mean_matched_dist": round(mean_d, 4),
                "median_matched_dist": round(float(np.median(dists)), 4),
                "max_matched_dist": round(float(dists.max()), 4),
            })
        print(f"{label:<30}" + "".join(f"{v:>6.2f}" for v in vals))

    # ==================================================================
    # 2. Method vs 10x native — centroid distances + ARI
    # ==================================================================
    tenx_label = "10x native"

    # Precompute spatial matches (once per method pair, reused across resolutions)
    print("\nPrecomputing spatial cell matches vs 10x native...")
    spatial_matches: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for key, label in available:
        if label == tenx_label:
            continue
        ma, mb = spatial_match(adatas[tenx_label], adatas[label])
        spatial_matches[label] = (ma, mb)
        print(f"  {label}: {len(ma)} matched cells")

    print("\n" + "=" * 95)
    print("2a. Centroid distance to 10x native (Hungarian matching)")
    print("=" * 95)

    vs_10x_rows = []
    print(header)
    print("-" * len(header))
    for key, label in available:
        vals = []
        for res in RESOLUTIONS:
            _, tenx_centers = cluster_centroids(projections[tenx_label], method_cl[tenx_label][res])
            _, m_centers = cluster_centroids(projections[label], method_cl[label][res])
            dists, _, _ = hungarian_centroid_dist(m_centers, tenx_centers)
            mean_d = float(dists.mean())
            vals.append(mean_d)
            vs_10x_rows.append({
                "method": label, "resolution": res,
                "matching": "hungarian",
                "n_10x_cl": len(tenx_centers), "n_method_cl": len(m_centers),
                "mean_matched_dist": round(mean_d, 4),
            })
        print(f"{label:<30}" + "".join(f"{v:>6.2f}" for v in vals))

    print(f"\n{'2b. Centroid distance to 10x native (argmax matching)':<60}")
    print(header)
    print("-" * len(header))
    for key, label in available:
        vals = []
        for res in RESOLUTIONS:
            if label == tenx_label:
                vals.append(0.0)
                vs_10x_rows.append({
                    "method": label, "resolution": res,
                    "matching": "argmax",
                    "n_10x_cl": len(set(method_cl[tenx_label][res])),
                    "n_method_cl": len(set(method_cl[label][res])),
                    "mean_matched_dist": 0.0,
                })
                continue
            ma, mb = spatial_matches[label]
            dists = argmax_centroid_dist(
                projections[tenx_label], method_cl[tenx_label][res],
                projections[label], method_cl[label][res],
                ma, mb,
            )
            valid = dists[~np.isnan(dists)]
            mean_d = float(valid.mean()) if len(valid) else float("nan")
            vals.append(mean_d)
            vs_10x_rows.append({
                "method": label, "resolution": res,
                "matching": "argmax",
                "n_10x_cl": len(set(method_cl[tenx_label][res])),
                "n_method_cl": len(set(method_cl[label][res])),
                "mean_matched_dist": round(mean_d, 4),
            })
        print(f"{label:<30}" + "".join(f"{v:>6.2f}" for v in vals))

    print(f"\n{'2c. ARI vs 10x native (cell-level spatial matching)':<60}")
    print(header)
    print("-" * len(header))
    for key, label in available:
        vals = []
        for res in RESOLUTIONS:
            if label == tenx_label:
                vals.append(1.0)
                continue
            ma, mb = spatial_matches[label]
            la = method_cl[tenx_label][res][ma]
            lb = method_cl[label][res][mb]
            ari = adjusted_rand_score(la, lb) if len(ma) >= 100 else float("nan")
            vals.append(ari)
        print(f"{label:<30}" + "".join(f"{v:>6.3f}" for v in vals))

    # ==================================================================
    # 3. Pairwise method centroid distances
    # ==================================================================
    print("\n" + "=" * 95)
    print("3. Pairwise centroid distances (resolution 1.0, Hungarian matching)")
    print("=" * 95)

    pairwise_rows = []
    res_show = 1.0
    pair_header = f"{'':>30}" + "".join(f"{l:>12}" for l in method_labels)
    print(pair_header)
    print("-" * len(pair_header))
    pair_matrix = np.zeros((len(method_labels), len(method_labels)))
    for i, la in enumerate(method_labels):
        for j, lb in enumerate(method_labels):
            if i == j:
                pair_matrix[i, j] = 0.0
                continue
            if j < i:
                pair_matrix[i, j] = pair_matrix[j, i]
                continue
            _, ca = cluster_centroids(projections[la], method_cl[la][res_show])
            _, cb = cluster_centroids(projections[lb], method_cl[lb][res_show])
            dists, _, _ = hungarian_centroid_dist(ca, cb)
            pair_matrix[i, j] = pair_matrix[j, i] = float(dists.mean())
            pairwise_rows.append({
                "method_a": la, "method_b": lb, "resolution": res_show,
                "mean_matched_dist": round(float(dists.mean()), 4),
            })

    for i, la in enumerate(method_labels):
        row = f"{la:>30}"
        for j in range(len(method_labels)):
            if i == j:
                row += f"{'—':>12}"
            else:
                row += f"{pair_matrix[i,j]:>12.2f}"
        print(row)

    # ---- save
    pd.DataFrame(vs_ref_rows).to_csv(TABLES / "centroid_distance.csv", index=False)
    pd.DataFrame(vs_10x_rows).to_csv(TABLES / "centroid_distance_vs_10x.csv", index=False)
    pd.DataFrame(pairwise_rows).to_csv(TABLES / "centroid_distance_pairwise.csv", index=False)

    # ==================================================================
    # 4. Cell type centroid distance — all methods vs reference
    # ==================================================================
    print("\n" + "=" * 80)
    print("4. Cell type centroid distance (marker-based, all methods)")
    print("=" * 80)

    ref_ct = score_and_assign(ref, set(ref.var_names))
    ct_rows = []
    for key, label in available:
        m_ct = score_and_assign(adatas[label], set(adatas[label].var_names))
        Z_m = projections[label]
        for ct in sorted(set(ref_ct.unique()) & set(m_ct.unique())):
            ref_mask = (ref_ct == ct).values
            m_mask = (m_ct == ct).values
            if ref_mask.sum() < 5 or m_mask.sum() < 5:
                continue
            dist = float(np.linalg.norm(
                Z_m[m_mask].mean(axis=0) - Z_ref[ref_mask].mean(axis=0)
            ))
            ct_rows.append({
                "method": label, "cell_type": ct,
                "n_ref": int(ref_mask.sum()), "n_method": int(m_mask.sum()),
                "distance": round(dist, 4),
            })

    ct_df = pd.DataFrame(ct_rows)
    ct_df.to_csv(TABLES / "celltype_centroid_distance.csv", index=False)

    print(f"\n{'Method':<30} {'Mean':>8} {'Median':>8} {'Max':>8}")
    print("-" * 58)
    for _, label in available:
        sub = ct_df[ct_df["method"] == label]
        print(f"{label:<30} {sub['distance'].mean():>8.2f} "
              f"{sub['distance'].median():>8.2f} {sub['distance'].max():>8.2f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
