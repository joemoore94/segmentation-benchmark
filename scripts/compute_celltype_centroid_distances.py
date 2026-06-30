"""Compute cell-type centroid distances for all methods vs scRNA-seq reference.

A fast version of the section 4 computation from run_tangent_deviation.py,
extended to all 29 methods. Skips the expensive Leiden-resolution sweep and
pairwise matrix. Writes (and overwrites) celltype_centroid_distance.csv.

Reads:  data/reference/scrna_3p_filtered_feature_bc_matrix.h5
        data/processed/roi/adata_*.h5ad
Writes: results/tables/celltype_centroid_distance.csv
        results/tables/celltype_centroid_distance_summary.csv

Usage:
    conda run -n segbench python scripts/compute_celltype_centroid_distances.py
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import PCA

REPO = Path(__file__).resolve().parents[1]
ROI_DIR = REPO / "data" / "processed" / "roi"
REF_PATH = REPO / "data" / "reference" / "scrna_3p_filtered_feature_bc_matrix.h5"
TABLES = REPO / "results" / "tables"

N_PCS = 30
RANDOM_STATE = 0

METHODS = [
    ("10x_native",                   "adata_10x.h5ad"),
    # Nuclear-only detectors
    ("cellpose",                     "adata_cellpose.h5ad"),
    ("stardist",                     "adata_stardist.h5ad"),
    ("mesmer",                       "adata_mesmer.h5ad"),
    ("10x_ranger",                   "adata_10x_ranger.h5ad"),
    # Voronoi expansion
    ("voronoi",                      "adata_voronoi.h5ad"),
    ("voronoi_stardist",             "adata_voronoi_stardist.h5ad"),
    ("voronoi_mesmer",               "adata_voronoi_mesmer.h5ad"),
    ("voronoi_10x_ranger",           "adata_voronoi_10x_ranger.h5ad"),
    # Geometric expansion 10µm
    ("cellpose_exp10um",             "adata_cellpose_exp10um.h5ad"),
    ("stardist_exp10um",             "adata_stardist_exp10um.h5ad"),
    ("mesmer_exp10um",               "adata_mesmer_exp10um.h5ad"),
    ("10x_ranger_exp10um",           "adata_10x_ranger_exp10um.h5ad"),
    # Geometric expansion 20µm
    ("cellpose_exp20um",             "adata_cellpose_exp20um.h5ad"),
    ("stardist_exp20um",             "adata_stardist_exp20um.h5ad"),
    ("mesmer_exp20um",               "adata_mesmer_exp20um.h5ad"),
    ("10x_ranger_exp20um",           "adata_10x_ranger_exp20um.h5ad"),
    # Watershed expansion
    ("watershed_10x",                "adata_watershed_10x.h5ad"),
    ("watershed_stardist",           "adata_watershed_stardist.h5ad"),
    ("watershed_mesmer",             "adata_watershed_mesmer.h5ad"),
    # Baysor PSC sweep (CellPose prior)
    ("baysor",                       "adata_baysor.h5ad"),
    ("baysor_prior",                 "adata_baysor_prior.h5ad"),
    ("baysor_prior_c05",             "adata_baysor_prior_c05.h5ad"),
    ("baysor_prior_c08",             "adata_baysor_prior_c08.h5ad"),
    ("baysor_prior_c10",             "adata_baysor_prior_c10.h5ad"),
    # Baysor other detectors at PSC=1.0
    ("baysor_stardist_prior_c10",    "adata_baysor_stardist_prior_c10.h5ad"),
    ("baysor_mesmer_prior_c10",      "adata_baysor_mesmer_prior_c10.h5ad"),
    ("baysor_10x_ranger_prior_c10",  "adata_baysor_10x_ranger_prior_c10.h5ad"),
    # Whole-cell NN
    ("cellpose_cyto3",               "adata_cellpose_cyto3.h5ad"),
    ("cellpose_cyto3_eosin",         "adata_cellpose_cyto3_eosin.h5ad"),
    ("cellpose_cyto3_density",       "adata_cellpose_cyto3_density.h5ad"),
    ("mesmer_wholecell_eosin",       "adata_mesmer_wholecell_eosin.h5ad"),
    ("mesmer_wholecell_density",     "adata_mesmer_wholecell_density.h5ad"),
]

try:
    from segbench.constants import METHOD_LABELS
except ImportError:
    METHOD_LABELS = {k: k for k, _ in METHODS}

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


def score_and_assign(adata_raw, gene_universe: set[str]) -> pd.Series:
    a = adata_raw.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    scores: dict[str, np.ndarray] = {}
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

    print("Loading scRNA-seq reference...")
    ref = sc.read_10x_h5(str(REF_PATH))
    ref.var_names_make_unique()
    sc.pp.filter_cells(ref, min_genes=200)
    sc.pp.filter_genes(ref, min_cells=3)
    ref = ref[ref.obs["n_genes"] < 6000].copy()
    print(f"  {ref.n_obs} cells, {ref.n_vars} genes after QC")

    print("\nLoading Xenium adatas...")
    adatas: dict[str, ad.AnnData] = {}
    available: list[tuple[str, str]] = []
    for key, fname in METHODS:
        path = ROI_DIR / fname
        if not path.exists():
            print(f"  {key}: skipped (not found)")
            continue
        label = METHOD_LABELS.get(key, key)
        adatas[label] = ad.read_h5ad(path)
        available.append((key, label))
        print(f"  {label}: {adatas[label].n_obs} cells")

    # Shared genes: reference ∩ all adatas (prevents KeyError for genes absent in new adatas)
    ref_genes = set(ref.var_names)
    all_adata_genes = set.intersection(*[set(adatas[lbl].var_names) for _, lbl in available])
    shared_genes = sorted(ref_genes & all_adata_genes)
    print(f"\nShared genes (ref ∩ all methods): {len(shared_genes)}")

    # Fit PCA on reference
    print(f"Fitting PCA on reference ({len(shared_genes)} genes, {N_PCS} PCs)...")
    X_ref = normalize_log(ref[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    Z_ref = pca.fit_transform(X_ref)
    print(f"  {pca.explained_variance_ratio_.sum():.1%} variance explained")

    # Project all methods
    print("\nProjecting methods into reference PCA space...")
    projections: dict[str, np.ndarray] = {}
    for key, label in available:
        projections[label] = pca.transform(
            normalize_log(adatas[label][:, shared_genes].X)
        )
        print(f"  {label}: done")

    # Annotate reference cell types
    print("\nAnnotating reference cell types...")
    ref_ct = score_and_assign(ref, set(ref.var_names))
    print(f"  Reference: {ref_ct.value_counts().to_dict()}")

    # Cell type centroid distances
    print("\nComputing cell-type centroid distances...")
    ct_rows = []
    for key, label in available:
        m_ct = score_and_assign(adatas[label], set(adatas[label].var_names))
        Z_m = projections[label]
        shared_cts = sorted(set(ref_ct.unique()) & set(m_ct.unique()))
        for ct in shared_cts:
            ref_mask = (ref_ct == ct).values
            m_mask = (m_ct == ct).values
            if ref_mask.sum() < 5 or m_mask.sum() < 5:
                continue
            dist = float(np.linalg.norm(
                Z_m[m_mask].mean(axis=0) - Z_ref[ref_mask].mean(axis=0)
            ))
            ct_rows.append({
                "method": label,
                "cell_type": ct,
                "n_ref": int(ref_mask.sum()),
                "n_method": int(m_mask.sum()),
                "distance": round(dist, 4),
            })

    ct_df = pd.DataFrame(ct_rows)
    ct_df.to_csv(TABLES / "celltype_centroid_distance.csv", index=False)
    print(f"Saved celltype_centroid_distance.csv ({len(ct_df)} rows)")

    # Summary: mean distance per method
    summary_rows = []
    print(f"\n{'Method':<35} {'Mean':>6} {'Median':>8} {'Max':>6}")
    print("-" * 60)
    for key, label in available:
        sub = ct_df[ct_df["method"] == label]
        if sub.empty:
            continue
        row = {
            "method_key": key,
            "method_label": label,
            "mean_dist": round(float(sub["distance"].mean()), 4),
            "median_dist": round(float(sub["distance"].median()), 4),
            "max_dist": round(float(sub["distance"].max()), 4),
        }
        summary_rows.append(row)
        print(f"{label:<35} {row['mean_dist']:>6.2f} {row['median_dist']:>8.2f} {row['max_dist']:>6.2f}")

    pd.DataFrame(summary_rows).to_csv(
        TABLES / "celltype_centroid_distance_summary.csv", index=False
    )
    print(f"\nSaved celltype_centroid_distance_summary.csv")
    print("Done.")


if __name__ == "__main__":
    main()
