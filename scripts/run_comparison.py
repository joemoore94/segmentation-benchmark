"""Cross-method comparison anchored to 10x native (Xenium Ranger) segmentation.

10x native is the platform's own reference segmentation, making it the most
defensible anchor for benchmarking other methods: the question becomes "how
closely does method X reproduce the platform's own cell calls?" rather than
implicitly treating CellPose as ground truth.

Usage::

    conda run -n segbench python scripts/run_comparison.py
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad

from segbench.compare import (
    cell_count_summary,
    cell_type_agreement,
    cluster_cell_types,
    cluster_embedding,
    expression_correlation,
    match_cells_by_centroid,
    match_cluster_labels,
    match_cluster_labels_argmax,
    size_summary,
)
from segbench.constants import TOTAL_TRANSCRIPTS_FULL_ROI
from segbench.spatial import disagreement_spatial_structure, disagreement_table

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

MAX_MATCH_DIST = 10.0


def _save_disagreement(method_key, matches, labels_10x, labels_comp, adata_10x):
    """Save disagreement tables for both Hungarian and argmax cluster alignment."""
    for suffix, matcher_fn in [("", match_cluster_labels), ("_argmax", match_cluster_labels_argmax)]:
        labels_aligned = matcher_fn(labels_10x, labels_comp, matches)
        dt = disagreement_table(matches, labels_10x, labels_aligned, adata_10x)
        dt.to_csv(TABLES_DIR / f"disagreement_table_10x_{method_key}{suffix}.csv", index=False)
        spatial = disagreement_spatial_structure(dt, n_perm=9999, seed=0)
        with open(TABLES_DIR / f"disagreement_spatial_10x_{method_key}{suffix}.json", "w") as f:
            json.dump(spatial, f, indent=2)
        tag = "hungarian" if not suffix else "argmax"
        print(f"  [{tag}] disagree={dt['disagree'].mean():.3f}  Moran's I={spatial.get('morans_i', 'N/A')}")


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adata_baysor_prior = ad.read_h5ad(ROI_DIR / "adata_baysor_prior.h5ad")
    adata_voronoi = ad.read_h5ad(ROI_DIR / "adata_voronoi.h5ad")
    adata_mesmer = ad.read_h5ad(ROI_DIR / "adata_mesmer.h5ad")
    adata_voronoi_mesmer = ad.read_h5ad(ROI_DIR / "adata_voronoi_mesmer.h5ad")
    adata_voronoi_stardist = ad.read_h5ad(ROI_DIR / "adata_voronoi_stardist.h5ad")
    adatas = {
        "10x_native": adata_10x,
        "cellpose": adata_cellpose,
        "baysor": adata_baysor,
        "stardist": adata_stardist,
        "baysor_prior": adata_baysor_prior,
        "voronoi": adata_voronoi,
        "mesmer": adata_mesmer,
        "voronoi_mesmer": adata_voronoi_mesmer,
        "voronoi_stardist": adata_voronoi_stardist,
    }

    for optional_key, optional_file in [
        ("baysor_prior_c05",          "adata_baysor_prior_c05.h5ad"),
        ("baysor_prior_c08",          "adata_baysor_prior_c08.h5ad"),
        ("baysor_prior_c10",          "adata_baysor_prior_c10.h5ad"),
        ("baysor_stardist_prior_c10", "adata_baysor_stardist_prior_c10.h5ad"),
        ("baysor_mesmer_prior_c10",      "adata_baysor_mesmer_prior_c10.h5ad"),
        ("baysor_10x_ranger_prior_c10", "adata_baysor_10x_ranger_prior_c10.h5ad"),
        ("voronoi_10x_ranger",          "adata_voronoi_10x_ranger.h5ad"),
        ("10x_ranger",                  "adata_10x_ranger.h5ad"),
        ("cellpose_exp10um",            "adata_cellpose_exp10um.h5ad"),
        ("cellpose_exp20um",            "adata_cellpose_exp20um.h5ad"),
        ("stardist_exp10um",            "adata_stardist_exp10um.h5ad"),
        ("stardist_exp20um",            "adata_stardist_exp20um.h5ad"),
        ("mesmer_exp10um",              "adata_mesmer_exp10um.h5ad"),
        ("mesmer_exp20um",              "adata_mesmer_exp20um.h5ad"),
        ("10x_ranger_exp10um",          "adata_10x_ranger_exp10um.h5ad"),
        ("10x_ranger_exp20um",          "adata_10x_ranger_exp20um.h5ad"),
        ("watershed_10x",              "adata_watershed_10x.h5ad"),
        ("watershed_stardist",         "adata_watershed_stardist.h5ad"),
        ("watershed_mesmer",           "adata_watershed_mesmer.h5ad"),
        ("cellpose_cyto3",             "adata_cellpose_cyto3.h5ad"),
        ("cellpose_cyto3_eosin",       "adata_cellpose_cyto3_eosin.h5ad"),
        ("cellpose_cyto3_density",     "adata_cellpose_cyto3_density.h5ad"),
        ("mesmer_wholecell_eosin",     "adata_mesmer_wholecell_eosin.h5ad"),
        ("mesmer_wholecell_density",   "adata_mesmer_wholecell_density.h5ad"),
        ("bidcell",                      "adata_bidcell.h5ad"),
        ("segger",                       "adata_segger.h5ad"),
    ]:
        p = ROI_DIR / optional_file
        if p.exists():
            adatas[optional_key] = ad.read_h5ad(p)

    counts = cell_count_summary(adatas)
    counts["transcript_capture_rate"] = counts["total_transcripts"] / TOTAL_TRANSCRIPTS_FULL_ROI
    counts.to_csv(TABLES_DIR / "cell_counts.csv")
    print("=== Cell counts + QC (full 2mm x 2mm ROI) ===")
    print(counts)

    sizes = size_summary(adatas)
    sizes.to_csv(TABLES_DIR / "size_summary.csv")
    print("\n=== Size summary ===")
    print(sizes)

    print("\n=== Clustering cell types (10x native, anchor) ===")
    labels_10x = cluster_cell_types(adata_10x, seed=0)
    print(f"10x_native: {labels_10x.nunique()} clusters")
    cluster_embedding(adata_10x, seed=0).to_csv(TABLES_DIR / "embedding_10x_native.csv")

    # 10x native vs. CellPose
    matches_cellpose = match_cells_by_centroid(adata_10x, adata_cellpose, max_dist=MAX_MATCH_DIST)
    matches_cellpose.to_csv(TABLES_DIR / "matches_10x_cellpose.csv", index=False)
    print(f"\n=== 10x native vs. CellPose ({len(matches_cellpose)} pairs) ===")

    corr_cellpose = expression_correlation(adata_10x, adata_cellpose, matches_cellpose)
    corr_cellpose.to_csv(TABLES_DIR / "expression_correlation_10x_cellpose.csv", index=False)
    print(corr_cellpose["correlation"].describe())

    print("\n=== Clustering cell types (CellPose) ===")
    labels_cellpose = cluster_cell_types(adata_cellpose, seed=0)
    print(f"cellpose: {labels_cellpose.nunique()} clusters")
    cluster_embedding(adata_cellpose, seed=0).to_csv(TABLES_DIR / "embedding_cellpose.csv")

    agreement_cellpose = cell_type_agreement(labels_10x, labels_cellpose, matches_cellpose)
    agreement_cellpose["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_cellpose.csv")
    print(f"ARI = {agreement_cellpose['ari']:.4f}, n_matched = {agreement_cellpose['n_matched']}")

    _save_disagreement("cellpose", matches_cellpose, labels_10x, labels_cellpose, adata_10x)

    # 10x native vs. Baysor
    matches_baysor = match_cells_by_centroid(adata_10x, adata_baysor, max_dist=MAX_MATCH_DIST)
    matches_baysor.to_csv(TABLES_DIR / "matches_10x_baysor.csv", index=False)
    print(f"\n=== 10x native vs. Baysor ({len(matches_baysor)} pairs) ===")

    corr_baysor = expression_correlation(adata_10x, adata_baysor, matches_baysor)
    corr_baysor.to_csv(TABLES_DIR / "expression_correlation_10x_baysor.csv", index=False)
    print(corr_baysor["correlation"].describe())

    print("\n=== Clustering cell types (Baysor) ===")
    labels_baysor = cluster_cell_types(adata_baysor, seed=0)
    print(f"baysor: {labels_baysor.nunique()} clusters")
    cluster_embedding(adata_baysor, seed=0).to_csv(TABLES_DIR / "embedding_baysor.csv")

    agreement_baysor = cell_type_agreement(labels_10x, labels_baysor, matches_baysor)
    agreement_baysor["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_baysor.csv")
    print(f"ARI = {agreement_baysor['ari']:.4f}, n_matched = {agreement_baysor['n_matched']}")

    _save_disagreement("baysor", matches_baysor, labels_10x, labels_baysor, adata_10x)

    # 10x native vs. StarDist
    matches_stardist = match_cells_by_centroid(adata_10x, adata_stardist, max_dist=MAX_MATCH_DIST)
    matches_stardist.to_csv(TABLES_DIR / "matches_10x_stardist.csv", index=False)
    print(f"\n=== 10x native vs. StarDist ({len(matches_stardist)} pairs) ===")

    corr_stardist = expression_correlation(adata_10x, adata_stardist, matches_stardist)
    corr_stardist.to_csv(TABLES_DIR / "expression_correlation_10x_stardist.csv", index=False)
    print(corr_stardist["correlation"].describe())

    print("\n=== Clustering cell types (StarDist) ===")
    labels_stardist = cluster_cell_types(adata_stardist, seed=0)
    print(f"stardist: {labels_stardist.nunique()} clusters")
    cluster_embedding(adata_stardist, seed=0).to_csv(TABLES_DIR / "embedding_stardist.csv")

    agreement_stardist = cell_type_agreement(labels_10x, labels_stardist, matches_stardist)
    agreement_stardist["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_stardist.csv")
    print(f"ARI = {agreement_stardist['ari']:.4f}, n_matched = {agreement_stardist['n_matched']}")

    _save_disagreement("stardist", matches_stardist, labels_10x, labels_stardist, adata_10x)

    # 10x native vs. Baysor (CellPose prior)
    matches_baysor_prior = match_cells_by_centroid(
        adata_10x, adata_baysor_prior, max_dist=MAX_MATCH_DIST
    )
    matches_baysor_prior.to_csv(TABLES_DIR / "matches_10x_baysor_prior.csv", index=False)
    print(f"\n=== 10x native vs. Baysor (prior) ({len(matches_baysor_prior)} pairs) ===")

    corr_baysor_prior = expression_correlation(adata_10x, adata_baysor_prior, matches_baysor_prior)
    corr_baysor_prior.to_csv(TABLES_DIR / "expression_correlation_10x_baysor_prior.csv", index=False)
    print(corr_baysor_prior["correlation"].describe())

    print("\n=== Clustering cell types (Baysor prior) ===")
    labels_baysor_prior = cluster_cell_types(adata_baysor_prior, seed=0)
    print(f"baysor_prior: {labels_baysor_prior.nunique()} clusters")
    cluster_embedding(adata_baysor_prior, seed=0).to_csv(TABLES_DIR / "embedding_baysor_prior.csv")

    agreement_baysor_prior = cell_type_agreement(labels_10x, labels_baysor_prior, matches_baysor_prior)
    agreement_baysor_prior["confusion"].to_csv(
        TABLES_DIR / "cell_type_confusion_10x_baysor_prior.csv"
    )
    print(f"ARI = {agreement_baysor_prior['ari']:.4f}, n_matched = {agreement_baysor_prior['n_matched']}")

    _save_disagreement("baysor_prior", matches_baysor_prior, labels_10x, labels_baysor_prior, adata_10x)

    # 10x native vs. Voronoi (CellPose nuclei, nearest-centroid transcript assignment)
    matches_voronoi = match_cells_by_centroid(adata_10x, adata_voronoi, max_dist=MAX_MATCH_DIST)
    matches_voronoi.to_csv(TABLES_DIR / "matches_10x_voronoi.csv", index=False)
    print(f"\n=== 10x native vs. Voronoi ({len(matches_voronoi)} pairs) ===")

    corr_voronoi = expression_correlation(adata_10x, adata_voronoi, matches_voronoi)
    corr_voronoi.to_csv(TABLES_DIR / "expression_correlation_10x_voronoi.csv", index=False)
    print(corr_voronoi["correlation"].describe())

    print("\n=== Clustering cell types (Voronoi) ===")
    labels_voronoi = cluster_cell_types(adata_voronoi, seed=0)
    print(f"voronoi: {labels_voronoi.nunique()} clusters")
    cluster_embedding(adata_voronoi, seed=0).to_csv(TABLES_DIR / "embedding_voronoi.csv")

    agreement_voronoi = cell_type_agreement(labels_10x, labels_voronoi, matches_voronoi)
    agreement_voronoi["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_voronoi.csv")
    print(f"ARI = {agreement_voronoi['ari']:.4f}, n_matched = {agreement_voronoi['n_matched']}")

    _save_disagreement("voronoi", matches_voronoi, labels_10x, labels_voronoi, adata_10x)

    # 10x native vs. Mesmer
    matches_mesmer = match_cells_by_centroid(adata_10x, adata_mesmer, max_dist=MAX_MATCH_DIST)
    matches_mesmer.to_csv(TABLES_DIR / "matches_10x_mesmer.csv", index=False)
    print(f"\n=== 10x native vs. Mesmer ({len(matches_mesmer)} pairs) ===")

    corr_mesmer = expression_correlation(adata_10x, adata_mesmer, matches_mesmer)
    corr_mesmer.to_csv(TABLES_DIR / "expression_correlation_10x_mesmer.csv", index=False)
    print(corr_mesmer["correlation"].describe())

    print("\n=== Clustering cell types (Mesmer) ===")
    labels_mesmer = cluster_cell_types(adata_mesmer, seed=0)
    print(f"mesmer: {labels_mesmer.nunique()} clusters")
    cluster_embedding(adata_mesmer, seed=0).to_csv(TABLES_DIR / "embedding_mesmer.csv")

    agreement_mesmer = cell_type_agreement(labels_10x, labels_mesmer, matches_mesmer)
    agreement_mesmer["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_mesmer.csv")
    print(f"ARI = {agreement_mesmer['ari']:.4f}, n_matched = {agreement_mesmer['n_matched']}")

    _save_disagreement("mesmer", matches_mesmer, labels_10x, labels_mesmer, adata_10x)

    # 10x native vs. Voronoi (Mesmer centroids)
    matches_voronoi_mesmer = match_cells_by_centroid(
        adata_10x, adata_voronoi_mesmer, max_dist=MAX_MATCH_DIST
    )
    matches_voronoi_mesmer.to_csv(TABLES_DIR / "matches_10x_voronoi_mesmer.csv", index=False)
    print(f"\n=== 10x native vs. Voronoi (Mesmer) ({len(matches_voronoi_mesmer)} pairs) ===")

    corr_voronoi_mesmer = expression_correlation(adata_10x, adata_voronoi_mesmer, matches_voronoi_mesmer)
    corr_voronoi_mesmer.to_csv(TABLES_DIR / "expression_correlation_10x_voronoi_mesmer.csv", index=False)
    print(corr_voronoi_mesmer["correlation"].describe())

    print("\n=== Clustering cell types (Voronoi-Mesmer) ===")
    labels_voronoi_mesmer = cluster_cell_types(adata_voronoi_mesmer, seed=0)
    print(f"voronoi_mesmer: {labels_voronoi_mesmer.nunique()} clusters")
    cluster_embedding(adata_voronoi_mesmer, seed=0).to_csv(
        TABLES_DIR / "embedding_voronoi_mesmer.csv"
    )

    agreement_voronoi_mesmer = cell_type_agreement(
        labels_10x, labels_voronoi_mesmer, matches_voronoi_mesmer
    )
    agreement_voronoi_mesmer["confusion"].to_csv(
        TABLES_DIR / "cell_type_confusion_10x_voronoi_mesmer.csv"
    )
    print(f"ARI = {agreement_voronoi_mesmer['ari']:.4f}, n_matched = {agreement_voronoi_mesmer['n_matched']}")

    _save_disagreement("voronoi_mesmer", matches_voronoi_mesmer, labels_10x, labels_voronoi_mesmer, adata_10x)

    # 10x native vs. Voronoi (StarDist centroids)
    matches_voronoi_stardist = match_cells_by_centroid(
        adata_10x, adata_voronoi_stardist, max_dist=MAX_MATCH_DIST
    )
    matches_voronoi_stardist.to_csv(TABLES_DIR / "matches_10x_voronoi_stardist.csv", index=False)
    print(f"\n=== 10x native vs. Voronoi (StarDist) ({len(matches_voronoi_stardist)} pairs) ===")

    corr_voronoi_stardist = expression_correlation(adata_10x, adata_voronoi_stardist, matches_voronoi_stardist)
    corr_voronoi_stardist.to_csv(TABLES_DIR / "expression_correlation_10x_voronoi_stardist.csv", index=False)
    print(corr_voronoi_stardist["correlation"].describe())

    print("\n=== Clustering cell types (Voronoi-StarDist) ===")
    labels_voronoi_stardist = cluster_cell_types(adata_voronoi_stardist, seed=0)
    print(f"voronoi_stardist: {labels_voronoi_stardist.nunique()} clusters")
    cluster_embedding(adata_voronoi_stardist, seed=0).to_csv(
        TABLES_DIR / "embedding_voronoi_stardist.csv"
    )

    agreement_voronoi_stardist = cell_type_agreement(
        labels_10x, labels_voronoi_stardist, matches_voronoi_stardist
    )
    agreement_voronoi_stardist["confusion"].to_csv(
        TABLES_DIR / "cell_type_confusion_10x_voronoi_stardist.csv"
    )
    print(f"ARI = {agreement_voronoi_stardist['ari']:.4f}, n_matched = {agreement_voronoi_stardist['n_matched']}")

    _save_disagreement("voronoi_stardist", matches_voronoi_stardist, labels_10x, labels_voronoi_stardist, adata_10x)

    # Generic comparison loop for optional methods
    from segbench.constants import METHOD_LABELS
    optional_methods = [
        "baysor_prior_c05", "baysor_prior_c08", "baysor_prior_c10",
        "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10",
        "baysor_10x_ranger_prior_c10",
        "voronoi_10x_ranger", "10x_ranger",
        "cellpose_exp10um", "cellpose_exp20um",
        "stardist_exp10um", "stardist_exp20um",
        "mesmer_exp10um", "mesmer_exp20um",
        "10x_ranger_exp10um", "10x_ranger_exp20um",
        "watershed_10x", "watershed_stardist", "watershed_mesmer",
        "cellpose_cyto3", "cellpose_cyto3_eosin", "cellpose_cyto3_density",
        "mesmer_wholecell_eosin", "mesmer_wholecell_density",
        "bidcell", "segger",
    ]
    for method_key in optional_methods:
        if method_key not in adatas:
            print(f"\n=== {METHOD_LABELS.get(method_key, method_key)}: skipped (not found) ===")
            continue

        label = METHOD_LABELS.get(method_key, method_key)
        adata_m = adatas[method_key]
        print(f"\n=== 10x native vs. {label} ===")

        matches = match_cells_by_centroid(adata_10x, adata_m, max_dist=MAX_MATCH_DIST)
        matches.to_csv(TABLES_DIR / f"matches_10x_{method_key}.csv", index=False)
        print(f"  {len(matches)} matched pairs")

        corr = expression_correlation(adata_10x, adata_m, matches)
        corr.to_csv(TABLES_DIR / f"expression_correlation_10x_{method_key}.csv", index=False)
        print(f"  median correlation: {corr['correlation'].median():.3f}")

        labels_m = cluster_cell_types(adata_m, seed=0)
        print(f"  {labels_m.nunique()} clusters")
        cluster_embedding(adata_m, seed=0).to_csv(
            TABLES_DIR / f"embedding_{method_key}.csv"
        )

        agreement = cell_type_agreement(labels_10x, labels_m, matches)
        agreement["confusion"].to_csv(
            TABLES_DIR / f"cell_type_confusion_10x_{method_key}.csv"
        )
        print(f"  ARI = {agreement['ari']:.4f}, n_matched = {agreement['n_matched']}")

        _save_disagreement(method_key, matches, labels_10x, labels_m, adata_10x)


if __name__ == "__main__":
    main()
