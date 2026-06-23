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
    size_summary,
)
from segbench.spatial import disagreement_spatial_structure, disagreement_table

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

MAX_MATCH_DIST = 10.0
TOTAL_TRANSCRIPTS_FULL_ROI = 3_392_051


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

    segger_path = ROI_DIR / "adata_segger.h5ad"
    if segger_path.exists():
        adata_segger = ad.read_h5ad(segger_path)
        adatas["segger"] = adata_segger

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

    labels_cellpose_aligned = match_cluster_labels(labels_10x, labels_cellpose, matches_cellpose)
    disagreement_cellpose = disagreement_table(
        matches_cellpose, labels_10x, labels_cellpose_aligned, adata_10x
    )
    disagreement_cellpose.to_csv(TABLES_DIR / "disagreement_table_10x_cellpose.csv", index=False)
    spatial_cellpose = disagreement_spatial_structure(disagreement_cellpose, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_10x_cellpose.json", "w") as f:
        json.dump(spatial_cellpose, f, indent=2)
    print(spatial_cellpose)

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

    labels_baysor_aligned = match_cluster_labels(labels_10x, labels_baysor, matches_baysor)
    disagreement_baysor = disagreement_table(
        matches_baysor, labels_10x, labels_baysor_aligned, adata_10x
    )
    disagreement_baysor.to_csv(TABLES_DIR / "disagreement_table_10x_baysor.csv", index=False)
    spatial_baysor = disagreement_spatial_structure(disagreement_baysor, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_10x_baysor.json", "w") as f:
        json.dump(spatial_baysor, f, indent=2)
    print(spatial_baysor)

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

    labels_stardist_aligned = match_cluster_labels(labels_10x, labels_stardist, matches_stardist)
    disagreement_stardist = disagreement_table(
        matches_stardist, labels_10x, labels_stardist_aligned, adata_10x
    )
    disagreement_stardist.to_csv(TABLES_DIR / "disagreement_table_10x_stardist.csv", index=False)
    spatial_stardist = disagreement_spatial_structure(disagreement_stardist, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_10x_stardist.json", "w") as f:
        json.dump(spatial_stardist, f, indent=2)
    print(spatial_stardist)

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

    labels_baysor_prior_aligned = match_cluster_labels(
        labels_10x, labels_baysor_prior, matches_baysor_prior
    )
    disagreement_baysor_prior = disagreement_table(
        matches_baysor_prior, labels_10x, labels_baysor_prior_aligned, adata_10x
    )
    disagreement_baysor_prior.to_csv(
        TABLES_DIR / "disagreement_table_10x_baysor_prior.csv", index=False
    )
    spatial_baysor_prior = disagreement_spatial_structure(
        disagreement_baysor_prior, n_perm=9999, seed=0
    )
    with open(TABLES_DIR / "disagreement_spatial_10x_baysor_prior.json", "w") as f:
        json.dump(spatial_baysor_prior, f, indent=2)
    print(spatial_baysor_prior)

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

    labels_voronoi_aligned = match_cluster_labels(labels_10x, labels_voronoi, matches_voronoi)
    disagreement_voronoi = disagreement_table(
        matches_voronoi, labels_10x, labels_voronoi_aligned, adata_10x
    )
    disagreement_voronoi.to_csv(TABLES_DIR / "disagreement_table_10x_voronoi.csv", index=False)
    spatial_voronoi = disagreement_spatial_structure(disagreement_voronoi, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_10x_voronoi.json", "w") as f:
        json.dump(spatial_voronoi, f, indent=2)
    print(spatial_voronoi)

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

    labels_mesmer_aligned = match_cluster_labels(labels_10x, labels_mesmer, matches_mesmer)
    disagreement_mesmer = disagreement_table(
        matches_mesmer, labels_10x, labels_mesmer_aligned, adata_10x
    )
    disagreement_mesmer.to_csv(TABLES_DIR / "disagreement_table_10x_mesmer.csv", index=False)
    spatial_mesmer = disagreement_spatial_structure(disagreement_mesmer, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_10x_mesmer.json", "w") as f:
        json.dump(spatial_mesmer, f, indent=2)
    print(spatial_mesmer)

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

    labels_voronoi_mesmer_aligned = match_cluster_labels(
        labels_10x, labels_voronoi_mesmer, matches_voronoi_mesmer
    )
    disagreement_voronoi_mesmer = disagreement_table(
        matches_voronoi_mesmer, labels_10x, labels_voronoi_mesmer_aligned, adata_10x
    )
    disagreement_voronoi_mesmer.to_csv(
        TABLES_DIR / "disagreement_table_10x_voronoi_mesmer.csv", index=False
    )
    spatial_voronoi_mesmer = disagreement_spatial_structure(
        disagreement_voronoi_mesmer, n_perm=9999, seed=0
    )
    with open(TABLES_DIR / "disagreement_spatial_10x_voronoi_mesmer.json", "w") as f:
        json.dump(spatial_voronoi_mesmer, f, indent=2)
    print(spatial_voronoi_mesmer)

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

    labels_voronoi_stardist_aligned = match_cluster_labels(
        labels_10x, labels_voronoi_stardist, matches_voronoi_stardist
    )
    disagreement_voronoi_stardist = disagreement_table(
        matches_voronoi_stardist, labels_10x, labels_voronoi_stardist_aligned, adata_10x
    )
    disagreement_voronoi_stardist.to_csv(
        TABLES_DIR / "disagreement_table_10x_voronoi_stardist.csv", index=False
    )
    spatial_voronoi_stardist = disagreement_spatial_structure(
        disagreement_voronoi_stardist, n_perm=9999, seed=0
    )
    with open(TABLES_DIR / "disagreement_spatial_10x_voronoi_stardist.json", "w") as f:
        json.dump(spatial_voronoi_stardist, f, indent=2)
    print(spatial_voronoi_stardist)

    # 10x native vs. Baysor (prior c=0.8), if available
    baysor_c08_path = ROI_DIR / "adata_baysor_prior_c08.h5ad"
    if baysor_c08_path.exists():
        adata_baysor_c08 = ad.read_h5ad(baysor_c08_path)
        adatas["baysor_prior_c08"] = adata_baysor_c08
        matches_baysor_c08 = match_cells_by_centroid(
            adata_10x, adata_baysor_c08, max_dist=MAX_MATCH_DIST
        )
        matches_baysor_c08.to_csv(TABLES_DIR / "matches_10x_baysor_prior_c08.csv", index=False)
        print(f"\n=== 10x native vs. Baysor (prior c=0.8) ({len(matches_baysor_c08)} pairs) ===")

        corr_baysor_c08 = expression_correlation(adata_10x, adata_baysor_c08, matches_baysor_c08)
        corr_baysor_c08.to_csv(
            TABLES_DIR / "expression_correlation_10x_baysor_prior_c08.csv", index=False
        )
        print(corr_baysor_c08["correlation"].describe())

        print("\n=== Clustering cell types (Baysor prior c=0.8) ===")
        labels_baysor_c08 = cluster_cell_types(adata_baysor_c08, seed=0)
        print(f"baysor_prior_c08: {labels_baysor_c08.nunique()} clusters")
        cluster_embedding(adata_baysor_c08, seed=0).to_csv(
            TABLES_DIR / "embedding_baysor_prior_c08.csv"
        )

        agreement_baysor_c08 = cell_type_agreement(
            labels_10x, labels_baysor_c08, matches_baysor_c08
        )
        agreement_baysor_c08["confusion"].to_csv(
            TABLES_DIR / "cell_type_confusion_10x_baysor_prior_c08.csv"
        )
        print(f"ARI = {agreement_baysor_c08['ari']:.4f}, n_matched = {agreement_baysor_c08['n_matched']}")

        labels_baysor_c08_aligned = match_cluster_labels(
            labels_10x, labels_baysor_c08, matches_baysor_c08
        )
        disagreement_baysor_c08 = disagreement_table(
            matches_baysor_c08, labels_10x, labels_baysor_c08_aligned, adata_10x
        )
        disagreement_baysor_c08.to_csv(
            TABLES_DIR / "disagreement_table_10x_baysor_prior_c08.csv", index=False
        )
        spatial_baysor_c08 = disagreement_spatial_structure(
            disagreement_baysor_c08, n_perm=9999, seed=0
        )
        with open(TABLES_DIR / "disagreement_spatial_10x_baysor_prior_c08.json", "w") as f:
            json.dump(spatial_baysor_c08, f, indent=2)
        print(spatial_baysor_c08)
    else:
        print("\n=== Baysor (prior c=0.8): skipped (file not found) ===")

    # 10x native vs. Baysor (prior c=1.0), if available
    baysor_c10_path = ROI_DIR / "adata_baysor_prior_c10.h5ad"
    if baysor_c10_path.exists():
        adata_baysor_c10 = ad.read_h5ad(baysor_c10_path)
        adatas["baysor_prior_c10"] = adata_baysor_c10
        matches_baysor_c10 = match_cells_by_centroid(
            adata_10x, adata_baysor_c10, max_dist=MAX_MATCH_DIST
        )
        matches_baysor_c10.to_csv(TABLES_DIR / "matches_10x_baysor_prior_c10.csv", index=False)
        print(f"\n=== 10x native vs. Baysor (prior c=1.0) ({len(matches_baysor_c10)} pairs) ===")

        corr_baysor_c10 = expression_correlation(adata_10x, adata_baysor_c10, matches_baysor_c10)
        corr_baysor_c10.to_csv(
            TABLES_DIR / "expression_correlation_10x_baysor_prior_c10.csv", index=False
        )
        print(corr_baysor_c10["correlation"].describe())

        print("\n=== Clustering cell types (Baysor prior c=1.0) ===")
        labels_baysor_c10 = cluster_cell_types(adata_baysor_c10, seed=0)
        print(f"baysor_prior_c10: {labels_baysor_c10.nunique()} clusters")
        cluster_embedding(adata_baysor_c10, seed=0).to_csv(
            TABLES_DIR / "embedding_baysor_prior_c10.csv"
        )

        agreement_baysor_c10 = cell_type_agreement(
            labels_10x, labels_baysor_c10, matches_baysor_c10
        )
        agreement_baysor_c10["confusion"].to_csv(
            TABLES_DIR / "cell_type_confusion_10x_baysor_prior_c10.csv"
        )
        print(f"ARI = {agreement_baysor_c10['ari']:.4f}, n_matched = {agreement_baysor_c10['n_matched']}")

        labels_baysor_c10_aligned = match_cluster_labels(
            labels_10x, labels_baysor_c10, matches_baysor_c10
        )
        disagreement_baysor_c10 = disagreement_table(
            matches_baysor_c10, labels_10x, labels_baysor_c10_aligned, adata_10x
        )
        disagreement_baysor_c10.to_csv(
            TABLES_DIR / "disagreement_table_10x_baysor_prior_c10.csv", index=False
        )
        spatial_baysor_c10 = disagreement_spatial_structure(
            disagreement_baysor_c10, n_perm=9999, seed=0
        )
        with open(TABLES_DIR / "disagreement_spatial_10x_baysor_prior_c10.json", "w") as f:
            json.dump(spatial_baysor_c10, f, indent=2)
        print(spatial_baysor_c10)
    else:
        print("\n=== Baysor (prior c=1.0): skipped (file not found) ===")

    # 10x native vs. BIDCell (multimodal, if available)
    bidcell_path = ROI_DIR / "adata_bidcell.h5ad"
    if bidcell_path.exists():
        adata_bidcell = ad.read_h5ad(bidcell_path)
        adatas["bidcell"] = adata_bidcell
        matches_bidcell = match_cells_by_centroid(
            adata_10x, adata_bidcell, max_dist=MAX_MATCH_DIST
        )
        matches_bidcell.to_csv(TABLES_DIR / "matches_10x_bidcell.csv", index=False)
        print(f"\n=== 10x native vs. BIDCell ({len(matches_bidcell)} pairs) ===")

        corr_bidcell = expression_correlation(adata_10x, adata_bidcell, matches_bidcell)
        corr_bidcell.to_csv(
            TABLES_DIR / "expression_correlation_10x_bidcell.csv", index=False
        )
        print(corr_bidcell["correlation"].describe())

        print("\n=== Clustering cell types (BIDCell) ===")
        labels_bidcell = cluster_cell_types(adata_bidcell, seed=0)
        print(f"bidcell: {labels_bidcell.nunique()} clusters")
        cluster_embedding(adata_bidcell, seed=0).to_csv(
            TABLES_DIR / "embedding_bidcell.csv"
        )

        agreement_bidcell = cell_type_agreement(
            labels_10x, labels_bidcell, matches_bidcell
        )
        agreement_bidcell["confusion"].to_csv(
            TABLES_DIR / "cell_type_confusion_10x_bidcell.csv"
        )
        print(f"ARI = {agreement_bidcell['ari']:.4f}, n_matched = {agreement_bidcell['n_matched']}")

        labels_bidcell_aligned = match_cluster_labels(
            labels_10x, labels_bidcell, matches_bidcell
        )
        disagreement_bidcell = disagreement_table(
            matches_bidcell, labels_10x, labels_bidcell_aligned, adata_10x
        )
        disagreement_bidcell.to_csv(
            TABLES_DIR / "disagreement_table_10x_bidcell.csv", index=False
        )
        spatial_bidcell = disagreement_spatial_structure(
            disagreement_bidcell, n_perm=9999, seed=0
        )
        with open(TABLES_DIR / "disagreement_spatial_10x_bidcell.json", "w") as f:
            json.dump(spatial_bidcell, f, indent=2)
        print(spatial_bidcell)
    else:
        print("\n=== BIDCell: skipped (adata_bidcell.h5ad not found) ===")

    # 10x native vs. Segger (GNN multimodal, if available)
    if "segger" in adatas:
        adata_segger = adatas["segger"]
        matches_segger = match_cells_by_centroid(
            adata_10x, adata_segger, max_dist=MAX_MATCH_DIST
        )
        matches_segger.to_csv(TABLES_DIR / "matches_10x_segger.csv", index=False)
        print(f"\n=== 10x native vs. Segger ({len(matches_segger)} pairs) ===")

        corr_segger = expression_correlation(adata_10x, adata_segger, matches_segger)
        corr_segger.to_csv(TABLES_DIR / "expression_correlation_10x_segger.csv", index=False)
        print(corr_segger["correlation"].describe())

        print("\n=== Clustering cell types (Segger) ===")
        labels_segger = cluster_cell_types(adata_segger, seed=0)
        print(f"segger: {labels_segger.nunique()} clusters")
        cluster_embedding(adata_segger, seed=0).to_csv(TABLES_DIR / "embedding_segger.csv")

        agreement_segger = cell_type_agreement(labels_10x, labels_segger, matches_segger)
        agreement_segger["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_segger.csv")
        print(f"ARI = {agreement_segger['ari']:.4f}, n_matched = {agreement_segger['n_matched']}")

        labels_segger_aligned = match_cluster_labels(labels_10x, labels_segger, matches_segger)
        disagreement_segger = disagreement_table(
            matches_segger, labels_10x, labels_segger_aligned, adata_10x
        )
        disagreement_segger.to_csv(TABLES_DIR / "disagreement_table_10x_segger.csv", index=False)
        spatial_segger = disagreement_spatial_structure(disagreement_segger, n_perm=9999, seed=0)
        with open(TABLES_DIR / "disagreement_spatial_10x_segger.json", "w") as f:
            json.dump(spatial_segger, f, indent=2)
        print(spatial_segger)
    else:
        print("\n=== Segger: skipped (adata_segger.h5ad not found) ===")


if __name__ == "__main__":
    main()
