"""Cross-method comparison of CellPose, Baysor, 10x native, and StarDist on the ROI.

Loads the per-method AnnData (CellPose/StarDist from ``quantify_masks.py``,
Baysor from ``build_baysor_adata.py``, 10x native from ``build_10x_adata.py``),
all covering the full 2mm x 2mm ROI, then runs the comparison and
spatial-structure metrics from ``segbench.compare`` / ``segbench.spatial`` and
writes summary tables to ``results/tables/``.

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

MAX_MATCH_DIST = 10.0  # microns

# qv>=20 non-control transcripts in the full 2mm x 2mm ROI (see
# docs/dataset.md) -- denominator for transcript capture rate.
TOTAL_TRANSCRIPTS_FULL_ROI = 3_392_051


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adatas = {
        "cellpose": adata_cellpose,
        "baysor": adata_baysor,
        "10x_native": adata_10x,
        "stardist": adata_stardist,
    }

    counts = cell_count_summary(adatas)
    counts["transcript_capture_rate"] = (
        counts["total_transcripts"] / TOTAL_TRANSCRIPTS_FULL_ROI
    )
    counts.to_csv(TABLES_DIR / "cell_counts.csv")
    print("=== Cell counts + QC (full 2mm x 2mm ROI, all methods) ===")
    print(counts)

    sizes = size_summary(adatas)
    sizes.to_csv(TABLES_DIR / "size_summary.csv")
    print("\n=== Size summary (full ROI per method) ===")
    print(sizes)

    matches = match_cells_by_centroid(adata_cellpose, adata_baysor, max_dist=MAX_MATCH_DIST)
    matches.to_csv(TABLES_DIR / "matches_cellpose_baysor.csv", index=False)
    print(f"\n=== Matched pairs (max_dist={MAX_MATCH_DIST}um) ===")
    print(f"{len(matches)} matched pairs out of {adata_cellpose.n_obs} cellpose / "
          f"{adata_baysor.n_obs} baysor cells")

    corr = expression_correlation(adata_cellpose, adata_baysor, matches)
    corr.to_csv(TABLES_DIR / "expression_correlation.csv", index=False)
    print("\n=== Expression correlation (matched pairs) ===")
    print(corr["correlation"].describe())

    print("\n=== Clustering cell types (CellPose) ===")
    labels_cellpose = cluster_cell_types(adata_cellpose, seed=0)
    print(f"cellpose: {labels_cellpose.nunique()} clusters")
    cluster_embedding(adata_cellpose, seed=0).to_csv(TABLES_DIR / "embedding_cellpose.csv")

    print("\n=== Clustering cell types (Baysor) ===")
    labels_baysor = cluster_cell_types(adata_baysor, seed=0)
    print(f"baysor: {labels_baysor.nunique()} clusters")
    cluster_embedding(adata_baysor, seed=0).to_csv(TABLES_DIR / "embedding_baysor.csv")

    agreement = cell_type_agreement(labels_cellpose, labels_baysor, matches)
    agreement["confusion"].to_csv(TABLES_DIR / "cell_type_confusion.csv")
    print("\n=== Cell type agreement ===")
    print(f"ARI = {agreement['ari']:.4f}, n_matched = {agreement['n_matched']}")
    print(agreement["confusion"])

    # Independent per-method Leiden runs assign arbitrary cluster ids, so
    # remap baysor's clusters onto cellpose's vocabulary before comparing
    # labels directly (ARI above is permutation-invariant and unaffected).
    labels_baysor_aligned = match_cluster_labels(labels_cellpose, labels_baysor, matches)
    disagreement = disagreement_table(
        matches, labels_cellpose, labels_baysor_aligned, adata_cellpose
    )
    disagreement.to_csv(TABLES_DIR / "disagreement_table.csv", index=False)

    spatial = disagreement_spatial_structure(disagreement, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial.json", "w") as f:
        json.dump(spatial, f, indent=2)
    print("\n=== Spatial structure of disagreement ===")
    print(spatial)

    # CellPose vs. 10x's own (Xenium Ranger) segmentation -- both cover the
    # full 2mm x 2mm ROI, so this is the "does our segmentation reproduce the
    # platform's reference segmentation" comparison, run the same way as
    # CellPose vs. Baysor above.
    matches_10x = match_cells_by_centroid(adata_cellpose, adata_10x, max_dist=MAX_MATCH_DIST)
    matches_10x.to_csv(TABLES_DIR / "matches_cellpose_10x.csv", index=False)
    print(f"\n=== CellPose vs 10x_native matched pairs (max_dist={MAX_MATCH_DIST}um) ===")
    print(f"{len(matches_10x)} matched pairs out of {adata_cellpose.n_obs} cellpose / "
          f"{adata_10x.n_obs} 10x_native cells")

    corr_10x = expression_correlation(adata_cellpose, adata_10x, matches_10x)
    corr_10x.to_csv(TABLES_DIR / "expression_correlation_cellpose_10x.csv", index=False)
    print("\n=== CellPose vs 10x_native expression correlation (matched pairs) ===")
    print(corr_10x["correlation"].describe())

    print("\n=== Clustering cell types (10x_native) ===")
    labels_10x = cluster_cell_types(adata_10x, seed=0)
    print(f"10x_native: {labels_10x.nunique()} clusters")
    cluster_embedding(adata_10x, seed=0).to_csv(TABLES_DIR / "embedding_10x_native.csv")

    agreement_10x = cell_type_agreement(labels_cellpose, labels_10x, matches_10x)
    agreement_10x["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_cellpose_10x.csv")
    print("\n=== CellPose vs 10x_native cell type agreement ===")
    print(f"ARI = {agreement_10x['ari']:.4f}, n_matched = {agreement_10x['n_matched']}")
    print(agreement_10x["confusion"])

    labels_10x_aligned = match_cluster_labels(labels_cellpose, labels_10x, matches_10x)
    disagreement_10x = disagreement_table(
        matches_10x, labels_cellpose, labels_10x_aligned, adata_cellpose
    )
    disagreement_10x.to_csv(TABLES_DIR / "disagreement_table_cellpose_10x.csv", index=False)

    spatial_10x = disagreement_spatial_structure(disagreement_10x, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_cellpose_10x.json", "w") as f:
        json.dump(spatial_10x, f, indent=2)
    print("\n=== CellPose vs 10x_native spatial structure of disagreement ===")
    print(spatial_10x)

    # CellPose vs. StarDist -- both are nuclear segmentations of the same
    # DAPI image (full 2mm x 2mm ROI), run the same way as the comparisons
    # above.
    matches_stardist = match_cells_by_centroid(
        adata_cellpose, adata_stardist, max_dist=MAX_MATCH_DIST
    )
    matches_stardist.to_csv(TABLES_DIR / "matches_cellpose_stardist.csv", index=False)
    print(f"\n=== CellPose vs StarDist matched pairs (max_dist={MAX_MATCH_DIST}um) ===")
    print(f"{len(matches_stardist)} matched pairs out of {adata_cellpose.n_obs} cellpose / "
          f"{adata_stardist.n_obs} stardist cells")

    corr_stardist = expression_correlation(adata_cellpose, adata_stardist, matches_stardist)
    corr_stardist.to_csv(TABLES_DIR / "expression_correlation_cellpose_stardist.csv", index=False)
    print("\n=== CellPose vs StarDist expression correlation (matched pairs) ===")
    print(corr_stardist["correlation"].describe())

    print("\n=== Clustering cell types (StarDist) ===")
    labels_stardist = cluster_cell_types(adata_stardist, seed=0)
    print(f"stardist: {labels_stardist.nunique()} clusters")
    cluster_embedding(adata_stardist, seed=0).to_csv(TABLES_DIR / "embedding_stardist.csv")

    agreement_stardist = cell_type_agreement(labels_cellpose, labels_stardist, matches_stardist)
    agreement_stardist["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_cellpose_stardist.csv")
    print("\n=== CellPose vs StarDist cell type agreement ===")
    print(f"ARI = {agreement_stardist['ari']:.4f}, n_matched = {agreement_stardist['n_matched']}")
    print(agreement_stardist["confusion"])

    labels_stardist_aligned = match_cluster_labels(
        labels_cellpose, labels_stardist, matches_stardist
    )
    disagreement_stardist = disagreement_table(
        matches_stardist, labels_cellpose, labels_stardist_aligned, adata_cellpose
    )
    disagreement_stardist.to_csv(TABLES_DIR / "disagreement_table_cellpose_stardist.csv", index=False)

    spatial_stardist = disagreement_spatial_structure(disagreement_stardist, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_cellpose_stardist.json", "w") as f:
        json.dump(spatial_stardist, f, indent=2)
    print("\n=== CellPose vs StarDist spatial structure of disagreement ===")
    print(spatial_stardist)


if __name__ == "__main__":
    main()
