"""Cross-method comparison of CellPose vs. Baysor segmentation on the ROI.

Loads the per-method AnnData produced by ``quantify_masks.py`` (CellPose) and
the Baysor quantification, then runs the comparison and spatial-structure
metrics from ``segbench.compare`` / ``segbench.spatial`` and writes summary
tables to ``results/tables/``.

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
    expression_correlation,
    match_cells_by_centroid,
    match_cluster_labels,
    size_summary,
    subset_to_region,
)
from segbench.spatial import disagreement_spatial_structure, disagreement_table

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

MAX_MATCH_DIST = 10.0  # microns

# Baysor only segmented the centered 1mm x 1mm sub-region of the 2mm x 2mm
# ROI (CPU-tractability, see docs/dataset.md). Subsetting CellPose to the
# same sub-region gives a direct, area-matched cell count/size comparison.
SUB_REGION = ((500.0, 1500.0), (500.0, 1500.0))  # (x_range, y_range), microns

# qv>=20 non-control transcripts in the 1mm x 1mm sub-region (input to
# Baysor; see docs/dataset.md) -- denominator for transcript capture rate.
TOTAL_TRANSCRIPTS_1MM2 = 770_748


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adatas = {"cellpose": adata_cellpose, "baysor": adata_baysor, "10x_native": adata_10x}

    counts = cell_count_summary(adatas)
    counts.to_csv(TABLES_DIR / "cell_counts.csv")
    print("=== Cell counts (full ROI per method) ===")
    print(counts)

    sizes = size_summary(adatas)
    sizes.to_csv(TABLES_DIR / "size_summary.csv")
    print("\n=== Size summary (full ROI per method) ===")
    print(sizes)

    # Direct, area-matched comparison: CellPose and 10x_native subset to
    # Baysor's 1mm x 1mm sub-region, so raw cell counts/sizes are comparable
    # without density normalization.
    x_range, y_range = SUB_REGION
    adata_cellpose_sub = subset_to_region(adata_cellpose, x_range, y_range)
    adata_10x_sub = subset_to_region(adata_10x, x_range, y_range)
    adatas_sub = {
        "cellpose": adata_cellpose_sub,
        "baysor": adata_baysor,
        "10x_native": adata_10x_sub,
    }

    counts_sub = cell_count_summary(adatas_sub)
    counts_sub["transcript_capture_rate"] = (
        counts_sub["total_transcripts"] / TOTAL_TRANSCRIPTS_1MM2
    )
    counts_sub.to_csv(TABLES_DIR / "cell_counts_1mm2.csv")
    print("\n=== Cell counts + QC (1mm x 1mm sub-region, all methods) ===")
    print(counts_sub)

    sizes_sub = size_summary(adatas_sub)
    sizes_sub.to_csv(TABLES_DIR / "size_summary_1mm2.csv")
    print("\n=== Size summary (1mm x 1mm sub-region, all methods) ===")
    print(sizes_sub)

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

    print("\n=== Clustering cell types (Baysor) ===")
    labels_baysor = cluster_cell_types(adata_baysor, seed=0)
    print(f"baysor: {labels_baysor.nunique()} clusters")

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

    spatial = disagreement_spatial_structure(disagreement, seed=0)
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

    spatial_10x = disagreement_spatial_structure(disagreement_10x, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_cellpose_10x.json", "w") as f:
        json.dump(spatial_10x, f, indent=2)
    print("\n=== CellPose vs 10x_native spatial structure of disagreement ===")
    print(spatial_10x)


if __name__ == "__main__":
    main()
