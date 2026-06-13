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
)
from segbench.spatial import disagreement_spatial_structure, disagreement_table

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

MAX_MATCH_DIST = 10.0  # microns


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adatas = {"cellpose": adata_cellpose, "baysor": adata_baysor}

    counts = cell_count_summary(adatas)
    counts.to_csv(TABLES_DIR / "cell_counts.csv")
    print("=== Cell counts ===")
    print(counts)

    sizes = size_summary(adatas)
    sizes.to_csv(TABLES_DIR / "size_summary.csv")
    print("\n=== Size summary ===")
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


if __name__ == "__main__":
    main()
