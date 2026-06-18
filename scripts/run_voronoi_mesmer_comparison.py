"""Run the 10x native vs. Voronoi (Mesmer centroids) comparison only.

Voronoi-Mesmer assigns every transcript to its nearest Mesmer nuclear
centroid, exactly as Voronoi-CellPose does for CellPose centroids.
Comparing the two Voronoi variants isolates the effect of nuclear-detector
choice (CellPose vs. Mesmer) on downstream cell-type calls, independent of
cytoplasmic-transcript coverage (both methods capture 100%).

Usage::

    conda run -n segbench python scripts/run_voronoi_mesmer_comparison.py
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad

from segbench.compare import (
    cell_type_agreement,
    cluster_cell_types,
    cluster_embedding,
    expression_correlation,
    match_cells_by_centroid,
    match_cluster_labels,
)
from segbench.spatial import disagreement_spatial_structure, disagreement_table

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
MAX_MATCH_DIST = 10.0


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_voronoi_mesmer = ad.read_h5ad(ROI_DIR / "adata_voronoi_mesmer.h5ad")

    labels_10x = cluster_cell_types(adata_10x, seed=0)

    matches = match_cells_by_centroid(adata_10x, adata_voronoi_mesmer, max_dist=MAX_MATCH_DIST)
    matches.to_csv(TABLES_DIR / "matches_10x_voronoi_mesmer.csv", index=False)
    print(f"=== 10x native vs. Voronoi (Mesmer) ({len(matches)} pairs) ===")

    corr = expression_correlation(adata_10x, adata_voronoi_mesmer, matches)
    corr.to_csv(TABLES_DIR / "expression_correlation_10x_voronoi_mesmer.csv", index=False)
    print(corr["correlation"].describe())

    print("\n=== Clustering cell types (Voronoi-Mesmer) ===")
    labels_vm = cluster_cell_types(adata_voronoi_mesmer, seed=0)
    print(f"voronoi_mesmer: {labels_vm.nunique()} clusters")
    cluster_embedding(adata_voronoi_mesmer, seed=0).to_csv(
        TABLES_DIR / "embedding_voronoi_mesmer.csv"
    )

    agreement = cell_type_agreement(labels_10x, labels_vm, matches)
    agreement["confusion"].to_csv(TABLES_DIR / "cell_type_confusion_10x_voronoi_mesmer.csv")
    print(f"ARI = {agreement['ari']:.4f}, n_matched = {agreement['n_matched']}")

    labels_vm_aligned = match_cluster_labels(labels_10x, labels_vm, matches)
    disagree = disagreement_table(matches, labels_10x, labels_vm_aligned, adata_10x)
    disagree.to_csv(TABLES_DIR / "disagreement_table_10x_voronoi_mesmer.csv", index=False)

    spatial = disagreement_spatial_structure(disagree, n_perm=9999, seed=0)
    with open(TABLES_DIR / "disagreement_spatial_10x_voronoi_mesmer.json", "w") as f:
        json.dump(spatial, f, indent=2)
    print(spatial)


if __name__ == "__main__":
    main()
