"""Quick ARI / disagreement-rate comparison of candidate methods vs. 10x native.

Loads whatever candidate h5ads are present in data/processed/roi/ and runs
the same match + cluster + agreement pipeline as run_comparison.py, printing
a summary table. Does NOT write CSVs or run spatial structure tests — this
is an assessment pass to decide which candidates to promote to the full pipeline.

Expected candidates (any subset that exists will be tested):
  adata_cellpose_exp10um.h5ad    nucleus expansion 10 µm
  adata_cellpose_exp20um.h5ad    nucleus expansion 20 µm
  adata_voronoi.h5ad             Voronoi from CellPose centroids
  adata_baysor_prior_c05.h5ad    Baysor prior confidence 0.5
  adata_baysor_prior_c08.h5ad    Baysor prior confidence 0.8

Usage::

    conda run -n segbench python scripts/run_candidate_comparison.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from segbench.compare import (
    cell_type_agreement,
    cluster_cell_types,
    match_cells_by_centroid,
    match_cluster_labels,
)
from segbench.spatial import disagreement_table, disagreement_spatial_structure

ROI_DIR = Path("data/processed/roi")
MAX_MATCH_DIST = 10.0
TOTAL_TX = 3_392_051

CANDIDATES = {
    "exp10um": "adata_cellpose_exp10um.h5ad",
    "exp20um": "adata_cellpose_exp20um.h5ad",
    "voronoi": "adata_voronoi.h5ad",
    "baysor_prior_c05": "adata_baysor_prior_c05.h5ad",
    "baysor_prior_c08": "adata_baysor_prior_c08.h5ad",
}


def main() -> None:
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    print(f"10x native: {adata_10x.n_obs} cells")

    print("\nClustering 10x native (anchor)...")
    labels_10x = cluster_cell_types(adata_10x, seed=0)

    rows = []
    for name, fname in CANDIDATES.items():
        path = ROI_DIR / fname
        if not path.exists():
            print(f"\n{name}: {fname} not found, skipping")
            continue

        print(f"\n=== {name} ===")
        adata = ad.read_h5ad(path)
        n_cells = adata.n_obs
        n_tx = int(adata.X.sum())
        capture = n_tx / TOTAL_TX

        matches = match_cells_by_centroid(adata_10x, adata, max_dist=MAX_MATCH_DIST)
        labels = cluster_cell_types(adata, seed=0)
        agreement = cell_type_agreement(labels_10x, labels, matches)
        labels_aligned = match_cluster_labels(labels_10x, labels, matches)
        dtable = disagreement_table(matches, labels_10x, labels_aligned, adata_10x)
        spatial = disagreement_spatial_structure(dtable, n_perm=999, seed=0)

        row = {
            "method": name,
            "n_cells": n_cells,
            "capture": f"{capture:.1%}",
            "matched_pairs": len(matches),
            "ari": round(agreement["ari"], 3),
            "disagree_rate": f"{spatial['disagreement_rate']:.1%}",
            "morans_i": round(spatial["morans_i"], 3),
        }
        rows.append(row)
        print(f"  n_cells={n_cells}, capture={capture:.1%}, "
              f"pairs={len(matches)}, ARI={agreement['ari']:.3f}, "
              f"disagree={spatial['disagreement_rate']:.1%}, "
              f"Moran's I={spatial['morans_i']:.3f}")

    print("\n\n=== CANDIDATE SUMMARY (vs. 10x native) ===")
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

    print("\n=== BASELINE (from run_comparison.py, for reference) ===")
    baselines = [
        ("cellpose",       "20,166", "35.4%", "18,966", "0.547", "30.8%", "0.178"),
        ("stardist",       "24,745", "40.8%", "21,429", "0.545", "33.5%", "0.215"),
        ("baysor",         "18,321", "98.6%", "10,953", "0.305", "51.7%", "0.033"),
        ("baysor_prior",   "19,061", "98.7%", "11,454", "0.318", "51.9%", "0.036"),
    ]
    bl_df = pd.DataFrame(
        baselines,
        columns=["method", "n_cells", "capture", "matched_pairs", "ari",
                 "disagree_rate", "morans_i"]
    )
    print(bl_df.to_string(index=False))


if __name__ == "__main__":
    main()
