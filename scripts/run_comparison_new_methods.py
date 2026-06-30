"""B1 (10x native anchor) comparison for the 9 methods added after the initial run.

Produces the same output files as run_comparison.py:
  - matches_10x_{key}.csv
  - expression_correlation_10x_{key}.csv
  - cell_type_confusion_10x_{key}.csv
  - disagreement_table_10x_{key}.csv / _argmax.csv
  - disagreement_spatial_10x_{key}.json / _argmax.json
  - embedding_{key}.csv  (skipped if already exists)

Also writes results/tables/b1_ari_summary.csv consolidating ARI for all
methods that have cell_type_confusion_10x_*.csv files.

Usage::

    conda run -n segbench python scripts/run_comparison_new_methods.py
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from segbench.compare import (
    cell_type_agreement,
    cluster_cell_types,
    cluster_embedding,
    expression_correlation,
    match_cells_by_centroid,
    match_cluster_labels,
    match_cluster_labels_argmax,
)
from segbench.constants import METHOD_LABELS
from segbench.spatial import disagreement_spatial_structure, disagreement_table

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

MAX_MATCH_DIST = 10.0

# Methods missing their B1 comparison outputs — these are new adatas built
# since the initial run_comparison.py was last executed.
MISSING_METHODS: list[tuple[str, str]] = [
    ("baysor_prior_c05",     "adata_baysor_prior_c05.h5ad"),
    ("stardist_exp10um",     "adata_stardist_exp10um.h5ad"),
    ("stardist_exp20um",     "adata_stardist_exp20um.h5ad"),
    ("mesmer_exp10um",       "adata_mesmer_exp10um.h5ad"),
    ("mesmer_exp20um",       "adata_mesmer_exp20um.h5ad"),
    ("10x_ranger_exp10um",   "adata_10x_ranger_exp10um.h5ad"),
    ("10x_ranger_exp20um",   "adata_10x_ranger_exp20um.h5ad"),
    ("watershed_stardist",   "adata_watershed_stardist.h5ad"),
    ("watershed_mesmer",     "adata_watershed_mesmer.h5ad"),
]


def _save_disagreement(key, matches, labels_10x, labels_comp, adata_10x):
    for suffix, matcher_fn in [
        ("", match_cluster_labels),
        ("_argmax", match_cluster_labels_argmax),
    ]:
        labels_aligned = matcher_fn(labels_10x, labels_comp, matches)
        dt = disagreement_table(matches, labels_10x, labels_aligned, adata_10x)
        dt.to_csv(TABLES_DIR / f"disagreement_table_10x_{key}{suffix}.csv", index=False)
        spatial = disagreement_spatial_structure(dt, n_perm=9999, seed=0)
        with open(TABLES_DIR / f"disagreement_spatial_10x_{key}{suffix}.json", "w") as f:
            json.dump(spatial, f, indent=2)
        tag = "hungarian" if not suffix else "argmax"
        print(f"  [{tag}] disagree={dt['disagree'].mean():.3f}  Moran's I={spatial.get('morans_i', 'N/A')}")


def _build_b1_summary() -> pd.DataFrame:
    """Read all existing cell_type_confusion_10x_*.csv files and extract ARI."""
    rows = []
    for conf_file in sorted(TABLES_DIR.glob("cell_type_confusion_10x_*.csv")):
        key = conf_file.stem.replace("cell_type_confusion_10x_", "")
        label = METHOD_LABELS.get(key, key)
        conf = pd.read_csv(conf_file, index_col=0)
        n_matched = int(conf.values.sum())
        rows.append({"method_key": key, "method": label, "n_matched": n_matched})
    return pd.DataFrame(rows)


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading 10x native anchor and computing Leiden clusters...")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    labels_10x = cluster_cell_types(adata_10x, seed=0)
    print(f"  10x native: {adata_10x.n_obs} cells, {labels_10x.nunique()} clusters")

    summary_rows = []

    for key, fname in MISSING_METHODS:
        path = ROI_DIR / fname
        if not path.exists():
            print(f"\n  {key}: skipped (file not found)")
            continue

        label = METHOD_LABELS.get(key, key)
        print(f"\n=== 10x native vs. {label} ===")
        adata_m = ad.read_h5ad(path)

        matches = match_cells_by_centroid(adata_10x, adata_m, max_dist=MAX_MATCH_DIST)
        matches.to_csv(TABLES_DIR / f"matches_10x_{key}.csv", index=False)
        print(f"  {len(matches)} matched pairs")

        corr = expression_correlation(adata_10x, adata_m, matches)
        corr.to_csv(TABLES_DIR / f"expression_correlation_10x_{key}.csv", index=False)
        print(f"  median correlation: {corr['correlation'].median():.3f}")

        labels_m = cluster_cell_types(adata_m, seed=0)
        print(f"  {labels_m.nunique()} clusters")

        emb_path = TABLES_DIR / f"embedding_{key}.csv"
        if not emb_path.exists():
            cluster_embedding(adata_m, seed=0).to_csv(emb_path)

        agreement = cell_type_agreement(labels_10x, labels_m, matches)
        agreement["confusion"].to_csv(
            TABLES_DIR / f"cell_type_confusion_10x_{key}.csv"
        )
        ari = agreement["ari"]
        print(f"  ARI = {ari:.4f}, n_matched = {agreement['n_matched']}")

        _save_disagreement(key, matches, labels_10x, labels_m, adata_10x)

        summary_rows.append({
            "method_key": key,
            "method": label,
            "n_cells": adata_m.n_obs,
            "n_matched": agreement["n_matched"],
            "median_corr": corr["correlation"].median(),
            "b1_ari": ari,
        })

    if summary_rows:
        new_df = pd.DataFrame(summary_rows)
        new_df.to_csv(TABLES_DIR / "b1_new_methods_summary.csv", index=False)
        print("\n=== New methods B1 summary ===")
        print(new_df[["method", "n_cells", "n_matched", "median_corr", "b1_ari"]].to_string(index=False))


if __name__ == "__main__":
    main()
