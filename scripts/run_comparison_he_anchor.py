"""Cross-method comparison anchored to the H&E morphology segmentation.

Uses Mesmer WC DAPI+eosin as the H&E reference anchor. Eosin stains cytoplasm
and stroma from the registered H&E image, making this the closest available
proxy for true whole-cell morphology. Comparing all methods against this anchor
reveals which methods approximate real cell boundaries vs. which reproduce the
DAPI-only expansion logic of 10x native.

Key difference from run_comparison.py: no spatial disagreement tables (those
are 10x-native-specific); this script focuses on the three benchmark metrics
that generalize across anchors — centroid matching, expression correlation,
and Leiden clustering ARI.

Outputs go to results/tables/ with _he_anchor suffix to avoid collisions with
the 10x-native-anchored equivalents.

Usage::

    conda run -n segbench python scripts/run_comparison_he_anchor.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from segbench.compare import (
    cell_count_summary,
    cell_type_agreement,
    cluster_cell_types,
    cluster_embedding,
    expression_correlation,
    match_cells_by_centroid,
)
from segbench.constants import METHOD_LABELS, TOTAL_TRANSCRIPTS_FULL_ROI

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
MAX_MATCH_DIST = 10.0


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    anchor_key = "mesmer_wholecell_eosin"
    anchor_label = METHOD_LABELS.get(anchor_key, anchor_key)
    adata_anchor = ad.read_h5ad(ROI_DIR / "adata_mesmer_wholecell_eosin.h5ad")
    print(f"H&E anchor: {anchor_label}  ({adata_anchor.n_obs} cells)")

    labels_anchor = cluster_cell_types(adata_anchor, seed=0)
    print(f"  {labels_anchor.nunique()} Leiden clusters")
    cluster_embedding(adata_anchor, seed=0).to_csv(
        TABLES_DIR / "embedding_he_anchor.csv"
    )

    # All methods to compare against the H&E anchor. This includes 10x native
    # (which is otherwise used as anchor in run_comparison.py) so it gets a
    # score against the H&E reference like every other method.
    all_methods: list[tuple[str, str]] = [
        ("10x_native",                  "adata_10x.h5ad"),
        ("cellpose",                    "adata_cellpose.h5ad"),
        ("stardist",                    "adata_stardist.h5ad"),
        ("mesmer",                      "adata_mesmer.h5ad"),
        ("10x_ranger",                  "adata_10x_ranger.h5ad"),
        ("voronoi",                     "adata_voronoi.h5ad"),
        ("voronoi_stardist",            "adata_voronoi_stardist.h5ad"),
        ("voronoi_mesmer",              "adata_voronoi_mesmer.h5ad"),
        ("voronoi_10x_ranger",          "adata_voronoi_10x_ranger.h5ad"),
        ("baysor",                      "adata_baysor.h5ad"),
        ("baysor_prior",                "adata_baysor_prior.h5ad"),
        ("baysor_prior_c05",            "adata_baysor_prior_c05.h5ad"),
        ("baysor_prior_c08",            "adata_baysor_prior_c08.h5ad"),
        ("baysor_prior_c10",            "adata_baysor_prior_c10.h5ad"),
        ("baysor_stardist_prior_c10",   "adata_baysor_stardist_prior_c10.h5ad"),
        ("baysor_mesmer_prior_c10",     "adata_baysor_mesmer_prior_c10.h5ad"),
        ("baysor_10x_ranger_prior_c10", "adata_baysor_10x_ranger_prior_c10.h5ad"),
        ("cellpose_exp10um",            "adata_cellpose_exp10um.h5ad"),
        ("cellpose_exp20um",            "adata_cellpose_exp20um.h5ad"),
        ("stardist_exp10um",            "adata_stardist_exp10um.h5ad"),
        ("stardist_exp20um",            "adata_stardist_exp20um.h5ad"),
        ("mesmer_exp10um",              "adata_mesmer_exp10um.h5ad"),
        ("mesmer_exp20um",              "adata_mesmer_exp20um.h5ad"),
        ("10x_ranger_exp10um",          "adata_10x_ranger_exp10um.h5ad"),
        ("10x_ranger_exp20um",          "adata_10x_ranger_exp20um.h5ad"),
        ("watershed_10x",               "adata_watershed_10x.h5ad"),
        ("watershed_stardist",          "adata_watershed_stardist.h5ad"),
        ("watershed_mesmer",            "adata_watershed_mesmer.h5ad"),
        ("cellpose_cyto3",              "adata_cellpose_cyto3.h5ad"),
        ("cellpose_cyto3_eosin",        "adata_cellpose_cyto3_eosin.h5ad"),
        ("cellpose_cyto3_density",      "adata_cellpose_cyto3_density.h5ad"),
        ("mesmer_wholecell_density",    "adata_mesmer_wholecell_density.h5ad"),
    ]

    rows = []
    for method_key, filename in all_methods:
        path = ROI_DIR / filename
        if not path.exists():
            print(f"  [skip] {method_key}: {filename} not found")
            continue

        label = METHOD_LABELS.get(method_key, method_key)
        print(f"\n=== H&E anchor vs. {label} ===")
        adata_m = ad.read_h5ad(path)

        matches = match_cells_by_centroid(adata_anchor, adata_m, max_dist=MAX_MATCH_DIST)
        matches.to_csv(TABLES_DIR / f"matches_he_anchor_{method_key}.csv", index=False)
        print(f"  {len(matches)} matched pairs")

        corr = expression_correlation(adata_anchor, adata_m, matches)
        corr.to_csv(
            TABLES_DIR / f"expression_correlation_he_anchor_{method_key}.csv", index=False
        )
        median_corr = corr["correlation"].median()
        print(f"  median correlation: {median_corr:.3f}")

        labels_m = cluster_cell_types(adata_m, seed=0)
        cluster_embedding(adata_m, seed=0).to_csv(
            TABLES_DIR / f"embedding_{method_key}.csv"
        )

        agreement = cell_type_agreement(labels_anchor, labels_m, matches)
        agreement["confusion"].to_csv(
            TABLES_DIR / f"cell_type_confusion_he_anchor_{method_key}.csv"
        )
        ari = agreement["ari"]
        n_matched = agreement["n_matched"]
        print(f"  ARI = {ari:.4f}, n_matched = {n_matched}")

        rows.append({
            "method": method_key,
            "label": label,
            "n_cells": adata_m.n_obs,
            "capture_pct": float(adata_m.X.sum()) / TOTAL_TRANSCRIPTS_FULL_ROI * 100,
            "n_matched": n_matched,
            "median_corr": median_corr,
            "ari": ari,
        })

    summary = pd.DataFrame(rows).sort_values("ari", ascending=False)
    summary.to_csv(TABLES_DIR / "comparison_he_anchor_summary.csv", index=False)
    print("\n=== H&E Anchor Summary (sorted by ARI) ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
