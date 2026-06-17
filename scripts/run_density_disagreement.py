"""Mellon cell-state density vs. segmentation-method disagreement.

For each 10x-native-anchored pairwise comparison, tests whether cells where the
two methods' cell-type calls disagree sit in lower-density regions of 10x
native's phenotypic (PCA) space than cells where they agree -- i.e. whether
segmentation disagreement concentrates on phenotypically ambiguous/transitional
cells (Mellon, Otto et al. 2024, Nature Methods).

Reads ``data/processed/roi/adata_10x.h5ad`` and the
``disagreement_table_10x_*.csv`` files produced by ``run_comparison.py``, and
writes ``results/tables/10x_log_density.csv`` and
``results/tables/density_disagreement_summary.csv``.

Usage::

    conda run -n segbench python scripts/run_density_disagreement.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mellon
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

COMPARISONS = {
    "10x native vs. CellPose": "disagreement_table_10x_cellpose.csv",
    "10x native vs. StarDist": "disagreement_table_10x_stardist.csv",
    "10x native vs. Voronoi": "disagreement_table_10x_voronoi.csv",
    "10x native vs. Baysor": "disagreement_table_10x_baysor.csv",
    "10x native vs. Baysor (prior)": "disagreement_table_10x_baysor_prior.csv",
}


def compute_log_density(adata: ad.AnnData, seed: int = 0) -> pd.Series:
    a = adata.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=30, random_state=seed)
    est = mellon.DensityEstimator(random_state=seed)
    log_density = est.fit_predict(a.obsm["X_pca"])
    return pd.Series(log_density, index=a.obs_names, name="log_density")


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    print(f"Estimating 10x native cell-state density with Mellon ({adata_10x.n_obs} cells)...")
    log_density = compute_log_density(adata_10x)
    log_density.to_csv(TABLES_DIR / "10x_log_density.csv")

    rows = []
    for label, fname in COMPARISONS.items():
        disagreement = pd.read_csv(TABLES_DIR / fname)
        disagreement["log_density"] = disagreement["id_a"].map(log_density)

        agree = disagreement.loc[disagreement["disagree"] == 0.0, "log_density"].dropna()
        disagree = disagreement.loc[disagreement["disagree"] == 1.0, "log_density"].dropna()

        stat, p = mannwhitneyu(disagree, agree, alternative="two-sided")
        row = {
            "comparison": label,
            "n_agree": len(agree),
            "n_disagree": len(disagree),
            "median_log_density_agree": agree.median(),
            "median_log_density_disagree": disagree.median(),
            "mannwhitney_u": stat,
            "p_value": p,
        }
        rows.append(row)
        print(f"\n=== {label}: density vs. disagreement ===")
        print(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES_DIR / "density_disagreement_summary.csv", index=False)


if __name__ == "__main__":
    main()
