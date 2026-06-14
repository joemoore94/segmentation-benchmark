"""Mellon cell-state density vs. segmentation-method disagreement.

For each CellPose-anchored pairwise comparison (Baysor, 10x native, StarDist),
tests whether cells where the two methods' cell-type calls disagree sit in
lower-density regions of CellPose's phenotypic (PCA) space than cells where
they agree -- i.e. whether segmentation disagreement concentrates on
phenotypically ambiguous/transitional cells (Mellon, Otto et al. 2024,
Nature Methods).

Reads ``data/processed/roi/adata_cellpose.h5ad`` and the
``disagreement_table*.csv`` files produced by ``run_comparison.py``, and
writes ``results/tables/cellpose_log_density.csv`` and
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

# Each table's `id_a` column indexes into adata_cellpose.obs_names.
COMPARISONS = {
    "CellPose vs. Baysor": "disagreement_table.csv",
    "CellPose vs. 10x native": "disagreement_table_cellpose_10x.csv",
    "CellPose vs. StarDist": "disagreement_table_cellpose_stardist.csv",
    "CellPose vs. Baysor (prior)": "disagreement_table_cellpose_baysor_prior.csv",
}


def compute_log_density(adata: ad.AnnData, seed: int = 0) -> pd.Series:
    """Normalize -> log1p -> PCA -> Mellon density estimate, per cell."""
    a = adata.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=30, random_state=seed)
    est = mellon.DensityEstimator(random_state=seed)
    log_density = est.fit_predict(a.obsm["X_pca"])
    return pd.Series(log_density, index=a.obs_names.astype(int), name="log_density")


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    print(f"Estimating CellPose cell-state density with Mellon ({adata_cellpose.n_obs} cells)...")
    log_density = compute_log_density(adata_cellpose)
    log_density.to_csv(TABLES_DIR / "cellpose_log_density.csv")

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
    print(f"\nwrote {TABLES_DIR / 'density_disagreement_summary.csv'}")


if __name__ == "__main__":
    main()
