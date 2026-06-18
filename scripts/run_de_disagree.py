"""Differential expression: cells that methods agree on vs. cells they disagree on.

For each pairwise comparison, splits the 10x-native matched cells into two
groups -- those where the two methods assign the same Leiden cluster (agree)
and those where they differ (disagree) -- and finds genes that distinguish
the two groups using a Wilcoxon rank-sum test on log-normalised counts.

The biological question is: what genes are characteristic of cells that are
robustly identifiable regardless of segmentation method vs. cells whose
identity is ambiguous?

Reads ``data/processed/roi/adata_10x.h5ad`` and the per-comparison
``disagreement_table_10x_*.csv`` files. Writes ranked gene tables to
``results/tables/de_disagree_10x_*.csv``.

Usage::

    conda run -n segbench python scripts/run_de_disagree.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")

COMPARISONS = {
    "10x native vs. CellPose": "disagreement_table_10x_cellpose.csv",
    "10x native vs. StarDist": "disagreement_table_10x_stardist.csv",
    "10x native vs. Mesmer": "disagreement_table_10x_mesmer.csv",
    "10x native vs. Voronoi": "disagreement_table_10x_voronoi.csv",
    "10x native vs. Baysor": "disagreement_table_10x_baysor.csv",
    "10x native vs. Baysor (prior)": "disagreement_table_10x_baysor_prior.csv",
}


def run_de(adata_full: ad.AnnData, disagreement: pd.DataFrame, label: str) -> pd.DataFrame:
    cell_ids = disagreement["id_a"].tolist()
    group_map = dict(zip(disagreement["id_a"], disagreement["disagree"].astype(int)))

    adata = adata_full[adata_full.obs_names.isin(cell_ids)].copy()
    adata.obs["group"] = pd.Categorical(
        [("disagree" if group_map[c] == 1 else "agree") for c in adata.obs_names]
    )

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.tl.rank_genes_groups(
        adata,
        groupby="group",
        groups=["disagree"],
        reference="agree",
        method="wilcoxon",
        key_added="de",
        pts=True,
    )

    result = sc.get.rank_genes_groups_df(adata, group="disagree", key="de")
    result.insert(0, "comparison", label)
    return result


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")

    for label, fname in COMPARISONS.items():
        print(f"\n=== DE: {label} ===")
        disagreement = pd.read_csv(TABLES_DIR / fname)
        n_agree = (disagreement["disagree"] == 0.0).sum()
        n_disagree = (disagreement["disagree"] == 1.0).sum()
        print(f"agree: {n_agree}, disagree: {n_disagree}")

        result = run_de(adata_10x, disagreement, label)

        out_name = fname.replace("disagreement_table_", "de_disagree_")
        result.to_csv(TABLES_DIR / out_name, index=False)

        top = result.nsmallest(10, "pvals_adj")[["names", "logfoldchanges", "pvals_adj", "scores"]]
        print("Top 10 upregulated in disagree group:")
        print(top.to_string(index=False))


if __name__ == "__main__":
    main()
