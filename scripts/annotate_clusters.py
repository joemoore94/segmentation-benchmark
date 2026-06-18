"""Annotate 10x-native Leiden clusters with cell type labels.

Reproduces the same normalize -> PCA -> neighbors -> Leiden pipeline used in
run_comparison.py (seed=0, resolution=1.0), computes per-cluster marker genes
(Wilcoxon), applies manual cell type annotations based on canonical breast
tissue markers (Janesick et al. 2023), then:

  1. Saves cluster -> cell type mapping to results/tables/cluster_annotations.csv
  2. Writes annotated UMAP figure (results/figures/umap_annotated.png)
  3. Joins cell type labels onto every disagreement table and saves
     per-cell-type disagreement rates per comparison to
     results/tables/celltype_disagreement.csv
  4. Writes per-comparison per-cell-type disagree breakdown figure
     (results/figures/celltype_disagreement.png)

Usage::

    conda run -n segbench python scripts/annotate_clusters.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

# Manual annotation based on top Wilcoxon marker genes per cluster.
# Clusters 0, 1, 3, 8 are luminal epithelial subtypes; grouped under
# "Luminal epithelial" for cross-method comparisons. Clusters 2 and 7 are
# both myeloid/macrophage populations with slightly different marker profiles.
CLUSTER_ANNOTATIONS = {
    "0":  "Luminal epithelial",   # ESR1, FOXA1, GATA3, ANKRD30A
    "1":  "Luminal epithelial",   # MYBPC1, MUC1, CLIC6, SERPINA3, GATA3, PGR
    "2":  "Macrophages",          # FCER1G, LYZ, FCGR3A, HAVCR2
    "3":  "Luminal epithelial",   # SERPINA3, MUC1, CLIC6, S100A14, FASN, PGR
    "4":  "Myoepithelial",        # ACTA2, KRT14, MYH11, MYLK
    "5":  "T cells",              # CD3E, TRAC, CD96, GZMA, IL7R
    "6":  "B cells",              # MS4A1, BANK1, CD52, SELL
    "7":  "Macrophages",          # CD14, AIF1, LYZ, FCER1G, FGL2
    "8":  "Luminal epithelial",   # TACSTD2, KRT7, GATA3, MUC1, CCND1
    "9":  "CAFs",                 # LUM, SFRP4, FBLN1, THBS2, MMP2
    "10": "Smooth muscle",        # ACTA2, MYH11, MYLK, RGS5, CAV1
    "11": "Endothelial",          # PECAM1, VWF, AQP1, CD93, CLEC14A
    "12": "Plasma cells",         # MZB1, SLAMF7, ITM2C
    "13": "CAFs",                 # POSTN, PDGFRB, CTHRC1, LUM
    "14": "Adipocytes",           # ADIPOQ, PLIN1, PPARG, LPL
}

CELLTYPE_COLORS = {
    "Luminal epithelial": "#E377C2",
    "Myoepithelial":      "#8172B2",
    "T cells":            "#2CA02C",
    "B cells":            "#17BECF",
    "Macrophages":        "#D62728",
    "Plasma cells":       "#FF7F0E",
    "CAFs":               "#7F7F7F",
    "Smooth muscle":      "#BCBD22",
    "Endothelial":        "#1F77B4",
    "Adipocytes":         "#8C564B",
}

COMPARISONS = [
    ("cellpose",       "CellPose"),
    ("stardist",       "StarDist"),
    ("mesmer",         "Mesmer"),
    ("voronoi",        "Voronoi (CP)"),
    ("voronoi_mesmer", "Voronoi (M)"),
    ("baysor",         "Baysor"),
]


def build_annotated_adata() -> ad.AnnData:
    adata = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    sc.settings.verbosity = 0
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph")
    sc.tl.umap(adata, random_state=0)
    adata.obs["cell_type"] = adata.obs["leiden"].map(CLUSTER_ANNOTATIONS).astype("category")
    return adata


def save_annotation_table(adata: ad.AnnData) -> None:
    rows = []
    for cluster, ct in CLUSTER_ANNOTATIONS.items():
        n = (adata.obs["leiden"] == cluster).sum()
        rows.append({"leiden_cluster": cluster, "cell_type": ct, "n_cells": n})
    df = pd.DataFrame(rows).sort_values("leiden_cluster", key=lambda x: x.astype(int))
    df.to_csv(TABLES_DIR / "cluster_annotations.csv", index=False)
    print("Cluster annotations:")
    print(df.to_string(index=False))


def fig_umap_annotated(adata: ad.AnnData) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Left: UMAP by cell type (annotated)
    cell_types = list(CELLTYPE_COLORS.keys())
    for ct in cell_types:
        mask = adata.obs["cell_type"] == ct
        axes[0].scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=CELLTYPE_COLORS[ct], s=4, alpha=0.6, label=ct,
        )
    axes[0].set_xlabel("UMAP1")
    axes[0].set_ylabel("UMAP2")
    axes[0].set_title("10x native: cell type annotation")
    axes[0].legend(markerscale=3, fontsize=9, loc="best")

    # Right: UMAP by Leiden cluster number (to cross-reference)
    n_clusters = adata.obs["leiden"].nunique()
    palette = sns.color_palette("tab20", n_clusters)
    for i, cluster in enumerate(sorted(adata.obs["leiden"].unique(), key=int)):
        mask = adata.obs["leiden"] == cluster
        axes[1].scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=[palette[i]], s=4, alpha=0.6, label=cluster,
        )
    axes[1].set_xlabel("UMAP1")
    axes[1].set_ylabel("UMAP2")
    axes[1].set_title("10x native: Leiden cluster IDs")
    axes[1].legend(markerscale=3, fontsize=8, ncol=2, loc="best", title="Cluster")

    fig.suptitle("10x native — 23,629 cells, 2mm × 2mm ROI")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "umap_annotated.png", dpi=150)
    plt.close(fig)
    print("Saved umap_annotated.png")


def build_celltype_disagreement(adata: ad.AnnData) -> pd.DataFrame:
    # Build cell_id -> cell_type map from 10x native
    id_to_ct = adata.obs["cell_type"].to_dict()

    rows = []
    for method, label in COMPARISONS:
        path = TABLES_DIR / f"disagreement_table_10x_{method}.csv"
        if not path.exists():
            print(f"  missing: {path.name}, skipping")
            continue
        df = pd.read_csv(path)
        df["cell_type"] = df["id_a"].map(id_to_ct)
        for ct in sorted(df["cell_type"].dropna().unique()):
            sub = df[df["cell_type"] == ct]
            n_total = len(sub)
            n_disagree = sub["disagree"].sum()
            rows.append({
                "comparison": label,
                "cell_type": ct,
                "n_matched": n_total,
                "n_disagree": int(n_disagree),
                "disagree_rate": n_disagree / n_total if n_total > 0 else float("nan"),
            })

    result = pd.DataFrame(rows)
    result.to_csv(TABLES_DIR / "celltype_disagreement.csv", index=False)
    return result


def fig_celltype_disagreement(df: pd.DataFrame) -> None:
    comparisons = [label for _, label in COMPARISONS]
    cell_types = list(CELLTYPE_COLORS.keys())

    # Pivot to cell_type × comparison matrix of disagree rates
    pivot = df.pivot(index="cell_type", columns="comparison", values="disagree_rate")
    pivot = pivot.reindex(index=cell_types, columns=comparisons)

    fig, axes = plt.subplots(1, 2, figsize=(18, 6))

    # Left: heatmap of disagree rates
    sns.heatmap(
        pivot * 100,
        annot=True, fmt=".0f", cmap="YlOrRd",
        linewidths=0.5, ax=axes[0],
        cbar_kws={"label": "Disagreement rate (%)"},
        vmin=0, vmax=80,
    )
    axes[0].set_title("Disagreement rate (%) by cell type and method")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    axes[0].tick_params(axis="x", rotation=30)

    # Right: grouped bar chart — each group is a cell type, bars are comparisons
    x = np.arange(len(cell_types))
    width = 0.13
    comparison_colors = ["#4C72B0", "#8172B2", "#D62728", "#17BECF", "#BCBD22", "#DD8452"]
    for i, (comp, color) in enumerate(zip(comparisons, comparison_colors)):
        vals = [pivot.loc[ct, comp] * 100 if ct in pivot.index else 0 for ct in cell_types]
        axes[1].bar(x + i * width, vals, width, label=comp, color=color, alpha=0.85)

    axes[1].set_xticks(x + width * (len(comparisons) - 1) / 2)
    axes[1].set_xticklabels(cell_types, rotation=35, ha="right")
    axes[1].set_ylabel("Disagreement rate (%)")
    axes[1].set_title("Disagreement rate by cell type and method")
    axes[1].legend(fontsize=9)
    axes[1].set_ylim(0, 85)

    fig.suptitle("Which cell types drive method disagreement with 10x native?")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "celltype_disagreement.png", dpi=150)
    plt.close(fig)
    print("Saved celltype_disagreement.png")


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print("Building annotated 10x native AnnData...")
    adata = build_annotated_adata()

    save_annotation_table(adata)
    fig_umap_annotated(adata)

    print("\nBuilding per-cell-type disagreement breakdown...")
    df = build_celltype_disagreement(adata)

    print("\nDisagreement rates by cell type and method:")
    pivot = df.pivot(index="cell_type", columns="comparison", values="disagree_rate")
    print((pivot * 100).round(1).to_string())

    fig_celltype_disagreement(df)
    print("\nDone.")


if __name__ == "__main__":
    main()
