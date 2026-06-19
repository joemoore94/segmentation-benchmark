"""Cell type vs. agreement/disagreement: full-ROI spatial comparison.

Three-panel figure using 10x-native cell centroids (full 2mm × 2mm ROI):
  A) Cells coloured by annotated cell type.
  B) Same cells coloured by agreement (blue) or disagreement (red) against
     CellPose; unmatched 10x cells shown in gray.
  C) Disagree rate (%) by cell type for all six comparison methods —
     quantifies which types drive the spatial pattern in panels A and B.

Usage::

    conda run -n segbench python scripts/make_agreement_explainer.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

ROI_DIR = Path("data/processed/roi")
TABLES  = Path("results/tables")
FIGURES = Path("results/figures")

CLUSTER_ANNOTATIONS = {
    "0":  "Luminal epithelial",
    "1":  "Luminal epithelial",
    "2":  "Macrophages",
    "3":  "Luminal epithelial",
    "4":  "Myoepithelial",
    "5":  "T cells",
    "6":  "B cells",
    "7":  "Macrophages",
    "8":  "Luminal epithelial",
    "9":  "CAFs",
    "10": "Smooth muscle",
    "11": "Endothelial",
    "12": "Plasma cells",
    "13": "CAFs",
    "14": "Adipocytes",
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


def recluster_10x() -> ad.AnnData:
    sc.settings.verbosity = 0
    adata = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph")
    adata.obs["cell_type"] = adata.obs["leiden"].map(CLUSTER_ANNOTATIONS)
    return adata


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="poster")

    print("Loading data...")
    adata = recluster_10x()
    disagree = pd.read_csv(TABLES / "disagreement_table_10x_cellpose.csv")

    obs = adata.obs[["centroid_x", "centroid_y", "cell_type"]].copy()
    obs.index = obs.index.astype(str)

    disagree_indexed = disagree.set_index("id_a")[["disagree"]]
    disagree_indexed.index = disagree_indexed.index.astype(str)
    obs = obs.join(disagree_indexed, how="left")

    # Panel C data — per-cell-type disagree rate, all methods
    ct_disagree = pd.read_csv(TABLES / "celltype_disagreement.csv")
    method_order = ["CellPose", "StarDist", "Mesmer", "Voronoi (CP)", "Voronoi (M)", "Baysor"]
    pivot = ct_disagree.pivot(index="cell_type", columns="comparison", values="disagree_rate")
    pivot = pivot.reindex(index=list(CELLTYPE_COLORS.keys()), columns=method_order) * 100

    fig, axes = plt.subplots(1, 3, figsize=(30, 9),
                             gridspec_kw={"width_ratios": [2, 2, 2.2]})

    # Panel A — cell types
    for ct, color in CELLTYPE_COLORS.items():
        sub = obs[obs["cell_type"] == ct]
        if len(sub):
            axes[0].scatter(sub["centroid_x"], sub["centroid_y"],
                            c=color, s=3, alpha=0.6, rasterized=True)
    axes[0].set_title("A · Cell type (10x native)", fontweight="bold")
    axes[0].set_xlabel("x (µm)")
    axes[0].set_ylabel("y (µm)")
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    handles = [mpatches.Patch(color=CELLTYPE_COLORS[ct], label=ct) for ct in CELLTYPE_COLORS]
    axes[0].legend(handles=handles, fontsize=8, loc="upper right",
                   title="Cell type", title_fontsize=9, markerscale=2)

    # Panel B — agree / disagree vs. CellPose
    unmatched      = obs[obs["disagree"].isna()]
    agree          = obs[obs["disagree"] == 0.0]
    disagree_cells = obs[obs["disagree"] == 1.0]
    axes[1].scatter(unmatched["centroid_x"], unmatched["centroid_y"],
                    c="#CCCCCC", s=3, alpha=0.4, rasterized=True,
                    label=f"Unmatched ({len(unmatched):,})")
    axes[1].scatter(agree["centroid_x"], agree["centroid_y"],
                    c="#4C72B0", s=3, alpha=0.6, rasterized=True,
                    label=f"Agree ({len(agree):,})")
    axes[1].scatter(disagree_cells["centroid_x"], disagree_cells["centroid_y"],
                    c="#C44E52", s=3, alpha=0.6, rasterized=True,
                    label=f"Disagree ({len(disagree_cells):,})")
    axes[1].set_title("B · Agreement vs. CellPose", fontweight="bold")
    axes[1].set_xlabel("x (µm)")
    axes[1].set_ylabel("")
    axes[1].set_aspect("equal")
    axes[1].invert_yaxis()
    axes[1].legend(fontsize=9, loc="upper right", markerscale=3)

    # Panel C — disagree rate heatmap, all methods × cell types
    sns.heatmap(
        pivot, ax=axes[2],
        cmap="YlOrRd", vmin=0, vmax=80,
        annot=True, fmt=".0f",
        linewidths=0.5, linecolor="#dddddd",
        cbar_kws={"label": "Disagree rate (%)"},
        yticklabels=[
            f"▌ {ct}" for ct in pivot.index
        ],
    )
    for tick, ct in zip(axes[2].get_yticklabels(), pivot.index):
        tick.set_color(CELLTYPE_COLORS.get(ct, "black"))
    axes[2].set_title("C · Disagree rate (%) by cell type and method", fontweight="bold")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("")
    axes[2].tick_params(axis="x", rotation=35, labelsize=9)
    axes[2].tick_params(axis="y", rotation=0, labelsize=9)

    fig.suptitle(
        "Cell type vs. agreement with CellPose (10x native, full 2mm × 2mm ROI)",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "agreement_explainer.png", dpi=150)
    plt.close(fig)
    print("Saved agreement_explainer.png")


if __name__ == "__main__":
    main()
