"""Generate a figure that illustrates how agreement/disagreement is determined.

Four-panel figure:
  A) 10x native cells in a tissue patch, coloured by annotated cell type.
  B) CellPose cells in the same patch, coloured by their Hungarian-aligned
     cell type (so colours mean the same thing as in A).
  C) Agreement map: same patch, 10x native cells coloured green (agree) or
     red (disagree) after Hungarian label alignment.
  D) Label-matching scatter across ALL matched pairs: x = 10x native cluster,
     y = CellPose aligned cluster, diagonal = agreement, off-diagonal =
     disagreement. Annotated with cell type names.

Usage::

    conda run -n segbench python scripts/make_agreement_explainer.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import linear_sum_assignment

ROI_DIR  = Path("data/processed/roi")
TABLES   = Path("results/tables")
FIGURES  = Path("results/figures")

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
    "Unmatched":          "#DDDDDD",
}

# 500 µm × 500 µm patch centred on the luminal epithelial disagreement
# hotspot (x≈691, y≈855 from build_agreement_explainer exploration).
PATCH_X = (450, 950)
PATCH_Y = (600, 1100)


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


def recompute_hungarian(raw_confusion: pd.DataFrame) -> dict[str, str]:
    """Return {comp_cluster_id: matched_10x_cluster_id} from raw confusion matrix."""
    mat = raw_confusion.to_numpy().astype(float)
    row_ind, col_ind = linear_sum_assignment(-mat)
    return {
        raw_confusion.columns[c]: raw_confusion.index[r]
        for r, c in zip(row_ind, col_ind)
    }


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    adata_10x = recluster_10x()
    adata_cp  = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    disagree  = pd.read_csv(TABLES / "disagreement_table_10x_cellpose.csv")

    # Hungarian mapping: CellPose original cluster → matched 10x cluster → cell type
    raw_conf = pd.read_csv(TABLES / "cell_type_confusion_10x_cellpose.csv", index_col=0)
    raw_conf.index   = raw_conf.index.astype(str)
    raw_conf.columns = raw_conf.columns.astype(str)
    comp_to_10x = recompute_hungarian(raw_conf)   # {cp_cluster: 10x_cluster}

    # --- build 10x patch data ---
    obs_10x = adata_10x.obs[["centroid_x", "centroid_y", "leiden", "cell_type"]].copy()
    obs_10x.index = obs_10x.index.astype(str)
    in_patch_10x = (
        obs_10x["centroid_x"].between(*PATCH_X) &
        obs_10x["centroid_y"].between(*PATCH_Y)
    )
    patch_10x = obs_10x[in_patch_10x].copy()

    # Merge disagree info onto patch cells
    disagree_sub = disagree[disagree["id_a"].isin(patch_10x.index)]
    patch_10x = patch_10x.join(
        disagree_sub.set_index("id_a")[["disagree", "id_b", "label_b"]],
        how="left",
    )

    # --- build CellPose patch data ---
    obs_cp = adata_cp.obs[["centroid_x", "centroid_y"]].copy()
    obs_cp.index = obs_cp.index.astype(str)

    # We want CellPose cells in the same spatial region (matched or not).
    in_patch_cp = (
        obs_cp["centroid_x"].between(*PATCH_X) &
        obs_cp["centroid_y"].between(*PATCH_Y)
    )
    patch_cp = obs_cp[in_patch_cp].copy()

    # Assign cell type via their Hungarian-aligned label.
    # Matched cells: look up their original CellPose cluster from the matches CSV.
    matches = pd.read_csv(TABLES / "matches_10x_cellpose.csv")
    # matches has id_a (10x), id_b (cellpose)
    # We need the original (pre-alignment) CellPose cluster for each matched cell.
    # That comes from the embedding CSV for cellpose — but we don't have unaligned labels saved.
    # Use label_b from disagree table (which IS the aligned label) to get cell type.
    id_b_to_aligned = dict(zip(disagree["id_b"].astype(str), disagree["label_b"].astype(str)))
    patch_cp["aligned_label"] = patch_cp.index.map(id_b_to_aligned)
    patch_cp["cell_type"] = patch_cp["aligned_label"].map(
        lambda x: CLUSTER_ANNOTATIONS.get(str(int(float(x))), "Unmatched") if pd.notna(x) else "Unmatched"
    )

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))

    def set_patch_axes(ax, title):
        ax.set_xlim(*PATCH_X)
        ax.set_ylim(*PATCH_Y)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_title(title, fontsize=11, fontweight="bold")

    # Panel A — 10x native cell types
    for ct, color in CELLTYPE_COLORS.items():
        sub = patch_10x[patch_10x["cell_type"] == ct]
        if len(sub):
            axes[0].scatter(sub["centroid_x"], sub["centroid_y"],
                            c=color, s=28, alpha=0.85, label=ct, zorder=2)
    set_patch_axes(axes[0], "A · 10x native\n(annotated cell type)")
    handles = [mpatches.Patch(color=CELLTYPE_COLORS[ct], label=ct)
               for ct in CELLTYPE_COLORS if ct in patch_10x["cell_type"].values]
    axes[0].legend(handles=handles, fontsize=7, loc="upper right",
                   title="Cell type", title_fontsize=8, markerscale=1.5)

    # Panel B — CellPose cells, same spatial region, coloured by aligned cell type
    for ct, color in CELLTYPE_COLORS.items():
        sub = patch_cp[patch_cp["cell_type"] == ct]
        if len(sub):
            axes[1].scatter(sub["centroid_x"], sub["centroid_y"],
                            c=color, s=28, alpha=0.85, label=ct, zorder=2,
                            marker="^")
    set_patch_axes(axes[1], "B · CellPose\n(Hungarian-aligned cell type)")
    axes[1].scatter([], [], c="gray", s=28, marker="^", label="Unmatched")
    axes[1].legend(handles=handles, fontsize=7, loc="upper right",
                   title="Cell type", title_fontsize=8, markerscale=1.5)
    axes[1].set_ylabel("")

    # Panel C — agreement / disagreement map (10x native centroids)
    agree    = patch_10x[patch_10x["disagree"] == 0.0]
    disagree_p = patch_10x[patch_10x["disagree"] == 1.0]
    unmatched  = patch_10x[patch_10x["disagree"].isna()]
    axes[2].scatter(agree["centroid_x"],    agree["centroid_y"],
                    c="#4C72B0", s=28, alpha=0.85, label=f"Agree ({len(agree)})", zorder=2)
    axes[2].scatter(disagree_p["centroid_x"], disagree_p["centroid_y"],
                    c="#C44E52", s=28, alpha=0.85, label=f"Disagree ({len(disagree_p)})", zorder=3)
    axes[2].scatter(unmatched["centroid_x"],  unmatched["centroid_y"],
                    c="#AAAAAA", s=14, alpha=0.5, label=f"Unmatched ({len(unmatched)})", zorder=1)
    set_patch_axes(axes[2], "C · Agreement map\n(10x native centroids)")
    axes[2].legend(fontsize=8, loc="upper right")
    axes[2].set_ylabel("")

    # Panel D — label-matching scatter across ALL pairs
    all_labels = disagree[["label_a", "label_b", "disagree"]].dropna()
    # Jitter to show density
    rng = np.random.default_rng(42)
    jitter = 0.35
    x_j = all_labels["label_a"] + rng.uniform(-jitter, jitter, len(all_labels))
    y_j = all_labels["label_b"] + rng.uniform(-jitter, jitter, len(all_labels))
    colors_d = np.where(all_labels["disagree"] == 0, "#4C72B0", "#C44E52")
    axes[3].scatter(x_j, y_j, c=colors_d, s=4, alpha=0.15, rasterized=True)

    # Annotate cluster IDs with cell type on both axes
    ct_order = list(dict.fromkeys(CLUSTER_ANNOTATIONS.values()))
    cluster_ids = sorted(CLUSTER_ANNOTATIONS.keys(), key=int)
    present_in_a = sorted(all_labels["label_a"].unique())
    present_in_b = sorted(all_labels["label_b"].unique())
    for cid in present_in_a:
        ct = CLUSTER_ANNOTATIONS[str(int(cid))]
        axes[3].annotate(f"{int(cid)}\n{ct[:6]}", xy=(cid, -1.8),
                         ha="center", va="top", fontsize=5.5, rotation=0)
    for cid in present_in_b:
        ct = CLUSTER_ANNOTATIONS[str(int(cid))]
        axes[3].annotate(f"{int(cid)} {ct[:8]}", xy=(-1.5, cid),
                         ha="right", va="center", fontsize=5.5)

    # Diagonal line = perfect agreement
    diag_max = max(all_labels["label_a"].max(), all_labels["label_b"].max())
    axes[3].plot([0, diag_max], [0, diag_max], "k--", lw=1, alpha=0.5, label="Perfect agreement")
    axes[3].set_xlabel("10x native cluster label", fontsize=9)
    axes[3].set_ylabel("CellPose cluster label\n(after Hungarian alignment)", fontsize=9)
    axes[3].set_title("D · Label matching: all pairs\ndiagonal = agree, off-diagonal = disagree",
                      fontsize=11, fontweight="bold")
    agree_n    = (all_labels["disagree"] == 0).sum()
    disagree_n = (all_labels["disagree"] == 1).sum()
    axes[3].legend(handles=[
        mpatches.Patch(color="#4C72B0", label=f"Agree ({agree_n:,})"),
        mpatches.Patch(color="#C44E52", label=f"Disagree ({disagree_n:,})"),
    ], fontsize=8, loc="upper left")
    axes[3].set_xlim(-2, diag_max + 0.5)
    axes[3].set_ylim(-2.5, diag_max + 0.5)
    axes[3].tick_params(labelbottom=False, labelleft=False)

    fig.suptitle(
        "How agreement/disagreement is determined  (10x native vs. CellPose, 500 µm × 500 µm patch + all pairs)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "agreement_explainer.png", dpi=150)
    plt.close(fig)
    print("Saved agreement_explainer.png")


if __name__ == "__main__":
    main()
