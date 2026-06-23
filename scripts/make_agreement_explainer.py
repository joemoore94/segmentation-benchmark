"""Cell type vs. agreement/disagreement: full-ROI spatial comparison.

Layout (7 rows × 2 columns):
  Row 1 — reference:
    A) 10x native cells coloured by annotated cell type.
    B) Average disagree rate by cell type across all six methods.
  Rows 2–7 — one row per comparison method:
    Left)  Spatial agree/disagree map (10x native centroids).
    Right) Disagree rate (%) by cell type for that method.

Usage::

    conda run -n segbench python scripts/make_agreement_explainer.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from segbench.constants import CELLTYPE_COLORS, CLUSTER_ANNOTATIONS, COMPARISON_ORDER
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES  = Path("results/tables")
FIGURES = Path("results/figures")


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
    apply_style()

    print("Loading data...")
    adata = recluster_10x()

    obs_full = adata.obs[["centroid_x", "centroid_y", "cell_type"]].copy()
    obs_full.index = obs_full.index.astype(str)

    ct_disagree = pd.read_csv(TABLES / "celltype_disagreement.csv")

    # Cell type order: sort by average disagree rate ascending so highest is at top of barh
    avg_rate = ct_disagree.groupby("cell_type")["disagree_rate"].mean()
    ct_sorted = avg_rate.sort_values(ascending=True).index.tolist()
    bar_colors = [CELLTYPE_COLORS[ct] for ct in ct_sorted]

    # ---------------------------------------------------------------- layout
    n_methods = len(COMPARISON_ORDER)
    fig = plt.figure(figsize=(18, 8 + 5.4 * n_methods))
    gs = gridspec.GridSpec(
        1 + n_methods, 2, figure=fig,
        height_ratios=[3] + [1.8] * n_methods,
        hspace=0.5, wspace=0.35,
    )
    ax_ct  = fig.add_subplot(gs[0, 0])   # Panel A: cell type map
    ax_avg = fig.add_subplot(gs[0, 1])   # Panel B: average disagree rate

    method_rows = [
        (fig.add_subplot(gs[1 + i, 0]), fig.add_subplot(gs[1 + i, 1]))
        for i in range(n_methods)
    ]

    # ---------------------------------------------------------------- Panel A: cell type map
    for ct, color in CELLTYPE_COLORS.items():
        sub = obs_full[obs_full["cell_type"] == ct]
        if len(sub):
            ax_ct.scatter(sub["centroid_x"], sub["centroid_y"],
                          c=color, s=2, alpha=0.5, rasterized=True)
    ax_ct.set_title("Cell type (10x native)", fontweight="bold")
    ax_ct.set_xlabel("x (µm)")
    ax_ct.set_ylabel("y (µm)")
    ax_ct.set_aspect("equal")
    ax_ct.invert_yaxis()
    handles = [mpatches.Patch(color=CELLTYPE_COLORS[ct], label=ct)
               for ct in CELLTYPE_COLORS]
    ax_ct.legend(handles=handles, fontsize=10, loc="upper right",
                 title="Cell type", title_fontsize=10, markerscale=2)

    # ---------------------------------------------------------------- Panel B: average disagree rate
    avg_vals = avg_rate.reindex(ct_sorted).fillna(0) * 100
    ax_avg.barh(ct_sorted, avg_vals.values, color=bar_colors, edgecolor="white", height=0.7)
    ax_avg.set_xlim(0, 85)
    ax_avg.axvline(50, color="black", linewidth=0.8, linestyle="--", alpha=0.35)
    ax_avg.set_xlabel("Disagree rate (%)")
    ax_avg.set_title("Average disagree rate by cell type\n(all methods)", fontweight="bold")
    ax_avg.tick_params(axis="y", labelsize=12)
    for val, ct in zip(avg_vals.values, ct_sorted):
        if val > 2:
            ax_avg.text(val + 1, ct, f"{val:.0f}%", va="center", fontsize=11)

    # ---------------------------------------------------------------- Method rows
    panel_letters = "CDEFGHIJKLMNOP"[:n_methods]
    for (ax_sp, ax_bar), (m, label), letter in zip(method_rows, COMPARISON_ORDER, panel_letters):
        # Load disagree table for this method
        df = pd.read_csv(TABLES / f"disagreement_table_10x_{m}.csv")
        df_indexed = df.set_index("id_a")[["disagree"]]
        df_indexed.index = df_indexed.index.astype(str)
        obs = obs_full.join(df_indexed, how="left")

        # Spatial map
        unmatched = obs[obs["disagree"].isna()]
        agree     = obs[obs["disagree"] == 0.0]
        dis       = obs[obs["disagree"] == 1.0]
        ax_sp.scatter(unmatched["centroid_x"], unmatched["centroid_y"],
                      c="#CCCCCC", s=2, alpha=0.35, rasterized=True)
        ax_sp.scatter(agree["centroid_x"], agree["centroid_y"],
                      c="#4C72B0", s=2, alpha=0.5, rasterized=True)
        ax_sp.scatter(dis["centroid_x"], dis["centroid_y"],
                      c="#C44E52", s=2, alpha=0.5, rasterized=True)
        ax_sp.set_title(label, fontweight="bold")
        ax_sp.set_xlabel("x (µm)")
        ax_sp.set_ylabel("y (µm)")
        ax_sp.invert_yaxis()
        dis_rate = df["disagree"].mean() * 100
        ax_sp.text(0.02, 0.97, f"Disagree: {dis_rate:.1f}%",
                   transform=ax_sp.transAxes, fontsize=11,
                   va="top", ha="left",
                   bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

        # Bar chart
        rates = (
            ct_disagree[ct_disagree["comparison"] == label]
            .set_index("cell_type")["disagree_rate"]
            .reindex(ct_sorted).fillna(0) * 100
        )
        ax_bar.barh(ct_sorted, rates.values, color=bar_colors, edgecolor="white", height=0.7)
        ax_bar.set_xlim(0, 85)
        ax_bar.axvline(50, color="black", linewidth=0.8, linestyle="--", alpha=0.35)
        ax_bar.set_xlabel("Disagree rate (%)")
        ax_bar.set_title("Disagree rate by cell type", fontweight="bold")
        ax_bar.tick_params(axis="y", labelsize=12)
        for val, ct in zip(rates.values, ct_sorted):
            if val > 2:
                ax_bar.text(val + 1, ct, f"{val:.0f}%", va="center", fontsize=11)

    fig.legend(handles=[
        mpatches.Patch(color="#CCCCCC", label="Unmatched"),
        mpatches.Patch(color="#4C72B0", label="Agree"),
        mpatches.Patch(color="#C44E52", label="Disagree"),
    ], loc="lower center", ncols=3, fontsize=11, framealpha=0.9)
    fig.subplots_adjust(top=0.99, bottom=0.04)
    fig.savefig(FIGURES / "agreement_explainer.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved agreement_explainer.png")


if __name__ == "__main__":
    main()
