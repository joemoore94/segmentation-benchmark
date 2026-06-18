"""Generate publication-quality figures from the cross-method comparison tables.

Reads ``results/tables/*`` (produced by ``run_comparison.py``) and the per-method
AnnData in ``data/processed/roi/`` and writes PNGs to ``results/figures/``.

Usage::

    conda run -n segbench python scripts/make_figures.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import numpy as np

from segbench.io import PIXEL_SIZE

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

# qv>=20 non-control transcripts in the full 2mm x 2mm ROI (see docs/dataset.md)
TOTAL_TRANSCRIPTS_FULL_ROI = 3_392_051

sns.set_theme(style="whitegrid", context="talk")


METHOD_COLORS = {
    "cellpose": "#4C72B0",
    "baysor": "#DD8452",
    "10x_native": "#55A868",
    "stardist": "#8172B2",
    "baysor_prior": "#937860",
    "voronoi": "#17BECF",
    "mesmer": "#D62728",
    "voronoi_mesmer": "#BCBD22",
}
METHOD_LABELS = {
    "cellpose": "CellPose",
    "baysor": "Baysor",
    "10x_native": "10x native",
    "stardist": "StarDist",
    "baysor_prior": "Baysor (prior)",
    "voronoi": "Voronoi (CP)",
    "mesmer": "Mesmer",
    "voronoi_mesmer": "Voronoi (M)",
}

# Methods included in the main analysis figures (baysor_prior is supplemental).
MAIN_METHODS = ["10x_native", "cellpose", "stardist", "mesmer", "voronoi", "voronoi_mesmer", "baysor"]

# The six pairwise comparisons shown in multi-panel figures, in 2×3 grid order.
COMPARISON_ORDER = [
    ("cellpose",       "CellPose"),
    ("stardist",       "StarDist"),
    ("mesmer",         "Mesmer"),
    ("voronoi",        "Voronoi (CP)"),
    ("voronoi_mesmer", "Voronoi (M)"),
    ("baysor",         "Baysor"),
]

# Keys used in the density_disagreement_summary.csv comparison column.
DENSITY_CSV_KEY = {
    "cellpose":       "10x native vs. CellPose",
    "stardist":       "10x native vs. StarDist",
    "mesmer":         "10x native vs. Mesmer",
    "voronoi":        "10x native vs. Voronoi",
    "voronoi_mesmer": "10x native vs. Voronoi (Mesmer)",
    "baysor":         "10x native vs. Baysor",
}


def fig_cell_counts_and_sizes() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts.csv", index_col="method")
    methods = [m for m in MAIN_METHODS if m in counts.index]

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adata_voronoi = ad.read_h5ad(ROI_DIR / "adata_voronoi.h5ad")
    adata_mesmer = ad.read_h5ad(ROI_DIR / "adata_mesmer.h5ad")
    adata_voronoi_mesmer = ad.read_h5ad(ROI_DIR / "adata_voronoi_mesmer.h5ad")

    transcripts_by_method = {
        "cellpose": np.asarray(adata_cellpose.X.sum(axis=1)).ravel(),
        "baysor": np.asarray(adata_baysor.X.sum(axis=1)).ravel(),
        "10x_native": np.asarray(adata_10x.X.sum(axis=1)).ravel(),
        "stardist": np.asarray(adata_stardist.X.sum(axis=1)).ravel(),
        "voronoi": np.asarray(adata_voronoi.X.sum(axis=1)).ravel(),
        "mesmer": np.asarray(adata_mesmer.X.sum(axis=1)).ravel(),
        "voronoi_mesmer": np.asarray(adata_voronoi_mesmer.X.sum(axis=1)).ravel(),
    }

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    axes[0].bar(
        [METHOD_LABELS[m] for m in methods],
        counts.loc[methods, "n_cells"].to_numpy(),
        color=[METHOD_COLORS[m] for m in methods],
    )
    axes[0].set_ylabel("Cell count (full 2mm × 2mm ROI)")
    axes[0].set_title("Cell count")
    axes[0].tick_params(axis="x", rotation=40)
    for label in axes[0].get_xticklabels():
        label.set_horizontalalignment("right")

    tx_long = pd.concat(
        [pd.DataFrame({"method": METHOD_LABELS[m], "log10_tx": np.log10(transcripts_by_method[m] + 1)})
         for m in methods],
        ignore_index=True,
    )
    sns.violinplot(
        data=tx_long, x="method", y="log10_tx",
        hue="method",
        palette={METHOD_LABELS[m]: METHOD_COLORS[m] for m in methods},
        order=[METHOD_LABELS[m] for m in methods],
        cut=0,
        inner="box",
        legend=False,
        ax=axes[1],
    )
    axes[1].set_xticks(range(len(methods)))
    axes[1].set_xticklabels([METHOD_LABELS[m] for m in methods], rotation=40, ha="right")
    axes[1].set_xlabel("")
    tick_vals = [1, 10, 100, 1000]
    axes[1].set_yticks([np.log10(v) for v in tick_vals])
    axes[1].set_yticklabels([str(v) for v in tick_vals])
    axes[1].set_ylabel("Transcripts per cell")
    axes[1].set_title("Transcripts/cell distribution")

    cellpose_area_um2 = adata_cellpose.obs["area"] * PIXEL_SIZE**2
    stardist_area_um2 = adata_stardist.obs["area"] * PIXEL_SIZE**2
    mesmer_area_um2 = adata_mesmer.obs["area"] * PIXEL_SIZE**2
    tenx_nucleus_area_um2 = adata_10x.obs["nucleus_area_um2"]
    sns.histplot(cellpose_area_um2, bins=50, ax=axes[2],
                 color=METHOD_COLORS["cellpose"], label=METHOD_LABELS["cellpose"], alpha=0.4)
    sns.histplot(stardist_area_um2, bins=50, ax=axes[2],
                 color=METHOD_COLORS["stardist"], label=METHOD_LABELS["stardist"], alpha=0.4)
    sns.histplot(mesmer_area_um2, bins=50, ax=axes[2],
                 color=METHOD_COLORS["mesmer"], label=METHOD_LABELS["mesmer"], alpha=0.4)
    sns.histplot(tenx_nucleus_area_um2, bins=50, ax=axes[2],
                 color=METHOD_COLORS["10x_native"], label=METHOD_LABELS["10x_native"], alpha=0.4)
    axes[2].set_xlabel("Nucleus area (µm²)")
    axes[2].set_title("Nuclear mask size (nuclear methods + 10x native)")
    axes[2].legend()

    fig.suptitle("Cell count and QC: all methods (full 2mm × 2mm ROI)")
    capture = counts.loc[methods, "transcript_capture_rate"]
    capture_str = ", ".join(f"{METHOD_LABELS[m]} {capture[m]:.0%}" for m in methods)
    fig.text(
        0.5, 0.01,
        f"Transcript capture: {capture_str}",
        ha="center", fontsize=10, style="italic",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIGURES_DIR / "cell_counts_and_sizes.png", dpi=150)
    plt.close(fig)


def fig_expression_correlation() -> None:
    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"expression_correlation_10x_{m}.csv"))
        for m, label in COMPARISON_ORDER
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=True)
    for ax, (m, label, corr) in zip(axes.flatten(), pairs):
        median = corr["correlation"].median()
        sns.histplot(corr["correlation"].dropna(), bins=40, ax=ax, color=METHOD_COLORS[m])
        ax.axvline(median, color="black", linestyle="--", label=f"median = {median:.3f}")
        ax.set_xlabel("Pearson correlation")
        ax.set_title(f"10x native vs. {label}")
        ax.legend(fontsize=9)

    axes[0, 0].set_ylabel("Number of pairs")
    axes[1, 0].set_ylabel("Number of pairs")
    fig.suptitle("Per-cell expression agreement vs. 10x native (matched cell pairs)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "expression_correlation.png", dpi=150)
    plt.close(fig)


def fig_disagreement_spatial_map() -> None:
    pairs = [
        (label, pd.read_csv(TABLES_DIR / f"disagreement_table_10x_{m}.csv"))
        for m, label in COMPARISON_ORDER
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 13))
    for ax, (label, df) in zip(axes.flatten(), pairs):
        sns.scatterplot(
            data=df, x="centroid_x", y="centroid_y",
            hue="disagree",
            palette={0.0: "#4C72B0", 1.0: "#C44E52"},
            s=4, alpha=0.6, ax=ax, legend=False,
        )
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(f"10x native vs. {label}")
        handles = [
            mpatches.Patch(color="#4C72B0", label="Agree"),
            mpatches.Patch(color="#C44E52", label="Disagree"),
        ]
        ax.legend(handles=handles, fontsize=9)

    fig.suptitle("Cell-type agreement (blue) vs. disagreement (red)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "disagreement_spatial_map.png", dpi=150)
    plt.close(fig)


def fig_cell_type_confusion() -> None:
    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"cell_type_confusion_10x_{m}.csv", index_col=0))
        for m, label in COMPARISON_ORDER
    ]

    fig, axes = plt.subplots(2, 3, figsize=(20, 14))
    for ax, (m, label, confusion) in zip(axes.flatten(), pairs):
        sns.heatmap(confusion, annot=False, cmap="viridis", ax=ax)
        ax.set_xlabel(f"{label} Leiden cluster")
        ax.set_ylabel("10x native Leiden cluster")
        ax.set_title(f"10x native vs. {label}")

    fig.suptitle("Cell-type cluster correspondence (matched pairs, Hungarian-aligned labels)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "cell_type_confusion.png", dpi=150)
    plt.close(fig)


def fig_density_vs_disagreement() -> None:
    log_density = pd.read_csv(TABLES_DIR / "10x_log_density.csv", index_col=0)["log_density"]
    summary = pd.read_csv(TABLES_DIR / "density_disagreement_summary.csv", index_col="comparison")

    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"disagreement_table_10x_{m}.csv"))
        for m, label in COMPARISON_ORDER
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
    for ax, (m, label, df) in zip(axes.flatten(), pairs):
        df = df.copy()
        df["log_density"] = df["id_a"].map(log_density)
        sns.kdeplot(
            data=df, x="log_density", hue="disagree",
            palette={0.0: "#4C72B0", 1.0: "#C44E52"}, common_norm=False,
            fill=True, alpha=0.3, ax=ax,
        )
        medians = df.groupby("disagree")["log_density"].median()
        ax.axvline(medians[0.0], color="#4C72B0", linestyle="--")
        ax.axvline(medians[1.0], color="#C44E52", linestyle="--")
        csv_key = DENSITY_CSV_KEY[m]
        p = summary.loc[csv_key, "p_value"] if csv_key in summary.index else float("nan")
        ax.set_title(f"10x native vs. {label}\n(p = {p:.1e})")
        ax.set_xlabel("Phenotypic log-density (Mellon)")
        ax.legend(title="Disagree", labels=["Yes", "No"], fontsize=9)

    axes[0, 0].set_ylabel("Density")
    axes[1, 0].set_ylabel("Density")
    fig.suptitle("10x native phenotypic density (Mellon) vs. cell-type call disagreement")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "density_vs_disagreement.png", dpi=150)
    plt.close(fig)


def fig_pca_umap() -> None:
    methods = [m for m in MAIN_METHODS if m != "10x_native"] + ["10x_native"]
    embeddings = {m: pd.read_csv(TABLES_DIR / f"embedding_{m}.csv", index_col=0) for m in methods}

    fig, axes = plt.subplots(2, len(methods), figsize=(42, 10))
    for col, method in enumerate(methods):
        emb = embeddings[method]
        n_clusters = emb["leiden"].nunique()
        palette = sns.color_palette("tab20", n_clusters)

        sns.scatterplot(
            data=emb, x="PC1", y="PC2", hue="leiden", palette=palette,
            s=10, alpha=0.6, ax=axes[0, col], legend=False,
        )
        axes[0, col].set_title(f"{METHOD_LABELS[method]}\n({n_clusters} clusters)")
        axes[0, col].set_xlabel("PC1")

        sns.scatterplot(
            data=emb, x="UMAP1", y="UMAP2", hue="leiden", palette=palette,
            s=10, alpha=0.6, ax=axes[1, col], legend=False,
        )
        axes[1, col].set_xlabel("UMAP1")

    axes[0, 0].set_ylabel("PC2")
    axes[1, 0].set_ylabel("UMAP2")
    fig.suptitle("Per-method Leiden clustering: PCA (top) and UMAP (bottom)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "pca_umap_clusters.png", dpi=150)
    plt.close(fig)


def fig_local_morans_map() -> None:
    pairs = [
        (label, pd.read_csv(TABLES_DIR / f"local_morans_10x_{m}.csv"))
        for m, label in COMPARISON_ORDER
    ]
    LISA_COLORS = {"HH": "#C44E52", "LL": "#4C72B0", "HL": "#DD8452", "LH": "#CCB974"}

    fig, axes = plt.subplots(2, 3, figsize=(18, 13))
    for ax, (label, df) in zip(axes.flatten(), pairs):
        for cluster, color in LISA_COLORS.items():
            sub = df[df["lisa_cluster"] == cluster]
            ax.scatter(sub["centroid_x"], sub["centroid_y"], c=color, s=4, alpha=0.6, label=cluster)
        ax.set_title(f"10x native vs. {label}")
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.legend(title="LISA", markerscale=2, fontsize=9)

    fig.suptitle("Local Moran's I: HH = disagreement hotspot, LL = agreement coldspot")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "local_morans_map.png", dpi=150)
    plt.close(fig)


def fig_de_volcano() -> None:
    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"de_disagree_10x_{m}.csv"))
        for m, label in COMPARISON_ORDER
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 11), sharey=True)
    for ax, (m, label, df) in zip(axes.flatten(), pairs):
        sig = df["pvals_adj"] < 0.05
        ax.scatter(
            df.loc[~sig, "logfoldchanges"], -np.log10(df.loc[~sig, "pvals_adj"] + 1e-300),
            c="#AAAAAA", s=8, alpha=0.5, label="n.s.",
        )
        ax.scatter(
            df.loc[sig, "logfoldchanges"], -np.log10(df.loc[sig, "pvals_adj"] + 1e-300),
            c="#C44E52", s=10, alpha=0.7, label="adj. p < 0.05",
        )
        for _, row in df[sig].nlargest(5, "scores").iterrows():
            ax.annotate(
                row["names"],
                xy=(row["logfoldchanges"], -np.log10(row["pvals_adj"] + 1e-300)),
                fontsize=7, ha="left",
            )
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel("log fold change (disagree vs. agree)")
        ax.set_title(f"10x native vs. {label}")
        ax.legend(fontsize=8)

    axes[0, 0].set_ylabel("-log10(adj. p)")
    axes[1, 0].set_ylabel("-log10(adj. p)")
    fig.suptitle("DE: disagree vs. agree cells (Wilcoxon, 10x native cells)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "de_volcano.png", dpi=150)
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_cell_counts_and_sizes()
    fig_expression_correlation()
    fig_disagreement_spatial_map()
    fig_cell_type_confusion()
    fig_density_vs_disagreement()
    fig_pca_umap()
    lisa_files = [f"local_morans_10x_{m}.csv" for m, _ in COMPARISON_ORDER]
    if all((TABLES_DIR / f).exists() for f in lisa_files):
        fig_local_morans_map()
    de_files = [f"de_disagree_10x_{m}.csv" for m, _ in COMPARISON_ORDER]
    if all((TABLES_DIR / f).exists() for f in de_files):
        fig_de_volcano()
    print(f"wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
