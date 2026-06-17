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
}
METHOD_LABELS = {
    "cellpose": "CellPose",
    "baysor": "Baysor",
    "10x_native": "10x native",
    "stardist": "StarDist",
    "baysor_prior": "Baysor (prior)",
    "voronoi": "Voronoi",
}


def fig_cell_counts_and_sizes() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts.csv", index_col="method")
    methods = list(counts.index)

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adata_baysor_prior = ad.read_h5ad(ROI_DIR / "adata_baysor_prior.h5ad")
    adata_voronoi = ad.read_h5ad(ROI_DIR / "adata_voronoi.h5ad")

    cellpose_area_um2 = adata_cellpose.obs["area"] * PIXEL_SIZE**2
    tenx_nucleus_area_um2 = adata_10x.obs["nucleus_area_um2"]
    stardist_area_um2 = adata_stardist.obs["area"] * PIXEL_SIZE**2
    cellpose_transcripts = np.asarray(adata_cellpose.X.sum(axis=1)).ravel()
    baysor_transcripts = np.asarray(adata_baysor.X.sum(axis=1)).ravel()
    tenx_transcripts = np.asarray(adata_10x.X.sum(axis=1)).ravel()
    stardist_transcripts = np.asarray(adata_stardist.X.sum(axis=1)).ravel()
    baysor_prior_transcripts = np.asarray(adata_baysor_prior.X.sum(axis=1)).ravel()
    voronoi_transcripts = np.asarray(adata_voronoi.X.sum(axis=1)).ravel()
    transcripts_by_method = {
        "cellpose": cellpose_transcripts,
        "baysor": baysor_transcripts,
        "10x_native": tenx_transcripts,
        "stardist": stardist_transcripts,
        "baysor_prior": baysor_prior_transcripts,
        "voronoi": voronoi_transcripts,
    }

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].bar(
        [METHOD_LABELS[m] for m in methods],
        counts.loc[methods, "n_cells"].to_numpy(),
        color=[METHOD_COLORS[m] for m in methods],
    )
    axes[0].set_ylabel("Cell count (full 2mm × 2mm ROI)")
    axes[0].set_title("Cell count")
    axes[0].tick_params(axis="x", rotation=30)
    for label in axes[0].get_xticklabels():
        label.set_horizontalalignment("right")

    # Transcripts/cell is computed identically for all methods (sum of the
    # per-cell gene-count matrix), so this panel is a true apples-to-apples
    # QC comparison -- unlike "cell size", which means different things
    # (nucleus pixel area vs. transcript count) per method.
    tx_long = pd.concat(
        [pd.DataFrame({"method": m, "log10_tx": np.log10(transcripts_by_method[m] + 1)})
         for m in methods],
        ignore_index=True,
    )
    sns.violinplot(
        data=tx_long, x="method", y="log10_tx",
        hue="method",
        palette={m: METHOD_COLORS[m] for m in methods},
        order=methods,
        cut=0,
        inner="box",
        legend=False,
        ax=axes[1],
    )
    axes[1].set_xticks(range(len(methods)))
    axes[1].set_xticklabels([METHOD_LABELS[m] for m in methods], rotation=30, ha="right")
    axes[1].set_xlabel("")
    tick_vals = [1, 10, 100, 1000]
    axes[1].set_yticks([np.log10(v) for v in tick_vals])
    axes[1].set_yticklabels([str(v) for v in tick_vals])
    axes[1].set_ylabel("Transcripts per cell")
    axes[1].set_title("Transcripts/cell distribution")

    # Nucleus area is directly comparable across CellPose and StarDist (pixel
    # area converted to µm²) and 10x native (already in µm²) -- all three are
    # nuclear segmentation areas, so this checks whether the two nuclear
    # U-Net/star-convex methods are similarly sized to the platform's own
    # reference segmentation.
    sns.histplot(
        cellpose_area_um2, bins=50, ax=axes[2],
        color=METHOD_COLORS["cellpose"], label=METHOD_LABELS["cellpose"], alpha=0.4,
    )
    sns.histplot(
        stardist_area_um2, bins=50, ax=axes[2],
        color=METHOD_COLORS["stardist"], label=METHOD_LABELS["stardist"], alpha=0.4,
    )
    sns.histplot(
        tenx_nucleus_area_um2, bins=50, ax=axes[2],
        color=METHOD_COLORS["10x_native"], label=METHOD_LABELS["10x_native"], alpha=0.4,
    )
    axes[2].set_xlabel("Nucleus area (µm²)")
    axes[2].set_title("Nucleus area: CellPose vs. StarDist vs. 10x native (QC)")
    axes[2].legend()

    fig.suptitle(
        "Cell count and QC: CellPose vs. Baysor vs. 10x native vs. StarDist (full 2mm × 2mm ROI)"
    )
    capture = counts.loc[methods, "transcript_capture_rate"]
    capture_str = ", ".join(
        f"{METHOD_LABELS[m]} {capture[m]:.0%}" for m in methods
    )
    fig.text(
        0.5, 0.01,
        f"Transcript capture rate (assigned / {TOTAL_TRANSCRIPTS_FULL_ROI:,} total"
        f" qv≥20 transcripts in ROI): {capture_str}. CellPose and StarDist are"
        " nuclear-only, so cytoplasmic transcripts are not assigned to any cell.",
        ha="center", fontsize=11, style="italic",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIGURES_DIR / "cell_counts_and_sizes.png", dpi=150)
    plt.close(fig)


def fig_expression_correlation() -> None:
    corr_cellpose = pd.read_csv(TABLES_DIR / "expression_correlation_10x_cellpose.csv")
    corr_stardist = pd.read_csv(TABLES_DIR / "expression_correlation_10x_stardist.csv")
    corr_voronoi = pd.read_csv(TABLES_DIR / "expression_correlation_10x_voronoi.csv")
    corr_baysor = pd.read_csv(TABLES_DIR / "expression_correlation_10x_baysor.csv")
    corr_baysor_prior = pd.read_csv(TABLES_DIR / "expression_correlation_10x_baysor_prior.csv")

    fig, axes = plt.subplots(1, 5, figsize=(32, 5.5), sharey=True)
    for ax, corr, label, color in [
        (axes[0], corr_cellpose, "10x native vs. CellPose", METHOD_COLORS["cellpose"]),
        (axes[1], corr_stardist, "10x native vs. StarDist", METHOD_COLORS["stardist"]),
        (axes[2], corr_voronoi, "10x native vs. Voronoi", METHOD_COLORS["voronoi"]),
        (axes[3], corr_baysor, "10x native vs. Baysor", METHOD_COLORS["baysor"]),
        (axes[4], corr_baysor_prior, "10x native vs. Baysor (prior)", METHOD_COLORS["baysor_prior"]),
    ]:
        median = corr["correlation"].median()
        sns.histplot(corr["correlation"].dropna(), bins=40, ax=ax, color=color)
        ax.axvline(median, color="black", linestyle="--", label=f"median = {median:.2f}")
        ax.set_xlabel("Pearson correlation (matched cell pair)")
        ax.set_title(label)
        ax.legend()

    axes[0].set_ylabel("Number of pairs")
    fig.suptitle("Per-cell expression agreement vs. 10x native")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "expression_correlation.png", dpi=150)
    plt.close(fig)


def fig_disagreement_spatial_map() -> None:
    disagreement_cellpose = pd.read_csv(TABLES_DIR / "disagreement_table_10x_cellpose.csv")
    disagreement_stardist = pd.read_csv(TABLES_DIR / "disagreement_table_10x_stardist.csv")
    disagreement_voronoi = pd.read_csv(TABLES_DIR / "disagreement_table_10x_voronoi.csv")
    disagreement_baysor = pd.read_csv(TABLES_DIR / "disagreement_table_10x_baysor.csv")
    disagreement_baysor_prior = pd.read_csv(TABLES_DIR / "disagreement_table_10x_baysor_prior.csv")

    fig, axes = plt.subplots(1, 5, figsize=(30, 6))
    for ax, disagreement, label in [
        (axes[0], disagreement_cellpose, "10x native vs. CellPose"),
        (axes[1], disagreement_stardist, "10x native vs. StarDist"),
        (axes[2], disagreement_voronoi, "10x native vs. Voronoi"),
        (axes[3], disagreement_baysor, "10x native vs. Baysor"),
        (axes[4], disagreement_baysor_prior, "10x native vs. Baysor (prior)"),
    ]:
        sns.scatterplot(
            data=disagreement,
            x="centroid_x",
            y="centroid_y",
            hue="disagree",
            palette={0.0: "#4C72B0", 1.0: "#C44E52"},
            s=20,
            alpha=0.7,
            ax=ax,
            legend=False,
        )
        ax.set_xlabel("x (µm, ROI coordinates)")
        ax.set_ylabel("y (µm, ROI coordinates)")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(label)
        handles = [
            mpatches.Patch(color="#4C72B0", label="No"),
            mpatches.Patch(color="#C44E52", label="Yes"),
        ]
        ax.legend(handles=handles, title="Disagree")

    fig.suptitle("Cell-type agreement (blue) vs. disagreement (red)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "disagreement_spatial_map.png", dpi=150)
    plt.close(fig)


def fig_cell_type_confusion() -> None:
    confusion_cellpose = pd.read_csv(TABLES_DIR / "cell_type_confusion_10x_cellpose.csv", index_col=0)
    confusion_stardist = pd.read_csv(TABLES_DIR / "cell_type_confusion_10x_stardist.csv", index_col=0)
    confusion_voronoi = pd.read_csv(TABLES_DIR / "cell_type_confusion_10x_voronoi.csv", index_col=0)
    confusion_baysor = pd.read_csv(TABLES_DIR / "cell_type_confusion_10x_baysor.csv", index_col=0)
    confusion_baysor_prior = pd.read_csv(
        TABLES_DIR / "cell_type_confusion_10x_baysor_prior.csv", index_col=0
    )

    fig, axes = plt.subplots(1, 5, figsize=(32, 6))

    sns.heatmap(confusion_cellpose, annot=False, cmap="viridis", ax=axes[0])
    axes[0].set_xlabel("CellPose Leiden cluster")
    axes[0].set_ylabel("10x native Leiden cluster")
    axes[0].set_title("10x native vs. CellPose")

    sns.heatmap(confusion_stardist, annot=False, cmap="viridis", ax=axes[1])
    axes[1].set_xlabel("StarDist Leiden cluster")
    axes[1].set_ylabel("10x native Leiden cluster")
    axes[1].set_title("10x native vs. StarDist")

    sns.heatmap(confusion_voronoi, annot=False, cmap="viridis", ax=axes[2])
    axes[2].set_xlabel("Voronoi Leiden cluster")
    axes[2].set_ylabel("10x native Leiden cluster")
    axes[2].set_title("10x native vs. Voronoi")

    sns.heatmap(confusion_baysor, annot=False, cmap="viridis", ax=axes[3])
    axes[3].set_xlabel("Baysor Leiden cluster")
    axes[3].set_ylabel("10x native Leiden cluster")
    axes[3].set_title("10x native vs. Baysor")

    sns.heatmap(confusion_baysor_prior, annot=False, cmap="viridis", ax=axes[4])
    axes[4].set_xlabel("Baysor (prior) Leiden cluster")
    axes[4].set_ylabel("10x native Leiden cluster")
    axes[4].set_title("10x native vs. Baysor (prior)")

    fig.suptitle("Cell-type cluster correspondence (matched pairs, Hungarian-aligned labels)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "cell_type_confusion.png", dpi=150)
    plt.close(fig)


def fig_density_vs_disagreement() -> None:
    log_density = pd.read_csv(TABLES_DIR / "10x_log_density.csv", index_col=0)["log_density"]
    summary = pd.read_csv(TABLES_DIR / "density_disagreement_summary.csv", index_col="comparison")

    disagreement_cellpose = pd.read_csv(TABLES_DIR / "disagreement_table_10x_cellpose.csv")
    disagreement_stardist = pd.read_csv(TABLES_DIR / "disagreement_table_10x_stardist.csv")
    disagreement_voronoi = pd.read_csv(TABLES_DIR / "disagreement_table_10x_voronoi.csv")
    disagreement_baysor = pd.read_csv(TABLES_DIR / "disagreement_table_10x_baysor.csv")
    disagreement_baysor_prior = pd.read_csv(TABLES_DIR / "disagreement_table_10x_baysor_prior.csv")

    fig, axes = plt.subplots(1, 5, figsize=(30, 5), sharex=True, sharey=True)
    for ax, disagreement, label in [
        (axes[0], disagreement_cellpose, "10x native vs. CellPose"),
        (axes[1], disagreement_stardist, "10x native vs. StarDist"),
        (axes[2], disagreement_voronoi, "10x native vs. Voronoi"),
        (axes[3], disagreement_baysor, "10x native vs. Baysor"),
        (axes[4], disagreement_baysor_prior, "10x native vs. Baysor (prior)"),
    ]:
        disagreement = disagreement.copy()
        disagreement["log_density"] = disagreement["id_a"].map(log_density)

        sns.kdeplot(
            data=disagreement, x="log_density", hue="disagree",
            palette={0.0: "#4C72B0", 1.0: "#C44E52"}, common_norm=False,
            fill=True, alpha=0.3, ax=ax,
        )
        medians = disagreement.groupby("disagree")["log_density"].median()
        ax.axvline(medians[0.0], color="#4C72B0", linestyle="--")
        ax.axvline(medians[1.0], color="#C44E52", linestyle="--")

        p = summary.loc[label, "p_value"]
        ax.set_title(f"{label}\n(Mann-Whitney p = {p:.1e})")
        ax.set_xlabel("10x native phenotypic log-density (Mellon)")
        ax.legend(title="Disagree", labels=["Yes", "No"])

    axes[0].set_ylabel("Density")
    fig.suptitle("10x native phenotypic density (Mellon) vs. cell-type call disagreement")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "density_vs_disagreement.png", dpi=150)
    plt.close(fig)


def fig_pca_umap() -> None:
    methods = ["cellpose", "stardist", "voronoi", "baysor", "baysor_prior", "10x_native"]
    embeddings = {m: pd.read_csv(TABLES_DIR / f"embedding_{m}.csv", index_col=0) for m in methods}

    fig, axes = plt.subplots(2, 6, figsize=(32, 10))
    for col, method in enumerate(methods):
        emb = embeddings[method]
        n_clusters = emb["leiden"].nunique()
        palette = sns.color_palette("tab20", n_clusters)

        sns.scatterplot(
            data=emb, x="PC1", y="PC2", hue="leiden", palette=palette,
            s=8, alpha=0.6, ax=axes[0, col], legend=False,
        )
        axes[0, col].set_title(f"{METHOD_LABELS[method]}: PCA")

        sns.scatterplot(
            data=emb, x="UMAP1", y="UMAP2", hue="leiden", palette=palette,
            s=8, alpha=0.6, ax=axes[1, col], legend=False,
        )
        axes[1, col].set_title(f"{METHOD_LABELS[method]}: UMAP")
        axes[1, col].set_xlabel(f"UMAP1 ({n_clusters} Leiden clusters)")

    fig.suptitle("Per-method Leiden clustering: PCA and UMAP embeddings")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "pca_umap_clusters.png", dpi=150)
    plt.close(fig)


def fig_local_morans_map() -> None:
    files = {
        "10x native vs. CellPose": "local_morans_10x_cellpose.csv",
        "10x native vs. StarDist": "local_morans_10x_stardist.csv",
        "10x native vs. Voronoi": "local_morans_10x_voronoi.csv",
        "10x native vs. Baysor": "local_morans_10x_baysor.csv",
        "10x native vs. Baysor (prior)": "local_morans_10x_baysor_prior.csv",
    }
    LISA_COLORS = {"HH": "#C44E52", "LL": "#4C72B0", "HL": "#DD8452", "LH": "#CCB974"}

    fig, axes = plt.subplots(1, 5, figsize=(30, 6))
    for ax, (label, fname) in zip(axes, files.items()):
        df = pd.read_csv(TABLES_DIR / fname)
        for cluster, color in LISA_COLORS.items():
            sub = df[df["lisa_cluster"] == cluster]
            ax.scatter(sub["centroid_x"], sub["centroid_y"], c=color, s=6, alpha=0.6, label=cluster)
        ax.set_title(label)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.legend(title="LISA", markerscale=2, fontsize=8)

    fig.suptitle(
        "Local Moran's I clusters: HH = disagreement hotspot, LL = agreement coldspot"
    )
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "local_morans_map.png", dpi=150)
    plt.close(fig)


def fig_de_volcano() -> None:
    files = {
        "10x native vs. CellPose": "de_disagree_10x_cellpose.csv",
        "10x native vs. StarDist": "de_disagree_10x_stardist.csv",
        "10x native vs. Voronoi": "de_disagree_10x_voronoi.csv",
        "10x native vs. Baysor": "de_disagree_10x_baysor.csv",
        "10x native vs. Baysor (prior)": "de_disagree_10x_baysor_prior.csv",
    }

    fig, axes = plt.subplots(1, 5, figsize=(32, 6), sharey=True)
    for ax, (label, fname) in zip(axes, files.items()):
        df = pd.read_csv(TABLES_DIR / fname)
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
        ax.set_title(label)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("-log10(adj. p)")
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
    if all((TABLES_DIR / f).exists() for f in [
        "local_morans_10x_cellpose.csv", "local_morans_10x_stardist.csv",
        "local_morans_10x_voronoi.csv", "local_morans_10x_baysor.csv",
        "local_morans_10x_baysor_prior.csv",
    ]):
        fig_local_morans_map()
    if all((TABLES_DIR / f).exists() for f in [
        "de_disagree_10x_cellpose.csv", "de_disagree_10x_stardist.csv",
        "de_disagree_10x_voronoi.csv", "de_disagree_10x_baysor.csv",
        "de_disagree_10x_baysor_prior.csv",
    ]):
        fig_de_volcano()
    print(f"wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
