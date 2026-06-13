"""Generate publication-quality figures from the cross-method comparison tables.

Reads ``results/tables/*`` (produced by ``run_comparison.py``) and the per-method
AnnData in ``data/processed/roi/`` and writes PNGs to ``results/figures/``.

Usage::

    conda run -n segbench python scripts/make_figures.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import numpy as np

from segbench.compare import subset_to_region
from segbench.io import PIXEL_SIZE

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

# Baysor only segmented the centered 1mm x 1mm sub-region of the 2mm x 2mm
# ROI (CPU-tractability, see docs/dataset.md). Subsetting CellPose to the
# same sub-region (see run_comparison.py) gives a direct, area-matched
# cell count/size comparison.
SUB_REGION = ((500.0, 1500.0), (500.0, 1500.0))  # (x_range, y_range), microns

# qv>=20 non-control transcripts in the 1mm x 1mm sub-region (see docs/dataset.md)
TOTAL_TRANSCRIPTS_1MM2 = 770_748

sns.set_theme(style="whitegrid", context="talk")


METHOD_COLORS = {
    "cellpose": "#4C72B0",
    "baysor": "#DD8452",
    "10x_native": "#55A868",
    "stardist": "#8172B2",
}
METHOD_LABELS = {
    "cellpose": "CellPose",
    "baysor": "Baysor",
    "10x_native": "10x native",
    "stardist": "StarDist",
}


def fig_cell_counts_and_sizes() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts_1mm2.csv", index_col="method")
    methods = list(counts.index)

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    x_range, y_range = SUB_REGION
    adata_cellpose_sub = subset_to_region(adata_cellpose, x_range, y_range)
    adata_10x_sub = subset_to_region(adata_10x, x_range, y_range)
    adata_stardist_sub = subset_to_region(adata_stardist, x_range, y_range)

    cellpose_area_um2 = adata_cellpose_sub.obs["area"] * PIXEL_SIZE**2
    tenx_nucleus_area_um2 = adata_10x_sub.obs["nucleus_area_um2"]
    stardist_area_um2 = adata_stardist_sub.obs["area"] * PIXEL_SIZE**2
    cellpose_transcripts = np.asarray(adata_cellpose_sub.X.sum(axis=1)).ravel()
    baysor_transcripts = np.asarray(adata_baysor.X.sum(axis=1)).ravel()
    tenx_transcripts = np.asarray(adata_10x_sub.X.sum(axis=1)).ravel()
    stardist_transcripts = np.asarray(adata_stardist_sub.X.sum(axis=1)).ravel()
    transcripts_by_method = {
        "cellpose": cellpose_transcripts,
        "baysor": baysor_transcripts,
        "10x_native": tenx_transcripts,
        "stardist": stardist_transcripts,
    }

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].bar(
        [METHOD_LABELS[m] for m in methods],
        counts.loc[methods, "n_cells"].to_numpy(),
        color=[METHOD_COLORS[m] for m in methods],
    )
    axes[0].set_ylabel("Cell count (1mm × 1mm sub-region)")
    axes[0].set_title("Cell count")

    # Transcripts/cell is computed identically for all methods (sum of the
    # per-cell gene-count matrix), so this panel is a true apples-to-apples
    # QC comparison -- unlike "cell size", which means different things
    # (nucleus pixel area vs. transcript count) per method.
    for method in methods:
        sns.histplot(
            transcripts_by_method[method], bins=40, log_scale=True, ax=axes[1],
            color=METHOD_COLORS[method], label=METHOD_LABELS[method], alpha=0.4,
        )
    axes[1].set_xlabel("Transcripts per cell")
    axes[1].set_title("Transcripts/cell (QC)")
    axes[1].legend()

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
        "Cell count and QC: CellPose vs. Baysor vs. 10x native vs. StarDist (1mm × 1mm sub-region)"
    )
    capture = counts.loc[methods, "transcript_capture_rate"]
    capture_str = ", ".join(
        f"{METHOD_LABELS[m]} {capture[m]:.0%}" for m in methods
    )
    fig.text(
        0.5, 0.01,
        f"Transcript capture rate (assigned / {TOTAL_TRANSCRIPTS_1MM2:,} total"
        f" qv≥20 transcripts in region): {capture_str}. CellPose and StarDist are"
        " nuclear-only, so cytoplasmic transcripts are not assigned to any cell.",
        ha="center", fontsize=11, style="italic",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIGURES_DIR / "cell_counts_and_sizes.png", dpi=150)
    plt.close(fig)


def fig_expression_correlation() -> None:
    corr_baysor = pd.read_csv(TABLES_DIR / "expression_correlation.csv")
    corr_10x = pd.read_csv(TABLES_DIR / "expression_correlation_cellpose_10x.csv")
    corr_stardist = pd.read_csv(TABLES_DIR / "expression_correlation_cellpose_stardist.csv")

    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5), sharey=True)
    for ax, corr, label, color in [
        (axes[0], corr_baysor, "CellPose vs. Baysor", METHOD_COLORS["baysor"]),
        (axes[1], corr_10x, "CellPose vs. 10x native", METHOD_COLORS["10x_native"]),
        (axes[2], corr_stardist, "CellPose vs. StarDist", METHOD_COLORS["stardist"]),
    ]:
        median = corr["correlation"].median()
        sns.histplot(corr["correlation"].dropna(), bins=40, ax=ax, color=color)
        ax.axvline(median, color="black", linestyle="--", label=f"median = {median:.2f}")
        ax.set_xlabel("Pearson correlation (matched cell pair)")
        ax.set_title(label)
        ax.legend()

    axes[0].set_ylabel("Number of pairs")
    fig.suptitle("Per-cell expression agreement vs. CellPose")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "expression_correlation.png", dpi=150)
    plt.close(fig)


def fig_disagreement_spatial_map() -> None:
    disagreement_baysor = pd.read_csv(TABLES_DIR / "disagreement_table.csv")
    disagreement_10x = pd.read_csv(TABLES_DIR / "disagreement_table_cellpose_10x.csv")
    disagreement_stardist = pd.read_csv(TABLES_DIR / "disagreement_table_cellpose_stardist.csv")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, disagreement, label in [
        (axes[0], disagreement_baysor, "CellPose vs. Baysor"),
        (axes[1], disagreement_10x, "CellPose vs. 10x native"),
        (axes[2], disagreement_stardist, "CellPose vs. StarDist"),
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
        )
        ax.set_xlabel("x (µm, ROI coordinates)")
        ax.set_ylabel("y (µm, ROI coordinates)")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(label)
        ax.legend(title="Disagree", labels=["No", "Yes"])

    fig.suptitle("Cell-type agreement (blue) vs. disagreement (red)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "disagreement_spatial_map.png", dpi=150)
    plt.close(fig)


def fig_cell_type_confusion() -> None:
    confusion_baysor = pd.read_csv(TABLES_DIR / "cell_type_confusion.csv", index_col=0)
    confusion_10x = pd.read_csv(TABLES_DIR / "cell_type_confusion_cellpose_10x.csv", index_col=0)
    confusion_stardist = pd.read_csv(
        TABLES_DIR / "cell_type_confusion_cellpose_stardist.csv", index_col=0
    )

    fig, axes = plt.subplots(1, 3, figsize=(22, 6))

    sns.heatmap(confusion_baysor, annot=False, cmap="viridis", ax=axes[0])
    axes[0].set_xlabel("Baysor Leiden cluster")
    axes[0].set_ylabel("CellPose Leiden cluster")
    axes[0].set_title("CellPose vs. Baysor")

    sns.heatmap(confusion_10x, annot=False, cmap="viridis", ax=axes[1])
    axes[1].set_xlabel("10x native Leiden cluster")
    axes[1].set_ylabel("CellPose Leiden cluster")
    axes[1].set_title("CellPose vs. 10x native")

    sns.heatmap(confusion_stardist, annot=False, cmap="viridis", ax=axes[2])
    axes[2].set_xlabel("StarDist Leiden cluster")
    axes[2].set_ylabel("CellPose Leiden cluster")
    axes[2].set_title("CellPose vs. StarDist")

    fig.suptitle("Cell-type cluster correspondence (matched pairs)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "cell_type_confusion.png", dpi=150)
    plt.close(fig)


def fig_density_vs_disagreement() -> None:
    log_density = pd.read_csv(TABLES_DIR / "cellpose_log_density.csv", index_col=0)["log_density"]
    summary = pd.read_csv(TABLES_DIR / "density_disagreement_summary.csv", index_col="comparison")

    disagreement_baysor = pd.read_csv(TABLES_DIR / "disagreement_table.csv")
    disagreement_10x = pd.read_csv(TABLES_DIR / "disagreement_table_cellpose_10x.csv")
    disagreement_stardist = pd.read_csv(TABLES_DIR / "disagreement_table_cellpose_stardist.csv")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)
    for ax, disagreement, label in [
        (axes[0], disagreement_baysor, "CellPose vs. Baysor"),
        (axes[1], disagreement_10x, "CellPose vs. 10x native"),
        (axes[2], disagreement_stardist, "CellPose vs. StarDist"),
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
        ax.set_xlabel("CellPose phenotypic log-density (Mellon)")
        ax.legend(title="Disagree", labels=["Yes", "No"])

    axes[0].set_ylabel("Density")
    fig.suptitle("CellPose phenotypic density (Mellon) vs. cell-type call disagreement")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "density_vs_disagreement.png", dpi=150)
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_cell_counts_and_sizes()
    fig_expression_correlation()
    fig_disagreement_spatial_map()
    fig_cell_type_confusion()
    fig_density_vs_disagreement()
    print(f"wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
