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


def fig_cell_counts_and_sizes() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts_1mm2.csv", index_col="method")

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    x_range, y_range = SUB_REGION
    adata_cellpose_sub = subset_to_region(adata_cellpose, x_range, y_range)
    cellpose_area_um2 = adata_cellpose_sub.obs["area"] * PIXEL_SIZE**2
    cellpose_transcripts = np.asarray(adata_cellpose_sub.X.sum(axis=1)).ravel()
    baysor_transcripts = np.asarray(adata_baysor.X.sum(axis=1)).ravel()

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].bar(
        counts.index, counts["n_cells"].to_numpy(), color=["#4C72B0", "#DD8452"]
    )
    axes[0].set_ylabel("Cell count (1mm × 1mm sub-region)")
    axes[0].set_title("Cell count")

    # Transcripts/cell is computed identically for both methods (sum of the
    # per-cell gene-count matrix), so this panel is a true apples-to-apples
    # QC comparison -- unlike "cell size", which means different things
    # (nucleus pixel area vs. transcript count) per method.
    sns.histplot(
        cellpose_transcripts, bins=40, log_scale=True, ax=axes[1],
        color="#4C72B0", label="CellPose", alpha=0.5,
    )
    sns.histplot(
        baysor_transcripts, bins=40, log_scale=True, ax=axes[1],
        color="#DD8452", label="Baysor", alpha=0.5,
    )
    axes[1].set_xlabel("Transcripts per cell")
    axes[1].set_title("Transcripts/cell (QC)")
    axes[1].legend()

    sns.histplot(cellpose_area_um2, bins=50, ax=axes[2], color="#4C72B0")
    axes[2].set_xlabel("Nucleus area (µm²)")
    axes[2].set_title("CellPose nucleus area (QC)")

    fig.suptitle("Cell count and QC: CellPose vs. Baysor (1mm × 1mm sub-region)")
    capture_cp = counts.loc["cellpose", "transcript_capture_rate"]
    capture_bs = counts.loc["baysor", "transcript_capture_rate"]
    fig.text(
        0.5, 0.01,
        f"Transcript capture rate (assigned / {TOTAL_TRANSCRIPTS_1MM2:,} total"
        f" qv≥20 transcripts in region): CellPose {capture_cp:.0%},"
        f" Baysor {capture_bs:.0%}. CellPose is nuclear-only, so cytoplasmic"
        " transcripts are not assigned to any cell.",
        ha="center", fontsize=11, style="italic",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIGURES_DIR / "cell_counts_and_sizes.png", dpi=150)
    plt.close(fig)


def fig_expression_correlation() -> None:
    corr = pd.read_csv(TABLES_DIR / "expression_correlation.csv")

    fig, ax = plt.subplots(figsize=(8, 5.5))
    sns.histplot(corr["correlation"].dropna(), bins=40, ax=ax, color="#55A868")
    ax.axvline(
        corr["correlation"].median(), color="black", linestyle="--",
        label=f"median = {corr['correlation'].median():.2f}",
    )
    ax.set_xlabel("Pearson correlation (matched cell pair)")
    ax.set_ylabel("Number of pairs")
    ax.set_title("Per-cell expression agreement: CellPose vs. Baysor")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "expression_correlation.png", dpi=150)
    plt.close(fig)


def fig_disagreement_spatial_map() -> None:
    disagreement = pd.read_csv(TABLES_DIR / "disagreement_table.csv")

    fig, ax = plt.subplots(figsize=(6, 6))
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
    ax.set_title("Cell-type agreement (blue) vs. disagreement (red)")
    ax.legend(title="Disagree", labels=["No", "Yes"])
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "disagreement_spatial_map.png", dpi=150)
    plt.close(fig)


def fig_cell_type_confusion() -> None:
    confusion = pd.read_csv(TABLES_DIR / "cell_type_confusion.csv", index_col=0)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(confusion, annot=False, cmap="viridis", ax=ax)
    ax.set_xlabel("Baysor Leiden cluster")
    ax.set_ylabel("CellPose Leiden cluster")
    ax.set_title("Cell-type cluster correspondence (matched pairs)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "cell_type_confusion.png", dpi=150)
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_cell_counts_and_sizes()
    fig_expression_correlation()
    fig_disagreement_spatial_map()
    fig_cell_type_confusion()
    print(f"wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
