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

from segbench.io import PIXEL_SIZE

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

# CellPose ran on the full 2mm x 2mm ROI; Baysor ran on the centered 1mm x 1mm
# sub-region (CPU-tractability). Raw cell counts aren't comparable across the
# two areas, so we report density (cells/mm^2) instead.
CELLPOSE_ROI_AREA_MM2 = 2.0 * 2.0
BAYSOR_ROI_AREA_MM2 = 1.0 * 1.0

sns.set_theme(style="whitegrid", context="talk")


def fig_cell_counts_and_sizes() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts.csv", index_col="method")
    density = pd.Series(
        {
            "cellpose": counts.loc["cellpose", "n_cells"] / CELLPOSE_ROI_AREA_MM2,
            "baysor": counts.loc["baysor", "n_cells"] / BAYSOR_ROI_AREA_MM2,
        }
    )

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_baysor = ad.read_h5ad(ROI_DIR / "adata_baysor.h5ad")
    cellpose_area_um2 = adata_cellpose.obs["area"] * PIXEL_SIZE**2
    baysor_n_transcripts = adata_baysor.obs["n_transcripts"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].bar(density.index, density.to_numpy(), color=["#4C72B0", "#DD8452"])
    axes[0].set_ylabel("Cell density (cells / mm²)")
    axes[0].set_title("Cell density")

    sns.histplot(cellpose_area_um2, bins=50, ax=axes[1], color="#4C72B0")
    axes[1].set_xlabel("Nucleus area (µm²)")
    axes[1].set_title("CellPose: cell size")

    sns.histplot(baysor_n_transcripts, bins=50, ax=axes[2], color="#DD8452")
    axes[2].set_xlabel("Transcripts per cell")
    axes[2].set_title("Baysor: cell size")

    fig.suptitle("Cell count and size distributions by segmentation method")
    fig.text(
        0.5, 0.01,
        "Densities are not directly comparable across panels: CellPose ran on the full"
        " 2mm × 2mm ROI; Baysor ran on the centered 1mm × 1mm sub-region.",
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
