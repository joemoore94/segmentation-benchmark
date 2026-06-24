"""Generate publication-quality figures from the cross-method comparison tables.

Reads ``results/tables/*`` (produced by ``run_comparison.py``) and the per-method
AnnData in ``data/processed/roi/`` and writes PNGs to ``results/figures/``.

Usage::

    conda run -n segbench python scripts/make_figures.py
"""

from __future__ import annotations

import math
from pathlib import Path

import anndata as ad
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment

from segbench.constants import (
    CLUSTER_ANNOTATIONS,
    COMPARISON_ORDER,
    MAIN_METHODS,
    METHOD_COLORS,
    METHOD_LABELS,
)
from segbench.io import PIXEL_SIZE
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

TOTAL_TRANSCRIPTS_FULL_ROI = 3_392_051

DPI = 200
PANEL_W = 9
PANEL_H = 7

DENSITY_CSV_KEY = {
    "voronoi":          "10x native vs. Voronoi",
    "voronoi_stardist": "10x native vs. Voronoi (StarDist)",
    "voronoi_mesmer":   "10x native vs. Voronoi (Mesmer)",
    "baysor":           "10x native vs. Baysor",
    "baysor_prior_c08": "10x native vs. Baysor (prior 0.8)",
    "baysor_prior_c10": "10x native vs. Baysor (prior 1.0)",
}

apply_style()


def _grid_dims(n: int) -> tuple[int, int]:
    ncols = min(n, 3)
    nrows = math.ceil(n / ncols)
    return nrows, ncols


def _available_comparisons() -> list[tuple[str, str]]:
    return [
        (m, label) for m, label in COMPARISON_ORDER
        if (TABLES_DIR / f"disagreement_table_10x_{m}.csv").exists()
    ]


def _make_grid(n: int, sharex: bool = False, sharey: bool = False):
    nrows, ncols = _grid_dims(n)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(PANEL_W * ncols, PANEL_H * nrows),
        sharex=sharex, sharey=sharey,
    )
    flat = np.array(axes).flatten()
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig, flat[:n], nrows, ncols


def fig_cell_counts_and_sizes() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts.csv", index_col="method")
    methods = [m for m in MAIN_METHODS if m in counts.index]

    FILE_KEY = {"10x_native": "10x"}
    adatas_for_violin = {}
    for m in methods:
        fname = f"adata_{FILE_KEY.get(m, m)}.h5ad"
        adatas_for_violin[m] = ad.read_h5ad(ROI_DIR / fname)

    transcripts_by_method = {
        m: np.asarray(adatas_for_violin[m].X.sum(axis=1)).ravel()
        for m in methods
    }

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adata_mesmer = ad.read_h5ad(ROI_DIR / "adata_mesmer.h5ad")
    adata_10x = adatas_for_violin.get("10x_native") or ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")

    fig, axes = plt.subplots(1, 3, figsize=(28, 8))

    axes[0].bar(
        [METHOD_LABELS[m] for m in methods],
        counts.loc[methods, "n_cells"].to_numpy(),
        color=[METHOD_COLORS[m] for m in methods],
    )
    axes[0].set_ylabel("Cell count")
    axes[0].set_title("Cell count (full 2mm × 2mm ROI)", fontweight="bold")
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
        cut=0, inner="box", legend=False, ax=axes[1],
    )
    axes[1].set_xticks(range(len(methods)))
    axes[1].set_xticklabels([METHOD_LABELS[m] for m in methods], rotation=40, ha="right")
    axes[1].set_xlabel("")
    tick_vals = [1, 10, 100, 1000]
    axes[1].set_yticks([np.log10(v) for v in tick_vals])
    axes[1].set_yticklabels([str(v) for v in tick_vals])
    axes[1].set_ylabel("Transcripts per cell")
    axes[1].set_title("Transcripts/cell distribution", fontweight="bold")

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
    axes[2].set_title("Nuclear mask size", fontweight="bold")
    axes[2].legend()

    fig.suptitle("Cell count and QC", fontweight="bold")
    capture = counts.loc[methods, "transcript_capture_rate"]
    capture_str = ", ".join(f"{METHOD_LABELS[m]} {capture[m]:.0%}" for m in methods)
    fig.text(0.5, 0.01, f"Transcript capture: {capture_str}",
             ha="center", style="italic")
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(FIGURES_DIR / "cell_counts_and_sizes.png", dpi=DPI)
    plt.close(fig)


def fig_expression_correlation() -> None:
    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"expression_correlation_10x_{m}.csv"))
        for m, label in _available_comparisons()
    ]

    fig, flat, nrows, ncols = _make_grid(len(pairs), sharey=True)
    for ax, (m, label, corr) in zip(flat, pairs):
        median = corr["correlation"].median()
        sns.histplot(corr["correlation"].dropna(), bins=40, ax=ax, color=METHOD_COLORS[m])
        ax.axvline(median, color="black", linestyle="--")
        ax.text(0.04, 0.94, f"median = {median:.3f}", transform=ax.transAxes,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        ax.set_xlabel("Pearson correlation")
        ax.set_title(f"10x native vs. {label}", fontweight="bold")

    flat[0].set_ylabel("Number of pairs")
    fig.suptitle("Per-cell expression agreement vs. 10x native", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES_DIR / "expression_correlation.png", dpi=DPI)
    plt.close(fig)


def fig_disagreement_spatial_map() -> None:
    pairs = [
        (label, pd.read_csv(TABLES_DIR / f"disagreement_table_10x_{m}.csv"))
        for m, label in _available_comparisons()
    ]

    fig, flat, _, _ = _make_grid(len(pairs))
    for ax, (label, df) in zip(flat, pairs):
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
        ax.set_title(f"10x native vs. {label}", fontweight="bold")

    fig.suptitle("Cell-type agreement vs. disagreement", fontweight="bold")
    fig.legend(handles=[
        mpatches.Patch(color="#4C72B0", label="Agree"),
        mpatches.Patch(color="#C44E52", label="Disagree"),
    ], loc="lower center", ncols=2, framealpha=0.9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(FIGURES_DIR / "disagreement_spatial_map.png", dpi=DPI)
    plt.close(fig)


def fig_density_vs_disagreement() -> None:
    log_density = pd.read_csv(TABLES_DIR / "10x_log_density.csv", index_col=0)["log_density"]
    summary = pd.read_csv(TABLES_DIR / "density_disagreement_summary.csv", index_col="comparison")

    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"disagreement_table_10x_{m}.csv"))
        for m, label in _available_comparisons()
    ]

    fig, flat, nrows, ncols = _make_grid(len(pairs), sharex=True, sharey=True)
    for ax, (m, label, df) in zip(flat, pairs):
        df = df.copy()
        df["log_density"] = df["id_a"].map(log_density)
        sns.kdeplot(
            data=df, x="log_density", hue="disagree",
            palette={0.0: "#4C72B0", 1.0: "#C44E52"}, common_norm=False,
            fill=True, alpha=0.3, ax=ax, legend=False,
        )
        medians = df.groupby("disagree")["log_density"].median()
        ax.axvline(medians[0.0], color="#4C72B0", linestyle="--")
        ax.axvline(medians[1.0], color="#C44E52", linestyle="--")
        csv_key = DENSITY_CSV_KEY.get(m)
        p = summary.loc[csv_key, "p_value"] if csv_key and csv_key in summary.index else float("nan")
        ax.set_title(f"10x native vs. {label}", fontweight="bold")
        ax.text(0.04, 0.94, f"p = {p:.1e}", transform=ax.transAxes,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        ax.set_xlabel("Phenotypic log-density (Mellon)")

    flat[0].set_ylabel("Density")
    fig.suptitle("Phenotypic density vs. cell-type disagreement", fontweight="bold")
    fig.legend(handles=[
        mpatches.Patch(color="#4C72B0", alpha=0.5, label="Agree"),
        mpatches.Patch(color="#C44E52", alpha=0.5, label="Disagree"),
    ], loc="lower center", ncols=2, framealpha=0.9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(FIGURES_DIR / "density_vs_disagreement.png", dpi=DPI)
    plt.close(fig)


def fig_pca_umap() -> None:
    methods = [m for m in MAIN_METHODS if m != "10x_native"] + ["10x_native"]
    methods = [m for m in methods if (TABLES_DIR / f"embedding_{m}.csv").exists()]
    embeddings = {m: pd.read_csv(TABLES_DIR / f"embedding_{m}.csv", index_col=0) for m in methods}

    fig, flat, _, _ = _make_grid(len(methods))
    for ax, method in zip(flat, methods):
        emb = embeddings[method]
        n_clusters = emb["leiden"].nunique()
        palette = sns.color_palette("tab20", n_clusters)
        sns.scatterplot(
            data=emb, x="UMAP1", y="UMAP2", hue="leiden", palette=palette,
            s=12, alpha=0.6, ax=ax, legend=False,
        )
        ax.set_title(f"{METHOD_LABELS[method]} ({n_clusters} clusters)", fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

    fig.suptitle("Per-method Leiden clustering (UMAP)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES_DIR / "pca_umap_clusters.png", dpi=DPI)
    plt.close(fig)


def fig_local_morans_map() -> None:
    pairs = [
        (label, pd.read_csv(TABLES_DIR / f"local_morans_10x_{m}.csv"))
        for m, label in _available_comparisons()
        if (TABLES_DIR / f"local_morans_10x_{m}.csv").exists()
    ]
    LISA_COLORS = {"HH": "#C44E52", "LL": "#4C72B0", "HL": "#DD8452", "LH": "#CCB974"}

    fig, flat, _, _ = _make_grid(len(pairs))
    for ax, (label, df) in zip(flat, pairs):
        for cluster, color in LISA_COLORS.items():
            sub = df[df["lisa_cluster"] == cluster]
            ax.scatter(sub["centroid_x"], sub["centroid_y"], c=color, s=4, alpha=0.6, label=cluster)
        ax.set_title(f"10x native vs. {label}", fontweight="bold")
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        ax.set_aspect("equal")
        ax.invert_yaxis()

    fig.suptitle("Local Moran's I: HH = disagreement hotspot, LL = agreement coldspot",
                 fontweight="bold")
    fig.legend(handles=[mpatches.Patch(color=color, label=cluster)
                        for cluster, color in LISA_COLORS.items()],
               title="LISA cluster", loc="lower center", ncols=4, framealpha=0.9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(FIGURES_DIR / "local_morans_map.png", dpi=DPI)
    plt.close(fig)


def fig_de_volcano() -> None:
    pairs = [
        (m, label, pd.read_csv(TABLES_DIR / f"de_disagree_10x_{m}.csv"))
        for m, label in _available_comparisons()
        if (TABLES_DIR / f"de_disagree_10x_{m}.csv").exists()
    ]

    fig, flat, _, _ = _make_grid(len(pairs), sharey=True)
    for ax, (m, label, df) in zip(flat, pairs):
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
                ha="left",
            )
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel("log fold change (disagree vs. agree)")
        ax.set_title(f"10x native vs. {label}", fontweight="bold")

    flat[0].set_ylabel("-log10(adj. p)")
    fig.suptitle("DE: disagree vs. agree cells (Wilcoxon, 10x native cells)", fontweight="bold")
    fig.legend(handles=[
        mpatches.Patch(color="#AAAAAA", label="n.s."),
        mpatches.Patch(color="#C44E52", label="adj. p < 0.05"),
    ], loc="lower center", ncols=2, framealpha=0.9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(FIGURES_DIR / "de_volcano.png", dpi=DPI)
    plt.close(fig)


def fig_cluster_confusion() -> None:
    from matplotlib.patches import Rectangle

    ct_short = {
        "Luminal epithelial": "Lum. epi.",
        "Macrophages": "Macro.",
        "T cells": "T cells",
        "B cells": "B cells",
        "Myoepithelial": "Myoepi.",
        "CAFs": "CAFs",
        "Smooth muscle": "Sm. muscle",
        "Endothelial": "Endoth.",
        "Plasma cells": "Plasma",
        "Adipocytes": "Adipo.",
    }

    avail = [
        (m, label) for m, label in _available_comparisons()
        if (TABLES_DIR / f"cell_type_confusion_10x_{m}.csv").exists()
    ]
    fig, flat, nrows, ncols = _make_grid(len(avail))

    for ax, (method, label) in zip(flat, avail):
        path = TABLES_DIR / f"cell_type_confusion_10x_{method}.csv"
        raw = pd.read_csv(path, index_col="label_a")
        raw.index = raw.index.astype(str)
        raw.columns = raw.columns.astype(str)

        ref_ids = sorted(raw.index, key=int)
        comp_ids = sorted(raw.columns, key=int)
        raw = raw.loc[ref_ids, comp_ids]

        row_ind, col_ind = linear_sum_assignment(-raw.to_numpy())
        matched = set(zip(row_ind.tolist(), col_ind.tolist()))

        row_sums = raw.sum(axis=1).replace(0, np.nan)
        norm = raw.div(row_sums, axis=0).fillna(0) * 100

        annot_text = norm.map(lambda v: f"{v:.0f}" if v >= 5 else "")

        row_labels = [f"{c}: {ct_short.get(CLUSTER_ANNOTATIONS.get(c, ''), c)}"
                      for c in ref_ids]

        sns.heatmap(
            norm.values, ax=ax,
            cmap="Blues", vmin=0, vmax=100,
            annot=np.array(annot_text), fmt="",
            annot_kws={"weight": "bold"},
            xticklabels=comp_ids, yticklabels=row_labels,
            linewidths=0.4, linecolor="#e0e0e0",
            cbar=False,
        )

        for r, c in matched:
            ax.add_patch(Rectangle((c, r), 1, 1, fill=False,
                                   edgecolor="red", linewidth=2.5))

        ax.set_title(f"{label}  ({len(comp_ids)} clusters)", fontweight="bold")
        ax.set_xlabel(f"{label} cluster")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=0)
        ax.tick_params(axis="y", rotation=0)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    sm = ScalarMappable(cmap="Blues", norm=Normalize(vmin=0, vmax=100))
    sm.set_array([])
    cbar_ax = fig.add_axes([0.92, 0.08, 0.015, 0.84])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("% of 10x native cluster")

    fig.suptitle(
        "Cluster-level confusion matrices (row-normalised)  ·  "
        "Red border = Hungarian-matched pair",
        fontstyle="italic", fontweight="bold",
    )
    fig.subplots_adjust(left=0.08, right=0.88, top=0.93, bottom=0.05,
                        hspace=0.4, wspace=0.4)
    fig.savefig(FIGURES_DIR / "confusion_clusters.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_cell_counts_and_sizes()
    fig_expression_correlation()
    fig_disagreement_spatial_map()
    fig_density_vs_disagreement()
    fig_pca_umap()
    lisa_files = [f"local_morans_10x_{m}.csv" for m, _ in _available_comparisons()]
    if any((TABLES_DIR / f).exists() for f in lisa_files):
        fig_local_morans_map()
    de_files = [f"de_disagree_10x_{m}.csv" for m, _ in _available_comparisons()]
    if any((TABLES_DIR / f).exists() for f in de_files):
        fig_de_volcano()
    confusion_csvs = [f"cell_type_confusion_10x_{m}.csv" for m, _ in _available_comparisons()]
    if any((TABLES_DIR / f).exists() for f in confusion_csvs):
        fig_cluster_confusion()
    print(f"wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
