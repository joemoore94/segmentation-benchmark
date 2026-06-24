"""Marker gene recovery by segmentation method and cell type.

For cells matched between 10x native and each comparison method, uses the
10x-native cell-type annotation as ground truth and asks: how much of each
cell type's canonical marker signal does each method recover?

For a given cell type T and marker gene M, the recovery score for method X is:
  mean(log-norm expression of M in X's matched cells of type T)
relative to the 10x-native baseline.

Left panel: heatmap of relative recovery (method / 10x native) per
            (cell type × marker gene) combination. Values near 1 = full
            recovery; <1 = marker expression reduced in that method.

Right panel: absolute mean log-norm expression per cell type for three
             selected markers, one line per method, to show the scale of
             differences.

Reads:  data/processed/roi/adata_*.h5ad
        results/tables/disagreement_table_10x_*.csv
Writes: results/figures/marker_recovery.png

Usage::

    conda run -n segbench python scripts/run_marker_recovery.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import scanpy as sc
from segbench.constants import (
    CLUSTER_ANNOTATIONS, METHOD_COLORS, METHOD_LABELS, NUCLEAR_ONLY,
)
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES  = Path("results/tables")
FIGURES = Path("results/figures")

CELL_TYPE_MARKERS: dict[str, list[str]] = {
    "Luminal epithelial": ["GATA3", "ESR1", "PGR", "MUC1"],
    "Macrophages":        ["CD14", "LYZ", "FCER1G"],
    "T cells":            ["CD3E", "TRAC", "IL7R"],
    "B cells":            ["MS4A1", "BANK1"],
    "Myoepithelial":      ["ACTA2", "KRT14", "MYH11"],
    "CAFs":               ["LUM", "FBLN1", "MMP2"],
    "Endothelial":        ["PECAM1", "VWF"],
}

_ALL_COMPARISONS = [
    "cellpose", "stardist", "mesmer",
    "voronoi", "voronoi_stardist", "voronoi_mesmer",
    "baysor", "baysor_prior_c08", "baysor_prior_c10", "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10", "bidcell", "segger",
]
COMPARISONS = [(k, METHOD_LABELS[k], METHOD_COLORS[k]) for k in _ALL_COMPARISONS]
PLOT_COMPARISONS = [c for c in COMPARISONS if c[0] not in NUCLEAR_ONLY]

# Markers to show in the right-panel line plots (one per cell type family)
HIGHLIGHT_MARKERS = ["GATA3", "CD3E", "LYZ"]
HIGHLIGHT_CELLTYPES = ["Luminal epithelial", "T cells", "Macrophages"]


def log_norm(adata: ad.AnnData) -> np.ndarray:
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    X = X / totals * 1e4
    return np.log1p(X)


def mean_marker_expression(
    adata: ad.AnnData, X_lognorm: np.ndarray, cell_ids: list[str],
    markers: list[str]
) -> dict[str, float]:
    id_to_idx = {name: i for i, name in enumerate(adata.obs_names)}
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}
    idxs = [id_to_idx[cid] for cid in cell_ids if cid in id_to_idx]
    result = {}
    for marker in markers:
        gene_idx = gene_to_idx.get(marker)
        if gene_idx is None:
            result[marker] = np.nan
            continue
        vals = X_lognorm[idxs, gene_idx] if idxs else np.array([])
        result[marker] = float(np.mean(vals)) if len(vals) > 0 else np.nan
    return result


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style()

    print("Loading 10x native and building cell type labels...")
    adata_10x_raw = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adata_10x = adata_10x_raw.copy()
    sc.pp.normalize_total(adata_10x)
    sc.pp.log1p(adata_10x)
    sc.pp.pca(adata_10x, n_comps=30, random_state=0)
    sc.pp.neighbors(adata_10x, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata_10x, resolution=1.0, random_state=0, flavor="igraph")
    cell_type_labels = adata_10x.obs["leiden"].map(CLUSTER_ANNOTATIONS)
    X_10x_lognorm = log_norm(adata_10x_raw)

    # Collect all markers
    all_markers = [m for markers in CELL_TYPE_MARKERS.values() for m in markers]

    # Build rows: for each method, cell type, marker → mean expression
    rows_ref = []
    for cell_type, markers in CELL_TYPE_MARKERS.items():
        ct_ids = list(cell_type_labels[cell_type_labels == cell_type].index)
        expr = mean_marker_expression(adata_10x_raw, X_10x_lognorm, ct_ids, markers)
        for marker, val in expr.items():
            rows_ref.append({"method": "10x native", "cell_type": cell_type,
                             "marker": marker, "mean_lognorm": val})
    df_ref = pd.DataFrame(rows_ref)

    print("Computing marker recovery per method...")
    rows_comp = []
    for method, label, _ in COMPARISONS:
        path = TABLES / f"disagreement_table_10x_{method}.csv"
        if not path.exists():
            print(f"  Missing {path.name}, skipping")
            continue
        dtable = pd.read_csv(path)
        id_b_from_a = dict(zip(dtable["id_a"].astype(str), dtable["id_b"].astype(str)))

        adata_comp = ad.read_h5ad(ROI_DIR / f"adata_{method}.h5ad")
        X_comp_lognorm = log_norm(adata_comp)

        for cell_type, markers in CELL_TYPE_MARKERS.items():
            ct_ids_10x = list(cell_type_labels[cell_type_labels == cell_type].index)
            # Get matched cells in the comparison method
            comp_ids = [id_b_from_a[cid] for cid in ct_ids_10x if cid in id_b_from_a]
            expr = mean_marker_expression(adata_comp, X_comp_lognorm, comp_ids, markers)
            for marker, val in expr.items():
                rows_comp.append({"method": label, "cell_type": cell_type,
                                  "marker": marker, "mean_lognorm": val})
        print(f"  {label}: done")

    df_comp = pd.DataFrame(rows_comp)
    df_all = pd.concat([df_ref, df_comp], ignore_index=True)

    # Relative recovery: method / 10x native per (cell_type, marker)
    baseline = df_ref.set_index(["cell_type", "marker"])["mean_lognorm"].rename("baseline")
    df_comp2 = df_comp.join(baseline, on=["cell_type", "marker"])
    df_comp2["relative_recovery"] = df_comp2["mean_lognorm"] / df_comp2["baseline"].replace(0, np.nan)

    # Pivot for heatmap: rows = cell_type + marker, cols = method (plot subset only)
    plot_labels = {l for _, l, _ in PLOT_COMPARISONS}
    df_plot = df_comp2[df_comp2["method"].isin(plot_labels)]
    df_plot = df_plot.copy()
    df_plot["ct_marker"] = df_plot["cell_type"] + " · " + df_plot["marker"]
    pivot = df_plot.pivot(index="ct_marker", columns="method", values="relative_recovery")
    row_order = [f"{ct} · {m}" for ct, markers in CELL_TYPE_MARKERS.items() for m in markers
                 if f"{ct} · {m}" in pivot.index]
    col_order = [l for _, l, _ in PLOT_COMPARISONS if l in pivot.columns]
    pivot = pivot.reindex(index=row_order, columns=col_order)

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(26, 13))
    gs_outer = gridspec.GridSpec(1, 2, figure=fig,
                                 width_ratios=[1.4, 1], wspace=0.38)
    ax_left = fig.add_subplot(gs_outer[0])
    gs_right = gridspec.GridSpecFromSubplotSpec(
        3, 1, subplot_spec=gs_outer[1], hspace=0.55)
    axes_right = [fig.add_subplot(gs_right[i]) for i in range(3)]

    # Left: recovery heatmap
    sns.heatmap(
        pivot,
        annot=True, fmt=".2f",
        cmap="RdYlGn", center=1.0, vmin=0.5, vmax=1.4,
        linewidths=0.4, linecolor="white",
        ax=ax_left,
        cbar_kws={"label": "Expression relative to 10x native", "shrink": 0.7},
        annot_kws={"size": 11},
    )
    ax_left.set_title(
        "Marker gene recovery relative to 10x native\n"
        "(matched cells, 10x-native cell type labels)",
        fontweight="bold", fontsize=13)
    ax_left.set_xlabel("")
    ax_left.set_ylabel("")
    ax_left.tick_params(axis="x", rotation=35, labelsize=11)
    ax_left.tick_params(axis="y", labelsize=11)
    plt.setp(ax_left.get_xticklabels(), ha="right")

    # Cell type group dividers
    marker_counts = [len(v) for v in CELL_TYPE_MARKERS.values()]
    cumulative = np.cumsum(marker_counts)
    for c in cumulative[:-1]:
        ax_left.axhline(c, color="black", linewidth=1.5)

    # Right: absolute expression for 3 highlight markers (one panel each)
    methods_all = ["10x native"] + [l for _, l, _ in PLOT_COMPARISONS]
    colors_all  = ["#55A868"] + [c for _, _, c in PLOT_COMPARISONS]
    x = np.arange(len(methods_all))

    for i, (ct, marker) in enumerate(zip(HIGHLIGHT_CELLTYPES, HIGHLIGHT_MARKERS)):
        ax_sub = axes_right[i]
        vals = [df_all[(df_all["method"] == m) & (df_all["cell_type"] == ct) &
                       (df_all["marker"] == marker)]["mean_lognorm"].values
                for m in methods_all]
        vals_flat = [v[0] if len(v) > 0 else np.nan for v in vals]
        ax_sub.bar(x, vals_flat, color=colors_all, edgecolor="white", alpha=0.88)
        ax_sub.set_xticks(x)
        if i == 2:
            ax_sub.set_xticklabels(methods_all, rotation=35, ha="right", fontsize=11)
        else:
            ax_sub.set_xticklabels([], fontsize=11)
        ax_sub.set_title(f"{marker}  ·  {ct}", fontweight="bold", fontsize=12)
        ax_sub.set_ylabel("Mean log-norm", fontsize=11)
        ax_sub.tick_params(axis="y", labelsize=11)

    fig.suptitle(
        "Does segmentation method affect recovery of canonical cell-type marker genes?",
        fontsize=14, fontweight="bold",
    )
    fig.savefig(FIGURES / "marker_recovery.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved marker_recovery.png")


if __name__ == "__main__":
    main()
