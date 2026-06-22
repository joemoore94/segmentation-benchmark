"""Pseudobulk correlation: do methods agree on cell-type-level expression?

Individual cell misclassification does not necessarily imply population-level
expression errors. This script sums counts across all matched cells of each
10x-native cell type for every method, normalises to CPM, and computes
Pearson correlation between each method's per-type pseudobulk and the 10x-
native reference. A high correlation (>0.95) at the pseudobulk level alongside
low ARI at the single-cell level would indicate that methods disagree on how
to assign individual cells to types but still produce coherent population
profiles.

Left panel: heatmap of per-cell-type pseudobulk Pearson correlation,
            method (rows) × 10x-native cell type (columns).
Right panel: global pseudobulk Pearson (all cell types pooled) vs. ARI,
             showing whether the two metrics are correlated or dissociated.

Reads:  data/processed/roi/adata_*.h5ad
        results/tables/disagreement_table_10x_*.csv
Writes: results/figures/pseudobulk_correlation.png
        results/tables/pseudobulk_correlation.csv

Usage::

    conda run -n segbench python scripts/run_pseudobulk_correlation.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import scanpy as sc
from scipy.stats import pearsonr
from segbench.style import apply_style

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

COMPARISONS = [
    ("cellpose",       "CellPose",     "#4C72B0", 0.547),
    ("stardist",       "StarDist",     "#8172B2", 0.545),
    ("mesmer",         "Mesmer",       "#D62728", 0.557),
    ("voronoi",        "Voronoi (CP)", "#17BECF", 0.630),
    ("voronoi_mesmer", "Voronoi (M)",  "#BCBD22", 0.686),
    ("baysor",         "Baysor",       "#DD8452", 0.305),
]

CELL_TYPES = [
    "Luminal epithelial", "Macrophages", "T cells", "B cells",
    "Myoepithelial", "CAFs", "Smooth muscle", "Endothelial",
    "Plasma cells", "Adipocytes",
]


def to_dense(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray(X.todense())
    return np.asarray(X)


def pseudobulk(adata: ad.AnnData, cell_ids: list[str]) -> np.ndarray:
    idx_map = {name: i for i, name in enumerate(adata.obs_names)}
    idxs = [idx_map[c] for c in cell_ids if c in idx_map]
    if not idxs:
        return np.zeros(adata.n_vars)
    X = to_dense(adata.X[idxs, :])
    raw = X.sum(axis=0)
    total = raw.sum()
    return raw / total * 1e6 if total > 0 else raw  # CPM


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style()

    print("Loading 10x native and building cell type labels...")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    sc.pp.normalize_total(adata_10x)
    sc.pp.log1p(adata_10x)
    sc.pp.pca(adata_10x, n_comps=30, random_state=0)
    sc.pp.neighbors(adata_10x, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata_10x, resolution=1.0, random_state=0, flavor="igraph")
    cell_type_labels = adata_10x.obs["leiden"].map(CLUSTER_ANNOTATIONS)

    # Reload raw counts for pseudobulk
    adata_10x_raw = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")

    # 10x-native pseudobulk per cell type
    ref_pseudobulk: dict[str, np.ndarray] = {}
    for ct in CELL_TYPES:
        ct_ids = list(cell_type_labels[cell_type_labels == ct].index)
        ref_pseudobulk[ct] = pseudobulk(adata_10x_raw, ct_ids)

    # Global reference pseudobulk (all cells)
    all_ct_ids = list(cell_type_labels.index)
    ref_global = pseudobulk(adata_10x_raw, all_ct_ids)

    rows = []
    for method, label, color, ari in COMPARISONS:
        path = TABLES / f"disagreement_table_10x_{method}.csv"
        if not path.exists():
            print(f"  Missing {path.name}, skipping")
            continue
        dtable = pd.read_csv(path)
        id_b_from_a = dict(zip(dtable["id_a"].astype(str), dtable["id_b"].astype(str)))

        adata_comp = ad.read_h5ad(ROI_DIR / f"adata_{method}.h5ad")

        # Shared genes for fair comparison
        shared = adata_10x_raw.var_names.intersection(adata_comp.var_names)

        corrs_per_ct: dict[str, float] = {}
        for ct in CELL_TYPES:
            ct_ids_10x = list(cell_type_labels[cell_type_labels == ct].index)
            comp_ids = [id_b_from_a[cid] for cid in ct_ids_10x if cid in id_b_from_a]
            if len(comp_ids) < 5:
                corrs_per_ct[ct] = np.nan
                continue
            ref_pb = ref_pseudobulk[ct][
                [list(adata_10x_raw.var_names).index(g) for g in shared]
            ]
            comp_pb_full = pseudobulk(adata_comp, comp_ids)
            comp_pb = comp_pb_full[
                [list(adata_comp.var_names).index(g) for g in shared]
            ]
            r, _ = pearsonr(np.log1p(ref_pb), np.log1p(comp_pb))
            corrs_per_ct[ct] = round(float(r), 4)

        # Global correlation (matched cells, all types pooled)
        all_comp_ids = [id_b_from_a[cid] for cid in all_ct_ids if cid in id_b_from_a]
        ref_global_shared = ref_global[
            [list(adata_10x_raw.var_names).index(g) for g in shared]
        ]
        comp_global_full = pseudobulk(adata_comp, all_comp_ids)
        comp_global = comp_global_full[
            [list(adata_comp.var_names).index(g) for g in shared]
        ]
        r_global, _ = pearsonr(np.log1p(ref_global_shared), np.log1p(comp_global))
        r_global = round(float(r_global), 4)
        print(f"  {label}: global pseudobulk r={r_global:.4f}, ARI={ari:.3f}")

        row = {"method": label, "ari": ari, "global_r": r_global, "color": color}
        row.update(corrs_per_ct)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "pseudobulk_correlation.csv", index=False)

    # Pivot for heatmap: method × cell type
    pivot = df.set_index("method")[CELL_TYPES]
    pivot = pivot.reindex(index=[l for _, l, _, _ in COMPARISONS if l in pivot.index])

    # ---------------------------------------------------------------- figure
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(26, 9),
                                             gridspec_kw={"width_ratios": [2, 1]})

    sns.heatmap(
        pivot,
        annot=True, fmt=".3f",
        cmap="YlOrRd",
        vmin=0.90, vmax=1.0,
        linewidths=0.4, linecolor="white",
        ax=ax_left,
        cbar_kws={"label": "Pseudobulk Pearson r (log CPM)", "shrink": 0.7},
        annot_kws={"size": 10, "weight": "bold"},
    )
    ax_left.set_title(
        "Per-cell-type pseudobulk Pearson correlation vs. 10x native\n"
        "(matched cells only, log CPM, shared genes)",
        fontweight="bold", fontsize=12,
    )
    ax_left.set_xlabel("")
    ax_left.set_ylabel("")
    ax_left.tick_params(axis="x", rotation=35, labelsize=10)
    plt.setp(ax_left.get_xticklabels(), ha="right")
    ax_left.tick_params(axis="y", rotation=0, labelsize=11)

    # Right: global pseudobulk r vs. ARI scatter
    method_order = [l for _, l, _, _ in COMPARISONS if l in df["method"].values]
    x_vals = df.set_index("method").loc[method_order, "ari"].values
    y_vals = df.set_index("method").loc[method_order, "global_r"].values
    colors = [c for _, l, c, _ in COMPARISONS if l in df["method"].values]

    ax_right.scatter(x_vals, y_vals, c=colors, s=220, zorder=5, edgecolors="white", linewidths=1.5)
    for label, x, y in zip(method_order, x_vals, y_vals):
        ax_right.annotate(label, (x, y), textcoords="offset points", xytext=(6, 4), fontsize=10)
    ax_right.set_xlabel("Single-cell ARI vs. 10x native", fontsize=11)
    ax_right.set_ylabel("Global pseudobulk Pearson r (log CPM)", fontsize=11)
    ax_right.set_title("Population-level vs.\nsingle-cell agreement", fontweight="bold", fontsize=12)
    ax_right.set_xlim(0.25, 0.75)

    fig.suptitle(
        "Population-level expression agreement with 10x native (pseudobulk CPM correlation)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "pseudobulk_correlation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved pseudobulk_correlation.png")


if __name__ == "__main__":
    main()
