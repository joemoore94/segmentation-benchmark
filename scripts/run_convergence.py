"""Method convergence: internal vs. reference-space clustering.

For the same set of expansion methods, computes pairwise ARI in two regimes:
  1. "Internal": each method clustered in its own PCA space
  2. "Reference": each method projected into scRNA-seq reference PCA, then
     clustered in that shared space

The delta (reference - internal) shows whether the reference projection makes
methods agree more (positive = convergence) or less (negative = divergence).

Reads:  data/reference/scrna_3p_filtered_feature_bc_matrix.h5
        data/processed/roi/adata_*.h5ad
Writes: results/figures/convergence_delta_ari.png
        results/tables/convergence_ari.csv

Usage::

    conda run -n segbench python scripts/run_convergence.py
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from sklearn.decomposition import PCA

from segbench.compare import (
    cell_type_agreement,
    cluster_cell_types,
    match_cells_by_centroid,
)
from segbench.constants import METHOD_FAMILIES, METHOD_LABELS
from segbench.style import apply_style

ROI_DIR  = Path("data/processed/roi")
REF_PATH = Path("data/reference/scrna_3p_filtered_feature_bc_matrix.h5")
TABLES   = Path("results/tables")
FIGURES  = Path("results/figures")

MAX_MATCH_DIST = 15.0
LEIDEN_RES = 1.0
N_PCS = 30
RANDOM_STATE = 0

METHODS = [
    ("10x_native",                  "adata_10x.h5ad"),
    ("voronoi",                     "adata_voronoi.h5ad"),
    ("voronoi_stardist",            "adata_voronoi_stardist.h5ad"),
    ("voronoi_mesmer",              "adata_voronoi_mesmer.h5ad"),
    ("voronoi_10x_ranger",          "adata_voronoi_10x_ranger.h5ad"),
    ("baysor",                      "adata_baysor.h5ad"),
    ("baysor_prior_c10",            "adata_baysor_prior_c10.h5ad"),
    ("baysor_stardist_prior_c10",   "adata_baysor_stardist_prior_c10.h5ad"),
    ("baysor_mesmer_prior_c10",     "adata_baysor_mesmer_prior_c10.h5ad"),
    ("baysor_10x_ranger_prior_c10", "adata_baysor_10x_ranger_prior_c10.h5ad"),
]

FAMILY_COLORS = {
    "Reference":          "#55A868",
    "Voronoi":            "#17BECF",
    "Transcript-density": "#DD8452",
}


def normalize_log(X) -> np.ndarray:
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    X = X / totals * np.median(totals)
    return np.log1p(X)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()
    sc.settings.verbosity = 0

    # ---- load adatas
    print("Loading AnnData files...")
    adatas: dict[str, ad.AnnData] = {}
    available: list[tuple[str, str]] = []
    for key, fname in METHODS:
        path = ROI_DIR / fname
        if not path.exists():
            print(f"  {METHOD_LABELS[key]}: skipped")
            continue
        label = METHOD_LABELS[key]
        adatas[label] = ad.read_h5ad(path)
        available.append((key, label))
        print(f"  {label}: {adatas[label].n_obs} cells")

    method_order = [label for _, label in available]
    n = len(method_order)

    # ---- internal clustering (each method in its own PCA)
    print("\nClustering in internal PCA space...")
    internal_labels: dict[str, pd.Series] = {}
    for _, label in available:
        internal_labels[label] = cluster_cell_types(
            adatas[label], resolution=LEIDEN_RES, seed=RANDOM_STATE,
        )
        internal_labels[label].index = internal_labels[label].index.astype(str)
        print(f"  {label}: {internal_labels[label].nunique()} clusters")

    # ---- reference PCA + clustering
    print("\nLoading scRNA-seq reference...")
    ref = sc.read_10x_h5(str(REF_PATH))
    ref.var_names_make_unique()
    sc.pp.filter_cells(ref, min_genes=200)
    sc.pp.filter_genes(ref, min_cells=3)
    ref = ref[ref.obs["n_genes"] < 6000].copy()

    xenium_genes = sorted(set(adatas[method_order[0]].var_names))
    shared_genes = sorted(set(xenium_genes) & set(ref.var_names))
    print(f"  Shared genes: {len(shared_genes)}")

    X_ref = normalize_log(ref[:, shared_genes].X)
    pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE)
    pca.fit(X_ref)
    print(f"  {N_PCS} PCs, {pca.explained_variance_ratio_.sum():.1%} variance")

    print("\nClustering in reference PCA space...")
    ref_labels: dict[str, pd.Series] = {}
    for _, label in available:
        X_norm = normalize_log(adatas[label][:, shared_genes].X)
        Z = pca.transform(X_norm)
        tmp = ad.AnnData(
            X=Z,
            obs=pd.DataFrame(index=adatas[label].obs_names.astype(str)),
        )
        sc.pp.neighbors(tmp, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)
        sc.tl.leiden(tmp, resolution=LEIDEN_RES, random_state=RANDOM_STATE, flavor="igraph")
        ref_labels[label] = tmp.obs["leiden"]
        print(f"  {label}: {ref_labels[label].nunique()} clusters")

    # ---- pairwise ARI: match once per pair, score in both spaces
    print("\nComputing pairwise ARI...")
    ari_internal = np.full((n, n), np.nan)
    ari_ref = np.full((n, n), np.nan)
    np.fill_diagonal(ari_internal, 1.0)
    np.fill_diagonal(ari_ref, 1.0)

    rows = []
    for i, j in combinations(range(n), 2):
        la, lb = method_order[i], method_order[j]
        matches = match_cells_by_centroid(adatas[la], adatas[lb], max_dist=MAX_MATCH_DIST)
        matches["id_a"] = matches["id_a"].astype(str)
        matches["id_b"] = matches["id_b"].astype(str)

        res_int = cell_type_agreement(internal_labels[la], internal_labels[lb], matches)
        res_ref = cell_type_agreement(ref_labels[la], ref_labels[lb], matches)

        ari_internal[i, j] = ari_internal[j, i] = res_int["ari"]
        ari_ref[i, j] = ari_ref[j, i] = res_ref["ari"]

        print(f"  {la} vs {lb}: internal={res_int['ari']:.3f}  "
              f"ref={res_ref['ari']:.3f}  Δ={res_ref['ari'] - res_int['ari']:+.3f}  "
              f"({res_int['n_matched']} matched)")

        rows.append({
            "method_a": la, "method_b": lb,
            "ari_internal": round(res_int["ari"], 4),
            "ari_reference": round(res_ref["ari"], 4),
            "delta": round(res_ref["ari"] - res_int["ari"], 4),
            "n_matched": res_int["n_matched"],
        })

    pd.DataFrame(rows).to_csv(TABLES / "convergence_ari.csv", index=False)

    # ---- summary stats
    delta = ari_ref - ari_internal
    mask = ~np.eye(n, dtype=bool)
    mean_int = np.nanmean(ari_internal[mask])
    mean_ref = np.nanmean(ari_ref[mask])
    mean_delta = np.nanmean(delta[mask])
    print(f"\nMean off-diagonal ARI — internal: {mean_int:.3f}, "
          f"reference: {mean_ref:.3f}, Δ: {mean_delta:+.3f}")

    # Per-family breakdown
    key_for_label = {METHOD_LABELS[k]: k for k, _ in available}
    family_for_label = {l: METHOD_FAMILIES[key_for_label[l]] for l in method_order}

    for tag, test_fn in [
        ("within Voronoi",
         lambda a, b: family_for_label[a] == "Voronoi" and family_for_label[b] == "Voronoi"),
        ("within Baysor",
         lambda a, b: family_for_label[a] == "Transcript-density" and family_for_label[b] == "Transcript-density"),
        ("cross-family",
         lambda a, b: family_for_label[a] != family_for_label[b]),
    ]:
        vals = [delta[i, j] for i, j in combinations(range(n), 2)
                if test_fn(method_order[i], method_order[j])]
        if vals:
            print(f"  {tag}: mean Δ = {np.mean(vals):+.3f} ({len(vals)} pairs)")

    # ---- figure
    fig, axes = plt.subplots(1, 3, figsize=(42, 13))

    off_vals = np.concatenate([
        ari_internal[mask & ~np.isnan(ari_internal)],
        ari_ref[mask & ~np.isnan(ari_ref)],
    ])
    ari_vmin = float(np.nanmin(off_vals))

    panels = [
        (ari_internal,
         f"Internal PCA\n(mean off-diag = {mean_int:.3f})",
         "YlOrRd", ari_vmin, 1.0, ".2f", "ARI"),
        (ari_ref,
         f"Reference PCA\n(mean off-diag = {mean_ref:.3f})",
         "YlOrRd", ari_vmin, 1.0, ".2f", "ARI"),
        (delta,
         f"Δ ARI (reference − internal)\n(mean = {mean_delta:+.3f})",
         "RdBu_r", None, None, "+.2f", "ΔARI"),
    ]

    for ax_idx, (mat, title, cmap, vmin, vmax, fmt, cbar_label) in enumerate(panels):
        ax = axes[ax_idx]
        plot_mat = mat.copy()
        diag_mask = np.eye(n, dtype=bool)
        np.fill_diagonal(plot_mat, np.nan)

        if vmin is None:
            abs_max = np.nanmax(np.abs(plot_mat))
            vmin, vmax = -abs_max, abs_max

        df_plot = pd.DataFrame(plot_mat, index=method_order, columns=method_order)
        sns.heatmap(
            df_plot, annot=True, fmt=fmt, cmap=cmap,
            vmin=vmin, vmax=vmax,
            linewidths=0.5, linecolor="white",
            mask=diag_mask, ax=ax,
            cbar_kws={"label": cbar_label, "shrink": 0.8},
            annot_kws={"size": 13, "weight": "bold"},
        )
        ax.set_title(title, fontweight="bold")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

        for tick in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            text = tick.get_text()
            family = family_for_label.get(text)
            if family and family in FAMILY_COLORS:
                tick.set_color(FAMILY_COLORS[family])
                tick.set_fontweight("bold")

    fig.suptitle(
        "Method convergence: does reference projection increase inter-method agreement?\n"
        f"(Leiden res = {LEIDEN_RES}, mutual nearest-centroid ≤ {MAX_MATCH_DIST:.0f} µm)",
        fontweight="bold", fontsize=22,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "convergence_delta_ari.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved convergence_delta_ari.png")


if __name__ == "__main__":
    main()
