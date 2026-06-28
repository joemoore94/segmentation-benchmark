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


def fig_nuclear_mask_sizes() -> None:
    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adata_mesmer = ad.read_h5ad(ROI_DIR / "adata_mesmer.h5ad")

    nuclear_methods = [
        ("cellpose", adata_cellpose),
        ("stardist", adata_stardist),
        ("mesmer",   adata_mesmer),
    ]
    ranger_path = ROI_DIR / "adata_10x_ranger.h5ad"
    if ranger_path.exists():
        nuclear_methods.append(("10x_ranger", ad.read_h5ad(ranger_path)))

    fig, ax = plt.subplots(figsize=(14, 10))
    for key, adata_nuc in nuclear_methods:
        area_um2 = adata_nuc.obs["area"] * PIXEL_SIZE**2
        sns.histplot(area_um2, bins=50, ax=ax,
                     color=METHOD_COLORS[key], label=METHOD_LABELS[key], alpha=0.4)
    ax.set_xlabel("Nucleus area (µm²)")
    ax.set_ylabel("Count")
    ax.set_title("Nuclear mask size comparison", fontweight="bold")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "nuclear_mask_sizes.png", dpi=DPI)
    plt.close(fig)


def fig_transcripts_per_cell() -> None:
    counts = pd.read_csv(TABLES_DIR / "cell_counts.csv", index_col="method")
    methods = [m for m in MAIN_METHODS if m in counts.index]

    FILE_KEY = {"10x_native": "10x"}
    adatas = {}
    for m in methods:
        fname = f"adata_{FILE_KEY.get(m, m)}.h5ad"
        adatas[m] = ad.read_h5ad(ROI_DIR / fname)

    transcripts_by_method = {
        m: np.asarray(adatas[m].X.sum(axis=1)).ravel()
        for m in methods
    }

    fig, ax = plt.subplots(figsize=(max(14, 2 * len(methods)), 10))
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
        cut=0, inner="box", legend=False, ax=ax,
    )
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels([METHOD_LABELS[m] for m in methods], rotation=40, ha="right")
    ax.set_xlabel("")
    tick_vals = [1, 10, 100, 1000]
    ax.set_yticks([np.log10(v) for v in tick_vals])
    ax.set_yticklabels([str(v) for v in tick_vals])
    ax.set_ylabel("Transcripts per cell")
    ax.set_ylim(bottom=ax.get_ylim()[0] - 0.15)
    ax.set_title("Transcripts/cell distribution (expansion methods)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "transcripts_per_cell.png", dpi=DPI)
    plt.close(fig)


def fig_expression_correlation() -> None:
    rows = []
    for m, label in _available_comparisons():
        corr = pd.read_csv(TABLES_DIR / f"expression_correlation_10x_{m}.csv")
        vals = corr["correlation"].dropna()
        for v in vals:
            rows.append({"method": m, "label": label, "correlation": v})
    df = pd.DataFrame(rows)

    medians = df.groupby("method")["correlation"].median().sort_values(ascending=False)
    order = [m for m in medians.index]
    labels = [METHOD_LABELS[m] for m in order]
    colors = [METHOD_COLORS[m] for m in order]

    fig, ax = plt.subplots(figsize=(PANEL_W, max(PANEL_H, 0.7 * len(order) + 2)))
    sns.violinplot(
        data=df, y="method", x="correlation", hue="method", order=order,
        hue_order=order, palette=dict(zip(order, colors)),
        inner=None, linewidth=0.8, cut=0, legend=False, ax=ax,
    )
    sns.boxplot(
        data=df, y="method", x="correlation", order=order,
        width=0.15, fliersize=0, linewidth=1.2,
        boxprops=dict(facecolor="white", zorder=3),
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="black"), capprops=dict(color="black"),
        ax=ax,
    )
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Pearson correlation with 10x native")
    ax.set_ylabel("")
    ax.set_title("Per-cell expression agreement vs. 10x native", fontweight="bold")

    for i, m in enumerate(order):
        med = medians[m]
        ax.annotate(f"{med:.3f}", xy=(med, i), xytext=(6, 0),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=8, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              ec="none", alpha=0.75))

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "expression_correlation.png", dpi=DPI,
                bbox_inches="tight")
    plt.close(fig)


def fig_disagreement_spatial_map() -> None:
    _MATCHER_SUFFIXES = {"hungarian": "", "argmax": "_argmax"}
    for matcher_name, suffix in _MATCHER_SUFFIXES.items():
        matcher_label = "Hungarian (one-to-one)" if matcher_name == "hungarian" else "Argmax (many-to-one)"
        comparisons = _available_comparisons()
        pairs = []
        for m, label in comparisons:
            path = TABLES_DIR / f"disagreement_table_10x_{m}{suffix}.csv"
            if path.exists():
                pairs.append((label, pd.read_csv(path)))
        if not pairs:
            continue

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

        fig.suptitle(f"Cell-type agreement vs. disagreement - {matcher_label}", fontweight="bold")
        fig.legend(handles=[
            mpatches.Patch(color="#4C72B0", label="Agree"),
            mpatches.Patch(color="#C44E52", label="Disagree"),
        ], loc="lower center", ncols=2, framealpha=0.9)
        fig.tight_layout(rect=[0, 0.04, 1, 0.95])
        out_name = "disagreement_spatial_map.png" if not suffix else f"disagreement_spatial_map{suffix}.png"
        fig.savefig(FIGURES_DIR / out_name, dpi=DPI)
        plt.close(fig)


def fig_density_vs_disagreement() -> None:
    log_density = pd.read_csv(TABLES_DIR / "10x_log_density.csv", index_col=0)["log_density"]

    _MATCHER_SUFFIXES = {"hungarian": "", "argmax": "_argmax"}
    for matcher_name, suffix in _MATCHER_SUFFIXES.items():
        matcher_label = "Hungarian (one-to-one)" if matcher_name == "hungarian" else "Argmax (many-to-one)"
        summary_path = TABLES_DIR / f"density_disagreement_summary{suffix}.csv"
        summary = pd.read_csv(summary_path, index_col="comparison") if summary_path.exists() else pd.DataFrame()

        comparisons = _available_comparisons()
        pairs = []
        for m, label in comparisons:
            path = TABLES_DIR / f"disagreement_table_10x_{m}{suffix}.csv"
            if path.exists():
                pairs.append((m, label, pd.read_csv(path)))
        if not pairs:
            continue

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
            p = summary.loc[csv_key, "p_value"] if not summary.empty and csv_key and csv_key in summary.index else float("nan")
            ax.set_title(f"10x native vs. {label}", fontweight="bold")
            ax.text(0.04, 0.94, f"p = {p:.1e}", transform=ax.transAxes,
                    va="top", ha="left",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
            ax.set_xlabel("Phenotypic log-density (Mellon)")

        flat[0].set_ylabel("Density")
        fig.suptitle(f"Phenotypic density vs. cell-type disagreement - {matcher_label}", fontweight="bold")
        fig.legend(handles=[
            mpatches.Patch(color="#4C72B0", alpha=0.5, label="Agree"),
            mpatches.Patch(color="#C44E52", alpha=0.5, label="Disagree"),
        ], loc="lower center", ncols=2, framealpha=0.9)
        fig.tight_layout(rect=[0, 0.04, 1, 0.95])
        out_name = "density_vs_disagreement.png" if not suffix else f"density_vs_disagreement{suffix}.png"
        fig.savefig(FIGURES_DIR / out_name, dpi=DPI)
        plt.close(fig)


def _save_single_umap(
    emb: pd.DataFrame,
    color_labels: pd.Series,
    palette: dict[str, object],
    title: str,
    out_path: Path,
) -> None:
    apply_style(scatter=True)

    present = sorted(
        [k for k in color_labels.unique() if k != "unmatched"],
        key=lambda x: (0, int(x)) if str(x).lstrip("-").isdigit() else (1, str(x)),
    )
    if "unmatched" in color_labels.values:
        present.append("unmatched")
    pal = {k: palette.get(k, "#AAAAAA") for k in present}

    fig, ax = plt.subplots(figsize=(10, 8))
    for label_val in present:
        mask = color_labels == label_val
        alpha = 0.15 if label_val == "unmatched" else 0.6
        s = 6 if label_val == "unmatched" else 12
        ax.scatter(
            emb.loc[mask, "UMAP1"], emb.loc[mask, "UMAP2"],
            c=[pal[label_val]], s=s, alpha=alpha, label=label_val, rasterized=True,
        )
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(markerscale=2, fontsize=7, ncol=3, loc="lower right",
              framealpha=0.8, handletextpad=0.3, columnspacing=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_pca_umap() -> None:
    umap_dir = FIGURES_DIR / "umap"
    umap_dir.mkdir(parents=True, exist_ok=True)

    methods = [m for m in MAIN_METHODS if m != "10x_native"]
    methods = [m for m in methods if (TABLES_DIR / f"embedding_{m}.csv").exists()]

    emb_ref = pd.read_csv(TABLES_DIR / "embedding_10x_native.csv", index_col=0)
    n_ref = emb_ref["leiden"].nunique()
    palette = {str(c): sns.color_palette("tab20", n_ref)[c] for c in range(n_ref)}
    palette["unmatched"] = "#DDDDDD"

    _save_single_umap(emb_ref, emb_ref["leiden"].astype(str), palette,
                       f"10x native ({n_ref} clusters)",
                       umap_dir / "umap_10x_native.png")

    for method in methods:
        emb = pd.read_csv(TABLES_DIR / f"embedding_{method}.csv", index_col=0)
        label = METHOD_LABELS[method]
        n_clusters_method = emb["leiden"].nunique()

        for alg, suffix in [("hungarian", ""), ("argmax", "_argmax")]:
            dt_path = TABLES_DIR / f"disagreement_table_10x_{method}{suffix}.csv"
            if not dt_path.exists():
                continue
            dt = pd.read_csv(dt_path)
            aligned = dt.set_index("id_b")["label_b"].astype(str)
            color_col = emb.index.astype(str).to_series().map(aligned).fillna("unmatched")
            out = umap_dir / f"umap_{method}_{alg}.png"
            _save_single_umap(emb, color_col, palette,
                              f"{label} ({alg}, {n_clusters_method} clusters)", out)

    print(f"  Individual UMAPs saved to {umap_dir}/")


def fig_local_morans_map() -> None:
    LISA_COLORS = {"HH": "#C44E52", "LL": "#4C72B0", "HL": "#DD8452", "LH": "#CCB974"}
    _MATCHER_SUFFIXES = {"hungarian": "", "argmax": "_argmax"}
    for matcher_name, suffix in _MATCHER_SUFFIXES.items():
        matcher_label = "Hungarian (one-to-one)" if matcher_name == "hungarian" else "Argmax (many-to-one)"
        pairs = [
            (label, pd.read_csv(TABLES_DIR / f"local_morans_10x_{m}{suffix}.csv"))
            for m, label in _available_comparisons()
            if (TABLES_DIR / f"local_morans_10x_{m}{suffix}.csv").exists()
        ]
        if not pairs:
            continue

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

        fig.suptitle(f"Local Moran's I - {matcher_label}",
                     fontweight="bold")
        fig.legend(handles=[mpatches.Patch(color=color, label=cluster)
                            for cluster, color in LISA_COLORS.items()],
                   title="LISA cluster", loc="lower center", ncols=4, framealpha=0.9)
        fig.tight_layout(rect=[0, 0.04, 1, 0.95])
        out_name = "local_morans_map.png" if not suffix else f"local_morans_map{suffix}.png"
        fig.savefig(FIGURES_DIR / out_name, dpi=DPI)
        plt.close(fig)


def fig_de_volcano() -> None:
    _MATCHER_SUFFIXES = {"hungarian": "", "argmax": "_argmax"}
    for matcher_name, suffix in _MATCHER_SUFFIXES.items():
        matcher_label = "Hungarian (one-to-one)" if matcher_name == "hungarian" else "Argmax (many-to-one)"
        pairs = [
            (m, label, pd.read_csv(TABLES_DIR / f"de_disagree_10x_{m}{suffix}.csv"))
            for m, label in _available_comparisons()
            if (TABLES_DIR / f"de_disagree_10x_{m}{suffix}.csv").exists()
        ]
        if not pairs:
            continue

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
        fig.suptitle(f"DE: disagree vs. agree cells - {matcher_label}", fontweight="bold")
        fig.legend(handles=[
            mpatches.Patch(color="#AAAAAA", label="n.s."),
            mpatches.Patch(color="#C44E52", label="adj. p < 0.05"),
        ], loc="lower center", ncols=2, framealpha=0.9)
        fig.tight_layout(rect=[0, 0.04, 1, 0.95])
        out_name = "de_volcano.png" if not suffix else f"de_volcano{suffix}.png"
        fig.savefig(FIGURES_DIR / out_name, dpi=DPI)
        plt.close(fig)


def fig_spatial_morans_dotplot() -> None:
    """Dumbbell plot of Global Moran's I (Hungarian vs Argmax) for all methods."""
    import json

    ALL_METHODS_ORDERED = [
        "cellpose", "stardist", "mesmer", "10x_ranger",
        "voronoi", "voronoi_stardist", "voronoi_mesmer", "voronoi_10x_ranger",
        "baysor", "baysor_prior",
        "baysor_prior_c08", "baysor_prior_c10",
        "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10",
        "baysor_10x_ranger_prior_c10",
    ]

    rows = []
    for m in ALL_METHODS_ORDERED:
        h_path = TABLES_DIR / f"disagreement_spatial_10x_{m}.json"
        a_path = TABLES_DIR / f"disagreement_spatial_10x_{m}_argmax.json"
        if not h_path.exists() or not a_path.exists():
            continue
        with open(h_path) as f:
            h = json.load(f)
        with open(a_path) as f:
            a = json.load(f)
        rows.append({
            "method": METHOD_LABELS.get(m, m),
            "hungarian": h["morans_i"],
            "argmax": a["morans_i"],
        })

    if not rows:
        return

    df = pd.DataFrame(rows)
    n = len(df)
    y = np.arange(n)

    fig, ax = plt.subplots(figsize=(10, max(6, n * 0.5)))

    for i in range(n):
        ax.plot([df.iloc[i]["hungarian"], df.iloc[i]["argmax"]], [y[i], y[i]],
                color="#8fbc8f", linewidth=2, zorder=1)

    ax.scatter(df["hungarian"], y, s=60, c="#55A868", marker="s",
               edgecolors="black", linewidth=0.5, zorder=2, label="Hungarian")
    ax.scatter(df["argmax"], y, s=60, c="#8fbc8f", marker="D",
               edgecolors="black", linewidth=0.5, zorder=2, label="Argmax")

    for i in range(n):
        h_val = df.iloc[i]["hungarian"]
        a_val = df.iloc[i]["argmax"]
        if abs(h_val - a_val) < 0.015:
            left, right = (h_val, a_val) if h_val < a_val else (a_val, h_val)
            left_is_h = h_val <= a_val
            ax.annotate(f"{left:.3f}", (left, y[i]),
                        textcoords="offset points", xytext=(-8, 0),
                        ha="right", va="center", fontsize=7, color="#555555")
            ax.annotate(f"{right:.3f}", (right, y[i]),
                        textcoords="offset points", xytext=(8, 0),
                        ha="left", va="center", fontsize=7, color="#555555")
        else:
            ax.annotate(f"{h_val:.3f}", (h_val, y[i]),
                        textcoords="offset points", xytext=(0, 8),
                        ha="center", fontsize=7, color="#555555")
            ax.annotate(f"{a_val:.3f}", (a_val, y[i]),
                        textcoords="offset points", xytext=(0, -12),
                        ha="center", fontsize=7, color="#555555")

    ax.set_yticks(y)
    ax.set_yticklabels(df["method"])
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlabel("Global Moran's I")
    ax.set_title("Spatial autocorrelation of disagreement vs. 10x native",
                 fontweight="bold", fontsize=13)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9, framealpha=0.9)
    ax.grid(True, axis="x", alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "spatial_morans_dotplot.png", dpi=DPI,
                bbox_inches="tight")
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

    hungarian_fn = lambda vals: set(zip(*[x.tolist() for x in linear_sum_assignment(-vals)]))
    argmax_fn = lambda vals: {(int(vals[:, c].argmax()), c) for c in range(vals.shape[1])}

    ncols = 2
    nrows = math.ceil(len(avail) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 12, nrows * 9))
    flat = np.array(axes).flatten()
    for ax in flat[len(avail):]:
        ax.set_visible(False)

    for ax, (method, label) in zip(flat, avail):
        path = TABLES_DIR / f"cell_type_confusion_10x_{method}.csv"
        raw = pd.read_csv(path, index_col="label_a")
        raw.index = raw.index.astype(str)
        raw.columns = raw.columns.astype(str)

        ref_ids = sorted(raw.index, key=int)
        comp_ids = sorted(raw.columns, key=int)
        raw = raw.loc[ref_ids, comp_ids]

        vals = raw.to_numpy()
        hungarian_matched = hungarian_fn(vals)
        argmax_matched = argmax_fn(vals)

        row_sums = raw.sum(axis=1).replace(0, np.nan)
        norm = raw.div(row_sums, axis=0).fillna(0) * 100

        annot_text = norm.map(lambda v: f"{v:.0f}" if v >= 5 else "")

        row_labels = [f"{c}: {ct_short.get(CLUSTER_ANNOTATIONS.get(c, ''), c)}"
                      for c in ref_ids]

        x_labels = [str(c) if i % 2 == 0 else "" for i, c in enumerate(comp_ids)]

        sns.heatmap(
            np.zeros_like(norm.values), ax=ax,
            cmap="Greys", vmin=0, vmax=1,
            annot=np.array(annot_text), fmt="",
            annot_kws={"weight": "bold"},
            xticklabels=x_labels, yticklabels=row_labels,
            linewidths=0.4, linecolor="#e0e0e0",
            cbar=False,
            square=True,
        )

        matched = set(argmax_matched) | set(hungarian_matched)
        for r in range(norm.shape[0]):
            for c in range(norm.shape[1]):
                if (r, c) in matched or norm.values[r, c] < 1:
                    continue
                alpha = np.clip(norm.values[r, c] / 100, 0.03, 0.3)
                ax.add_patch(Rectangle((c, r), 1, 1, fill=True,
                                       facecolor="#888888", alpha=alpha,
                                       linewidth=0))

        for r, c in sorted(matched):
            both = (r, c) in argmax_matched and (r, c) in hungarian_matched
            color = "purple" if both else "red" if (r, c) in hungarian_matched else "royalblue"
            alpha = np.clip(norm.values[r, c] / 100, 0.08, 0.9)
            ax.add_patch(Rectangle((c, r), 1, 1, fill=True,
                                   facecolor=color, edgecolor=color,
                                   alpha=alpha, linewidth=1.5))

        ax.set_title(f"{label}  ({len(comp_ids)} clusters)", fontweight="bold")
        ax.set_xlabel(f"{label} cluster")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=0)
        ax.tick_params(axis="y", rotation=0)

    legend_handles = [
        mpatches.Patch(edgecolor="red", facecolor="red", alpha=0.5, label="Hungarian (one-to-one)"),
        mpatches.Patch(edgecolor="royalblue", facecolor="royalblue", alpha=0.5, label="Argmax (many-to-one)"),
        mpatches.Patch(edgecolor="purple", facecolor="purple", alpha=0.5, label="Overlap (both match)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               frameon=True, fontsize=16, bbox_to_anchor=(0.48, -0.005),
               handlelength=3, handleheight=2)

    fig.suptitle(
        "Cluster-level confusion matrices (row-normalised)",
        fontstyle="italic", fontweight="bold",
    )
    fig.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.04,
                        hspace=0.45, wspace=0.35)
    fig.savefig(FIGURES_DIR / "confusion_clusters.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig_nuclear_mask_sizes()
    fig_transcripts_per_cell()
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
    fig_spatial_morans_dotplot()
    print(f"wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
