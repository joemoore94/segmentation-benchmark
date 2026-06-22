"""Cell size vs. disagreement: does nucleus or cell area predict method mismatch?

For each matched pair (10x native vs. comparison method), this script joins the
10x-native cell area (nucleus_area_um2 and cell_area_um2) and tests whether
smaller or larger cells are more likely to disagree.

Expected signal for nuclear-only methods: larger 10x-native cells have more
cytoplasmic volume — transcripts that nuclear masks miss — so they should show
higher disagreement rates. Voronoi methods capture all transcripts regardless of
cell size, so this relationship should be flat.

Left panel: binned disagree probability vs. 10x-native cell area for all methods.
Right panel: KDE of 10x-native cell area, coloured by agree/disagree, for the
             three nuclear methods combined.

Reads:  data/processed/roi/adata_10x.h5ad
        results/tables/disagreement_table_10x_*.csv
Writes: results/figures/cell_size_disagreement.png

Usage::

    conda run -n segbench python scripts/run_cell_size_disagreement.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES  = Path("results/tables")
FIGURES = Path("results/figures")

COMPARISONS = [
    ("cellpose",          "CellPose",     "#4C72B0", "Nuclear"),
    ("stardist",          "StarDist",     "#8172B2", "Nuclear"),
    ("mesmer",            "Mesmer",       "#D62728", "Nuclear"),
    ("voronoi",           "Voronoi (CP)", "#17BECF", "Voronoi"),
    ("voronoi_stardist",  "Voronoi (SD)", "#9467BD", "Voronoi"),
    ("voronoi_mesmer",    "Voronoi (M)",  "#BCBD22", "Voronoi"),
    ("baysor",            "Baysor",       "#DD8452", "Transcript-density"),
]

PLOT_COMPARISONS = [c for c in COMPARISONS if c[0] not in ("cellpose", "stardist", "mesmer")]

N_BINS = 10


def load_disagree_with_area(adata_10x: ad.AnnData, method: str) -> pd.DataFrame:
    path = TABLES / f"disagreement_table_10x_{method}.csv"
    df = pd.read_csv(path)
    area = adata_10x.obs[["cell_area_um2", "nucleus_area_um2"]]
    df = df.join(area, on="id_a")
    return df


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    print("Loading 10x native AnnData...")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")

    all_dfs: dict[str, pd.DataFrame] = {}
    for method, label, _, _ in COMPARISONS:
        path = TABLES / f"disagreement_table_10x_{method}.csv"
        if not path.exists():
            print(f"  Missing {path.name}, skipping")
            continue
        df = load_disagree_with_area(adata_10x, method)
        all_dfs[label] = df
        stat, p = mannwhitneyu(
            df.loc[df["disagree"] == 0, "cell_area_um2"].dropna(),
            df.loc[df["disagree"] == 1, "cell_area_um2"].dropna(),
            alternative="two-sided",
        )
        agree_med = df.loc[df["disagree"] == 0, "cell_area_um2"].median()
        disagree_med = df.loc[df["disagree"] == 1, "cell_area_um2"].median()
        print(f"  {label}: median area agree={agree_med:.1f} disagree={disagree_med:.1f} µm², p={p:.2e}")

    # Common area range for binning (use 10x native full distribution)
    area_col = "cell_area_um2"
    all_areas = np.concatenate([df[area_col].dropna().values for df in all_dfs.values()])
    bin_edges = np.percentile(all_areas, np.linspace(0, 100, N_BINS + 1))
    bin_edges = np.unique(bin_edges)
    bin_mids  = (bin_edges[:-1] + bin_edges[1:]) / 2

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(22, 9))

    # Left: binned disagree probability vs. cell area (plot subset)
    for method, label, color, family in PLOT_COMPARISONS:
        if label not in all_dfs:
            continue
        df = all_dfs[label]
        df = df.dropna(subset=[area_col])
        df["area_bin"] = pd.cut(df[area_col], bins=bin_edges, labels=bin_mids, include_lowest=True)
        binned = df.groupby("area_bin", observed=True)["disagree"].mean()
        ls = "--" if family == "Voronoi" else ("-." if family == "Transcript-density" else "-")
        ax_left.plot(binned.index.astype(float), binned.values,
                     color=color, lw=2.5, ls=ls, marker="o", markersize=6,
                     label=label)

    ax_left.set_xlabel("10x-native cell area (µm²)")
    ax_left.set_ylabel("Disagree probability")
    ax_left.set_title("Disagree probability vs. cell area", fontweight="bold")
    ax_left.legend(fontsize=10)
    ax_left.set_ylim(0, None)

    # Right: KDE of cell area for agree vs. disagree — nuclear methods pooled
    nuclear_methods = [l for _, l, _, fam in COMPARISONS if fam == "Nuclear" and l in all_dfs]
    agree_areas    = np.concatenate([all_dfs[l].loc[all_dfs[l]["disagree"] == 0, area_col].dropna().values
                                     for l in nuclear_methods])
    disagree_areas = np.concatenate([all_dfs[l].loc[all_dfs[l]["disagree"] == 1, area_col].dropna().values
                                     for l in nuclear_methods])

    sns.kdeplot(agree_areas,    ax=ax_right, label="Agree",    color="#2CA02C", fill=True, alpha=0.35, log_scale=True)
    sns.kdeplot(disagree_areas, ax=ax_right, label="Disagree", color="#D62728", fill=True, alpha=0.35, log_scale=True)

    ax_right.set_xlabel("10x-native cell area (µm², log scale)")
    ax_right.set_ylabel("Density")
    ax_right.set_title("Cell area: agree vs. disagree\n(nuclear methods pooled)", fontweight="bold")
    ax_right.legend(fontsize=11)

    fig.suptitle(
        "Does 10x-native cell size predict disagreement with each segmentation method?",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "cell_size_disagreement.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved cell_size_disagreement.png")


if __name__ == "__main__":
    main()
