"""Clustering agreement vs 10x native — grouped vertical bar chart.

Shows ARI (×100), Hungarian disagreement %, and argmax disagreement % for
every method compared to 10x native.

Reads:  results/tables/disagreement_table_10x_*.csv
        results/tables/pairwise_ari.csv
Writes: results/figures/clustering_agreement.png

Usage::

    conda run -n segbench python scripts/make_cluster_comparison.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from segbench.constants import METHOD_LABELS
from segbench.style import apply_style

TABLES_DIR = Path("results/tables")
FIGURES = Path("results/figures")

ALL_METHODS = [
    "cellpose", "stardist", "mesmer", "10x_ranger",
    "voronoi", "voronoi_stardist", "voronoi_mesmer", "voronoi_10x_ranger",
    "baysor", "baysor_prior", "baysor_prior_c08", "baysor_prior_c10",
    "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10",
    "baysor_10x_ranger_prior_c10",
]


def _load_ari_lookup() -> dict[str, float]:
    lookup: dict[str, float] = {}
    res_path = TABLES_DIR / "resolution_ari.csv"
    if res_path.exists():
        res = pd.read_csv(res_path)
        at_1 = res[res["resolution"] == 1.0]
        for _, row in at_1.iterrows():
            lookup[row["method"]] = row["ari"]
    pw_path = TABLES_DIR / "pairwise_ari.csv"
    if pw_path.exists():
        pw = pd.read_csv(pw_path)
        for _, row in pw[pw["method_a"] == "10x native"].iterrows():
            if row["method_b"] not in lookup:
                lookup[row["method_b"]] = row["ari"]
    return lookup


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    ari_lookup = _load_ari_lookup()
    rows: list[dict] = []

    for method in ALL_METHODS:
        label = METHOD_LABELS[method]

        hungarian_path = TABLES_DIR / f"disagreement_table_10x_{method}.csv"
        argmax_path = TABLES_DIR / f"disagreement_table_10x_{method}_argmax.csv"
        if not hungarian_path.exists() or not argmax_path.exists():
            print(f"  {label}: skipped (missing files)")
            continue

        df_h = pd.read_csv(hungarian_path)
        df_a = pd.read_csv(argmax_path)

        rate_h = df_h["disagree"].mean() * 100
        rate_a = df_a["disagree"].mean() * 100
        ari = ari_lookup.get(label, float("nan")) * 100

        rows.append({
            "method": label,
            "ari": ari,
            "hungarian": rate_h,
            "argmax": rate_a,
        })
        print(f"  {label}: ARI={ari:.1f}  hungarian={rate_h:.1f}%  argmax={rate_a:.1f}%")

    df = pd.DataFrame(rows)

    x = np.arange(len(df))
    bar_w = 0.25

    fig, ax = plt.subplots(figsize=(max(20, len(df) * 2), 14))

    bars_ari = ax.bar(
        x - bar_w, df["ari"].values, bar_w,
        label="ARI (×100)", color="#4C72B0", edgecolor="white", linewidth=0.5,
    )
    bars_h = ax.bar(
        x, df["hungarian"].values, bar_w,
        label="Disagree % (Hungarian)", color="#C44E52", edgecolor="white", linewidth=0.5,
    )
    bars_a = ax.bar(
        x + bar_w, df["argmax"].values, bar_w,
        label="Disagree % (Argmax)", color="#DD8452", edgecolor="white", linewidth=0.5,
    )

    for bars in [bars_ari, bars_h, bars_a]:
        for bar in bars:
            h = bar.get_height()
            if np.isnan(h):
                continue
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.8,
                    f"{h:.1f}", ha="center", va="bottom", fontsize=12,
                    rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(df["method"].tolist(), rotation=45, ha="right")
    ax.set_ylabel("Percent")
    ax.legend(loc="upper right", fontsize=18)
    ax.set_ylim(0, df[["ari", "hungarian", "argmax"]].max(skipna=True).max() + 8)

    fig.suptitle(
        "Clustering agreement vs. 10x native",
        fontsize=28, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "clustering_agreement.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {FIGURES / 'clustering_agreement.png'}")


if __name__ == "__main__":
    main()
