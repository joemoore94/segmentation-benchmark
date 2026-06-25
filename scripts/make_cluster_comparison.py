"""Compare Hungarian vs argmax cluster alignment across segmentation methods.

Grouped bar chart showing per-method disagreement rate under each alignment
strategy (Hungarian one-to-one vs argmax many-to-one).

Reads:  results/tables/disagreement_table_10x_*.csv
Writes: results/figures/cluster_comparison.png

Usage::

    conda run -n segbench python scripts/make_cluster_comparison.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from segbench.constants import (
    MAIN_METHODS,
    METHOD_FAMILIES,
    METHOD_LABELS,
)
from segbench.style import apply_style

TABLES_DIR = Path("results/tables")
FIGURES = Path("results/figures")

METHODS = [m for m in MAIN_METHODS if m != "10x_native"]

FAMILY_COLORS = {
    "Reference": "#333333",
    "Voronoi": "#1B9E77",
    "Transcript-density": "#D95F02",
}


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    rows: list[dict] = []

    for method in METHODS:
        label = METHOD_LABELS[method]
        family = METHOD_FAMILIES[method]

        hungarian_path = TABLES_DIR / f"disagreement_table_10x_{method}.csv"
        argmax_path = TABLES_DIR / f"disagreement_table_10x_{method}_argmax.csv"
        if not hungarian_path.exists() or not argmax_path.exists():
            print(f"  {label}: skipped (missing files)")
            continue

        df_h = pd.read_csv(hungarian_path)
        df_a = pd.read_csv(argmax_path)

        rate_h = df_h["disagree"].mean()
        rate_a = df_a["disagree"].mean()

        rows.append({
            "method": label,
            "method_key": method,
            "family": family,
            "hungarian": rate_h,
            "argmax": rate_a,
        })
        print(f"  {label}: hungarian={rate_h:.3f}  argmax={rate_a:.3f}")

    df = pd.DataFrame(rows)
    method_order = df["method"].tolist()

    fig, ax = plt.subplots(figsize=(18, 12))

    y_pos = np.arange(len(method_order))
    bar_h = 0.35

    vals_hungarian = df["hungarian"].values * 100
    vals_argmax = df["argmax"].values * 100

    bars_h = ax.barh(
        y_pos - bar_h / 2, vals_hungarian, bar_h,
        label="Hungarian (1-to-1)", color="#4C72B0", edgecolor="white", linewidth=0.5,
    )
    bars_a = ax.barh(
        y_pos + bar_h / 2, vals_argmax, bar_h,
        label="Argmax (many-to-1)", color="#DD8452", edgecolor="white", linewidth=0.5,
    )

    for bar_set in [bars_h, bars_a]:
        for bar in bar_set:
            w = bar.get_width()
            ax.text(w + 0.3, bar.get_y() + bar.get_height() / 2,
                    f"{w:.1f}", va="center", fontsize=16)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(method_order)
    ax.set_xlabel("Disagreement rate (%)")
    ax.invert_yaxis()
    ax.legend(loc="lower right", fontsize=20)
    ax.set_xlim(0, max(vals_hungarian.max(), vals_argmax.max()) + 5)

    fig.suptitle(
        "Cluster alignment strategy: Hungarian vs argmax disagreement",
        fontsize=28, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "cluster_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {FIGURES / 'cluster_comparison.png'}")


if __name__ == "__main__":
    main()
