"""Compare Hungarian vs argmax cluster alignment across segmentation methods.

Two-panel figure: per-method disagreement rate under each alignment strategy
(left, grouped bar chart) and the difference (argmax minus Hungarian) per
method (right, diverging bar chart).

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

    families = df.set_index("method")["family"].reindex(method_order)
    bar_colors = [FAMILY_COLORS[f] for f in families]

    fig, (ax_left, ax_right) = plt.subplots(
        1, 2, figsize=(28, 10),
        gridspec_kw={"width_ratios": [1.8, 1]},
    )

    # --- Left panel: grouped bars, Hungarian vs argmax ---
    y_pos = np.arange(len(method_order))
    bar_h = 0.35

    vals_hungarian = df["hungarian"].values * 100
    vals_argmax = df["argmax"].values * 100

    bars_h = ax_left.barh(
        y_pos - bar_h / 2, vals_hungarian, bar_h,
        label="Hungarian (1-to-1)", color="#4C72B0", edgecolor="white", linewidth=0.5,
    )
    bars_a = ax_left.barh(
        y_pos + bar_h / 2, vals_argmax, bar_h,
        label="Argmax (many-to-1)", color="#DD8452", edgecolor="white", linewidth=0.5,
    )

    for bar_set in [bars_h, bars_a]:
        for bar in bar_set:
            w = bar.get_width()
            ax_left.text(w + 0.3, bar.get_y() + bar.get_height() / 2,
                         f"{w:.1f}", va="center", fontsize=14)

    ax_left.set_yticks(y_pos)
    ax_left.set_yticklabels(method_order)
    ax_left.set_xlabel("Disagreement rate (%)")
    ax_left.set_title("Hungarian vs argmax alignment", fontweight="bold")
    ax_left.invert_yaxis()
    ax_left.legend(loc="lower right", fontsize=18)
    ax_left.set_xlim(0, max(vals_hungarian.max(), vals_argmax.max()) + 5)

    # --- Right panel: difference (argmax - hungarian) ---
    diff = vals_argmax - vals_hungarian
    diff_colors = ["#2CA02C" if d < 0 else "#D62728" for d in diff]

    ax_right.barh(y_pos, diff, 0.6, color=diff_colors, edgecolor="white", linewidth=0.5)
    for i, d in enumerate(diff):
        ha = "left" if d >= 0 else "right"
        offset = 0.3 if d >= 0 else -0.3
        ax_right.text(d + offset, i, f"{d:+.1f}", va="center", ha=ha, fontsize=14)

    ax_right.set_yticks(y_pos)
    ax_right.set_yticklabels([])
    ax_right.set_xlabel("Difference (pp)")
    ax_right.set_title("Argmax − Hungarian", fontweight="bold")
    ax_right.axvline(0, color="black", linewidth=0.8, linestyle="-")
    ax_right.invert_yaxis()

    from matplotlib.patches import Patch
    ax_right.legend(
        handles=[
            Patch(facecolor="#2CA02C", label="Argmax lower"),
            Patch(facecolor="#D62728", label="Argmax higher"),
        ],
        loc="lower right", fontsize=16,
    )

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
