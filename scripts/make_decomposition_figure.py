"""ARI decomposition figure: nuclear detector quality × cytoplasmic coverage.

The two Voronoi controls isolate two independent sources of improvement over
nuclear-only segmentation:

  CellPose nuclear  ──(+cytoplasm)──▶  Voronoi (CP)  ──(+nuclei)──▶  Voronoi (M)
       0.547                                0.630                          0.686

Left panel: slope chart of the three-step path (ARI on Y axis).
Right panel: all three key metrics (ARI, disagree rate, transcript capture)
             at each step in a grouped bar chart.

Usage::

    conda run -n segbench python scripts/make_decomposition_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

TABLES = Path("results/tables")
FIGURES = Path("results/figures")

# The three steps in the decomposition path.
STEPS = [
    {
        "label":    "CellPose\nNuclear",
        "ari":      0.547,
        "disagree": 0.308,
        "capture":  0.354,
        "color":    "#4C72B0",
        "short":    "CP Nuclear",
    },
    {
        "label":    "Voronoi\n(CellPose)",
        "ari":      0.630,
        "disagree": 0.219,
        "capture":  1.000,
        "color":    "#17BECF",
        "short":    "Voronoi (CP)",
    },
    {
        "label":    "Voronoi\n(Mesmer)",
        "ari":      0.686,
        "disagree": 0.188,
        "capture":  1.000,
        "color":    "#D62728",
        "short":    "Voronoi (M)",
    },
]

# Context: other methods shown as reference dots on the slope chart.
OTHER = [
    {"label": "StarDist\nNuclear", "ari": 0.545, "color": "#8172B2"},
    {"label": "Mesmer\nNuclear",   "ari": 0.557, "color": "#C44E52"},
    {"label": "Baysor",            "ari": 0.305, "color": "#DD8452"},
    {"label": "10x native\n(ref)", "ari": 1.000, "color": "#55A868"},
]


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="poster")

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(20, 9),
                                             gridspec_kw={"width_ratios": [1.1, 1]})

    # ---------------------------------------------------------------- Left: slope / path chart
    xs = [0, 1, 2]
    aris = [s["ari"] for s in STEPS]
    colors = [s["color"] for s in STEPS]

    # Reference method dots at x=-0.5 (not on the path)
    for other in OTHER:
        ax_left.plot(-0.5, other["ari"], "D", color=other["color"], ms=10, zorder=3)
        ax_left.text(-0.42, other["ari"], other["label"],
                     va="center", ha="left", fontsize=9, color=other["color"])

    # Path line
    ax_left.plot(xs, aris, "-", color="#444444", linewidth=2.5, zorder=2)

    # Step nodes
    for x, step in zip(xs, STEPS):
        ax_left.plot(x, step["ari"], "o", color=step["color"], ms=18, zorder=4,
                     markeredgecolor="white", markeredgewidth=2)
        ax_left.text(x, step["ari"] - 0.035, f"ARI {step['ari']:.3f}",
                     ha="center", va="top", fontsize=10, fontweight="bold")
        ax_left.text(x, step["ari"] + 0.025, step["label"],
                     ha="center", va="bottom", fontsize=10)

    # Gain annotations on the arrows
    gains = [
        (0.5, (aris[0] + aris[1]) / 2, "+cytoplasmic\ncoverage", f"+{aris[1]-aris[0]:.3f} ARI"),
        (1.5, (aris[1] + aris[2]) / 2, "+better\nnuclei", f"+{aris[2]-aris[1]:.3f} ARI"),
    ]
    for gx, gy, label, delta in gains:
        ax_left.annotate(
            f"{label}\n{delta}",
            xy=(gx, gy), ha="center", va="center", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="#bbbb88", alpha=0.95),
        )

    ax_left.set_xlim(-0.85, 2.55)
    ax_left.set_ylim(0.22, 1.08)
    ax_left.set_xticks(xs)
    ax_left.set_xticklabels(["Step 1", "Step 2", "Step 3"], fontsize=11)
    ax_left.set_ylabel("ARI vs. 10x native")
    ax_left.set_title("ARI decomposition: what drives\nagreement with 10x native?",
                       fontweight="bold")
    ax_left.axhline(1.0, color="#55A868", linewidth=1, linestyle="--", alpha=0.5)
    ax_left.axvline(-0.25, color="#cccccc", linewidth=1, linestyle=":")

    # ---------------------------------------------------------------- Right: grouped bar chart
    metrics = ["ARI", "Agreement rate\n(1 − disagree)", "Transcript\ncapture"]
    step_labels = [s["short"] for s in STEPS]
    step_colors = [s["color"] for s in STEPS]

    vals = np.array([
        [s["ari"] for s in STEPS],
        [1 - s["disagree"] for s in STEPS],
        [s["capture"] for s in STEPS],
    ])

    x = np.arange(len(metrics))
    width = 0.22
    offsets = [-width, 0, width]

    for i, (step, color, offset) in enumerate(zip(step_labels, step_colors, offsets)):
        bars = ax_right.bar(x + offset, vals[:, i], width, label=step, color=color,
                            alpha=0.88, edgecolor="white")
        for bar, val in zip(bars, vals[:, i]):
            ax_right.text(bar.get_x() + bar.get_width() / 2,
                          bar.get_height() + 0.012,
                          f"{val:.0%}", ha="center", va="bottom", fontsize=9,
                          fontweight="bold")

    ax_right.set_xticks(x)
    ax_right.set_xticklabels(metrics, fontsize=11)
    ax_right.set_ylim(0, 1.18)
    ax_right.set_ylabel("Value (0–1 scale)")
    ax_right.set_title("Three metrics at each step",
                        fontweight="bold")
    ax_right.legend(title="Method", fontsize=10, title_fontsize=10)
    ax_right.axhline(1.0, color="#55A868", linewidth=1, linestyle="--", alpha=0.4,
                     label="10x native (ref)")

    fig.suptitle(
        "Decomposing agreement with 10x native: cytoplasmic coverage vs. nuclear detector quality",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "decomposition.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved decomposition.png")


if __name__ == "__main__":
    main()
