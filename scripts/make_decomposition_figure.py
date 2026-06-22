"""ARI decomposition figure: nuclear detector quality × cytoplasmic coverage.

The two Voronoi controls isolate two independent sources of improvement over
nuclear-only segmentation:

  CellPose nuclear  ──(+cytoplasm)──▶  Voronoi (CP)  ──(+nuclei)──▶  Voronoi (M)
       0.547                                0.630                          0.686

Left panel: horizontal bar chart of all methods' ARI, with bracketed gain
            annotations showing the decomposition path.
Right panel: ARI, agreement rate, and transcript capture for the three
             decomposition steps.

Usage::

    conda run -n segbench python scripts/make_decomposition_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from segbench.style import apply_style

FIGURES = Path("results/figures")

# All methods sorted by ARI (bottom to top in a horizontal bar chart).
ALL_METHODS = [
    {"label": "Baysor",        "ari": 0.305, "color": "#DD8452", "family": "Transcript-density"},
    {"label": "StarDist",      "ari": 0.545, "color": "#8172B2", "family": "Nuclear"},
    {"label": "CellPose",      "ari": 0.547, "color": "#4C72B0", "family": "Nuclear"},
    {"label": "Mesmer",        "ari": 0.557, "color": "#D62728", "family": "Nuclear"},
    {"label": "Voronoi (CP)",  "ari": 0.630, "color": "#17BECF", "family": "Voronoi"},
    {"label": "Voronoi (M)",   "ari": 0.686, "color": "#BCBD22", "family": "Voronoi"},
]

# The three decomposition steps (subset of ALL_METHODS).
STEPS = [
    {"short": "CP Nuclear",   "ari": 0.547, "disagree": 0.308, "capture": 0.354, "color": "#4C72B0"},
    {"short": "Voronoi (CP)", "ari": 0.630, "disagree": 0.219, "capture": 1.000, "color": "#17BECF"},
    {"short": "Voronoi (M)",  "ari": 0.686, "disagree": 0.188, "capture": 1.000, "color": "#D62728"},
]

# Y-position indices in ALL_METHODS for the decomposition path.
PATH_INDICES = {m["label"]: i for i, m in enumerate(ALL_METHODS)}
CP_IDX       = PATH_INDICES["CellPose"]
VCP_IDX      = PATH_INDICES["Voronoi (CP)"]
VM_IDX       = PATH_INDICES["Voronoi (M)"]


def bracket(ax, y0, y1, x, label, color="#444444", fontsize=10):
    """Draw a right-side bracket from y0 to y1 at x, labelled with `label`."""
    pad = 0.008
    ax.annotate("", xy=(x + pad, y1), xytext=(x + pad, y0),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.8))
    ax.text(x + pad + 0.012, (y0 + y1) / 2, label,
            va="center", ha="left", fontsize=fontsize,
            color=color, fontweight="bold")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_style()

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(22, 9),
                                             gridspec_kw={"width_ratios": [1.2, 1]})

    # ---------------------------------------------------------------- Left: all-methods ARI bar chart
    labels  = [m["label"] for m in ALL_METHODS]
    aris    = [m["ari"]   for m in ALL_METHODS]
    colors  = [m["color"] for m in ALL_METHODS]
    ys      = np.arange(len(ALL_METHODS))

    bars = ax_left.barh(ys, aris, color=colors, height=0.6, edgecolor="white", zorder=3)

    # Value labels at the end of each bar
    for y, ari, m in zip(ys, aris, ALL_METHODS):
        ax_left.text(ari + 0.005, y, f"{ari:.3f}",
                     va="center", ha="left", fontsize=10, fontweight="bold",
                     color=m["color"])

    # Highlight decomposition path bars with a thicker border
    for idx in [CP_IDX, VCP_IDX, VM_IDX]:
        bar = bars[idx]
        bar.set_linewidth(2.5)
        bar.set_edgecolor("#333333")

    # Gain brackets to the right of the bars
    x_bracket = 0.72
    bracket(ax_left, CP_IDX, VCP_IDX, x_bracket,
            f"+{aris[VCP_IDX] - aris[CP_IDX]:.3f}\n+cytoplasmic\ncoverage")
    bracket(ax_left, VCP_IDX, VM_IDX, x_bracket,
            f"+{aris[VM_IDX] - aris[VCP_IDX]:.3f}\n+better nuclei")

    ax_left.set_yticks(ys)
    ax_left.set_yticklabels(labels, fontsize=11)
    ax_left.set_xlabel("ARI vs. 10x native")
    ax_left.set_xlim(0, 0.88)
    ax_left.set_title("ARI vs. 10x native — all methods", fontweight="bold")

    # Family legend
    family_colors = {"Nuclear": "#4C72B0", "Voronoi": "#17BECF", "Transcript-density": "#DD8452"}
    handles = [mpatches.Patch(color=c, label=f) for f, c in family_colors.items()]
    ax_left.legend(handles=handles, fontsize=10, loc="lower right")

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
