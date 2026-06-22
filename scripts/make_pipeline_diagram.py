"""Pipeline overview diagram.

Four rows, one per segmentation paradigm, converging on a shared evaluation
pipeline. Color-coded by method family. No external dependencies beyond
matplotlib.

Writes: results/figures/pipeline_diagram.png

Usage::

    conda run -n segbench python scripts/make_pipeline_diagram.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

FIGURES = Path("results/figures")

# ---------- layout constants ------------------------------------------------
W, H = 17, 6          # figure size in inches
XL, XR = 0, 17        # data x range
YB, YT = 0, 6         # data y range

# x breakpoints (all in data units 0–17)
X_APPROACH_L = 0.3    # left edge of approach box
X_APPROACH_R = 5.8    # right edge of approach box
X_GAP1       = 6.1    # arrow start
X_MATRIX_L   = 6.4    # left edge of "cell×gene" box
X_MATRIX_R   = 8.8    # right edge
X_GAP2       = 9.1    # arrow start
X_EVAL_L     = 9.4    # left edge of eval pipeline
X_EVAL_R     = 16.7   # right edge

# y centers for each paradigm row (top to bottom)
ROWS = {
    "nuclear":  5.0,
    "voronoi":  3.5,
    "baysor":   2.1,
    "ref":      0.7,
}
BOX_H = 0.85           # box height

FAMILY_COLORS = {
    "nuclear":  "#4C72B0",
    "voronoi":  "#17BECF",
    "baysor":   "#DD8452",
    "ref":      "#55A868",
}
FAMILY_ALPHA = 0.18    # background fill alpha

EVAL_STEPS = [
    ("Normalize\n+ log1p",        "#e0e0e0"),
    ("PCA\n(30 PCs)",             "#e0e0e0"),
    ("Leiden\nclustering",        "#e0e0e0"),
    ("Hungarian\nalignment",      "#e0e0e0"),
    ("ARI  ·  Moran's I\nPearson r (pseudobulk)", "#d4edda"),
]


# ---------- helpers ---------------------------------------------------------

def box(ax, xl, xr, yc, h, label, color, alpha=1.0,
        fontsize=10, bold=False, text_color="white", sublabel=None):
    rect = FancyBboxPatch(
        (xl, yc - h / 2), xr - xl, h,
        boxstyle="round,pad=0.05",
        facecolor=color, edgecolor="white",
        linewidth=1.2, alpha=alpha, zorder=3,
    )
    ax.add_patch(rect)
    weight = "bold" if bold else "normal"
    y_text = yc if sublabel is None else yc + h * 0.14
    ax.text((xl + xr) / 2, y_text, label,
            ha="center", va="center", fontsize=fontsize,
            color=text_color, fontweight=weight, zorder=4)
    if sublabel:
        ax.text((xl + xr) / 2, yc - h * 0.22, sublabel,
                ha="center", va="center", fontsize=fontsize - 1.5,
                color=text_color, alpha=0.85, zorder=4, style="italic")


def arrow(ax, x0, x1, y, color="#555555", lw=1.5):
    ax.annotate(
        "", xy=(x1, y), xytext=(x0, y),
        arrowprops=dict(
            arrowstyle="-|>", color=color, lw=lw,
            mutation_scale=12,
            connectionstyle="arc3,rad=0.0",
        ),
        zorder=5,
    )


def brace_arrow(ax, x0, x1, y_from, y_to, color="#555555"):
    """Arrow that starts at (x0, y_from) and ends at (x1, y_to)."""
    ax.annotate(
        "", xy=(x1, y_to), xytext=(x0, y_from),
        arrowprops=dict(
            arrowstyle="-|>", color=color, lw=1.3,
            mutation_scale=10,
            connectionstyle="arc3,rad=0.0",
        ),
        zorder=5,
    )


# ---------- main ------------------------------------------------------------

def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_xlim(XL, XR)
    ax.set_ylim(YB, YT)
    ax.axis("off")

    # -------- background family bands --------
    band_ranges = {
        "nuclear": (ROWS["nuclear"] - BOX_H * 0.6, ROWS["nuclear"] + BOX_H * 0.6),
        "voronoi": (ROWS["voronoi"] - BOX_H * 0.6, ROWS["voronoi"] + BOX_H * 0.6),
        "baysor":  (ROWS["baysor"]  - BOX_H * 0.6, ROWS["baysor"]  + BOX_H * 0.6),
        "ref":     (ROWS["ref"]     - BOX_H * 0.6, ROWS["ref"]     + BOX_H * 0.6),
    }
    for family, (ylo, yhi) in band_ranges.items():
        ax.axhspan(ylo, yhi, xmin=0.0, xmax=1.0,
                   color=FAMILY_COLORS[family], alpha=0.06, zorder=0)

    # -------- family labels on the left --------
    label_x = 0.12
    ax.text(label_x, ROWS["nuclear"], "Nuclear\n(pixel mask)", ha="center", va="center",
            fontsize=8.5, color=FAMILY_COLORS["nuclear"], fontweight="bold",
            rotation=90, zorder=4)
    ax.text(label_x, ROWS["voronoi"], "Voronoi\n(centroid)", ha="center", va="center",
            fontsize=8.5, color=FAMILY_COLORS["voronoi"], fontweight="bold",
            rotation=90, zorder=4)
    ax.text(label_x, ROWS["baysor"], "Transcript\ndensity", ha="center", va="center",
            fontsize=8.5, color=FAMILY_COLORS["baysor"], fontweight="bold",
            rotation=90, zorder=4)
    ax.text(label_x, ROWS["ref"], "Reference", ha="center", va="center",
            fontsize=8.5, color=FAMILY_COLORS["ref"], fontweight="bold",
            rotation=90, zorder=4)

    # -------- approach boxes --------
    box(ax, X_APPROACH_L, X_APPROACH_R, ROWS["nuclear"], BOX_H,
        "DAPI → CellPose / StarDist / Mesmer",
        FAMILY_COLORS["nuclear"], bold=True, fontsize=10.5,
        sublabel="nuclear masks · transcripts inside mask only (35–52% capture)")

    box(ax, X_APPROACH_L, X_APPROACH_R, ROWS["voronoi"], BOX_H,
        "Nuclear centroids → Voronoi (CellPose) · Voronoi (Mesmer)",
        FAMILY_COLORS["voronoi"], bold=True, fontsize=10.5,
        sublabel="nearest-centroid transcript assignment · 100% capture")

    box(ax, X_APPROACH_L, X_APPROACH_R, ROWS["baysor"], BOX_H,
        "Transcripts → Baysor EM",
        FAMILY_COLORS["baysor"], bold=True, fontsize=10.5,
        sublabel="transcript-density EM · ~99% capture · 4 tiles")

    box(ax, X_APPROACH_L, X_APPROACH_R, ROWS["ref"], BOX_H,
        "10x native (Xenium Ranger output)",
        FAMILY_COLORS["ref"], bold=True, fontsize=10.5,
        sublabel="provided whole-cell segmentation · reference anchor")

    # -------- arrows from approach → cell×gene --------
    for row_key in ("nuclear", "voronoi", "baysor", "ref"):
        y = ROWS[row_key]
        brace_arrow(ax, X_APPROACH_R + 0.05, X_MATRIX_L - 0.05,
                    y, y, color="#777777")

    # -------- cell×gene matrix box (shared) --------
    cy = (ROWS["nuclear"] + ROWS["ref"]) / 2 + 0.1
    cx_mid = (X_MATRIX_L + X_MATRIX_R) / 2
    cell_h = ROWS["nuclear"] - ROWS["ref"] + BOX_H * 0.5
    rect = FancyBboxPatch(
        (X_MATRIX_L, cy - cell_h / 2), X_MATRIX_R - X_MATRIX_L, cell_h,
        boxstyle="round,pad=0.1",
        facecolor="#f0f4f8", edgecolor="#aaaaaa",
        linewidth=1.5, zorder=3,
    )
    ax.add_patch(rect)
    ax.text(cx_mid, cy + 0.18, "Cell × gene matrix", ha="center", va="center",
            fontsize=11, fontweight="bold", color="#333333", zorder=4)
    ax.text(cx_mid, cy - 0.22, "(per method)", ha="center", va="center",
            fontsize=9, color="#666666", style="italic", zorder=4)

    # -------- arrow from cell×gene → eval pipeline --------
    arrow(ax, X_MATRIX_R + 0.05, X_EVAL_L - 0.05, cy, color="#555555", lw=2)

    # -------- evaluation pipeline (horizontal chain) --------
    n_steps = len(EVAL_STEPS)
    step_w = (X_EVAL_R - X_EVAL_L) / n_steps
    eval_h = 0.72
    for i, (step_label, fill) in enumerate(EVAL_STEPS):
        xl = X_EVAL_L + i * step_w + 0.06
        xr = X_EVAL_L + (i + 1) * step_w - 0.06
        is_last = i == n_steps - 1
        edge_color = "#27ae60" if is_last else "#aaaaaa"
        text_color = "#1a5c35" if is_last else "#333333"
        rect = FancyBboxPatch(
            (xl, cy - eval_h / 2), xr - xl, eval_h,
            boxstyle="round,pad=0.08",
            facecolor=fill, edgecolor=edge_color,
            linewidth=1.5 if is_last else 1.0, zorder=3,
        )
        ax.add_patch(rect)
        ax.text((xl + xr) / 2, cy, step_label,
                ha="center", va="center", fontsize=9,
                color=text_color, fontweight="bold" if is_last else "normal", zorder=4)
        if i < n_steps - 1:
            arrow(ax, xr + 0.01, X_EVAL_L + (i + 1) * step_w + 0.05,
                  cy, color="#888888", lw=1.2)

    # -------- column headers --------
    hdr_y = YT - 0.18
    ax.text((X_APPROACH_L + X_APPROACH_R) / 2, hdr_y,
            "Segmentation", ha="center", va="top",
            fontsize=11, color="#444444", fontweight="bold")
    ax.text((X_MATRIX_L + X_MATRIX_R) / 2, hdr_y,
            "Quantification", ha="center", va="top",
            fontsize=11, color="#444444", fontweight="bold")
    ax.text((X_EVAL_L + X_EVAL_R) / 2, hdr_y,
            "Evaluation (shared)", ha="center", va="top",
            fontsize=11, color="#444444", fontweight="bold")

    # column header underlines
    for xl, xr in [(X_APPROACH_L, X_APPROACH_R),
                   (X_MATRIX_L, X_MATRIX_R),
                   (X_EVAL_L, X_EVAL_R)]:
        ax.plot([xl, xr], [hdr_y - 0.08, hdr_y - 0.08],
                color="#cccccc", lw=1, zorder=2)

    fig.tight_layout(pad=0.3)
    fig.savefig(FIGURES / "pipeline_diagram.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("Saved pipeline_diagram.png")


if __name__ == "__main__":
    main()
