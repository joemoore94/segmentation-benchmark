"""Generate a comprehensive PDF report of all segmentation benchmark results.

Includes narrative prose, embedded figures, and data tables. Pages are sized
dynamically to fit their content.

Usage::

    conda run -n segbench python scripts/make_pdf_report.py
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from segbench.constants import (
    CLUSTER_ANNOTATIONS,
    METHOD_FAMILIES,
    METHOD_LABELS,
    NEGATIVE_PAIRS_TIER1,
    NEGATIVE_PAIRS_TIER2,
    TOTAL_TRANSCRIPTS_FULL_ROI,
)

TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")
OUTPUT_PDF = Path("results/segmentation_benchmark_complete_results.pdf")

PAGE_W = 8.5

FAMILY_COLORS = {
    "Reference":          "#d4edda",
    "Nuclear":            "#cce5ff",
    "Voronoi":            "#fff3cd",
    "Transcript-density": "#f8d7da",
}

METHOD_ORDER = [
    "10x_native",
    "cellpose", "stardist", "mesmer", "10x_ranger",
    "voronoi", "voronoi_stardist", "voronoi_mesmer", "voronoi_10x_ranger",
    "baysor", "baysor_prior", "baysor_prior_c08", "baysor_prior_c10",
    "baysor_stardist_prior_c10", "baysor_mesmer_prior_c10",
    "baysor_10x_ranger_prior_c10",
]


def label(key: str) -> str:
    return METHOD_LABELS.get(key, key)


def family(key: str) -> str:
    return METHOD_FAMILIES.get(key, "")


def sort_methods(df: pd.DataFrame, col: str = "method") -> pd.DataFrame:
    order = {k: i for i, k in enumerate(METHOD_ORDER)}
    df = df.copy()
    df["_sort"] = df[col].map(order).fillna(99)
    return df.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)


def method_row_colors(methods: list[str]) -> list[str]:
    return [FAMILY_COLORS.get(family(m), "white") for m in methods]


# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------

def fmt_pct(v):
    try:
        return f"{float(v) * 100:.1f}%"
    except (ValueError, TypeError):
        return str(v)

def fmt_pct_raw(v):
    try:
        return f"{float(v):.1f}%"
    except (ValueError, TypeError):
        return str(v)

def fmt_corr(v):
    try:
        return f"{float(v):.3f}"
    except (ValueError, TypeError):
        return str(v)

def fmt_int(v):
    try:
        return f"{int(float(v)):,}"
    except (ValueError, TypeError):
        return str(v)

def fmt_sci(v):
    try:
        f = float(v)
        if f < 0.001:
            return f"{f:.2e}"
        return f"{f:.4f}"
    except (ValueError, TypeError):
        return str(v)

def fmt_float2(v):
    try:
        return f"{float(v):.2f}"
    except (ValueError, TypeError):
        return str(v)


# ---------------------------------------------------------------------------
# Rendering helpers
# ---------------------------------------------------------------------------

def render_title_page(pdf):
    fig, ax = plt.subplots(figsize=(PAGE_W, 11))
    ax.axis("off")
    ax.text(0.5, 0.70,
            "Segmentation Benchmarking on\nXenium Spatial Transcriptomics",
            ha="center", va="center", fontsize=28, fontweight="bold",
            transform=ax.transAxes)
    ax.text(0.5, 0.55,
            "Complete Results",
            ha="center", va="center", fontsize=20, fontstyle="italic",
            transform=ax.transAxes)
    ax.text(0.5, 0.40,
            "Dataset: Xenium FFPE Human Breast (Janesick et al. 2023)\n"
            "ROI: 2 mm x 2 mm  |  ~23,600 cells  |  ~3.4 M transcripts  |  380-gene panel",
            ha="center", va="center", fontsize=12, transform=ax.transAxes)
    ax.text(0.5, 0.28,
            "17 segmentation methods: 4 nuclear detectors x 2 expansion strategies\n"
            "+ Baysor prior confidence sweep + reference-free negative marker validation",
            ha="center", va="center", fontsize=11, color="#555555",
            transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_toc(pdf):
    fig, ax = plt.subplots(figsize=(PAGE_W, 11))
    ax.axis("off")
    ax.set_title("Contents", fontsize=22, fontweight="bold", pad=20, loc="left")
    toc = [
        ("1", "Dataset & Methods Overview",
         "Dataset summary, methods inventory, cluster annotations"),
        ("2", "Cell & Transcript Recovery",
         "Cell counts, transcript capture, size distributions"),
        ("3", "Clustering Comparison",
         "Resolution sensitivity, cluster alignment, expression correlation, pseudobulk"),
        ("4", "Spatial Structure of Disagreement",
         "Global Moran's I, LISA hotspot/coldspot maps"),
        ("5", "Cell-Type Sensitivity",
         "Per-cell-type disagreement rates"),
        ("6", "Disagreement Drivers: Cell State vs Geometry",
         "Phenotypic density, cell size, differential expression"),
        ("7", "Phenotypic Landscape Distortion",
         "Shared PCA/UMAP manifold, density ratio maps"),
        ("8", "Pairwise Method Agreement",
         "Inter-method ARI matrix"),
        ("9", "Marker Gene Recovery",
         "Cell-type marker detection relative to 10x native"),
        ("10", "Population-Level Convergence",
         "Per-cell-type and per-cluster pseudobulk correlation"),
        ("11", "Negative Marker Analysis",
         "Cross-lineage contamination, violation rates"),
        ("12", "Confusion Matrices",
         "15 cluster-level confusion matrices"),
        ("13", "Gaps & Coverage Matrix",
         "Analysis coverage across methods, known data issues"),
    ]
    y = 0.90
    for num, title, desc in toc:
        ax.text(0.05, y, f"Section {num}", fontsize=9, color="#888888",
                fontweight="bold", transform=ax.transAxes, va="center")
        ax.text(0.18, y, title, fontsize=11, fontweight="bold",
                transform=ax.transAxes, va="center")
        ax.text(0.18, y - 0.022, desc, fontsize=8, color="#666666",
                transform=ax.transAxes, va="center")
        y -= 0.065
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_section_title(pdf, num, title):
    fig, ax = plt.subplots(figsize=(PAGE_W, 4))
    ax.axis("off")
    ax.text(0.5, 0.65, f"Section {num}",
            ha="center", va="center", fontsize=16, color="#888888",
            transform=ax.transAxes)
    ax.text(0.5, 0.30, title,
            ha="center", va="center", fontsize=24, fontweight="bold",
            transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_text(pdf, text: str, title: str | None = None):
    wrapped = textwrap.fill(text, width=95)
    n_lines = wrapped.count("\n") + 1
    height = max(3, 1.5 + n_lines * 0.22)
    fig, ax = plt.subplots(figsize=(PAGE_W, height))
    ax.axis("off")
    y_top = 0.95
    if title:
        ax.text(0.05, y_top, title, fontsize=14, fontweight="bold",
                transform=ax.transAxes, va="top")
        y_top -= 0.08
    ax.text(0.05, y_top, wrapped, fontsize=9.5, transform=ax.transAxes,
            va="top", ha="left", linespacing=1.5,
            fontfamily="serif")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_figure(pdf, path: Path, caption: str | None = None):
    if not path.exists():
        return
    img = mpimg.imread(str(path))
    h, w = img.shape[:2]
    aspect = h / w
    fig_w = PAGE_W - 0.5
    fig_h = fig_w * aspect
    total_h = fig_h + (0.6 if caption else 0.1)
    fig = plt.figure(figsize=(PAGE_W, total_h))
    ax = fig.add_axes([0.025, (0.5 if caption else 0.05) / total_h,
                       0.95, fig_h / total_h])
    ax.imshow(img)
    ax.axis("off")
    if caption:
        fig.text(0.5, 0.02, caption, ha="center", fontsize=8, color="#555555",
                 style="italic", wrap=True)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_table(pdf, df: pd.DataFrame, title: str,
                 footnote: str | None = None,
                 font_size: float = 8,
                 row_colors: list[str] | None = None,
                 landscape: bool = False):
    n_rows, n_cols = df.shape
    row_h = 0.30
    header_h = 0.35
    title_h = 0.6
    footnote_h = 0.5 if footnote else 0.1
    table_h = header_h + n_rows * row_h
    total_h = title_h + table_h + footnote_h

    if landscape:
        fig_w, fig_h = 11, max(5, total_h)
    else:
        fig_w = PAGE_W
        fig_h = max(3, total_h)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=15, loc="left")

    tbl = ax.table(
        cellText=df.values.tolist(),
        colLabels=list(df.columns),
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_facecolor("#4472C4")
            cell.set_text_props(color="white", fontweight="bold", fontsize=font_size)
        else:
            if row_colors and r - 1 < len(row_colors):
                cell.set_facecolor(row_colors[r - 1])
            else:
                cell.set_facecolor("#f9f9f9" if r % 2 == 0 else "white")
        if c == 0:
            cell.set_text_props(ha="left")

    tbl.auto_set_column_width(list(range(n_cols)))

    if footnote:
        fig.text(0.05, 0.01, footnote, fontsize=7, color="#666666",
                 va="bottom", wrap=True)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def render_gap_note(pdf, text: str):
    fig, ax = plt.subplots(figsize=(PAGE_W, 2.5))
    ax.axis("off")
    ax.text(0.5, 0.5, f"Gap Note: {text}",
            ha="center", va="center", fontsize=10,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#fff3cd",
                      edgecolor="#ffc107"))
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_spatial_jsons(pattern: str) -> pd.DataFrame:
    rows = []
    for f in sorted(TABLES_DIR.glob(pattern)):
        with open(f) as fh:
            d = json.load(fh)
        name = f.stem.replace("disagreement_spatial_", "")
        rows.append({
            "comparison": name,
            "n_cells": d.get("n_cells", ""),
            "disagreement_rate": d.get("disagreement_rate", ""),
            "morans_i": d.get("morans_i", ""),
            "p_value": d.get("p_value", ""),
        })
    return pd.DataFrame(rows)


def compute_expression_stats(filepath: Path) -> dict:
    df = pd.read_csv(filepath)
    corrs = df["correlation"].dropna()
    return {
        "n_pairs": len(corrs),
        "median_r": corrs.median(),
        "p25_r": corrs.quantile(0.25),
        "p75_r": corrs.quantile(0.75),
    }


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def section_1(pdf):
    render_section_title(pdf, 1, "Dataset & Methods Overview")

    render_text(pdf,
        "Xenium FFPE Human Breast (Custom Add-on Panel), Janesick et al. 2023, Nature Communications. "
        "Invasive ductal carcinoma; matched scRNA-seq + Visium from the same tissue blocks: GEO GSE243275. "
        "All analysis runs on a 2 mm x 2 mm ROI (~23,600 cells, ~3.4 M transcripts, 380-gene panel) "
        "with a mix of tumor, stroma, and immune-infiltrated regions.",
        title="Dataset")

    ds = pd.DataFrame([
        {"Property": "Tissue", "Value": "FFPE Human Breast, invasive ductal carcinoma"},
        {"Property": "Source", "Value": "Janesick et al. 2023, Nature Communications"},
        {"Property": "Platform", "Value": "10x Xenium (Custom Add-on Panel)"},
        {"Property": "ROI Dimensions", "Value": "2 mm x 2 mm"},
        {"Property": "Total Cells (10x native)", "Value": "23,629"},
        {"Property": "Total Transcripts (ROI)", "Value": f"{TOTAL_TRANSCRIPTS_FULL_ROI:,}"},
        {"Property": "Gene Panel", "Value": "380 genes"},
        {"Property": "Companion scRNA-seq", "Value": "GEO GSE243275"},
    ])
    render_table(pdf, ds, "Table 1.1 - Dataset Summary", font_size=10)

    render_text(pdf,
        "Four nuclear detectors (CellPose, StarDist, Mesmer, 10x Ranger) are tested with two expansion "
        "strategies (Voronoi nearest-centroid and Baysor transcript-density EM), producing a 4-detector x "
        "2-expansion comparison that cleanly separates nuclear detection quality from expansion strategy. "
        "10x native (Xenium Ranger's full segmentation with proprietary expansion) serves as the reference "
        "anchor. Baysor is additionally tested across a prior_segmentation_confidence sweep (0.0, 0.2, 0.8, "
        "1.0) with CellPose prior, and at PSC=1.0 with all four nuclear detectors. At PSC=1.0, nuclear "
        "transcripts are hard-locked to the detector's assignments and only cytoplasmic transcripts use "
        "density-adaptive expansion.",
        title="Methods")

    methods_data = []
    method_details = {
        "10x_native":       ("10x Ranger", "Proprietary", "-"),
        "cellpose":         ("CellPose 3.x", "None (nuclear only)", "-"),
        "stardist":         ("StarDist 2D_versatile_fluo", "None (nuclear only)", "-"),
        "mesmer":           ("DeepCell/Mesmer (Docker)", "None (nuclear only)", "-"),
        "10x_ranger":       ("10x Ranger (nuclei)", "None (nuclear only)", "-"),
        "voronoi":          ("CellPose", "Voronoi (nearest centroid)", "-"),
        "voronoi_stardist": ("StarDist", "Voronoi (nearest centroid)", "-"),
        "voronoi_mesmer":   ("Mesmer", "Voronoi (nearest centroid)", "-"),
        "voronoi_10x_ranger": ("10x Ranger", "Voronoi (nearest centroid)", "-"),
        "baysor":           ("-", "Baysor transcript-density EM", "0.0"),
        "baysor_prior":     ("CellPose", "Baysor + nuclear prior", "0.2"),
        "baysor_prior_c08": ("CellPose", "Baysor + nuclear prior", "0.8"),
        "baysor_prior_c10": ("CellPose", "Baysor + nuclear prior", "1.0"),
        "baysor_stardist_prior_c10": ("StarDist", "Baysor + nuclear prior", "1.0"),
        "baysor_mesmer_prior_c10":   ("Mesmer", "Baysor + nuclear prior", "1.0"),
        "baysor_10x_ranger_prior_c10": ("10x Ranger", "Baysor + nuclear prior", "1.0"),
    }
    for key in METHOD_ORDER:
        det, exp, psc = method_details.get(key, ("-", "-", "-"))
        methods_data.append({
            "Method": label(key), "Family": family(key),
            "Detector": det, "Expansion": exp, "PSC": psc,
        })
    colors = method_row_colors(METHOD_ORDER)
    render_table(pdf, pd.DataFrame(methods_data),
                 "Table 1.2 - Methods Inventory", font_size=8, row_colors=colors)

    ann = pd.read_csv(TABLES_DIR / "cluster_annotations.csv")
    ann["leiden_cluster"] = ann["leiden_cluster"].astype(str)
    ann.columns = ["Leiden Cluster", "Cell Type", "n_cells"]
    ann["n_cells"] = ann["n_cells"].apply(fmt_int)
    render_table(pdf, ann, "Table 1.3 - Cluster Annotations (10x Native Reference)",
                 font_size=9)


def section_2(pdf):
    render_section_title(pdf, 2, "Cell & Transcript Recovery")

    # Nuclear detectors
    render_text(pdf,
        "All four detectors operate on the same DAPI image but produce substantially different masks. "
        "Mesmer detects the largest nuclei (median ~45 um^2, long tail to 200 um^2), capturing 51.8% of "
        "transcripts - nearly double CellPose's 35.4%. StarDist and 10x Ranger produce similar-sized masks "
        "but StarDist finds more nuclei (24,745 vs 23,624). 10x Ranger captures only 38% of transcripts "
        "despite detecting nearly as many cells as the 10x native whole-cell segmentation (23,624 vs "
        "23,629), confirming that 10x native's 99% capture comes from its proprietary expansion, not "
        "from larger nuclei.",
        title="Nuclear Detectors")

    render_figure(pdf, FIGURES_DIR / "nuclear_mask_sizes.png",
                  "Nuclear mask size comparison across all four detectors")

    # Expansion methods
    render_text(pdf,
        "Voronoi expansion captures 100% of transcripts by construction regardless of detector. Median "
        "transcripts per cell under Voronoi varies with detector quality: Voronoi (CP) leads at 149 "
        "tx/cell because CellPose detects fewer, larger nuclei, concentrating more transcripts per cell. "
        "Voronoi (10x) at 128 tx/cell is closest to 10x native (124), consistent with using the same "
        "nuclear seeds. Baysor without a prior captures 98.6% but detects fewer cells (18,321) - the "
        "density model merges adjacent cells freely. At PSC 0.2, the prior barely changes behavior "
        "(19,061 cells, 53 tx/cell). At PSC 0.8-1.0, the hard-locked nuclear seeds prevent merging and "
        "cell count jumps to ~30,000+. The four PSC=1.0 variants reveal how detector choice propagates "
        "through density-adaptive expansion: Baysor (SD prior) produces the most cells (34,230) because "
        "StarDist detects the most nuclei, while Baysor (M prior) has the highest median tx/cell (74) "
        "because Mesmer's larger nuclei anchor more cytoplasmic transcripts per cell.",
        title="Expansion Methods")

    render_figure(pdf, FIGURES_DIR / "transcripts_per_cell.png",
                  "Transcripts/cell distribution across expansion methods")

    # Table: Cell counts
    cc = pd.read_csv(TABLES_DIR / "cell_counts.csv")
    cc = sort_methods(cc)
    colors = method_row_colors(cc["method"].tolist())
    display = pd.DataFrame({
        "Method": cc["method"].map(label), "Family": cc["method"].map(family),
        "n_cells": cc["n_cells"].apply(fmt_int),
        "Total Transcripts": cc["total_transcripts"].apply(fmt_int),
        "Median tx/cell": cc["median_transcripts_per_cell"].apply(fmt_int),
        "Mean tx/cell": cc["mean_transcripts_per_cell"].apply(fmt_float2),
        "Capture Rate": cc["transcript_capture_rate"].apply(fmt_pct),
    })
    render_table(pdf, display, "Table 2.1 - Cell Counts & Transcript Recovery",
                 font_size=8, row_colors=colors)

    # Size distributions
    ss = pd.read_csv(TABLES_DIR / "size_summary.csv")
    ss = sort_methods(ss)
    colors_ss = method_row_colors(ss["method"].tolist())
    display_ss = pd.DataFrame({
        "Method": ss["method"].map(label), "Size Metric": ss["size_col"],
        "P10": ss["p10"].apply(fmt_float2), "P25": ss["p25"].apply(fmt_float2),
        "Median": ss["median"].apply(fmt_float2), "P75": ss["p75"].apply(fmt_float2),
        "P90": ss["p90"].apply(fmt_float2), "Mean": ss["mean"].apply(fmt_float2),
        "Std": ss["std"].apply(fmt_float2),
    })
    render_table(pdf, display_ss, "Table 2.2 - Cell/Nuclear Size Distributions",
                 font_size=7.5, row_colors=colors_ss,
                 footnote="Size metric is 'area' (pixels^2) for morphological methods, "
                          "'n_transcripts' for Baysor/10x native.")


def section_3(pdf):
    render_section_title(pdf, 3, "Clustering Comparison")

    render_text(pdf,
        "Leiden clustering runs independently on each method's cells (normalize -> PCA -> neighbors -> "
        "Leiden at resolution 1.0). Cluster labels are aligned across methods before computing confusion "
        "matrices and disagreement, using two algorithms: Hungarian (one-to-one) and argmax (many-to-one).",
        title="Approach")

    # Resolution sensitivity
    render_text(pdf,
        "The method ordering is stable across Leiden resolutions 0.3-2.0 under both alignment algorithms. "
        "Voronoi (Mesmer) leads at most resolutions (0.3, 0.6, 0.8-1.2); at resolutions 0.5 and 0.7, "
        "Voronoi (StarDist) briefly takes the lead, and at 1.5+ StarDist's higher cell count gives it a "
        "durable advantage as finer clustering demands more cells per cluster. Baysor without a prior is "
        "consistently lowest. The Hungarian alignment forces unmatched clusters into poor pairings when "
        "cluster counts differ, inflating disagreement for methods that produce more clusters. The argmax "
        "alignment lets multiple clusters map to the same reference cluster, reducing this artifact.",
        title="Resolution Stability")

    rs = pd.read_csv(TABLES_DIR / "resolution_sensitivity.csv")

    # ARI is alignment-agnostic: deduplicate and emit as standalone table + CSV
    ari_standalone = (rs[rs["matcher"] == "hungarian"]
                      [["resolution", "method", "ari", "n_clusters_10x", "n_clusters_comp"]]
                      .copy()
                      .sort_values(["resolution", "method"])
                      .reset_index(drop=True))
    ari_csv_path = TABLES_DIR / "resolution_ari.csv"
    ari_standalone.to_csv(ari_csv_path, index=False)
    print(f"  Wrote {ari_csv_path}")

    # Standalone ARI plot
    fig_ari, ax_ari = plt.subplots(figsize=(10, 6))
    for method_name in ari_standalone["method"].unique():
        mdata = ari_standalone[ari_standalone["method"] == method_name]
        ax_ari.plot(mdata["resolution"], mdata["ari"], "o-", label=method_name, markersize=4)
    ax_ari.set_xlabel("Leiden resolution")
    ax_ari.set_ylabel("ARI vs 10x native")
    ax_ari.set_title("ARI across Leiden resolutions (alignment-agnostic)")
    ax_ari.legend(fontsize=7, ncol=2, loc="lower left")
    ax_ari.grid(True, alpha=0.3)
    ari_plot_path = FIGURES_DIR / "resolution_ari.png"
    fig_ari.savefig(ari_plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig_ari)
    print(f"  Wrote {ari_plot_path}")

    render_figure(pdf, ari_plot_path,
                  "ARI vs Leiden resolution (alignment-agnostic)")

    ari_display = pd.DataFrame({
        "Res.": ari_standalone["resolution"],
        "Method": ari_standalone["method"],
        "ARI": ari_standalone["ari"].apply(fmt_corr),
        "Clusters (10x)": ari_standalone["n_clusters_10x"].astype(int),
        "Clusters (method)": ari_standalone["n_clusters_comp"].astype(int),
    })
    render_table(pdf, ari_display,
                 "Table 3.1a - ARI vs Resolution (alignment-agnostic)",
                 font_size=6.5,
                 footnote="ARI is partition-based and does not depend on cluster alignment. "
                          "Saved to results/tables/resolution_ari.csv.")

    # Generate 2-panel plots (disagreement + Moran's I only, no ARI)
    from segbench.constants import METHOD_COLORS as _MC
    resolutions = sorted(rs["resolution"].unique())
    tick_labels = {r: str(r) for r in [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]}

    for matcher_name, matcher_val in [("Hungarian", "hungarian"), ("Argmax", "argmax")]:
        sub = rs[rs["matcher"] == matcher_val]
        matcher_label = ("Hungarian (one-to-one)" if matcher_val == "hungarian"
                         else "Argmax (many-to-one)")
        fig2, (ax_d, ax_m) = plt.subplots(1, 2, figsize=(20, 9))
        methods_in_data = sub["method"].unique()
        label_to_key = {v: k for k, v in METHOD_LABELS.items()}
        for method_name in methods_in_data:
            s = sub[sub["method"] == method_name]
            mkey = label_to_key.get(method_name, method_name)
            color = _MC.get(mkey, None)
            ax_d.plot(s["resolution"], s["disagree_pct"], "o-", color=color,
                      label=method_name, linewidth=2.5, markersize=8)
            ax_m.plot(s["resolution"], s["morans_i"], "o-", color=color,
                      label=method_name, linewidth=2.5, markersize=8)
        for ax in (ax_d, ax_m):
            ax.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.4)
            ax.set_xlabel("Leiden resolution")
            ax.set_xticks(resolutions)
            ax.set_xticklabels([tick_labels.get(r, "") for r in resolutions])
            ax.legend(fontsize=10)
        ax_d.set_ylabel("Disagreement (%)")
        ax_d.set_title("Cell-type disagreement across resolutions", fontweight="bold")
        ax_m.set_ylabel("Global Moran's I of disagreement")
        ax_m.set_title("Spatial structure of disagreement across resolutions", fontweight="bold")
        fig2.suptitle(f"Leiden resolution sensitivity - {matcher_label} cluster alignment",
                      fontsize=18, fontweight="bold")
        fig2.tight_layout()
        out_path = FIGURES_DIR / f"resolution_disagree_morans_{matcher_val}.png"
        fig2.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig2)
        print(f"  Wrote {out_path}")

    render_figure(pdf, FIGURES_DIR / "resolution_disagree_morans_hungarian.png",
                  "Disagreement and Moran's I across Leiden resolutions - Hungarian alignment")
    render_figure(pdf, FIGURES_DIR / "resolution_disagree_morans_argmax.png",
                  "Disagreement and Moran's I across Leiden resolutions - argmax alignment")

    for matcher_name, matcher_val in [("Hungarian", "hungarian"), ("Argmax", "argmax")]:
        subset = rs[rs["matcher"] == matcher_val].copy()
        subset = subset.sort_values(["resolution", "method"]).reset_index(drop=True)

        n_10x = subset["n_clusters_10x"].astype(int)
        n_comp = subset["n_clusters_comp"].astype(int)

        if matcher_val == "hungarian":
            matched_col = np.minimum(n_10x, n_comp)
            unmatched_col = np.abs(n_10x - n_comp)
            display = pd.DataFrame({
                "Res.": subset["resolution"], "Method": subset["method"],
                "Clust. (10x)": n_10x, "Clust. (method)": n_comp,
                "Matched": matched_col, "Unmatched": unmatched_col,
                "Disagree %": subset["disagree_pct"].apply(fmt_pct_raw),
                "Moran's I": subset["morans_i"].apply(fmt_corr),
            })
            fn = ("Matched = min(n_10x, n_method). Unmatched clusters are forced "
                  "into poor pairings, inflating disagreement.")
        else:
            # argmax: every method cluster maps to its best 10x cluster (all mapped),
            # but some 10x clusters may receive no mapping
            display = pd.DataFrame({
                "Res.": subset["resolution"], "Method": subset["method"],
                "Clust. (10x)": n_10x, "Clust. (method)": n_comp,
                "10x receiving": np.minimum(n_10x, n_comp),
                "Method mapped": n_comp,
                "Disagree %": subset["disagree_pct"].apply(fmt_pct_raw),
                "Moran's I": subset["morans_i"].apply(fmt_corr),
            })
            fn = ("All method clusters map to their best-overlap 10x cluster. "
                  "'10x receiving' = at most n_10x clusters get at least one mapping; "
                  "'Method mapped' = all method clusters are assigned (= n_method).")

        render_table(pdf, display,
                     f"Table 3.1{'b' if matcher_val == 'hungarian' else 'c'} - "
                     f"Resolution Sensitivity: Disagreement & Moran's I ({matcher_name})",
                     font_size=6.5, footnote=fn)

    # Cluster comparison
    render_text(pdf,
        "10x native and Voronoi methods converge on 14-16 clusters with median sizes above 1,000 cells. "
        "Baysor without a prior and at PSC 0.2 produce 21-24 smaller clusters, consistent with "
        "over-segmentation. At PSC 0.8-1.0, Baysor prior variants produce 20-23 clusters with higher cell "
        "counts (29,000-34,000) because the hard-locked nuclear seeds prevent merging; their median cluster "
        "sizes approach the Voronoi range.",
        title="Cluster Distributions")

    render_figure(pdf, FIGURES_DIR / "cluster_comparison.png",
                  "Leiden clustering comparison across methods")

    # Expression correlation
    render_text(pdf,
        "Per-cell expression correlation is high for all methods (median 0.79-0.96), but cluster-label "
        "agreement tells a different story. Voronoi methods disagree with 10x native on 19-32% of matched "
        "cells. Baysor without a prior disagrees on 52% (Hungarian) / 44% (argmax); at PSC 0.2, "
        "disagreement is unchanged (52% / 39%), but PSC 1.0 variants drop to 32-38% (Hungarian) / "
        "31-33% (argmax), approaching the Voronoi range.",
        title="Expression Correlation vs Cluster Agreement")

    render_figure(pdf, FIGURES_DIR / "expression_correlation.png",
                  "Per-cell-pair expression correlation")

    # Master agreement table
    hungarian_jsons = load_spatial_jsons("disagreement_spatial_10x_*.json")
    hungarian_jsons = hungarian_jsons[~hungarian_jsons["comparison"].str.contains("argmax")]
    argmax_jsons = load_spatial_jsons("disagreement_spatial_10x_*_argmax.json")
    ari_df = pd.read_csv(TABLES_DIR / "pairwise_ari.csv")
    ari_10x = ari_df[ari_df["method_a"] == "10x native"].copy()

    rows = []
    for _, hrow in hungarian_jsons.iterrows():
        comp = hrow["comparison"]
        method_key = comp.replace("10x_", "", 1) if comp.startswith("10x_") else comp
        if method_key == comp:
            continue
        arow = argmax_jsons[argmax_jsons["comparison"].str.replace("_argmax", "") == comp]
        ari_match = ari_10x[ari_10x["method_b"] == label(method_key)]
        ec_file = TABLES_DIR / f"expression_correlation_10x_{method_key}.csv"
        med_r = ""
        if ec_file.exists():
            stats = compute_expression_stats(ec_file)
            med_r = fmt_corr(stats["median_r"])
        rows.append({
            "Method": label(method_key), "method_key": method_key,
            "Matched": fmt_int(hrow["n_cells"]), "Med. r": med_r,
            "ARI": fmt_corr(ari_match["ari"].values[0]) if len(ari_match) else "",
            "H. Disagree": fmt_pct(hrow["disagreement_rate"]),
            "H. Moran's I": fmt_corr(hrow["morans_i"]),
            "A. Disagree": fmt_pct(arow["disagreement_rate"].values[0]) if len(arow) else "",
            "A. Moran's I": fmt_corr(arow["morans_i"].values[0]) if len(arow) else "",
        })
    if rows:
        master = pd.DataFrame(rows)
        order = {k: i for i, k in enumerate(METHOD_ORDER)}
        master["_sort"] = master["method_key"].map(order).fillna(99)
        master = master.sort_values("_sort").drop(columns=["_sort", "method_key"]).reset_index(drop=True)
        render_table(pdf, master,
                     "Table 3.2 - Clustering Agreement vs 10x Native (All Methods)",
                     font_size=7.5)

    # Per-cluster pseudobulk
    render_text(pdf,
        "To test whether cluster-level expression profiles agree, matched cells are grouped by 10x "
        "native's 15 Leiden clusters and pseudobulked per method. Nuclear methods drop to r = 0.86-0.87 "
        "on luminal epithelial clusters (0, 1, 3, 8) - the same populations driving single-cell "
        "disagreement - while Voronoi variants stay above 0.99 across all clusters. Baysor shows a "
        "comparable luminal dip plus reduced correlation on macrophage clusters (2, 7), consistent with "
        "transcript-density boundaries partitioning those populations differently.",
        title="Per-Cluster Pseudobulk")

    render_figure(pdf, FIGURES_DIR / "pseudobulk_by_cluster.png",
                  "Per-cluster pseudobulk correlation vs 10x native")

    pbc = pd.read_csv(TABLES_DIR / "pseudobulk_by_cluster.csv")
    pbc_display = pbc.copy()
    for c in pbc_display.columns[1:]:
        pbc_display[c] = pbc_display[c].apply(fmt_corr)
    render_table(pdf, pbc_display,
                 "Table 3.3 - Per-Cluster Pseudobulk Correlation", font_size=6.5,
                 landscape=True)

    # Confusion matrices figure
    render_text(pdf,
        "Each row is one 10x native cluster; columns are the comparison method's clusters. Red cells mark "
        "Hungarian (one-to-one) matched pairs, blue cells mark argmax (many-to-one) matches, and purple "
        "cells mark pairs selected by both algorithms. Voronoi methods produce clean matches under both "
        "algorithms. Baysor's 15x21 matrix shows the key difference: under Hungarian, 6 clusters are "
        "forced into empty pairings; under argmax, every column maps to the highest-overlap reference "
        "cluster with no wasted assignments.",
        title="Cluster Alignment")

    render_figure(pdf, FIGURES_DIR / "confusion_clusters.png",
                  "Confusion matrices with Hungarian and argmax alignment")
    if (FIGURES_DIR / "confusion_clusters_argmax.png").exists():
        render_figure(pdf, FIGURES_DIR / "confusion_clusters_argmax.png",
                      "Confusion matrices - argmax alignment")

    # UMAP examples
    render_figure(pdf, FIGURES_DIR / "umap" / "umap_baysor_hungarian.png",
                  "Baysor UMAP - Hungarian alignment")
    render_figure(pdf, FIGURES_DIR / "umap" / "umap_baysor_argmax.png",
                  "Baysor UMAP - argmax alignment")


def section_4(pdf):
    render_section_title(pdf, 4, "Spatial Structure of Disagreement")

    render_text(pdf,
        "Nuclear and Voronoi disagreements are spatially structured (Moran's I 0.076-0.215 under "
        "Hungarian), concentrated in luminal epithelial territory. Mesmer has the most agreement coldspots "
        "(32.5% LL); Voronoi (Mesmer) has the fewest disagreement hotspots (9.5% HH), consistent with "
        "residual errors being diffuse boundary noise. Under Hungarian alignment Baysor's near-zero "
        "Moran's I (0.033) reflects noise from forced cluster mismatches; under argmax alignment Moran's I "
        "increases to 0.079, revealing that Baysor's genuine disagreements are spatially structured - "
        "just less so than morphological methods.",
        title="Hungarian Alignment")

    render_figure(pdf, FIGURES_DIR / "disagreement_spatial_map.png",
                  "Disagreement mapped spatially - Hungarian")
    render_figure(pdf, FIGURES_DIR / "local_morans_map.png",
                  "LISA hotspot/coldspot maps - Hungarian")

    # Hungarian table
    hung = load_spatial_jsons("disagreement_spatial_10x_*.json")
    hung = hung[~hung["comparison"].str.contains("argmax")]
    hung_display = pd.DataFrame({
        "Comparison": hung["comparison"].str.replace("10x_", "10x vs "),
        "n_cells": hung["n_cells"].apply(fmt_int),
        "Disagree Rate": hung["disagreement_rate"].apply(fmt_pct),
        "Moran's I": hung["morans_i"].apply(fmt_corr),
        "p_value": hung["p_value"].apply(fmt_sci),
    })
    render_table(pdf, hung_display,
                 "Table 4.1 - Global Spatial Autocorrelation (Hungarian)", font_size=8)

    # Argmax
    render_figure(pdf, FIGURES_DIR / "disagreement_spatial_map_argmax.png",
                  "Disagreement mapped spatially - argmax")

    argm = load_spatial_jsons("disagreement_spatial_10x_*_argmax.json")
    argm_display = pd.DataFrame({
        "Comparison": argm["comparison"].str.replace("_argmax", "").str.replace("10x_", "10x vs "),
        "n_cells": argm["n_cells"].apply(fmt_int),
        "Disagree Rate": argm["disagreement_rate"].apply(fmt_pct),
        "Moran's I": argm["morans_i"].apply(fmt_corr),
        "p_value": argm["p_value"].apply(fmt_sci),
    })
    render_table(pdf, argm_display,
                 "Table 4.2 - Global Spatial Autocorrelation (Argmax)", font_size=8)

    # LISA distribution
    lisa_rows = []
    for f in sorted(TABLES_DIR.glob("local_morans_10x_*.csv")):
        comp_name = f.stem.replace("local_morans_10x_", "")
        lm = pd.read_csv(f)
        n = len(lm)
        counts = lm["lisa_cluster"].value_counts()
        lisa_rows.append({
            "Comparison": f"10x vs {label(comp_name)}",
            "n_cells": fmt_int(n),
            "HH (hotspot)": f"{counts.get('HH', 0)} ({counts.get('HH', 0)/n*100:.1f}%)",
            "HL": f"{counts.get('HL', 0)} ({counts.get('HL', 0)/n*100:.1f}%)",
            "LH": f"{counts.get('LH', 0)} ({counts.get('LH', 0)/n*100:.1f}%)",
            "LL (coldspot)": f"{counts.get('LL', 0)} ({counts.get('LL', 0)/n*100:.1f}%)",
        })
    if lisa_rows:
        render_table(pdf, pd.DataFrame(lisa_rows),
                     "Table 4.3 - LISA Cluster Distribution", font_size=8,
                     footnote="HH = disagreement hotspot, LL = agreement coldspot. "
                              "Data for 8 of 17 comparisons.")


def section_5(pdf):
    render_section_title(pdf, 5, "Cell-Type Sensitivity")

    render_text(pdf,
        "Adipocytes and myoepithelial cells have the highest per-cell disagreement (~50-68% and ~40-47%) "
        "but are rare. Luminal epithelial cells dominate by volume: ~35% disagreement across ~8,500 cells "
        "drives the majority of total disagreement events. These clusters likely encompass malignant and "
        "normal epithelial cells; both share canonical markers (GATA3, PGR, ESR1, MUC1) and are "
        "inseparable by nuclear morphology alone. T cells and B cells are robustly identified regardless "
        "of method or alignment algorithm.")

    render_figure(pdf, FIGURES_DIR / "agreement_explainer.png",
                  "Cell type vs agreement - Hungarian")

    ct = pd.read_csv(TABLES_DIR / "celltype_disagreement.csv")
    pivot = ct.pivot(index="comparison", columns="cell_type", values="disagree_rate")
    pivot = pivot.map(fmt_pct)
    pivot.insert(0, "Method", pivot.index)
    pivot = pivot.reset_index(drop=True)
    render_table(pdf, pivot,
                 "Table 5.1 - Per-Cell-Type Disagreement Rate", font_size=7.5,
                 footnote="7 of 17 methods. Missing: 10x Ranger, Voronoi(10x), all Baysor prior variants.")

    # Full detail
    ct_display = pd.DataFrame({
        "Method": ct["comparison"], "Cell Type": ct["cell_type"],
        "n_matched": ct["n_matched"].apply(fmt_int),
        "n_disagree": ct["n_disagree"].apply(fmt_int),
        "Disagree Rate": ct["disagree_rate"].apply(fmt_pct),
    })
    render_table(pdf, ct_display,
                 "Table 5.2 - Cell-Type Disagreement (Full Detail)", font_size=7.5)


def section_6(pdf):
    render_section_title(pdf, 6, "Disagreement Drivers: Cell State vs Geometry")

    # Density
    render_text(pdf,
        "Nuclear methods disagree on cells in higher-density phenotypic regions (Mann-Whitney p << 0.001). "
        "The DE volcano confirms this: disagreeing cells are enriched for luminal epithelial markers "
        "(MYBPC1, SERPINA3, CLIC6, PGR, GATA3, MUC1), cytoplasmic transcripts underrepresented in "
        "nuclear-only masks. Voronoi (CellPose) disagreement is density-neutral (p = 0.19) with few DE "
        "genes, indicating residual errors are geometric. Baysor disagreement is also density-neutral but "
        "enriched for macrophage markers (CD14, MRC1, CD163), consistent with transcript-density "
        "boundaries partitioning macrophage-rich regions differently.",
        title="Phenotypic Density vs Disagreement")

    render_figure(pdf, FIGURES_DIR / "density_vs_disagreement.png",
                  "Phenotypic density vs disagreement - Hungarian")
    if (FIGURES_DIR / "density_vs_disagreement_argmax.png").exists():
        render_figure(pdf, FIGURES_DIR / "density_vs_disagreement_argmax.png",
                      "Phenotypic density vs disagreement - argmax")

    dd = pd.read_csv(TABLES_DIR / "density_disagreement_summary.csv")
    dd_display = pd.DataFrame({
        "Comparison": dd["comparison"],
        "n_agree": dd["n_agree"].apply(fmt_int),
        "n_disagree": dd["n_disagree"].apply(fmt_int),
        "Med. Density (agree)": dd["median_log_density_agree"].apply(fmt_float2),
        "Med. Density (disagree)": dd["median_log_density_disagree"].apply(fmt_float2),
        "p_value": dd["p_value"].apply(fmt_sci),
    })
    render_table(pdf, dd_display,
                 "Table 6.1 - Phenotypic Density vs Disagreement", font_size=8)

    # DE
    render_figure(pdf, FIGURES_DIR / "de_volcano.png",
                  "DE: disagree vs agree cells - Hungarian")

    de_files = sorted(TABLES_DIR.glob("de_disagree_10x_*.csv"))
    de_summary_rows = []
    for f in de_files:
        comp_name = f.stem.replace("de_disagree_", "")
        de = pd.read_csv(f)
        n_sig = (de["pvals_adj"] < 0.05).sum()
        top = de.iloc[0] if len(de) > 0 else None
        de_summary_rows.append({
            "Comparison": comp_name.replace("_", " vs ", 1),
            "n_sig (p<0.05)": n_sig,
            "Top Gene": top["names"] if top is not None else "",
            "Top log2FC": fmt_float2(top["logfoldchanges"]) if top is not None else "",
            "Top Score": fmt_float2(top["scores"]) if top is not None else "",
        })
    if de_summary_rows:
        render_table(pdf, pd.DataFrame(de_summary_rows),
                     "Table 6.2 - DE Summary (Disagree vs Agree Cells)", font_size=8)

    # Top 20 DE per comparison
    for f in de_files:
        comp_name = f.stem.replace("de_disagree_", "")
        de = pd.read_csv(f)
        top20 = de.head(20).copy()
        de_display = pd.DataFrame({
            "Gene": top20["names"],
            "Score": top20["scores"].apply(fmt_float2),
            "log2FC": top20["logfoldchanges"].apply(fmt_float2),
            "adj p-value": top20["pvals_adj"].apply(fmt_sci),
            "% non-zero": top20["pct_nz_group"].apply(fmt_pct),
        })
        render_table(pdf, de_display,
                     f"Table 6.3 - Top 20 DE Genes: {comp_name.replace('_', ' vs ', 1)}",
                     font_size=8)

    # Cell size
    render_text(pdf,
        "Smaller 10x-native cells are significantly more likely to disagree with every morphological "
        "method (p << 0.001). The direction is counter-intuitive: larger cells have more cytoplasm, yet "
        "it is smaller cells that disagree more. The pattern holds for Voronoi methods too, ruling out "
        "transcript capture as the cause. Smaller cells likely correspond to densely packed regions where "
        "any method's cluster assignment is noisier. Baysor shows no size dependence (p = 0.28); its "
        "boundaries are insensitive to morphologically defined cell area.",
        title="Cell Size vs Disagreement")

    render_figure(pdf, FIGURES_DIR / "cell_size_disagreement.png",
                  "Cell size vs disagreement - Hungarian")

    size_data = pd.DataFrame([
        {"Comparison": "10x vs CellPose",     "Med. Area (agree)": "123.9 um^2", "Med. Area (disagree)": "121.4 um^2", "p": "7.2e-07"},
        {"Comparison": "10x vs StarDist",     "Med. Area (agree)": "121.2 um^2", "Med. Area (disagree)": "116.9 um^2", "p": "2.4e-12"},
        {"Comparison": "10x vs Mesmer",       "Med. Area (agree)": "126.2 um^2", "Med. Area (disagree)": "119.1 um^2", "p": "3.2e-07"},
        {"Comparison": "10x vs Voronoi (CP)", "Med. Area (agree)": "126.8 um^2", "Med. Area (disagree)": "111.5 um^2", "p": "4.9e-32"},
        {"Comparison": "10x vs Voronoi (SD)", "Med. Area (agree)": "123.8 um^2", "Med. Area (disagree)": "112.4 um^2", "p": "8.8e-29"},
        {"Comparison": "10x vs Voronoi (M)",  "Med. Area (agree)": "125.9 um^2", "Med. Area (disagree)": "117.7 um^2", "p": "3.7e-19"},
        {"Comparison": "10x vs Baysor",       "Med. Area (agree)": "167.0 um^2", "Med. Area (disagree)": "173.4 um^2", "p": "0.28 n.s."},
    ])
    render_table(pdf, size_data, "Table 6.4 - Cell Size vs Disagreement", font_size=9)


def section_7(pdf):
    render_section_title(pdf, 7, "Phenotypic Landscape Distortion")

    render_text(pdf,
        "All methods are projected into a shared PCA space fit on 10x native (30 PCs, 55% variance "
        "explained) and embedded in a joint UMAP. Density ratio maps (log2 method/10x) show which "
        "phenotypic regions each method enriches or depletes. Nuclear methods show depleted regions in "
        "high-density luminal epithelial areas, consistent with missed cytoplasmic transcripts pulling "
        "cells toward lower-expression PCA states. Voronoi methods track 10x native closely. Baysor "
        "shows enrichment in a distinct region corresponding to its finer resolution of macrophage and "
        "stromal subtypes.")

    render_figure(pdf, FIGURES_DIR / "manifold_shared_umap.png",
                  "All methods in shared PCA/UMAP space")
    render_figure(pdf, FIGURES_DIR / "manifold_distortion.png",
                  "Phenotypic landscape distortion vs 10x native")


def section_8(pdf):
    render_section_title(pdf, 8, "Pairwise Method Agreement")

    render_text(pdf,
        "CellPose and StarDist agree with each other at ARI 0.764 (higher than the Voronoi pair at "
        "0.661) because both are nuclear-morphology methods on the same DAPI image. Switching to Voronoi "
        "assignment lowers within-paradigm agreement because the two Voronoi variants use different "
        "centroids, shifting boundaries even where centroids are close. What Voronoi does raise is "
        "agreement with the 10x-native whole-cell reference (0.63-0.69): compatibility with the "
        "platform's own segmentation, not cross-method reproducibility. Baysor remains isolated from all "
        "morphological methods (ARI 0.30-0.46 regardless of partner).")

    render_figure(pdf, FIGURES_DIR / "pairwise_consensus.png",
                  "Pairwise ARI between all segmentation methods")

    ari_all = pd.read_csv(TABLES_DIR / "pairwise_ari.csv")
    ari_display = pd.DataFrame({
        "Method A": ari_all["method_a"], "Method B": ari_all["method_b"],
        "ARI": ari_all["ari"].apply(fmt_corr),
        "n_matched": ari_all["n_matched"].apply(fmt_int),
    })
    render_table(pdf, ari_display, "Table 8.1 - Pairwise ARI (All Pairs)", font_size=9)

    ari_mat = pd.read_csv(TABLES_DIR / "pairwise_ari_matrix.csv", index_col=0)
    ari_mat_display = ari_mat.map(fmt_corr)
    ari_mat_display.insert(0, "Method", ari_mat_display.index)
    ari_mat_display = ari_mat_display.reset_index(drop=True)
    render_table(pdf, ari_mat_display, "Table 8.2 - Pairwise ARI Matrix", font_size=7)


def section_9(pdf):
    render_section_title(pdf, 9, "Marker Gene Recovery")

    render_text(pdf,
        "Using 10x-native cell-type annotations as ground truth, nuclear methods recover 75-92% of "
        "cytoplasmic marker expression relative to 10x native, with the largest deficits for extranuclear "
        "markers like MUC1, SERPINA3, and LYZ. Voronoi methods recover near-100% across all cell types. "
        "Baysor recovers macrophage markers (LYZ, CD14) at or above 10x-native levels while showing "
        "slightly reduced T cell marker (CD3E) recovery.")

    render_figure(pdf, FIGURES_DIR / "marker_recovery.png",
                  "Marker gene recovery relative to 10x native")

    render_gap_note(pdf,
        "Marker recovery data exists only as a figure. No CSV output; data computed on-the-fly from "
        "AnnData objects in scripts/run_marker_recovery.py.")


def section_10(pdf):
    render_section_title(pdf, 10, "Population-Level Convergence")

    render_text(pdf,
        "Pseudobulk is computed within each of 10 annotated cell types (not as a whole-ROI sum), so the "
        "correlation tests whether each method's cell-type compartments recover the same expression "
        "programs as 10x native. Baysor's per-cell-type correlations range from 0.94 (plasma cells) to "
        "0.997 (CAFs), degrading predictably on rare populations with fewer cells. Despite its low "
        "single-cell ARI of 0.305, Baysor is competitive with nuclear methods at the cell-type level - "
        "its aggregate r of 0.999 sits above CellPose (0.970) and StarDist (0.975). Nuclear methods show "
        "reduced pseudobulk r (0.97-0.98) because missing cytoplasmic transcripts suppress marker signal "
        "systematically across all cells of a type. Voronoi methods achieve both high single-cell ARI and "
        "near-perfect pseudobulk agreement.")

    render_figure(pdf, FIGURES_DIR / "pseudobulk_correlation.png",
                  "Pseudobulk correlation vs 10x native")

    pb = pd.read_csv(TABLES_DIR / "pseudobulk_correlation.csv")
    cols_to_show = ["method", "ari", "global_r"]
    ct_cols = [c for c in pb.columns if c not in ("method", "ari", "global_r", "color")]
    cols_to_show.extend(ct_cols)
    pb_display = pb[cols_to_show].copy()
    pb_display.columns = ["Method", "ARI", "Global r"] + ct_cols
    for c in ["ARI", "Global r"] + ct_cols:
        pb_display[c] = pb_display[c].apply(fmt_corr)
    render_table(pdf, pb_display,
                 "Table 10.1 - Per-Cell-Type Pseudobulk Correlation", font_size=7,
                 landscape=True)


def section_11(pdf):
    render_section_title(pdf, 11, "Negative Marker Analysis")

    render_text(pdf,
        "Negative marker analysis provides a reference-free segmentation quality metric using biologically "
        "impossible co-expression (e.g., CD3E + GATA3 = T cell + luminal epithelial in the same cell). "
        "Each method is scored by violation rate without relying on 10x native as ground truth. 11 Tier 1 "
        "pairs and 2 Tier 2 pairs are defined from the 380-gene panel. Across both expansion strategies, "
        "Voronoi consistently outperforms Baysor PSC=1.0 on ARI (0.58-0.69 vs 0.50-0.53), but Baysor "
        "prior variants have lower negative marker violation rates (0.31 per 1000 tx vs 0.37-0.43 for "
        "Voronoi), indicating that density-adaptive boundaries produce fewer cross-lineage contamination "
        "artifacts even when they disagree with the 10x reference.")

    render_figure(pdf, FIGURES_DIR / "negative_marker.png",
                  "Negative marker violation rates")

    nm = pd.read_csv(TABLES_DIR / "negative_marker_summary.csv")
    nm_display = pd.DataFrame({
        "Method": nm["method"], "n_cells": nm["n_cells"].apply(fmt_int),
        "Med. tx/cell": nm["median_tx_per_cell"].apply(fmt_int),
        "Tier1 Viol. Rate": nm["tier1_violation_rate"].apply(fmt_pct),
        "Tier1 Violations": nm["tier1_violations"].apply(fmt_int),
        "Mean Contam.": nm["mean_contamination_score"].apply(fmt_sci),
        "Viol./1000 tx": nm["violations_per_1000tx"].apply(fmt_float2),
    })
    render_table(pdf, nm_display, "Table 11.1 - Negative Marker Summary", font_size=8,
                 footnote="10 of 17 methods. Missing: Voronoi(10x), Baysor(SD/M/10x prior 1.0).")

    # Per-pair pivoted
    nv = pd.read_csv(TABLES_DIR / "negative_marker_violations.csv")
    pivot = nv.pivot_table(index="pair", columns="method", values="violation_rate", aggfunc="first")
    pivot = pivot.map(lambda x: fmt_sci(x) if pd.notna(x) else "-")
    pivot.insert(0, "Marker Pair", pivot.index)
    pivot = pivot.reset_index(drop=True)
    render_table(pdf, pivot, "Table 11.2 - Per-Pair Violation Rates by Method",
                 font_size=6.5, landscape=True)

    # Pair definitions
    pairs_data = []
    for ga, gb, ta, tb in NEGATIVE_PAIRS_TIER1:
        pairs_data.append({"Gene A": ga, "Gene B": gb, "Type A": ta, "Type B": tb, "Tier": 1})
    for ga, gb, ta, tb in NEGATIVE_PAIRS_TIER2:
        pairs_data.append({"Gene A": ga, "Gene B": gb, "Type A": ta, "Type B": tb, "Tier": 2})
    render_table(pdf, pd.DataFrame(pairs_data),
                 "Table 11.3 - Negative Marker Pair Definitions", font_size=9)


def section_12(pdf):
    render_section_title(pdf, 12, "Confusion Matrices")

    render_text(pdf,
        "One confusion matrix per method comparison. Each row is a 10x native Leiden cluster (annotated "
        "with cell type); columns are the comparison method's clusters. These matrices show how cluster "
        "identity maps between methods at the single-cell level.")

    confusion_files = sorted(TABLES_DIR.glob("cell_type_confusion_10x_*.csv"))
    for i, f in enumerate(confusion_files, 1):
        method_key = f.stem.replace("cell_type_confusion_10x_", "")
        cm = pd.read_csv(f)
        row_labels = []
        for _, row in cm.iterrows():
            cid = str(row.iloc[0])
            ct = CLUSTER_ANNOTATIONS.get(cid, "?")
            row_labels.append(f"{cid}: {ct}")
        first_col = cm.columns[0]
        cm[first_col] = cm[first_col].astype(str)
        cm[first_col] = row_labels
        cm.columns = ["10x Cluster"] + [str(c) for c in cm.columns[1:]]
        render_table(pdf, cm,
                     f"Table 12.{i} - Confusion: 10x vs {label(method_key)}",
                     font_size=7, landscape=True)


def section_13(pdf):
    render_section_title(pdf, 13, "Gaps & Coverage Matrix")

    analyses = [
        "Cell Counts", "Expr. Corr.", "ARI", "Confusion",
        "Pseudobulk", "Res. Sens.", "CellType Dis.",
        "Density Dis.", "DE", "Neg. Marker", "LISA",
    ]

    ec_methods = {f.stem.replace("expression_correlation_10x_", "")
                  for f in TABLES_DIR.glob("expression_correlation_10x_*.csv")}
    cm_methods = {f.stem.replace("cell_type_confusion_10x_", "")
                  for f in TABLES_DIR.glob("cell_type_confusion_10x_*.csv")}
    rs_df = pd.read_csv(TABLES_DIR / "resolution_sensitivity.csv")
    label_to_key = {v: k for k, v in METHOD_LABELS.items()}
    rs_methods = {label_to_key.get(m, m) for m in rs_df["method"].unique()}
    pb = pd.read_csv(TABLES_DIR / "pseudobulk_correlation.csv")
    pb_keys = {label_to_key.get(m, m) for m in pb["method"].unique()}
    ct = pd.read_csv(TABLES_DIR / "celltype_disagreement.csv")
    ct_keys = {label_to_key.get(m, m) for m in ct["comparison"].unique()}
    dd = pd.read_csv(TABLES_DIR / "density_disagreement_summary.csv")
    dd_comps = set(dd["comparison"].unique())
    de_methods = {f.stem.replace("de_disagree_10x_", "")
                  for f in TABLES_DIR.glob("de_disagree_10x_*.csv")}
    nm = pd.read_csv(TABLES_DIR / "negative_marker_summary.csv")
    nm_methods = set(nm["method_key"].unique())
    lisa_methods = {f.stem.replace("local_morans_10x_", "")
                    for f in TABLES_DIR.glob("local_morans_10x_*.csv")}
    cc_methods = set(pd.read_csv(TABLES_DIR / "cell_counts.csv")["method"].unique())
    ari_df = pd.read_csv(TABLES_DIR / "pairwise_ari.csv")
    ari_methods = set()
    for _, row in ari_df.iterrows():
        if row["method_a"] == "10x native":
            ari_methods.add(label_to_key.get(row["method_b"], row["method_b"]))

    rows = []
    for key in METHOD_ORDER:
        ck = lambda s: "Y" if s else "-"
        lbl = label(key)
        dd_match = any(lbl in c for c in dd_comps)
        rows.append({
            "Method": lbl, "Family": family(key),
            analyses[0]: ck(key in cc_methods),
            analyses[1]: ck(key in ec_methods),
            analyses[2]: ck(key in ari_methods),
            analyses[3]: ck(key in cm_methods),
            analyses[4]: ck(key in pb_keys),
            analyses[5]: ck(key in rs_methods),
            analyses[6]: ck(key in ct_keys),
            analyses[7]: ck(dd_match),
            analyses[8]: ck(key in de_methods),
            analyses[9]: ck(key in nm_methods),
            analyses[10]: ck(key in lisa_methods),
        })

    colors = method_row_colors(METHOD_ORDER)
    render_table(pdf, pd.DataFrame(rows),
                 "Table 13.1 - Analysis Coverage Matrix",
                 font_size=7, row_colors=colors, landscape=True)

    issues = pd.DataFrame([
        {"Issue": "Voronoi(10x) median expression correlation",
         "Status": "Pending", "Detail": "Gene-name encoding mismatch; needs re-run"},
        {"Issue": "Baysor PSC sweep limited to CellPose prior",
         "Status": "By design", "Detail": "Only CellPose at PSC 0.2/0.8/1.0; others only 1.0"},
        {"Issue": "LISA limited to 8 comparisons",
         "Status": "Partial", "Detail": "Missing Baysor PSC 0.8/1.0, 10x Ranger, Voronoi(10x)"},
        {"Issue": "Cell-type disagreement: 7 methods",
         "Status": "Partial", "Detail": "Missing Baysor prior variants, 10x Ranger, Voronoi(10x)"},
        {"Issue": "DE analysis: 8 methods",
         "Status": "Partial", "Detail": "Missing Voronoi(10x), Baysor PSC 0.8/1.0, 10x Ranger"},
        {"Issue": "Negative marker: 10 methods",
         "Status": "Partial", "Detail": "Missing Voronoi(10x), Baysor(SD/M/10x prior 1.0)"},
        {"Issue": "Pseudobulk ARI=0.0 for Baysor(SD/M prior 1.0)",
         "Status": "Data issue", "Detail": "ARI not computed at generation time"},
        {"Issue": "Marker gene recovery",
         "Status": "Figure only", "Detail": "No CSV; data computed on-the-fly from AnnData"},
    ])
    render_table(pdf, issues, "Table 13.2 - Known Data Issues & Gaps", font_size=8)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
    })

    print(f"Generating PDF report -> {OUTPUT_PDF}")
    with PdfPages(str(OUTPUT_PDF)) as pdf:
        render_title_page(pdf)
        render_toc(pdf)
        section_1(pdf)
        section_2(pdf)
        section_3(pdf)
        section_4(pdf)
        section_5(pdf)
        section_6(pdf)
        section_7(pdf)
        section_8(pdf)
        section_9(pdf)
        section_10(pdf)
        section_11(pdf)
        section_12(pdf)
        section_13(pdf)

    print(f"Done. {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
