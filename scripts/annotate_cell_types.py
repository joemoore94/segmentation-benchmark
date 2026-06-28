"""Cell type annotation via Leiden clustering and differential expression.

Clusters the 10x native cells (Leiden, resolution 1.0), runs Wilcoxon DE to
find the top marker genes per cluster, and assigns cell type labels based on
canonical breast tissue markers. Produces an evidence table showing the DE
genes that justify each annotation.

Reads:  data/processed/roi/adata_10x.h5ad
Writes: results/tables/annotation_evidence.csv
        results/figures/annotation_dotplot.png
        results/figures/annotation_evidence_heatmap.png

Usage::

    conda run -n segbench python scripts/annotate_cell_types.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from segbench.constants import CLUSTER_ANNOTATIONS
from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
TABLES = Path("results/tables")
FIGURES = Path("results/figures")
DPI = 200

# Canonical markers used to justify each annotation. These are the genes a
# biologist would look for to confirm cell type identity in breast tissue.
# The DE results should surface these (or close relatives) as top hits.
CANONICAL_MARKERS: dict[str, list[str]] = {
    "Luminal epithelial": [
        "EPCAM", "KRT7", "KRT8", "GATA3", "ESR1", "PGR", "FOXA1",
        "MUC1", "ANKRD30A", "TACSTD2", "CCND1",
    ],
    "Myoepithelial": ["KRT14", "KRT5", "ACTA2", "MYLK", "DST"],
    "T cells": ["CD3E", "CD3G", "TRAC", "TRBC1", "CD96", "IL7R", "CCL5"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19"],
    "Plasma cells": ["MZB1", "SLAMF7", "TENT5C", "TNFRSF17"],
    "Macrophages": ["CD14", "CD68", "CD163", "AIF1", "LYZ", "FCER1G"],
    "CAFs": ["LUM", "SFRP4", "FBLN1", "CCDC80", "THBS2", "MMP2", "PDGFRA"],
    "Smooth muscle": ["MYH11", "ACTA2", "MYLK", "RGS5", "CAV1"],
    "Endothelial": ["PECAM1", "VWF", "AQP1", "CD93", "CLEC14A"],
    "Adipocytes": ["ADIPOQ", "PLIN1", "PPARG", "LPL", "G0S2"],
}


def load_and_cluster() -> ad.AnnData:
    adata = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    sc.settings.verbosity = 0
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph")
    return adata


def build_evidence_table(adata: ad.AnnData, n_top: int = 8) -> pd.DataFrame:
    """Run DE per cluster and build a table of top markers + annotation."""
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    rows = []
    for cluster in sorted(adata.obs["leiden"].unique(), key=int):
        n_cells = int((adata.obs["leiden"] == cluster).sum())
        de = sc.get.rank_genes_groups_df(adata, group=cluster)
        top = de[de["pvals_adj"] < 0.05].head(n_top)

        top_genes = top["names"].tolist()
        top_lfc = top["logfoldchanges"].tolist()

        annotation = CLUSTER_ANNOTATIONS.get(cluster, "Unknown")

        canonical = CANONICAL_MARKERS.get(annotation, [])
        canonical_in_top = [g for g in canonical if g in top_genes]
        canonical_in_de = [
            g for g in canonical
            if g in de[de["pvals_adj"] < 0.05]["names"].values
        ]

        rows.append({
            "cluster": cluster,
            "n_cells": n_cells,
            "annotation": annotation,
            "top_de_genes": ", ".join(
                f"{g} (+{lfc:.1f})" for g, lfc in zip(top_genes, top_lfc)
            ),
            "canonical_markers_in_top8": ", ".join(canonical_in_top),
            "canonical_markers_in_de": ", ".join(canonical_in_de),
            "n_canonical_in_de": len(canonical_in_de),
            "n_canonical_total": len(canonical),
        })

    return pd.DataFrame(rows)


def fig_marker_detection(adata: ad.AnnData) -> None:
    """Stacked bar chart: fraction of cells in each cluster expressing canonical markers."""
    apply_style()

    panel_genes = set(adata.var_names)

    # Pick top 2-3 most specific markers per cell type
    KEY_MARKERS: dict[str, list[str]] = {
        "Luminal epithelial": ["EPCAM", "GATA3", "ESR1"],
        "Myoepithelial": ["KRT14", "ACTA2"],
        "T cells": ["CD3E", "TRAC", "CCL5"],
        "B cells": ["MS4A1", "CD79A"],
        "Plasma cells": ["MZB1", "SLAMF7"],
        "Macrophages": ["CD14", "LYZ", "FCER1G"],
        "CAFs": ["LUM", "SFRP4", "PDGFRA"],
        "Smooth muscle": ["MYH11", "RGS5"],
        "Endothelial": ["PECAM1", "VWF"],
        "Adipocytes": ["ADIPOQ", "PLIN1"],
    }

    all_markers = []
    marker_to_ct = {}
    for ct, genes in KEY_MARKERS.items():
        for g in genes:
            if g in panel_genes:
                all_markers.append(g)
                marker_to_ct[g] = ct

    X = adata[:, all_markers].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    detected = pd.DataFrame(
        (X > 0).astype(float), columns=all_markers, index=adata.obs_names,
    )
    detected["cluster"] = adata.obs["leiden"].values

    pct = detected.groupby("cluster").mean() * 100
    pct = pct.loc[sorted(pct.index, key=int)]

    labels = [
        f"{c}: {CLUSTER_ANNOTATIONS.get(c, '?')}" for c in pct.index
    ]

    fig, ax = plt.subplots(figsize=(22, 12))
    sns.heatmap(
        pct, annot=True, fmt=".0f", cmap="YlOrRd",
        linewidths=0.5, ax=ax, vmin=0, vmax=100,
        yticklabels=labels, annot_kws={"size": 11},
        cbar_kws={"label": "% cells expressing marker"},
    )
    ax.set_title("Canonical marker detection rate by Leiden cluster (10x native)",
                 fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Add cell type group labels along the top
    prev_ct = None
    start = 0
    for i, g in enumerate(all_markers):
        ct = marker_to_ct[g]
        if ct != prev_ct:
            if prev_ct is not None:
                mid = (start + i) / 2
                ax.text(mid, -0.8, prev_ct, ha="center", va="bottom",
                        fontsize=11, fontweight="bold", rotation=30)
            start = i
            prev_ct = ct
    mid = (start + len(all_markers)) / 2
    ax.text(mid, -0.8, prev_ct, ha="center", va="bottom",
            fontsize=11, fontweight="bold", rotation=30)

    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(FIGURES / "annotation_dotplot.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_dotplot_canonical(adata: ad.AnnData) -> None:
    """Custom dotplot: markers on y-axis (grouped by cell type on the left),
    clusters on top x-axis (rotated labels)."""
    apply_style()

    panel_genes = set(adata.var_names)

    seen = set()
    marker_order = []
    marker_ct = []
    for ct, genes in CANONICAL_MARKERS.items():
        for g in genes:
            if g in panel_genes and g not in seen:
                marker_order.append(g)
                marker_ct.append(ct)
                seen.add(g)

    cluster_order = sorted(adata.obs["leiden"].unique(), key=int)
    cluster_labels = [
        f"{c}: {CLUSTER_ANNOTATIONS.get(c, '?')}" for c in cluster_order
    ]

    X = adata[:, marker_order].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    expr_df = pd.DataFrame(X, columns=marker_order, index=adata.obs_names)
    expr_df["cluster"] = adata.obs["leiden"].values

    detected_df = expr_df[marker_order] > 0
    detected_df["cluster"] = expr_df["cluster"]
    pct = detected_df.groupby("cluster", observed=False).mean() * 100
    pct = pct.loc[cluster_order]

    mean_expr = expr_df.groupby("cluster", observed=False)[marker_order].mean()
    mean_expr = mean_expr.loc[cluster_order]
    col_max = mean_expr.max().replace(0, 1)
    mean_scaled = mean_expr / col_max

    n_markers = len(marker_order)
    n_clusters = len(cluster_order)

    fig, ax = plt.subplots(figsize=(max(14, n_clusters * 1.1), max(16, n_markers * 0.45)))

    for i, marker in enumerate(marker_order):
        for j, cluster in enumerate(cluster_order):
            size = float(pct.loc[cluster, marker])
            color_val = float(mean_scaled.loc[cluster, marker])
            ax.scatter(
                j, n_markers - 1 - i,
                s=max(size * 8, 1),
                c=[plt.cm.Reds(color_val)],
                edgecolors="black", linewidth=0.3,
            )

    ax.set_xlim(-0.5, n_clusters - 0.5)
    ax.set_ylim(-0.5, n_markers - 0.5)

    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xticks(range(n_clusters))
    ax.set_xticklabels(cluster_labels, rotation=45, ha="left", fontsize=12)

    ax.set_yticks(range(n_markers))
    ax.set_yticklabels(list(reversed(marker_order)), fontsize=11)

    # Cell type brackets on the left
    ct_groups = []
    prev = None
    start = 0
    for i, ct in enumerate(marker_ct):
        if ct != prev:
            if prev is not None:
                ct_groups.append((prev, start, i - 1))
            start = i
            prev = ct
    ct_groups.append((prev, start, len(marker_ct) - 1))

    for ct_name, idx_start, idx_end in ct_groups:
        y_top = n_markers - 1 - idx_start + 0.35
        y_bot = n_markers - 1 - idx_end - 0.35
        y_mid = (y_top + y_bot) / 2
        bracket_x = -1.8

        ax.plot([bracket_x, bracket_x], [y_bot, y_top],
                color="black", lw=2, clip_on=False,
                transform=ax.transData)
        ax.plot([bracket_x, bracket_x + 0.15], [y_top, y_top],
                color="black", lw=2, clip_on=False)
        ax.plot([bracket_x, bracket_x + 0.15], [y_bot, y_bot],
                color="black", lw=2, clip_on=False)

        ax.text(bracket_x - 0.2, y_mid, ct_name, ha="right", va="center",
                fontsize=12, fontweight="bold", clip_on=False)

    # Size legend
    for pct_val, label in [(20, "20%"), (50, "50%"), (80, "80%"), (100, "100%")]:
        ax.scatter([], [], s=pct_val * 8, c="gray", edgecolors="black",
                   linewidth=0.3, label=label)
    ax.legend(title="% expressing", loc="lower right", fontsize=11,
              title_fontsize=12, framealpha=0.9, labelspacing=1.2)

    ax.set_title("Canonical marker expression by Leiden cluster (10x native)",
                 fontweight="bold", pad=60)
    ax.grid(True, alpha=0.15)

    fig.tight_layout()
    fig.subplots_adjust(left=0.32)
    fig.savefig(FIGURES / "annotation_dotplot_flipped.png", dpi=DPI,
                bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    print("Loading and clustering 10x native...")
    adata = load_and_cluster()

    print("Running DE and building evidence table...")
    evidence = build_evidence_table(adata)
    evidence.to_csv(TABLES / "annotation_evidence.csv", index=False)

    print("\nAnnotation evidence:")
    print("-" * 100)
    for _, row in evidence.iterrows():
        print(
            f"Cluster {row['cluster']:>2s} ({row['n_cells']:5d} cells) "
            f"→ {row['annotation']}"
        )
        print(f"  Top DE: {row['top_de_genes']}")
        print(
            f"  Canonical markers found: {row['n_canonical_in_de']}"
            f"/{row['n_canonical_total']}"
        )
        if row["canonical_markers_in_de"]:
            print(f"  Matched: {row['canonical_markers_in_de']}")
        print()

    print("Generating figures...")
    fig_marker_detection(adata)
    print("  Saved annotation_dotplot.png")

    fig_dotplot_canonical(adata)
    print("  Saved annotation_dotplot_flipped.png")

    print("Done.")


if __name__ == "__main__":
    main()
