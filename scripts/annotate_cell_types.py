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


def fig_evidence_heatmap(adata: ad.AnnData, evidence: pd.DataFrame) -> None:
    """Heatmap of mean expression of top DE genes per cluster."""
    apply_style()

    all_top_genes: list[str] = []
    for _, row in evidence.iterrows():
        genes = [g.split(" (")[0] for g in row["top_de_genes"].split(", ")]
        for g in genes:
            if g not in all_top_genes:
                all_top_genes.append(g)

    X = adata[:, all_top_genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    expr = pd.DataFrame(X, columns=all_top_genes, index=adata.obs_names)
    expr["cluster"] = adata.obs["leiden"].values

    cluster_means = expr.groupby("cluster").mean()
    cluster_means = cluster_means.loc[
        sorted(cluster_means.index, key=int)
    ]

    labels = [
        f"{c}: {evidence.loc[evidence['cluster'] == c, 'annotation'].values[0]}"
        for c in cluster_means.index
    ]

    zscore = (cluster_means - cluster_means.mean()) / cluster_means.std().replace(0, 1)

    fig, ax = plt.subplots(figsize=(max(20, len(all_top_genes) * 0.5), 12))
    sns.heatmap(
        zscore, cmap="RdBu_r", center=0, linewidths=0.3,
        yticklabels=labels, ax=ax,
        cbar_kws={"label": "Z-score (mean expression)"},
    )
    ax.set_title("Top DE genes per cluster (z-scored mean expression)",
                 fontweight="bold")
    ax.set_xlabel("Gene")
    ax.set_ylabel("")
    plt.xticks(rotation=60, ha="right", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "annotation_evidence_heatmap.png", dpi=DPI,
                bbox_inches="tight")
    plt.close(fig)


def fig_dotplot(adata: ad.AnnData) -> None:
    """Dotplot of canonical markers grouped by cell type."""
    apply_style()

    panel_genes = set(adata.var_names)
    marker_dict = {}
    for ct, genes in CANONICAL_MARKERS.items():
        present = [g for g in genes if g in panel_genes]
        if present:
            marker_dict[ct] = present

    fig = sc.pl.dotplot(
        adata, var_names=marker_dict, groupby="leiden",
        standard_scale="var", show=False, return_fig=True,
        figsize=(22, 10),
    )
    fig.savefig(FIGURES / "annotation_dotplot.png", dpi=DPI, bbox_inches="tight")
    plt.close()


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
    fig_evidence_heatmap(adata, evidence)
    print("  Saved annotation_evidence_heatmap.png")

    fig_dotplot(adata)
    print("  Saved annotation_dotplot.png")

    print("Done.")


if __name__ == "__main__":
    main()
