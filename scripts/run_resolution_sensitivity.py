"""Leiden resolution sensitivity: does the ARI ordering survive resolution changes?

Runs the normalize → PCA → neighbors → Leiden pipeline at five resolutions
(0.5, 0.8, 1.0, 1.5, 2.0) for each method and the 10x-native reference,
computes ARI vs. 10x native, and checks whether the method ordering is stable.

Reads:  data/processed/roi/adata_*.h5ad
Writes: results/tables/resolution_sensitivity.csv
        results/figures/resolution_sensitivity.png

Usage::

    conda run -n segbench python scripts/run_resolution_sensitivity.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from segbench.compare import cell_type_agreement, cluster_cell_types
from segbench.spatial import disagreement_table, morans_i
from segbench.style import apply_style

ROI_DIR   = Path("data/processed/roi")
TABLES    = Path("results/tables")
FIGURES   = Path("results/figures")

RESOLUTIONS = [0.5, 0.8, 1.0, 1.5, 2.0]

METHODS = [
    ("cellpose",       "CellPose"),
    ("stardist",       "StarDist"),
    ("mesmer",         "Mesmer"),
    ("voronoi",        "Voronoi (CP)"),
    ("voronoi_mesmer", "Voronoi (M)"),
    ("baysor",         "Baysor"),
]

METHOD_COLORS = {
    "CellPose":    "#4C72B0",
    "StarDist":    "#8172B2",
    "Mesmer":      "#D62728",
    "Voronoi (CP)":"#17BECF",
    "Voronoi (M)": "#BCBD22",
    "Baysor":      "#DD8452",
}


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 0
    apply_style()

    print("Loading AnnData files...")
    adata_10x = ad.read_h5ad(ROI_DIR / "adata_10x.h5ad")
    adatas = {m: ad.read_h5ad(ROI_DIR / f"adata_{m}.h5ad") for m, _ in METHODS}

    print("Loading match tables...")
    matches = {
        m: pd.read_csv(TABLES / f"disagreement_table_10x_{m}.csv")
        for m, _ in METHODS
    }

    rows = []
    for res in RESOLUTIONS:
        print(f"\nResolution {res}:")
        labels_10x = cluster_cell_types(adata_10x, resolution=res)
        labels_10x.index = labels_10x.index.astype(str)
        n_10x = labels_10x.nunique()
        for method, label in METHODS:
            labels_comp = cluster_cell_types(adatas[method], resolution=res)
            labels_comp.index = labels_comp.index.astype(str)
            n_comp = labels_comp.nunique()
            # Normalise id types to string to avoid int/float/str mismatches
            m = matches[method].copy()
            m["id_a"] = m["id_a"].astype(str)
            m["id_b"] = m["id_b"].astype(str)
            result = cell_type_agreement(labels_10x, labels_comp, m)
            ari = result["ari"]
            dt = disagreement_table(m, labels_10x, labels_comp, adata_10x)
            coords = dt[["centroid_x", "centroid_y"]].to_numpy()
            mi = morans_i(coords, dt["disagree"].to_numpy(dtype=float))
            print(f"  {label:18s}: ARI={ari:.4f}  I={mi:.4f}  clusters: 10x={n_10x}, {label}={n_comp}")
            rows.append({
                "resolution": res,
                "method": label,
                "ari": round(ari, 4),
                "morans_i": round(mi, 4),
                "n_clusters_10x": n_10x,
                "n_clusters_comp": n_comp,
            })

    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "resolution_sensitivity.csv", index=False)
    print("\nSaved resolution_sensitivity.csv")

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(22, 9))

    # Left: ARI vs resolution, one line per method
    ax = axes[0]
    for _, label in METHODS:
        sub = df[df["method"] == label]
        ax.plot(sub["resolution"], sub["ari"], "o-", color=METHOD_COLORS[label],
                label=label, linewidth=2.5, markersize=8)
    ax.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.4, label="default (1.0)")
    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("ARI vs. 10x native")
    ax.set_title("ARI across Leiden resolutions", fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_xticks(RESOLUTIONS)

    # Right: Moran's I of disagreement vs resolution
    ax2 = axes[1]
    for _, label in METHODS:
        sub = df[df["method"] == label]
        ax2.plot(sub["resolution"], sub["morans_i"], "o-", color=METHOD_COLORS[label],
                 label=label, linewidth=2.5, markersize=8)
    ax2.axvline(1.0, color="black", linewidth=1, linestyle="--", alpha=0.4)
    ax2.set_xlabel("Leiden resolution")
    ax2.set_ylabel("Global Moran's I of disagreement")
    ax2.set_title("Spatial structure of disagreement across resolutions", fontweight="bold")
    ax2.set_xticks(RESOLUTIONS)
    ax2.legend(fontsize=10)

    fig.suptitle(
        "Leiden resolution sensitivity: ARI and spatial autocorrelation vs. 10x native",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "resolution_sensitivity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved resolution_sensitivity.png")


if __name__ == "__main__":
    main()
