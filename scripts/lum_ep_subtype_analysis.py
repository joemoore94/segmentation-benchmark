"""Luminal epithelial subtype analysis on 10x native segmentation.

Two questions:
1. What markers differentiate the lum ep subclusters (0, 1, 3, 8 at res=1.0)?
2. At what Leiden resolution do these subclusters first split?
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from segbench.constants import CLUSTER_ANNOTATIONS

DATA = Path(__file__).resolve().parents[1] / "data" / "processed" / "roi"
RESULTS = Path(__file__).resolve().parents[1] / "results"
FIG_DIR = RESULTS / "figures" / "lum_ep_subtypes"
TABLE_DIR = RESULTS / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

LUM_EP_CLUSTERS = {"0", "1", "3", "8"}


def preprocess(adata: ad.AnnData, resolution: float = 1.0) -> ad.AnnData:
    a = adata.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=30, random_state=0)
    sc.pp.neighbors(a, n_neighbors=15, random_state=0)
    sc.tl.leiden(a, resolution=resolution, random_state=0, flavor="igraph")
    return a


# ── Load and cluster at reference resolution ─────────────────────────────
raw = ad.read_h5ad(DATA / "adata_10x.h5ad")
adata = preprocess(raw, resolution=1.0)
sc.tl.umap(adata, random_state=0)

adata.obs["cell_type"] = adata.obs["leiden"].map(CLUSTER_ANNOTATIONS).astype("category")
print(f"Total cells: {adata.n_obs}")
print(f"Leiden clusters at res=1.0: {adata.obs['leiden'].nunique()}")

lum_mask = adata.obs["leiden"].isin(LUM_EP_CLUSTERS)
n_lum = lum_mask.sum()
print(f"\nLuminal epithelial cells (clusters {sorted(LUM_EP_CLUSTERS)}): {n_lum}")
print(adata.obs.loc[lum_mask, "leiden"].value_counts().sort_index())

# ── 1. DE between lum ep subclusters ─────────────────────────────────────
lum = adata[lum_mask].copy()
lum.obs["subcluster"] = lum.obs["leiden"].copy()

sc.tl.rank_genes_groups(lum, groupby="subcluster", method="wilcoxon",
                        use_raw=False)

print("\n=== Top 10 markers per luminal epithelial subcluster ===")
top_n = 10
result_df = sc.get.rank_genes_groups_df(lum, group=None)
for cluster_id in sorted(LUM_EP_CLUSTERS):
    sub = result_df[result_df["group"] == cluster_id].head(top_n)
    print(f"\nCluster {cluster_id} ({lum.obs['subcluster'].value_counts()[cluster_id]} cells):")
    for _, row in sub.iterrows():
        print(f"  {row['names']:12s}  logFC={row['logfoldchanges']:+.2f}  "
              f"pval_adj={row['pvals_adj']:.2e}  score={row['scores']:.1f}")

result_df.to_csv(TABLE_DIR / "lum_ep_subcluster_markers_10x.csv", index=False)

# Dotplot of top markers
top_genes = {}
for c in sorted(LUM_EP_CLUSTERS):
    genes = result_df[result_df["group"] == c].head(5)["names"].tolist()
    top_genes[c] = genes

all_top_genes = []
for c in sorted(LUM_EP_CLUSTERS):
    for g in top_genes[c]:
        if g not in all_top_genes:
            all_top_genes.append(g)

fig, ax = plt.subplots(figsize=(12, 5))
sc.pl.dotplot(lum, var_names=all_top_genes, groupby="subcluster",
              ax=ax, show=False)
plt.suptitle("Top markers per luminal epithelial subcluster (10x native, res=1.0)",
             y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "lum_ep_subcluster_dotplot.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"\nSaved dotplot → {FIG_DIR / 'lum_ep_subcluster_dotplot.png'}")

# Heatmap of top markers
sc.pl.rank_genes_groups_heatmap(lum, n_genes=8, groupby="subcluster",
                                show=False, use_raw=False)
plt.savefig(FIG_DIR / "lum_ep_subcluster_heatmap.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved heatmap → {FIG_DIR / 'lum_ep_subcluster_heatmap.png'}")

# UMAP colored by subcluster (lum ep cells only)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.umap(adata, color="cell_type", ax=axes[0], show=False, title="All cells by type")
sc.pl.umap(adata, color="leiden", ax=axes[1], show=False,
           title="All cells by Leiden cluster (res=1.0)",
           groups=list(LUM_EP_CLUSTERS),
           na_color="lightgrey", na_in_legend=False)
plt.tight_layout()
plt.savefig(FIG_DIR / "lum_ep_umap_context.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved UMAP context → {FIG_DIR / 'lum_ep_umap_context.png'}")

# ── 2. Resolution sweep: when do lum ep clusters split? ──────────────────
resolutions = np.arange(0.1, 2.05, 0.1)
sweep_results = []

for res in resolutions:
    a = preprocess(raw, resolution=res)
    a.obs["cell_type"] = a.obs["leiden"].map(
        lambda x: CLUSTER_ANNOTATIONS.get(x, "Unknown")
    )
    n_clusters = a.obs["leiden"].nunique()

    lum_labels = a.obs.loc[lum_mask, "leiden"]
    n_lum_clusters = lum_labels.nunique()
    lum_cluster_sizes = lum_labels.value_counts().to_dict()

    sweep_results.append({
        "resolution": round(res, 1),
        "total_clusters": n_clusters,
        "lum_ep_clusters": n_lum_clusters,
        "lum_ep_cluster_sizes": lum_cluster_sizes,
    })
    print(f"res={res:.1f}: {n_clusters} total clusters, "
          f"{n_lum_clusters} contain lum ep cells "
          f"(sizes: {dict(sorted(lum_cluster_sizes.items(), key=lambda x: -x[1]))})")

sweep_df = pd.DataFrame([
    {"resolution": r["resolution"],
     "total_clusters": r["total_clusters"],
     "lum_ep_clusters": r["lum_ep_clusters"]}
    for r in sweep_results
])
sweep_df.to_csv(TABLE_DIR / "lum_ep_resolution_sweep_10x.csv", index=False)

# Plot resolution sweep
fig, ax1 = plt.subplots(figsize=(8, 4))
ax1.plot(sweep_df["resolution"], sweep_df["total_clusters"],
         "o-", color="steelblue", label="Total clusters")
ax1.set_xlabel("Leiden resolution")
ax1.set_ylabel("Total clusters", color="steelblue")
ax1.tick_params(axis="y", labelcolor="steelblue")

ax2 = ax1.twinx()
ax2.plot(sweep_df["resolution"], sweep_df["lum_ep_clusters"],
         "s-", color="orchid", label="Lum ep subclusters")
ax2.set_ylabel("Luminal epithelial subclusters", color="orchid")
ax2.tick_params(axis="y", labelcolor="orchid")

ax1.axvline(1.0, color="grey", ls="--", alpha=0.5, label="res=1.0 (default)")
ax1.set_title("Leiden resolution sweep — luminal epithelial splitting (10x native)")
fig.legend(loc="upper left", bbox_to_anchor=(0.12, 0.88))
plt.tight_layout()
plt.savefig(FIG_DIR / "lum_ep_resolution_sweep.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"\nSaved resolution sweep → {FIG_DIR / 'lum_ep_resolution_sweep.png'}")

print("\nDone.")
