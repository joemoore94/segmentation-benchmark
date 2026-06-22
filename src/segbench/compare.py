"""Cross-method comparisons of per-cell AnnData produced by ``quantify``."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from numpy.typing import NDArray
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score


def subset_to_region(
    adata: AnnData, x_range: tuple[float, float], y_range: tuple[float, float]
) -> AnnData:
    """Subset ``adata`` to cells whose centroid falls within a rectangular region.

    ``x_range``/``y_range`` are ``(min, max)`` bounds in the same physical
    units as ``obs["centroid_x"]``/``obs["centroid_y"]`` (microns, ROI-local).
    Used to compare methods on the exact same physical area when their
    full segmentation runs covered different-sized ROIs.
    """
    x = adata.obs["centroid_x"]
    y = adata.obs["centroid_y"]
    mask = x.between(*x_range) & y.between(*y_range)
    return adata[mask].copy()


def cell_count_summary(adatas: dict[str, AnnData]) -> pd.DataFrame:
    """Per-method cell counts and transcripts-per-cell summary statistics.

    ``X.sum(axis=1)`` (total assigned transcripts per cell) is comparable
    across methods regardless of how each defines cell "size" (pixel area
    for CellPose/Mesmer, ``n_transcripts`` for Baysor).
    """
    rows = []
    for method, adata in adatas.items():
        counts = np.asarray(adata.X.sum(axis=1)).ravel()
        rows.append(
            {
                "method": method,
                "n_cells": adata.n_obs,
                "total_transcripts": counts.sum(),
                "median_transcripts_per_cell": np.median(counts),
                "mean_transcripts_per_cell": counts.mean(),
            }
        )
    return pd.DataFrame(rows).set_index("method")


def match_cells_by_centroid(
    adata_a: AnnData, adata_b: AnnData, max_dist: float
) -> pd.DataFrame:
    """Pair up cells from two methods by nearest-neighbor centroid distance.

    ``adata_a.obs``/``adata_b.obs`` must have ``centroid_x``/``centroid_y`` in
    the same physical units (e.g. microns from a common ROI origin). Matching
    is mutual nearest-neighbor: cell ``a`` is paired with cell ``b`` only if
    each is the other's closest centroid and the distance is ``<= max_dist``.

    Returns
    -------
    DataFrame with columns ``id_a``, ``id_b``, ``distance``, one row per
    matched pair. Cells with no match within ``max_dist`` are omitted.
    """
    xy_a = adata_a.obs[["centroid_x", "centroid_y"]].to_numpy()
    xy_b = adata_b.obs[["centroid_x", "centroid_y"]].to_numpy()

    tree_a = cKDTree(xy_a)
    tree_b = cKDTree(xy_b)

    dist_ab, idx_b = tree_b.query(xy_a)
    dist_ba, idx_a = tree_a.query(xy_b)

    pairs = []
    for i_a, (i_b, d) in enumerate(zip(idx_b, dist_ab)):
        if d <= max_dist and idx_a[i_b] == i_a:
            pairs.append((adata_a.obs_names[i_a], adata_b.obs_names[i_b], d))

    return pd.DataFrame(pairs, columns=["id_a", "id_b", "distance"])


def size_summary(adatas: dict[str, AnnData]) -> pd.DataFrame:
    """Per-method cell "size" distribution.

    Size is pixel area (``obs["area"]``, for mask-based methods like CellPose
    and Mesmer) or transcript count (``obs["n_transcripts"]``, for Baysor),
    whichever is present. These are not directly comparable across methods,
    but the percentile spread within each method shows how broad/narrow each
    method's size distribution is.
    """
    rows = []
    for method, adata in adatas.items():
        size_col = "area" if "area" in adata.obs.columns else "n_transcripts"
        sizes = adata.obs[size_col].to_numpy(dtype=float)
        rows.append(
            {
                "method": method,
                "size_col": size_col,
                "p10": np.percentile(sizes, 10),
                "p25": np.percentile(sizes, 25),
                "median": np.median(sizes),
                "p75": np.percentile(sizes, 75),
                "p90": np.percentile(sizes, 90),
                "mean": sizes.mean(),
                "std": sizes.std(),
            }
        )
    return pd.DataFrame(rows).set_index("method")


def _dense(adata: AnnData) -> NDArray[np.floating]:
    x = adata.X
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float64)


def expression_correlation(
    adata_a: AnnData, adata_b: AnnData, matches: pd.DataFrame
) -> pd.DataFrame:
    """Per-cell-pair Pearson correlation of expression profiles, over shared genes.

    For each row in ``matches`` (as returned by :func:`match_cells_by_centroid`),
    correlates ``adata_a``'s and ``adata_b``'s count vectors for that cell pair
    across the genes present in both ``adata_a.var_names`` and
    ``adata_b.var_names``. Pairs where either cell has zero variance (e.g. all
    zero counts) get ``NaN`` correlation.
    """
    shared_genes = adata_a.var_names.intersection(adata_b.var_names)
    xa = _dense(adata_a[:, shared_genes])
    xb = _dense(adata_b[:, shared_genes])

    idx_a = {name: i for i, name in enumerate(adata_a.obs_names)}
    idx_b = {name: i for i, name in enumerate(adata_b.obs_names)}

    correlations = []
    for id_a, id_b in zip(matches["id_a"], matches["id_b"]):
        va = xa[idx_a[id_a]]
        vb = xb[idx_b[id_b]]
        if va.std() == 0 or vb.std() == 0:
            correlations.append(np.nan)
        else:
            correlations.append(float(np.corrcoef(va, vb)[0, 1]))

    out = matches.copy()
    out["correlation"] = correlations
    return out


def _cluster_pipeline(
    adata: AnnData,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
    seed: int,
    compute_umap: bool,
) -> AnnData:
    a = adata.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    n_comps = min(n_pcs, a.n_vars - 1, a.n_obs - 1)
    sc.pp.pca(a, n_comps=n_comps, random_state=seed)
    sc.pp.neighbors(a, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(a, resolution=resolution, random_state=seed, flavor="igraph")
    if compute_umap:
        sc.tl.umap(a, random_state=seed)
    return a


def cluster_cell_types(
    adata: AnnData,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    seed: int = 0,
) -> pd.Series:
    """Standard normalize -> PCA -> neighbors -> Leiden pipeline.

    Returns a ``pd.Series`` of Leiden cluster labels indexed by
    ``adata.obs_names``. Operates on a copy; ``adata`` itself is unmodified.
    """
    a = _cluster_pipeline(adata, n_pcs, n_neighbors, resolution, seed, compute_umap=False)
    return a.obs["leiden"].copy()


def cluster_embedding(
    adata: AnnData,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    seed: int = 0,
) -> pd.DataFrame:
    """Same pipeline as :func:`cluster_cell_types`, plus a 2D UMAP embedding.

    Returns a ``pd.DataFrame`` indexed by ``adata.obs_names`` with columns
    ``PC1``/``PC2`` (first two PCA components), ``UMAP1``/``UMAP2``, and
    ``leiden`` (cluster label) -- for visualizing how distinct each method's
    Leiden clusters are in PCA/UMAP space. Operates on a copy; ``adata``
    itself is unmodified.
    """
    a = _cluster_pipeline(adata, n_pcs, n_neighbors, resolution, seed, compute_umap=True)
    return pd.DataFrame(
        {
            "PC1": a.obsm["X_pca"][:, 0],
            "PC2": a.obsm["X_pca"][:, 1],
            "UMAP1": a.obsm["X_umap"][:, 0],
            "UMAP2": a.obsm["X_umap"][:, 1],
            "leiden": a.obs["leiden"].to_numpy(),
        },
        index=a.obs_names,
    )


def cell_type_agreement(
    labels_a: pd.Series, labels_b: pd.Series, matches: pd.DataFrame
) -> dict[str, object]:
    """Agreement between two methods' cell-type labels, for matched cell pairs.

    Returns a dict with ``"ari"`` (adjusted Rand index over matched pairs
    where both cells have a label), ``"confusion"`` (a ``label_a`` x
    ``label_b`` contingency table), and ``"n_matched"`` (number of pairs used).
    """
    label_a = matches["id_a"].map(labels_a)
    label_b = matches["id_b"].map(labels_b)
    valid = label_a.notna() & label_b.notna()
    label_a = label_a[valid]
    label_b = label_b[valid]

    ari = adjusted_rand_score(label_a, label_b) if len(label_a) > 0 else float("nan")
    confusion = pd.crosstab(
        label_a.rename("label_a"), label_b.rename("label_b")
    )
    return {"ari": ari, "confusion": confusion, "n_matched": int(valid.sum())}


def match_cluster_labels(
    labels_a: pd.Series, labels_b: pd.Series, matches: pd.DataFrame
) -> pd.Series:
    """Relabel ``labels_b``'s clusters onto ``labels_a``'s cluster vocabulary.

    Independent per-method clusterings (e.g. two separate
    :func:`cluster_cell_types` runs) assign arbitrary cluster ids, so
    "cluster 3" in method A has no relationship to "cluster 3" in method B.
    This finds a one-to-one assignment between A's and B's clusters that
    maximizes overlap (the Hungarian algorithm on the ``label_a`` x
    ``label_b`` contingency table over ``matches``), then renames B's
    clusters to their matched A cluster name. B clusters with no counterpart
    (e.g. if B has more clusters than A) are renamed ``"b_only_<original>"``.

    Returns ``labels_b`` with values renamed accordingly (same index).
    """
    label_a = matches["id_a"].map(labels_a)
    label_b = matches["id_b"].map(labels_b)
    valid = label_a.notna() & label_b.notna()
    confusion = pd.crosstab(label_a[valid], label_b[valid])

    row_ind, col_ind = linear_sum_assignment(-confusion.to_numpy())
    rename = {confusion.columns[c]: confusion.index[r] for r, c in zip(row_ind, col_ind)}

    for name in labels_b.unique():
        rename.setdefault(name, f"b_only_{name}")

    return labels_b.map(rename)
