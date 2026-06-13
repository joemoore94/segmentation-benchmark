"""Cross-method comparisons of per-cell AnnData produced by ``quantify``."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree


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
