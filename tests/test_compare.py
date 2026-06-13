from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from segbench.compare import cell_count_summary, match_cells_by_centroid


def _make_adata(counts: list[list[float]], centroids: list[tuple[float, float]]) -> AnnData:
    obs = pd.DataFrame(
        {
            "centroid_x": [c[0] for c in centroids],
            "centroid_y": [c[1] for c in centroids],
        },
        index=[str(i) for i in range(len(centroids))],
    )
    var = pd.DataFrame(index=["geneA", "geneB"])
    return AnnData(X=np.array(counts, dtype=np.float32), obs=obs, var=var)


def test_cell_count_summary() -> None:
    a = _make_adata([[1, 1], [2, 0]], [(0, 0), (10, 10)])
    b = _make_adata([[3, 1], [1, 1], [0, 2]], [(0, 0), (10, 10), (20, 20)])

    summary = cell_count_summary({"a": a, "b": b})

    assert summary.loc["a", "n_cells"] == 2
    assert summary.loc["a", "total_transcripts"] == 4
    assert summary.loc["a", "median_transcripts_per_cell"] == 2
    assert summary.loc["b", "n_cells"] == 3
    assert summary.loc["b", "total_transcripts"] == 8


def test_match_cells_by_centroid_mutual_nearest_neighbor() -> None:
    # a has 3 cells; b has 2 cells close to a's cells 0 and 1, and one far away.
    a = _make_adata(
        [[1, 0], [1, 0], [1, 0]], [(0.0, 0.0), (10.0, 0.0), (100.0, 100.0)]
    )
    b = _make_adata([[1, 0], [1, 0]], [(0.5, 0.0), (10.5, 0.0)])

    matches = match_cells_by_centroid(a, b, max_dist=2.0)

    assert len(matches) == 2
    assert set(matches["id_a"]) == {"0", "1"}
    assert set(matches["id_b"]) == {"0", "1"}
    assert (matches["distance"] <= 2.0).all()


def test_match_cells_by_centroid_respects_max_dist() -> None:
    a = _make_adata([[1, 0]], [(0.0, 0.0)])
    b = _make_adata([[1, 0]], [(5.0, 0.0)])

    matches = match_cells_by_centroid(a, b, max_dist=2.0)

    assert len(matches) == 0
