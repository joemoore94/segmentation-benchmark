from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from segbench.compare import (
    cell_count_summary,
    cell_type_agreement,
    cluster_cell_types,
    cluster_embedding,
    expression_correlation,
    match_cells_by_centroid,
    match_cluster_labels,
    size_summary,
)


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


def test_size_summary_uses_area_or_n_transcripts() -> None:
    var = pd.DataFrame(index=["g1"])

    a = AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"area": [10.0, 20.0, 30.0]}, index=["1", "2", "3"]),
        var=var,
    )
    b = AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame({"n_transcripts": [5.0, 15.0]}, index=["c1", "c2"]),
        var=var,
    )

    summary = size_summary({"a": a, "b": b})

    assert summary.loc["a", "size_col"] == "area"
    assert summary.loc["a", "median"] == 20.0
    assert summary.loc["b", "size_col"] == "n_transcripts"
    assert summary.loc["b", "median"] == 10.0


def test_expression_correlation_perfect_and_zero_variance() -> None:
    var = pd.DataFrame(index=["g1", "g2", "g3"])
    a = AnnData(
        X=np.array([[1.0, 2.0, 3.0]], dtype=np.float32),
        obs=pd.DataFrame({"centroid_x": [0.0], "centroid_y": [0.0]}, index=["a1"]),
        var=var,
    )
    b = AnnData(
        X=np.array([[2.0, 4.0, 6.0], [5.0, 5.0, 5.0]], dtype=np.float32),
        obs=pd.DataFrame(
            {"centroid_x": [0.0, 0.0], "centroid_y": [0.0, 0.0]}, index=["b1", "b2"]
        ),
        var=var,
    )
    matches = pd.DataFrame(
        {"id_a": ["a1", "a1"], "id_b": ["b1", "b2"], "distance": [0.0, 0.0]}
    )

    result = expression_correlation(a, b, matches)

    perfect = result.loc[result["id_b"] == "b1", "correlation"].item()
    zero_var = result.loc[result["id_b"] == "b2", "correlation"].item()
    assert perfect == pytest.approx(1.0)
    assert np.isnan(zero_var)


def test_cluster_cell_types_separates_distinct_groups() -> None:
    rng = np.random.default_rng(0)
    n_per_group = 15
    group_a = rng.poisson(lam=[20, 20, 1, 1], size=(n_per_group, 4)).astype(np.float32)
    group_b = rng.poisson(lam=[1, 1, 20, 20], size=(n_per_group, 4)).astype(np.float32)
    x = np.vstack([group_a, group_b])

    adata = AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(2 * n_per_group)]),
        var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
    )

    labels = cluster_cell_types(adata, n_neighbors=10, seed=0)

    assert list(labels.index) == list(adata.obs_names)
    assert labels.nunique() >= 2


def test_cluster_embedding_has_pca_and_umap_columns() -> None:
    rng = np.random.default_rng(0)
    n_per_group = 15
    group_a = rng.poisson(lam=[20, 20, 1, 1], size=(n_per_group, 4)).astype(np.float32)
    group_b = rng.poisson(lam=[1, 1, 20, 20], size=(n_per_group, 4)).astype(np.float32)
    x = np.vstack([group_a, group_b])

    adata = AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(2 * n_per_group)]),
        var=pd.DataFrame(index=["g1", "g2", "g3", "g4"]),
    )

    embedding = cluster_embedding(adata, n_neighbors=10, seed=0)

    assert list(embedding.index) == list(adata.obs_names)
    assert list(embedding.columns) == ["PC1", "PC2", "UMAP1", "UMAP2", "leiden"]
    assert embedding["leiden"].nunique() >= 2


def test_cell_type_agreement_ari_and_confusion() -> None:
    labels_a = pd.Series({"0": "A", "1": "A", "2": "B"})
    labels_b = pd.Series({"x": "A", "y": "A", "z": "B"})
    matches = pd.DataFrame(
        {"id_a": ["0", "1", "2"], "id_b": ["x", "y", "z"], "distance": [0.0, 0.0, 0.0]}
    )

    result = cell_type_agreement(labels_a, labels_b, matches)

    assert result["ari"] == pytest.approx(1.0)
    assert result["n_matched"] == 3
    assert result["confusion"].loc["A", "A"] == 2
    assert result["confusion"].loc["B", "B"] == 1


def test_match_cluster_labels_aligns_independent_cluster_ids() -> None:
    # method A calls these clusters "0"/"1"; method B calls the *same*
    # underlying populations "5"/"9" (arbitrary, unrelated ids).
    labels_a = pd.Series({"0": "0", "1": "0", "2": "1"})
    labels_b = pd.Series({"x": "9", "y": "9", "z": "5"})
    matches = pd.DataFrame(
        {"id_a": ["0", "1", "2"], "id_b": ["x", "y", "z"], "distance": [0.0, 0.0, 0.0]}
    )

    remapped = match_cluster_labels(labels_a, labels_b, matches)

    assert remapped["x"] == "0"
    assert remapped["y"] == "0"
    assert remapped["z"] == "1"


def test_match_cluster_labels_extra_b_clusters_get_b_only_prefix() -> None:
    labels_a = pd.Series({"0": "0", "1": "0"})
    labels_b = pd.Series({"x": "5", "y": "5", "z": "7"})
    matches = pd.DataFrame(
        {"id_a": ["0", "1"], "id_b": ["x", "y"], "distance": [0.0, 0.0]}
    )

    remapped = match_cluster_labels(labels_a, labels_b, matches)

    assert remapped["x"] == "0"
    assert remapped["y"] == "0"
    assert remapped["z"] == "b_only_7"
