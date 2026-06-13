from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from segbench.spatial import (
    disagreement_spatial_structure,
    disagreement_table,
    morans_i,
    morans_i_permutation_test,
)


def _grid_coords(n_side: int) -> np.ndarray:
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(float)


def test_morans_i_high_for_spatially_clustered_values() -> None:
    coords = _grid_coords(10)
    # left half = 1, right half = 0: strong spatial clustering.
    values = (coords[:, 0] < 5).astype(float)

    i = morans_i(coords, values, k=4)

    assert i > 0.5


def test_morans_i_near_zero_for_checkerboard() -> None:
    coords = _grid_coords(10)
    # checkerboard pattern: neighbors always differ.
    values = ((coords[:, 0].astype(int) + coords[:, 1].astype(int)) % 2).astype(float)

    i = morans_i(coords, values, k=4)

    assert i < 0


def test_morans_i_permutation_test_significant_for_clustered_pattern() -> None:
    coords = _grid_coords(10)
    values = (coords[:, 0] < 5).astype(float)

    observed, p_value = morans_i_permutation_test(coords, values, k=4, n_perm=199, seed=0)

    assert observed > 0.5
    assert p_value < 0.05


def test_disagreement_table_merges_labels_and_coords() -> None:
    matches = pd.DataFrame(
        {"id_a": ["a1", "a2", "a3"], "id_b": ["b1", "b2", "b3"], "distance": [0.0, 0.0, 0.0]}
    )
    labels_a = pd.Series({"a1": "X", "a2": "X", "a3": "Y"})
    labels_b = pd.Series({"b1": "X", "b2": "Y", "b3": "Y"})

    adata_a = AnnData(
        X=np.zeros((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"centroid_x": [0.0, 1.0, 2.0], "centroid_y": [0.0, 1.0, 2.0]},
            index=["a1", "a2", "a3"],
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    table = disagreement_table(matches, labels_a, labels_b, adata_a)

    assert len(table) == 3
    assert table.loc[table["id_a"] == "a1", "disagree"].item() == 0.0
    assert table.loc[table["id_a"] == "a2", "disagree"].item() == 1.0
    assert table.loc[table["id_a"] == "a3", "disagree"].item() == 0.0
    assert table.loc[table["id_a"] == "a2", "centroid_x"].item() == 1.0


def test_disagreement_table_handles_mismatched_categoricals() -> None:
    """labels_a/labels_b may be independent Leiden clusterings (Categorical
    dtype with different category sets), which can't be compared with !=
    directly."""
    matches = pd.DataFrame(
        {"id_a": ["a1", "a2"], "id_b": ["b1", "b2"], "distance": [0.0, 0.0]}
    )
    labels_a = pd.Series({"a1": "0", "a2": "1"}, dtype="category")
    labels_b = pd.Series({"b1": "0", "b2": "0"}, dtype="category")

    adata_a = AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"centroid_x": [0.0, 1.0], "centroid_y": [0.0, 1.0]}, index=["a1", "a2"]
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    table = disagreement_table(matches, labels_a, labels_b, adata_a)

    assert table.loc[table["id_a"] == "a1", "disagree"].item() == 0.0
    assert table.loc[table["id_a"] == "a2", "disagree"].item() == 1.0


def test_disagreement_table_drops_unlabeled_pairs() -> None:
    matches = pd.DataFrame({"id_a": ["a1", "a2"], "id_b": ["b1", "b2"], "distance": [0.0, 0.0]})
    labels_a = pd.Series({"a1": "X"})  # a2 has no label
    labels_b = pd.Series({"b1": "X", "b2": "Y"})

    adata_a = AnnData(
        X=np.zeros((2, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"centroid_x": [0.0, 1.0], "centroid_y": [0.0, 1.0]}, index=["a1", "a2"]
        ),
        var=pd.DataFrame(index=["g1"]),
    )

    table = disagreement_table(matches, labels_a, labels_b, adata_a)

    assert len(table) == 1
    assert table.loc[0, "id_a"] == "a1"


def test_disagreement_spatial_structure_returns_expected_keys() -> None:
    n_side = 10
    coords = _grid_coords(n_side)
    disagree = (coords[:, 0] < 5).astype(float)

    table = pd.DataFrame(
        {
            "centroid_x": coords[:, 0],
            "centroid_y": coords[:, 1],
            "disagree": disagree,
        }
    )

    result = disagreement_spatial_structure(table, k=4, n_perm=199, seed=0)

    assert result["n_cells"] == n_side * n_side
    assert result["disagreement_rate"] == 0.5
    assert result["morans_i"] > 0.5
    assert result["p_value"] < 0.05
