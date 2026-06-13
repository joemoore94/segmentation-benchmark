from __future__ import annotations

import numpy as np
import pandas as pd

from segbench.quantify import quantify_baysor, quantify_cells


def test_quantify_cells_counts_and_obs() -> None:
    masks = np.zeros((10, 10), dtype=np.int32)
    masks[1:4, 1:4] = 1  # cell 1: rows/cols 1-3
    masks[6:9, 6:9] = 2  # cell 2: rows/cols 6-8

    transcripts = pd.DataFrame(
        {
            "x_location": [1.0, 2.0, 7.0, 7.0, 8.0, 50.0],
            "y_location": [1.0, 2.0, 7.0, 7.0, 7.0, 50.0],
            "feature_name": ["geneA", "geneB", "geneA", "geneA", "geneB", "geneA"],
        }
    )

    adata = quantify_cells(masks, transcripts)

    assert set(adata.obs_names) == {"1", "2"}
    assert set(adata.var_names) == {"geneA", "geneB"}

    cell1 = adata["1"]
    cell2 = adata["2"]
    assert cell1[:, "geneA"].X.item() == 1
    assert cell1[:, "geneB"].X.item() == 1
    assert cell2[:, "geneA"].X.item() == 2
    assert cell2[:, "geneB"].X.item() == 1

    assert adata.obs.loc["1", "area"] == 9
    assert adata.obs.loc["2", "area"] == 9


def test_quantify_cells_centroid_in_physical_units() -> None:
    masks = np.zeros((10, 10), dtype=np.int32)
    masks[0:2, 0:2] = 1  # cell 1: rows/cols 0-1, pixel centroid (0.5, 0.5)

    transcripts = pd.DataFrame(
        {
            "x_location": [100.0],
            "y_location": [200.0],
            "feature_name": ["geneA"],
        }
    )

    adata = quantify_cells(masks, transcripts, pixel_size=2.0, origin=(100.0, 200.0))

    assert adata.obs.loc["1", "centroid_x"] == 0.5 * 2.0 + 100.0
    assert adata.obs.loc["1", "centroid_y"] == 0.5 * 2.0 + 200.0


def test_quantify_cells_zero_transcript_cell_is_zero_filled() -> None:
    masks = np.zeros((10, 10), dtype=np.int32)
    masks[1:4, 1:4] = 1
    masks[6:9, 6:9] = 2

    transcripts = pd.DataFrame(
        {
            "x_location": [2.0],
            "y_location": [2.0],
            "feature_name": ["geneA"],
        }
    )

    adata = quantify_cells(masks, transcripts)

    assert set(adata.obs_names) == {"1", "2"}
    assert adata["2"].X.sum() == 0


def test_quantify_baysor_counts_and_centroids() -> None:
    segmentation = pd.DataFrame(
        {
            "x": [1.0, 3.0, 10.0, 12.0, 5.0],
            "y": [1.0, 3.0, 10.0, 12.0, 50.0],
            "gene": ["geneA", "geneB", "geneA", "geneA", "geneB"],
            "cell": ["c1", "c1", "c2", "c2", np.nan],
        }
    )

    adata = quantify_baysor(segmentation)

    assert set(adata.obs_names) == {"c1", "c2"}
    assert set(adata.var_names) == {"geneA", "geneB"}

    c1 = adata["c1"]
    c2 = adata["c2"]
    assert c1[:, "geneA"].X.item() == 1
    assert c1[:, "geneB"].X.item() == 1
    assert c2[:, "geneA"].X.item() == 2
    assert c2[:, "geneB"].X.item() == 0

    assert adata.obs.loc["c1", "centroid_x"] == 2.0
    assert adata.obs.loc["c1", "centroid_y"] == 2.0
    assert adata.obs.loc["c1", "n_transcripts"] == 2
    assert adata.obs.loc["c2", "centroid_x"] == 11.0
    assert adata.obs.loc["c2", "n_transcripts"] == 2
