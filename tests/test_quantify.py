from __future__ import annotations

import numpy as np
import pandas as pd

from segbench.quantify import quantify_cells


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
