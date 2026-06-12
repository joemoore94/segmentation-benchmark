from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel


def test_load_xenium_calls_spatialdata_io_with_expected_kwargs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench import io

    mock_xenium = MagicMock(return_value="sdata")
    monkeypatch.setattr("spatialdata_io.xenium", mock_xenium)

    result = io.load_xenium(tmp_path)

    mock_xenium.assert_called_once_with(
        tmp_path, morphology_mip=False, aligned_images=False
    )
    assert result == "sdata"


@pytest.fixture
def synthetic_sdata() -> SpatialData:
    image = Image2DModel.parse(
        np.arange(20 * 30, dtype=np.uint8).reshape(1, 20, 30), dims=("c", "y", "x")
    )
    points = PointsModel.parse(
        pd.DataFrame(
            {"x": [1.0, 15.0, 25.0], "y": [1.0, 15.0, 18.0], "gene": ["a", "b", "c"]}
        ),
        coordinates={"x": "x", "y": "y"},
    )
    return SpatialData(images={"morphology": image}, points={"transcripts": points})


def test_crop_roi_subsets_image_and_points(synthetic_sdata: SpatialData) -> None:
    from segbench import io

    cropped = io.crop_roi(synthetic_sdata, xmin=0, ymin=0, xmax=20, ymax=20)

    assert cropped["morphology"].shape == (1, 20, 20)
    transcripts = cropped["transcripts"].compute()
    assert sorted(transcripts["gene"]) == ["a", "b"]


def test_crop_roi_outside_data_raises(synthetic_sdata: SpatialData) -> None:
    from segbench import io

    with pytest.raises(ValueError, match="does not intersect"):
        io.crop_roi(synthetic_sdata, xmin=1000, ymin=1000, xmax=2000, ymax=2000)
