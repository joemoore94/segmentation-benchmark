from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    """A miniature Xenium bundle, in the same layout as ``data/raw``."""
    # DAPI: 20x30 px, written uncompressed (compression isn't what's under test).
    dapi = np.arange(20 * 30, dtype=np.uint16).reshape(20, 30)
    tifffile.imwrite(tmp_path / "morphology_focus.ome.tif", dapi)

    # H&E: 3 separate single-band pages (Hematoxylin/Eosin/Residual).
    he_pages = np.stack(
        [np.full((20, 30), v, dtype=np.uint8) for v in (10, 20, 30)]
    )
    tifffile.imwrite(tmp_path / "he_image.ome.tif", he_pages)

    transcripts = pd.DataFrame(
        {
            "x_location": [1.0, 5.0, 25.0],
            "y_location": [1.0, 5.0, 18.0],
            "feature_name": [b"geneA", b"geneB", b"geneC"],
        }
    )
    transcripts.to_parquet(tmp_path / "transcripts.parquet")

    cells = pd.DataFrame(
        {
            "cell_id": [b"cell-1", b"cell-2", b"cell-3"],
            "x_centroid": [2.0, 6.0, 25.0],
            "y_centroid": [2.0, 6.0, 18.0],
        }
    )
    cells.to_parquet(tmp_path / "cells.parquet")

    for fname in ("cell_boundaries.parquet", "nucleus_boundaries.parquet"):
        boundaries = pd.DataFrame(
            {
                "cell_id": [b"cell-1", b"cell-1", b"cell-3", b"cell-3"],
                "vertex_x": [1.5, 2.5, 24.5, 25.5],
                "vertex_y": [1.5, 2.5, 17.5, 18.5],
            }
        )
        boundaries.to_parquet(tmp_path / fname)

    return tmp_path


def test_load_roi_dapi_crops_in_pixel_space(data_dir: Path) -> None:
    from segbench import io

    crop = io.load_roi_dapi(data_dir, xmin=0, ymin=0, xmax=10, ymax=10, pixel_size=1.0)

    assert crop.shape == (10, 10)
    full = tifffile.imread(data_dir / "morphology_focus.ome.tif")
    np.testing.assert_array_equal(crop, full[0:10, 0:10])


def test_load_roi_he_bandjoins_pages(data_dir: Path) -> None:
    from segbench import io

    crop = io.load_roi_he(data_dir, xmin=0, ymin=0, xmax=10, ymax=10, pixel_size=1.0)

    assert crop.shape == (10, 10, 3)
    assert np.all(crop[:, :, 0] == 10)
    assert np.all(crop[:, :, 1] == 20)
    assert np.all(crop[:, :, 2] == 30)


def test_load_roi_transcripts_crops_and_shifts_origin(data_dir: Path) -> None:
    from segbench import io

    cropped = io.load_roi_transcripts(data_dir, xmin=0, ymin=0, xmax=10, ymax=10)

    assert sorted(cropped["feature_name"]) == [b"geneA", b"geneB"]
    assert cropped["x_location"].tolist() == [1.0, 5.0]
    assert cropped["y_location"].tolist() == [1.0, 5.0]


def test_load_roi_cells_crops_and_shifts_origin(data_dir: Path) -> None:
    from segbench import io

    cropped = io.load_roi_cells(data_dir, xmin=0, ymin=0, xmax=10, ymax=10)

    assert sorted(cropped["cell_id"]) == [b"cell-1", b"cell-2"]
    assert cropped["x_centroid"].tolist() == [2.0, 6.0]
    assert cropped["y_centroid"].tolist() == [2.0, 6.0]


def test_load_roi_boundaries_filters_by_cell_id_and_shifts_origin(
    data_dir: Path,
) -> None:
    from segbench import io

    cell_ids = pd.Series([b"cell-1"])
    cropped = io.load_roi_boundaries(data_dir, cell_ids, xmin=1.0, ymin=1.0, kind="cell")

    assert set(cropped["cell_id"]) == {b"cell-1"}
    assert cropped["vertex_x"].tolist() == [0.5, 1.5]
    assert cropped["vertex_y"].tolist() == [0.5, 1.5]


def test_extract_roi_writes_expected_files(data_dir: Path, tmp_path: Path) -> None:
    from segbench import io

    out_dir = tmp_path / "roi"
    captured: dict[str, object] = {}

    def fake_load_roi_counts(data_dir: Path, cell_ids: pd.Series) -> object:
        captured["cell_ids"] = list(cell_ids)

        class FakeAdata:
            def write_h5ad(self, path: Path) -> None:
                Path(path).write_bytes(b"fake-h5ad")

        return FakeAdata()

    import segbench.io as io_module

    orig = io_module.load_roi_counts
    io_module.load_roi_counts = fake_load_roi_counts
    try:
        io.extract_roi(data_dir, out_dir, xmin=0, ymin=0, xmax=10, ymax=10, pixel_size=1.0)
    finally:
        io_module.load_roi_counts = orig

    for fname in (
        "dapi.tif",
        "he.tif",
        "transcripts.parquet",
        "cells.parquet",
        "cell_boundaries.parquet",
        "nucleus_boundaries.parquet",
        "cell_feature_matrix.h5ad",
    ):
        assert (out_dir / fname).exists(), fname

    assert captured["cell_ids"] == [b"cell-1", b"cell-2"]
