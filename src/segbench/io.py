"""Extract a rectangular ROI from a 10x Xenium output bundle.

We work directly with the per-file downloads (OME-TIFFs, parquet tables,
``cell_feature_matrix.h5``) rather than ``spatialdata_io.xenium()``, which
expects the Xenium Explorer ``*.zarr.zip`` bundles (``cells.zarr.zip``,
``transcripts.zarr.zip``, ...) that aren't part of this download.

Two image-format gotchas drive the implementation:

- ``morphology_focus.ome.tif`` (DAPI) is JPEG2000-compressed. pyvips/libtiff
  silently decodes it to all-zeros, so it's read via tifffile's zarr store
  instead.
- ``he_image.ome.tif`` stores its 3 channels (Hematoxylin/Eosin/Residual) as
  separate pages, not as bands of one page. pyvips reads each page
  separately and they're band-joined into an RGB-like array.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyvips
import scanpy as sc
import tifffile
import zarr
from anndata import AnnData
from numpy.typing import NDArray

PIXEL_SIZE = 0.2125  # microns per pixel (experiment.xenium)


def _to_px(coord_um: float, pixel_size: float) -> int:
    return int(round(coord_um / pixel_size))


def load_roi_dapi(
    data_dir: str | Path,
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    pixel_size: float = PIXEL_SIZE,
) -> NDArray[np.uint16]:
    """Crop the DAPI ``morphology_focus`` image to an ROI given in microns."""
    x0, y0 = _to_px(xmin, pixel_size), _to_px(ymin, pixel_size)
    x1, y1 = _to_px(xmax, pixel_size), _to_px(ymax, pixel_size)
    store = tifffile.imread(Path(data_dir) / "morphology_focus.ome.tif", aszarr=True)
    try:
        z = zarr.open(store, mode="r")
        # Pyramidal OME-TIFFs open as a Group of per-level arrays ("0" = full res);
        # a single-resolution TIFF opens as an Array directly.
        if isinstance(z, zarr.Group):
            arr: zarr.Array = z["0"]  # type: ignore[assignment]
        else:
            arr = z
        return np.asarray(arr[y0:y1, x0:x1])
    finally:
        store.close()


def load_roi_he(
    data_dir: str | Path,
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    pixel_size: float = PIXEL_SIZE,
) -> NDArray[np.uint8]:
    """Crop the deconvolved H&E image to an ROI given in microns.

    Returns an ``(h, w, 3)`` array with Hematoxylin/Eosin/Residual channels.
    """
    x0, y0 = _to_px(xmin, pixel_size), _to_px(ymin, pixel_size)
    x1, y1 = _to_px(xmax, pixel_size), _to_px(ymax, pixel_size)
    path = str(Path(data_dir) / "he_image.ome.tif")
    pages = [
        pyvips.Image.new_from_file(path, page=i).crop(x0, y0, x1 - x0, y1 - y0)
        for i in range(3)
    ]
    rgb = pages[0].bandjoin(pages[1:])
    return np.ndarray(
        buffer=rgb.write_to_memory(),
        dtype=np.uint8,
        shape=(rgb.height, rgb.width, rgb.bands),
    )


def _crop_points(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
) -> pd.DataFrame:
    mask = df[x_col].between(xmin, xmax) & df[y_col].between(ymin, ymax)
    cropped = df.loc[mask].copy()
    cropped[x_col] -= xmin
    cropped[y_col] -= ymin
    return cropped


def load_roi_transcripts(
    data_dir: str | Path, xmin: float, ymin: float, xmax: float, ymax: float
) -> pd.DataFrame:
    """Transcripts within the ROI, with coordinates shifted to the ROI origin."""
    df = pd.read_parquet(Path(data_dir) / "transcripts.parquet")
    return _crop_points(df, "x_location", "y_location", xmin, ymin, xmax, ymax)


def load_roi_cells(
    data_dir: str | Path, xmin: float, ymin: float, xmax: float, ymax: float
) -> pd.DataFrame:
    """10x cell table within the ROI, with centroids shifted to the ROI origin."""
    df = pd.read_parquet(Path(data_dir) / "cells.parquet")
    return _crop_points(df, "x_centroid", "y_centroid", xmin, ymin, xmax, ymax)


def load_roi_boundaries(
    data_dir: str | Path,
    cell_ids: pd.Series,
    xmin: float,
    ymin: float,
    kind: str = "cell",
) -> pd.DataFrame:
    """Cell or nucleus boundary polygons for ``cell_ids``, shifted to the ROI origin."""
    fname = "cell_boundaries.parquet" if kind == "cell" else "nucleus_boundaries.parquet"
    df = pd.read_parquet(Path(data_dir) / fname)
    cropped = df.loc[df["cell_id"].isin(cell_ids)].copy()
    cropped["vertex_x"] -= xmin
    cropped["vertex_y"] -= ymin
    return cropped


def load_roi_counts(data_dir: str | Path, cell_ids: pd.Series) -> AnnData:
    """10x reference cell x gene counts, subset to ``cell_ids``."""
    adata = sc.read_10x_h5(Path(data_dir) / "cell_feature_matrix.h5")
    barcodes = cell_ids.apply(lambda b: b.decode() if isinstance(b, bytes) else b)
    return adata[adata.obs_names.isin(barcodes)].copy()


def extract_roi(
    data_dir: str | Path,
    out_dir: str | Path,
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    pixel_size: float = PIXEL_SIZE,
) -> None:
    """Crop every input the segmentation methods need to an ROI and write to ``out_dir``.

    ROI bounds (``xmin``, ``ymin``, ``xmax``, ``ymax``) are in microns, in the
    Xenium "global" coordinate system. All written tables have coordinates
    shifted so the ROI's top-left corner is the origin.
    """
    data_dir = Path(data_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dapi = load_roi_dapi(data_dir, xmin, ymin, xmax, ymax, pixel_size)
    tifffile.imwrite(out_dir / "dapi.tif", dapi)

    he = load_roi_he(data_dir, xmin, ymin, xmax, ymax, pixel_size)
    tifffile.imwrite(out_dir / "he.tif", he)

    transcripts = load_roi_transcripts(data_dir, xmin, ymin, xmax, ymax)
    transcripts.to_parquet(out_dir / "transcripts.parquet")

    cells = load_roi_cells(data_dir, xmin, ymin, xmax, ymax)
    cells.to_parquet(out_dir / "cells.parquet")

    cell_ids = cells["cell_id"]
    load_roi_boundaries(data_dir, cell_ids, xmin, ymin, kind="cell").to_parquet(
        out_dir / "cell_boundaries.parquet"
    )
    load_roi_boundaries(data_dir, cell_ids, xmin, ymin, kind="nucleus").to_parquet(
        out_dir / "nucleus_boundaries.parquet"
    )

    counts = load_roi_counts(data_dir, cell_ids)
    counts.write_h5ad(out_dir / "cell_feature_matrix.h5ad")
