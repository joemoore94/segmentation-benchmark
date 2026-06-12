"""Aggregate per-transcript detections into a per-cell x gene AnnData."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from skimage.measure import regionprops_table


def quantify_cells(
    masks: NDArray[np.integer],
    transcripts: pd.DataFrame,
    x_col: str = "x_location",
    y_col: str = "y_location",
    gene_col: str = "feature_name",
    pixel_size: float = 1.0,
    origin: tuple[float, float] = (0.0, 0.0),
) -> ad.AnnData:
    """Aggregate transcripts into a cell x gene count matrix using a label mask.

    Parameters
    ----------
    masks:
        2D integer label image (0 = background, >0 = cell id), as produced by
        a segmentation method.
    transcripts:
        Per-transcript table with physical (x, y) coordinates and a gene/feature
        name column.
    x_col, y_col:
        Column names for the transcript's physical x/y coordinates.
    gene_col:
        Column name for the transcript's gene/feature identity.
    pixel_size:
        Physical units per pixel (e.g. microns/pixel), used to convert
        transcript coordinates to mask pixel indices.
    origin:
        Physical (x, y) coordinate of the mask's pixel ``(0, 0)``, used to
        convert transcript coordinates to mask pixel indices.

    Returns
    -------
    AnnData with one observation per labeled cell (``obs_names`` = label ids
    as strings), ``var_names`` = gene names, ``X`` = transcript counts, and
    ``obs["area"]``, ``obs["centroid_y"]``, ``obs["centroid_x"]`` (in pixel
    coordinates) from :func:`skimage.measure.regionprops_table`. Cells with no
    assigned transcripts get all-zero rows.
    """
    col = ((transcripts[x_col].to_numpy() - origin[0]) / pixel_size).astype(np.int64)
    row = ((transcripts[y_col].to_numpy() - origin[1]) / pixel_size).astype(np.int64)

    in_bounds = (
        (row >= 0) & (row < masks.shape[0]) & (col >= 0) & (col < masks.shape[1])
    )
    cell_id = np.zeros(len(transcripts), dtype=masks.dtype)
    cell_id[in_bounds] = masks[row[in_bounds], col[in_bounds]]

    assigned = transcripts.loc[cell_id > 0].copy()
    assigned["cell_id"] = cell_id[cell_id > 0]

    counts = assigned.groupby(["cell_id", gene_col]).size().unstack(fill_value=0)

    props = pd.DataFrame(
        regionprops_table(masks, properties=("label", "area", "centroid"))
    ).set_index("label")
    props = props.rename(columns={"centroid-0": "centroid_y", "centroid-1": "centroid_x"})

    obs = props.reindex(props.index.union(counts.index))
    counts = counts.reindex(obs.index, fill_value=0)

    obs.index = obs.index.astype(str)
    counts.index = counts.index.astype(str)

    return ad.AnnData(
        X=counts.to_numpy(dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )
