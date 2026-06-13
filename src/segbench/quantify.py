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
    ``obs["area"]`` (in pixel\N{SUPERSCRIPT TWO}) and ``obs["centroid_x"]``/
    ``obs["centroid_y"]`` converted to the same physical coordinate system as
    ``transcripts`` (via ``pixel_size``/``origin``), so cells from different
    methods can be matched by centroid. Cells with no assigned transcripts get
    all-zero rows.
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
    props["centroid_x"] = props["centroid_x"] * pixel_size + origin[0]
    props["centroid_y"] = props["centroid_y"] * pixel_size + origin[1]

    obs = props.reindex(props.index.union(counts.index))
    counts = counts.reindex(obs.index, fill_value=0)

    obs.index = obs.index.astype(str)
    counts.index = counts.index.astype(str)

    return ad.AnnData(
        X=counts.to_numpy(dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )


def quantify_baysor(
    segmentation: pd.DataFrame,
    gene_col: str = "gene",
    cell_col: str = "cell",
    x_col: str = "x",
    y_col: str = "y",
) -> ad.AnnData:
    """Aggregate Baysor's per-molecule segmentation output into a cell x gene AnnData.

    Unlike :func:`quantify_cells`, Baysor assigns each transcript to a cell
    directly (no pixel mask): cells are defined by ``segmentation[cell_col]``,
    and molecules with a missing (NaN) cell id are unassigned/noise and dropped.

    Returns
    -------
    AnnData with one observation per Baysor cell id, ``var_names`` = gene
    names, ``X`` = transcript counts, and ``obs["n_transcripts"]``,
    ``obs["centroid_x"]``, ``obs["centroid_y"]`` (the mean molecule position
    per cell, in the same physical units as ``x_col``/``y_col``). Baysor has
    no pixel mask, so cell "size" is given as ``n_transcripts`` rather than
    a pixel area.
    """
    assigned = segmentation.dropna(subset=[cell_col])

    counts = assigned.groupby([cell_col, gene_col]).size().unstack(fill_value=0)

    obs = assigned.groupby(cell_col)[[x_col, y_col]].mean().rename(
        columns={x_col: "centroid_x", y_col: "centroid_y"}
    )
    obs["n_transcripts"] = assigned.groupby(cell_col).size()
    obs = obs.reindex(counts.index)

    obs.index = obs.index.astype(str)
    counts.index = counts.index.astype(str)

    return ad.AnnData(
        X=counts.to_numpy(dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )
