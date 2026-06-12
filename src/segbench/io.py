"""Loading the Xenium bundle and cropping it down to an ROI."""

from __future__ import annotations

from pathlib import Path

from spatialdata import SpatialData, bounding_box_query


def load_xenium(data_dir: str | Path) -> SpatialData:
    """Load a 10x Xenium output bundle (``data/raw``) as a SpatialData object.

    Only the ``morphology_focus`` image is read (not ``morphology_mip``, which
    we don't download), and aligned H&E/IF images are skipped since this
    bundle has no alignment matrix for them.
    """
    import spatialdata_io

    return spatialdata_io.xenium(
        data_dir,
        morphology_mip=False,
        aligned_images=False,
    )


def crop_roi(
    sdata: SpatialData,
    xmin: float,
    ymin: float,
    xmax: float,
    ymax: float,
    target_coordinate_system: str = "global",
) -> SpatialData:
    """Crop ``sdata`` to the rectangular ROI ``[xmin, xmax] x [ymin, ymax]``.

    Coordinates are in the units of ``target_coordinate_system`` (microns,
    for the Xenium "global" coordinate system).
    """
    cropped = bounding_box_query(
        sdata,
        axes=("x", "y"),
        min_coordinate=[xmin, ymin],
        max_coordinate=[xmax, ymax],
        target_coordinate_system=target_coordinate_system,
    )
    if cropped is None or not cropped.elements_paths_in_memory():
        raise ValueError(
            f"ROI [{xmin}, {xmax}] x [{ymin}, {ymax}] does not intersect any data"
        )
    return cropped
