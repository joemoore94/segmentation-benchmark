"""CellPose segmentation wrapper."""

from __future__ import annotations

from typing import cast

import numpy as np
from numpy.typing import NDArray


def run_cellpose(
    image: NDArray[np.floating],
    model_type: str = "nuclei",
    channels: list[int] | None = None,
    diameter: float | None = None,
    gpu: bool = True,
) -> NDArray[np.int32]:
    """Segment cells in ``image`` with CellPose and return an integer label mask.

    Parameters
    ----------
    image:
        2D (Y, X) or multi-channel (Y, X, C) image array.
    model_type:
        CellPose classical pretrained model name (e.g. ``"nuclei"``, ``"cyto2"``,
        ``"cyto3"``). Default is ``"nuclei"``, appropriate for DAPI input. These
        are the lightweight U-Net models from CellPose <4 (the CellPose 4.x
        Cellpose-SAM models are CPU-prohibitive for whole-slide ROIs).
    channels:
        CellPose channel spec, e.g. ``[0, 0]`` for grayscale, or
        ``[cytoplasm_channel, nucleus_channel]`` for two-channel input.
        Defaults to ``[0, 0]``.
    diameter:
        Expected cell diameter in pixels. ``None`` lets CellPose estimate it.
    gpu:
        Whether to request GPU acceleration (CUDA) (CellPose falls back to CPU
        if unavailable).

    Returns
    -------
    Integer label mask, same (Y, X) shape as ``image``, where 0 is background
    and each cell has a unique positive integer id.
    """
    from cellpose import models

    if channels is None:
        channels = [0, 0]

    model = models.CellposeModel(gpu=gpu, model_type=model_type)
    masks, _flows, _styles = model.eval(image, diameter=diameter, channels=channels)
    return cast(NDArray[np.int32], masks.astype(np.int32))
