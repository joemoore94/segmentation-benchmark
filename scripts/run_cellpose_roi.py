"""Run CellPose on the extracted ROI's DAPI image and save the label mask.

Usage::

    conda run -n segbench python scripts/run_cellpose_roi.py
"""

from __future__ import annotations

import time
from pathlib import Path

import tifffile

from segbench.segmentation.cellpose_run import run_cellpose

ROI_DIR = Path("data/processed/roi")


def main() -> None:
    dapi = tifffile.imread(ROI_DIR / "dapi.tif")
    print(f"DAPI image: {dapi.shape} {dapi.dtype}")

    t0 = time.time()
    masks = run_cellpose(dapi, gpu=True)
    print(f"CellPose done in {time.time() - t0:.1f}s, {masks.max()} cells")

    out_path = ROI_DIR / "masks_cellpose.tif"
    tifffile.imwrite(out_path, masks)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
