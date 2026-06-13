"""Run StarDist on the extracted ROI's DAPI image and save the label mask.

Usage::

    conda run -n segbench python scripts/run_stardist_roi.py
"""

from __future__ import annotations

import time
from pathlib import Path

import tifffile

from segbench.segmentation.stardist_run import run_stardist

ROI_DIR = Path("data/processed/roi")


def main() -> None:
    t0 = time.time()
    out_dir = run_stardist(ROI_DIR, "dapi.tif", "stardist_out")
    print(f"StarDist done in {time.time() - t0:.1f}s")

    masks = tifffile.imread(out_dir / "mask.tif")
    out_path = ROI_DIR / "masks_stardist.tif"
    tifffile.imwrite(out_path, masks)
    print(f"wrote {out_path}, {masks.max()} cells")


if __name__ == "__main__":
    main()
