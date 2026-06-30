"""Build a transcript density heatmap on the same pixel grid as the DAPI image.

Bins transcript coordinates into a 2D histogram matching the DAPI resolution
(0.2125 µm/pixel, 9412 × 9412 px), then smooths with a Gaussian kernel to
produce a continuous density surface suitable as a pseudo-membrane channel.

Reads:  data/processed/roi/transcripts_baysor.parquet
        data/processed/roi/dapi.tif (for grid dimensions)
Writes: data/processed/roi/density_heatmap.tif

Usage::

    conda run -n segbench python scripts/build_density_heatmap.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter

from segbench.io import PIXEL_SIZE

ROI_DIR = Path("data/processed/roi")
SIGMA_PX = 3  # ~0.6 µm, preserves subcellular resolution


def main() -> None:
    dapi = tifffile.imread(ROI_DIR / "dapi.tif")
    h, w = dapi.shape[:2]
    print(f"Grid: {h} × {w} px (PIXEL_SIZE={PIXEL_SIZE} µm)")

    tx = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    print(f"Transcripts: {len(tx)}")

    col = (tx["x_location"].values / PIXEL_SIZE).astype(np.int32)
    row = (tx["y_location"].values / PIXEL_SIZE).astype(np.int32)

    valid = (row >= 0) & (row < h) & (col >= 0) & (col < w)
    row, col = row[valid], col[valid]
    print(f"  {valid.sum()} transcripts within grid ({(~valid).sum()} out of bounds)")

    density = np.zeros((h, w), dtype=np.float32)
    np.add.at(density, (row, col), 1)
    print(f"  Raw counts: max={density.max():.0f}, nonzero={np.count_nonzero(density)}")

    density = gaussian_filter(density, sigma=SIGMA_PX)
    print(f"  Smoothed (σ={SIGMA_PX} px): max={density.max():.2f}")

    density_uint16 = (density / density.max() * 65535).astype(np.uint16)
    out = ROI_DIR / "density_heatmap.tif"
    tifffile.imwrite(out, density_uint16)
    print(f"Saved {out.name} ({density_uint16.shape}, uint16)")


if __name__ == "__main__":
    main()
