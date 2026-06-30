"""Watershed expansion from 10x Ranger nuclear seeds on DAPI gradient.

Uses the 10x Ranger nuclear mask as seeds and the inverted DAPI as the
landscape. Watershed basins expand from each nucleus until they meet at
intensity valleys between cells.

Reads:  data/processed/roi/masks_10x_ranger.tif (or rasterized equivalent)
        data/processed/roi/dapi.tif
        data/processed/roi/transcripts_baysor.parquet
Writes: data/processed/roi/masks_watershed_10x.tif
        data/processed/roi/adata_watershed_10x.h5ad

Usage::

    conda run -n segbench python scripts/run_watershed_expansion.py
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter
from skimage.segmentation import watershed

from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells

ROI_DIR = Path("data/processed/roi")


def main() -> None:
    print("Loading DAPI and 10x Ranger nuclear mask...")
    dapi = tifffile.imread(ROI_DIR / "dapi.tif").astype(np.float32)

    ranger_mask_path = ROI_DIR / "masks_10x_ranger.tif"
    if not ranger_mask_path.exists():
        print(f"  {ranger_mask_path.name} not found, trying to rasterize from adata...")
        import anndata as ad
        from skimage.draw import polygon as draw_polygon

        adata_ranger = ad.read_h5ad(ROI_DIR / "adata_10x_ranger.h5ad")
        mask = np.zeros(dapi.shape[:2], dtype=np.int32)
        cx = adata_ranger.obs["centroid_x"].values
        cy = adata_ranger.obs["centroid_y"].values
        for i, (x, y) in enumerate(zip(cx, cy), 1):
            r, c = int(y / PIXEL_SIZE), int(x / PIXEL_SIZE)
            if 0 <= r < mask.shape[0] and 0 <= c < mask.shape[1]:
                mask[max(0, r-3):r+4, max(0, c-3):c+4] = i
        tifffile.imwrite(ranger_mask_path, mask)
        print(f"  Rasterized {mask.max()} centroids as 7×7 px seeds")
    else:
        mask = tifffile.imread(ranger_mask_path)

    print(f"  DAPI: {dapi.shape}, Ranger mask: {mask.shape}, {mask.max()} seeds")

    dapi_smooth = gaussian_filter(dapi, sigma=2)
    landscape = -dapi_smooth

    print("Running watershed...")
    t0 = time.time()
    ws_mask = watershed(landscape, markers=mask, compactness=0.01)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s — {ws_mask.max()} cells")

    out_mask = ROI_DIR / "masks_watershed_10x.tif"
    tifffile.imwrite(out_mask, ws_mask.astype(np.int32))
    print(f"  Saved {out_mask.name}")

    print("Quantifying transcripts...")
    transcripts = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    adata = quantify_cells(ws_mask.astype(np.int32), transcripts, pixel_size=PIXEL_SIZE)
    adata_path = ROI_DIR / "adata_watershed_10x.h5ad"
    adata.write_h5ad(adata_path)

    total_tx = int(adata.X.sum())
    median_tx = float(np.median(np.asarray(adata.X.sum(axis=1)).ravel()))
    print(f"  {adata.n_obs} cells, median {median_tx:.0f} tx/cell, "
          f"capture {total_tx / len(transcripts) * 100:.1f}%")
    print(f"  Saved: {adata_path.name}")


if __name__ == "__main__":
    main()
