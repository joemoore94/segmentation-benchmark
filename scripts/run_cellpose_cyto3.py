"""Run CellPose cyto3 whole-cell model on DAPI and optional second channel.

Three modes:
  --dapi-only     DAPI as single grayscale input (no cytoplasm channel)
  --eosin         DAPI (nuclear) + eosin from H&E (cytoplasm)
  --density       DAPI (nuclear) + transcript density heatmap (cytoplasm)

Each mode produces a whole-cell label mask and a quantified AnnData.

Usage::

    conda run -n segbench python scripts/run_cellpose_cyto3.py --dapi-only
    conda run -n segbench python scripts/run_cellpose_cyto3.py --eosin
    conda run -n segbench python scripts/run_cellpose_cyto3.py --density
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells
from segbench.segmentation.cellpose_run import run_cellpose

ROI_DIR = Path("data/processed/roi")


def load_and_normalize(path: Path) -> np.ndarray:
    img = tifffile.imread(path).astype(np.float32)
    p995 = np.percentile(img[img > 0], 99.5) if (img > 0).any() else 1.0
    return np.clip(img, 0, p995) / p995


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dapi-only", action="store_true")
    group.add_argument("--eosin", action="store_true")
    group.add_argument("--density", action="store_true")
    args = parser.parse_args()

    dapi = load_and_normalize(ROI_DIR / "dapi.tif")
    print(f"DAPI: {dapi.shape}, range [{dapi.min():.3f}, {dapi.max():.3f}]")

    if args.dapi_only:
        image = dapi
        channels = [0, 0]
        suffix = "cyto3"
        key = "cellpose_cyto3"
    elif args.eosin:
        eosin = load_and_normalize(ROI_DIR / "eosin.tif")
        print(f"Eosin: {eosin.shape}, range [{eosin.min():.3f}, {eosin.max():.3f}]")
        image = np.stack([eosin, dapi], axis=-1)
        channels = [1, 2]
        suffix = "cyto3_eosin"
        key = "cellpose_cyto3_eosin"
    else:
        density = load_and_normalize(ROI_DIR / "density_heatmap.tif")
        print(f"Density: {density.shape}, range [{density.min():.3f}, {density.max():.3f}]")
        image = np.stack([density, dapi], axis=-1)
        channels = [1, 2]
        suffix = "cyto3_density"
        key = "cellpose_cyto3_density"

    print(f"\nRunning CellPose cyto3 ({suffix})...")
    print(f"  Input shape: {image.shape}, channels={channels}")

    t0 = time.time()
    masks = run_cellpose(image, model_type="cyto3", channels=channels, gpu=False)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s — {masks.max()} cells detected")

    mask_path = ROI_DIR / f"masks_cellpose_{suffix}.tif"
    tifffile.imwrite(mask_path, masks)
    print(f"  Saved mask: {mask_path.name}")

    print("Quantifying transcripts...")
    transcripts = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    adata = quantify_cells(masks, transcripts, pixel_size=PIXEL_SIZE)
    adata_path = ROI_DIR / f"adata_{key}.h5ad"
    adata.write_h5ad(adata_path)

    total_tx = int(adata.X.sum())
    median_tx = float(np.median(np.asarray(adata.X.sum(axis=1)).ravel()))
    print(f"  {adata.n_obs} cells, median {median_tx:.0f} tx/cell, "
          f"capture {total_tx / len(transcripts) * 100:.1f}%")
    print(f"  Saved: {adata_path.name}")


if __name__ == "__main__":
    main()
