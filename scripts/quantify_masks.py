"""Quantify each segmentation method's label mask into a cell x gene AnnData.

All three methods are quantified against the same qv>=20, non-control
transcript set (``transcripts_baysor.parquet``) so cell counts and
per-cell expression are directly comparable.

Usage::

    conda run -n segbench python scripts/quantify_masks.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import tifffile

from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells

ROI_DIR = Path("data/processed/roi")


def main() -> None:
    transcripts = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    print(f"transcripts: {len(transcripts)} rows")

    cellpose_masks = tifffile.imread(ROI_DIR / "masks_cellpose.tif")
    adata_cellpose = quantify_cells(cellpose_masks, transcripts, pixel_size=PIXEL_SIZE)
    adata_cellpose.write_h5ad(ROI_DIR / "adata_cellpose.h5ad")
    print(f"cellpose: {adata_cellpose.n_obs} cells -> adata_cellpose.h5ad")

    stardist_masks = tifffile.imread(ROI_DIR / "masks_stardist.tif")
    adata_stardist = quantify_cells(stardist_masks, transcripts, pixel_size=PIXEL_SIZE)
    adata_stardist.write_h5ad(ROI_DIR / "adata_stardist.h5ad")
    print(f"stardist: {adata_stardist.n_obs} cells -> adata_stardist.h5ad")

    mesmer_masks = tifffile.imread(ROI_DIR / "mesmer_out" / "mask.tif")
    adata_mesmer = quantify_cells(mesmer_masks, transcripts, pixel_size=PIXEL_SIZE)
    adata_mesmer.write_h5ad(ROI_DIR / "adata_mesmer.h5ad")
    print(f"mesmer: {adata_mesmer.n_obs} cells -> adata_mesmer.h5ad")


if __name__ == "__main__":
    main()
