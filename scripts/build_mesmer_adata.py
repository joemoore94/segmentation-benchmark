"""Build AnnData for Mesmer (DeepCell) segmentation from the Docker output mask.

The vanvalenlab/deepcell-applications container writes a 4-D label image
(1, H, W, 1); this script squeezes it to (H, W) before quantifying transcripts.

Usage::

    conda run -n segbench python scripts/build_mesmer_adata.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells

ROI_DIR = Path("data/processed/roi")
TOTAL_TRANSCRIPTS_FULL_ROI = 3_392_051


def main() -> None:
    mask_raw = tifffile.imread(ROI_DIR / "mesmer_out" / "mask.tif")
    print(f"raw mask shape: {mask_raw.shape}, dtype: {mask_raw.dtype}")
    mask = mask_raw[0, ..., 0] if mask_raw.ndim == 4 else mask_raw
    print(f"squeezed mask: {mask.shape}, {mask.max()} cells")

    transcripts = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    print(f"loaded {len(transcripts)} transcripts")

    adata = quantify_cells(mask, transcripts, pixel_size=PIXEL_SIZE)

    n_tx = int(np.asarray(adata.X).sum())
    capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI
    print(f"Mesmer: {adata.n_obs} cells, {n_tx} transcripts ({capture:.1%} capture)")

    out = ROI_DIR / "adata_mesmer.h5ad"
    adata.write_h5ad(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
