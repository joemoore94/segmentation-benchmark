"""Convert BIDCell segmentation mask to a per-cell AnnData.

BIDCell outputs a pixel-level label mask (*_connected.tif). This script
uses the same quantify_cells() function as CellPose/StarDist/Mesmer to
assign transcripts to cells and build a gene x cell count matrix.

Usage::

    conda run -n segbench python scripts/build_bidcell_adata.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import tifffile

from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells

ROI_DIR = Path("data/processed/roi")
BIDCELL_OUT = ROI_DIR / "bidcell_out"


def main() -> None:
    mask_candidates = list(BIDCELL_OUT.rglob("*_connected.tif"))
    if not mask_candidates:
        raise FileNotFoundError(
            f"No *_connected.tif mask found in {BIDCELL_OUT}. "
            "Run scripts/run_bidcell.py first."
        )
    mask_path = mask_candidates[0]
    print(f"Loading BIDCell mask: {mask_path}")
    masks = tifffile.imread(mask_path)
    print(f"  Mask shape: {masks.shape}, max label: {masks.max()}")

    print("Loading transcripts...")
    transcripts = pd.read_parquet(ROI_DIR / "transcripts.parquet")
    print(f"  {len(transcripts)} transcripts")

    print("Quantifying cells...")
    adata = quantify_cells(
        masks, transcripts,
        pixel_size=PIXEL_SIZE,
        origin=(0.0, 0.0),
    )

    out_path = ROI_DIR / "adata_bidcell.h5ad"
    adata.write_h5ad(out_path)
    print(f"Wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{int(adata.X.sum())} assigned transcripts")


if __name__ == "__main__":
    main()
