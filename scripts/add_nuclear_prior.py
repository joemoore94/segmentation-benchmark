"""Add a nuclear-mask prior column to the Baysor transcript table.

Generalises add_cellpose_prior.py to work with any nuclear label mask.
For each transcript, looks up the label at that transcript's pixel location
and writes it as a new column. Transcripts outside all nuclei get 0.

Usage::

    conda run -n segbench python scripts/add_nuclear_prior.py stardist
    conda run -n segbench python scripts/add_nuclear_prior.py mesmer
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from segbench.io import PIXEL_SIZE

ROI_DIR = Path("data/processed/roi")

MASK_PATHS = {
    "cellpose": ROI_DIR / "masks_cellpose.tif",
    "stardist": ROI_DIR / "masks_stardist.tif",
    "mesmer":   ROI_DIR / "mesmer_out" / "mask.tif",
}


def main(method: str) -> None:
    if method not in MASK_PATHS:
        print(f"Unknown method: {method}. Choose from: {list(MASK_PATHS)}")
        sys.exit(1)

    transcripts = pd.read_csv(ROI_DIR / "transcripts_baysor.csv")
    print(f"loaded {len(transcripts)} transcripts")

    mask_path = MASK_PATHS[method]
    masks = tifffile.imread(mask_path)
    print(f"{mask_path.name}: {masks.shape}, {masks.max()} labels")

    col = (transcripts["x_location"].to_numpy() / PIXEL_SIZE).astype(np.int64)
    row = (transcripts["y_location"].to_numpy() / PIXEL_SIZE).astype(np.int64)
    in_bounds = (row >= 0) & (row < masks.shape[0]) & (col >= 0) & (col < masks.shape[1])

    prior_col = f"{method}_prior"
    prior = np.zeros(len(transcripts), dtype=np.int32)
    prior[in_bounds] = masks[row[in_bounds], col[in_bounds]]
    transcripts[prior_col] = prior

    n_assigned = (prior > 0).sum()
    print(f"{n_assigned} / {len(transcripts)} transcripts ({n_assigned / len(transcripts):.1%}) "
          f"fall within a {method} nucleus")

    out_path = ROI_DIR / f"transcripts_baysor_{method}_prior.csv"
    transcripts.to_csv(out_path, index=False)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python add_nuclear_prior.py <cellpose|stardist|mesmer>")
        sys.exit(1)
    main(sys.argv[1])
