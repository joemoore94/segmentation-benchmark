"""Add a CellPose-nucleus prior-segmentation column to the Baysor transcript table.

For each transcript in ``transcripts_baysor.csv``, looks up the CellPose label
(``masks_cellpose.tif``, 0 = background) at that transcript's pixel location and
writes it as a new ``cellpose_prior`` column. This is the "hybrid" input for
Baysor's ``--prior-segmentation`` mode (``:cellpose_prior`` column, see
``run_baysor.sh``): Baysor's transcript-density model expands/refines cells
seeded from CellPose's nuclear segmentation, with 0 meaning "no prior label"
(Baysor's default ``unassigned-prior-label``).

Usage::

    conda run -n segbench python scripts/add_cellpose_prior.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from segbench.io import PIXEL_SIZE

ROI_DIR = Path("data/processed/roi")


def main() -> None:
    transcripts = pd.read_csv(ROI_DIR / "transcripts_baysor.csv")
    print(f"loaded {len(transcripts)} transcripts")

    masks = tifffile.imread(ROI_DIR / "masks_cellpose.tif")
    print(f"masks_cellpose.tif: {masks.shape}, {masks.max()} labels")

    col = (transcripts["x_location"].to_numpy() / PIXEL_SIZE).astype(np.int64)
    row = (transcripts["y_location"].to_numpy() / PIXEL_SIZE).astype(np.int64)
    in_bounds = (row >= 0) & (row < masks.shape[0]) & (col >= 0) & (col < masks.shape[1])

    prior = np.zeros(len(transcripts), dtype=np.int32)
    prior[in_bounds] = masks[row[in_bounds], col[in_bounds]]
    transcripts["cellpose_prior"] = prior

    n_assigned = (prior > 0).sum()
    print(f"{n_assigned} / {len(transcripts)} transcripts ({n_assigned / len(transcripts):.1%}) "
          f"fall within a CellPose nucleus")

    out_path = ROI_DIR / "transcripts_baysor_prior.csv"
    transcripts.to_csv(out_path, index=False)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
