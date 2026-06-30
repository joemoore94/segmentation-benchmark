"""Run Mesmer whole-cell segmentation with eosin or density as membrane channel.

Wraps run_mesmer.py with the correct arguments for whole-cell compartment
and quantifies the resulting mask.

Usage::

    conda run -n mesmer python scripts/run_mesmer_wholecell.py --eosin
    conda run -n mesmer python scripts/run_mesmer_wholecell.py --density
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import tifffile

DATA_DIR = Path("data/processed/roi")
IMAGE_MPP = 0.2125


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--eosin", action="store_true")
    group.add_argument("--density", action="store_true")
    args = parser.parse_args()

    if args.eosin:
        membrane_file = "eosin.tif"
        out_dir = "mesmer_wholecell_eosin"
        key = "mesmer_wholecell_eosin"
    else:
        membrane_file = "density_heatmap.tif"
        out_dir = "mesmer_wholecell_density"
        key = "mesmer_wholecell_density"

    nuclear = tifffile.imread(DATA_DIR / "dapi.tif").astype(np.float32)
    membrane = tifffile.imread(DATA_DIR / membrane_file).astype(np.float32)
    print(f"Nuclear: {nuclear.shape}, Membrane ({membrane_file}): {membrane.shape}")

    from deepcell.applications import Mesmer

    image = np.stack([nuclear, membrane], axis=-1)[np.newaxis, ...]
    print(f"Input shape: {image.shape}, compartment=whole-cell, mpp={IMAGE_MPP}")

    app = Mesmer()
    labels = app.predict(image, image_mpp=IMAGE_MPP, compartment="whole-cell")
    mask = labels[0, ..., 0].astype(np.int32)
    print(f"Mesmer whole-cell done — {mask.max()} cells detected")

    out_path = DATA_DIR / out_dir
    out_path.mkdir(parents=True, exist_ok=True)
    mask_file = out_path / "mask.tif"
    tifffile.imwrite(mask_file, mask)
    print(f"Saved mask: {mask_file}")

    # Quantify — import segbench only if available (may not be in mesmer env)
    try:
        import pandas as pd
        from segbench.quantify import quantify_cells
        transcripts = pd.read_parquet(DATA_DIR / "transcripts_baysor.parquet")
        adata = quantify_cells(mask, transcripts, pixel_size=IMAGE_MPP)
        adata_path = DATA_DIR / f"adata_{key}.h5ad"
        adata.write_h5ad(adata_path)
        total_tx = int(adata.X.sum())
        median_tx = float(np.median(np.asarray(adata.X.sum(axis=1)).ravel()))
        print(f"  {adata.n_obs} cells, median {median_tx:.0f} tx/cell, "
              f"capture {total_tx / len(transcripts) * 100:.1f}%")
        print(f"  Saved: {adata_path.name}")
    except ImportError:
        print("  segbench not available in this env — run quantification separately")
        print(f"  Mask saved at {mask_file}, quantify with segbench env")


if __name__ == "__main__":
    main()
