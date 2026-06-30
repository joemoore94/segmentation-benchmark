"""Build the missing factorial adata files for StarDist, Mesmer, and 10x Ranger.

Completes the N-detector × N-expansion factorial grid by adding geometric
expansion (10µm, 20µm) and watershed for the three detectors that were missing
these strategies:

  StarDist   + geometric 10µm  → adata_stardist_exp10um.h5ad
  StarDist   + geometric 20µm  → adata_stardist_exp20um.h5ad
  StarDist   + watershed       → adata_watershed_stardist.h5ad
  Mesmer     + geometric 10µm  → adata_mesmer_exp10um.h5ad
  Mesmer     + geometric 20µm  → adata_mesmer_exp20um.h5ad
  Mesmer     + watershed       → adata_watershed_mesmer.h5ad
  10x Ranger + geometric 10µm  → adata_10x_ranger_exp10um.h5ad
  10x Ranger + geometric 20µm  → adata_10x_ranger_exp20um.h5ad

Nuclear masks used as seeds / bases:
  StarDist   → data/processed/roi/masks_stardist.tif
  Mesmer     → data/processed/roi/mesmer_out/mask.tif  (may be 4-D, auto-squeezed)
  10x Ranger → data/processed/roi/masks_10x_ranger.tif

Watershed uses DAPI (data/processed/roi/dapi.tif) as the landscape, identical
to run_watershed_expansion.py but parameterised over detector.

Usage::

    conda run -n segbench python scripts/build_missing_factorial_adatas.py
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter
from skimage.segmentation import expand_labels, watershed

from segbench.constants import TOTAL_TRANSCRIPTS_FULL_ROI
from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells

ROI_DIR = Path("data/processed/roi")
EXPANSION_DISTANCES_UM = [10.0, 20.0]


def load_mask(path: Path) -> np.ndarray:
    raw = tifffile.imread(path)
    if raw.ndim == 4:
        raw = raw[0, ..., 0]
    return raw.astype(np.int32)


def build_expansion(
    mask: np.ndarray,
    transcripts: pd.DataFrame,
    distance_um: float,
) -> tuple[int, float, object]:
    radius_px = distance_um / PIXEL_SIZE
    expanded = expand_labels(mask, distance=radius_px)
    adata = quantify_cells(expanded, transcripts, pixel_size=PIXEL_SIZE)
    n_tx = int(np.asarray(adata.X).sum())
    capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI * 100
    return adata, n_tx, capture


def build_watershed(
    mask: np.ndarray,
    dapi: np.ndarray,
    transcripts: pd.DataFrame,
    out_mask_path: Path,
) -> tuple[object, int, float]:
    dapi_smooth = gaussian_filter(dapi.astype(np.float32), sigma=2)
    landscape = -dapi_smooth

    t0 = time.time()
    ws_mask = watershed(landscape, markers=mask, compactness=0.01).astype(np.int32)
    print(f"    watershed done in {time.time()-t0:.1f}s — {ws_mask.max()} cells")

    tifffile.imwrite(out_mask_path, ws_mask)

    adata = quantify_cells(ws_mask, transcripts, pixel_size=PIXEL_SIZE)
    n_tx = int(np.asarray(adata.X).sum())
    capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI * 100
    return adata, n_tx, capture


def main() -> None:
    transcripts = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    print(f"Loaded {len(transcripts)} transcripts")

    dapi = tifffile.imread(ROI_DIR / "dapi.tif").astype(np.float32)

    detectors = [
        ("stardist",   ROI_DIR / "masks_stardist.tif"),
        ("mesmer",     ROI_DIR / "mesmer_out" / "mask.tif"),
        ("10x_ranger", ROI_DIR / "masks_10x_ranger.tif"),
    ]

    for det_name, mask_path in detectors:
        print(f"\n{'='*60}")
        print(f"Detector: {det_name}  (mask: {mask_path})")
        mask = load_mask(mask_path)
        print(f"  mask shape: {mask.shape}, {mask.max()} seeds")

        for dist_um in EXPANSION_DISTANCES_UM:
            key = f"{det_name}_exp{dist_um:.0f}um"
            out = ROI_DIR / f"adata_{key}.h5ad"
            if out.exists():
                print(f"  [skip] {out.name} already exists")
                continue
            print(f"\n  --- Geometric {dist_um:.0f}µm ---")
            adata, n_tx, capture = build_expansion(mask, transcripts, dist_um)
            adata.write_h5ad(out)
            print(f"  {adata.n_obs} cells, {n_tx} transcripts ({capture:.1f}% capture) → {out.name}")

        if det_name == "10x_ranger":
            continue

        ws_key = f"watershed_{det_name}"
        ws_out = ROI_DIR / f"adata_{ws_key}.h5ad"
        ws_mask_out = ROI_DIR / f"masks_{ws_key}.tif"
        if ws_out.exists():
            print(f"\n  [skip] {ws_out.name} already exists")
        else:
            print(f"\n  --- Watershed ---")
            adata, n_tx, capture = build_watershed(mask, dapi, transcripts, ws_mask_out)
            adata.write_h5ad(ws_out)
            print(f"  {adata.n_obs} cells, {n_tx} transcripts ({capture:.1f}% capture) → {ws_out.name}")

    print("\nAll missing factorial adatas built.")


if __name__ == "__main__":
    main()
