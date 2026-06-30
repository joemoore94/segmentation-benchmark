"""Extract eosin and hematoxylin channels from registered H&E image.

Color-deconvolves the H&E using Ruifrok & Johnston stain vectors (skimage
rgb2hed) and saves each channel as a standalone TIFF for downstream use as
input to whole-cell segmentation models (Cellpose cyto3, Mesmer).

Also generates a registration-verification overlay of DAPI vs. hematoxylin
to check for systematic shifts at single-cell resolution.

Reads:  data/processed/roi/he.tif
        data/processed/roi/dapi.tif
Writes: data/processed/roi/eosin.tif
        data/processed/roi/hematoxylin.tif
        results/figures/registration_overlay.png

Usage::

    conda run -n segbench python scripts/build_eosin_channel.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.color import rgb2hed

from segbench.style import apply_style

ROI_DIR = Path("data/processed/roi")
FIGURES = Path("results/figures")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)

    print("Loading H&E and DAPI...")
    he = tifffile.imread(ROI_DIR / "he.tif")
    dapi = tifffile.imread(ROI_DIR / "dapi.tif")
    print(f"  H&E shape: {he.shape}, dtype: {he.dtype}")
    print(f"  DAPI shape: {dapi.shape}, dtype: {dapi.dtype}")

    print("Color deconvolution (rgb2hed)...")
    hed = rgb2hed(he)
    hem_raw = hed[:, :, 0]
    eos_raw = hed[:, :, 1]

    hem_clipped = np.clip(hem_raw, 0, np.percentile(hem_raw[hem_raw > 0], 99.5))
    eos_clipped = np.clip(eos_raw, 0, np.percentile(eos_raw[eos_raw > 0], 99.5))

    hem_uint16 = (hem_clipped / hem_clipped.max() * 65535).astype(np.uint16)
    eos_uint16 = (eos_clipped / eos_clipped.max() * 65535).astype(np.uint16)

    tifffile.imwrite(ROI_DIR / "hematoxylin.tif", hem_uint16)
    tifffile.imwrite(ROI_DIR / "eosin.tif", eos_uint16)
    print(f"  Saved hematoxylin.tif ({hem_uint16.shape}, uint16)")
    print(f"  Saved eosin.tif ({eos_uint16.shape}, uint16)")

    # --- Registration overlay: DAPI (green) + hematoxylin (magenta) ---
    print("Generating registration overlay...")
    dapi_norm = dapi.astype(np.float32)
    dapi_norm = np.clip(dapi_norm, 0, np.percentile(dapi_norm[dapi_norm > 0], 99.5))
    dapi_norm /= dapi_norm.max()

    hem_norm = hem_clipped / hem_clipped.max()

    # Full ROI overview
    apply_style(scatter=True)
    fig, axes = plt.subplots(1, 3, figsize=(36, 12))

    # Left: full overlay
    overlay = np.zeros((*dapi_norm.shape, 3), dtype=np.float32)
    overlay[:, :, 0] = hem_norm  # magenta = R
    overlay[:, :, 1] = dapi_norm  # green = G
    overlay[:, :, 2] = hem_norm  # magenta = B
    axes[0].imshow(np.clip(overlay, 0, 1))
    axes[0].set_title("Full ROI: DAPI (green) + hematoxylin (magenta)", fontweight="bold")
    axes[0].axis("off")

    # Center + right: zoomed crops to check cell-level alignment
    crop_size = 500  # pixels
    cy, cx = dapi_norm.shape[0] // 3, dapi_norm.shape[1] // 3
    for i, (y0, x0, label) in enumerate([
        (cy, cx, "Crop A (tumor nest)"),
        (cy * 2, cx * 2, "Crop B (stroma)"),
    ]):
        crop = overlay[y0:y0 + crop_size, x0:x0 + crop_size]
        axes[i + 1].imshow(np.clip(crop, 0, 1))
        axes[i + 1].set_title(label, fontweight="bold")
        axes[i + 1].axis("off")

    fig.suptitle("H&E–DAPI registration verification", fontweight="bold", fontsize=22)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIGURES / "registration_overlay.png", dpi=200,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved registration_overlay.png")

    # --- Stats ---
    from scipy.stats import pearsonr
    tissue = he.mean(axis=2) < 220
    n_sample = min(100_000, tissue.sum())
    idx = np.where(tissue.ravel())[0]
    rng = np.random.default_rng(42)
    sample = rng.choice(idx, size=n_sample, replace=False)
    r, _ = pearsonr(hem_norm.ravel()[sample], dapi_norm.ravel()[sample])
    print(f"\n  Hematoxylin vs DAPI correlation (tissue, n={n_sample}): r={r:.3f}")
    print(f"  Eosin range: [{eos_raw.min():.4f}, {eos_raw.max():.4f}]")
    print(f"  Hematoxylin range: [{hem_raw.min():.4f}, {hem_raw.max():.4f}]")


if __name__ == "__main__":
    main()
