"""Split the full-ROI transcript table into 4 overlapping quadrant tiles for Baysor.

Baysor's main-EM runtime scales roughly as N^1.8 with transcript count
(docs/dataset.md), so the full 2mm x 2mm ROI (3.39M transcripts) is run as 4
separate ~1mm x 1mm quadrants instead of one job. Each quadrant is padded by
``MARGIN`` on its interior edges so Baysor sees full local context near tile
boundaries; ``build_baysor_adata.py`` later keeps only cells whose centroid
falls in each tile's non-padded "core" region, so the padding is discarded
rather than double-counted.

Coordinates are left in the original ROI-global frame (no re-origin): Baysor's
``scale``/``scale_std`` are physical and origin-independent, and keeping global
coords means the merged output's centroids are directly usable for
``match_cells_by_centroid`` against the other (full-ROI) methods.

Usage::

    conda run -n segbench python scripts/tile_baysor_transcripts.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROI_DIR = Path("data/processed/roi")
OUT_DIR = ROI_DIR / "baysor_tiles"

ROI_EXTENT = (0.0, 2000.0)  # microns, both x and y
MARGIN = 50.0  # microns of padding on each interior tile edge

TILES = {
    "x0_y0": ((0.0, 1000.0), (0.0, 1000.0)),
    "x1_y0": ((1000.0, 2000.0), (0.0, 1000.0)),
    "x0_y1": ((0.0, 1000.0), (1000.0, 2000.0)),
    "x1_y1": ((1000.0, 2000.0), (1000.0, 2000.0)),
}


def _pad(bounds: tuple[float, float]) -> tuple[float, float]:
    lo, hi = bounds
    return (max(ROI_EXTENT[0], lo - MARGIN), min(ROI_EXTENT[1], hi + MARGIN))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    transcripts = pd.read_csv(ROI_DIR / "transcripts_baysor.csv")
    print(f"loaded {len(transcripts)} transcripts")

    manifest = {}
    for name, (x_core, y_core) in TILES.items():
        x_pad, y_pad = _pad(x_core), _pad(y_core)
        mask = (
            transcripts["x_location"].between(*x_pad)
            & transcripts["y_location"].between(*y_pad)
        )
        tile = transcripts.loc[mask]
        out_path = OUT_DIR / f"transcripts_{name}.csv"
        tile.to_csv(out_path, index=False)
        manifest[name] = {
            "core_x": x_core,
            "core_y": y_core,
            "padded_x": x_pad,
            "padded_y": y_pad,
            "n_transcripts": len(tile),
        }
        print(f"{name}: core x={x_core} y={y_core}, padded x={x_pad} y={y_pad}, "
              f"{len(tile)} transcripts -> {out_path}")

    with open(OUT_DIR / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)


if __name__ == "__main__":
    main()
