"""Split the CellPose-prior transcript table into the same 4 tiles as ``tile_baysor_transcripts.py``.

Identical tiling scheme (2x2 grid of ~1mm x 1mm quadrants, each padded by
``MARGIN`` on interior edges) applied to ``transcripts_baysor_prior.csv``
(produced by ``add_cellpose_prior.py``), so the CellPose-prior Baysor run is
directly comparable to the unprimed run -- same tiles, same merge logic in
``build_baysor_prior_adata.py``.

Usage::

    conda run -n segbench python scripts/tile_baysor_prior_transcripts.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROI_DIR = Path("data/processed/roi")
OUT_DIR = ROI_DIR / "baysor_prior_tiles"

ROI_EXTENT = (0.0, 2000.0)  # microns, both x and y
TILE_SIZE = 1000.0  # microns
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

    transcripts = pd.read_csv(ROI_DIR / "transcripts_baysor_prior.csv")
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
