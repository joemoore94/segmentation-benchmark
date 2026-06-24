"""Tile transcripts with a nuclear prior column for Baysor.

Generalises tile_baysor_prior_transcripts.py for any nuclear method.

Usage::

    conda run -n segbench python scripts/tile_nuclear_prior_transcripts.py stardist
    conda run -n segbench python scripts/tile_nuclear_prior_transcripts.py mesmer
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROI_DIR = Path("data/processed/roi")

ROI_EXTENT = (0.0, 2000.0)
MARGIN = 50.0

TILES = {
    "x0_y0": ((0.0, 1000.0), (0.0, 1000.0)),
    "x1_y0": ((1000.0, 2000.0), (0.0, 1000.0)),
    "x0_y1": ((0.0, 1000.0), (1000.0, 2000.0)),
    "x1_y1": ((1000.0, 2000.0), (1000.0, 2000.0)),
}


def _pad(bounds: tuple[float, float]) -> tuple[float, float]:
    lo, hi = bounds
    return (max(ROI_EXTENT[0], lo - MARGIN), min(ROI_EXTENT[1], hi + MARGIN))


def main(method: str) -> None:
    in_path = ROI_DIR / f"transcripts_baysor_{method}_prior.csv"
    out_dir = ROI_DIR / f"baysor_{method}_prior_c10_tiles"
    out_dir.mkdir(parents=True, exist_ok=True)

    transcripts = pd.read_csv(in_path)
    print(f"loaded {len(transcripts)} transcripts from {in_path.name}")

    manifest = {}
    for name, (x_core, y_core) in TILES.items():
        x_pad, y_pad = _pad(x_core), _pad(y_core)
        mask = (
            transcripts["x_location"].between(*x_pad)
            & transcripts["y_location"].between(*y_pad)
        )
        tile = transcripts.loc[mask]
        out_path = out_dir / f"transcripts_{name}.csv"
        tile.to_csv(out_path, index=False)
        manifest[name] = {
            "core_x": x_core,
            "core_y": y_core,
            "padded_x": x_pad,
            "padded_y": y_pad,
            "n_transcripts": len(tile),
        }
        print(f"{name}: {len(tile)} transcripts -> {out_path}")

    with open(out_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"manifest -> {out_dir / 'manifest.json'}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python tile_nuclear_prior_transcripts.py <stardist|mesmer>")
        sys.exit(1)
    main(sys.argv[1])
