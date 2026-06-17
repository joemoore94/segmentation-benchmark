"""Merge tiled Baysor output into a per-cell AnnData for a given prior confidence.

Identical merge logic to build_baysor_prior_adata.py, parameterised for the
confidence sensitivity experiment (0.5 and 0.8). Tile output directories are
baysor_prior_c05_tiles/ and baysor_prior_c08_tiles/.

Usage::

    conda run -n segbench python scripts/build_baysor_conf_adata.py 0.5
    conda run -n segbench python scripts/build_baysor_conf_adata.py 0.8
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from segbench.quantify import quantify_baysor

ROI_DIR = Path("data/processed/roi")
GENE_RE = r"^b'(.*)'$"


def conf_label(confidence: float) -> str:
    return f"c{int(confidence * 10):02d}"


def main(confidence: float) -> None:
    label = conf_label(confidence)
    tiles_dir = ROI_DIR / f"baysor_prior_{label}_tiles"

    with open(tiles_dir / "manifest.json") as f:
        manifest = json.load(f)

    tables = []
    for name, info in manifest.items():
        seg = pd.read_csv(tiles_dir / name / "segmentation.csv")
        seg["gene"] = seg["gene"].str.replace(GENE_RE, r"\1", regex=True)

        assigned = seg.dropna(subset=["cell"])
        centroids = assigned.groupby("cell")[["x", "y"]].mean()
        core_x, core_y = info["core_x"], info["core_y"]
        core_cells = centroids.index[
            centroids["x"].between(*core_x) & centroids["y"].between(*core_y)
        ]

        tile_table = assigned[assigned["cell"].isin(core_cells)].copy()
        tile_table["cell"] = name + "_" + tile_table["cell"].astype(str)
        tables.append(tile_table)
        print(f"{name}: {len(core_cells)} core cells, "
              f"{len(tile_table)} assigned molecules")

    merged = pd.concat(tables, ignore_index=True)
    before = len(merged)
    merged = merged.drop_duplicates(subset="transcript_id", keep="first")
    print(f"dropped {before - len(merged)} duplicate transcripts from overlap bands")

    adata = quantify_baysor(merged)

    out_path = ROI_DIR / f"adata_baysor_prior_{label}.h5ad"
    adata.write_h5ad(out_path)
    print(f"wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{int(adata.X.sum())} assigned transcripts")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python build_baysor_conf_adata.py <confidence>")
        sys.exit(1)
    main(float(sys.argv[1]))
