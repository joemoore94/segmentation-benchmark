"""Merge tiled Baysor output into a per-cell AnnData for a nuclear prior variant.

Usage::

    conda run -n segbench python scripts/build_baysor_nuclear_prior_adata.py mesmer
    conda run -n segbench python scripts/build_baysor_nuclear_prior_adata.py 10x_ranger
    conda run -n segbench python scripts/build_baysor_nuclear_prior_adata.py stardist
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from segbench.quantify import quantify_baysor

ROI_DIR = Path("data/processed/roi")
GENE_RE = r"^b'(.*)'$"


def main(method: str) -> None:
    tiles_dir = ROI_DIR / f"baysor_{method}_prior_c10_tiles"

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

    out_path = ROI_DIR / f"adata_baysor_{method}_prior_c10.h5ad"
    adata.write_h5ad(out_path)
    print(f"wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{int(adata.X.sum())} assigned transcripts")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python build_baysor_nuclear_prior_adata.py <method>")
        sys.exit(1)
    main(sys.argv[1])
