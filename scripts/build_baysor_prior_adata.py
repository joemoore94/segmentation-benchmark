"""Merge the 4 CellPose-prior Baysor tile outputs into ``adata_baysor_prior.h5ad``.

Identical merge logic to ``build_baysor_adata.py`` (core-region partitioning +
``transcript_id`` dedup across overlapping tile padding), applied to
``data/processed/roi/baysor_prior_tiles/<tile>/segmentation.csv`` -- the
CellPose-nucleus-prior hybrid run (see ``add_cellpose_prior.py``,
``tile_baysor_prior_transcripts.py``, ``configs/baysor_prior_config.toml``).

Usage::

    conda run -n segbench python scripts/build_baysor_prior_adata.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from segbench.quantify import quantify_baysor

ROI_DIR = Path("data/processed/roi")
TILES_DIR = ROI_DIR / "baysor_prior_tiles"

GENE_RE = r"^b'(.*)'$"


def main() -> None:
    with open(TILES_DIR / "manifest.json") as f:
        manifest = json.load(f)

    tables = []
    for name, info in manifest.items():
        seg = pd.read_csv(TILES_DIR / name / "segmentation.csv")
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

    out_path = ROI_DIR / "adata_baysor_prior.h5ad"
    adata.write_h5ad(out_path)
    print(f"wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{int(adata.X.sum())} assigned transcripts")


if __name__ == "__main__":
    main()
