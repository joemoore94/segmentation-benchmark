"""Merge the 4 Baysor tile outputs into one full-ROI ``adata_baysor.h5ad``.

Each tile in ``data/processed/roi/baysor_tiles/<name>/segmentation.csv`` was
run on a padded ~1.05mm x 1.05mm region (see ``tile_baysor_transcripts.py``).
For each tile, this keeps only cells whose centroid (mean molecule position)
falls within that tile's non-padded "core" region -- the 4 cores exactly
partition the full 2mm x 2mm ROI, so this gives one cell per physical
location with no double-counting, while each kept cell still has its full set
of assigned molecules (including any in the tile's padding).

Also applies the ``b'GENENAME'`` -> ``GENENAME`` gene-name cleanup (the
``feature_name``/``gene`` columns are written from a bytes-typed pandas
column, see docs/dataset.md) before aggregating with
:func:`segbench.quantify.quantify_baysor`.

Usage::

    conda run -n segbench python scripts/build_baysor_adata.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from segbench.quantify import quantify_baysor

ROI_DIR = Path("data/processed/roi")
TILES_DIR = ROI_DIR / "baysor_tiles"

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
    # A transcript in the overlap band between two tiles can be assigned to a
    # core cell in both tiles' outputs; keep only the first occurrence so it's
    # not double-counted in the merged totals.
    before = len(merged)
    merged = merged.drop_duplicates(subset="transcript_id", keep="first")
    print(f"dropped {before - len(merged)} duplicate transcripts from overlap bands")

    adata = quantify_baysor(merged)

    out_path = ROI_DIR / "adata_baysor.h5ad"
    adata.write_h5ad(out_path)
    print(f"wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{int(adata.X.sum())} assigned transcripts")


if __name__ == "__main__":
    main()
