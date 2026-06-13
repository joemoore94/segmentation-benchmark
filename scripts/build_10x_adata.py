"""Build an AnnData for 10x's native (Xenium Ranger) cell segmentation on the ROI.

10x's own segmentation is already computed -- nothing to run here. This just
merges its quantified cell x gene matrix (``cell_feature_matrix.h5ad``) with
per-cell spatial metadata (``cells.parquet``) into the same ``obs`` schema
(``centroid_x``/``centroid_y``/``n_transcripts``) used by ``adata_cellpose.h5ad``
and ``adata_baysor.h5ad``, so it can be compared via ``segbench.compare``.

Usage::

    conda run -n segbench python scripts/build_10x_adata.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

ROI_DIR = Path("data/processed/roi")


def main() -> None:
    adata = ad.read_h5ad(ROI_DIR / "cell_feature_matrix.h5ad")

    cells = pd.read_parquet(ROI_DIR / "cells.parquet")
    cells["cell_id"] = cells["cell_id"].apply(lambda b: b.decode())
    cells = cells.set_index("cell_id").reindex(adata.obs_names)

    adata.obs["centroid_x"] = cells["x_centroid"].to_numpy()
    adata.obs["centroid_y"] = cells["y_centroid"].to_numpy()
    adata.obs["n_transcripts"] = cells["transcript_counts"].to_numpy()
    # Already in um^2 (unlike CellPose's `area`, which is in px^2) -- kept
    # under distinct names so they don't get picked up by compare.size_summary's
    # `area`/`n_transcripts` branch and mixed with CellPose's pixel areas.
    adata.obs["nucleus_area_um2"] = cells["nucleus_area"].to_numpy()
    adata.obs["cell_area_um2"] = cells["cell_area"].to_numpy()

    out_path = ROI_DIR / "adata_10x.h5ad"
    adata.write_h5ad(out_path)
    print(f"wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
