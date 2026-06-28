"""Build expansion and Voronoi whole-cell AnnDatas from CellPose nuclei.

Three methods that bridge the nuclear-only vs. whole-cell gap:

- **Expansion 10µm**: each CellPose nucleus mask dilated outward by 10 µm
  (skimage expand_labels); overlapping expansions meet at the midpoint.
- **Expansion 20µm**: same with 20 µm radius, approaching full Voronoi coverage.
- **Voronoi**: every transcript assigned to the nearest CellPose nuclear
  centroid (scipy cKDTree), covering the entire ROI with no gaps.

Usage::

    conda run -n segbench python scripts/build_expanded_adatas.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import tifffile
from scipy.spatial import cKDTree
from skimage.segmentation import expand_labels

from segbench.constants import TOTAL_TRANSCRIPTS_FULL_ROI
from segbench.io import PIXEL_SIZE
from segbench.quantify import quantify_cells

ROI_DIR = Path("data/processed/roi")
EXPANSION_DISTANCES_UM = [10.0, 20.0]


def build_expansion(
    masks: np.ndarray,
    transcripts: pd.DataFrame,
    distance_um: float,
) -> ad.AnnData:
    radius_px = distance_um / PIXEL_SIZE
    expanded = expand_labels(masks, distance=radius_px)
    print(f"  expand_labels done: {(expanded > 0).sum()} foreground pixels "
          f"({(expanded > 0).mean():.1%} of ROI)")
    return quantify_cells(expanded, transcripts, pixel_size=PIXEL_SIZE)


def build_voronoi(
    adata_cellpose: ad.AnnData,
    transcripts: pd.DataFrame,
) -> ad.AnnData:
    centroids = adata_cellpose.obs[["centroid_x", "centroid_y"]].to_numpy()
    cell_ids = adata_cellpose.obs_names.to_numpy()

    xy_tx = transcripts[["x_location", "y_location"]].to_numpy()
    tree = cKDTree(centroids)
    _, idx = tree.query(xy_tx, workers=-1)

    tx = transcripts.copy()
    tx["cell_id"] = cell_ids[idx]

    counts = (
        tx.groupby(["cell_id", "feature_name"])
        .size()
        .unstack(fill_value=0)
    )

    obs = adata_cellpose.obs[["centroid_x", "centroid_y", "area"]].copy()
    obs = obs.reindex(counts.index)
    obs.index = obs.index.astype(str)
    counts.index = counts.index.astype(str)

    return ad.AnnData(
        X=counts.to_numpy(dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )


def main() -> None:
    transcripts = pd.read_parquet(ROI_DIR / "transcripts_baysor.parquet")
    print(f"loaded {len(transcripts)} transcripts")

    masks = tifffile.imread(ROI_DIR / "masks_cellpose.tif")
    print(f"cellpose mask shape: {masks.shape}, labels: {masks.max()}")

    adata_cellpose = ad.read_h5ad(ROI_DIR / "adata_cellpose.h5ad")

    for dist_um in EXPANSION_DISTANCES_UM:
        print(f"\n=== Expansion {dist_um:.0f}µm ===")
        adata = build_expansion(masks, transcripts, dist_um)
        out = ROI_DIR / f"adata_cellpose_exp{dist_um:.0f}um.h5ad"
        adata.write_h5ad(out)
        n_tx = int(adata.X.sum())
        capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI
        print(f"  {adata.n_obs} cells, {n_tx} transcripts ({capture:.1%} capture) -> {out.name}")

    print("\n=== Voronoi (CellPose centroids) ===")
    adata_vor = build_voronoi(adata_cellpose, transcripts)
    out = ROI_DIR / "adata_voronoi.h5ad"
    adata_vor.write_h5ad(out)
    n_tx = int(adata_vor.X.sum())
    capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI
    print(f"  {adata_vor.n_obs} cells, {n_tx} transcripts ({capture:.1%} capture) -> {out.name}")

    print("\n=== Voronoi (Mesmer centroids) ===")
    adata_mesmer = ad.read_h5ad(ROI_DIR / "adata_mesmer.h5ad")
    adata_vor_mesmer = build_voronoi(adata_mesmer, transcripts)
    out = ROI_DIR / "adata_voronoi_mesmer.h5ad"
    adata_vor_mesmer.write_h5ad(out)
    n_tx = int(adata_vor_mesmer.X.sum())
    capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI
    print(f"  {adata_vor_mesmer.n_obs} cells, {n_tx} transcripts ({capture:.1%} capture) -> {out.name}")

    print("\n=== Voronoi (StarDist centroids) ===")
    adata_stardist = ad.read_h5ad(ROI_DIR / "adata_stardist.h5ad")
    adata_vor_stardist = build_voronoi(adata_stardist, transcripts)
    out = ROI_DIR / "adata_voronoi_stardist.h5ad"
    adata_vor_stardist.write_h5ad(out)
    n_tx = int(adata_vor_stardist.X.sum())
    capture = n_tx / TOTAL_TRANSCRIPTS_FULL_ROI
    print(f"  {adata_vor_stardist.n_obs} cells, {n_tx} transcripts ({capture:.1%} capture) -> {out.name}")


if __name__ == "__main__":
    main()
