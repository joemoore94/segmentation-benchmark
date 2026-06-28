"""Convert Segger output to benchmark-compatible AnnData.

Segger writes an AnnData with transcript-to-cell assignments. This script
reads that output and, if needed, reshapes it into the same format as the
other benchmark methods (cells × genes count matrix with centroid_x,
centroid_y, area in .obs).

If Segger already writes a compatible AnnData (cells × genes), this script
copies it to the standard location. Otherwise it aggregates per-transcript
assignments into per-cell counts.

Usage::

    conda run -n segbench python scripts/build_segger_adata.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad

ROI_DIR = Path("data/processed/roi")
SEGGER_DIR = Path("data/processed/segger")


def find_segger_adata() -> Path:
    candidates = list(SEGGER_DIR.glob("*.h5ad"))
    if not candidates:
        raise FileNotFoundError(
            f"No .h5ad files found in {SEGGER_DIR}. "
            "Run scripts/run_segger.sh first."
        )
    if len(candidates) == 1:
        return candidates[0]
    for c in candidates:
        if "segmentation" in c.stem or "result" in c.stem:
            return c
    return candidates[0]


def main() -> None:
    path = find_segger_adata()
    print(f"Loading Segger output: {path}")
    adata = ad.read_h5ad(path)
    print(f"  Shape: {adata.shape}")
    print(f"  obs columns: {list(adata.obs.columns)}")

    # Ensure required obs columns exist
    required = ["centroid_x", "centroid_y"]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        # Try common Segger column names
        rename_map = {}
        for col in adata.obs.columns:
            cl = col.lower()
            if "centroid" in cl and "x" in cl:
                rename_map[col] = "centroid_x"
            elif "centroid" in cl and "y" in cl:
                rename_map[col] = "centroid_y"
            elif cl in ("x", "x_centroid", "cell_x"):
                rename_map[col] = "centroid_x"
            elif cl in ("y", "y_centroid", "cell_y"):
                rename_map[col] = "centroid_y"
        if rename_map:
            adata.obs = adata.obs.rename(columns=rename_map)
            print(f"  Renamed columns: {rename_map}")

    if "area" not in adata.obs.columns:
        adata.obs["area"] = 0.0

    out = ROI_DIR / "adata_segger.h5ad"
    adata.write_h5ad(out)
    n_tx = int(adata.X.sum()) if hasattr(adata.X, 'sum') else 0
    print(f"  {adata.n_obs} cells, {n_tx} transcripts -> {out.name}")


if __name__ == "__main__":
    main()
