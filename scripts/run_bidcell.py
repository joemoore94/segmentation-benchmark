"""Run BIDCell segmentation on the extracted ROI.

BIDCell (Fu et al. 2024, Nature Communications) is a self-supervised deep
learning method that jointly uses DAPI morphology and transcript positions
with biologically-informed loss functions (including negative marker
co-expression penalties) to learn cell boundaries.

Prerequisites:
    pip install bidcell
    Download sc references to data/bidcell_ref/ (see configs/bidcell_config.yaml)

Usage::

    conda run -n segbench python scripts/run_bidcell.py

The output mask is saved to data/processed/roi/bidcell_out/ and converted
to an AnnData by build_bidcell_adata.py.
"""

from __future__ import annotations

import time
from pathlib import Path


def main() -> None:
    from bidcell import BIDCellModel

    config_path = "configs/bidcell_config.yaml"
    print(f"Loading BIDCell config: {config_path}")

    model = BIDCellModel(config_path)

    t0 = time.time()
    print("Running BIDCell pipeline (preprocess -> train -> predict)...")
    print("  This will take several hours on GPU.")
    model.run_pipeline()
    elapsed = time.time() - t0
    print(f"BIDCell pipeline complete in {elapsed / 60:.1f} minutes")

    out_dir = Path("data/processed/roi/bidcell_out")
    masks = list(out_dir.rglob("*_connected.tif"))
    if masks:
        print(f"Output mask: {masks[0]}")
    else:
        print("Warning: no *_connected.tif mask found in output directory")


if __name__ == "__main__":
    main()
