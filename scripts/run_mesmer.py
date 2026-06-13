"""Run Mesmer (DeepCell) segmentation on a nuclear (+ optional membrane) image.

DeepCell pins TensorFlow 2.8, which conflicts with the rest of the
``segbench`` stack, so this script runs in its own ``mesmer`` conda env::

    conda run -n mesmer python scripts/run_mesmer.py \\
        <data_dir> <nuclear_file> <output_dir> <compartment> <image_mpp> [membrane_file]

``nuclear_file``, ``membrane_file``, and ``output_dir`` are paths relative to
``data_dir``. ``compartment`` is one of ``nuclear``, ``whole-cell``, ``both``.
If ``membrane_file`` is omitted, the membrane channel is filled with zeros
(our ``morphology_focus`` image is DAPI-only).

Writes ``<data_dir>/<output_dir>/mask.tif`` (int32 label image).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tifffile
from deepcell.applications import Mesmer


def main() -> None:
    data_dir = Path(sys.argv[1])
    nuclear_file = sys.argv[2]
    output_dir = sys.argv[3]
    compartment = sys.argv[4]
    image_mpp = float(sys.argv[5])
    membrane_file = sys.argv[6] if len(sys.argv) > 6 else None

    nuclear = tifffile.imread(data_dir / nuclear_file).astype(np.float32)
    membrane = (
        tifffile.imread(data_dir / membrane_file).astype(np.float32)
        if membrane_file is not None
        else np.zeros_like(nuclear)
    )
    image = np.stack([nuclear, membrane], axis=-1)[np.newaxis, ...]  # (1, H, W, 2)

    app = Mesmer()
    labels = app.predict(image, image_mpp=image_mpp, compartment=compartment)

    out_dir = data_dir / output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    mask = labels[0, ..., 0].astype(np.int32)
    tifffile.imwrite(out_dir / "mask.tif", mask)
    print(f"Mesmer done, {mask.max()} cells -> {out_dir / 'mask.tif'}")


if __name__ == "__main__":
    main()
