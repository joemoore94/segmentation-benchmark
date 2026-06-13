"""Run StarDist segmentation on a nuclear (DAPI) image.

StarDist depends on TensorFlow, which conflicts with the PyTorch-based
``segbench`` stack, so this script runs in its own ``stardist`` conda env::

    conda run -n stardist python scripts/run_stardist.py \\
        <data_dir> <nuclear_file> <output_dir> [model_name]

``nuclear_file`` and ``output_dir`` are paths relative to ``data_dir``.
``model_name`` is a StarDist2D pretrained model name, default
``2D_versatile_fluo`` (suited to DAPI/fluorescence nuclear images).

Writes ``<data_dir>/<output_dir>/mask.tif`` (int32 label image).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tifffile
from csbdeep.utils import normalize
from stardist.models import StarDist2D


def main() -> None:
    data_dir = Path(sys.argv[1])
    nuclear_file = sys.argv[2]
    output_dir = sys.argv[3]
    model_name = sys.argv[4] if len(sys.argv) > 4 else "2D_versatile_fluo"

    image = tifffile.imread(data_dir / nuclear_file)
    image = normalize(image, 1, 99.8)

    model = StarDist2D.from_pretrained(model_name)
    n_tiles = model._guess_n_tiles(image)
    labels, _details = model.predict_instances(image, n_tiles=n_tiles)

    out_dir = data_dir / output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    mask = labels.astype(np.int32)
    tifffile.imwrite(out_dir / "mask.tif", mask)
    print(f"StarDist done, {mask.max()} cells -> {out_dir / 'mask.tif'}")


if __name__ == "__main__":
    main()
