"""Mesmer (DeepCell) segmentation wrapper, run in a separate conda env.

DeepCell pins TensorFlow 2.8, which conflicts with the rest of the
``segbench`` stack, so Mesmer runs in its own ``mesmer`` conda env
(``conda create -n mesmer python=3.10 && conda run -n mesmer pip install
deepcell``), via ``conda run -n mesmer python scripts/run_mesmer.py``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from segbench.io import PIXEL_SIZE

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "run_mesmer.py"
_CONDA_ENV = "mesmer"


def run_mesmer(
    data_dir: Path,
    nuclear_file: str,
    output_dir: str,
    compartment: str = "nuclear",
    image_mpp: float = PIXEL_SIZE,
    membrane_file: str | None = None,
) -> Path:
    """Segment cells with Mesmer, run via ``conda run -n mesmer``.

    ``nuclear_file``, ``membrane_file``, and ``output_dir`` are paths relative
    to ``data_dir``. ``compartment`` defaults to ``"nuclear"`` since our
    ``morphology_focus`` image is DAPI-only (no membrane channel).

    Returns ``data_dir / output_dir``, where ``scripts/run_mesmer.py`` writes
    the predicted label mask (``mask.tif``).
    """
    args = [
        "conda", "run", "-n", _CONDA_ENV, "python", str(_SCRIPT),
        str(data_dir), nuclear_file, output_dir, compartment, str(image_mpp),
    ]
    if membrane_file is not None:
        args.append(membrane_file)

    subprocess.run(args, check=True)
    return data_dir / output_dir
