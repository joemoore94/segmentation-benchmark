"""Mesmer (DeepCell) segmentation wrapper, run via Docker.

deepcell-tf has no native Apple Silicon build, so Mesmer runs inside the
``vanvalenlab/deepcell-applications`` container instead of in-process.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from segbench.io import PIXEL_SIZE

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "run_mesmer.sh"


def run_mesmer(
    data_dir: Path,
    nuclear_file: str,
    output_dir: str,
    compartment: str = "nuclear",
    image_mpp: float = PIXEL_SIZE,
    membrane_file: str | None = None,
) -> Path:
    """Segment cells with Mesmer via the deepcell-applications Docker image.

    ``nuclear_file``, ``membrane_file``, and ``output_dir`` are paths relative
    to ``data_dir``, which is bind-mounted into the container. ``compartment``
    defaults to ``"nuclear"`` since our ``morphology_focus`` image is
    DAPI-only (no membrane channel).

    Returns ``data_dir / output_dir``, where the container writes the
    predicted label mask.
    """
    args = [
        str(_SCRIPT),
        str(data_dir),
        nuclear_file,
        output_dir,
        compartment,
        str(image_mpp),
    ]
    if membrane_file is not None:
        args.append(membrane_file)

    subprocess.run(args, check=True)
    return data_dir / output_dir
