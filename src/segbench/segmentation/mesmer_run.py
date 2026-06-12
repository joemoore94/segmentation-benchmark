"""Mesmer (DeepCell) segmentation wrapper, run via Docker.

deepcell-tf has no native Apple Silicon build, so Mesmer runs inside the
``vanvalenlab/deepcell-applications`` container instead of in-process.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "run_mesmer.sh"


def run_mesmer(
    data_dir: Path,
    nuclear_file: str,
    membrane_file: str,
    output_dir: str,
    compartment: str = "whole-cell",
) -> Path:
    """Segment cells with Mesmer via the deepcell-applications Docker image.

    ``nuclear_file``, ``membrane_file``, and ``output_dir`` are paths relative
    to ``data_dir``, which is bind-mounted into the container.

    Returns ``data_dir / output_dir``, where the container writes the
    predicted label mask.
    """
    subprocess.run(
        [
            str(_SCRIPT),
            str(data_dir),
            nuclear_file,
            membrane_file,
            output_dir,
            compartment,
        ],
        check=True,
    )
    return data_dir / output_dir
