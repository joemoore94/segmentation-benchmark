"""Mesmer (DeepCell) segmentation wrapper, run via Docker.

The ``vanvalenlab/deepcell-applications`` image bundles pretrained model
weights, bypassing the DEEPCELL_ACCESS_TOKEN requirement added in the pip
package. Requires Docker to be installed and the current user in the
``docker`` group (``sudo usermod -aG docker $USER``).
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
    to ``data_dir``. ``compartment`` defaults to ``"nuclear"`` since our
    ``morphology_focus`` image is DAPI-only (no membrane channel).
    ``membrane_file`` is accepted for API compatibility but ignored (Mesmer
    fills the membrane channel with zeros when only a nuclear image is given).

    Returns ``data_dir / output_dir``, where the script writes ``mask.tif``.
    """
    subprocess.run(
        ["bash", str(_SCRIPT), str(data_dir), nuclear_file, output_dir,
         compartment, str(image_mpp)],
        check=True,
    )
    return data_dir / output_dir
