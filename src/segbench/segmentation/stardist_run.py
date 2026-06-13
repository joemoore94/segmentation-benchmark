"""StarDist segmentation wrapper, run in a separate conda env.

StarDist depends on TensorFlow, which conflicts with the PyTorch-based
``segbench`` stack, so StarDist runs in its own ``stardist`` conda env
(``conda create -n stardist python=3.10 && conda run -n stardist pip install
stardist tensorflow-cpu``), via ``conda run -n stardist python
scripts/run_stardist.py``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "run_stardist.py"
_CONDA_ENV = "stardist"


def run_stardist(
    data_dir: Path,
    nuclear_file: str,
    output_dir: str,
    model_name: str = "2D_versatile_fluo",
) -> Path:
    """Segment nuclei with StarDist, run via ``conda run -n stardist``.

    ``nuclear_file`` and ``output_dir`` are paths relative to ``data_dir``.
    ``model_name`` is a StarDist2D pretrained model name, default
    ``"2D_versatile_fluo"`` (suited to DAPI/fluorescence nuclear images).

    Returns ``data_dir / output_dir``, where ``scripts/run_stardist.py``
    writes the predicted label mask (``mask.tif``).
    """
    args = [
        "conda", "run", "-n", _CONDA_ENV, "python", str(_SCRIPT),
        str(data_dir), nuclear_file, output_dir, model_name,
    ]
    subprocess.run(args, check=True)
    return data_dir / output_dir
