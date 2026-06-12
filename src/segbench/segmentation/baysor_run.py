"""Baysor transcript-based segmentation wrapper, run via the Julia CLI."""

from __future__ import annotations

import subprocess
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "run_baysor.sh"


def run_baysor(
    transcripts: Path,
    config: Path,
    output_dir: Path,
    prior_segmentation_column: str | None = None,
) -> Path:
    """Segment cells from transcript coordinates with Baysor.

    ``prior_segmentation_column``, if given, names a column in ``transcripts``
    holding a prior cell assignment (e.g. from 10x or a DAPI-based
    segmentation) used to seed Baysor.

    Returns ``output_dir``, containing Baysor's segmentation results.
    """
    args = [str(_SCRIPT), str(transcripts), str(config), str(output_dir)]
    if prior_segmentation_column is not None:
        args.append(prior_segmentation_column)

    subprocess.run(args, check=True)
    return output_dir
