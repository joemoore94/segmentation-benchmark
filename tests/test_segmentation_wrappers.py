from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest


def test_run_mesmer_invokes_script_with_expected_args(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.segmentation import mesmer_run

    mock_run = MagicMock()
    monkeypatch.setattr(mesmer_run.subprocess, "run", mock_run)

    out = mesmer_run.run_mesmer(tmp_path, "dapi.tif", "membrane.tif", "mesmer_out")

    mock_run.assert_called_once()
    args = mock_run.call_args[0][0]
    assert args[1:] == [str(tmp_path), "dapi.tif", "membrane.tif", "mesmer_out", "whole-cell"]
    assert out == tmp_path / "mesmer_out"


def test_run_baysor_with_prior_column(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.segmentation import baysor_run

    mock_run = MagicMock()
    monkeypatch.setattr(baysor_run.subprocess, "run", mock_run)

    transcripts = tmp_path / "transcripts.csv"
    config = tmp_path / "config.toml"
    out_dir = tmp_path / "baysor_out"

    result = baysor_run.run_baysor(
        transcripts, config, out_dir, prior_segmentation_column="cell_id"
    )

    mock_run.assert_called_once()
    args = mock_run.call_args[0][0]
    assert args[1:] == [str(transcripts), str(config), str(out_dir), "cell_id"]
    assert result == out_dir


def test_run_baysor_without_prior_column(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.segmentation import baysor_run

    mock_run = MagicMock()
    monkeypatch.setattr(baysor_run.subprocess, "run", mock_run)

    transcripts = tmp_path / "transcripts.csv"
    config = tmp_path / "config.toml"
    out_dir = tmp_path / "baysor_out"

    baysor_run.run_baysor(transcripts, config, out_dir)

    args = mock_run.call_args[0][0]
    assert args[1:] == [str(transcripts), str(config), str(out_dir)]
