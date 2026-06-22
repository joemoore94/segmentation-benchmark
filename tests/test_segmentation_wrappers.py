from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest


def test_run_mesmer_nuclear_only_invokes_script_with_expected_args(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.io import PIXEL_SIZE
    from segbench.segmentation import mesmer_run

    mock_run = MagicMock()
    monkeypatch.setattr(mesmer_run.subprocess, "run", mock_run)

    out = mesmer_run.run_mesmer(tmp_path, "dapi.tif", "mesmer_out")

    mock_run.assert_called_once()
    args = mock_run.call_args[0][0]
    assert args == ["bash", str(mesmer_run._SCRIPT),
                    str(tmp_path), "dapi.tif", "mesmer_out", "nuclear", str(PIXEL_SIZE)]
    assert out == tmp_path / "mesmer_out"


def test_run_mesmer_membrane_file_accepted_but_ignored(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.io import PIXEL_SIZE
    from segbench.segmentation import mesmer_run

    mock_run = MagicMock()
    monkeypatch.setattr(mesmer_run.subprocess, "run", mock_run)

    mesmer_run.run_mesmer(
        tmp_path, "dapi.tif", "mesmer_out", compartment="whole-cell", membrane_file="membrane.tif"
    )

    args = mock_run.call_args[0][0]
    assert args == ["bash", str(mesmer_run._SCRIPT),
                    str(tmp_path), "dapi.tif", "mesmer_out", "whole-cell", str(PIXEL_SIZE)]


def test_run_stardist_invokes_script_with_expected_args(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.segmentation import stardist_run

    mock_run = MagicMock()
    monkeypatch.setattr(stardist_run.subprocess, "run", mock_run)

    out = stardist_run.run_stardist(tmp_path, "dapi.tif", "stardist_out")

    mock_run.assert_called_once()
    args = mock_run.call_args[0][0]
    assert args[:5] == ["conda", "run", "-n", "stardist", "python"]
    assert args[5] == str(stardist_run._SCRIPT)
    assert args[6:] == [str(tmp_path), "dapi.tif", "stardist_out", "2D_versatile_fluo"]
    assert out == tmp_path / "stardist_out"


def test_run_stardist_with_custom_model_name(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from segbench.segmentation import stardist_run

    mock_run = MagicMock()
    monkeypatch.setattr(stardist_run.subprocess, "run", mock_run)

    stardist_run.run_stardist(tmp_path, "dapi.tif", "stardist_out", model_name="2D_paper_dsb2018")

    args = mock_run.call_args[0][0]
    assert args[6:] == [str(tmp_path), "dapi.tif", "stardist_out", "2D_paper_dsb2018"]


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
