from __future__ import annotations

import numpy as np
import pytest


def test_run_cellpose_default_channels_and_dtype(monkeypatch: pytest.MonkeyPatch) -> None:
    """run_cellpose should default to single-channel input and return int32 masks."""
    from segbench.segmentation import cellpose_run

    image = np.zeros((32, 32), dtype=np.float32)
    fake_masks = np.zeros((32, 32), dtype=np.uint16)
    fake_masks[10:20, 10:20] = 1

    captured: dict[str, object] = {}

    class FakeCellposeModel:
        def __init__(self, gpu: bool, pretrained_model: str) -> None:
            captured["gpu"] = gpu
            captured["pretrained_model"] = pretrained_model

        def eval(
            self, img: np.ndarray, diameter: float | None, channels: list[int]
        ) -> tuple[np.ndarray, None, None]:
            captured["channels"] = channels
            captured["diameter"] = diameter
            return fake_masks, None, None

    monkeypatch.setattr("cellpose.models.CellposeModel", FakeCellposeModel)

    result = cellpose_run.run_cellpose(image, model_type="cpsam_v2", gpu=False)

    assert result.dtype == np.int32
    assert captured["channels"] == [0, 0]
    assert captured["pretrained_model"] == "cpsam_v2"
    assert captured["gpu"] is False
