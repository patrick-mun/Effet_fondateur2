from __future__ import annotations

from pathlib import Path

import pytest

from effet_fondateur.calibration.r11b_resources import (
    IBD_DENSITY_STEP,
    IBD_REPLICATES,
    IBD_TOOLS,
    PROBE_DIMENSIONS,
    R11bResourceError,
    _one_sided_prediction_bound,
    _shapeit_probe_root,
    _shapeit_fit_bound,
    _variant_positions,
    status,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r11b_design_covers_all_prescribed_resource_cells() -> None:
    assert PROBE_DIMENSIONS == {
        "dlow": (250, 250, "D_LOW"),
        "bridge": (2_000, 1_500, "D_MAX"),
        "target": (4_000, 3_000, "D_MAX"),
    }
    assert set(IBD_DENSITY_STEP) == {"D_MAX", "D_MID", "D_LOW"}
    assert IBD_TOOLS == ("HAP_IBD", "REFINED_IBD")
    assert len(IBD_DENSITY_STEP) * len(IBD_TOOLS) * IBD_REPLICATES == 18
    positions = _variant_positions()
    assert len(positions) == 26_000
    assert (positions[1:] > positions[:-1]).all()


def test_r11b_prediction_bounds_are_conservative() -> None:
    repeated_bound = _one_sided_prediction_bound([10.0, 11.0, 12.0])
    scaling_bound = _shapeit_fit_bound(
        [250, 500, 1_000, 3_500, 7_000],
        [300.0, 1_200.0, 3_900.0, 36_000.0, 130_000.0],
        7_000,
    )

    assert repeated_bound > 11.0
    assert scaling_bound > 130_000.0


def test_r11b_prediction_rejects_insufficient_replicates() -> None:
    with pytest.raises(R11bResourceError, match="r11b_prediction_values_invalid"):
        _one_sided_prediction_bound([10.0, 11.0])


def test_r11b_interrupted_bridge_uses_new_attempt_without_moving_other_probes(
    tmp_path: Path,
) -> None:
    assert _shapeit_probe_root(tmp_path, "bridge").name == "shapeit_bridge_attempt02"
    assert _shapeit_probe_root(tmp_path, "dlow").name == "shapeit_dlow_attempt01"
    assert _shapeit_probe_root(tmp_path, "target").name == "shapeit_target_attempt01"


def test_r11b_status_is_read_only_for_empty_workspace(tmp_path: Path) -> None:
    workspace = tmp_path / "r11b"

    result = status(REPOSITORY_ROOT, workspace)

    assert result["prepared"] is False
    assert result["ibd_completed"] == 0
    assert result["decision"] is None
    assert not workspace.exists()
