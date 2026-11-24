from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from effet_fondateur.calibration.r07_revision import (
    CONFIGURATION_PATH,
    R07RevisionError,
    assess_r07_revision,
    validate_tracked_receipt,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r07_v2_preserves_scientific_minima_and_closes_known_shapeit_budget() -> None:
    receipt = assess_r07_revision(REPOSITORY_ROOT)

    assert receipt["status"] == "REVISED_READY_FOR_R11B_RESOURCE_CONFIRMATION"
    assert receipt["scientific_minima_preserved"] is True
    assert receipt["shapeit_schedule"]["total_calls"] == 640
    assert receipt["shapeit_schedule"]["projected_core_hours_with_25_percent_reserve"] < 75_000
    assert receipt["ibd_schedule"] == {
        "independent_worlds": 90,
        "tool_density_calls": 540,
        "allocation_core_hours": 30_000,
        "resource_projection_complete": False,
    }
    assert receipt["scientific_g5_authorized"] is False
    assert receipt["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_r07_v2_tracked_receipt_matches_contract() -> None:
    assert validate_tracked_receipt(REPOSITORY_ROOT)["r07_revision_complete"] is True


def test_r07_v2_refuses_reduced_primary_effective_size(
    tmp_path: Path,
) -> None:
    configuration = yaml.safe_load((REPOSITORY_ROOT / CONFIGURATION_PATH).read_text())
    configuration["scientific_minima"]["primary_independent_units"] = 499
    copied_root = tmp_path / "repository"
    copied_configuration = copied_root / CONFIGURATION_PATH
    copied_projection = copied_root / (
        "calibration_workspace/resource_probes/"
        "r11_shapeit_sample_scaling_v1/projection.json"
    )
    copied_configuration.parent.mkdir(parents=True)
    copied_projection.parent.mkdir(parents=True)
    copied_configuration.write_text(yaml.safe_dump(configuration), encoding="utf-8")
    copied_projection.write_bytes(
        (
            REPOSITORY_ROOT
            / "calibration_workspace/resource_probes/"
            "r11_shapeit_sample_scaling_v1/projection.json"
        ).read_bytes()
    )

    with pytest.raises(R07RevisionError, match="r07_scientific_minimum_mismatch"):
        assess_r07_revision(copied_root)
