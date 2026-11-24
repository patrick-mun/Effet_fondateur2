from __future__ import annotations

from pathlib import Path

from effet_fondateur.calibration.contracts import load_json
from effet_fondateur.calibration.r12c_scientific import (
    PREFLIGHT_PATH,
    preflight_document,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r12c_preflight_preserves_scientific_scope_without_opening_g5() -> None:
    receipt = preflight_document(REPOSITORY_ROOT)

    assert receipt["scope"] == "R13_SCIENTIFIC_BLIND_VALIDATION"
    assert receipt["partition_id"] == "VALIDATION_BLIND"
    assert receipt["status"] == "READY_FOR_G5_AUDIT"
    assert receipt["gate_g5_passed"] is False
    assert receipt["scientific_validation_authorized"] is False
    assert receipt["missing_before_g5"] == []
    assert receipt["metric_rule_bundle"]["bundle_state"] == "FULL_RULE_SET_FROZEN"
    assert receipt["metric_rule_bundle"]["stage_count"] == 25
    assert receipt["metric_rule_bundle"]["scientific_evaluation_authorized"] is False
    assert receipt["metric_producer_methods"]["rule_binding_count"] == 135
    assert receipt["metric_producer_methods"]["stage_producers_complete"] is True
    assert receipt["stage_metric_producers"]["implemented_stages"] == [
        "00",
        "01",
        "02",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
        "10A",
        "11",
        "12",
        "12A",
        "13",
        "14",
        "15",
        "16",
        "16A",
        "16B",
        "16C",
        "17",
        "18",
        "19",
    ]
    assert receipt["stage_metric_producers"]["missing_stage_count"] == 0
    assert receipt["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_r12c_tracked_preflight_is_reproducible() -> None:
    assert load_json(REPOSITORY_ROOT / PREFLIGHT_PATH) == preflight_document(
        REPOSITORY_ROOT
    )
