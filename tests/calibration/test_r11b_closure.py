from __future__ import annotations

from pathlib import Path

import pytest

import effet_fondateur.calibration.r11b_closure as r11b_closure
from effet_fondateur.calibration.contracts import load_json
from effet_fondateur.calibration.r11b_closure import (
    RECEIPT_RELATIVE,
    R11bClosureError,
    receipt_document,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r11b_tracked_closure_matches_completed_resource_evidence() -> None:
    receipt = receipt_document(REPOSITORY_ROOT)

    assert receipt["status"] == "RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE"
    assert receipt["shapeit_probe_count"] == 3
    assert receipt["ibd_probe_count"] == 18
    assert len(receipt["checks"]) == 5
    assert all(check["status"] == "PASS" for check in receipt["checks"])
    assert receipt["new_scientific_release_construction_authorized"] is True
    assert receipt["scientific_g5_passed"] is False
    assert receipt["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_r11b_tracked_receipt_is_reproducible() -> None:
    assert load_json(REPOSITORY_ROOT / RECEIPT_RELATIVE) == receipt_document(
        REPOSITORY_ROOT
    )


def test_r11b_closure_rejects_missing_decision(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        r11b_closure,
        "DECISION_RELATIVE",
        Path("calibration_workspace/resource_probes/missing/decision.json"),
    )

    with pytest.raises(R11bClosureError, match="r11b_decision_missing"):
        receipt_document(REPOSITORY_ROOT)
