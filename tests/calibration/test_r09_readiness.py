"""Reçu statique R09 et restrictions de la porte R10."""

from pathlib import Path

from effet_fondateur.calibration.contracts import load_json, load_yaml
from effet_fondateur.calibration.readiness import (
    assess_r09_repository,
    validate_readiness_receipt,
)


ROOT = Path(__file__).resolve().parents[2]
RECEIPT = ROOT / "calibration/readiness/R09_READY_TO_RUN.json"


def test_static_assessment_opens_only_r10_micro_construction() -> None:
    assessment = assess_r09_repository(ROOT)
    assert assessment.ready_to_run is True
    assert assessment.scope == "R10_MICRO_PANEL_CONSTRUCTION"
    assert len(assessment.terms) == 8
    assert len(assessment.evidence_sha256) == 15


def test_versioned_receipt_matches_all_current_evidence_hashes() -> None:
    assessment = validate_readiness_receipt(ROOT, RECEIPT)
    assert "NO_PANEL_GENERATED_IN_R09" in assessment.checks
    assert "NO_EXTERNAL_TOOL_EXECUTION_IN_R09" in assessment.checks


def test_receipt_keeps_all_qualification_flags_false() -> None:
    receipt = load_json(RECEIPT)
    assert receipt["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_micro_configuration_is_technical_and_not_generated() -> None:
    config = load_yaml(ROOT / "calibration/configurations/micro_panel_v1.yaml")
    assert config["generated"] is False
    assert config["scientific_performance_decision_allowed"] is False
    assert config["allowed_verdict_scope"] == "TECHNICAL_INTEGRITY_ONLY"
    assert config["autosomes"] == 22


def test_r09_tool_harness_cannot_execute_external_tools() -> None:
    harness = load_yaml(ROOT / "calibration/configurations/tool_harness_v1.yaml")
    assert harness["r09_external_execution_allowed"] is False
    assert all(
        tool["r09_mode"] == "SIMULATED_CONTRACT_ONLY"
        for tool in harness["external_tools"].values()
    )
