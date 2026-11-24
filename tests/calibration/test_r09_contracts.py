"""Fermeture statique des contrats de la release unifiée."""

from pathlib import Path

import pytest

from effet_fondateur.calibration.contracts import (
    SCHEMA_FILES,
    STAGE_IDS,
    TRUTH_CONSUMERS,
    VIEW_CONSUMERS,
    CalibrationContractError,
    load_yaml,
    r08_section_sha256,
    validate_contract_bundle,
    validate_schema_set,
)
from effet_fondateur.calibration.evaluation import validate_criteria_document


ROOT = Path(__file__).resolve().parents[2]
CALIBRATION = ROOT / "calibration"


def test_schema_set_is_closed_and_valid() -> None:
    hashes = validate_schema_set(CALIBRATION / "schemas")
    assert tuple(hashes) == SCHEMA_FILES
    assert all(len(digest) == 64 for digest in hashes.values())


def test_bundle_reproduces_r03_to_r08_closure() -> None:
    bundle = load_yaml(CALIBRATION / "configurations/contract_bundle_v1.yaml")
    validate_contract_bundle(bundle, CALIBRATION / "PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md")
    assert len(bundle["governance_roots"]) == 8
    assert tuple(bundle["stage_ids"]) == STAGE_IDS
    assert set(bundle["truth_consumers"]) == set(TRUTH_CONSUMERS)
    assert set(bundle["view_consumers"]) == set(VIEW_CONSUMERS)
    assert bundle["dynamic_outputs"]["discovery_by_glob"] == "forbidden"


def test_r08_hash_is_bound_to_exact_plan_section() -> None:
    bundle = load_yaml(CALIBRATION / "configurations/contract_bundle_v1.yaml")
    assert bundle["r08_prespecification_sha256"] == r08_section_sha256(
        CALIBRATION / "PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md"
    )


def test_all_stage_criteria_are_frozen() -> None:
    criteria = load_yaml(CALIBRATION / "configurations/stage_criteria_v1.yaml")
    validate_criteria_document(criteria)
    assert tuple(criteria["stages"]) == STAGE_IDS
    assert "p_MC <= 0,05" in criteria["stages"]["16B"]["pass_contract"]
    assert "max(0,5 cM" in criteria["stages"]["13"]["pass_contract"]
    assert "HTML ou PDF" in criteria["stages"]["18"]["pass_contract"]


def test_bundle_rejects_a_truth_policy_drift() -> None:
    bundle = load_yaml(CALIBRATION / "configurations/contract_bundle_v1.yaml")
    bundle["truth_consumers"]["T-IBD"] = ["16B"]
    with pytest.raises(CalibrationContractError, match="truth_policy_mismatch"):
        validate_contract_bundle(bundle, CALIBRATION / "PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md")
