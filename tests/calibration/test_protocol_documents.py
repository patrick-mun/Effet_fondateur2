"""Contrôles structurels du protocole et des critères de calibration."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any

import pytest
import yaml
from yaml.constructor import ConstructorError
from yaml.resolver import BaseResolver


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PROTOCOL_PATH = (
    REPOSITORY_ROOT / "calibration" / "protocoles" / "PROTOCOLE_GENERAL.md"
)
CRITERIA_PATH = (
    REPOSITORY_ROOT / "calibration" / "protocoles" / "CRITERES_ACCEPTATION.yaml"
)
DECISION_REGISTER_PATH = (
    REPOSITORY_ROOT / "calibration" / "protocoles" / "REGISTRE_DECISIONS.md"
)


class UniqueKeyLoader(yaml.SafeLoader):
    """Charge un YAML sûr en refusant les clés répétées silencieusement."""


def _construct_unique_mapping(
    loader: UniqueKeyLoader,
    node: yaml.MappingNode,
    deep: bool = False,
) -> dict[Any, Any]:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise ConstructorError(
                "while constructing a mapping",
                node.start_mark,
                f"duplicate key: {key!r}",
                key_node.start_mark,
            )
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


UniqueKeyLoader.add_constructor(
    BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def _criteria_document() -> dict[str, Any]:
    document = yaml.load(
        CRITERIA_PATH.read_text(encoding="utf-8"),
        Loader=UniqueKeyLoader,
    )
    assert isinstance(document, dict)
    return document


def test_protocol_is_approved() -> None:
    protocol = PROTOCOL_PATH.read_text(encoding="utf-8")

    assert "**Statut :** `APPROVED`" in protocol
    assert "## 9. Aveuglement et déroulement d'un run" in protocol
    assert "## 11. Métriques et formules" in protocol
    assert "## 13. Cycle des critères d'acceptation" in protocol
    assert "## 20. Approbation" in protocol


def test_criteria_yaml_has_no_duplicate_keys() -> None:
    assert _criteria_document()["document"]["id"] == "CAL-ACCEPTANCE"


def test_every_criterion_has_a_complete_decision_contract() -> None:
    criteria = _criteria_document()["criteria"]
    required_fields = {
        "id",
        "reports",
        "metric",
        "unit",
        "direction",
        "acceptance",
        "missing",
    }
    allowed_missing_statuses = {
        "FAIL",
        "NOT_EVALUABLE",
        "REVIEW_REQUIRED",
    }

    assert len(criteria) == 70
    assert len({criterion["id"] for criterion in criteria}) == len(criteria)
    for criterion in criteria:
        assert required_fields <= criterion.keys()
        assert criterion["missing"]["status"] in allowed_missing_statuses
        assert criterion["missing"]["reason"]


def test_all_independent_stage_reports_are_covered() -> None:
    criteria = _criteria_document()["criteria"]
    covered_reports = {
        report_id
        for criterion in criteria
        for report_id in criterion["reports"]
    }
    expected_reports = {f"CAL-REP-{stage:02d}" for stage in range(20)} | {
        "CAL-REP-10A",
        "CAL-REP-12A",
        "CAL-REP-16A",
        "CAL-REP-16B",
        "CAL-REP-16C",
    }

    assert expected_reports <= covered_reports


def test_pilot_thresholds_are_explicitly_unsealed() -> None:
    criteria = _criteria_document()["criteria"]
    pilot_criteria = [
        criterion
        for criterion in criteria
        if criterion["acceptance"]["lifecycle"] == "pilot_then_seal"
    ]

    assert pilot_criteria
    for criterion in pilot_criteria:
        acceptance = criterion["acceptance"]
        threshold_fields = [
            field for field in ("value", "tolerance") if field in acceptance
        ]
        assert len(threshold_fields) == 1
        assert acceptance[threshold_fields[0]] is None


def test_fixed_contracts_have_a_non_null_value() -> None:
    criteria = _criteria_document()["criteria"]
    fixed_criteria = [
        criterion
        for criterion in criteria
        if criterion["acceptance"]["lifecycle"] == "fixed_contract"
    ]

    assert fixed_criteria
    for criterion in fixed_criteria:
        assert criterion["acceptance"]["value"] is not None


def test_d03_and_pre_pilot_d04_are_human_approved() -> None:
    document = _criteria_document()

    assert document["document"]["status"] == "APPROVED_PRE_PILOT"
    assert document["approval"]["d03_protocol_documented"] is True
    assert document["approval"]["d03_protocol_human_approved"] is True
    assert document["approval"]["human_approval"] is True
    assert document["approval"]["pilot_thresholds_sealed"] is False
    assert document["approval"]["validation_panel_opened"] is False


def test_approved_documents_match_registered_hashes() -> None:
    register = DECISION_REGISTER_PATH.read_text(encoding="utf-8")
    protocol_hash = sha256(PROTOCOL_PATH.read_bytes()).hexdigest()
    criteria_hash = sha256(CRITERIA_PATH.read_bytes()).hexdigest()

    assert "CAL-DEC-002 — Approbation de D03 et D04 avant pilote" in register
    assert protocol_hash in register
    assert criteria_hash in register


@pytest.mark.parametrize(
    "definition",
    [
        "sensitivity",
        "specificity",
        "false_positive_rate",
        "signed_bias",
        "mean_absolute_error",
        "rmse",
        "interval_coverage",
    ],
)
def test_core_metric_definitions_are_explicit(definition: str) -> None:
    assert _criteria_document()["metric_definitions"][definition]
