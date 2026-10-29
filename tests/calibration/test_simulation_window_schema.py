"""Contrats JSON D14b des fenêtres de simulation."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.simulation.windows import (
    SimulationWindowError,
    parse_window_configuration,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = (
    ROOT / "calibration/configurations/d14b_micro_chr5_age100_calibration_v1.json"
)
SCHEMA_PATH = ROOT / "calibration/schemas/simulation_window_config.schema.json"


def _document() -> dict[str, object]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_versioned_d14b_configuration_validates_and_parses() -> None:
    document = _document()
    jsonschema.Draft202012Validator(
        json.loads(SCHEMA_PATH.read_text(encoding="utf-8")),
        format_checker=jsonschema.FormatChecker(),
    ).validate(document)
    parsed = parse_window_configuration(document)
    assert parsed.panel_role == "CALIBRATION"
    assert {spec.purpose for spec in parsed.window_specs} == {
        "PRIMARY",
        "CENSORING_CHALLENGE",
    }
    assert len(parsed.ibd_requirements) == 5


def test_primary_cannot_override_formula_width() -> None:
    document = _document()
    document["window_specs"][0]["challenge_half_width_cm"] = 1.0
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(
            json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        ).validate(document)
    with pytest.raises(
        SimulationWindowError,
        match="simulation_window_primary_width_override_forbidden",
    ):
        parse_window_configuration(document)


def test_duplicate_ibd_length_is_refused_by_cross_validation() -> None:
    document = _document()
    duplicate = copy.deepcopy(document["marker_requirements"]["ibd_scenarios"][0])
    duplicate["scenario_id"] = "duplicate_length"
    document["marker_requirements"]["ibd_scenarios"].append(duplicate)
    with pytest.raises(
        SimulationWindowError, match="simulation_window_ibd_requirement_invalid"
    ):
        parse_window_configuration(document)


def test_configuration_rejects_probability_above_prespecified_limit() -> None:
    document = _document()
    document["method"]["maximum_arm_censor_probability"] = 0.1
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(
            json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        ).validate(document)
