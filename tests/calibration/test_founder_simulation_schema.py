"""Tests des contrats JSON D14."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema


ROOT = Path(__file__).resolve().parents[2]


def test_d14_configuration_is_valid_and_schemas_are_well_formed() -> None:
    config_schema = json.loads(
        (
            ROOT / "calibration/schemas/founder_simulation_config.schema.json"
        ).read_text()
    )
    manifest_schema = json.loads(
        (
            ROOT / "calibration/schemas/founder_simulation_manifest.schema.json"
        ).read_text()
    )
    configuration = json.loads(
        (
            ROOT
            / "calibration/configurations/"
            "d14_diagnostic_100kb_age150_calibration_v1.json"
        ).read_text()
    )
    jsonschema.Draft202012Validator.check_schema(config_schema)
    jsonschema.Draft202012Validator.check_schema(manifest_schema)
    jsonschema.Draft202012Validator(config_schema).validate(configuration)
