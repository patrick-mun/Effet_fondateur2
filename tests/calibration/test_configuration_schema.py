"""Tests du contrat de configuration des panels de calibration."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = (
    REPOSITORY_ROOT / "calibration" / "schemas" / "calibration_config.schema.json"
)


def _schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def _validator() -> jsonschema.Draft202012Validator:
    schema = _schema()
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(schema)


def _micro_config() -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "config_id": "micro_panel_v1",
        "phase": "micro_panel",
        "panel": {
            "panel_id": "micro_technical_001",
            "role": "TECHNICAL",
            "individual_count": 24,
            "genome_scope": "REDUCED_AUTOSOMAL",
            "chromosomes": [19],
        },
        "random_seeds": {
            "master": 1001,
            "haplotypes": 1002,
            "pedigrees": 1003,
            "recombination": 1004,
            "measurement_error": 1005,
            "blinding": 1006,
        },
        "target": {
            "mode": "SIMULATED_ONLY",
            "variant_id": "synthetic_target_001",
            "ref": "A",
            "alt": "G",
        },
        "references": {
            "assembly": "GRCh38",
            "haplotype_source": {
                "source_id": "synthetic_haplotypes_v1",
                "mode": "SYNTHETIC_FIXTURE",
                "fixture_id": "tiny_haplotype_fixture_v1",
            },
            "genetic_map_source": {
                "source_id": "public_map_v1",
                "mode": "PUBLIC_CATALOG",
                "catalog_path": "config/references/genetic_maps.yaml",
                "public_reference_id": "public_map_release_v1",
            },
        },
        "scenarios": {
            "founder": {
                "ages_generations": [5, 20],
                "independent_family_counts": [3, 5],
            },
            "ibd": {
                "lengths_cm": [1.0, 2.0],
                "include_negative_controls": True,
            },
            "measurement_error": {
                "missingness_rates": [0.0, 0.01],
                "genotyping_error_rates": [0.0, 0.001],
                "phase_error_rates": [0.0, 0.01],
                "marker_density_fractions": [0.5, 1.0],
            },
            "population": {"stratum_count": 2, "include_admixed": False},
            "pedigree": {
                "trio_count": 1,
                "kinship_classes": ["UNRELATED", "FIRST_DEGREE"],
            },
            "autosomal_positions": [{"chromosome": 19, "position_bp": 100000}],
            "replications": {
                "kinship_pairs_per_class": 1,
                "ibd_negative_pairs": 1,
                "ibd_positive_per_length": 1,
                "dating_per_age": 1,
            },
        },
        "safety": {
            "allow_real_data": False,
            "allow_dock6_inputs": False,
            "expose_private_truth_to_pipeline": False,
            "allow_identifier_reuse": False,
            "allow_run_overwrite": False,
        },
        "reporting": {
            "language": "fr",
            "primary_format": "HTML",
            "include_json_summary": True,
            "include_sha256": True,
        },
    }


def _errors(config: dict[str, Any]) -> list[jsonschema.ValidationError]:
    return sorted(_validator().iter_errors(config), key=lambda error: list(error.path))


def test_schema_and_micro_panel_are_valid() -> None:
    assert _errors(_micro_config()) == []


@pytest.mark.parametrize(
    "missing_seed",
    [
        "master",
        "haplotypes",
        "pedigrees",
        "recombination",
        "measurement_error",
        "blinding",
    ],
)
def test_every_random_seed_is_required(missing_seed: str) -> None:
    config = _micro_config()
    del config["random_seeds"][missing_seed]

    assert _errors(config)


@pytest.mark.parametrize(
    ("phase", "role", "individual_count"),
    [
        ("micro_panel", "TECHNICAL", 23),
        ("micro_panel", "TECHNICAL", 61),
        ("pilot_calibration", "CALIBRATION", 499),
        ("pilot_validation", "VALIDATION", 501),
        ("qualification_calibration", "CALIBRATION", 2499),
        ("qualification_validation", "VALIDATION", 2501),
        ("pilot_calibration", "VALIDATION", 500),
        ("pilot_validation", "CALIBRATION", 500),
    ],
)
def test_inconsistent_phase_role_or_individual_count_is_rejected(
    phase: str,
    role: str,
    individual_count: int,
) -> None:
    config = _config_for_phase(phase, role, individual_count)

    assert _errors(config)


@pytest.mark.parametrize(
    ("phase", "role", "individual_count"),
    [
        ("micro_panel", "TECHNICAL", 24),
        ("pilot_calibration", "CALIBRATION", 500),
        ("pilot_validation", "VALIDATION", 500),
        ("qualification_calibration", "CALIBRATION", 2500),
        ("qualification_validation", "VALIDATION", 2500),
    ],
)
def test_each_declared_phase_has_a_consistent_profile(
    phase: str,
    role: str,
    individual_count: int,
) -> None:
    config = _config_for_phase(phase, role, individual_count)

    assert _errors(config) == []


@pytest.mark.parametrize(
    ("section", "field", "invalid_value"),
    [
        ("population", "stratum_count", 4),
        ("population", "include_admixed", False),
        ("pedigree", "trio_count", 99),
        ("replications", "kinship_pairs_per_class", 199),
        ("replications", "ibd_negative_pairs", 2999),
        ("replications", "ibd_positive_per_length", 399),
        ("replications", "dating_per_age", 499),
    ],
)
def test_qualification_minima_are_enforced(
    section: str,
    field: str,
    invalid_value: int | bool,
) -> None:
    config = _config_for_phase(
        "qualification_calibration",
        "CALIBRATION",
        2500,
    )
    config["scenarios"][section][field] = invalid_value

    assert _errors(config)


def test_scientific_panels_require_all_autosomes() -> None:
    config = _config_for_phase("pilot_calibration", "CALIBRATION", 500)
    config["panel"]["chromosomes"] = list(range(1, 22))

    assert _errors(config)


@pytest.mark.parametrize(
    "unsafe_path",
    [
        "/tmp/catalog.json",
        "~/catalog.json",
        "C:\\private\\catalog.json",
        "\\\\server\\share\\catalog.json",
        "file:///tmp/catalog.json",
        "../catalog.json",
        "config/../catalog.json",
        "data/input/dock6_real/catalog.json",
        "DATA/INPUT/DOCK6_REAL/catalog.json",
        "calibration_workspace/run/catalog.json",
        "Calibration_Workspace/run/catalog.json",
        "calibration/runs/catalog.json",
        "calibration/local/catalog.json",
    ],
)
def test_forbidden_catalog_paths_are_rejected(unsafe_path: str) -> None:
    config = _micro_config()
    config["references"]["genetic_map_source"]["catalog_path"] = unsafe_path

    assert _errors(config)


@pytest.mark.parametrize(
    "unsafe_flag",
    [
        "allow_real_data",
        "allow_dock6_inputs",
        "expose_private_truth_to_pipeline",
        "allow_identifier_reuse",
        "allow_run_overwrite",
    ],
)
def test_safety_flags_cannot_be_enabled(unsafe_flag: str) -> None:
    config = _micro_config()
    config["safety"][unsafe_flag] = True

    assert _errors(config)


def test_unknown_configuration_key_is_rejected() -> None:
    config = _micro_config()
    config["unexpected"] = True

    assert _errors(config)


def _config_for_phase(
    phase: str,
    role: str,
    individual_count: int,
) -> dict[str, Any]:
    config = deepcopy(_micro_config())
    config["phase"] = phase
    config["panel"]["role"] = role
    config["panel"]["individual_count"] = individual_count
    if phase != "micro_panel":
        config["panel"]["genome_scope"] = "AUTOSOMAL_GENOMEWIDE"
        config["panel"]["chromosomes"] = list(range(1, 23))
    if phase.startswith("qualification_"):
        config["scenarios"]["population"] = {
            "stratum_count": 5,
            "include_admixed": True,
        }
        config["scenarios"]["pedigree"]["trio_count"] = 100
        config["scenarios"]["replications"] = {
            "kinship_pairs_per_class": 200,
            "ibd_negative_pairs": 3000,
            "ibd_positive_per_length": 400,
            "dating_per_age": 500,
        }
    return config
