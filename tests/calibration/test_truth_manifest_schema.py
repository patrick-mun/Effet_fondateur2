"""Tests du contrat de vérité privée des panels synthétiques."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = (
    REPOSITORY_ROOT / "calibration" / "schemas" / "truth_manifest.schema.json"
)


def _field(
    data_type: str,
    unit: str = "none",
    *,
    nullable: bool = False,
) -> dict[str, Any]:
    return {"data_type": data_type, "unit": unit, "nullable": nullable}


TRUTH_CONTRACTS: dict[str, tuple[dict[str, dict[str, Any]], list[str]]] = {
    "SAMPLE_REGISTRY": (
        {
            "sample_id": _field("string"),
            "panel_id": _field("string"),
            "family_id": _field("string", nullable=True),
            "is_synthetic": _field("boolean"),
        },
        ["CAL-REP-02"],
    ),
    "HAPLOTYPES_AND_GENOTYPES": (
        {
            "sample_id": _field("string"),
            "chromosome": _field("integer", "chromosome"),
            "position_bp": _field("integer", "bp"),
            "genotype": _field("categorical", "genotype"),
            "haplotype_1_allele": _field("string", "allele"),
            "haplotype_2_allele": _field("string", "allele"),
        },
        ["CAL-REP-03", "CAL-REP-12", "CAL-REP-12A"],
    ),
    "ANCESTRY": (
        {
            "sample_id": _field("string"),
            "stratum_id": _field("string"),
            "mixture_proportions": _field("object", "proportion"),
        },
        ["CAL-REP-08", "CAL-REP-16A"],
    ),
    "PEDIGREE_KINSHIP": (
        {
            "sample_1_id": _field("string"),
            "sample_2_id": _field("string"),
            "relationship_class": _field("categorical"),
            "kinship_coefficient": _field("number", "proportion"),
        },
        ["CAL-REP-07"],
    ),
    "TARGET_VARIANT": (
        {
            "sample_id": _field("string"),
            "variant_id": _field("string"),
            "genotype": _field("categorical", "genotype"),
            "mutant_haplotype": _field("categorical", "haplotype_index"),
        },
        ["CAL-REP-04", "CAL-REP-12"],
    ),
    "FOUNDER_TRANSMISSION": (
        {
            "replication_id": _field("string"),
            "family_id": _field("string"),
            "sample_id": _field("string"),
            "age_generations": _field("integer", "generations"),
            "left_boundary_bp": _field("integer", "bp"),
            "right_boundary_bp": _field("integer", "bp"),
            "left_length_cm": _field("number", "cM"),
            "right_length_cm": _field("number", "cM"),
            "total_length_cm": _field("number", "cM"),
            "recombination_events": _field("object"),
            "founder_haplotype_id": _field("string"),
            "carrier_haplotype_index": _field("integer", "haplotype_index"),
        },
        ["CAL-REP-13", "CAL-REP-14"],
    ),
    "IBD_SEGMENTS": (
        {
            "replication_id": _field("string"),
            "sample_1_id": _field("string"),
            "sample_2_id": _field("string"),
            "chromosome": _field("integer", "chromosome"),
            "start_bp": _field("integer", "bp", nullable=True),
            "end_bp": _field("integer", "bp", nullable=True),
            "start_cm": _field("number", "cM", nullable=True),
            "end_cm": _field("number", "cM", nullable=True),
            "scenario_id": _field("string"),
            "is_positive": _field("boolean"),
        },
        ["CAL-REP-16C"],
    ),
    "ROH_SEGMENTS": (
        {
            "replication_id": _field("string"),
            "sample_id": _field("string"),
            "chromosome": _field("integer", "chromosome"),
            "start_bp": _field("integer", "bp", nullable=True),
            "end_bp": _field("integer", "bp", nullable=True),
            "start_cm": _field("number", "cM", nullable=True),
            "end_cm": _field("number", "cM", nullable=True),
            "scenario_id": _field("string"),
            "is_positive": _field("boolean"),
        },
        ["CAL-REP-16"],
    ),
    "MEASUREMENT_ERRORS": (
        {
            "error_event_id": _field("string"),
            "error_profile_id": _field("string"),
            "sample_id": _field("string", nullable=True),
            "error_type": _field("categorical"),
            "affected_field": _field("string"),
            "original_value": _field("string", nullable=True),
            "injected_value": _field("string", nullable=True),
            "chromosome": _field("integer", "chromosome", nullable=True),
            "position_bp": _field("integer", "bp", nullable=True),
            "injected": _field("boolean"),
        },
        ["CAL-REP-05", "CAL-REP-10"],
    ),
}


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Refuse les clés JSON dupliquées, sinon masquées par ``json.loads``."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Clé JSON dupliquée dans le schéma D06 : {key}")
        result[key] = value
    return result


def _schema() -> dict[str, Any]:
    return json.loads(
        SCHEMA_PATH.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_keys,
    )


def _validator() -> jsonschema.Draft202012Validator:
    schema = _schema()
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(
        schema,
        format_checker=jsonschema.Draft202012Validator.FORMAT_CHECKER,
    )


def _artifact(truth_type: str, index: int) -> dict[str, Any]:
    fields, consumers = TRUTH_CONTRACTS[truth_type]
    slug = truth_type.lower()
    return {
        "artifact_id": f"truth_{slug}_v1",
        "truth_type": truth_type,
        "relative_path": f"private_truth/{slug}.tsv",
        "format": "TSV",
        "sha256": f"{index:x}" * 64,
        "record_count": 24,
        "sensitivity": "SYNTHETIC_GENETIC_PRIVATE",
        "blinded_export_allowed": False,
        "record_contract": {
            "schema_id": f"truth.{slug}.v1",
            "schema_version": "1.0.0",
            "schema_sha256": f"{index + 1:x}" * 64,
            "fields": deepcopy(fields),
        },
        "consumed_by_evaluators": list(consumers),
    }


def _valid_manifest() -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "truth_manifest_id": "truth_manifest_micro_001",
        "calibration_run_id": "cal_run_micro_001",
        "config_id": "micro_panel_v1",
        "panel": {
            "panel_id": "micro_technical_001",
            "role": "TECHNICAL",
            "individual_count": 24,
        },
        "created_at": "2026-10-28T12:00:00Z",
        "disclosure_state": "SEALED_PRIVATE",
        "coordinate_system": {
            "assembly": "GRCh38",
            "base_position": "BP_1_BASED_CLOSED",
            "genetic_distance": "cM",
            "age": "generations",
            "proportion": "UNIT_INTERVAL",
            "haplotype_index": "ONE_OR_TWO",
        },
        "provenance": {
            "generator_method_id": "synthetic_panel_generator_v1",
            "generator_version": "1.0.0",
            "code_commit_sha": "a" * 40,
            "configuration_sha256": "b" * 64,
            "seed_registry_sha256": "c" * 64,
            "public_reference_set_sha256": "d" * 64,
        },
        "target": {
            "mode": "SIMULATED_ONLY",
            "variant_id": "synthetic_target_001",
            "chromosome": 19,
            "position_bp": 100000,
            "ref": "A",
            "alt": "G",
        },
        "dimensions": {
            "individual_count": 24,
            "family_count": 3,
            "variant_count": 500,
            "haplotype_count": 48,
            "simulation_repetition_count": 1,
            "truth_artifact_count": 9,
        },
        "separation": {
            "real_data_source_count": 0,
            "real_identifier_mapping_count": 0,
            "calibration_validation_identifier_overlap_count": 0,
            "private_truth_reference_in_blinded_manifest_count": 0,
        },
        "truth_artifacts": [
            _artifact(truth_type, index)
            for index, truth_type in enumerate(TRUTH_CONTRACTS, start=1)
        ],
    }


def _errors(manifest: dict[str, Any]) -> list[jsonschema.ValidationError]:
    return sorted(
        _validator().iter_errors(manifest),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )


def _artifact_by_type(manifest: dict[str, Any], truth_type: str) -> dict[str, Any]:
    return next(
        artifact
        for artifact in manifest["truth_artifacts"]
        if artifact["truth_type"] == truth_type
    )


def test_schema_and_complete_truth_manifest_are_valid() -> None:
    assert _errors(_valid_manifest()) == []


@pytest.mark.parametrize("missing_truth_type", list(TRUTH_CONTRACTS))
def test_every_truth_category_is_required(missing_truth_type: str) -> None:
    manifest = _valid_manifest()
    artifacts = manifest["truth_artifacts"]
    missing_index = next(
        index
        for index, artifact in enumerate(artifacts)
        if artifact["truth_type"] == missing_truth_type
    )
    replacement = deepcopy(artifacts[(missing_index + 1) % len(artifacts)])
    replacement["artifact_id"] = f"duplicate_for_{missing_truth_type.lower()}"
    replacement["relative_path"] = (
        f"private_truth/duplicate_for_{missing_truth_type.lower()}.tsv"
    )
    artifacts[missing_index] = replacement

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("truth_type", "required_field"),
    [
        ("SAMPLE_REGISTRY", "sample_id"),
        ("HAPLOTYPES_AND_GENOTYPES", "haplotype_1_allele"),
        ("ANCESTRY", "mixture_proportions"),
        ("PEDIGREE_KINSHIP", "kinship_coefficient"),
        ("TARGET_VARIANT", "mutant_haplotype"),
        ("FOUNDER_TRANSMISSION", "age_generations"),
        ("IBD_SEGMENTS", "is_positive"),
        ("ROH_SEGMENTS", "is_positive"),
        ("MEASUREMENT_ERRORS", "error_type"),
    ],
)
def test_required_typed_fields_cannot_be_omitted(
    truth_type: str,
    required_field: str,
) -> None:
    manifest = _valid_manifest()
    artifact = _artifact_by_type(manifest, truth_type)
    del artifact["record_contract"]["fields"][required_field]

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("truth_type", "field_name", "invalid_definition"),
    [
        ("HAPLOTYPES_AND_GENOTYPES", "position_bp", _field("number", "cM")),
        ("FOUNDER_TRANSMISSION", "age_generations", _field("number", "cM")),
        ("IBD_SEGMENTS", "start_cm", _field("number", "cM")),
        ("ROH_SEGMENTS", "start_bp", _field("integer", "bp")),
        ("MEASUREMENT_ERRORS", "sample_id", _field("string", nullable=False)),
    ],
)
def test_critical_field_types_and_units_are_enforced(
    truth_type: str,
    field_name: str,
    invalid_definition: dict[str, Any],
) -> None:
    manifest = _valid_manifest()
    artifact = _artifact_by_type(manifest, truth_type)
    artifact["record_contract"]["fields"][field_name] = invalid_definition

    assert _errors(manifest)


@pytest.mark.parametrize(
    "unsafe_path",
    [
        "/private_truth/sample_registry.tsv",
        "../private_truth/sample_registry.tsv",
        "private_truth/../blinded_inputs/sample_registry.tsv",
        "blinded_inputs/sample_registry.tsv",
        "private_truth\\sample_registry.tsv",
        "file:///private_truth/sample_registry.tsv",
    ],
)
def test_truth_artifacts_must_use_safe_private_paths(unsafe_path: str) -> None:
    manifest = _valid_manifest()
    manifest["truth_artifacts"][0]["relative_path"] = unsafe_path

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("truth_type", "required_consumer"),
    [
        (truth_type, consumer)
        for truth_type, (_, consumers) in TRUTH_CONTRACTS.items()
        for consumer in consumers
    ],
)
def test_required_evaluator_consumers_are_declared(
    truth_type: str,
    required_consumer: str,
) -> None:
    manifest = _valid_manifest()
    artifact = _artifact_by_type(manifest, truth_type)
    artifact["consumed_by_evaluators"].remove(required_consumer)

    assert _errors(manifest)


@pytest.mark.parametrize(
    "separation_field",
    [
        "real_data_source_count",
        "real_identifier_mapping_count",
        "calibration_validation_identifier_overlap_count",
        "private_truth_reference_in_blinded_manifest_count",
    ],
)
def test_every_separation_violation_is_rejected(separation_field: str) -> None:
    manifest = _valid_manifest()
    manifest["separation"][separation_field] = 1

    assert _errors(manifest)


def test_truth_cannot_be_exported_to_blinded_inputs() -> None:
    manifest = _valid_manifest()
    manifest["truth_artifacts"][0]["blinded_export_allowed"] = True

    assert _errors(manifest)


def test_truth_manifest_must_be_sealed_and_private() -> None:
    manifest = _valid_manifest()
    manifest["disclosure_state"] = "OPENED_FOR_EVALUATION"

    assert _errors(manifest)


@pytest.mark.parametrize(
    "invalid_sha256",
    ["", "a" * 63, "A" * 64, "not-a-sha256"],
)
def test_artifact_checksums_are_strict(invalid_sha256: str) -> None:
    manifest = _valid_manifest()
    manifest["truth_artifacts"][0]["sha256"] = invalid_sha256

    assert _errors(manifest)


def test_unknown_manifest_key_is_rejected() -> None:
    manifest = _valid_manifest()
    manifest["unexpected"] = True

    assert _errors(manifest)
