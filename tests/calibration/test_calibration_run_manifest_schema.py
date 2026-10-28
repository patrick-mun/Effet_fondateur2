"""Tests du contrat de reproductibilité des runs de calibration."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "schemas"
    / "calibration_run_manifest.schema.json"
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Refuse les clés JSON dupliquées dans le contrat D07."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Clé JSON dupliquée dans le schéma D07 : {key}")
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
        format_checker=jsonschema.FormatChecker(),
    )


def _sha(character: str) -> str:
    return character * 64


def _schema_contract(schema_id: str, hash_character: str) -> dict[str, Any]:
    return {
        "schema_id": schema_id,
        "schema_version": "1.0.0",
        "schema_sha256": _sha(hash_character),
    }


def _artifact(
    artifact_id: str,
    artifact_type: str,
    relative_path: str,
    hash_character: str,
    *,
    sensitivity: str,
    pipeline_access_allowed: bool,
    schema: dict[str, Any] | None = None,
    producer_command_id: str | None = "generate_panel",
) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "artifact_type": artifact_type,
        "relative_path": relative_path,
        "sha256": _sha(hash_character),
        "size_bytes": 128,
        "media_type": "application/json",
        "producer_command_id": producer_command_id,
        "schema": schema,
        "sensitivity": sensitivity,
        "immutable": True,
        "pipeline_access_allowed": pipeline_access_allowed,
    }


def _valid_manifest() -> dict[str, Any]:
    config_schema = _schema_contract("calibration_config.schema.json", "1")
    truth_schema = _schema_contract("truth_manifest.schema.json", "2")
    artifacts = [
        _artifact(
            "calibration_config",
            "CONFIGURATION",
            "control/config.resolved.json",
            "3",
            sensitivity="INTERNAL_CONTROL",
            pipeline_access_allowed=False,
            schema=config_schema,
            producer_command_id=None,
        ),
        _artifact(
            "truth_manifest",
            "TRUTH_MANIFEST",
            "private_truth/truth_manifest.json",
            "4",
            sensitivity="SYNTHETIC_GENETIC_PRIVATE",
            pipeline_access_allowed=False,
            schema=truth_schema,
        ),
        _artifact(
            "seed_registry",
            "PRIVATE_TRUTH",
            "private_truth/seed_registry.json",
            "8",
            sensitivity="SYNTHETIC_GENETIC_PRIVATE",
            pipeline_access_allowed=False,
        ),
        _artifact(
            "blinded_panel",
            "BLINDED_INPUT",
            "blinded_inputs/panel.vcf.gz",
            "5",
            sensitivity="SYNTHETIC_BLINDED",
            pipeline_access_allowed=True,
            producer_command_id="generate_panel",
        ),
        _artifact(
            "environment_lock",
            "ENVIRONMENT_LOCK",
            "control/environment.freeze.txt",
            "6",
            sensitivity="INTERNAL_CONTROL",
            pipeline_access_allowed=False,
            producer_command_id=None,
        ),
        _artifact(
            "checksum_index",
            "CHECKSUM_INDEX",
            "checksums.sha256",
            "7",
            sensitivity="INTERNAL_CONTROL",
            pipeline_access_allowed=False,
            producer_command_id="generate_panel",
        ),
    ]
    return {
        "schema_version": "1.0.0",
        "run_manifest_id": "run_manifest_micro_001",
        "calibration_run_id": "cal_run_micro_001",
        "config_id": "micro_panel_v1",
        "phase": "micro_panel",
        "panel": {
            "panel_id": "micro_technical_001",
            "role": "TECHNICAL",
            "individual_count": 24,
        },
        "manifest_state": "FINALIZED",
        "access_scope": "CALIBRATION_CONTROL_PRIVATE",
        "created_at": "2026-10-28T12:00:00Z",
        "completed_at": "2026-10-28T12:01:00Z",
        "run_status": "SUCCEEDED",
        "attempt": {
            "attempt_number": 1,
            "resumed_from_manifest_sha256": None,
            "inputs_unchanged": True,
            "parameters_unchanged": True,
            "code_unchanged": True,
        },
        "provenance": {
            "repository_id": "Effet_fondateur2",
            "code_commit_sha": "a" * 40,
            "dirty_worktree": False,
            "pipeline_version": "2.0.0",
            "python_implementation": "CPython",
            "python_version": "3.12.7",
            "platform_system": "Darwin",
            "platform_release": "25.0.0",
            "machine": "arm64",
            "timezone": "Indian/Reunion",
            "locale": "fr_FR.UTF-8",
            "dependency_lock_path": "control/environment.freeze.txt",
            "dependency_lock_sha256": _sha("6"),
        },
        "contracts": {
            "configuration": {
                "relative_path": "control/config.resolved.json",
                "sha256": _sha("3"),
                "schema": config_schema,
            },
            "truth_manifest": {
                "relative_path": "private_truth/truth_manifest.json",
                "sha256": _sha("4"),
                "schema": truth_schema,
                "pipeline_access_allowed": False,
            },
        },
        "random_seeds": {
            "registry_path": "private_truth/seed_registry.json",
            "registry_sha256": _sha("8"),
            "values": {
                "master": 1001,
                "haplotypes": 1002,
                "pedigrees": 1003,
                "recombination": 1004,
                "measurement_error": 1005,
                "blinding": 1006,
            },
        },
        "sources": [
            {
                "source_id": "tiny_haplotype_fixture_v1",
                "source_type": "SYNTHETIC_FIXTURE",
                "version": "1.0.0",
                "origin_id": "synthetic_fixture_catalog_v1",
                "license_id": "INTERNAL_SYNTHETIC",
                "assembly": "GRCh38",
                "local_artifact_id": "blinded_panel",
                "sha256": _sha("5"),
                "size_bytes": 128,
                "acquired_at": "2026-10-28T11:59:00Z",
                "shared_between_panels": False,
            }
        ],
        "tools": [
            {
                "tool_id": "python",
                "tool_type": "EXECUTABLE",
                "version": "3.12.7",
                "invocation": "python",
                "sha256": _sha("9"),
                "version_probe_argv": ["python", "--version"],
                "version_probe_output_sha256": _sha("a"),
            }
        ],
        "commands": [
            {
                "command_id": "generate_panel",
                "scope": "GENERATOR",
                "stage_id": "GEN",
                "tool_id": "python",
                "normalized_argv": [
                    ".venv/bin/python",
                    "-m",
                    "effet_fondateur.calibration.generate",
                    "--config",
                    "control/config.resolved.json",
                ],
                "working_directory": "RUN_ROOT",
                "environment": [
                    {"name": "TZ", "value": "Indian/Reunion"},
                    {"name": "PYTHONHASHSEED", "value": "1001"},
                ],
                "seed_names": [
                    "master",
                    "haplotypes",
                    "pedigrees",
                    "recombination",
                    "measurement_error",
                    "blinding",
                ],
                "input_artifact_ids": ["calibration_config"],
                "output_artifact_ids": [
                    "truth_manifest",
                    "seed_registry",
                    "blinded_panel",
                    "checksum_index",
                ],
                "resource_limits": {
                    "threads": 1,
                    "timeout_seconds": 300,
                    "memory_mb": 1024,
                },
                "started_at": "2026-10-28T12:00:00Z",
                "completed_at": "2026-10-28T12:01:00Z",
                "duration_seconds": 60.0,
                "exit_code": 0,
                "status": "SUCCEEDED",
                "stdout_artifact_id": None,
                "stderr_artifact_id": None,
            }
        ],
        "artifacts": artifacts,
        "replay": {
            "workspace_layout_version": "1.0.0",
            "start_from_clean_workspace": True,
            "command_order": ["generate_panel"],
            "verification_argv": [
                ".venv/bin/python",
                "-m",
                "effet_fondateur.calibration.verify_run",
                "--run-root",
                "{RUN_ROOT}",
            ],
            "expected_checksum_index_sha256": _sha("7"),
            "determinism_expected": True,
            "nondeterministic_components": [],
        },
        "integrity": {
            "checksum_index_artifact_id": "checksum_index",
            "artifact_count": 6,
            "all_sources_hashed": True,
            "all_tools_hashed": True,
            "all_artifacts_hashed": True,
            "commands_normalized": True,
            "seeds_complete": True,
            "real_data_source_count": 0,
            "private_truth_reference_in_pipeline_command_count": 0,
            "calibration_validation_identifier_overlap_count": 0,
            "run_overwrite_detected": False,
        },
        "deviations": [],
    }


def _errors(manifest: dict[str, Any]) -> list[jsonschema.ValidationError]:
    return sorted(
        _validator().iter_errors(manifest),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )


def _artifact_by_type(
    manifest: dict[str, Any],
    artifact_type: str,
) -> dict[str, Any]:
    return next(
        artifact
        for artifact in manifest["artifacts"]
        if artifact["artifact_type"] == artifact_type
    )


def test_schema_and_complete_run_manifest_are_valid() -> None:
    assert _errors(_valid_manifest()) == []


@pytest.mark.parametrize(
    "required_section",
    [
        "attempt",
        "provenance",
        "contracts",
        "random_seeds",
        "sources",
        "tools",
        "commands",
        "artifacts",
        "replay",
        "integrity",
        "deviations",
    ],
)
def test_every_reproducibility_section_is_required(required_section: str) -> None:
    manifest = _valid_manifest()
    del manifest[required_section]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "provenance_field",
    [
        "code_commit_sha",
        "pipeline_version",
        "python_version",
        "platform_system",
        "platform_release",
        "machine",
        "timezone",
        "locale",
        "dependency_lock_path",
        "dependency_lock_sha256",
    ],
)
def test_environment_provenance_is_complete(provenance_field: str) -> None:
    manifest = _valid_manifest()
    del manifest["provenance"][provenance_field]

    assert _errors(manifest)


def test_dirty_worktree_is_not_reproducible() -> None:
    manifest = _valid_manifest()
    manifest["provenance"]["dirty_worktree"] = True

    assert _errors(manifest)


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
    manifest = _valid_manifest()
    del manifest["random_seeds"]["values"][missing_seed]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "source_field",
    [
        "version",
        "origin_id",
        "license_id",
        "local_artifact_id",
        "sha256",
        "size_bytes",
        "acquired_at",
    ],
)
def test_sources_are_versioned_and_hashed(source_field: str) -> None:
    manifest = _valid_manifest()
    del manifest["sources"][0][source_field]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "forbidden_origin",
    [
        "data/input/dock6_real/panel.vcf.gz",
        "DATA/INPUT/reference.vcf.gz",
        "dock6_real_metadata",
    ],
)
def test_real_data_sources_are_rejected(forbidden_origin: str) -> None:
    manifest = _valid_manifest()
    manifest["sources"][0]["origin_id"] = forbidden_origin

    assert _errors(manifest)


def test_only_public_references_may_be_shared_between_panels() -> None:
    manifest = _valid_manifest()
    manifest["sources"][0]["shared_between_panels"] = True

    assert _errors(manifest)


@pytest.mark.parametrize(
    "tool_field",
    [
        "version",
        "invocation",
        "sha256",
        "version_probe_argv",
        "version_probe_output_sha256",
    ],
)
def test_tools_have_exact_versions_and_binary_provenance(tool_field: str) -> None:
    manifest = _valid_manifest()
    del manifest["tools"][0][tool_field]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "command_field",
    [
        "tool_id",
        "normalized_argv",
        "working_directory",
        "environment",
        "seed_names",
        "input_artifact_ids",
        "output_artifact_ids",
        "resource_limits",
        "started_at",
        "completed_at",
        "duration_seconds",
        "exit_code",
        "status",
    ],
)
def test_commands_contain_their_complete_execution_contract(
    command_field: str,
) -> None:
    manifest = _valid_manifest()
    del manifest["commands"][0][command_field]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "unsafe_argument",
    [
        "/tmp/panel.vcf.gz",
        "--input=/tmp/panel.vcf.gz",
        "~/panel.vcf.gz",
        "--input=~/panel.vcf.gz",
        "file:///tmp/panel.vcf.gz",
        "--input=file:///tmp/panel.vcf.gz",
        "C:\\temp\\panel.vcf.gz",
        "--input=C:\\temp\\panel.vcf.gz",
        "../panel.vcf.gz",
        "--input=../panel.vcf.gz",
        "data/input/reference.vcf.gz",
        "dock6_real/panel.vcf.gz",
    ],
)
def test_commands_reject_local_or_real_data_paths(unsafe_argument: str) -> None:
    manifest = _valid_manifest()
    manifest["commands"][0]["normalized_argv"].append(unsafe_argument)

    assert _errors(manifest)


@pytest.mark.parametrize(
    "truth_argument",
    [
        "private_truth/truth_manifest.json",
        "PRIVATE_TRUTH/segments.tsv",
        "--truth=private_truth/segments.tsv",
    ],
)
def test_pipeline_commands_cannot_reference_private_truth(
    truth_argument: str,
) -> None:
    manifest = _valid_manifest()
    command = manifest["commands"][0]
    command["scope"] = "PIPELINE"
    command["stage_id"] = "00"
    command["normalized_argv"].append(truth_argument)

    assert _errors(manifest)


@pytest.mark.parametrize(
    "secret_variable",
    ["TOKEN", "PASSWORD", "AWS_SECRET_ACCESS_KEY", "GITHUB_TOKEN"],
)
def test_command_environment_cannot_record_secret_variables(
    secret_variable: str,
) -> None:
    manifest = _valid_manifest()
    manifest["commands"][0]["environment"].append(
        {"name": secret_variable, "value": "synthetic-secret"}
    )

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("status", "exit_code"),
    [("SUCCEEDED", 1), ("FAILED", 0)],
)
def test_command_status_and_exit_code_must_agree(
    status: str,
    exit_code: int,
) -> None:
    manifest = _valid_manifest()
    manifest["commands"][0]["status"] = status
    manifest["commands"][0]["exit_code"] = exit_code

    assert _errors(manifest)


@pytest.mark.parametrize(
    "required_artifact_type",
    [
        "CONFIGURATION",
        "TRUTH_MANIFEST",
        "PRIVATE_TRUTH",
        "BLINDED_INPUT",
        "ENVIRONMENT_LOCK",
        "CHECKSUM_INDEX",
    ],
)
def test_minimum_reproduction_artifacts_are_required(
    required_artifact_type: str,
) -> None:
    manifest = _valid_manifest()
    artifact = _artifact_by_type(manifest, required_artifact_type)
    artifact["artifact_type"] = "LOG"

    assert _errors(manifest)


@pytest.mark.parametrize(
    "unsafe_path",
    [
        "/tmp/artifact.json",
        "../artifact.json",
        "control/../artifact.json",
        "file:///tmp/artifact.json",
        "data/input/artifact.json",
        "dock6_real/artifact.json",
        "calibration_workspace/run/artifact.json",
    ],
)
def test_artifacts_use_portable_relative_paths(unsafe_path: str) -> None:
    manifest = _valid_manifest()
    manifest["artifacts"][0]["relative_path"] = unsafe_path

    assert _errors(manifest)


def test_truth_artifacts_remain_private_and_inaccessible_to_pipeline() -> None:
    manifest = _valid_manifest()
    truth = _artifact_by_type(manifest, "TRUTH_MANIFEST")
    truth["relative_path"] = "blinded_inputs/truth_manifest.json"
    truth["pipeline_access_allowed"] = True

    assert _errors(manifest)


def test_blinded_inputs_must_be_explicitly_accessible_to_pipeline() -> None:
    manifest = _valid_manifest()
    blinded = _artifact_by_type(manifest, "BLINDED_INPUT")
    blinded["pipeline_access_allowed"] = False

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("attempt_field", "invalid_value"),
    [
        ("resumed_from_manifest_sha256", _sha("b")),
        ("inputs_unchanged", False),
        ("parameters_unchanged", False),
        ("code_unchanged", False),
    ],
)
def test_first_attempt_cannot_claim_a_changed_or_resumed_run(
    attempt_field: str,
    invalid_value: str | bool,
) -> None:
    manifest = _valid_manifest()
    manifest["attempt"][attempt_field] = invalid_value

    assert _errors(manifest)


def test_resume_requires_previous_manifest_hash_and_unchanged_contracts() -> None:
    manifest = _valid_manifest()
    manifest["attempt"]["attempt_number"] = 2
    manifest["attempt"]["resumed_from_manifest_sha256"] = None

    assert _errors(manifest)

    manifest["attempt"]["resumed_from_manifest_sha256"] = _sha("b")
    manifest["attempt"]["parameters_unchanged"] = False

    assert _errors(manifest)


@pytest.mark.parametrize(
    "replay_field",
    [
        "workspace_layout_version",
        "start_from_clean_workspace",
        "command_order",
        "verification_argv",
        "expected_checksum_index_sha256",
        "determinism_expected",
        "nondeterministic_components",
    ],
)
def test_replay_instructions_are_complete(replay_field: str) -> None:
    manifest = _valid_manifest()
    del manifest["replay"][replay_field]

    assert _errors(manifest)


def test_nondeterminism_must_be_named() -> None:
    manifest = _valid_manifest()
    manifest["replay"]["determinism_expected"] = False

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("integrity_field", "invalid_value"),
    [
        ("all_sources_hashed", False),
        ("all_tools_hashed", False),
        ("all_artifacts_hashed", False),
        ("commands_normalized", False),
        ("seeds_complete", False),
        ("real_data_source_count", 1),
        ("private_truth_reference_in_pipeline_command_count", 1),
        ("calibration_validation_identifier_overlap_count", 1),
        ("run_overwrite_detected", True),
    ],
)
def test_integrity_violations_are_rejected(
    integrity_field: str,
    invalid_value: bool | int,
) -> None:
    manifest = _valid_manifest()
    manifest["integrity"][integrity_field] = invalid_value

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("phase", "role", "individual_count"),
    [
        ("micro_panel", "TECHNICAL", 23),
        ("pilot_calibration", "VALIDATION", 500),
        ("pilot_validation", "VALIDATION", 499),
        ("qualification_calibration", "CALIBRATION", 2499),
        ("qualification_validation", "CALIBRATION", 2500),
    ],
)
def test_phase_panel_role_and_size_must_match(
    phase: str,
    role: str,
    individual_count: int,
) -> None:
    manifest = _valid_manifest()
    manifest["phase"] = phase
    manifest["panel"]["role"] = role
    manifest["panel"]["individual_count"] = individual_count

    assert _errors(manifest)


def test_invalidated_run_requires_an_invalidating_deviation() -> None:
    manifest = _valid_manifest()
    manifest["run_status"] = "INVALIDATED"

    assert _errors(manifest)


def test_invalidating_deviation_cannot_be_retained() -> None:
    manifest = _valid_manifest()
    manifest["deviations"] = [
        {
            "deviation_id": "deviation_001",
            "recorded_at": "2026-10-28T12:01:00Z",
            "justification": "Écart synthétique de test.",
            "impact": "INVALIDATING",
            "decision": "RETAIN",
        }
    ]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "invalid_timestamp",
    ["2026-10-28", "not-a-date", ""],
)
def test_timestamps_use_date_time_format(invalid_timestamp: str) -> None:
    manifest = _valid_manifest()
    manifest["completed_at"] = invalid_timestamp

    assert _errors(manifest)


def test_unknown_manifest_key_is_rejected() -> None:
    manifest = _valid_manifest()
    manifest["oral_instruction"] = "non reproductible"

    assert _errors(manifest)
