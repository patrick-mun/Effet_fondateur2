"""Tests des gardes physiques et référentielles de calibration D08."""

from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from effet_fondateur.calibration.isolation import (
    CalibrationIsolationError,
    IsolationAudit,
    validate_panel_separation,
    validate_run_isolation,
)


def _sha(character: str) -> str:
    return character * 64


def _artifact(
    artifact_id: str,
    artifact_type: str,
    relative_path: str,
    hash_character: str,
    content: str,
    *,
    producer_command_id: str | None,
    pipeline_access_allowed: bool,
) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "artifact_type": artifact_type,
        "relative_path": relative_path,
        "sha256": _sha(hash_character),
        "size_bytes": len(content.encode("utf-8")),
        "producer_command_id": producer_command_id,
        "pipeline_access_allowed": pipeline_access_allowed,
    }


def _build_run(
    tmp_path: Path,
) -> tuple[dict[str, Any], Path, Path, Path, dict[str, str]]:
    repository_root = tmp_path / "repository"
    workspace_root = tmp_path / "workspace"
    run_id = "cal_run_micro_001"
    run_root = workspace_root / run_id
    repository_root.mkdir()
    run_root.mkdir(parents=True)

    contents = {
        "control/config.resolved.json": "{}\n",
        "private_truth/truth_manifest.json": "{}\n",
        "private_truth/seed_registry.json": "{}\n",
        "blinded_inputs/panel.vcf.gz": "synthetic-blinded\n",
        "control/environment.freeze.txt": "python==3.12.7\n",
        "checksums.sha256": "synthetic-checksums\n",
    }
    for relative_path, content in contents.items():
        path = run_root / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    artifacts = [
        _artifact(
            "calibration_config",
            "CONFIGURATION",
            "control/config.resolved.json",
            "1",
            contents["control/config.resolved.json"],
            producer_command_id=None,
            pipeline_access_allowed=False,
        ),
        _artifact(
            "truth_manifest",
            "TRUTH_MANIFEST",
            "private_truth/truth_manifest.json",
            "2",
            contents["private_truth/truth_manifest.json"],
            producer_command_id="generate_panel",
            pipeline_access_allowed=False,
        ),
        _artifact(
            "seed_registry",
            "PRIVATE_TRUTH",
            "private_truth/seed_registry.json",
            "3",
            contents["private_truth/seed_registry.json"],
            producer_command_id="generate_panel",
            pipeline_access_allowed=False,
        ),
        _artifact(
            "blinded_panel",
            "BLINDED_INPUT",
            "blinded_inputs/panel.vcf.gz",
            "4",
            contents["blinded_inputs/panel.vcf.gz"],
            producer_command_id="generate_panel",
            pipeline_access_allowed=True,
        ),
        _artifact(
            "environment_lock",
            "ENVIRONMENT_LOCK",
            "control/environment.freeze.txt",
            "5",
            contents["control/environment.freeze.txt"],
            producer_command_id=None,
            pipeline_access_allowed=False,
        ),
        _artifact(
            "checksum_index",
            "CHECKSUM_INDEX",
            "checksums.sha256",
            "6",
            contents["checksums.sha256"],
            producer_command_id="generate_panel",
            pipeline_access_allowed=False,
        ),
    ]
    generator_outputs = [
        "truth_manifest",
        "seed_registry",
        "blinded_panel",
        "checksum_index",
    ]
    manifest: dict[str, Any] = {
        "calibration_run_id": run_id,
        "panel": {
            "panel_id": "micro_technical_001",
            "role": "TECHNICAL",
        },
        "provenance": {
            "dependency_lock_path": "control/environment.freeze.txt",
            "dependency_lock_sha256": _sha("5"),
        },
        "contracts": {
            "configuration": {
                "relative_path": "control/config.resolved.json",
                "sha256": _sha("1"),
            },
            "truth_manifest": {
                "relative_path": "private_truth/truth_manifest.json",
                "sha256": _sha("2"),
            },
        },
        "random_seeds": {
            "registry_path": "private_truth/seed_registry.json",
            "registry_sha256": _sha("3"),
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
                "source_id": "synthetic_fixture_v1",
                "source_type": "SYNTHETIC_FIXTURE",
                "local_artifact_id": "blinded_panel",
                "sha256": _sha("4"),
                "size_bytes": len(
                    contents["blinded_inputs/panel.vcf.gz"].encode("utf-8")
                ),
                "shared_between_panels": False,
            }
        ],
        "tools": [{"tool_id": "python"}],
        "commands": [
            {
                "command_id": "generate_panel",
                "scope": "GENERATOR",
                "tool_id": "python",
                "normalized_argv": [
                    ".venv/bin/python",
                    "-m",
                    "effet_fondateur.calibration.generate",
                ],
                "input_artifact_ids": ["calibration_config"],
                "output_artifact_ids": generator_outputs,
                "stdout_artifact_id": None,
                "stderr_artifact_id": None,
            }
        ],
        "artifacts": artifacts,
        "replay": {
            "command_order": ["generate_panel"],
            "expected_checksum_index_sha256": _sha("6"),
        },
        "integrity": {
            "artifact_count": len(artifacts),
            "checksum_index_artifact_id": "checksum_index",
        },
        "deviations": [],
    }
    return manifest, run_root, workspace_root, repository_root, contents


def _validate(
    manifest: dict[str, Any],
    run_root: Path,
    workspace_root: Path,
    repository_root: Path,
) -> IsolationAudit:
    return validate_run_isolation(
        manifest,
        run_root=run_root,
        workspace_root=workspace_root,
        repository_root=repository_root,
    )


def _assert_error(expected_code: str, action: Callable[[], Any]) -> None:
    with pytest.raises(CalibrationIsolationError) as captured:
        action()
    assert captured.value.code == expected_code


def _add_pipeline_command(
    manifest: dict[str, Any],
    run_root: Path,
) -> dict[str, Any]:
    content = "synthetic-result\n"
    relative_path = "pipeline_runs/result.json"
    output_path = run_root / relative_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content, encoding="utf-8")
    artifact = _artifact(
        "pipeline_result",
        "PIPELINE_OUTPUT",
        relative_path,
        "7",
        content,
        producer_command_id="run_pipeline",
        pipeline_access_allowed=False,
    )
    command = {
        "command_id": "run_pipeline",
        "scope": "PIPELINE",
        "tool_id": "python",
        "normalized_argv": [
            ".venv/bin/effet-fondateur",
            "run",
            "--config",
            "blinded_inputs/pipeline.yaml",
        ],
        "input_artifact_ids": ["blinded_panel"],
        "output_artifact_ids": ["pipeline_result"],
        "stdout_artifact_id": None,
        "stderr_artifact_id": None,
    }
    manifest["artifacts"].append(artifact)
    manifest["commands"].append(command)
    manifest["replay"]["command_order"].append("run_pipeline")
    manifest["integrity"]["artifact_count"] += 1
    return command


def test_valid_run_is_confined_without_opening_artifact_content(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)

    def forbidden_open(*args: Any, **kwargs: Any) -> None:
        raise AssertionError("Le garde D08 ne doit ouvrir aucun artefact.")

    monkeypatch.setattr(Path, "open", forbidden_open)
    audit = _validate(manifest, run_root, workspace_root, repository_root)

    assert audit.artifact_count == 6
    assert audit.command_count == 1
    assert audit.pipeline_input_count == 0
    assert "NO_REAL_DATA_PATH" in audit.checks


@pytest.mark.parametrize(
    "forbidden_path",
    [
        "../../repository/data/input/dock6_real/secret.vcf",
        "data/input/dock6_real/secret.vcf",
        "dock6_real/secret.vcf",
        "/tmp/secret.vcf",
        "file:///tmp/secret.vcf",
        "C:\\temp\\secret.vcf",
    ],
)
def test_real_or_nonportable_artifact_paths_are_rejected_before_open(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    forbidden_path: str,
) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    manifest["artifacts"][3]["relative_path"] = forbidden_path
    opened: list[Path] = []

    def record_open(path: Path, *args: Any, **kwargs: Any) -> None:
        opened.append(path)
        raise AssertionError("Un chemin interdit ne doit jamais être ouvert.")

    monkeypatch.setattr(Path, "open", record_open)
    expected = "unsafe_relative_path"
    if not forbidden_path.startswith("../") and (
        "dock6_real" in forbidden_path.casefold()
        or "data/input" in forbidden_path.casefold()
    ):
        expected = "forbidden_real_data_path"

    _assert_error(
        expected,
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )
    assert opened == []


def test_symlink_to_dock6_real_is_rejected_before_open(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    forbidden_root = repository_root / "data" / "input" / "dock6_real"
    forbidden_root.mkdir(parents=True)
    real_file = forbidden_root / "secret.vcf"
    real_file.write_text("synthetic-probe\n", encoding="utf-8")
    blinded_file = run_root / "blinded_inputs" / "panel.vcf.gz"
    blinded_file.unlink()
    blinded_file.symlink_to(real_file)
    opened: list[Path] = []

    def record_open(path: Path, *args: Any, **kwargs: Any) -> None:
        opened.append(path)
        raise AssertionError("Le lien interdit ne doit pas être ouvert.")

    monkeypatch.setattr(Path, "open", record_open)

    _assert_error(
        "symlink_forbidden",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )
    assert opened == []


def test_symlinked_parent_directory_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    forbidden_root = repository_root / "data" / "input" / "dock6_real"
    forbidden_root.mkdir(parents=True)
    (forbidden_root / "panel.vcf.gz").write_text("probe\n", encoding="utf-8")
    blinded_file = run_root / "blinded_inputs" / "panel.vcf.gz"
    blinded_file.unlink()
    blinded_file.parent.rmdir()
    blinded_file.parent.symlink_to(forbidden_root, target_is_directory=True)

    _assert_error(
        "symlink_forbidden",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_symlinked_run_root_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    real_run_root = tmp_path / "real_run"
    run_root.rename(real_run_root)
    run_root.symlink_to(real_run_root, target_is_directory=True)

    _assert_error(
        "run_root_invalid",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_symlinked_workspace_root_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    real_workspace = tmp_path / "real_workspace"
    workspace_root.rename(real_workspace)
    workspace_root.symlink_to(real_workspace, target_is_directory=True)

    _assert_error(
        "workspace_root_invalid",
        lambda: _validate(
            manifest,
            workspace_root / run_root.name,
            workspace_root,
            repository_root,
        ),
    )


def test_hardlink_to_real_input_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    forbidden_root = repository_root / "data" / "input" / "dock6_real"
    forbidden_root.mkdir(parents=True)
    real_file = forbidden_root / "secret.vcf"
    real_file.write_text("synthetic-probe\n", encoding="utf-8")
    blinded_file = run_root / "blinded_inputs" / "panel.vcf.gz"
    blinded_file.unlink()
    blinded_file.hardlink_to(real_file)

    _assert_error(
        "hardlink_forbidden",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_run_outside_declared_workspace_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, _, repository_root, _ = _build_run(tmp_path)
    other_workspace = tmp_path / "other_workspace"
    other_workspace.mkdir()

    _assert_error(
        "run_outside_workspace",
        lambda: _validate(manifest, run_root, other_workspace, repository_root),
    )


def test_run_directory_must_match_run_identifier(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    manifest["calibration_run_id"] = "another_run"

    _assert_error(
        "run_directory_id_mismatch",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_missing_artifact_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    (run_root / "blinded_inputs" / "panel.vcf.gz").unlink()

    _assert_error(
        "declared_artifact_missing",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_directory_cannot_replace_a_declared_file(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    checksum_path = run_root / "checksums.sha256"
    checksum_path.unlink()
    checksum_path.mkdir()

    _assert_error(
        "declared_artifact_not_file",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_duplicate_artifact_identifier_is_rejected(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    duplicate = deepcopy(manifest["artifacts"][0])
    duplicate["relative_path"] = "control/duplicate.json"
    (run_root / "control" / "duplicate.json").write_text("{}\n", encoding="utf-8")
    manifest["artifacts"].append(duplicate)

    _assert_error(
        "duplicate_artifact_id",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_artifact_path_collision_is_case_insensitive(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    collision = deepcopy(manifest["artifacts"][0])
    collision["artifact_id"] = "case_collision"
    collision["relative_path"] = "CONTROL/CONFIG.RESOLVED.JSON"
    manifest["artifacts"].append(collision)

    _assert_error(
        "duplicate_artifact_path",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


@pytest.mark.parametrize(
    ("collection", "identifier_field", "expected_code"),
    [
        ("sources", "source_id", "duplicate_source_id"),
        ("tools", "tool_id", "duplicate_tool_id"),
        ("commands", "command_id", "duplicate_command_id"),
    ],
)
def test_manifest_identifiers_are_unique_within_each_namespace(
    tmp_path: Path,
    collection: str,
    identifier_field: str,
    expected_code: str,
) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    duplicate = deepcopy(manifest[collection][0])
    duplicate[identifier_field] = duplicate[identifier_field].upper()
    manifest[collection].append(duplicate)

    _assert_error(
        expected_code,
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        ("missing_source_artifact", "source_artifact_reference_missing"),
        ("source_hash_mismatch", "source_artifact_metadata_mismatch"),
        ("missing_command_tool", "command_tool_reference_missing"),
        ("missing_command_input", "command_artifact_reference_missing"),
        ("missing_command_output", "command_artifact_reference_missing"),
        ("missing_producer", "artifact_producer_reference_missing"),
        ("producer_mismatch", "artifact_producer_mismatch"),
        ("incomplete_replay", "replay_command_order_incomplete"),
        ("artifact_count", "artifact_count_mismatch"),
        ("checksum_reference", "checksum_index_reference_invalid"),
        ("config_contract", "contract_artifact_mismatch"),
        ("truth_contract", "contract_artifact_mismatch"),
        ("dependency_lock", "dependency_lock_artifact_mismatch"),
        ("seed_registry", "seed_registry_artifact_mismatch"),
    ],
)
def test_cross_reference_failures_are_blocking(
    tmp_path: Path,
    mutation: str,
    expected_code: str,
) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    if mutation == "missing_source_artifact":
        manifest["sources"][0]["local_artifact_id"] = "missing"
    elif mutation == "source_hash_mismatch":
        manifest["sources"][0]["sha256"] = _sha("f")
    elif mutation == "missing_command_tool":
        manifest["commands"][0]["tool_id"] = "missing"
    elif mutation == "missing_command_input":
        manifest["commands"][0]["input_artifact_ids"] = ["missing"]
    elif mutation == "missing_command_output":
        manifest["commands"][0]["output_artifact_ids"].append("missing")
    elif mutation == "missing_producer":
        manifest["artifacts"][1]["producer_command_id"] = "missing"
    elif mutation == "producer_mismatch":
        manifest["commands"][0]["output_artifact_ids"].remove("truth_manifest")
    elif mutation == "incomplete_replay":
        manifest["replay"]["command_order"] = []
    elif mutation == "artifact_count":
        manifest["integrity"]["artifact_count"] = 99
    elif mutation == "checksum_reference":
        manifest["integrity"]["checksum_index_artifact_id"] = "truth_manifest"
    elif mutation == "config_contract":
        manifest["contracts"]["configuration"]["sha256"] = _sha("f")
    elif mutation == "truth_contract":
        manifest["contracts"]["truth_manifest"]["relative_path"] = (
            "private_truth/other.json"
        )
    elif mutation == "dependency_lock":
        manifest["provenance"]["dependency_lock_sha256"] = _sha("f")
    elif mutation == "seed_registry":
        manifest["random_seeds"]["registry_sha256"] = _sha("f")

    _assert_error(
        expected_code,
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_pipeline_may_only_consume_blinded_inputs(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    command = _add_pipeline_command(manifest, run_root)
    command["input_artifact_ids"] = ["truth_manifest"]

    _assert_error(
        "pipeline_input_not_blinded",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_pipeline_command_with_only_blinded_inputs_is_valid(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    _add_pipeline_command(manifest, run_root)

    audit = _validate(manifest, run_root, workspace_root, repository_root)

    assert audit.pipeline_input_count == 1
    assert audit.command_count == 2


def test_pipeline_command_cannot_name_private_truth(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    command = _add_pipeline_command(manifest, run_root)
    command["normalized_argv"].append("private_truth/truth_manifest.json")

    _assert_error(
        "private_truth_exposed_to_pipeline",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


@pytest.mark.parametrize(
    "forbidden_argument",
    [
        "data/input/dock6_real/panel.vcf",
        "--input=/tmp/panel.vcf",
        "--input=../panel.vcf",
        "file:///tmp/panel.vcf",
    ],
)
def test_command_arguments_cannot_bypass_path_isolation(
    tmp_path: Path,
    forbidden_argument: str,
) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    manifest["commands"][0]["normalized_argv"].append(forbidden_argument)
    expected = (
        "forbidden_real_data_path"
        if "data/input" in forbidden_argument
        else "unsafe_command_argument"
    )

    _assert_error(
        expected,
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_private_truth_cannot_be_marked_for_pipeline_access(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    manifest["artifacts"][1]["pipeline_access_allowed"] = True

    _assert_error(
        "private_truth_exposed_to_pipeline",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def test_seeds_must_be_distinct_within_a_run(tmp_path: Path) -> None:
    manifest, run_root, workspace_root, repository_root, _ = _build_run(tmp_path)
    manifest["random_seeds"]["values"]["haplotypes"] = manifest[
        "random_seeds"
    ]["values"]["master"]

    _assert_error(
        "seed_collision_within_panel",
        lambda: _validate(manifest, run_root, workspace_root, repository_root),
    )


def _panel_pair(tmp_path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    calibration, _, _, _, _ = _build_run(tmp_path)
    validation = deepcopy(calibration)
    calibration["panel"] = {
        "panel_id": "pilot_calibration_001",
        "role": "CALIBRATION",
    }
    calibration["calibration_run_id"] = "calibration_run_001"
    validation["panel"] = {
        "panel_id": "pilot_validation_001",
        "role": "VALIDATION",
    }
    validation["calibration_run_id"] = "validation_run_001"
    validation["random_seeds"]["values"] = {
        name: value + 100
        for name, value in calibration["random_seeds"]["values"].items()
    }
    validation["sources"][0]["source_id"] = "validation_fixture_v1"
    validation["sources"][0]["sha256"] = _sha("a")
    return calibration, validation


def test_independent_panels_pass_separation_controls(tmp_path: Path) -> None:
    calibration, validation = _panel_pair(tmp_path)

    audit = validate_panel_separation(
        calibration,
        validation,
        calibration_identifiers=["CAL_001", "CAL_002"],
        validation_identifiers=["VAL_001", "VAL_002"],
    )

    assert audit.calibration_identifier_count == 2
    assert audit.validation_identifier_count == 2
    assert audit.shared_public_source_count == 0


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        ("calibration_role", "panel_role_invalid"),
        ("validation_role", "panel_role_invalid"),
        ("run_id", "calibration_validation_run_collision"),
        ("panel_id", "calibration_validation_panel_collision"),
        ("seed_cross", "calibration_validation_seed_collision"),
        ("seed_within", "seed_collision_within_panel"),
        ("source_id", "calibration_validation_source_collision"),
        ("source_content", "calibration_validation_source_content_collision"),
    ],
)
def test_panel_separation_rejects_collisions(
    tmp_path: Path,
    mutation: str,
    expected_code: str,
) -> None:
    calibration, validation = _panel_pair(tmp_path)
    if mutation == "calibration_role":
        calibration["panel"]["role"] = "TECHNICAL"
    elif mutation == "validation_role":
        validation["panel"]["role"] = "CALIBRATION"
    elif mutation == "run_id":
        validation["calibration_run_id"] = "CALIBRATION_RUN_001"
    elif mutation == "panel_id":
        validation["panel"]["panel_id"] = "PILOT_CALIBRATION_001"
    elif mutation == "seed_cross":
        validation["random_seeds"]["values"]["master"] = 1001
    elif mutation == "seed_within":
        validation["random_seeds"]["values"]["haplotypes"] = validation[
            "random_seeds"
        ]["values"]["master"]
    elif mutation == "source_id":
        validation["sources"][0]["source_id"] = "SYNTHETIC_FIXTURE_V1"
    elif mutation == "source_content":
        validation["sources"][0]["sha256"] = _sha("4")

    _assert_error(
        expected_code,
        lambda: validate_panel_separation(
            calibration,
            validation,
            calibration_identifiers=["CAL_001"],
            validation_identifiers=["VAL_001"],
        ),
    )


def test_panel_identifiers_are_disjoint_case_insensitively(tmp_path: Path) -> None:
    calibration, validation = _panel_pair(tmp_path)

    _assert_error(
        "calibration_validation_identifier_collision",
        lambda: validate_panel_separation(
            calibration,
            validation,
            calibration_identifiers=["SYNTHETIC_001"],
            validation_identifiers=["synthetic_001"],
        ),
    )


def test_declared_public_reference_may_be_shared(tmp_path: Path) -> None:
    calibration, validation = _panel_pair(tmp_path)
    public_source = {
        "source_id": "public_reference_v1",
        "source_type": "PUBLIC_REFERENCE",
        "sha256": _sha("b"),
        "shared_between_panels": True,
    }
    calibration["sources"].append(deepcopy(public_source))
    validation["sources"].append(deepcopy(public_source))

    audit = validate_panel_separation(
        calibration,
        validation,
        calibration_identifiers=["CAL_001"],
        validation_identifiers=["VAL_001"],
    )

    assert audit.shared_public_source_count == 1
