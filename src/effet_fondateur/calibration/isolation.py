"""Gardes d'isolation des panels et runs de calibration synthétiques.

Le module valide uniquement la structure, les métadonnées de fichiers et les
références du manifeste. Il n'ouvre jamais le contenu d'un artefact : un chemin
interdit ou un lien symbolique est donc rejeté avant toute lecture de données.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
import re
from typing import Any


_FORBIDDEN_REAL_COMPONENT = "dock6_real"
_PRIVATE_TRUTH_COMPONENT = "private_truth"
_URI_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
_WINDOWS_DRIVE = re.compile(r"^[A-Za-z]:[\\/]")
_PRIVATE_ARTIFACT_TYPES = {"TRUTH_MANIFEST", "PRIVATE_TRUTH"}


class CalibrationIsolationError(RuntimeError):
    """Signale une violation bloquante avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        message = code if detail is None else f"{code}: {detail}"
        super().__init__(message)


@dataclass(frozen=True)
class IsolationAudit:
    """Résumé non sensible d'un contrôle réussi de l'isolation d'un run."""

    calibration_run_id: str
    artifact_count: int
    command_count: int
    pipeline_input_count: int
    checks: tuple[str, ...]


@dataclass(frozen=True)
class PanelSeparationAudit:
    """Résumé d'une vérification d'indépendance calibration–validation."""

    calibration_run_id: str
    validation_run_id: str
    calibration_identifier_count: int
    validation_identifier_count: int
    shared_public_source_count: int


def _fail(code: str, detail: str | None = None) -> None:
    raise CalibrationIsolationError(code, detail)


def _as_records(manifest: Mapping[str, Any], field: str) -> list[Mapping[str, Any]]:
    records = manifest.get(field)
    if not isinstance(records, list) or any(
        not isinstance(record, Mapping) for record in records
    ):
        _fail("manifest_structure_invalid", field)
    return records


def _index_records(
    records: Sequence[Mapping[str, Any]],
    identifier_field: str,
    namespace: str,
) -> dict[str, Mapping[str, Any]]:
    indexed: dict[str, Mapping[str, Any]] = {}
    normalized: set[str] = set()
    for record in records:
        identifier = record.get(identifier_field)
        if not isinstance(identifier, str) or not identifier:
            _fail("manifest_identifier_invalid", namespace)
        collision_key = identifier.casefold()
        if collision_key in normalized:
            _fail(f"duplicate_{namespace}_id", identifier)
        normalized.add(collision_key)
        indexed[identifier] = record
    return indexed


def _relative_parts(raw_path: Any, field: str) -> tuple[str, ...]:
    if not isinstance(raw_path, str) or not raw_path or "\x00" in raw_path:
        _fail("unsafe_relative_path", field)
    if "\\" in raw_path or raw_path.startswith(("/", "~")):
        _fail("unsafe_relative_path", field)
    if _URI_SCHEME.match(raw_path) or _WINDOWS_DRIVE.match(raw_path):
        _fail("unsafe_relative_path", field)
    raw_parts = raw_path.split("/")
    if any(part in {"", ".", ".."} for part in raw_parts):
        _fail("unsafe_relative_path", field)
    declared_path = PurePosixPath(raw_path)
    if declared_path.is_absolute():
        _fail("unsafe_relative_path", field)
    folded_parts = tuple(part.casefold() for part in declared_path.parts)
    if any(_FORBIDDEN_REAL_COMPONENT in part for part in folded_parts):
        _fail("forbidden_real_data_path", field)
    if any(
        folded_parts[index : index + 2] == ("data", "input")
        for index in range(max(0, len(folded_parts) - 1))
    ):
        _fail("forbidden_real_data_path", field)
    return declared_path.parts


def _is_within(candidate: Path, root: Path) -> bool:
    try:
        candidate.relative_to(root)
    except ValueError:
        return False
    return True


def _reject_symlink_chain(candidate: Path, root: Path, detail: str) -> None:
    relative = candidate.relative_to(root)
    current = root
    if current.is_symlink():
        _fail("symlink_forbidden", detail)
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            _fail("symlink_forbidden", detail)


def _resolve_declared_file(
    run_root: Path,
    raw_path: Any,
    detail: str,
    forbidden_roots: Sequence[Path],
) -> Path:
    parts = _relative_parts(raw_path, detail)
    candidate = run_root.joinpath(*parts)
    _reject_symlink_chain(candidate, run_root, detail)
    resolved_candidate = candidate.resolve(strict=False)
    resolved_run_root = run_root.resolve(strict=True)
    if not _is_within(resolved_candidate, resolved_run_root):
        _fail("path_escapes_run_root", detail)
    if any(_is_within(resolved_candidate, root) for root in forbidden_roots):
        _fail("forbidden_real_data_path", detail)
    if not candidate.exists():
        _fail("declared_artifact_missing", detail)
    if not candidate.is_file():
        _fail("declared_artifact_not_file", detail)
    if candidate.stat(follow_symlinks=False).st_nlink != 1:
        _fail("hardlink_forbidden", detail)
    return resolved_candidate


def _single_artifact_by_type(
    artifacts: Sequence[Mapping[str, Any]],
    artifact_type: str,
) -> Mapping[str, Any]:
    matches = [
        artifact
        for artifact in artifacts
        if artifact.get("artifact_type") == artifact_type
    ]
    if len(matches) != 1:
        _fail("artifact_type_cardinality_invalid", artifact_type)
    return matches[0]


def _assert_contract_matches_artifact(
    contract: Any,
    artifact: Mapping[str, Any],
    detail: str,
) -> None:
    if not isinstance(contract, Mapping):
        _fail("manifest_structure_invalid", detail)
    if (
        contract.get("relative_path") != artifact.get("relative_path")
        or contract.get("sha256") != artifact.get("sha256")
    ):
        _fail("contract_artifact_mismatch", detail)


def _argument_values(argument: str) -> tuple[str, ...]:
    if argument.startswith("-") and "=" in argument:
        return argument, argument.split("=", 1)[1]
    return (argument,)


def _validate_command_argument(argument: Any, scope: str, command_id: str) -> None:
    if not isinstance(argument, str) or not argument or "\x00" in argument:
        _fail("unsafe_command_argument", command_id)
    for value in _argument_values(argument):
        folded = value.replace("\\", "/").casefold()
        if _FORBIDDEN_REAL_COMPONENT in folded or "data/input" in folded:
            _fail("forbidden_real_data_path", command_id)
        if scope == "PIPELINE" and _PRIVATE_TRUTH_COMPONENT in folded:
            _fail("private_truth_exposed_to_pipeline", command_id)
        if (
            value.startswith(("/", "~"))
            or _URI_SCHEME.match(value)
            or _WINDOWS_DRIVE.match(value)
            or ".." in PurePosixPath(folded).parts
        ):
            _fail("unsafe_command_argument", command_id)


def _validate_workspace_roots(
    manifest: Mapping[str, Any],
    run_root: Path,
    workspace_root: Path,
    repository_root: Path,
) -> tuple[Path, tuple[Path, ...]]:
    for path, code in (
        (workspace_root, "workspace_root_invalid"),
        (run_root, "run_root_invalid"),
        (repository_root, "repository_root_invalid"),
    ):
        if path.is_symlink() or not path.is_dir():
            _fail(code)

    resolved_workspace = workspace_root.resolve(strict=True)
    resolved_run = run_root.resolve(strict=True)
    resolved_repository = repository_root.resolve(strict=True)
    if resolved_run.parent != resolved_workspace:
        _fail("run_outside_workspace")
    run_id = manifest.get("calibration_run_id")
    if not isinstance(run_id, str) or run_root.name != run_id:
        _fail("run_directory_id_mismatch")

    forbidden_roots = (
        (resolved_repository / "data" / "input").resolve(strict=False),
        (resolved_repository / "data" / "output").resolve(strict=False),
        (resolved_repository / "data" / "runs").resolve(strict=False),
    )
    if any(_is_within(resolved_workspace, root) for root in forbidden_roots):
        _fail("workspace_under_forbidden_data_root")
    return resolved_run, forbidden_roots


def _validate_artifact_paths(
    artifacts: Sequence[Mapping[str, Any]],
    run_root: Path,
    forbidden_roots: Sequence[Path],
) -> None:
    normalized_paths: set[str] = set()
    for artifact in artifacts:
        artifact_id = str(artifact.get("artifact_id", "artifact"))
        raw_path = artifact.get("relative_path")
        parts = _relative_parts(raw_path, artifact_id)
        collision_key = "/".join(parts).casefold()
        if collision_key in normalized_paths:
            _fail("duplicate_artifact_path", artifact_id)
        normalized_paths.add(collision_key)
        _resolve_declared_file(
            run_root,
            raw_path,
            artifact_id,
            forbidden_roots,
        )


def _validate_cross_references(
    manifest: Mapping[str, Any],
    artifacts: Sequence[Mapping[str, Any]],
    artifacts_by_id: Mapping[str, Mapping[str, Any]],
    sources_by_id: Mapping[str, Mapping[str, Any]],
    tools_by_id: Mapping[str, Mapping[str, Any]],
    commands_by_id: Mapping[str, Mapping[str, Any]],
) -> int:
    pipeline_input_ids: set[str] = set()

    for source_id, source in sources_by_id.items():
        artifact_id = source.get("local_artifact_id")
        artifact = artifacts_by_id.get(artifact_id)
        if artifact is None:
            _fail("source_artifact_reference_missing", source_id)
        if (
            source.get("sha256") != artifact.get("sha256")
            or source.get("size_bytes") != artifact.get("size_bytes")
        ):
            _fail("source_artifact_metadata_mismatch", source_id)

    for command_id, command in commands_by_id.items():
        tool_id = command.get("tool_id")
        if tool_id not in tools_by_id:
            _fail("command_tool_reference_missing", command_id)
        scope = command.get("scope")
        if scope not in {"GENERATOR", "PIPELINE", "EVALUATOR", "REPORTER"}:
            _fail("command_scope_invalid", command_id)
        arguments = command.get("normalized_argv")
        if not isinstance(arguments, list):
            _fail("manifest_structure_invalid", command_id)
        for argument in arguments:
            _validate_command_argument(argument, str(scope), command_id)

        for field in ("input_artifact_ids", "output_artifact_ids"):
            references = command.get(field)
            if not isinstance(references, list):
                _fail("manifest_structure_invalid", command_id)
            for artifact_id in references:
                if artifact_id not in artifacts_by_id:
                    _fail("command_artifact_reference_missing", command_id)

        input_ids = set(command["input_artifact_ids"])
        output_ids = set(command["output_artifact_ids"])
        if input_ids & output_ids:
            _fail("command_input_output_collision", command_id)

        if scope == "PIPELINE":
            for artifact_id in input_ids:
                artifact = artifacts_by_id[artifact_id]
                if (
                    artifact.get("artifact_type") != "BLINDED_INPUT"
                    or artifact.get("pipeline_access_allowed") is not True
                ):
                    _fail("pipeline_input_not_blinded", command_id)
                pipeline_input_ids.add(artifact_id)
        if scope == "REPORTER":
            if any(
                artifacts_by_id[artifact_id].get("artifact_type")
                in _PRIVATE_ARTIFACT_TYPES
                for artifact_id in input_ids
            ):
                _fail("private_truth_exposed_to_reporter", command_id)

        for artifact_id in output_ids:
            producer = artifacts_by_id[artifact_id].get("producer_command_id")
            if producer is not None and producer not in commands_by_id:
                _fail("artifact_producer_reference_missing", artifact_id)
            if producer != command_id:
                _fail("artifact_producer_mismatch", artifact_id)
        for log_field in ("stdout_artifact_id", "stderr_artifact_id"):
            log_id = command.get(log_field)
            if log_id is not None and log_id not in artifacts_by_id:
                _fail("command_log_reference_missing", command_id)

    for artifact in artifacts:
        artifact_id = str(artifact.get("artifact_id", "artifact"))
        producer = artifact.get("producer_command_id")
        if producer is not None:
            command = commands_by_id.get(producer)
            if command is None:
                _fail("artifact_producer_reference_missing", artifact_id)
            if artifact_id not in command.get("output_artifact_ids", []):
                _fail("artifact_producer_mismatch", artifact_id)
        if (
            artifact.get("artifact_type") in _PRIVATE_ARTIFACT_TYPES
            and artifact.get("pipeline_access_allowed") is not False
        ):
            _fail("private_truth_exposed_to_pipeline", artifact_id)

    replay = manifest.get("replay")
    if not isinstance(replay, Mapping):
        _fail("manifest_structure_invalid", "replay")
    replay_order = replay.get("command_order")
    if not isinstance(replay_order, list) or set(replay_order) != set(commands_by_id):
        _fail("replay_command_order_incomplete")

    integrity = manifest.get("integrity")
    if not isinstance(integrity, Mapping):
        _fail("manifest_structure_invalid", "integrity")
    if integrity.get("artifact_count") != len(artifacts):
        _fail("artifact_count_mismatch")
    checksum_id = integrity.get("checksum_index_artifact_id")
    checksum_artifact = artifacts_by_id.get(checksum_id)
    if (
        checksum_artifact is None
        or checksum_artifact.get("artifact_type") != "CHECKSUM_INDEX"
        or replay.get("expected_checksum_index_sha256")
        != checksum_artifact.get("sha256")
    ):
        _fail("checksum_index_reference_invalid")
    return len(pipeline_input_ids)


def _validate_anchored_files(
    manifest: Mapping[str, Any],
    artifacts: Sequence[Mapping[str, Any]],
) -> None:
    config_artifact = _single_artifact_by_type(artifacts, "CONFIGURATION")
    truth_manifest_artifact = _single_artifact_by_type(artifacts, "TRUTH_MANIFEST")
    environment_artifact = _single_artifact_by_type(artifacts, "ENVIRONMENT_LOCK")
    _single_artifact_by_type(artifacts, "CHECKSUM_INDEX")

    contracts = manifest.get("contracts")
    if not isinstance(contracts, Mapping):
        _fail("manifest_structure_invalid", "contracts")
    _assert_contract_matches_artifact(
        contracts.get("configuration"),
        config_artifact,
        "configuration",
    )
    _assert_contract_matches_artifact(
        contracts.get("truth_manifest"),
        truth_manifest_artifact,
        "truth_manifest",
    )

    provenance = manifest.get("provenance")
    if not isinstance(provenance, Mapping):
        _fail("manifest_structure_invalid", "provenance")
    if (
        provenance.get("dependency_lock_path")
        != environment_artifact.get("relative_path")
        or provenance.get("dependency_lock_sha256")
        != environment_artifact.get("sha256")
    ):
        _fail("dependency_lock_artifact_mismatch")

    seed_registry = manifest.get("random_seeds")
    if not isinstance(seed_registry, Mapping):
        _fail("manifest_structure_invalid", "random_seeds")
    seed_matches = [
        artifact
        for artifact in artifacts
        if artifact.get("artifact_type") == "PRIVATE_TRUTH"
        and artifact.get("relative_path") == seed_registry.get("registry_path")
        and artifact.get("sha256") == seed_registry.get("registry_sha256")
    ]
    if len(seed_matches) != 1:
        _fail("seed_registry_artifact_mismatch")


def validate_run_isolation(
    manifest: Mapping[str, Any],
    *,
    run_root: Path,
    workspace_root: Path,
    repository_root: Path,
) -> IsolationAudit:
    """Valide l'isolation physique et référentielle d'un run finalisé.

    Les chemins sont contrôlés lexicalement puis physiquement sans ouvrir leur
    contenu. Toute référence à ``data/input`` ou ``dock6_real``, directe ou via
    un lien symbolique, provoque une erreur bloquante.
    """
    resolved_run_root, forbidden_roots = _validate_workspace_roots(
        manifest,
        run_root,
        workspace_root,
        repository_root,
    )
    artifacts = _as_records(manifest, "artifacts")
    sources = _as_records(manifest, "sources")
    tools = _as_records(manifest, "tools")
    commands = _as_records(manifest, "commands")
    deviations = _as_records(manifest, "deviations")

    artifacts_by_id = _index_records(artifacts, "artifact_id", "artifact")
    sources_by_id = _index_records(sources, "source_id", "source")
    tools_by_id = _index_records(tools, "tool_id", "tool")
    commands_by_id = _index_records(commands, "command_id", "command")
    _index_records(deviations, "deviation_id", "deviation")
    _seed_values(manifest, "run")

    _validate_artifact_paths(
        artifacts,
        resolved_run_root,
        forbidden_roots,
    )
    _validate_anchored_files(manifest, artifacts)
    pipeline_input_count = _validate_cross_references(
        manifest,
        artifacts,
        artifacts_by_id,
        sources_by_id,
        tools_by_id,
        commands_by_id,
    )
    run_id = manifest.get("calibration_run_id")
    if not isinstance(run_id, str):
        _fail("manifest_identifier_invalid", "calibration_run")
    return IsolationAudit(
        calibration_run_id=run_id,
        artifact_count=len(artifacts),
        command_count=len(commands),
        pipeline_input_count=pipeline_input_count,
        checks=(
            "WORKSPACE_CONFINED",
            "NO_SYMLINK",
            "NO_REAL_DATA_PATH",
            "IDENTIFIERS_UNIQUE",
            "REFERENCES_RESOLVED",
            "PRIVATE_TRUTH_ISOLATED",
        ),
    )


def _normalized_identifiers(values: Iterable[str], namespace: str) -> set[str]:
    normalized: set[str] = set()
    for value in values:
        if not isinstance(value, str) or not value:
            _fail("panel_identifier_invalid", namespace)
        folded = value.casefold()
        if folded in normalized:
            _fail("duplicate_identifier_within_panel", namespace)
        normalized.add(folded)
    return normalized


def _seed_values(manifest: Mapping[str, Any], namespace: str) -> set[int]:
    random_seeds = manifest.get("random_seeds")
    values = random_seeds.get("values") if isinstance(random_seeds, Mapping) else None
    if not isinstance(values, Mapping) or any(
        type(value) is not int for value in values.values()
    ):
        _fail("seed_registry_invalid", namespace)
    seed_values = set(values.values())
    if len(seed_values) != len(values):
        _fail("seed_collision_within_panel", namespace)
    return seed_values


def _source_index(manifest: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    sources = _as_records(manifest, "sources")
    return _index_records(sources, "source_id", "source")


def validate_panel_separation(
    calibration_manifest: Mapping[str, Any],
    validation_manifest: Mapping[str, Any],
    *,
    calibration_identifiers: Iterable[str],
    validation_identifiers: Iterable[str],
) -> PanelSeparationAudit:
    """Vérifie l'absence de collisions entre deux panels indépendants."""
    calibration_panel = calibration_manifest.get("panel")
    validation_panel = validation_manifest.get("panel")
    if not isinstance(calibration_panel, Mapping) or not isinstance(
        validation_panel, Mapping
    ):
        _fail("manifest_structure_invalid", "panel")
    if calibration_panel.get("role") != "CALIBRATION":
        _fail("panel_role_invalid", "calibration")
    if validation_panel.get("role") != "VALIDATION":
        _fail("panel_role_invalid", "validation")

    calibration_run_id = calibration_manifest.get("calibration_run_id")
    validation_run_id = validation_manifest.get("calibration_run_id")
    if not isinstance(calibration_run_id, str) or not isinstance(
        validation_run_id, str
    ):
        _fail("manifest_identifier_invalid", "calibration_run")
    if calibration_run_id.casefold() == validation_run_id.casefold():
        _fail("calibration_validation_run_collision")
    calibration_panel_id = calibration_panel.get("panel_id")
    validation_panel_id = validation_panel.get("panel_id")
    if not isinstance(calibration_panel_id, str) or not isinstance(
        validation_panel_id, str
    ):
        _fail("manifest_identifier_invalid", "panel")
    if calibration_panel_id.casefold() == validation_panel_id.casefold():
        _fail("calibration_validation_panel_collision")

    calibration_ids = _normalized_identifiers(
        calibration_identifiers,
        "calibration",
    )
    validation_ids = _normalized_identifiers(
        validation_identifiers,
        "validation",
    )
    if calibration_ids & validation_ids:
        _fail("calibration_validation_identifier_collision")

    calibration_seeds = _seed_values(calibration_manifest, "calibration")
    validation_seeds = _seed_values(validation_manifest, "validation")
    if calibration_seeds & validation_seeds:
        _fail("calibration_validation_seed_collision")

    calibration_sources = _source_index(calibration_manifest)
    validation_sources = _source_index(validation_manifest)
    normalized_calibration_sources = {
        source_id.casefold(): source
        for source_id, source in calibration_sources.items()
    }
    normalized_validation_sources = {
        source_id.casefold(): source
        for source_id, source in validation_sources.items()
    }
    shared_public_source_count = 0
    for source_id in (
        normalized_calibration_sources.keys()
        & normalized_validation_sources.keys()
    ):
        calibration_source = normalized_calibration_sources[source_id]
        validation_source = normalized_validation_sources[source_id]
        public_shared = (
            calibration_source.get("source_type") == "PUBLIC_REFERENCE"
            and validation_source.get("source_type") == "PUBLIC_REFERENCE"
            and calibration_source.get("shared_between_panels") is True
            and validation_source.get("shared_between_panels") is True
            and calibration_source.get("sha256") == validation_source.get("sha256")
        )
        if not public_shared:
            _fail("calibration_validation_source_collision", source_id)
        shared_public_source_count += 1

    calibration_private_hashes = {
        source.get("sha256")
        for source in calibration_sources.values()
        if source.get("source_type") != "PUBLIC_REFERENCE"
    }
    validation_private_hashes = {
        source.get("sha256")
        for source in validation_sources.values()
        if source.get("source_type") != "PUBLIC_REFERENCE"
    }
    if calibration_private_hashes & validation_private_hashes:
        _fail("calibration_validation_source_content_collision")

    return PanelSeparationAudit(
        calibration_run_id=calibration_run_id,
        validation_run_id=validation_run_id,
        calibration_identifier_count=len(calibration_ids),
        validation_identifier_count=len(validation_ids),
        shared_public_source_count=shared_public_source_count,
    )
