"""Routeur à privilège minimal pour les vérités privées après clôture."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import Any

from effet_fondateur.calibration.contracts import STAGE_IDS, TRUTH_CONSUMERS


class TruthRoutingError(ValueError):
    """Refus d'accès à une vérité privée."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class TruthGrant:
    """Autorisation bornée à une étape, une release et des artefacts exacts."""

    evaluation_id: str
    stage_id: str
    release_id: str
    closed_run_id: str
    artifacts: tuple[tuple[str, str, str], ...]


def allowed_truth_for_stage(stage_id: str) -> frozenset[str]:
    """Retourne les familles explicitement accordées à un évaluateur."""
    if stage_id not in STAGE_IDS:
        raise TruthRoutingError("unknown_stage", stage_id)
    return frozenset(
        truth_id
        for truth_id, consumers in TRUTH_CONSUMERS.items()
        if stage_id in consumers
    )


def _artifact_index(manifest: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    artifacts = manifest.get("truth_artifacts")
    if not isinstance(artifacts, list):
        raise TruthRoutingError("truth_manifest_invalid")
    indexed: dict[str, Mapping[str, Any]] = {}
    for artifact in artifacts:
        if not isinstance(artifact, Mapping):
            raise TruthRoutingError("truth_manifest_invalid")
        truth_id = artifact.get("truth_id")
        if not isinstance(truth_id, str) or truth_id in indexed:
            raise TruthRoutingError("truth_artifact_id_invalid", str(truth_id))
        indexed[truth_id] = artifact
    return indexed


def authorize_truth_access(
    *,
    evaluation_id: str,
    stage_id: str,
    requested_truth_ids: Sequence[str],
    truth_manifest: Mapping[str, Any],
    run_manifest: Mapping[str, Any],
) -> TruthGrant:
    """Autorise l'accès seulement après clôture et pour les familles nécessaires."""
    if not evaluation_id:
        raise TruthRoutingError("evaluation_id_invalid")
    if run_manifest.get("run_state") != "CLOSED":
        raise TruthRoutingError("run_not_closed")
    if not run_manifest.get("closed_at") or not run_manifest.get("manifest_sha256"):
        raise TruthRoutingError("run_closure_evidence_missing")
    release_id = run_manifest.get("release_id")
    run_id = run_manifest.get("run_id")
    if not isinstance(release_id, str) or not isinstance(run_id, str):
        raise TruthRoutingError("run_identity_invalid")
    if truth_manifest.get("release_id") != release_id:
        raise TruthRoutingError("truth_release_mismatch")
    if truth_manifest.get("disclosure_state") != "SEALED_PRIVATE":
        raise TruthRoutingError("truth_not_sealed")

    requested = tuple(requested_truth_ids)
    if len(set(requested)) != len(requested):
        raise TruthRoutingError("duplicate_truth_request")
    allowed = allowed_truth_for_stage(stage_id)
    if not set(requested).issubset(allowed):
        denied = sorted(set(requested) - allowed)
        raise TruthRoutingError("truth_scope_exceeded", ",".join(denied))
    indexed = _artifact_index(truth_manifest)
    grants: list[tuple[str, str, str]] = []
    for truth_id in requested:
        artifact = indexed.get(truth_id)
        if artifact is None:
            raise TruthRoutingError("truth_artifact_missing", truth_id)
        path = artifact.get("relative_path")
        digest = artifact.get("sha256")
        if not isinstance(path, str) or not isinstance(digest, str) or len(digest) != 64:
            raise TruthRoutingError("truth_artifact_contract_invalid", truth_id)
        parts = PurePosixPath(path).parts
        if not parts or parts[0] != "private_truth" or ".." in parts:
            raise TruthRoutingError("truth_path_invalid", truth_id)
        grants.append((truth_id, path, digest))
    return TruthGrant(
        evaluation_id=evaluation_id,
        stage_id=stage_id,
        release_id=release_id,
        closed_run_id=run_id,
        artifacts=tuple(grants),
    )


def validate_truth_policy(policy: Mapping[str, Any]) -> None:
    """Vérifie qu'une politique sérialisée reproduit exactement R04/R06."""
    expected = {
        stage_id: sorted(allowed_truth_for_stage(stage_id))
        for stage_id in STAGE_IDS
    }
    actual = policy.get("stage_truth_access")
    if actual != expected:
        raise TruthRoutingError("serialized_truth_policy_mismatch")
