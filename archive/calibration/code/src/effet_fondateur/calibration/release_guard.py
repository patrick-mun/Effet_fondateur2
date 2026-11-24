"""Gardes génériques d'isolement et d'aveuglement de la release unifiée."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from hashlib import sha256
import hmac
from pathlib import Path, PurePosixPath
import re
from typing import Any


PRIVATE_COMPONENTS = frozenset({"private_truth", "governance"})
FORBIDDEN_BLINDED_FIELDS = frozenset(
    {
        "master_id",
        "founder_age_true",
        "true_phase",
        "true_segment",
        "founder_class",
        "ibd_class",
        "roh_class",
        "injected_error_mask",
        "expected_verdict",
        "expected_direction",
        "latent_population",
        "master_seed",
        "private_path",
    }
)
_URI_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")


class ReleaseGuardError(ValueError):
    """Violation bloquante d'un invariant R06/R09."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class BlindedInventoryAudit:
    """Résumé non sensible d'une validation de liste blanche."""

    inventory_id: str
    artifact_count: int
    declared_bytes: int
    verified_files: bool


@dataclass(frozen=True)
class SeparationAudit:
    """Résumé des contrôles de disjonction, sans identifiant biologique."""

    partition_count: int
    isolation_unit_count: int
    reference_member_count: int
    negative_control_count: int


def _fail(code: str, detail: str | None = None) -> None:
    raise ReleaseGuardError(code, detail)


def _relative_parts(value: Any, detail: str) -> tuple[str, ...]:
    if not isinstance(value, str) or not value or "\x00" in value:
        _fail("unsafe_relative_path", detail)
    if value.startswith(("/", "~")) or "\\" in value or _URI_SCHEME.match(value):
        _fail("unsafe_relative_path", detail)
    path = PurePosixPath(value)
    if any(part in {"", ".", ".."} for part in path.parts):
        _fail("unsafe_relative_path", detail)
    folded = tuple(part.casefold() for part in path.parts)
    if PRIVATE_COMPONENTS.intersection(folded):
        _fail("private_path_in_blinded_inventory", detail)
    if any(
        folded[index : index + 2] == ("data", "input")
        for index in range(max(0, len(folded) - 1))
    ):
        _fail("real_input_path_forbidden", detail)
    return path.parts


def _within(candidate: Path, root: Path) -> bool:
    try:
        candidate.relative_to(root)
    except ValueError:
        return False
    return True


def _reject_symlink_chain(candidate: Path, root: Path, detail: str) -> None:
    current = root
    if current.is_symlink():
        _fail("symlink_forbidden", detail)
    for part in candidate.relative_to(root).parts:
        current = current / part
        if current.is_symlink():
            _fail("symlink_forbidden", detail)


def scan_blinded_document(value: Any, location: str = "$") -> None:
    """Refuse récursivement les champs et chemins privés dans un document aveugle."""
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str):
                _fail("non_string_blinded_key", location)
            folded_key = key.casefold()
            if folded_key in FORBIDDEN_BLINDED_FIELDS:
                _fail("truth_field_in_blinded_document", f"{location}.{key}")
            if folded_key.endswith(("_path", "_relative_path")) and isinstance(child, str):
                _relative_parts(child, f"{location}.{key}")
            scan_blinded_document(child, f"{location}.{key}")
        return
    if isinstance(value, list):
        for index, child in enumerate(value):
            scan_blinded_document(child, f"{location}[{index}]")
        return
    if isinstance(value, str):
        folded_parts = tuple(part.casefold() for part in value.replace("\\", "/").split("/"))
        if PRIVATE_COMPONENTS.intersection(folded_parts):
            _fail("private_reference_in_blinded_document", location)


def validate_blinded_inventory(
    inventory: Mapping[str, Any],
    *,
    view_root: Path,
    repository_root: Path,
    verify_files: bool,
) -> BlindedInventoryAudit:
    """Valide l'inventaire fermé remis au processus aveugle.

    ``verify_files=False`` est le mode statique R09. R10 et les étapes
    suivantes doivent activer la vérification physique et des empreintes.
    """
    scan_blinded_document(inventory)
    inventory_id = inventory.get("inventory_id")
    artifacts = inventory.get("artifacts")
    if not isinstance(inventory_id, str) or not inventory_id:
        _fail("inventory_id_invalid")
    if not isinstance(artifacts, list) or any(not isinstance(item, Mapping) for item in artifacts):
        _fail("inventory_artifacts_invalid")
    if view_root.is_symlink() or repository_root.is_symlink():
        _fail("root_symlink_forbidden")
    resolved_view = view_root.resolve(strict=verify_files)
    resolved_repository = repository_root.resolve(strict=True)
    real_input_root = (resolved_repository / "data" / "input").resolve(strict=False)
    paths: set[str] = set()
    artifact_ids: set[str] = set()
    declared_bytes = 0
    for artifact in artifacts:
        artifact_id = artifact.get("artifact_id")
        if not isinstance(artifact_id, str) or not artifact_id or artifact_id in artifact_ids:
            _fail("artifact_id_invalid", str(artifact_id))
        artifact_ids.add(artifact_id)
        parts = _relative_parts(artifact.get("relative_path"), artifact_id)
        normalized = "/".join(parts).casefold()
        if normalized in paths:
            _fail("duplicate_blinded_path", artifact_id)
        paths.add(normalized)
        candidate = view_root.joinpath(*parts)
        _reject_symlink_chain(candidate, view_root, artifact_id)
        resolved = candidate.resolve(strict=False)
        if not _within(resolved, resolved_view):
            _fail("path_escapes_view_root", artifact_id)
        if _within(resolved, real_input_root):
            _fail("real_input_path_forbidden", artifact_id)
        size = artifact.get("size_bytes")
        if not isinstance(size, int) or size < 0:
            _fail("artifact_size_invalid", artifact_id)
        declared_bytes += size
        if verify_files:
            if not candidate.is_file():
                _fail("blinded_artifact_missing", artifact_id)
            if candidate.stat(follow_symlinks=False).st_nlink != 1:
                _fail("hardlink_forbidden", artifact_id)
            payload = candidate.read_bytes()
            if len(payload) != size:
                _fail("artifact_size_mismatch", artifact_id)
            if sha256(payload).hexdigest() != artifact.get("sha256"):
                _fail("artifact_hash_mismatch", artifact_id)
    return BlindedInventoryAudit(
        inventory_id=inventory_id,
        artifact_count=len(artifacts),
        declared_bytes=declared_bytes,
        verified_files=verify_files,
    )


def partition_pseudonym(secret: bytes, partition_id: str, master_id: str) -> str:
    """Produit un pseudonyme stable dans une partition et non raccordable entre partitions."""
    if len(secret) < 32:
        _fail("pseudonym_secret_too_short")
    if partition_id not in {"CALIBRATION", "VALIDATION_BLIND", "RESOURCE_ONLY"}:
        _fail("partition_id_invalid", partition_id)
    if not master_id:
        _fail("master_id_invalid")
    digest = hmac.new(
        secret,
        f"Effet_fondateur2|{partition_id}|{master_id}".encode("utf-8"),
        "sha256",
    ).hexdigest()
    return f"SYN-{digest[:24]}"


def _sets_by_partition(
    assignments: Sequence[Mapping[str, Any]],
    field: str,
) -> dict[str, set[str]]:
    values: dict[str, set[str]] = {}
    for assignment in assignments:
        partition = assignment.get("partition_id")
        raw_values = assignment.get(field, [])
        if not isinstance(partition, str) or not isinstance(raw_values, list) or any(
            not isinstance(value, str) or not value for value in raw_values
        ):
            _fail("separation_record_invalid", field)
        partition_values = values.setdefault(partition, set())
        if partition_values.intersection(raw_values):
            _fail("duplicate_identifier_within_partition", field)
        partition_values.update(raw_values)
    return values


def _assert_pairwise_disjoint(values: Mapping[str, set[str]], field: str) -> None:
    partitions = sorted(values)
    for index, left in enumerate(partitions):
        for right in partitions[index + 1 :]:
            if values[left].intersection(values[right]):
                _fail("cross_partition_overlap", f"{field}:{left}:{right}")


def validate_partition_separation(
    assignments: Sequence[Mapping[str, Any]],
    reference_pools: Sequence[Mapping[str, Any]],
) -> SeparationAudit:
    """Vérifie les disjonctions R06 sans ouvrir de génotype."""
    required_fields = (
        "isolation_unit_ids",
        "master_identity_ids",
        "pedigree_component_ids",
        "event_ids",
        "segment_ids",
        "duplicate_group_ids",
        "seed_prefixes",
        "negative_control_ids",
    )
    sets = {field: _sets_by_partition(assignments, field) for field in required_fields}
    for field, values in sets.items():
        _assert_pairwise_disjoint(values, field)

    pool_members: dict[str, set[str]] = {}
    study_ids = set().union(*sets["master_identity_ids"].values()) if assignments else set()
    for pool in reference_pools:
        pool_id = pool.get("pool_id")
        members = pool.get("member_ids")
        if not isinstance(pool_id, str) or not isinstance(members, list) or any(
            not isinstance(member, str) or not member for member in members
        ):
            _fail("reference_pool_invalid")
        member_set = set(members)
        if len(member_set) != len(members):
            _fail("reference_member_duplicate", pool_id)
        if study_ids.intersection(member_set):
            _fail("reference_study_overlap", pool_id)
        pool_members[pool_id] = member_set
    _assert_pairwise_disjoint(pool_members, "reference_member_ids")
    if len(pool_members) != 6:
        _fail("reference_pool_count_invalid")
    return SeparationAudit(
        partition_count=len({item.get("partition_id") for item in assignments}),
        isolation_unit_count=sum(len(value) for value in sets["isolation_unit_ids"].values()),
        reference_member_count=sum(len(value) for value in pool_members.values()),
        negative_control_count=sum(len(value) for value in sets["negative_control_ids"].values()),
    )


def _flatten(value: Any, prefix: str = "") -> dict[str, Any]:
    if isinstance(value, Mapping):
        result: dict[str, Any] = {}
        for key in sorted(value):
            path = f"{prefix}.{key}" if prefix else str(key)
            result.update(_flatten(value[key], path))
        return result
    if isinstance(value, list):
        result = {}
        for index, child in enumerate(value):
            result.update(_flatten(child, f"{prefix}[{index}]"))
        return result
    return {prefix: value}


def validate_single_factor_pair(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    *,
    allowed_path: str,
    ignored_paths: Iterable[str] = ("run_id", "closed_at", "manifest_sha256"),
) -> tuple[str, Any, Any]:
    """Exige un diff canonique réduit exactement au facteur signé déclaré."""
    left_flat = _flatten(left)
    right_flat = _flatten(right)
    ignored = set(ignored_paths)
    changed = [
        path
        for path in sorted(set(left_flat) | set(right_flat))
        if path not in ignored and left_flat.get(path) != right_flat.get(path)
    ]
    if changed != [allowed_path]:
        _fail("single_factor_diff_invalid", ",".join(changed) or "no_difference")
    return allowed_path, left_flat.get(allowed_path), right_flat.get(allowed_path)


def sanitized_environment(
    environment: Mapping[str, str],
    allowed_names: Iterable[str],
) -> dict[str, str]:
    """Construit l'environnement minimal d'un processus aveugle."""
    allowed = set(allowed_names)
    result: dict[str, str] = {}
    for name, value in environment.items():
        if name not in allowed:
            continue
        folded = value.replace("\\", "/").casefold()
        if "private_truth" in folded or "/governance" in folded or "data/input" in folded:
            _fail("private_environment_value", name)
        result[name] = value
    return result
