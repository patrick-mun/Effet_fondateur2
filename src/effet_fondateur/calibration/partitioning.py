"""Partitionnement D10 des donneurs et de la référence de phasage.

Les identifiants publics ne sont jamais écrits par ce module. Les appelants
doivent conserver les listes de membres sous ``private_truth/partitions`` dans
le workspace ignoré par Git. Le manifest ne contient que des effectifs et des
empreintes SHA-256 des listes canoniques.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import csv
from dataclasses import dataclass
import hashlib
import re
from typing import Any
import unicodedata


PARTITION_METHOD_ID = "family_isolated_reference_partition_v1"
PARTITION_METHOD_VERSION = "1.0.0"
DONOR_CALIBRATION = "DONOR_CALIBRATION"
DONOR_VALIDATION = "DONOR_VALIDATION"
PHASING_REFERENCE_SHARED = "PHASING_REFERENCE_SHARED"
PARTITION_ROLES = (
    DONOR_CALIBRATION,
    DONOR_VALIDATION,
    PHASING_REFERENCE_SHARED,
)
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_EXPECTED_METHOD = {
    "method_id": PARTITION_METHOD_ID,
    "method_version": PARTITION_METHOD_VERSION,
    "identifier_normalization": "UTF8_NFC_CASEFOLD_COLLISION",
    "independence_unit": "PEDIGREE_OR_KINSHIP_CONNECTED_COMPONENT",
    "mixed_group_policy": "REMOVE_WHOLE_GROUP_FROM_PHASING_REFERENCE",
    "donor_balance_policy": "GREEDY_BALANCE_BY_GROUP_SIZE",
    "partition_order_key": "SHA256_METHOD_SOURCE_FINGERPRINT_GROUP_ID",
}
_EXPECTED_CONFIDENTIALITY = {
    "raw_identifiers_in_manifest": False,
    "membership_artifacts_private": True,
    "membership_artifacts_version_control_allowed": False,
    "aggregate_counts_publishable": True,
}


class CalibrationPartitionError(RuntimeError):
    """Signale une violation D10 avec un code stable et non sensible."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        message = code if detail is None else f"{code}: {detail}"
        super().__init__(message)


@dataclass(frozen=True)
class ReferencePartitions:
    """Listes privées produites par le partitionnement déterministe D10."""

    donor_calibration: tuple[str, ...]
    donor_validation: tuple[str, ...]
    phasing_reference_shared: tuple[str, ...]

    def by_role(self) -> dict[str, tuple[str, ...]]:
        """Retourne les membres indexés par rôle contractuel D10."""
        return {
            DONOR_CALIBRATION: self.donor_calibration,
            DONOR_VALIDATION: self.donor_validation,
            PHASING_REFERENCE_SHARED: self.phasing_reference_shared,
        }


@dataclass(frozen=True)
class PartitionAudit:
    """Résumé agrégé d'un manifest D10 validé sans exposer ses membres."""

    partition_manifest_id: str
    source_identifier_count: int
    donor_calibration_count: int
    donor_validation_count: int
    phasing_reference_count: int
    independence_group_count: int
    checks: tuple[str, ...]


def _fail(code: str, detail: str | None = None) -> None:
    raise CalibrationPartitionError(code, detail)


def _normalize_identifier(value: Any, namespace: str) -> tuple[str, str]:
    if not isinstance(value, str) or not value:
        _fail("partition_identifier_invalid", namespace)
    normalized = unicodedata.normalize("NFC", value)
    if normalized != normalized.strip() or any(
        character in normalized for character in ("\x00", "\n", "\r")
    ):
        _fail("partition_identifier_invalid", namespace)
    if any(ord(character) < 32 for character in normalized):
        _fail("partition_identifier_invalid", namespace)
    return normalized, normalized.casefold()


def _index_identifiers(
    values: Iterable[str],
    namespace: str,
) -> dict[str, str]:
    indexed: dict[str, str] = {}
    if isinstance(values, (str, bytes)):
        _fail("partition_identifier_collection_invalid", namespace)
    try:
        iterator = iter(values)
    except TypeError:
        _fail("partition_identifier_collection_invalid", namespace)
    for value in iterator:
        normalized, collision_key = _normalize_identifier(value, namespace)
        if collision_key in indexed:
            _fail("partition_identifier_duplicate", namespace)
        indexed[collision_key] = normalized
    return indexed


def _validate_sha256(value: Any, namespace: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        _fail("partition_sha256_invalid", namespace)
    return value


def _normalized_group_assignments(
    independence_group_by_identifier: Mapping[str, str],
    source_identifiers: Mapping[str, str],
) -> dict[str, str]:
    if not isinstance(independence_group_by_identifier, Mapping):
        _fail("independence_group_mapping_invalid")
    assignments: dict[str, str] = {}
    for raw_identifier, raw_group_id in independence_group_by_identifier.items():
        _, identifier_key = _normalize_identifier(raw_identifier, "group_member")
        _, group_key = _normalize_identifier(raw_group_id, "independence_group")
        if identifier_key in assignments:
            _fail("independence_group_member_duplicate")
        assignments[identifier_key] = group_key
    unknown = assignments.keys() - source_identifiers.keys()
    missing = source_identifiers.keys() - assignments.keys()
    if unknown:
        _fail("independence_group_member_unknown")
    if missing:
        _fail("independence_group_member_missing")
    return assignments


def canonical_membership_bytes(identifiers: Iterable[str]) -> bytes:
    """Encode une liste privée sous une forme canonique empreintable.

    Les identifiants sont normalisés en NFC, contrôlés contre les collisions
    insensibles à la casse, triés par leur clé de collision et terminés par un
    saut de ligne LF. Une liste vide est refusée.
    """
    indexed = _index_identifiers(identifiers, "membership")
    if not indexed:
        _fail("partition_membership_empty")
    payload = "".join(f"{indexed[key]}\n" for key in sorted(indexed))
    return payload.encode("utf-8")


def membership_sha256(identifiers: Iterable[str]) -> str:
    """Calcule le SHA-256 de la représentation canonique D10."""
    return hashlib.sha256(canonical_membership_bytes(identifiers)).hexdigest()


def parse_population_pedigree(
    lines: Iterable[str],
) -> tuple[tuple[str, ...], dict[str, str]]:
    """Extrait les sujets publics et leurs composantes familiales.

    Les groupes sont les composantes connexes créées par ``FamilyID`` et par
    les liens père/mère lorsque le parent appartient au même panel. Les noms de
    groupes retournés sont des empreintes non réversibles à longueur fixe ; ils
    ne recopient aucun identifiant source.
    """
    materialized_lines = [line for line in lines if line.strip()]
    if not materialized_lines:
        _fail("population_metadata_empty")
    fieldnames = materialized_lines[0].strip().split()
    required_columns = {
        "FamilyID",
        "SampleID",
        "FatherID",
        "MotherID",
    }
    if not required_columns <= set(fieldnames):
        _fail("population_metadata_columns_invalid")
    records: list[dict[str, str]] = []
    for line in materialized_lines[1:]:
        values = line.strip().split()
        if len(values) != len(fieldnames):
            _fail("population_metadata_row_invalid")
        records.append(dict(zip(fieldnames, values, strict=True)))
    sample_index = _index_identifiers(
        (record.get("SampleID", "") for record in records),
        "population_metadata",
    )
    if not sample_index:
        _fail("population_metadata_empty")

    parent: dict[str, str] = {key: key for key in sample_index}

    def find(key: str) -> str:
        root = key
        while parent[root] != root:
            root = parent[root]
        while parent[key] != key:
            next_key = parent[key]
            parent[key] = root
            key = next_key
        return root

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    family_anchor: dict[str, str] = {}
    for record in records:
        _, sample_key = _normalize_identifier(
            record.get("SampleID"),
            "population_metadata",
        )
        family_id = record.get("FamilyID", "")
        if isinstance(family_id, str) and family_id not in {"", "0", ".", "NA"}:
            _, family_key = _normalize_identifier(family_id, "family_id")
            anchor = family_anchor.setdefault(family_key, sample_key)
            union(sample_key, anchor)
        for parent_field in ("FatherID", "MotherID"):
            parent_id = record.get(parent_field, "")
            if not isinstance(parent_id, str) or parent_id in {"", "0", ".", "NA"}:
                continue
            _, parent_key = _normalize_identifier(parent_id, parent_field)
            if parent_key in sample_index:
                union(sample_key, parent_key)

    members_by_root: dict[str, list[str]] = {}
    for sample_key in sample_index:
        members_by_root.setdefault(find(sample_key), []).append(sample_key)
    group_by_root = {
        root: "PEDIGREE_"
        + hashlib.sha256("\n".join(sorted(members)).encode("utf-8")).hexdigest()[:24]
        for root, members in members_by_root.items()
    }
    group_by_identifier = {
        sample_index[sample_key]: group_by_root[find(sample_key)]
        for sample_key in sample_index
    }
    source_identifiers = tuple(
        sample_index[key] for key in sorted(sample_index)
    )
    return source_identifiers, group_by_identifier


def parse_unrelated_sequence_index(lines: Iterable[str]) -> tuple[str, ...]:
    """Extrait sans duplication les noms d'échantillons de l'index IGSR."""
    materialized_lines = list(lines)
    header_index = next(
        (
            index
            for index, line in enumerate(materialized_lines)
            if line.startswith("#ENA_FILE_PATH\t")
        ),
        None,
    )
    if header_index is None:
        _fail("sequence_index_header_missing")
    header = materialized_lines[header_index].lstrip("#")
    reader = csv.DictReader(
        [header, *materialized_lines[header_index + 1 :]],
        delimiter="\t",
    )
    if reader.fieldnames is None or "SAMPLE_NAME" not in reader.fieldnames:
        _fail("sequence_index_columns_invalid")
    samples: dict[str, str] = {}
    for record in reader:
        sample_name = record.get("SAMPLE_NAME", "")
        if isinstance(sample_name, str):
            sample_name = sample_name.strip()
        normalized, collision_key = _normalize_identifier(
            sample_name,
            "sequence_index",
        )
        if collision_key in samples and samples[collision_key] != normalized:
            _fail("partition_identifier_duplicate", "sequence_index")
        samples.setdefault(collision_key, normalized)
    if not samples:
        _fail("sequence_index_empty")
    return tuple(samples[key] for key in sorted(samples))


def partition_reference_identifiers(
    *,
    source_identifiers: Iterable[str],
    phasing_reference_candidate_identifiers: Iterable[str],
    independence_group_by_identifier: Mapping[str, str],
    source_fingerprint_sha256: str,
) -> ReferencePartitions:
    """Construit trois partitions disjointes par unité d'indépendance.

    Un groupe n'est conservé dans la référence de phasage que si tous ses
    membres appartiennent à la liste de candidats. Les groupes mixtes sont
    déplacés en bloc vers les donneurs, puis les groupes donneurs sont répartis
    de façon déterministe et équilibrée entre calibration et validation.
    """
    fingerprint = _validate_sha256(source_fingerprint_sha256, "source")
    source = _index_identifiers(source_identifiers, "source")
    candidates = _index_identifiers(
        phasing_reference_candidate_identifiers,
        "phasing_reference_candidates",
    )
    if len(source) < 3:
        _fail("partition_source_too_small")
    if candidates.keys() - source.keys():
        _fail("phasing_reference_candidate_unknown")
    assignments = _normalized_group_assignments(
        independence_group_by_identifier,
        source,
    )

    members_by_group: dict[str, list[str]] = {}
    keys_by_group: dict[str, set[str]] = {}
    for identifier_key, identifier in source.items():
        group_key = assignments[identifier_key]
        members_by_group.setdefault(group_key, []).append(identifier)
        keys_by_group.setdefault(group_key, set()).add(identifier_key)

    reference_groups: list[str] = []
    donor_groups: list[str] = []
    for group_key, member_keys in keys_by_group.items():
        if member_keys <= candidates.keys():
            reference_groups.append(group_key)
        else:
            donor_groups.append(group_key)
    if not reference_groups:
        _fail("phasing_reference_partition_empty")
    if len(donor_groups) < 2:
        _fail("donor_independence_groups_insufficient")

    domain = f"{PARTITION_METHOD_ID}\0{fingerprint}\0"
    donor_groups.sort(
        key=lambda group_key: hashlib.sha256(
            f"{domain}{group_key}".encode("utf-8")
        ).hexdigest()
    )
    calibration_members: list[str] = []
    validation_members: list[str] = []
    calibration_count = 0
    validation_count = 0
    for group_key in donor_groups:
        group_members = members_by_group[group_key]
        if calibration_count <= validation_count:
            calibration_members.extend(group_members)
            calibration_count += len(group_members)
        else:
            validation_members.extend(group_members)
            validation_count += len(group_members)

    if not calibration_members or not validation_members:
        _fail("donor_partition_empty")
    reference_members = [
        member
        for group_key in reference_groups
        for member in members_by_group[group_key]
    ]

    def ordered(values: Sequence[str]) -> tuple[str, ...]:
        return tuple(sorted(values, key=lambda value: value.casefold()))

    return ReferencePartitions(
        donor_calibration=ordered(calibration_members),
        donor_validation=ordered(validation_members),
        phasing_reference_shared=ordered(reference_members),
    )


def _records(value: Any, field: str) -> list[Mapping[str, Any]]:
    if not isinstance(value, list) or any(
        not isinstance(record, Mapping) for record in value
    ):
        _fail("partition_manifest_structure_invalid", field)
    return value


def _catalog_source(
    public_reference_catalog: Mapping[str, Any],
    source_id: str,
) -> Mapping[str, Any]:
    sources = _records(public_reference_catalog.get("sources"), "catalog.sources")
    matches = [source for source in sources if source.get("source_id") == source_id]
    if len(matches) != 1:
        _fail("partition_catalog_source_missing", source_id)
    return matches[0]


def _require_equal(actual: Any, expected: Any, code: str) -> None:
    if actual != expected:
        _fail(code)


def _validate_private_membership_path(value: Any) -> str:
    if (
        not isinstance(value, str)
        or "\x00" in value
        or "\\" in value
        or not value.startswith("private_truth/partitions/")
        or not value.endswith(".ids.txt")
        or any(part in {"", ".", ".."} for part in value.split("/"))
    ):
        _fail("partition_membership_path_invalid")
    return value


def validate_partition_manifest(
    manifest: Mapping[str, Any],
    *,
    memberships_by_partition_id: Mapping[str, Iterable[str]],
    source_identifiers: Iterable[str],
    independence_group_by_identifier: Mapping[str, str],
    public_reference_catalog: Mapping[str, Any],
    public_reference_catalog_sha256: str,
) -> PartitionAudit:
    """Valide le contenu privé et les engagements agrégés d'un manifest D10."""
    if not isinstance(manifest, Mapping):
        _fail("partition_manifest_structure_invalid", "manifest")
    catalog_hash = _validate_sha256(
        public_reference_catalog_sha256,
        "public_reference_catalog",
    )
    manifest_id = manifest.get("partition_manifest_id")
    if not isinstance(manifest_id, str) or not manifest_id:
        _fail("partition_manifest_identifier_invalid")
    _require_equal(
        manifest.get("schema_version"),
        "1.0.0",
        "partition_schema_version_mismatch",
    )
    _require_equal(
        manifest.get("manifest_state"),
        "FINALIZED",
        "partition_manifest_not_finalized",
    )

    catalog_contract = manifest.get("source_catalog")
    source_contract = manifest.get("source")
    method = manifest.get("method")
    confidentiality = manifest.get("confidentiality")
    integrity = manifest.get("integrity")
    if not isinstance(catalog_contract, Mapping):
        _fail("partition_manifest_structure_invalid", "source_catalog")
    if not isinstance(source_contract, Mapping):
        _fail("partition_manifest_structure_invalid", "source")
    if not isinstance(method, Mapping):
        _fail("partition_manifest_structure_invalid", "method")
    if not isinstance(confidentiality, Mapping):
        _fail("partition_manifest_structure_invalid", "confidentiality")
    if not isinstance(integrity, Mapping):
        _fail("partition_manifest_structure_invalid", "integrity")
    for field, expected in _EXPECTED_METHOD.items():
        _require_equal(
            method.get(field),
            expected,
            "partition_method_mismatch",
        )
    for field, expected in _EXPECTED_CONFIDENTIALITY.items():
        _require_equal(
            confidentiality.get(field),
            expected,
            "partition_confidentiality_violation",
        )
    _require_equal(
        catalog_contract.get("catalog_id"),
        public_reference_catalog.get("catalog_id"),
        "partition_catalog_id_mismatch",
    )
    _require_equal(
        catalog_contract.get("catalog_version"),
        public_reference_catalog.get("catalog_version"),
        "partition_catalog_version_mismatch",
    )
    _require_equal(
        catalog_contract.get("sha256"),
        catalog_hash,
        "partition_catalog_sha256_mismatch",
    )
    _require_equal(
        catalog_contract.get("repository_path"),
        "calibration/references/public_reference_catalog.json",
        "partition_catalog_path_mismatch",
    )

    source_id = source_contract.get("source_id")
    if not isinstance(source_id, str):
        _fail("partition_manifest_structure_invalid", "source.source_id")
    catalog_source = _catalog_source(public_reference_catalog, source_id)
    catalog_assembly = catalog_source.get("assembly")
    catalog_fingerprint = catalog_source.get("fingerprint")
    if not isinstance(catalog_assembly, Mapping) or not isinstance(
        catalog_fingerprint, Mapping
    ):
        _fail("partition_catalog_source_invalid", source_id)
    for manifest_field, catalog_value in (
        ("source_type", catalog_source.get("source_type")),
        ("version", catalog_source.get("version")),
        ("assembly", catalog_assembly.get("name")),
        ("fingerprint_sha256", catalog_fingerprint.get("sha256")),
    ):
        _require_equal(
            source_contract.get(manifest_field),
            catalog_value,
            "partition_source_catalog_mismatch",
        )

    partition_inputs = _records(
        manifest.get("partition_inputs"),
        "partition_inputs",
    )
    expected_input_ids = {
        "INDEPENDENCE_GROUPS": "igsr_1kg_3202_population_metadata_20130606",
        "PHASING_REFERENCE_CANDIDATES": "igsr_1kg_2504_unrelated_sequence_index",
    }
    if len(partition_inputs) != len(expected_input_ids):
        _fail("partition_input_cardinality_invalid")
    seen_purposes: set[str] = set()
    for partition_input in partition_inputs:
        purpose = partition_input.get("purpose")
        source_input_id = partition_input.get("source_id")
        if purpose not in expected_input_ids or purpose in seen_purposes:
            _fail("partition_input_cardinality_invalid")
        _require_equal(
            source_input_id,
            expected_input_ids[str(purpose)],
            "partition_input_source_mismatch",
        )
        if not isinstance(source_input_id, str):
            _fail("partition_input_source_mismatch")
        catalog_input = _catalog_source(
            public_reference_catalog,
            source_input_id,
        )
        catalog_input_fingerprint = catalog_input.get("fingerprint")
        if not isinstance(catalog_input_fingerprint, Mapping):
            _fail("partition_catalog_source_invalid", source_input_id)
        _require_equal(
            partition_input.get("version"),
            catalog_input.get("version"),
            "partition_input_catalog_mismatch",
        )
        _require_equal(
            partition_input.get("fingerprint_sha256"),
            catalog_input_fingerprint.get("sha256"),
            "partition_input_catalog_mismatch",
        )
        seen_purposes.add(str(purpose))

    source = _index_identifiers(source_identifiers, "source")
    assignments = _normalized_group_assignments(
        independence_group_by_identifier,
        source,
    )
    partitions = _records(manifest.get("partitions"), "partitions")
    if len(partitions) != 3:
        _fail("partition_count_invalid")
    partitions_by_id: dict[str, Mapping[str, Any]] = {}
    roles: set[str] = set()
    normalized_paths: set[str] = set()
    for partition in partitions:
        partition_id = partition.get("partition_id")
        role = partition.get("role")
        path = _validate_private_membership_path(
            partition.get("membership_artifact_path")
        )
        if not isinstance(partition_id, str) or not partition_id:
            _fail("partition_manifest_identifier_invalid")
        partition_key = partition_id.casefold()
        if partition_key in {key.casefold() for key in partitions_by_id}:
            _fail("partition_manifest_identifier_duplicate")
        if role not in PARTITION_ROLES or role in roles:
            _fail("partition_role_cardinality_invalid")
        if path.casefold() in normalized_paths:
            _fail("partition_membership_path_duplicate")
        _require_equal(
            partition.get("shared_between_panels"),
            role == PHASING_REFERENCE_SHARED,
            "partition_sharing_policy_violation",
        )
        partitions_by_id[partition_id] = partition
        roles.add(str(role))
        normalized_paths.add(path.casefold())
    if roles != set(PARTITION_ROLES):
        _fail("partition_role_cardinality_invalid")
    if not isinstance(memberships_by_partition_id, Mapping) or set(
        memberships_by_partition_id
    ) != set(partitions_by_id):
        _fail("partition_membership_registry_mismatch")

    group_sets: dict[str, set[str]] = {}
    role_counts: dict[str, int] = {}
    union: set[str] = set()
    for partition_id, partition in partitions_by_id.items():
        membership = _index_identifiers(
            memberships_by_partition_id[partition_id],
            partition_id,
        )
        if not membership:
            _fail("partition_membership_empty", partition_id)
        member_keys = set(membership)
        if union & member_keys:
            _fail("forbidden_partition_identifier_overlap")
        if member_keys - source.keys():
            _fail("partition_member_unknown")
        union.update(member_keys)
        partition_groups = {assignments[key] for key in member_keys}
        group_sets[partition_id] = partition_groups
        role = str(partition["role"])
        role_counts[role] = len(member_keys)
        _require_equal(
            partition.get("membership_sha256"),
            membership_sha256(membership.values()),
            "partition_membership_sha256_mismatch",
        )
        _require_equal(
            partition.get("member_count"),
            len(member_keys),
            "partition_member_count_mismatch",
        )
        _require_equal(
            partition.get("independence_group_count"),
            len(partition_groups),
            "partition_group_count_mismatch",
        )

    unknown = union - source.keys()
    unassigned = source.keys() - union
    if unknown:
        _fail("partition_member_unknown")
    if unassigned:
        _fail("partition_member_unassigned")
    partition_ids = tuple(partitions_by_id)
    for index, left_id in enumerate(partition_ids):
        for right_id in partition_ids[index + 1 :]:
            if group_sets[left_id] & group_sets[right_id]:
                _fail("forbidden_partition_independence_group_overlap")

    source_count = len(source)
    donor_difference = abs(
        role_counts[DONOR_CALIBRATION] - role_counts[DONOR_VALIDATION]
    )
    expected_integrity = {
        "partition_count": 3,
        "assigned_identifier_count": len(union),
        "unique_identifier_count": len(union),
        "source_identifier_count": source_count,
        "unassigned_identifier_count": 0,
        "unknown_identifier_count": 0,
        "forbidden_identifier_overlap_count": 0,
        "forbidden_independence_group_overlap_count": 0,
        "all_source_identifiers_assigned_once": True,
        "donor_member_count_difference": donor_difference,
    }
    for field, expected in expected_integrity.items():
        _require_equal(
            integrity.get(field),
            expected,
            "partition_integrity_mismatch",
        )
    _require_equal(
        source_contract.get("source_identifier_count"),
        source_count,
        "partition_source_count_mismatch",
    )

    return PartitionAudit(
        partition_manifest_id=manifest_id,
        source_identifier_count=source_count,
        donor_calibration_count=role_counts[DONOR_CALIBRATION],
        donor_validation_count=role_counts[DONOR_VALIDATION],
        phasing_reference_count=role_counts[PHASING_REFERENCE_SHARED],
        independence_group_count=len(set(assignments.values())),
        checks=(
            "SOURCE_CATALOG_MATCHED",
            "ALL_IDENTIFIERS_ASSIGNED_ONCE",
            "IDENTIFIER_PARTITIONS_DISJOINT",
            "INDEPENDENCE_GROUPS_DISJOINT",
            "MEMBERSHIP_SHA256_MATCHED",
            "RAW_IDENTIFIERS_NOT_IN_MANIFEST",
        ),
    )


def build_partition_manifest(
    *,
    partition_manifest_id: str,
    created_at: str,
    partitions: ReferencePartitions,
    source_identifiers: Iterable[str],
    independence_group_by_identifier: Mapping[str, str],
    public_reference_catalog: Mapping[str, Any],
    public_reference_catalog_sha256: str,
) -> tuple[dict[str, Any], dict[str, tuple[str, ...]]]:
    """Construit puis audite un manifest D10 sans identifiant brut intégré."""
    source = _index_identifiers(source_identifiers, "source")
    assignments = _normalized_group_assignments(
        independence_group_by_identifier,
        source,
    )
    panel = _catalog_source(
        public_reference_catalog,
        "igsr_1kg_3202_phased_grch38_20220422",
    )
    metadata = _catalog_source(
        public_reference_catalog,
        "igsr_1kg_3202_population_metadata_20130606",
    )
    unrelated_index = _catalog_source(
        public_reference_catalog,
        "igsr_1kg_2504_unrelated_sequence_index",
    )

    def fingerprint(record: Mapping[str, Any]) -> str:
        value = record.get("fingerprint")
        if not isinstance(value, Mapping):
            _fail("partition_catalog_source_invalid")
        return _validate_sha256(value.get("sha256"), "catalog_source")

    role_to_partition_id = {
        DONOR_CALIBRATION: "donor_calibration",
        DONOR_VALIDATION: "donor_validation",
        PHASING_REFERENCE_SHARED: "phasing_reference",
    }
    memberships = {
        role_to_partition_id[role]: members
        for role, members in partitions.by_role().items()
    }
    partition_records: list[dict[str, Any]] = []
    role_counts: dict[str, int] = {}
    for role, members in partitions.by_role().items():
        partition_id = role_to_partition_id[role]
        normalized_members = _index_identifiers(members, partition_id)
        group_count = len(
            {assignments[member_key] for member_key in normalized_members}
        )
        role_counts[role] = len(normalized_members)
        partition_records.append(
            {
                "partition_id": partition_id,
                "role": role,
                "membership_artifact_path": (
                    f"private_truth/partitions/{partition_id}.ids.txt"
                ),
                "membership_sha256": membership_sha256(members),
                "member_count": len(normalized_members),
                "independence_group_count": group_count,
                "shared_between_panels": role == PHASING_REFERENCE_SHARED,
            }
        )

    panel_assembly = panel.get("assembly")
    if not isinstance(panel_assembly, Mapping):
        _fail("partition_catalog_source_invalid")
    source_count = len(source)
    manifest: dict[str, Any] = {
        "schema_version": "1.0.0",
        "partition_manifest_id": partition_manifest_id,
        "manifest_state": "FINALIZED",
        "created_at": created_at,
        "source_catalog": {
            "catalog_id": public_reference_catalog.get("catalog_id"),
            "catalog_version": public_reference_catalog.get("catalog_version"),
            "repository_path": "calibration/references/public_reference_catalog.json",
            "sha256": public_reference_catalog_sha256,
        },
        "source": {
            "source_id": panel.get("source_id"),
            "source_type": panel.get("source_type"),
            "version": panel.get("version"),
            "assembly": panel_assembly.get("name"),
            "fingerprint_sha256": fingerprint(panel),
            "source_identifier_count": source_count,
        },
        "partition_inputs": [
            {
                "purpose": "INDEPENDENCE_GROUPS",
                "source_id": metadata.get("source_id"),
                "version": metadata.get("version"),
                "fingerprint_sha256": fingerprint(metadata),
            },
            {
                "purpose": "PHASING_REFERENCE_CANDIDATES",
                "source_id": unrelated_index.get("source_id"),
                "version": unrelated_index.get("version"),
                "fingerprint_sha256": fingerprint(unrelated_index),
            },
        ],
        "method": dict(_EXPECTED_METHOD),
        "confidentiality": dict(_EXPECTED_CONFIDENTIALITY),
        "partitions": partition_records,
        "integrity": {
            "partition_count": 3,
            "assigned_identifier_count": source_count,
            "unique_identifier_count": source_count,
            "source_identifier_count": source_count,
            "unassigned_identifier_count": 0,
            "unknown_identifier_count": 0,
            "forbidden_identifier_overlap_count": 0,
            "forbidden_independence_group_overlap_count": 0,
            "all_source_identifiers_assigned_once": True,
            "donor_member_count_difference": abs(
                role_counts[DONOR_CALIBRATION]
                - role_counts[DONOR_VALIDATION]
            ),
        },
    }
    validate_partition_manifest(
        manifest,
        memberships_by_partition_id=memberships,
        source_identifiers=source.values(),
        independence_group_by_identifier=independence_group_by_identifier,
        public_reference_catalog=public_reference_catalog,
        public_reference_catalog_sha256=public_reference_catalog_sha256,
    )
    return manifest, memberships
