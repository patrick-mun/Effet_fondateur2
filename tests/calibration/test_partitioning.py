"""Tests du partitionnement familial D10 et de son audit privé."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from effet_fondateur.calibration.partitioning import (
    CalibrationPartitionError,
    DONOR_CALIBRATION,
    DONOR_VALIDATION,
    PHASING_REFERENCE_SHARED,
    ReferencePartitions,
    membership_sha256,
    parse_population_pedigree,
    parse_unrelated_sequence_index,
    partition_reference_identifiers,
    validate_partition_manifest,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PUBLIC_CATALOG_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "references"
    / "public_reference_catalog.json"
)


def _public_catalog() -> dict[str, Any]:
    return json.loads(PUBLIC_CATALOG_PATH.read_text(encoding="utf-8"))


def _source_identifiers() -> list[str]:
    return [f"SYN_SOURCE_{index:03d}" for index in range(1, 13)]


def _phasing_candidates() -> list[str]:
    return [f"SYN_SOURCE_{index:03d}" for index in range(1, 8)]


def _independence_groups() -> dict[str, str]:
    groups = {
        identifier: f"GROUP_{index:03d}"
        for index, identifier in enumerate(_source_identifiers(), start=1)
    }
    groups["SYN_SOURCE_001"] = "MIXED_FAMILY_001"
    groups["SYN_SOURCE_008"] = "MIXED_FAMILY_001"
    return groups


def _source_fingerprint() -> str:
    catalog = _public_catalog()
    panel = next(
        source
        for source in catalog["sources"]
        if source["source_type"] == "PHASED_HAPLOTYPE_PANEL"
    )
    return str(panel["fingerprint"]["sha256"])


def _partitions() -> ReferencePartitions:
    return partition_reference_identifiers(
        source_identifiers=_source_identifiers(),
        phasing_reference_candidate_identifiers=_phasing_candidates(),
        independence_group_by_identifier=_independence_groups(),
        source_fingerprint_sha256=_source_fingerprint(),
    )


def _partition_id_by_role() -> dict[str, str]:
    return {
        DONOR_CALIBRATION: "donor_calibration",
        DONOR_VALIDATION: "donor_validation",
        PHASING_REFERENCE_SHARED: "phasing_reference",
    }


def _memberships_by_partition_id(
    partitions: ReferencePartitions | None = None,
) -> dict[str, list[str]]:
    selected = _partitions() if partitions is None else partitions
    ids = _partition_id_by_role()
    return {
        ids[role]: list(members)
        for role, members in selected.by_role().items()
    }


def _group_count(members: list[str], groups: dict[str, str]) -> int:
    return len({groups[member].casefold() for member in members})


def _valid_manifest() -> tuple[dict[str, Any], dict[str, list[str]]]:
    catalog = _public_catalog()
    catalog_hash = hashlib.sha256(PUBLIC_CATALOG_PATH.read_bytes()).hexdigest()
    panel = next(
        source
        for source in catalog["sources"]
        if source["source_type"] == "PHASED_HAPLOTYPE_PANEL"
    )
    memberships = _memberships_by_partition_id()
    groups = _independence_groups()
    role_by_id = {
        partition_id: role
        for role, partition_id in _partition_id_by_role().items()
    }
    partitions = []
    for partition_id, members in memberships.items():
        role = role_by_id[partition_id]
        partitions.append(
            {
                "partition_id": partition_id,
                "role": role,
                "membership_artifact_path": (
                    f"private_truth/partitions/{partition_id}.ids.txt"
                ),
                "membership_sha256": membership_sha256(members),
                "member_count": len(members),
                "independence_group_count": _group_count(members, groups),
                "shared_between_panels": role == PHASING_REFERENCE_SHARED,
            }
        )
    source_count = len(_source_identifiers())
    donor_difference = abs(
        len(memberships["donor_calibration"])
        - len(memberships["donor_validation"])
    )
    manifest = {
        "schema_version": "1.0.0",
        "partition_manifest_id": "synthetic_partition_manifest_001",
        "manifest_state": "FINALIZED",
        "created_at": "2026-10-28T12:00:00Z",
        "source_catalog": {
            "catalog_id": catalog["catalog_id"],
            "catalog_version": catalog["catalog_version"],
            "repository_path": "calibration/references/public_reference_catalog.json",
            "sha256": catalog_hash,
        },
        "source": {
            "source_id": panel["source_id"],
            "source_type": panel["source_type"],
            "version": panel["version"],
            "assembly": panel["assembly"]["name"],
            "fingerprint_sha256": panel["fingerprint"]["sha256"],
            "source_identifier_count": source_count,
        },
        "partition_inputs": [
            {
                "purpose": "INDEPENDENCE_GROUPS",
                "source_id": "igsr_1kg_3202_population_metadata_20130606",
                "version": "20130606_g1k_3202_samples_ped_population",
                "fingerprint_sha256": "4e164b717433fc9bfc16e29b20e5b4659b4216185b3275322d421fe3d323c132",
            },
            {
                "purpose": "PHASING_REFERENCE_CANDIDATES",
                "source_id": "igsr_1kg_2504_unrelated_sequence_index",
                "version": "1000G_2504_high_coverage.sequence.index",
                "fingerprint_sha256": "0ac5fdda0b60b6575b1867827f53e39d723f915b325b10bf023b2dc3d71c08f2",
            },
        ],
        "method": {
            "method_id": "family_isolated_reference_partition_v1",
            "method_version": "1.0.0",
            "identifier_normalization": "UTF8_NFC_CASEFOLD_COLLISION",
            "independence_unit": "PEDIGREE_OR_KINSHIP_CONNECTED_COMPONENT",
            "mixed_group_policy": "REMOVE_WHOLE_GROUP_FROM_PHASING_REFERENCE",
            "donor_balance_policy": "GREEDY_BALANCE_BY_GROUP_SIZE",
            "partition_order_key": "SHA256_METHOD_SOURCE_FINGERPRINT_GROUP_ID",
        },
        "confidentiality": {
            "raw_identifiers_in_manifest": False,
            "membership_artifacts_private": True,
            "membership_artifacts_version_control_allowed": False,
            "aggregate_counts_publishable": True,
        },
        "partitions": partitions,
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
            "donor_member_count_difference": donor_difference,
        },
    }
    return manifest, memberships


def _validate(
    manifest: dict[str, Any],
    memberships: dict[str, list[str]],
    *,
    source_identifiers: list[str] | None = None,
    groups: dict[str, str] | None = None,
    catalog: dict[str, Any] | None = None,
    catalog_sha256: str | None = None,
):
    return validate_partition_manifest(
        manifest,
        memberships_by_partition_id=memberships,
        source_identifiers=(
            _source_identifiers() if source_identifiers is None else source_identifiers
        ),
        independence_group_by_identifier=(
            _independence_groups() if groups is None else groups
        ),
        public_reference_catalog=_public_catalog() if catalog is None else catalog,
        public_reference_catalog_sha256=(
            hashlib.sha256(PUBLIC_CATALOG_PATH.read_bytes()).hexdigest()
            if catalog_sha256 is None
            else catalog_sha256
        ),
    )


def _assert_error(expected_code: str, action: Any) -> None:
    with pytest.raises(CalibrationPartitionError) as captured:
        action()
    assert captured.value.code == expected_code


def test_partitioning_is_deterministic_and_assigns_every_source_once() -> None:
    first = _partitions()
    second = partition_reference_identifiers(
        source_identifiers=reversed(_source_identifiers()),
        phasing_reference_candidate_identifiers=reversed(_phasing_candidates()),
        independence_group_by_identifier=dict(
            reversed(list(_independence_groups().items()))
        ),
        source_fingerprint_sha256=_source_fingerprint(),
    )

    assert first == second
    member_sets = [set(members) for members in first.by_role().values()]
    assert set.union(*member_sets) == set(_source_identifiers())
    assert sum(len(members) for members in member_sets) == len(
        set.union(*member_sets)
    )


def test_mixed_independence_group_is_removed_from_phasing_reference() -> None:
    partitions = _partitions()
    reference = set(partitions.phasing_reference_shared)
    donor_calibration = set(partitions.donor_calibration)
    donor_validation = set(partitions.donor_validation)

    assert "SYN_SOURCE_001" not in reference
    assert "SYN_SOURCE_008" not in reference
    assert {"SYN_SOURCE_001", "SYN_SOURCE_008"} <= (
        donor_calibration | donor_validation
    )
    assert not (
        {"SYN_SOURCE_001", "SYN_SOURCE_008"} & donor_calibration
        and {"SYN_SOURCE_001", "SYN_SOURCE_008"} & donor_validation
    )


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        ("candidate_unknown", "phasing_reference_candidate_unknown"),
        ("group_missing", "independence_group_member_missing"),
        ("source_duplicate", "partition_identifier_duplicate"),
        ("invalid_hash", "partition_sha256_invalid"),
        ("donor_groups_insufficient", "donor_independence_groups_insufficient"),
        ("source_string", "partition_identifier_collection_invalid"),
    ],
)
def test_partitioning_rejects_invalid_inputs(
    mutation: str,
    expected_code: str,
) -> None:
    source: Any = _source_identifiers()
    candidates = _phasing_candidates()
    groups = _independence_groups()
    fingerprint = _source_fingerprint()
    if mutation == "candidate_unknown":
        candidates.append("UNKNOWN_SOURCE")
    elif mutation == "group_missing":
        del groups[source[0]]
    elif mutation == "source_duplicate":
        source.append(source[0].lower())
    elif mutation == "invalid_hash":
        fingerprint = "not-a-sha256"
    elif mutation == "donor_groups_insufficient":
        candidates = source[:-1]
    elif mutation == "source_string":
        source = "SYN_SOURCE_001"

    _assert_error(
        expected_code,
        lambda: partition_reference_identifiers(
            source_identifiers=source,
            phasing_reference_candidate_identifiers=candidates,
            independence_group_by_identifier=groups,
            source_fingerprint_sha256=fingerprint,
        ),
    )


def test_valid_partition_manifest_passes_aggregate_audit() -> None:
    manifest, memberships = _valid_manifest()

    audit = _validate(manifest, memberships)

    assert audit.source_identifier_count == 12
    assert audit.donor_calibration_count + audit.donor_validation_count == 6
    assert audit.phasing_reference_count == 6
    assert audit.independence_group_count == 11
    assert "INDEPENDENCE_GROUPS_DISJOINT" in audit.checks


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        ("identifier_overlap", "forbidden_partition_identifier_overlap"),
        ("membership_hash", "partition_membership_sha256_mismatch"),
        ("member_count", "partition_member_count_mismatch"),
        ("catalog_hash", "partition_catalog_sha256_mismatch"),
        ("source_version", "partition_source_catalog_mismatch"),
        ("partition_input", "partition_input_catalog_mismatch"),
        ("method", "partition_method_mismatch"),
        ("confidentiality", "partition_confidentiality_violation"),
        ("private_path", "partition_membership_path_invalid"),
        ("sharing", "partition_sharing_policy_violation"),
        ("unassigned", "partition_member_unassigned"),
        ("unknown", "partition_member_unknown"),
    ],
)
def test_partition_manifest_audit_rejects_inconsistencies(
    mutation: str,
    expected_code: str,
) -> None:
    manifest, memberships = _valid_manifest()
    catalog_hash: str | None = None
    if mutation == "identifier_overlap":
        memberships["donor_validation"].append(memberships["donor_calibration"][0])
    elif mutation == "membership_hash":
        manifest["partitions"][0]["membership_sha256"] = "0" * 64
    elif mutation == "member_count":
        manifest["partitions"][0]["member_count"] += 1
    elif mutation == "catalog_hash":
        catalog_hash = "0" * 64
    elif mutation == "source_version":
        manifest["source"]["version"] = "unreviewed-version"
    elif mutation == "partition_input":
        manifest["partition_inputs"][0]["fingerprint_sha256"] = "0" * 64
    elif mutation == "method":
        manifest["method"]["mixed_group_policy"] = "ALLOW_GROUP_SPLIT"
    elif mutation == "confidentiality":
        manifest["confidentiality"]["raw_identifiers_in_manifest"] = True
    elif mutation == "private_path":
        manifest["partitions"][0]["membership_artifact_path"] = (
            "blinded_inputs/donor.ids.txt"
        )
    elif mutation == "sharing":
        manifest["partitions"][0]["shared_between_panels"] = True
    elif mutation == "unassigned":
        removed = memberships["phasing_reference"].pop()
        partition = next(
            item
            for item in manifest["partitions"]
            if item["partition_id"] == "phasing_reference"
        )
        partition["membership_sha256"] = membership_sha256(
            memberships["phasing_reference"]
        )
        partition["member_count"] -= 1
        partition["independence_group_count"] -= 1
        assert removed in _source_identifiers()
    elif mutation == "unknown":
        memberships["phasing_reference"].append("UNKNOWN_SOURCE")
        partition = next(
            item
            for item in manifest["partitions"]
            if item["partition_id"] == "phasing_reference"
        )
        partition["membership_sha256"] = membership_sha256(
            memberships["phasing_reference"]
        )
        partition["member_count"] += 1
        partition["independence_group_count"] += 1

    _assert_error(
        expected_code,
        lambda: _validate(
            manifest,
            memberships,
            catalog_sha256=catalog_hash,
        ),
    )


def test_partition_manifest_rejects_independence_group_crossing() -> None:
    manifest, memberships = _valid_manifest()
    groups = _independence_groups()
    calibration_member = memberships["donor_calibration"][0]
    validation_member = memberships["donor_validation"][0]
    groups[validation_member] = groups[calibration_member]
    validation_partition = next(
        item
        for item in manifest["partitions"]
        if item["partition_id"] == "donor_validation"
    )
    validation_partition["independence_group_count"] = _group_count(
        memberships["donor_validation"],
        groups,
    )

    _assert_error(
        "forbidden_partition_independence_group_overlap",
        lambda: _validate(manifest, memberships, groups=groups),
    )


def test_membership_hash_is_order_independent_but_case_preserving() -> None:
    assert membership_sha256(["ÉCHANTILLON_2", "ECHANTILLON_1"]) == (
        membership_sha256(["ECHANTILLON_1", "ÉCHANTILLON_2"])
    )
    assert membership_sha256(["SYN_A"]) != membership_sha256(["syn_a"])


def test_population_metadata_builds_connected_family_components() -> None:
    lines = [
        "FamilyID\tSampleID\tFatherID\tMotherID\tSex\tPopulation\tSuperpopulation\n",
        "FAM_A\tSYN_PARENT\t0\t0\t1\tP1\tSP1\n",
        "FAM_A\tSYN_CHILD\tSYN_PARENT\t0\t2\tP1\tSP1\n",
        "FAM_B\tSYN_SINGLETON\t0\t0\t1\tP2\tSP2\n",
    ]

    identifiers, groups = parse_population_pedigree(lines)

    assert identifiers == ("SYN_CHILD", "SYN_PARENT", "SYN_SINGLETON")
    assert groups["SYN_CHILD"] == groups["SYN_PARENT"]
    assert groups["SYN_SINGLETON"] != groups["SYN_CHILD"]
    assert all(group.startswith("PEDIGREE_") for group in groups.values())
    assert all("SYN_" not in group for group in groups.values())


def test_sequence_index_extracts_unique_sample_names() -> None:
    lines = [
        "##FileDate=20260101\n",
        "#ENA_FILE_PATH\tSAMPLE_NAME\tPOPULATION\n",
        "path/a\tSYN_A\tP1\n",
        "path/b\tSYN_A\tP1\n",
        "path/c\t SYN_B \tP2\n",
    ]

    assert parse_unrelated_sequence_index(lines) == ("SYN_A", "SYN_B")


@pytest.mark.parametrize(
    ("parser", "lines", "expected_code"),
    [
        (
            parse_population_pedigree,
            ["SampleID\tPopulation\n", "SYN_A\tP1\n"],
            "population_metadata_columns_invalid",
        ),
        (
            parse_unrelated_sequence_index,
            ["##FileDate=20260101\n"],
            "sequence_index_header_missing",
        ),
    ],
)
def test_public_metadata_parsers_reject_incomplete_contracts(
    parser: Any,
    lines: list[str],
    expected_code: str,
) -> None:
    _assert_error(expected_code, lambda: parser(lines))
