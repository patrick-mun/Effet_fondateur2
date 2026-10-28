"""Tests du schéma D10 du manifest privé de partitions."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest

from effet_fondateur.calibration.partitioning import membership_sha256


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "schemas"
    / "reference_partition_manifest.schema.json"
)
PUBLIC_CATALOG_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "references"
    / "public_reference_catalog.json"
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Refuse les clés JSON dupliquées dans les contrats D10."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Clé JSON dupliquée dans le contrat D10 : {key}")
        result[key] = value
    return result


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_keys,
    )


def _validator() -> jsonschema.Draft202012Validator:
    schema = _load_json(SCHEMA_PATH)
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(
        schema,
        format_checker=jsonschema.FormatChecker(),
    )


def _valid_manifest() -> dict[str, Any]:
    catalog = _load_json(PUBLIC_CATALOG_PATH)
    panel = next(
        source
        for source in catalog["sources"]
        if source["source_type"] == "PHASED_HAPLOTYPE_PANEL"
    )
    memberships = {
        "donor_calibration": ["SYN_CAL_001", "SYN_CAL_002"],
        "donor_validation": ["SYN_VAL_001", "SYN_VAL_002"],
        "phasing_reference": ["SYN_REF_001", "SYN_REF_002", "SYN_REF_003"],
    }
    partitions = [
        {
            "partition_id": "donor_calibration",
            "role": "DONOR_CALIBRATION",
            "membership_artifact_path": "private_truth/partitions/donor_calibration.ids.txt",
            "membership_sha256": membership_sha256(
                memberships["donor_calibration"]
            ),
            "member_count": 2,
            "independence_group_count": 2,
            "shared_between_panels": False,
        },
        {
            "partition_id": "donor_validation",
            "role": "DONOR_VALIDATION",
            "membership_artifact_path": "private_truth/partitions/donor_validation.ids.txt",
            "membership_sha256": membership_sha256(
                memberships["donor_validation"]
            ),
            "member_count": 2,
            "independence_group_count": 2,
            "shared_between_panels": False,
        },
        {
            "partition_id": "phasing_reference",
            "role": "PHASING_REFERENCE_SHARED",
            "membership_artifact_path": "private_truth/partitions/phasing_reference.ids.txt",
            "membership_sha256": membership_sha256(
                memberships["phasing_reference"]
            ),
            "member_count": 3,
            "independence_group_count": 3,
            "shared_between_panels": True,
        },
    ]
    return {
        "schema_version": "1.0.0",
        "partition_manifest_id": "synthetic_partition_manifest_001",
        "manifest_state": "FINALIZED",
        "created_at": "2026-10-28T12:00:00Z",
        "source_catalog": {
            "catalog_id": catalog["catalog_id"],
            "catalog_version": catalog["catalog_version"],
            "repository_path": "calibration/references/public_reference_catalog.json",
            "sha256": hashlib.sha256(PUBLIC_CATALOG_PATH.read_bytes()).hexdigest(),
        },
        "source": {
            "source_id": panel["source_id"],
            "source_type": panel["source_type"],
            "version": panel["version"],
            "assembly": panel["assembly"]["name"],
            "fingerprint_sha256": panel["fingerprint"]["sha256"],
            "source_identifier_count": 7,
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
            "assigned_identifier_count": 7,
            "unique_identifier_count": 7,
            "source_identifier_count": 7,
            "unassigned_identifier_count": 0,
            "unknown_identifier_count": 0,
            "forbidden_identifier_overlap_count": 0,
            "forbidden_independence_group_overlap_count": 0,
            "all_source_identifiers_assigned_once": True,
            "donor_member_count_difference": 0,
        },
    }


def _errors(manifest: dict[str, Any]) -> list[jsonschema.ValidationError]:
    return sorted(
        _validator().iter_errors(manifest),
        key=lambda error: list(error.absolute_path),
    )


def test_schema_and_synthetic_manifest_are_valid() -> None:
    assert _errors(_valid_manifest()) == []


@pytest.mark.parametrize(
    "missing_field",
    [
        "source_catalog",
        "source",
        "partition_inputs",
        "method",
        "confidentiality",
        "partitions",
        "integrity",
    ],
)
def test_each_partition_contract_section_is_required(missing_field: str) -> None:
    manifest = _valid_manifest()
    del manifest[missing_field]

    assert _errors(manifest)


@pytest.mark.parametrize(
    "unsafe_path",
    [
        "partitions/donor.ids.txt",
        "private_truth/partitions/../donor.ids.txt",
        "/private_truth/partitions/donor.ids.txt",
        "data/input/dock6_real/donor.ids.txt",
        "private_truth/partitions/donor.txt",
    ],
)
def test_membership_artifacts_must_remain_private_and_canonical(
    unsafe_path: str,
) -> None:
    manifest = _valid_manifest()
    manifest["partitions"][0]["membership_artifact_path"] = unsafe_path

    assert _errors(manifest)


def test_raw_identifiers_cannot_be_added_to_manifest() -> None:
    manifest = _valid_manifest()
    manifest["partitions"][0]["identifiers"] = ["PUBLIC_SAMPLE_001"]

    assert _errors(manifest)


@pytest.mark.parametrize(
    ("role", "shared_between_panels"),
    [
        ("DONOR_CALIBRATION", True),
        ("DONOR_VALIDATION", True),
        ("PHASING_REFERENCE_SHARED", False),
    ],
)
def test_only_the_phasing_reference_may_be_shared(
    role: str,
    shared_between_panels: bool,
) -> None:
    manifest = _valid_manifest()
    partition = next(
        partition for partition in manifest["partitions"] if partition["role"] == role
    )
    partition["shared_between_panels"] = shared_between_panels

    assert _errors(manifest)


@pytest.mark.parametrize(
    "integrity_field",
    [
        "unassigned_identifier_count",
        "unknown_identifier_count",
        "forbidden_identifier_overlap_count",
        "forbidden_independence_group_overlap_count",
    ],
)
def test_any_declared_partition_violation_is_rejected(
    integrity_field: str,
) -> None:
    manifest = _valid_manifest()
    manifest["integrity"][integrity_field] = 1

    assert _errors(manifest)


def test_duplicate_partition_role_is_rejected() -> None:
    manifest = _valid_manifest()
    manifest["partitions"][1]["role"] = "DONOR_CALIBRATION"

    assert _errors(manifest)


def test_unknown_reference_source_is_rejected() -> None:
    manifest = _valid_manifest()
    manifest["source"]["source_id"] = "unreviewed_public_panel"

    assert _errors(manifest)


def test_partition_input_purposes_cannot_be_swapped() -> None:
    manifest = _valid_manifest()
    manifest["partition_inputs"][0]["source_id"] = (
        "igsr_1kg_2504_unrelated_sequence_index"
    )

    assert _errors(manifest)


def test_duplicate_json_keys_are_rejected() -> None:
    with pytest.raises(ValueError, match="Clé JSON dupliquée"):
        json.loads(
            '{"schema_version": "1.0.0", "schema_version": "2.0.0"}',
            object_pairs_hook=_reject_duplicate_keys,
        )
