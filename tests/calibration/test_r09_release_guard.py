"""Chemins, aveuglement, pseudonymes et disjonctions R06/R09."""

from copy import deepcopy
from hashlib import sha256
from pathlib import Path

import pytest

from effet_fondateur.calibration.release_guard import (
    ReleaseGuardError,
    partition_pseudonym,
    sanitized_environment,
    scan_blinded_document,
    validate_blinded_inventory,
    validate_partition_separation,
    validate_single_factor_pair,
)


def _inventory(content: bytes = b"synthetic\n") -> dict[str, object]:
    return {
        "inventory_id": "INV-RESOURCE-001",
        "artifacts": [
            {
                "artifact_id": "synthetic_input",
                "relative_path": "inputs/panel.vcf.gz",
                "sha256": sha256(content).hexdigest(),
                "size_bytes": len(content),
                "media_type": "application/gzip",
            }
        ],
    }


def test_static_blinded_inventory_is_a_closed_relative_whitelist(tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    view = tmp_path / "view"
    repository.mkdir()
    view.mkdir()
    audit = validate_blinded_inventory(
        _inventory(), view_root=view, repository_root=repository, verify_files=False
    )
    assert audit.artifact_count == 1
    assert audit.verified_files is False


def test_physical_inventory_checks_size_and_hash(tmp_path: Path) -> None:
    content = b"synthetic\n"
    repository = tmp_path / "repository"
    view = tmp_path / "view"
    (view / "inputs").mkdir(parents=True)
    repository.mkdir()
    (view / "inputs/panel.vcf.gz").write_bytes(content)
    audit = validate_blinded_inventory(
        _inventory(content), view_root=view, repository_root=repository, verify_files=True
    )
    assert audit.declared_bytes == len(content)
    assert audit.verified_files is True


@pytest.mark.parametrize(
    ("mutation", "code"),
    [
        (("relative_path", "../private_truth/truth.tsv"), "unsafe_relative_path"),
        (("relative_path", "data/input/real.vcf"), "real_input_path_forbidden"),
        (("private_path", "inputs/panel.vcf.gz"), "truth_field_in_blinded_document"),
    ],
)
def test_blinded_inventory_rejects_truth_and_real_paths(
    tmp_path: Path,
    mutation: tuple[str, str],
    code: str,
) -> None:
    repository = tmp_path / "repository"
    view = tmp_path / "view"
    repository.mkdir()
    view.mkdir()
    inventory = _inventory()
    if mutation[0] == "private_path":
        inventory[mutation[0]] = mutation[1]
    else:
        inventory["artifacts"][0][mutation[0]] = mutation[1]  # type: ignore[index]
    with pytest.raises(ReleaseGuardError) as captured:
        validate_blinded_inventory(
            inventory, view_root=view, repository_root=repository, verify_files=False
        )
    assert captured.value.code == code


def test_partition_pseudonyms_are_stable_but_not_cross_partition_linkable() -> None:
    secret = b"x" * 32
    first = partition_pseudonym(secret, "CALIBRATION", "MASTER-001")
    assert first == partition_pseudonym(secret, "CALIBRATION", "MASTER-001")
    assert first != partition_pseudonym(secret, "VALIDATION_BLIND", "MASTER-001")
    assert "MASTER" not in first


def _assignments() -> list[dict[str, object]]:
    result = []
    for prefix, partition in (("C", "CALIBRATION"), ("V", "VALIDATION_BLIND")):
        result.append(
            {
                "partition_id": partition,
                "isolation_unit_ids": [f"{prefix}-ISO"],
                "master_identity_ids": [f"{prefix}-ID"],
                "pedigree_component_ids": [f"{prefix}-PED"],
                "event_ids": [f"{prefix}-EVENT"],
                "segment_ids": [f"{prefix}-SEG"],
                "duplicate_group_ids": [f"{prefix}-DUP"],
                "seed_prefixes": [f"{prefix}-SEED"],
                "negative_control_ids": [f"{prefix}-NEG"],
            }
        )
    return result


def _reference_pools() -> list[dict[str, object]]:
    return [
        {"pool_id": f"POOL-{index}", "member_ids": [f"REF-{index}"]}
        for index in range(6)
    ]


def test_partition_and_six_reference_pools_are_disjoint() -> None:
    audit = validate_partition_separation(_assignments(), _reference_pools())
    assert audit.partition_count == 2
    assert audit.reference_member_count == 6


def test_cross_partition_overlap_is_rejected() -> None:
    assignments = _assignments()
    assignments[1]["event_ids"] = ["C-EVENT"]
    with pytest.raises(ReleaseGuardError, match="cross_partition_overlap"):
        validate_partition_separation(assignments, _reference_pools())


def test_single_factor_pair_accepts_exactly_one_declared_change() -> None:
    left = {"run_id": "A", "scenario": {"density": "D_MAX"}, "seed": 7}
    right = {"run_id": "B", "scenario": {"density": "D_LOW"}, "seed": 7}
    assert validate_single_factor_pair(
        left, right, allowed_path="scenario.density"
    ) == ("scenario.density", "D_MAX", "D_LOW")
    invalid = deepcopy(right)
    invalid["seed"] = 8
    with pytest.raises(ReleaseGuardError, match="single_factor_diff_invalid"):
        validate_single_factor_pair(left, invalid, allowed_path="scenario.density")


def test_environment_is_reduced_and_private_values_are_rejected() -> None:
    assert sanitized_environment(
        {"PATH": "/usr/bin", "TOKEN": "secret"}, {"PATH"}
    ) == {"PATH": "/usr/bin"}
    with pytest.raises(ReleaseGuardError, match="private_environment_value"):
        sanitized_environment({"TMPDIR": "private_truth/tmp"}, {"TMPDIR"})


def test_recursive_scan_rejects_a_hidden_truth_field() -> None:
    with pytest.raises(ReleaseGuardError, match="truth_field_in_blinded_document"):
        scan_blinded_document({"nested": [{"expected_direction": "UP"}]})
