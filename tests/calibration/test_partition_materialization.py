"""Tests de matérialisation privée et atomique des partitions D10."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.materialize_partitions import (
    materialize_reference_partitions,
)
from effet_fondateur.calibration.partitioning import CalibrationPartitionError


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PUBLIC_CATALOG_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "references"
    / "public_reference_catalog.json"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _synthetic_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    metadata_path = tmp_path / "population_metadata.txt"
    metadata_path.write_text(
        "FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
        "FAM_A SYN_001 0 0 1 P1 SP1\n"
        "FAM_B SYN_002 0 0 2 P1 SP1\n"
        "FAM_C SYN_003 0 0 1 P2 SP2\n"
        "FAM_A SYN_004 SYN_001 0 2 P1 SP1\n"
        "FAM_D SYN_005 0 0 1 P2 SP2\n"
        "FAM_E SYN_006 0 0 2 P2 SP2\n",
        encoding="utf-8",
    )
    index_path = tmp_path / "unrelated.sequence.index"
    index_path.write_text(
        "##FileDate=20261028\n"
        "#ENA_FILE_PATH\tSAMPLE_NAME\tPOPULATION\n"
        "path/1\tSYN_001\tP1\n"
        "path/2\tSYN_002\tP1\n"
        "path/3\tSYN_003\tP2\n",
        encoding="utf-8",
    )
    catalog = deepcopy(json.loads(PUBLIC_CATALOG_PATH.read_text(encoding="utf-8")))
    for source in catalog["sources"]:
        if source["source_id"] == "igsr_1kg_3202_population_metadata_20130606":
            source["fingerprint"]["sha256"] = _sha256(metadata_path)
        elif source["source_id"] == "igsr_1kg_2504_unrelated_sequence_index":
            source["fingerprint"]["sha256"] = _sha256(index_path)
    catalog_path = tmp_path / "public_reference_catalog.json"
    catalog_path.write_text(
        json.dumps(catalog, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    workspace_root = tmp_path / "calibration_workspace"
    return catalog_path, metadata_path, index_path, workspace_root


def test_materialization_writes_only_private_memberships_and_aggregate_manifest(
    tmp_path: Path,
) -> None:
    catalog_path, metadata_path, index_path, workspace_root = _synthetic_inputs(
        tmp_path
    )

    audit = materialize_reference_partitions(
        public_catalog_path=catalog_path,
        population_metadata_path=metadata_path,
        unrelated_index_path=index_path,
        workspace_root=workspace_root,
        partition_run_id="d10_synthetic_001",
        created_at="2026-10-28T12:00:00Z",
    )

    run_root = workspace_root / "d10_synthetic_001"
    manifest_path = run_root / "private_truth" / "reference_partition_manifest.json"
    manifest_text = manifest_path.read_text(encoding="utf-8")
    assert audit.source_identifier_count == 6
    assert audit.phasing_reference_count == 2
    assert "SYN_" not in manifest_text
    assert len(list((run_root / "private_truth" / "partitions").glob("*.ids.txt"))) == 3
    assert manifest_path.stat().st_mode & 0o077 == 0


def test_materialization_refuses_overwrite_and_corrupt_input(
    tmp_path: Path,
) -> None:
    catalog_path, metadata_path, index_path, workspace_root = _synthetic_inputs(
        tmp_path
    )
    arguments = {
        "public_catalog_path": catalog_path,
        "population_metadata_path": metadata_path,
        "unrelated_index_path": index_path,
        "workspace_root": workspace_root,
        "partition_run_id": "d10_synthetic_001",
        "created_at": "2026-10-28T12:00:00Z",
    }
    materialize_reference_partitions(**arguments)

    with pytest.raises(CalibrationPartitionError) as overwrite:
        materialize_reference_partitions(**arguments)
    assert overwrite.value.code == "partition_run_exists"

    corrupt_arguments = dict(arguments)
    corrupt_arguments["partition_run_id"] = "d10_synthetic_corrupt"
    metadata_path.write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(CalibrationPartitionError) as corrupt:
        materialize_reference_partitions(**corrupt_arguments)
    assert corrupt.value.code == "partition_input_sha256_mismatch"


def test_materialization_rejects_duplicate_catalog_keys(tmp_path: Path) -> None:
    catalog_path, metadata_path, index_path, workspace_root = _synthetic_inputs(
        tmp_path
    )
    catalog_text = catalog_path.read_text(encoding="utf-8")
    catalog_path.write_text(
        '{"catalog_id": "duplicate",' + catalog_text.lstrip()[1:],
        encoding="utf-8",
    )

    with pytest.raises(CalibrationPartitionError) as captured:
        materialize_reference_partitions(
            public_catalog_path=catalog_path,
            population_metadata_path=metadata_path,
            unrelated_index_path=index_path,
            workspace_root=workspace_root,
            partition_run_id="d10_duplicate_catalog",
            created_at="2026-10-28T12:00:00Z",
        )
    assert captured.value.code == "partition_catalog_duplicate_key"
