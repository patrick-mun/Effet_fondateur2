"""Tests de matérialisation atomique de l'observation sur puce D11b."""

from __future__ import annotations

import gzip
import hashlib
import json
import sqlite3
from copy import deepcopy
from pathlib import Path
import zipfile

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_chip_observation import (
    materialize_chip_observation,
)
from effet_fondateur.calibration.simulation.chip_markers import ChipMarkerError


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "calibration/schemas/chip_marker_harmonization_manifest.schema.json"
ANNOTATION_COLUMNS = (
    "ProbeSet_ID", "Chr_id", "Start", "Stop", "Strand", "dbSNP_RS_ID",
    "Strand_Vs_dbSNP", "Probe_Count", "Cytoband", "ChrX_PAR", "Flank",
    "Allele_A", "Allele_B", "Associated_Gene", "Genetic_Map",
    "Microsatellite", "Enzyme_Fragment", "CNV_Record", "SNP_Interference",
    "Percent_GC", "Heterozygous_Allele_Frequencies",
    "Allele_Frequency_Count", "Allele_Frequencies", "Minor_Allele",
    "Minor_Allele_Frequency", "OMIM", "ProbeSet_Type", "Probe_GC_A",
    "Process_Flag", "Fragment_GC",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _official_package_fixture(
    tmp_path: Path,
    *,
    array_type: str = "CytoScan750k_Accel",
) -> Path:
    database_path = tmp_path / "CytoScan750k_Accel_Array.na36.annot.db"
    with sqlite3.connect(database_path) as connection:
        connection.execute(
            "CREATE TABLE Information (key TEXT, value TEXT)"
        )
        information = {
            "array_type": array_type,
            "genome-version": "hg38",
            "genome-version-ncbi": "38",
            "genome-version-ensembl-url-parameter": "GRCh38",
            "netaffx-annotation-netaffx-build": "36",
            "annotation_row_count": "4",
        }
        connection.executemany(
            "INSERT INTO Information VALUES (?, ?)", information.items()
        )
        connection.execute(
            "CREATE TABLE Annotations ("
            + ", ".join(f'"{column}"' for column in ANNOTATION_COLUMNS)
            + ")"
        )
        rows = []
        for marker_id, chromosome, position, allele_a, allele_b, strand in (
            ("AX-1", 19, 100, "A", "G", "+"),
            ("AX-2", 19, 200, "A", "G", "+"),
            ("AX-3", 19, 300, "C", "G", "+"),
            ("AX-X", 23, 400, "C", "T", "+"),
        ):
            row: list[object | None] = [None] * len(ANNOTATION_COLUMNS)
            for column, value in {
                "ProbeSet_ID": marker_id,
                "Chr_id": chromosome,
                "Start": position,
                "Stop": position,
                "Strand": strand,
                "Allele_A": allele_a,
                "Allele_B": allele_b,
                "ProbeSet_Type": "SNP",
            }.items():
                row[ANNOTATION_COLUMNS.index(column)] = value
            rows.append(row)
        connection.executemany(
            "INSERT INTO Annotations VALUES ("
            + ",".join("?" for _ in ANNOTATION_COLUMNS)
            + ")",
            rows,
        )
    package_path = tmp_path / "CytoScan750K_Accel_NA36.r1.zip"
    with zipfile.ZipFile(package_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(database_path, arcname=database_path.name)
        archive.writestr(
            "CytoScan750K_Accel_Array.snplist.r3.txt",
            "probeset_id\nAX-1\n",
        )
    database_path.unlink()
    return package_path


def _fixture(tmp_path: Path) -> dict[str, object]:
    dense_root = tmp_path / "calibration_workspace/dense_d11a"
    private_root = dense_root / "private_truth/d11"
    private_root.mkdir(parents=True)
    dense_path = private_root / "base_genotypes.unphased.vcf.gz"
    with gzip.open(dense_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            "SYN_CAL_000001\tSYN_CAL_000002\n"
            "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n"
            "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1\n"
            "chr19\t300\tv3\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/1\n"
        )
    manifest = {
        "base_haplotype_manifest_id": "dense_test.d11",
        "panel": {"role": "CALIBRATION", "individual_count": 2},
        "source": {
            "source_id": "igsr_1kg_3202_phased_grch38_20220422",
            "assembly": "GRCh38",
        },
        "method": {"variant_policy": "BIALLELIC_SNV_COMPLETE_PHASED"},
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 300},
        "artifacts": [
            {
                "artifact_id": "dense_test.unphased",
                "role": "UNPHASED_GENOTYPES",
                "path": "private_truth/d11/base_genotypes.unphased.vcf.gz",
                "sha256": _sha256(dense_path),
                "record_count": 3,
            }
        ],
        "integrity": {"variant_count": 3},
    }
    (private_root / "base_haplotype_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    package_path = _official_package_fixture(tmp_path)
    public_catalog_path = tmp_path / "public_reference_catalog.json"
    public_catalog_path.write_text(
        json.dumps(
            {
                "catalog_state": "FROZEN",
                "scope": {"real_study_data_allowed": False},
                "sources": [
                    {
                        "source_id": "thermofisher_cytoscan750k_accel_na36_r1",
                        "source_type": "ARRAY_DESIGN_MANIFEST",
                        "version": "NA36.r1",
                        "allowed_roles": ["OBSERVATION_MARKER_UNIVERSE"],
                        "origin": {
                            "provider": "THERMO_FISHER",
                            "canonical_url": (
                                "https://www.thermofisher.com/order/catalog/"
                                "product/952539"
                            ),
                        },
                        "assembly": {"name": "GRCh38"},
                        "rights": {
                            "review_status": "TERMS_VERIFIED",
                            "use_decision": "ALLOWED_LOCAL_CALIBRATION",
                            "redistribution_allowed_by_project": False,
                        },
                        "acquisition": {
                            "local_sha256_required_before_use": True
                        },
                        "fingerprint": {
                            "algorithm": "SHA-256",
                            "sha256": _sha256(package_path),
                            "verification": "PROJECT_VERIFIED",
                        },
                    }
                ]
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "public_catalog_path": public_catalog_path,
        "dense_root": dense_root,
        "official_chip_package_path": package_path,
        "expected_official_package_sha256": _sha256(package_path),
        "workspace_root": tmp_path / "calibration_workspace",
        "run_id": "d11b_synthetic_test",
        "created_at": "2026-10-28T10:00:00Z",
    }


def test_materialization_binds_dense_truth_and_is_schema_valid(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    manifest = materialize_chip_observation(**arguments)
    run_root = arguments["workspace_root"] / arguments["run_id"]
    private_root = run_root / "private_truth/d11b"
    jsonschema.Draft202012Validator(
        json.loads(SCHEMA.read_text(encoding="utf-8")),
        format_checker=jsonschema.FormatChecker(),
    ).validate(manifest)

    assert manifest["integrity"]["matched_marker_count"] == 1
    assert manifest["integrity"]["excluded_marker_count"] == 2
    assert manifest["integrity"]["normalized_autosomal_marker_count"] == 3
    assert manifest["integrity"]["package_non_autosomal_snp_count"] == 1
    assert manifest["confidentiality"]["real_sample_files_read"] is False
    assert all(path.stat().st_mode & 0o077 == 0 for path in private_root.iterdir())
    with gzip.open(
        private_root / "observed_chip_genotypes.vcf.gz", "rt", encoding="utf-8"
    ) as handle:
        observed = handle.read()
    assert "AX-1" in observed
    assert "\t200\t" not in observed
    assert "\t300\t" not in observed
    serialized = json.dumps(manifest)
    assert str(tmp_path) not in serialized


def test_materialization_refuses_wrong_package_hash_atomically(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["expected_official_package_sha256"] = "0" * 64
    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)
    assert captured.value.code == "chip_official_package_sha256_mismatch"
    assert not (arguments["workspace_root"] / arguments["run_id"]).exists()


def test_materialization_refuses_wrong_official_database_identity(
    tmp_path: Path,
) -> None:
    arguments = _fixture(tmp_path)
    invalid_root = tmp_path / "invalid_package"
    invalid_root.mkdir()
    invalid_package = _official_package_fixture(
        invalid_root,
        array_type="CytoScan750k_Array",
    )
    arguments["official_chip_package_path"] = invalid_package
    arguments["expected_official_package_sha256"] = _sha256(invalid_package)

    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)

    assert captured.value.code == "chip_official_annotation_metadata_mismatch"
    assert not (arguments["workspace_root"] / arguments["run_id"]).exists()


def test_materialization_refuses_sample_derived_package_before_read(
    tmp_path: Path,
) -> None:
    arguments = _fixture(tmp_path)
    forbidden_root = tmp_path / "RES_ACPA"
    forbidden_root.mkdir()
    forbidden_path = forbidden_root / "package.zip"
    forbidden_path.write_bytes(arguments["official_chip_package_path"].read_bytes())
    arguments["official_chip_package_path"] = forbidden_path
    arguments["expected_official_package_sha256"] = _sha256(forbidden_path)
    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)
    assert captured.value.code == "chip_official_package_source_forbidden"


def test_materialization_refuses_existing_run(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    materialize_chip_observation(**arguments)
    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)
    assert captured.value.code == "chip_harmonization_run_exists"


def test_materialization_refuses_non_synthetic_vcf_identifier(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    dense_path = (
        arguments["dense_root"]
        / "private_truth/d11/base_genotypes.unphased.vcf.gz"
    )
    with gzip.open(dense_path, "rt", encoding="utf-8") as handle:
        payload = handle.read().replace("SYN_CAL_000001", "REAL_SAMPLE")
    with gzip.open(dense_path, "wt", encoding="utf-8") as handle:
        handle.write(payload)
    dense_manifest_path = (
        arguments["dense_root"]
        / "private_truth/d11/base_haplotype_manifest.json"
    )
    dense_manifest = json.loads(dense_manifest_path.read_text(encoding="utf-8"))
    dense_manifest["artifacts"][0]["sha256"] = _sha256(dense_path)
    dense_manifest_path.write_text(
        json.dumps(dense_manifest, indent=2) + "\n", encoding="utf-8"
    )
    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)
    assert captured.value.code == "chip_dense_non_synthetic_identifier_forbidden"


def test_materialization_refuses_dense_variant_count_mismatch(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    dense_manifest_path = (
        arguments["dense_root"]
        / "private_truth/d11/base_haplotype_manifest.json"
    )
    dense_manifest = json.loads(dense_manifest_path.read_text(encoding="utf-8"))
    dense_manifest["integrity"]["variant_count"] = 4
    dense_manifest_path.write_text(
        json.dumps(dense_manifest, indent=2) + "\n", encoding="utf-8"
    )
    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)
    assert captured.value.code == "chip_dense_variant_count_mismatch"


def test_materialization_requires_frozen_chip_catalog_entry(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["public_catalog_path"].write_text(
        json.dumps(
            {
                "catalog_state": "FROZEN",
                "scope": {"real_study_data_allowed": False},
                "sources": [],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(ChipMarkerError) as captured:
        materialize_chip_observation(**arguments)
    assert captured.value.code == "chip_public_catalog_source_missing"


@pytest.mark.parametrize(
    ("section", "field", "invalid_value"),
    [
        ("chip_design", "design_version", "NA36.r7"),
        ("chip_design", "sample_derived_marker_source", True),
        ("confidentiality", "real_sample_files_read", True),
        ("policy", "pipeline_marker_mode", "UNION"),
    ],
)
def test_schema_rejects_scientifically_invalid_contracts(
    tmp_path: Path,
    section: str,
    field: str,
    invalid_value: object,
) -> None:
    arguments = _fixture(tmp_path)
    manifest = materialize_chip_observation(**arguments)
    invalid = deepcopy(manifest)
    invalid[section][field] = invalid_value
    validator = jsonschema.Draft202012Validator(
        json.loads(SCHEMA.read_text(encoding="utf-8")),
        format_checker=jsonschema.FormatChecker(),
    )
    assert list(validator.iter_errors(invalid))
