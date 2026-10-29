"""Tests de matérialisation atomique D11 avec bcftools simulé."""

from __future__ import annotations

from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import shutil

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_base_haplotypes import (
    materialize_base_haplotypes,
)
from effet_fondateur.calibration.partitioning import membership_sha256
from effet_fondateur.calibration.simulation.base_haplotypes import (
    BaseHaplotypeError,
)


ROOT = Path(__file__).resolve().parents[2]
PUBLIC_CATALOG = ROOT / "calibration/references/public_reference_catalog.json"
TECHNICAL_CATALOG = ROOT / "config/references/1000g_high_coverage_grch38_phased.json"
SCHEMA = ROOT / "calibration/schemas/base_haplotype_manifest.schema.json"
DONORS = ("SYN_SOURCE_A", "SYN_SOURCE_B", "SYN_SOURCE_C", "SYN_SOURCE_D")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()


def _fixture(tmp_path: Path) -> dict[str, object]:
    source_vcf = (
        tmp_path
        / "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    )
    source_vcf.write_bytes(b"synthetic compressed placeholder\n")
    source_index = Path(f"{source_vcf}.tbi")
    source_index.write_bytes(b"synthetic index placeholder\n")
    source_text = tmp_path / "source_fixture.vcf"
    source_text.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(DONORS)
        + "\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0\t1|0\n"
        "chr19\t300\tv3\tG\tA\t.\tPASS\t.\tGT\t1|1\t0|0\t1|0\t0|1\n",
        encoding="utf-8",
    )

    metadata = tmp_path / "population_metadata.txt"
    metadata.write_text(
        "FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
        "F1 SYN_SOURCE_A 0 0 1 P1 SP1\n"
        "F2 SYN_SOURCE_B 0 0 2 P1 SP1\n"
        "F3 SYN_SOURCE_C 0 0 1 P2 SP2\n"
        "F4 SYN_SOURCE_D 0 0 2 P2 SP2\n",
        encoding="utf-8",
    )
    genetic_map = tmp_path / "chr19.grch38.tsv"
    genetic_map.write_text(
        "MAP_ID\tASSEMBLY\tCHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
        "map\tGRCh38\t19\t50\t0\n"
        "map\tGRCh38\t19\t400\t3\n",
        encoding="utf-8",
    )

    technical_catalog = deepcopy(
        json.loads(TECHNICAL_CATALOG.read_text(encoding="utf-8"))
    )
    chromosome_record = next(
        record
        for record in technical_catalog["panels"][0]["chromosomes"]
        if record["chromosome"] == 19
    )
    chromosome_record["vcf_md5"] = _md5(source_vcf)
    chromosome_record["index_md5"] = _md5(source_index)
    technical_catalog_path = tmp_path / "technical_catalog.json"
    technical_catalog_path.write_text(
        json.dumps(technical_catalog, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    public_catalog = deepcopy(json.loads(PUBLIC_CATALOG.read_text(encoding="utf-8")))
    map_source: dict[str, object] | None = None
    for source in public_catalog["sources"]:
        if source["source_id"] == "igsr_1kg_3202_population_metadata_20130606":
            source["fingerprint"]["sha256"] = _sha256(metadata)
        elif source["source_id"] == "igsr_1kg_3202_phased_grch38_20220422":
            source["catalog_binding"]["sha256"] = _sha256(technical_catalog_path)
        elif source["source_id"] == "shapeit4_hapmap_grch38_git_7a0cab7":
            map_source = source
    assert map_source is not None
    public_catalog_path = tmp_path / "public_reference_catalog.json"
    public_catalog_path.write_text(
        json.dumps(public_catalog, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    genetic_map_manifest = tmp_path / "genetic_map_manifest.json"
    genetic_map_manifest.write_text(
        json.dumps(
            {
                "map_id": "shapeit4_hapmap_grch38_7a0cab7",
                "catalog_sha256": map_source["catalog_binding"]["sha256"],
                "archive": {"sha256": map_source["fingerprint"]["sha256"]},
                "maps": {
                    "19": {
                        "path": "normalized/chr19.grch38.tsv",
                        "sha256": _sha256(genetic_map),
                    }
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    partition_root = tmp_path / "d10_partition"
    membership_root = partition_root / "private_truth/partitions"
    membership_root.mkdir(parents=True)
    membership_path = membership_root / "donor_calibration.ids.txt"
    membership_path.write_text("".join(f"{donor}\n" for donor in DONORS), encoding="utf-8")
    partition_manifest = {
        "partition_manifest_id": "d10.synthetic.manifest",
        "partitions": [
            {
                "role": "DONOR_CALIBRATION",
                "membership_artifact_path": "private_truth/partitions/donor_calibration.ids.txt",
                "membership_sha256": membership_sha256(DONORS),
                "member_count": len(DONORS),
            }
        ],
    }
    (partition_root / "private_truth/reference_partition_manifest.json").write_text(
        json.dumps(partition_manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    workspace_root = tmp_path / "calibration_workspace"

    def fake_bcftools(command: list[str]) -> None:
        output_path = Path(command[command.index("--output-file") + 1])
        if "--samples-file" in command:
            shutil.copyfile(source_text, output_path)
            return
        input_path = Path(command[-1])
        with output_path.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_handle, mtime=0
            ) as output_handle:
                output_handle.write(input_path.read_bytes())

    return {
        "public_catalog_path": public_catalog_path,
        "technical_catalog_path": technical_catalog_path,
        "source_vcf_path": source_vcf,
        "source_index_path": source_index,
        "population_metadata_path": metadata,
        "genetic_map_path": genetic_map,
        "genetic_map_manifest_path": genetic_map_manifest,
        "partition_root": partition_root,
        "workspace_root": workspace_root,
        "run_id": "d11_synthetic_test",
        "panel_role": "CALIBRATION",
        "chromosome": 19,
        "start_bp": 100,
        "end_bp": 300,
        "random_seed": 1002,
        "created_at": "2026-10-28T10:00:00Z",
        "command_runner": fake_bcftools,
    }


def test_materialization_is_private_atomic_and_schema_valid(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    audit = materialize_base_haplotypes(**arguments)
    run_root = arguments["workspace_root"] / arguments["run_id"]
    private_root = run_root / "private_truth/d11"
    manifest_path = private_root / "base_haplotype_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(manifest)

    assert audit.individual_count == 4
    assert audit.haplotype_count == 8
    assert audit.variant_count == 3
    assert not (run_root / ".work").exists()
    assert manifest_path.stat().st_mode & 0o077 == 0
    serialized_manifest = json.dumps(manifest)
    assert all(donor not in serialized_manifest for donor in DONORS)
    assert "/Users/" not in serialized_manifest
    assert "Desktop" not in serialized_manifest

    with gzip.open(private_root / "donor_provenance.tsv.gz", "rt") as handle:
        provenance = handle.read()
    assert all(donor in provenance for donor in DONORS)
    with gzip.open(private_root / "base_haplotypes.phased.vcf.gz", "rt") as handle:
        phased = handle.read()
    assert all(donor not in phased for donor in DONORS)
    assert "SYN_CAL_000001" in phased


def test_materialization_refuses_overwrite(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    materialize_base_haplotypes(**arguments)
    with pytest.raises(BaseHaplotypeError) as captured:
        materialize_base_haplotypes(**arguments)
    assert captured.value.code == "base_haplotype_run_exists"


def test_materialization_refuses_corrupt_source_hash(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["source_vcf_path"].write_bytes(b"changed")
    with pytest.raises(BaseHaplotypeError) as captured:
        materialize_base_haplotypes(**arguments)
    assert captured.value.code == "base_haplotype_source_vcf_md5_mismatch"
    assert not (arguments["workspace_root"] / arguments["run_id"]).exists()


def test_materialization_refuses_unbound_technical_catalog(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["technical_catalog_path"].write_text("{}\n", encoding="utf-8")
    with pytest.raises(BaseHaplotypeError) as captured:
        materialize_base_haplotypes(**arguments)
    assert captured.value.code == "base_haplotype_technical_catalog_sha256_mismatch"


def test_materialization_refuses_corrupt_normalized_map(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["genetic_map_path"].write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(BaseHaplotypeError) as captured:
        materialize_base_haplotypes(**arguments)
    assert captured.value.code == "base_haplotype_genetic_map_sha256_mismatch"


def test_materialization_cleans_temporary_run_after_external_failure(
    tmp_path: Path,
) -> None:
    arguments = _fixture(tmp_path)

    def fail_command(_command: list[str]) -> None:
        raise BaseHaplotypeError("base_haplotype_external_command_failed")

    arguments["command_runner"] = fail_command
    with pytest.raises(BaseHaplotypeError):
        materialize_base_haplotypes(**arguments)
    workspace = arguments["workspace_root"]
    assert not (workspace / arguments["run_id"]).exists()
    assert not list(workspace.glob(".d11_synthetic_test.*"))


def test_materialization_rejects_duplicate_public_catalog_key(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    catalog_path = arguments["public_catalog_path"]
    original = catalog_path.read_text(encoding="utf-8")
    catalog_path.write_text(
        '{"catalog_id":"duplicate",' + original.lstrip()[1:],
        encoding="utf-8",
    )
    with pytest.raises(BaseHaplotypeError) as captured:
        materialize_base_haplotypes(**arguments)
    assert captured.value.code == "base_haplotype_json_duplicate_key"


def test_materialization_rejects_validation_role_with_calibration_partition(
    tmp_path: Path,
) -> None:
    arguments = _fixture(tmp_path)
    arguments["panel_role"] = "VALIDATION"
    with pytest.raises(BaseHaplotypeError) as captured:
        materialize_base_haplotypes(**arguments)
    assert captured.value.code == "base_haplotype_partition_role_missing"
