import csv
import json
import subprocess
from pathlib import Path
from typing import Sequence

import pytest
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.references import (
    ExtractedReferenceWindow,
    ReferenceHarmonizationIntegrityError,
    harmonize_reference_window,
)


def _reference_window(tmp_path: Path) -> ExtractedReferenceWindow:
    output_dir = tmp_path / "window"
    output_dir.mkdir()
    vcf_path = output_dir / "reference.chr19.1000-2000.vcf.gz"
    index_path = output_dir / "reference.chr19.1000-2000.vcf.gz.tbi"
    vcf_path.write_bytes(b"reference window")
    index_path.write_bytes(b"reference index")
    manifest = {
        "schema_version": "1.0.0",
        "created_at": "2026-08-06T10:00:00Z",
        "method_id": "bcftools_all_samples_region_extract_v1",
        "assembly": "GRCh38",
        "chromosome": 19,
        "region": {"start_bp": 1_000, "end_bp": 2_000},
        "sample_count": 2,
        "variant_count": 3,
        "sample_scope": "ALL_RELEASE_SAMPLES",
        "population_scope": "ALL_AVAILABLE_1KGP_POPULATIONS",
        "source": {
            "catalog_id": "synthetic_catalog_v1",
            "catalog_sha256": "a" * 64,
            "panel_id": "synthetic_panel_v1",
            "provider": "IGSR_1000_GENOMES",
            "release_id": "synthetic_release_v1",
            "cache_entry_key": "b" * 64,
            "cache_manifest_sha256": "c" * 64,
            "vcf_sha256": "d" * 64,
            "index_sha256": "e" * 64,
        },
        "tool": {"name": "bcftools", "version": "1.21"},
        "files": {
            "vcf": {
                "filename": vcf_path.name,
                "sha256": sha256_file(vcf_path),
                "size_bytes": vcf_path.stat().st_size,
            },
            "index": {
                "filename": index_path.name,
                "sha256": sha256_file(index_path),
                "size_bytes": index_path.stat().st_size,
            },
        },
        "checks": {
            "assembly": "PASS",
            "chromosome": "PASS",
            "sample_count": "PASS",
            "sample_set_and_order": "PASS",
            "phased_genotypes": "PASS",
            "tabix_index": "PASS",
        },
    }
    manifest_path = output_dir / "reference_window_manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return ExtractedReferenceWindow(
        output_dir=output_dir,
        manifest_path=manifest_path,
        vcf_path=vcf_path,
        index_path=index_path,
        chromosome=19,
        start_bp=1_000,
        end_bp=2_000,
        sample_count=2,
        variant_count=3,
        vcf_sha256=sha256_file(vcf_path),
        index_sha256=sha256_file(index_path),
    )


def _study_inputs(tmp_path: Path, bim_lines: list[str] | None = None):
    bim_path = tmp_path / "target_region.bim"
    lines = bim_lines or [
        "19\tprobe_direct\t0.1\t1100\tA\tG",
        "19\tprobe_swapped\t0.2\t1200\tC\tT",
        "19\ttarget_variant\t0.3\t1500\tG\tA",
    ]
    bim_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    phasing_manifest = {
        "schema_version": "1.0.0",
        "format": "PLINK_BED_BIM_FAM_WITH_CM",
        "assembly": "GRCh38",
        "chromosome": 19,
        "dataset_id": "target_region",
        "sample_set_id": "study_samples",
        "variant_set_id": "study_variants",
        "sample_count": 2,
        "variant_count": len(lines),
        "target_variant_id": "target_variant",
        "target_variant_order": len(lines),
        "map_id": "synthetic_map",
        "files": {
            "bed": "target_region.bed",
            "bim": bim_path.name,
            "fam": "target_region.fam",
            "genetic_map": "target_genetic_map.tsv",
        },
    }
    phasing_manifest_path = tmp_path / "phasing_input_manifest.json"
    phasing_manifest_path.write_text(json.dumps(phasing_manifest), encoding="utf-8")
    target_metadata = {
        "schema_version": "1.0.0",
        "assembly": "GRCh38",
        "gene": "GENE1",
        "chromosome": 19,
        "position_bp": 1_500,
        "ref": "G",
        "alt": "A",
        "transcript": "NM_TEST.1",
        "hgvs_c": "c.1G>A",
        "hgvs_p": "p.Gly1Arg",
        "project_variant_id": "target_variant",
        "rsid": None,
        "reference_validation": {
            "status": "CONFIRMED",
            "source_id": "synthetic_reference",
            "confirmed_at": "2026-08-06T10:00:00Z",
        },
        "annotation_validation": {
            "status": "CONFIRMED",
            "source_id": "synthetic_annotation",
            "confirmed_at": "2026-08-06T10:00:00Z",
        },
    }
    target_metadata_path = tmp_path / "target_variant.yaml"
    target_metadata_path.write_text(
        yaml.safe_dump(target_metadata, sort_keys=False), encoding="utf-8"
    )
    return bim_path, phasing_manifest_path, target_metadata_path


class FakeBcftools:
    def __init__(
        self,
        variants: tuple[str, ...] = (
            "chr19\t1100\trs1\tA\tG",
            "chr19\t1200\trs2\tT\tC",
            "chr19\t1300\trs3\tC\tG",
        ),
        output_samples: tuple[str, ...] = ("HG00001", "HG00002"),
    ) -> None:
        self.variants = variants
        self.output_samples = output_samples
        self.commands: list[list[str]] = []

    def __call__(
        self, command: Sequence[str], timeout_seconds: float
    ) -> subprocess.CompletedProcess[str]:
        arguments = list(command)
        self.commands.append(arguments)
        if arguments[1:] == ["--version"]:
            return subprocess.CompletedProcess(arguments, 0, "bcftools 1.21\n", "")
        if arguments[1:3] == ["query", "--list-samples"]:
            samples = (
                self.output_samples
                if "harmonized" in Path(arguments[-1]).name
                else ("HG00001", "HG00002")
            )
            return subprocess.CompletedProcess(
                arguments, 0, "".join(f"{sample}\n" for sample in samples), ""
            )
        if arguments[1] == "norm":
            Path(arguments[arguments.index("--output") + 1]).write_bytes(b"bcf")
            return subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[1] == "view":
            Path(arguments[arguments.index("--output") + 1]).write_bytes(
                b"normalized vcf"
            )
            return subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[1:3] == ["index", "--tbi"]:
            Path(f"{arguments[-1]}.tbi").write_bytes(b"tabix")
            return subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[1:3] == ["index", "--nrecords"]:
            return subprocess.CompletedProcess(
                arguments, 0, f"{len(self.variants)}\n", ""
            )
        if arguments[1] == "query" and "--include" in arguments:
            return subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[1] == "query" and "--format" in arguments:
            return subprocess.CompletedProcess(
                arguments, 0, "".join(f"{variant}\n" for variant in self.variants), ""
            )
        raise AssertionError(f"unexpected command: {arguments}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_harmonizes_direct_swapped_and_study_only_variants(tmp_path: Path) -> None:
    reference = _reference_window(tmp_path)
    bim_path, phasing_manifest_path, target_metadata_path = _study_inputs(tmp_path)
    runner = FakeBcftools()

    result = harmonize_reference_window(
        reference,
        bim_path,
        phasing_manifest_path,
        target_metadata_path,
        tmp_path / "harmonized",
        command_runner=runner,
    )

    assert result.reference_variant_count == 3
    assert result.study_variant_count == 3
    assert result.common_variant_count == 2
    rows = _read_rows(result.harmonization_path)
    assert [row["HARMONIZATION_STATUS"] for row in rows] == [
        "MATCHED_DIRECT",
        "MATCHED_SWAPPED",
        "STUDY_ONLY",
    ]
    assert rows[2]["IS_TARGET_VARIANT"] == "true"
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["checks"]["target_variant"] == "STUDY_ONLY"
    assert not (result.output_dir / "reference.decomposed.bcf").exists()
    norm_command = next(command for command in runner.commands if command[1] == "norm")
    assert norm_command[norm_command.index("--multiallelics") + 1] == "-any"


def test_marks_duplicate_study_mapping_as_ineligible(tmp_path: Path) -> None:
    reference = _reference_window(tmp_path)
    bim_path, phasing_manifest_path, target_metadata_path = _study_inputs(
        tmp_path,
        [
            "19\tprobe_a\t0.1\t1100\tA\tG",
            "19\tprobe_b\t0.1\t1100\tG\tA",
            "19\tprobe_common\t0.2\t1200\tC\tT",
            "19\ttarget_variant\t0.3\t1500\tG\tA",
        ],
    )

    result = harmonize_reference_window(
        reference,
        bim_path,
        phasing_manifest_path,
        target_metadata_path,
        tmp_path / "harmonized",
        command_runner=FakeBcftools(),
    )

    rows = _read_rows(result.harmonization_path)
    assert [row["HARMONIZATION_STATUS"] for row in rows[:2]] == [
        "DUPLICATE_STUDY_MATCH",
        "DUPLICATE_STUDY_MATCH",
    ]
    assert result.common_variant_count == 1


def test_completes_missing_plink_allele_from_unique_reference(tmp_path: Path) -> None:
    reference = _reference_window(tmp_path)
    bim_path, phasing_manifest_path, target_metadata_path = _study_inputs(
        tmp_path,
        [
            "19\tmonomorphic_probe\t0.1\t1100\tA\t0",
            "19\tprobe_common\t0.2\t1200\tC\tT",
            "19\ttarget_variant\t0.3\t1500\tG\tA",
        ],
    )

    result = harmonize_reference_window(
        reference,
        bim_path,
        phasing_manifest_path,
        target_metadata_path,
        tmp_path / "harmonized",
        command_runner=FakeBcftools(),
    )

    first_row = _read_rows(result.harmonization_path)[0]
    assert first_row["STUDY_A2"] == ""
    assert first_row["REFERENCE_REF"] == "A"
    assert first_row["REFERENCE_ALT"] == "G"
    assert first_row["HARMONIZATION_STATUS"] == "MATCHED_REFERENCE_COMPLETED"
    assert first_row["COMMON_PHASE_ELIGIBLE"] == "true"


def test_blocks_target_allele_collision_without_publishing(tmp_path: Path) -> None:
    reference = _reference_window(tmp_path)
    bim_path, phasing_manifest_path, target_metadata_path = _study_inputs(tmp_path)
    output_dir = tmp_path / "harmonized"

    with pytest.raises(
        ReferenceHarmonizationIntegrityError,
        match="target_variant_harmonization_ambiguous",
    ):
        harmonize_reference_window(
            reference,
            bim_path,
            phasing_manifest_path,
            target_metadata_path,
            output_dir,
            command_runner=FakeBcftools(
                variants=(
                    "chr19\t1100\trs1\tA\tG",
                    "chr19\t1200\trs2\tT\tC",
                    "chr19\t1500\trs_target_other\tG\tC",
                )
            ),
        )

    assert not output_dir.exists()


def test_blocks_sample_reordering_without_publishing(tmp_path: Path) -> None:
    reference = _reference_window(tmp_path)
    bim_path, phasing_manifest_path, target_metadata_path = _study_inputs(tmp_path)
    output_dir = tmp_path / "harmonized"

    with pytest.raises(
        ReferenceHarmonizationIntegrityError,
        match="sample_set_or_order_mismatch",
    ):
        harmonize_reference_window(
            reference,
            bim_path,
            phasing_manifest_path,
            target_metadata_path,
            output_dir,
            command_runner=FakeBcftools(output_samples=("HG00002", "HG00001")),
        )

    assert not output_dir.exists()


def test_detects_modified_reference_before_bcftools(tmp_path: Path) -> None:
    reference = _reference_window(tmp_path)
    bim_path, phasing_manifest_path, target_metadata_path = _study_inputs(tmp_path)
    reference.vcf_path.write_bytes(b"modified reference")
    runner = FakeBcftools()

    with pytest.raises(
        ReferenceHarmonizationIntegrityError, match="window_size_mismatch"
    ):
        harmonize_reference_window(
            reference,
            bim_path,
            phasing_manifest_path,
            target_metadata_path,
            tmp_path / "harmonized",
            command_runner=runner,
        )

    assert runner.commands == []
