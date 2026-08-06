"""Normalisation biallélique et harmonisation de la référence avec l'étude."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.references.window import ExtractedReferenceWindow


CommandRunner = Callable[[Sequence[str], float], subprocess.CompletedProcess[str]]
DNA_ALLELE = re.compile(r"^[ACGT]+$")
HARMONIZATION_COLUMNS = (
    "STUDY_VARIANT_ORDER",
    "STUDY_VARIANT_ID",
    "ASSEMBLY",
    "CHROMOSOME",
    "POSITION_BP",
    "STUDY_A1",
    "STUDY_A2",
    "REFERENCE_ID",
    "REFERENCE_REF",
    "REFERENCE_ALT",
    "HARMONIZATION_STATUS",
    "COMMON_PHASE_ELIGIBLE",
    "IS_TARGET_VARIANT",
    "DETAIL_CODE",
)
COMPLEMENT = str.maketrans("ACGT", "TGCA")


class ReferenceHarmonizationError(RuntimeError):
    """Signale une configuration ou une exécution d'harmonisation invalide."""


class ReferenceHarmonizationIntegrityError(ReferenceHarmonizationError):
    """Signale des entrées modifiées ou génétiquement incohérentes."""


@dataclass(frozen=True)
class HarmonizedReferenceWindow:
    """Décrit la référence normalisée et son audit d'harmonisation publiés."""

    output_dir: Path
    manifest_path: Path
    vcf_path: Path
    index_path: Path
    harmonization_path: Path
    reference_variant_count: int
    study_variant_count: int
    common_variant_count: int


@dataclass(frozen=True)
class _StudyVariant:
    order: int
    variant_id: str
    chromosome: int
    position_bp: int
    allele_1: str
    allele_2: str


@dataclass(frozen=True)
class _ReferenceVariant:
    chromosome: int
    position_bp: int
    reference_id: str | None
    ref: str
    alt: str

    @property
    def canonical_key(self) -> tuple[int, int, str, str]:
        return (self.chromosome, self.position_bp, self.ref, self.alt)


def _default_runner(
    command: Sequence[str], timeout_seconds: float
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command), check=False, capture_output=True, text=True, timeout=timeout_seconds
    )


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ReferenceHarmonizationError(f"invalid_harmonization_parameter:{name}")
    return value


def _positive_number(value: float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ReferenceHarmonizationError(f"invalid_harmonization_parameter:{name}")
    return float(value)


def _resolve_executable(command: str) -> str:
    if not isinstance(command, str) or not command.strip():
        raise ReferenceHarmonizationError(
            "invalid_harmonization_parameter:bcftools_command"
        )
    executable = shutil.which(command)
    if executable is None:
        raise ReferenceHarmonizationError("bcftools_not_found")
    return executable


def _run_checked(
    runner: CommandRunner,
    command: Sequence[str],
    timeout_seconds: float,
    operation: str,
) -> subprocess.CompletedProcess[str]:
    try:
        completed = runner(command, timeout_seconds)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ReferenceHarmonizationError(f"bcftools_{operation}_failed") from error
    if completed.returncode != 0:
        raise ReferenceHarmonizationError(
            f"bcftools_{operation}_failed:{completed.returncode}"
        )
    return completed


def _validated_reference_manifest(
    reference: ExtractedReferenceWindow,
) -> dict[str, Any]:
    expected_paths = (reference.manifest_path, reference.vcf_path, reference.index_path)
    if any(path.is_symlink() or not path.is_file() for path in expected_paths):
        raise ReferenceHarmonizationIntegrityError("reference_window_file_invalid")
    try:
        manifest = read_json(reference.manifest_path)
        validate_json_document(manifest, "reference_window_manifest.schema.json")
    except (ValueError, DocumentValidationError) as error:
        raise ReferenceHarmonizationIntegrityError(
            "reference_window_manifest_invalid"
        ) from error
    expected_identity = {
        "chromosome": reference.chromosome,
        "sample_count": reference.sample_count,
        "variant_count": reference.variant_count,
    }
    for field, expected in expected_identity.items():
        if manifest[field] != expected:
            raise ReferenceHarmonizationIntegrityError(
                f"reference_window_identity_mismatch:{field}"
            )
    for role, path in (("vcf", reference.vcf_path), ("index", reference.index_path)):
        record = manifest["files"][role]
        if path.parent != reference.output_dir or path.name != record["filename"]:
            raise ReferenceHarmonizationIntegrityError(
                f"reference_window_path_mismatch:{role}"
            )
        if path.stat().st_size != record["size_bytes"]:
            raise ReferenceHarmonizationIntegrityError(
                f"reference_window_size_mismatch:{role}"
            )
        if sha256_file(path) != record["sha256"]:
            raise ReferenceHarmonizationIntegrityError(
                f"reference_window_sha256_mismatch:{role}"
            )
    return manifest


def _load_phasing_manifest(path: Path, study_bim_path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ReferenceHarmonizationIntegrityError("phasing_manifest_file_invalid")
    try:
        manifest = read_json(path)
        validate_json_document(manifest, "phasing_input_manifest.schema.json")
    except (ValueError, DocumentValidationError) as error:
        raise ReferenceHarmonizationIntegrityError("phasing_manifest_invalid") from error
    if manifest["files"]["bim"] != study_bim_path.name:
        raise ReferenceHarmonizationIntegrityError("phasing_manifest_bim_mismatch")
    return manifest


def _load_target_metadata(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ReferenceHarmonizationIntegrityError("target_metadata_file_invalid")
    try:
        metadata = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            raise ValueError("target metadata is not an object")
        validate_json_document(metadata, "target_variant_metadata.schema.json")
    except (OSError, ValueError, yaml.YAMLError, DocumentValidationError) as error:
        raise ReferenceHarmonizationIntegrityError("target_metadata_invalid") from error
    return metadata


def _load_study_variants(path: Path) -> list[_StudyVariant]:
    if path.is_symlink() or not path.is_file():
        raise ReferenceHarmonizationIntegrityError("study_bim_file_invalid")
    variants: list[_StudyVariant] = []
    for order, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        fields = line.split()
        if len(fields) != 6:
            raise ReferenceHarmonizationIntegrityError("study_bim_malformed")
        chromosome_text, variant_id, _, position_text, allele_1, allele_2 = fields
        try:
            chromosome = int(chromosome_text.removeprefix("chr"))
            position_bp = int(position_text)
        except ValueError as error:
            raise ReferenceHarmonizationIntegrityError(
                "study_bim_coordinate_invalid"
            ) from error
        normalized_alleles = (allele_1.upper(), allele_2.upper())
        if (
            not variant_id
            or position_bp <= 0
            or any(DNA_ALLELE.fullmatch(allele) is None for allele in normalized_alleles)
            or normalized_alleles[0] == normalized_alleles[1]
        ):
            raise ReferenceHarmonizationIntegrityError("study_bim_variant_invalid")
        variants.append(
            _StudyVariant(
                order,
                variant_id,
                chromosome,
                position_bp,
                normalized_alleles[0],
                normalized_alleles[1],
            )
        )
    if not variants or len({variant.variant_id for variant in variants}) != len(variants):
        raise ReferenceHarmonizationIntegrityError("study_bim_variant_ids_invalid")
    return variants


def _sample_ids(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
    operation: str,
) -> list[str]:
    completed = _run_checked(
        runner,
        [executable, "query", "--list-samples", str(vcf_path)],
        timeout_seconds,
        operation,
    )
    sample_ids = [line for line in completed.stdout.splitlines() if line]
    if len(sample_ids) != len(set(sample_ids)):
        raise ReferenceHarmonizationIntegrityError("reference_duplicate_sample_id")
    return sample_ids


def _reference_variants(
    executable: str,
    vcf_path: Path,
    chromosome: int,
    runner: CommandRunner,
    timeout_seconds: float,
) -> list[_ReferenceVariant]:
    completed = _run_checked(
        runner,
        [
            executable,
            "query",
            "--format",
            "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n",
            str(vcf_path),
        ],
        timeout_seconds,
        "query_normalized_variants",
    )
    variants: list[_ReferenceVariant] = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 5:
            raise ReferenceHarmonizationIntegrityError(
                "normalized_reference_variant_malformed"
            )
        chromosome_text, position_text, reference_id, ref, alt = fields
        try:
            observed_chromosome = int(chromosome_text.removeprefix("chr"))
            position_bp = int(position_text)
        except ValueError as error:
            raise ReferenceHarmonizationIntegrityError(
                "normalized_reference_coordinate_invalid"
            ) from error
        ref = ref.upper()
        alt = alt.upper()
        if (
            observed_chromosome != chromosome
            or position_bp <= 0
            or DNA_ALLELE.fullmatch(ref) is None
            or DNA_ALLELE.fullmatch(alt) is None
            or ref == alt
            or "," in alt
            or (
                len(ref) > 1
                and len(alt) > 1
                and (ref[0] == alt[0] or ref[-1] == alt[-1])
            )
        ):
            raise ReferenceHarmonizationIntegrityError(
                "normalized_reference_variant_invalid"
            )
        variants.append(
            _ReferenceVariant(
                observed_chromosome,
                position_bp,
                None if reference_id == "." else reference_id,
                ref,
                alt,
            )
        )
    if not variants:
        raise ReferenceHarmonizationIntegrityError("normalized_reference_empty")
    return variants


def _indexed_variant_count(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
) -> int:
    completed = _run_checked(
        runner,
        [executable, "index", "--nrecords", str(vcf_path)],
        timeout_seconds,
        "count_harmonized_variants",
    )
    try:
        count = int(completed.stdout.strip())
    except ValueError as error:
        raise ReferenceHarmonizationIntegrityError(
            "harmonized_reference_variant_count_invalid"
        ) from error
    if count <= 0:
        raise ReferenceHarmonizationIntegrityError("normalized_reference_empty")
    return count


def _validate_phased_genotypes(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
) -> None:
    completed = _run_checked(
        runner,
        [
            executable,
            "query",
            "--include",
            'GT!="mis" && GT~"/"',
            "--format",
            "%CHROM\\t%POS\\n",
            str(vcf_path),
        ],
        timeout_seconds,
        "check_harmonized_phase",
    )
    if completed.stdout.strip():
        raise ReferenceHarmonizationIntegrityError(
            "harmonized_reference_unphased_genotype_detected"
        )


def _complement_pair(allele_1: str, allele_2: str) -> frozenset[str]:
    return frozenset((allele_1.translate(COMPLEMENT), allele_2.translate(COMPLEMENT)))


def _row_for_variant(
    study: _StudyVariant,
    candidates: list[_ReferenceVariant],
    assembly: str,
    target_variant_id: str,
) -> dict[str, Any]:
    study_alleles = frozenset((study.allele_1, study.allele_2))
    exact = [
        candidate
        for candidate in candidates
        if frozenset((candidate.ref, candidate.alt)) == study_alleles
    ]
    detail_code: str | None = None
    matched: _ReferenceVariant | None = None
    if len(exact) == 1:
        matched = exact[0]
        status = (
            "MATCHED_DIRECT"
            if (study.allele_1, study.allele_2) == (matched.ref, matched.alt)
            else "MATCHED_SWAPPED"
        )
        eligible = True
    elif len(exact) > 1:
        status = "AMBIGUOUS_REFERENCE_MATCH"
        eligible = False
        detail_code = "MULTIPLE_REFERENCE_RECORDS_FOR_CANONICAL_KEY"
    elif not candidates:
        status = "STUDY_ONLY"
        eligible = False
        detail_code = "NO_REFERENCE_VARIANT_AT_POSITION"
    else:
        status = "ALLELE_MISMATCH"
        eligible = False
        reference_pairs = [frozenset((candidate.ref, candidate.alt)) for candidate in candidates]
        if _complement_pair(study.allele_1, study.allele_2) in reference_pairs:
            detail_code = "STRAND_COMPLEMENT_NOT_APPLIED"
        else:
            detail_code = "REFERENCE_POSITION_WITH_DIFFERENT_ALLELES"
    return {
        "STUDY_VARIANT_ORDER": str(study.order),
        "STUDY_VARIANT_ID": study.variant_id,
        "ASSEMBLY": assembly,
        "CHROMOSOME": str(study.chromosome),
        "POSITION_BP": str(study.position_bp),
        "STUDY_A1": study.allele_1,
        "STUDY_A2": study.allele_2,
        "REFERENCE_ID": matched.reference_id if matched is not None else None,
        "REFERENCE_REF": matched.ref if matched is not None else None,
        "REFERENCE_ALT": matched.alt if matched is not None else None,
        "HARMONIZATION_STATUS": status,
        "COMMON_PHASE_ELIGIBLE": eligible,
        "IS_TARGET_VARIANT": study.variant_id == target_variant_id,
        "DETAIL_CODE": detail_code,
    }


def _mark_duplicate_study_matches(rows: list[dict[str, Any]]) -> None:
    rows_by_key: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        if row["COMMON_PHASE_ELIGIBLE"]:
            key = (
                row["CHROMOSOME"],
                row["POSITION_BP"],
                row["REFERENCE_REF"],
                row["REFERENCE_ALT"],
            )
            rows_by_key.setdefault(key, []).append(row)
    for duplicates in rows_by_key.values():
        if len(duplicates) > 1:
            for row in duplicates:
                row["HARMONIZATION_STATUS"] = "DUPLICATE_STUDY_MATCH"
                row["COMMON_PHASE_ELIGIBLE"] = False
                row["DETAIL_CODE"] = "MULTIPLE_STUDY_VARIANTS_FOR_CANONICAL_KEY"


def _write_harmonization(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HARMONIZATION_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    column: (
                        "true"
                        if row[column] is True
                        else "false"
                        if row[column] is False
                        else ""
                        if row[column] is None
                        else row[column]
                    )
                    for column in HARMONIZATION_COLUMNS
                }
            )


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "filename": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _result_from_manifest(
    output_dir: Path, manifest: dict[str, Any]
) -> HarmonizedReferenceWindow:
    return HarmonizedReferenceWindow(
        output_dir=output_dir,
        manifest_path=output_dir / "reference_harmonization_manifest.json",
        vcf_path=output_dir / manifest["files"]["vcf"]["filename"],
        index_path=output_dir / manifest["files"]["index"]["filename"],
        harmonization_path=output_dir
        / manifest["files"]["harmonization"]["filename"],
        reference_variant_count=manifest["reference_variant_count"],
        study_variant_count=manifest["study_variant_count"],
        common_variant_count=manifest["common_variant_count"],
    )


def harmonize_reference_window(
    reference: ExtractedReferenceWindow,
    study_bim_path: Path,
    phasing_manifest_path: Path,
    target_metadata_path: Path,
    output_dir: Path,
    *,
    bcftools_command: str = "bcftools",
    threads: int = 1,
    timeout_seconds: float = 7_200,
    minimum_common_variants: int = 1,
    command_runner: CommandRunner = _default_runner,
) -> HarmonizedReferenceWindow:
    """Décompose la référence et audite la jointure chr+position+REF+ALT."""
    resolved_threads = _positive_integer(threads, "threads")
    minimum_common = _positive_integer(
        minimum_common_variants, "minimum_common_variants"
    )
    timeout = _positive_number(timeout_seconds, "timeout_seconds")
    if output_dir.exists() or output_dir.is_symlink():
        raise ReferenceHarmonizationError("harmonization_output_exists")
    reference_manifest = _validated_reference_manifest(reference)
    phasing_manifest = _load_phasing_manifest(phasing_manifest_path, study_bim_path)
    target_metadata = _load_target_metadata(target_metadata_path)
    study_variants = _load_study_variants(study_bim_path)
    expected_identity = (
        reference_manifest["assembly"],
        reference.chromosome,
    )
    if (
        (phasing_manifest["assembly"], phasing_manifest["chromosome"])
        != expected_identity
        or (target_metadata["assembly"], target_metadata["chromosome"])
        != expected_identity
    ):
        raise ReferenceHarmonizationIntegrityError("assembly_or_chromosome_mismatch")
    if len(study_variants) != phasing_manifest["variant_count"]:
        raise ReferenceHarmonizationIntegrityError("study_variant_count_mismatch")
    if any(variant.chromosome != reference.chromosome for variant in study_variants):
        raise ReferenceHarmonizationIntegrityError("study_chromosome_mismatch")
    start_bp = reference_manifest["region"]["start_bp"]
    end_bp = reference_manifest["region"]["end_bp"]
    if any(
        not start_bp <= variant.position_bp <= end_bp for variant in study_variants
    ):
        raise ReferenceHarmonizationIntegrityError("study_variant_outside_reference_window")
    target_matches = [
        variant
        for variant in study_variants
        if variant.variant_id == phasing_manifest["target_variant_id"]
    ]
    if len(target_matches) != 1:
        raise ReferenceHarmonizationIntegrityError("target_variant_missing_or_ambiguous")
    target = target_matches[0]
    if (
        target.order != phasing_manifest["target_variant_order"]
        or target.variant_id != target_metadata["project_variant_id"]
        or target.position_bp != target_metadata["position_bp"]
        or frozenset((target.allele_1, target.allele_2))
        != frozenset((target_metadata["ref"], target_metadata["alt"]))
    ):
        raise ReferenceHarmonizationIntegrityError("target_variant_identity_mismatch")

    executable = _resolve_executable(bcftools_command)
    version_result = _run_checked(
        command_runner, [executable, "--version"], timeout, "version"
    )
    version_line = version_result.stdout.splitlines()[0] if version_result.stdout else ""
    if not version_line.startswith("bcftools "):
        raise ReferenceHarmonizationIntegrityError("bcftools_version_invalid")
    source_samples = _sample_ids(
        executable, reference.vcf_path, command_runner, timeout, "read_source_samples"
    )
    if len(source_samples) != reference.sample_count:
        raise ReferenceHarmonizationIntegrityError("reference_sample_count_mismatch")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    published = False
    try:
        decomposed_path = staging_dir / "reference.decomposed.bcf"
        vcf_name = (
            f"reference.harmonized.chr{reference.chromosome}.{start_bp}-{end_bp}.vcf.gz"
        )
        vcf_path = staging_dir / vcf_name
        index_path = staging_dir / f"{vcf_name}.tbi"
        _run_checked(
            command_runner,
            [
                executable,
                "norm",
                "--multiallelics",
                "-any",
                "--no-version",
                "--output-type",
                "b",
                "--threads",
                str(resolved_threads),
                "--output",
                str(decomposed_path),
                str(reference.vcf_path),
            ],
            timeout,
            "decompose_reference",
        )
        if decomposed_path.is_symlink() or not decomposed_path.is_file():
            raise ReferenceHarmonizationIntegrityError(
                "decomposed_reference_missing"
            )
        _run_checked(
            command_runner,
            [
                executable,
                "view",
                "--types",
                "snps,indels",
                "--min-alleles",
                "2",
                "--max-alleles",
                "2",
                "--no-version",
                "--output-type",
                "z",
                "--threads",
                str(resolved_threads),
                "--output",
                str(vcf_path),
                str(decomposed_path),
            ],
            timeout,
            "filter_reference_variants",
        )
        if vcf_path.is_symlink() or not vcf_path.is_file():
            raise ReferenceHarmonizationIntegrityError(
                "harmonized_reference_vcf_missing"
            )
        _run_checked(
            command_runner,
            [
                executable,
                "index",
                "--tbi",
                "--threads",
                str(resolved_threads),
                str(vcf_path),
            ],
            timeout,
            "index_harmonized_reference",
        )
        if index_path.is_symlink() or not index_path.is_file():
            raise ReferenceHarmonizationIntegrityError(
                "harmonized_reference_index_missing"
            )
        output_samples = _sample_ids(
            executable, vcf_path, command_runner, timeout, "read_harmonized_samples"
        )
        if output_samples != source_samples:
            raise ReferenceHarmonizationIntegrityError(
                "harmonized_reference_sample_set_or_order_mismatch"
            )
        reference_variants = _reference_variants(
            executable,
            vcf_path,
            reference.chromosome,
            command_runner,
            timeout,
        )
        indexed_variant_count = _indexed_variant_count(
            executable, vcf_path, command_runner, timeout
        )
        if indexed_variant_count != len(reference_variants):
            raise ReferenceHarmonizationIntegrityError(
                "harmonized_reference_variant_count_mismatch"
            )
        _validate_phased_genotypes(
            executable, vcf_path, command_runner, timeout
        )
        variants_by_position: dict[int, list[_ReferenceVariant]] = {}
        for variant in reference_variants:
            variants_by_position.setdefault(variant.position_bp, []).append(variant)
        rows = [
            _row_for_variant(
                study,
                variants_by_position.get(study.position_bp, []),
                reference_manifest["assembly"],
                phasing_manifest["target_variant_id"],
            )
            for study in study_variants
        ]
        _mark_duplicate_study_matches(rows)
        target_row = next(row for row in rows if row["IS_TARGET_VARIANT"])
        if target_row["HARMONIZATION_STATUS"] not in {
            "MATCHED_DIRECT",
            "MATCHED_SWAPPED",
            "STUDY_ONLY",
        }:
            raise ReferenceHarmonizationIntegrityError(
                "target_variant_harmonization_ambiguous"
            )
        common_variant_count = sum(row["COMMON_PHASE_ELIGIBLE"] for row in rows)
        if common_variant_count < minimum_common:
            raise ReferenceHarmonizationIntegrityError(
                "insufficient_common_variants_for_phasing"
            )
        harmonization_path = staging_dir / "variant_harmonization.tsv"
        _write_harmonization(harmonization_path, rows)
        validate_tsv_table(harmonization_path, "variant_harmonization.schema.json")
        manifest = {
            "schema_version": "1.0.0",
            "created_at": utc_now(),
            "method_id": "canonical_coordinate_allele_harmonization_v1",
            "assembly": reference_manifest["assembly"],
            "chromosome": reference.chromosome,
            "region": {"start_bp": start_bp, "end_bp": end_bp},
            "target_variant_id": phasing_manifest["target_variant_id"],
            "sample_count": len(output_samples),
            "reference_variant_count": len(reference_variants),
            "study_variant_count": len(study_variants),
            "common_variant_count": common_variant_count,
            "source": {
                "window_manifest_sha256": sha256_file(reference.manifest_path),
                "vcf_sha256": reference_manifest["files"]["vcf"]["sha256"],
                "index_sha256": reference_manifest["files"]["index"]["sha256"],
            },
            "study": {
                "phasing_manifest_sha256": sha256_file(phasing_manifest_path),
                "bim_sha256": sha256_file(study_bim_path),
                "variant_set_id": phasing_manifest["variant_set_id"],
            },
            "tool": {
                "name": "bcftools",
                "version": version_line.removeprefix("bcftools "),
            },
            "files": {
                "vcf": _file_record(vcf_path),
                "index": _file_record(index_path),
                "harmonization": _file_record(harmonization_path),
            },
            "checks": {
                "input_integrity": "PASS",
                "assembly_and_region": "PASS",
                "sample_set_and_order": "PASS",
                "phased_genotypes": "PASS",
                "biallelic_sequence_variants": "PASS",
                "canonical_variant_identity": "PASS",
                "target_variant": (
                    "STUDY_ONLY"
                    if target_row["HARMONIZATION_STATUS"] == "STUDY_ONLY"
                    else "MATCHED"
                ),
                "minimum_common_variants": "PASS",
                "tabix_index": "PASS",
            },
        }
        validate_json_document(
            manifest, "reference_harmonization_manifest.schema.json"
        )
        atomic_write_json(
            staging_dir / "reference_harmonization_manifest.json", manifest
        )
        decomposed_path.unlink()
        os.replace(staging_dir, output_dir)
        published = True
        return _result_from_manifest(output_dir, manifest)
    except (DocumentValidationError, TableValidationError, OSError) as error:
        raise ReferenceHarmonizationError("harmonization_publication_failed") from error
    finally:
        if not published:
            shutil.rmtree(staging_dir, ignore_errors=True)
