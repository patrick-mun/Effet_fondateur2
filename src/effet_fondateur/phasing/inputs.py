"""Construction contrôlée des entrées étude, carte et pedigree pour SHAPEIT5."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
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
from effet_fondateur.phasing.shapeit5 import SHAPEIT5_CONTRACT


CommandRunner = Callable[[Sequence[str], float], subprocess.CompletedProcess[str]]
COMMON_STATUSES = {
    "MATCHED_DIRECT",
    "MATCHED_SWAPPED",
    "MATCHED_REFERENCE_COMPLETED",
}
VARIANT_SELECTION_COLUMNS = (
    "VARIANT_ORDER",
    "VARIANT_ID",
    "CHROMOSOME",
    "POSITION_BP",
    "HARMONIZATION_STATUS",
    "SHAPEIT_ROLE",
    "INCLUDED",
    "CANONICAL_REF",
    "CANONICAL_ALT",
    "EXCLUSION_CODE",
    "IS_TARGET_VARIANT",
)
SAMPLE_MAPPING_COLUMNS = (
    "VCF_ORDER",
    "SAMPLE_ID",
    "FID",
    "IID",
    "VCF_SAMPLE_ID",
    "PEDIGREE_INCLUDED",
    "VCF_FATHER_ID",
    "VCF_MOTHER_ID",
)


class Shapeit5InputError(ValueError):
    """Signale des entrées ou paramètres incompatibles avec le contrat 12.5."""


class Shapeit5InputExternalError(RuntimeError):
    """Signale un échec PLINK ou bcftools pendant la préparation."""


class Shapeit5InputBlockError(RuntimeError):
    """Signale des données scientifiquement insuffisantes ou ambiguës."""


@dataclass(frozen=True)
class PreparedShapeit5Inputs:
    """Décrit les fichiers prêts à être consommés par l'adaptateur SHAPEIT5."""

    output_dir: Path
    manifest_path: Path
    study_vcf_path: Path
    study_index_path: Path
    genetic_map_path: Path
    pedigree_path: Path
    variant_selection_path: Path
    sample_mapping_path: Path
    sample_count: int
    study_variant_count: int
    common_variant_count: int
    pedigree_record_count: int


def _default_runner(
    command: Sequence[str], timeout_seconds: float
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command), check=False, capture_output=True, text=True, timeout=timeout_seconds
    )


def _positive_number(value: float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise Shapeit5InputError(f"invalid_shapeit5_input_parameter:{name}")
    return float(value)


def _resolve_executable(command: str, tool: str) -> str:
    if not isinstance(command, str) or not command.strip():
        raise Shapeit5InputError(f"{tool}_not_configured")
    executable = shutil.which(command)
    if executable is None:
        raise Shapeit5InputExternalError(f"{tool}_not_found")
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
        raise Shapeit5InputExternalError(f"{operation}_failed") from error
    if completed.returncode != 0:
        raise Shapeit5InputExternalError(
            f"{operation}_failed:{completed.returncode}"
        )
    return completed


def _read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise Shapeit5InputError("study_fam_malformed")
    if len({(row[0], row[1]) for row in rows}) != len(rows):
        raise Shapeit5InputError("study_fam_duplicate_sample")
    return rows


def _load_target_metadata(path: Path) -> dict[str, Any]:
    try:
        metadata = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(metadata, dict):
            raise ValueError("target metadata is not an object")
        validate_json_document(metadata, "target_variant_metadata.schema.json")
    except (OSError, ValueError, yaml.YAMLError, DocumentValidationError) as error:
        raise Shapeit5InputError("target_metadata_invalid") from error
    return metadata


def _load_manifest(path: Path, schema_name: str, error_code: str) -> dict[str, Any]:
    try:
        manifest = read_json(path)
        validate_json_document(manifest, schema_name)
    except (OSError, ValueError, json.JSONDecodeError, DocumentValidationError) as error:
        raise Shapeit5InputError(error_code) from error
    return manifest


def _selection_rows(
    harmonization_rows: list[dict[str, Any]], target_metadata: dict[str, Any]
) -> list[dict[str, Any]]:
    target_variant_id = target_metadata["project_variant_id"]
    rows: list[dict[str, Any]] = []
    for source in harmonization_rows:
        is_target = source["STUDY_VARIANT_ID"] == target_variant_id
        status = source["HARMONIZATION_STATUS"]
        if status in COMMON_STATUSES and source["COMMON_PHASE_ELIGIBLE"]:
            included = True
            role = "COMMON_TARGET" if is_target else "COMMON"
            canonical_ref = source["REFERENCE_REF"]
            canonical_alt = source["REFERENCE_ALT"]
            exclusion_code = None
        elif is_target and status == "STUDY_ONLY":
            included = True
            role = "RARE_TARGET"
            canonical_ref = target_metadata["ref"]
            canonical_alt = target_metadata["alt"]
            exclusion_code = None
        else:
            included = False
            role = "EXCLUDED"
            canonical_ref = None
            canonical_alt = None
            exclusion_code = (
                "study_only_without_confirmed_ref_alt"
                if status == "STUDY_ONLY"
                else f"harmonization_{status.lower()}"
            )
        rows.append(
            {
                "VARIANT_ORDER": source["STUDY_VARIANT_ORDER"],
                "VARIANT_ID": source["STUDY_VARIANT_ID"],
                "CHROMOSOME": source["CHROMOSOME"],
                "POSITION_BP": source["POSITION_BP"],
                "HARMONIZATION_STATUS": status,
                "SHAPEIT_ROLE": role,
                "INCLUDED": included,
                "CANONICAL_REF": canonical_ref,
                "CANONICAL_ALT": canonical_alt,
                "EXCLUSION_CODE": exclusion_code,
                "IS_TARGET_VARIANT": is_target,
            }
        )
    target_rows = [row for row in rows if row["IS_TARGET_VARIANT"]]
    if len(target_rows) != 1 or not target_rows[0]["INCLUDED"]:
        raise Shapeit5InputBlockError("target_variant_not_retained_for_shapeit5")
    included_rows = [row for row in rows if row["INCLUDED"]]
    common_count = sum(row["SHAPEIT_ROLE"] in {"COMMON", "COMMON_TARGET"} for row in rows)
    if len(included_rows) < 2:
        raise Shapeit5InputBlockError("fewer_than_two_shapeit5_study_variants")
    if common_count < 1:
        raise Shapeit5InputBlockError("no_common_variant_for_phase_common")
    return rows


def _sample_rows(
    fam_rows: list[list[str]], master_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[tuple[str, str, str]]]:
    master_by_plink = {(row["FID"], row["IID"]): row for row in master_rows}
    fam_keys = [(row[0], row[1]) for row in fam_rows]
    if any(key not in master_by_plink for key in fam_keys):
        raise Shapeit5InputError("study_fam_sample_missing_from_master")
    sample_ids = [master_by_plink[key]["SAMPLE_ID"] for key in fam_keys]
    if len(sample_ids) != len(set(sample_ids)):
        raise Shapeit5InputError("shapeit5_sample_id_collision")
    fam_key_set = set(fam_keys)
    mapping_rows: list[dict[str, Any]] = []
    pedigree_rows: list[tuple[str, str, str]] = []
    for order, fam_row in enumerate(fam_rows, start=1):
        fid, iid, pid, mid = fam_row[:4]
        master = master_by_plink[(fid, iid)]
        if (master["PID"], master["MID"]) != (pid, mid):
            raise Shapeit5InputError("study_fam_parent_mismatch")
        father = (
            master_by_plink[(fid, pid)]["SAMPLE_ID"]
            if pid != "0" and (fid, pid) in fam_key_set
            else None
        )
        mother = (
            master_by_plink[(fid, mid)]["SAMPLE_ID"]
            if mid != "0" and (fid, mid) in fam_key_set
            else None
        )
        pedigree_included = father is not None or mother is not None
        if pedigree_included:
            pedigree_rows.append(
                (master["SAMPLE_ID"], father or "NA", mother or "NA")
            )
        mapping_rows.append(
            {
                "VCF_ORDER": str(order),
                "SAMPLE_ID": master["SAMPLE_ID"],
                "FID": fid,
                "IID": iid,
                "VCF_SAMPLE_ID": master["SAMPLE_ID"],
                "PEDIGREE_INCLUDED": pedigree_included,
                "VCF_FATHER_ID": father,
                "VCF_MOTHER_ID": mother,
            }
        )
    return mapping_rows, pedigree_rows


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
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
                    for column in columns
                }
            )


def _write_shapeit_map(path: Path, map_rows: list[dict[str, Any]], chromosome: int) -> None:
    chromosome_rows = [
        row for row in map_rows if int(row["CHROMOSOME"]) == chromosome
    ]
    if len(chromosome_rows) < 2:
        raise Shapeit5InputBlockError("shapeit5_map_has_fewer_than_two_points")
    positions = [int(row["POSITION_BP"]) for row in chromosome_rows]
    positions_cm = [float(row["POSITION_CM"]) for row in chromosome_rows]
    if any(current <= previous for previous, current in zip(positions, positions[1:])):
        raise Shapeit5InputError("shapeit5_map_positions_not_increasing")
    if any(current < previous for previous, current in zip(positions_cm, positions_cm[1:])):
        raise Shapeit5InputBlockError("shapeit5_map_cm_not_monotonic")
    payload = "pos\tchr\tcM\n" + "".join(
        f"{row['POSITION_BP']}\t{chromosome}\t{row['POSITION_CM']}\n"
        for row in chromosome_rows
    )
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as gzip_handle:
            gzip_handle.write(payload.encode("utf-8"))


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "filename": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _stage_file_record(path: Path, stage_root: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(stage_root).as_posix(),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _bcftools_samples(
    executable: str,
    path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
) -> list[str]:
    result = _run_checked(
        runner,
        [executable, "query", "--list-samples", str(path)],
        timeout_seconds,
        "bcftools_list_samples",
    )
    samples = [line for line in result.stdout.splitlines() if line]
    if len(samples) != len(set(samples)):
        raise Shapeit5InputBlockError("vcf_duplicate_sample_id")
    return samples


def _bcftools_variants(
    executable: str,
    path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
) -> list[tuple[str, int, str, str, str]]:
    result = _run_checked(
        runner,
        [
            executable,
            "query",
            "--format",
            "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n",
            str(path),
        ],
        timeout_seconds,
        "bcftools_query_variants",
    )
    variants: list[tuple[str, int, str, str, str]] = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 5:
            raise Shapeit5InputBlockError("vcf_variant_record_malformed")
        try:
            position_bp = int(fields[1])
        except ValueError as error:
            raise Shapeit5InputBlockError("vcf_variant_position_invalid") from error
        variants.append((fields[0], position_bp, fields[2], fields[3], fields[4]))
    return variants


def _result(output_dir: Path, manifest: dict[str, Any]) -> PreparedShapeit5Inputs:
    files = manifest["files"]
    return PreparedShapeit5Inputs(
        output_dir=output_dir,
        manifest_path=output_dir / "shapeit5_inputs_manifest.json",
        study_vcf_path=output_dir / files["study_vcf"]["filename"],
        study_index_path=output_dir / files["study_index"]["filename"],
        genetic_map_path=output_dir / files["genetic_map"]["filename"],
        pedigree_path=output_dir / files["pedigree"]["filename"],
        variant_selection_path=output_dir / files["variant_selection"]["filename"],
        sample_mapping_path=output_dir / files["sample_mapping"]["filename"],
        sample_count=manifest["sample_count"],
        study_variant_count=manifest["study_variant_count"],
        common_variant_count=manifest["common_variant_count"],
        pedigree_record_count=manifest["pedigree_record_count"],
    )


def build_shapeit5_inputs(
    *,
    study_bed_path: Path,
    study_bim_path: Path,
    study_fam_path: Path,
    phasing_manifest_path: Path,
    target_genetic_map_path: Path,
    samples_master_path: Path,
    target_metadata_path: Path,
    harmonization_path: Path,
    harmonization_manifest_path: Path,
    reference_vcf_path: Path,
    reference_index_path: Path,
    reference_variant_set_id: str,
    stage_root: Path,
    output_dir: Path,
    plink_command: str = "plink",
    bcftools_command: str = "bcftools",
    timeout_seconds: float = 300,
    command_runner: CommandRunner = _default_runner,
) -> PreparedShapeit5Inputs:
    """Construit les fichiers 12.5 sans exécuter le logiciel de phasage."""
    timeout = _positive_number(timeout_seconds, "timeout_seconds")
    if output_dir.exists() or output_dir.is_symlink():
        raise Shapeit5InputError("shapeit5_inputs_output_exists")
    source_paths = (
        study_bed_path,
        study_bim_path,
        study_fam_path,
        phasing_manifest_path,
        target_genetic_map_path,
        samples_master_path,
        target_metadata_path,
        harmonization_path,
        harmonization_manifest_path,
        reference_vcf_path,
        reference_index_path,
    )
    if any(path.is_symlink() or not path.is_file() for path in source_paths):
        raise Shapeit5InputError("shapeit5_input_file_invalid")
    phasing_manifest = _load_manifest(
        phasing_manifest_path,
        "phasing_input_manifest.schema.json",
        "phasing_input_manifest_invalid",
    )
    harmonization_manifest = _load_manifest(
        harmonization_manifest_path,
        "reference_harmonization_manifest.schema.json",
        "reference_harmonization_manifest_invalid",
    )
    target_metadata = _load_target_metadata(target_metadata_path)
    if (
        len(
            {
                path.with_suffix("")
                for path in (study_bed_path, study_bim_path, study_fam_path)
            }
        )
        != 1
    ):
        raise Shapeit5InputError("study_plink_files_do_not_share_prefix")
    if (
        phasing_manifest["files"]["bed"] != study_bed_path.name
        or phasing_manifest["files"]["bim"] != study_bim_path.name
        or phasing_manifest["files"]["fam"] != study_fam_path.name
        or phasing_manifest["files"]["genetic_map"] != target_genetic_map_path.name
    ):
        raise Shapeit5InputError("phasing_manifest_file_mismatch")
    if (
        phasing_manifest["assembly"] != "GRCh38"
        or harmonization_manifest["assembly"] != "GRCh38"
        or target_metadata["assembly"] != "GRCh38"
    ):
        raise Shapeit5InputError("shapeit5_input_assembly_mismatch")
    chromosome = phasing_manifest["chromosome"]
    if (
        harmonization_manifest["chromosome"] != chromosome
        or target_metadata["chromosome"] != chromosome
    ):
        raise Shapeit5InputError("shapeit5_input_chromosome_mismatch")
    if (
        phasing_manifest["target_variant_id"]
        != target_metadata["project_variant_id"]
    ):
        raise Shapeit5InputError("shapeit5_target_identity_mismatch")
    if (
        harmonization_manifest["files"]["vcf"]["sha256"]
        != sha256_file(reference_vcf_path)
        or harmonization_manifest["files"]["index"]["sha256"]
        != sha256_file(reference_index_path)
        or harmonization_manifest["files"]["harmonization"]["sha256"]
        != sha256_file(harmonization_path)
    ):
        raise Shapeit5InputError("shapeit5_harmonization_artifact_mismatch")
    harmonization_table = validate_tsv_table(
        harmonization_path, "variant_harmonization.schema.json"
    )
    selection_rows = _selection_rows(harmonization_table.rows, target_metadata)
    included_rows = [row for row in selection_rows if row["INCLUDED"]]
    common_rows = [
        row
        for row in included_rows
        if row["SHAPEIT_ROLE"] in {"COMMON", "COMMON_TARGET"}
    ]
    target_row = next(row for row in included_rows if row["IS_TARGET_VARIANT"])
    fam_rows = _read_fam(study_fam_path)
    samples_table = validate_tsv_table(samples_master_path, "samples_master.schema.json")
    sample_rows, pedigree_rows = _sample_rows(fam_rows, samples_table.rows)
    map_table = validate_tsv_table(
        target_genetic_map_path, "target_genetic_map.schema.json"
    )
    plink_executable = _resolve_executable(plink_command, "plink")
    bcftools_executable = _resolve_executable(bcftools_command, "bcftools")
    plink_version_result = _run_checked(
        command_runner,
        [plink_executable, "--version"],
        timeout,
        "plink_version",
    )
    bcftools_version_result = _run_checked(
        command_runner,
        [bcftools_executable, "--version"],
        timeout,
        "bcftools_version",
    )
    plink_version = (
        plink_version_result.stdout.splitlines()[0]
        if plink_version_result.stdout
        else ""
    )
    bcftools_version = (
        bcftools_version_result.stdout.splitlines()[0]
        if bcftools_version_result.stdout
        else ""
    )
    if not plink_version or not bcftools_version.startswith("bcftools "):
        raise Shapeit5InputExternalError("shapeit5_input_tool_version_invalid")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    published = False
    try:
        selection_path = staging_dir / "shapeit5_variant_selection.tsv"
        sample_mapping_path = staging_dir / "shapeit5_sample_mapping.tsv"
        pedigree_path = staging_dir / "shapeit5.pedigree.tsv"
        genetic_map_path = staging_dir / "shapeit5.genetic_map.tsv.gz"
        _write_tsv(selection_path, VARIANT_SELECTION_COLUMNS, selection_rows)
        _write_tsv(sample_mapping_path, SAMPLE_MAPPING_COLUMNS, sample_rows)
        validate_tsv_table(selection_path, "shapeit5_variant_selection.schema.json")
        validate_tsv_table(sample_mapping_path, "shapeit5_sample_mapping.schema.json")
        pedigree_path.write_text(
            "".join("\t".join(row) + "\n" for row in pedigree_rows),
            encoding="utf-8",
        )
        _write_shapeit_map(genetic_map_path, map_table.rows, chromosome)

        extract_path = staging_dir / "variants.extract.txt"
        reference_alleles_path = staging_dir / "reference_alleles.tsv"
        update_ids_path = staging_dir / "sample_ids.update.tsv"
        extract_path.write_text(
            "".join(f"{row['VARIANT_ID']}\n" for row in included_rows),
            encoding="utf-8",
        )
        reference_alleles_path.write_text(
            "".join(
                f"{row['CANONICAL_REF']}\t{row['VARIANT_ID']}\n"
                for row in included_rows
            ),
            encoding="utf-8",
        )
        update_ids_path.write_text(
            "".join(
                f"{row['FID']}\t{row['IID']}\t{row['SAMPLE_ID']}\t{row['SAMPLE_ID']}\n"
                for row in sample_rows
            ),
            encoding="utf-8",
        )
        study_prefix = staging_dir / "study.shapeit5"
        _run_checked(
            command_runner,
            [
                plink_executable,
                "--bfile",
                str(study_bed_path.with_suffix("")),
                "--extract",
                str(extract_path),
                "--update-ids",
                str(update_ids_path),
                "--a2-allele",
                str(reference_alleles_path),
                "1",
                "2",
                "--keep-allele-order",
                "--real-ref-alleles",
                "--output-chr",
                "chrM",
                "--recode",
                "vcf-iid",
                "bgz",
                "--out",
                str(study_prefix),
            ],
            timeout,
            "plink_shapeit5_vcf_export",
        )
        study_vcf_path = Path(f"{study_prefix}.vcf.gz")
        study_index_path = Path(f"{study_vcf_path}.tbi")
        if study_vcf_path.is_symlink() or not study_vcf_path.is_file():
            raise Shapeit5InputExternalError("plink_shapeit5_vcf_missing")
        _run_checked(
            command_runner,
            [bcftools_executable, "index", "--tbi", str(study_vcf_path)],
            timeout,
            "bcftools_index_shapeit5_study",
        )
        if study_index_path.is_symlink() or not study_index_path.is_file():
            raise Shapeit5InputExternalError("shapeit5_study_index_missing")
        expected_samples = [row["SAMPLE_ID"] for row in sample_rows]
        observed_samples = _bcftools_samples(
            bcftools_executable,
            study_vcf_path,
            command_runner,
            timeout,
        )
        if observed_samples != expected_samples:
            raise Shapeit5InputBlockError("shapeit5_study_sample_set_or_order_mismatch")
        expected_variants = [
            (
                f"chr{chromosome}",
                int(row["POSITION_BP"]),
                row["VARIANT_ID"],
                row["CANONICAL_REF"],
                row["CANONICAL_ALT"],
            )
            for row in included_rows
        ]
        observed_variants = _bcftools_variants(
            bcftools_executable,
            study_vcf_path,
            command_runner,
            timeout,
        )
        if observed_variants != expected_variants:
            raise Shapeit5InputBlockError("shapeit5_study_variant_identity_or_order_mismatch")
        reference_samples = _bcftools_samples(
            bcftools_executable,
            reference_vcf_path,
            command_runner,
            timeout,
        )
        if set(reference_samples) & set(observed_samples):
            raise Shapeit5InputBlockError("study_reference_sample_id_overlap")
        reference_variants = set(
            _bcftools_variants(
                bcftools_executable,
                reference_vcf_path,
                command_runner,
                timeout,
            )
        )
        expected_common = {
            (
                f"chr{chromosome}",
                int(row["POSITION_BP"]),
                row["REFERENCE_ID"] or ".",
                row["REFERENCE_REF"],
                row["REFERENCE_ALT"],
            )
            for row in harmonization_table.rows
            if row["STUDY_VARIANT_ID"]
            in {common_row["VARIANT_ID"] for common_row in common_rows}
        }
        reference_keys = {
            (chromosome_name, position_bp, ref, alt)
            for chromosome_name, position_bp, _, ref, alt in reference_variants
        }
        if any(
            (chromosome_name, position_bp, ref, alt) not in reference_keys
            for chromosome_name, position_bp, _, ref, alt in expected_common
        ):
            raise Shapeit5InputBlockError("shapeit5_common_variant_missing_from_reference")
        if not (
            int(map_table.rows[0]["POSITION_BP"])
            <= int(included_rows[0]["POSITION_BP"])
            and int(map_table.rows[-1]["POSITION_BP"])
            >= int(included_rows[-1]["POSITION_BP"])
        ):
            raise Shapeit5InputBlockError("shapeit5_map_does_not_cover_study_variants")

        for temporary_path in (
            extract_path,
            reference_alleles_path,
            update_ids_path,
            Path(f"{study_prefix}.log"),
            Path(f"{study_prefix}.nosex"),
        ):
            if temporary_path.exists():
                temporary_path.unlink()
        input_region = (
            f"chr{chromosome}:{included_rows[0]['POSITION_BP']}-{included_rows[-1]['POSITION_BP']}"
        )
        scaffold_region = input_region
        manifest = {
            "schema_version": "1.0.0",
            "created_at": utc_now(),
            "method_id": "plink_vcf_shapeit5_inputs_v1",
            "adapter_id": SHAPEIT5_CONTRACT.adapter_id,
            "assembly": "GRCh38",
            "chromosome": chromosome,
            "input_region": input_region,
            "scaffold_region": scaffold_region,
            "sample_set_id": phasing_manifest["sample_set_id"],
            "study_variant_set_id": f"shapeit5_study_{sha256_file(study_vcf_path)[:16]}",
            "reference_variant_set_id": reference_variant_set_id,
            "sample_count": len(observed_samples),
            "study_variant_count": len(observed_variants),
            "common_variant_count": len(common_rows),
            "pedigree_record_count": len(pedigree_rows),
            "target_variant_id": target_metadata["project_variant_id"],
            "target_variant_role": target_row["SHAPEIT_ROLE"],
            "tools": {
                "plink": plink_version,
                "bcftools": bcftools_version.removeprefix("bcftools "),
            },
            "files": {
                "study_vcf": _file_record(study_vcf_path),
                "study_index": _file_record(study_index_path),
                "reference_vcf": _stage_file_record(reference_vcf_path, stage_root),
                "reference_index": _stage_file_record(reference_index_path, stage_root),
                "genetic_map": _file_record(genetic_map_path),
                "pedigree": _file_record(pedigree_path),
                "variant_selection": _file_record(selection_path),
                "sample_mapping": _file_record(sample_mapping_path),
            },
            "checks": {
                "sample_identity_and_order": "PASS",
                "pedigree_identity": "PASS",
                "canonical_ref_alt": "PASS",
                "reference_common_variants": "PASS",
                "target_variant_retained": "PASS",
                "genetic_map": "PASS",
                "tabix_index": "PASS",
            },
        }
        validate_json_document(manifest, "shapeit5_inputs_manifest.schema.json")
        atomic_write_json(staging_dir / "shapeit5_inputs_manifest.json", manifest)
        os.replace(staging_dir, output_dir)
        published = True
        return _result(output_dir, manifest)
    except (DocumentValidationError, TableValidationError, OSError) as error:
        raise Shapeit5InputError("shapeit5_inputs_publication_failed") from error
    finally:
        if not published:
            shutil.rmtree(staging_dir, ignore_errors=True)
