"""Étape 04 : injection auditée du variant cible dans le jeu PLINK dédié."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


GENOTYPE_AUDIT_COLUMNS = (
    "SAMPLE_ID",
    "FID",
    "IID",
    "GENOTYPE",
    "GENOTYPE_SOURCE",
    "MEASUREMENT_METHOD",
    "STATUS",
)
FORBIDDEN_GENOTYPE_SOURCES = frozenset(
    {
        "clinical",
        "clinical_status",
        "group",
        "group_label",
        "inferred_from_clinical_status",
        "inferred_from_group",
    }
)


class TargetVariantInputError(ValueError):
    """Signale une entrée technique invalide ou incohérente."""


class ExternalToolError(RuntimeError):
    """Signale un échec de PLINK."""


class TargetVariantBlockError(RuntimeError):
    """Signale une incompatibilité scientifique qui interdit l'injection."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="prepare-target-variant-dataset")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _resolve_input_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    physical_path = _resolve_input_path(artifact, run_dir)
    if not physical_path.is_file() or sha256_file(physical_path) != artifact["sha256"]:
        raise TargetVariantInputError("declared_input_modified")
    return physical_path


def _artifact_by_id(
    stage_inputs: dict[str, Any], artifact_id: str
) -> dict[str, Any]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise TargetVariantInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _load_variant_metadata(path: Path) -> dict[str, Any]:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        raise TargetVariantInputError("invalid_target_variant_metadata_yaml") from error
    if not isinstance(document, dict):
        raise TargetVariantInputError("invalid_target_variant_metadata_document")
    validate_json_document(document, "target_variant_metadata.schema.json")
    if document["ref"] == document["alt"]:
        raise TargetVariantBlockError("target_ref_equals_alt")
    return document


def _validate_metadata_against_config(
    metadata: dict[str, Any], config: dict[str, Any]
) -> None:
    expected_values = {
        "assembly": config["project"]["assembly"],
        "gene": config["target"]["gene"],
        "chromosome": config["target"]["chromosome"],
        "position_bp": config["target"]["position_bp"],
        "ref": config["target"]["ref"],
        "alt": config["target"]["alt"],
        "transcript": config["target"]["transcript"],
        "project_variant_id": config["target"]["project_variant_id"],
        "rsid": config["target"]["rsid"],
    }
    missing_config = [
        key
        for key, value in expected_values.items()
        if value is None and key != "rsid"
    ]
    if missing_config:
        raise TargetVariantBlockError(
            f"incomplete_target_configuration:{','.join(sorted(missing_config))}"
        )
    mismatches = [
        key for key, expected in expected_values.items() if metadata[key] != expected
    ]
    if mismatches:
        raise TargetVariantBlockError(
            f"target_metadata_config_mismatch:{','.join(sorted(mismatches))}"
        )


def _normalize_genotype(genotype: str, ref: str, alt: str) -> tuple[str, str, str]:
    normalized = genotype.strip().upper().replace("|", "/")
    alleles = normalized.split("/")
    if len(alleles) != 2 or any(allele not in {ref, alt} for allele in alleles):
        raise TargetVariantBlockError("genotype_outside_declared_ref_alt")
    ordered = sorted(alleles, key=lambda allele: (allele != ref, allele))
    return ordered[0], ordered[1], "/".join(ordered)


def _validated_genotype_rows(
    genotype_path: Path,
    sample_rows: tuple[dict[str, Any], ...],
    ref: str,
    alt: str,
) -> list[dict[str, str]]:
    genotype_table = validate_tsv_table(
        genotype_path, "target_variant_genotypes.schema.json"
    )
    genotypes_by_sample = {row["SAMPLE_ID"]: row for row in genotype_table.rows}
    expected_ids = {row["SAMPLE_ID"] for row in sample_rows}
    observed_ids = set(genotypes_by_sample)
    if expected_ids != observed_ids:
        missing_count = len(expected_ids - observed_ids)
        extra_count = len(observed_ids - expected_ids)
        raise TargetVariantBlockError(
            "target_genotype_sample_set_mismatch:"
            f"missing={missing_count},extra={extra_count}"
        )

    validated_rows: list[dict[str, str]] = []
    for sample_row in sample_rows:
        genotype_row = genotypes_by_sample[sample_row["SAMPLE_ID"]]
        source = str(genotype_row["GENOTYPE_SOURCE"])
        if source.casefold() in FORBIDDEN_GENOTYPE_SOURCES:
            raise TargetVariantBlockError("forbidden_genotype_inference_source")
        allele_1, allele_2, normalized = _normalize_genotype(
            str(genotype_row["GENOTYPE"]), ref, alt
        )
        registry_genotype = sample_row["TARGET_GENOTYPE"]
        registry_source = sample_row["TARGET_GENOTYPE_SOURCE"]
        if registry_genotype is not None:
            _, _, normalized_registry = _normalize_genotype(
                str(registry_genotype), ref, alt
            )
            if normalized_registry != normalized or registry_source != source:
                raise TargetVariantBlockError("sample_registry_genotype_mismatch")
        validated_rows.append(
            {
                "SAMPLE_ID": sample_row["SAMPLE_ID"],
                "FID": sample_row["FID"],
                "IID": sample_row["IID"],
                "GENOTYPE": normalized,
                "GENOTYPE_SOURCE": source,
                "MEASUREMENT_METHOD": str(genotype_row["MEASUREMENT_METHOD"]),
                "STATUS": "ACCEPTED",
                "ALLELE_1": allele_1,
                "ALLELE_2": allele_2,
            }
        )
    return validated_rows


def _read_fam(path: Path) -> list[list[str]]:
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    if not rows or any(len(row) != 6 for row in rows):
        raise TargetVariantInputError("malformed_target_base_fam")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    if not rows or any(len(row) != 6 for row in rows):
        raise TargetVariantInputError("malformed_target_base_bim")
    return rows


def _validate_base_dataset(
    descriptor: dict[str, Any],
    descriptor_path: Path,
    bed_path: Path,
    bim_path: Path,
    fam_path: Path,
    sample_rows: tuple[dict[str, Any], ...],
    metadata: dict[str, Any],
) -> tuple[list[list[str]], list[list[str]]]:
    validate_json_document(descriptor, "plink_dataset.schema.json")
    expected_files = {"bed": bed_path, "bim": bim_path, "fam": fam_path}
    for suffix, physical_path in expected_files.items():
        declared = descriptor["files"][suffix]
        if declared["name"] != physical_path.name or declared["sha256"] != sha256_file(
            physical_path
        ):
            raise TargetVariantInputError("base_dataset_descriptor_mismatch")
    if bed_path.read_bytes()[:3] != b"\x6c\x1b\x01":
        raise TargetVariantInputError("invalid_target_base_bed_header")
    if descriptor["assembly"] != metadata["assembly"]:
        raise TargetVariantBlockError("target_base_assembly_mismatch")
    if descriptor["target_chromosome"] != metadata["chromosome"]:
        raise TargetVariantBlockError("target_base_chromosome_mismatch")

    fam_rows = _read_fam(fam_path)
    bim_rows = _read_bim(bim_path)
    expected_sample_pairs = [(row["FID"], row["IID"]) for row in sample_rows]
    observed_sample_pairs = [(row[0], row[1]) for row in fam_rows]
    if observed_sample_pairs != expected_sample_pairs:
        raise TargetVariantBlockError("target_base_sample_order_mismatch")
    if (
        descriptor["sample_count"] != len(fam_rows)
        or descriptor["variant_count"] != len(bim_rows)
    ):
        raise TargetVariantInputError("target_base_descriptor_count_mismatch")
    variant_ids = [row[1] for row in bim_rows]
    if metadata["project_variant_id"] in variant_ids:
        raise TargetVariantBlockError("target_variant_already_present")
    if any(
        row[0] == str(metadata["chromosome"])
        and int(row[3]) == metadata["position_bp"]
        for row in bim_rows
    ):
        raise TargetVariantBlockError("target_position_already_present")
    return fam_rows, bim_rows


def _write_single_variant_ped_map(
    prefix: Path,
    fam_rows: list[list[str]],
    genotype_rows: list[dict[str, str]],
    metadata: dict[str, Any],
) -> None:
    with prefix.with_suffix(".map").open("w", encoding="utf-8", newline="") as map_file:
        writer = csv.writer(map_file, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                metadata["chromosome"],
                metadata["project_variant_id"],
                "0",
                metadata["position_bp"],
            )
        )
    with prefix.with_suffix(".ped").open("w", encoding="utf-8", newline="") as ped_file:
        writer = csv.writer(ped_file, delimiter="\t", lineterminator="\n")
        for fam_row, genotype_row in zip(fam_rows, genotype_rows, strict=True):
            writer.writerow(
                fam_row + [genotype_row["ALLELE_1"], genotype_row["ALLELE_2"]]
            )


def _resolve_executable(command: str | None) -> str:
    if command is None:
        raise ExternalToolError("plink_not_configured")
    resolved = shutil.which(command)
    if resolved is None:
        raise ExternalToolError("plink_not_available")
    return resolved


def _plink_version(executable: str) -> str | None:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            check=False,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    output = completed.stdout.strip() or completed.stderr.strip()
    return output.splitlines()[0][:500] if output else None


def _run_plink(executable: str, arguments: list[str], timeout_seconds: int) -> None:
    try:
        completed = subprocess.run(
            [executable, *arguments],
            capture_output=True,
            check=False,
            text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExternalToolError("plink_execution_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError("plink_command_failed")


def _complete_single_variant_alleles(prefix: Path, metadata: dict[str, Any]) -> None:
    bim_rows = _read_bim(prefix.with_suffix(".bim"))
    if len(bim_rows) != 1 or bim_rows[0][1] != metadata["project_variant_id"]:
        raise ExternalToolError("invalid_single_variant_bim")
    row = bim_rows[0]
    declared = {metadata["ref"], metadata["alt"]}
    observed = {allele for allele in row[4:6] if allele != "0"}
    if not observed or not observed <= declared:
        raise ExternalToolError("single_variant_allele_mismatch")
    if observed != declared:
        missing_allele = (declared - observed).pop()
        missing_index = 4 if row[4] == "0" else 5
        row[missing_index] = missing_allele
        prefix.with_suffix(".bim").write_text("\t".join(row) + "\n", encoding="utf-8")


def _validate_enriched_dataset(
    prefix: Path,
    sample_count: int,
    variant_count: int,
    target_variant_id: str,
    expected_fam_rows: list[list[str]],
    ref: str,
    alt: str,
) -> list[list[str]]:
    try:
        bed_header = prefix.with_suffix(".bed").read_bytes()[:3]
        fam_rows = _read_fam(prefix.with_suffix(".fam"))
        bim_rows = _read_bim(prefix.with_suffix(".bim"))
    except OSError as error:
        raise ExternalToolError("missing_enriched_plink_output") from error
    if bed_header != b"\x6c\x1b\x01":
        raise ExternalToolError("invalid_enriched_bed_header")
    if len(fam_rows) != sample_count or fam_rows != expected_fam_rows:
        raise ExternalToolError("unexpected_enriched_sample_set")
    if len(bim_rows) != variant_count:
        raise ExternalToolError("unexpected_enriched_variant_count")
    target_rows = [row for row in bim_rows if row[1] == target_variant_id]
    if len(target_rows) != 1 or set(target_rows[0][4:6]) != {ref, alt}:
        raise ExternalToolError("invalid_enriched_target_variant")
    return bim_rows


def _write_tsv(
    path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _evaluable_trio_count(fam_rows: list[list[str]]) -> int:
    sample_pairs = {(row[0], row[1]) for row in fam_rows}
    return sum(
        row[2] != "0"
        and row[3] != "0"
        and (row[0], row[2]) in sample_pairs
        and (row[0], row[3]) in sample_pairs
        for row in fam_rows
    )


def _mendel_error_count(
    path: Path, target_variant_id: str, evaluable_trio_count: int
) -> int:
    if not path.is_file() and evaluable_trio_count == 0:
        path.write_text("FID KID CHR SNP CODE ERROR\n", encoding="utf-8")
    try:
        lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines()]
    except OSError as error:
        raise ExternalToolError("missing_plink_mendel_report") from error
    rows = [row for row in lines[1:] if row]
    if any(target_variant_id not in row for row in rows):
        raise ExternalToolError("unexpected_variant_in_mendel_report")
    return len(rows)


def _variant_set_id(parent_variant_set_id: str, metadata: dict[str, Any]) -> str:
    digest = hashlib.sha256(
        (
            f"{parent_variant_set_id}\n{metadata['chromosome']}:"
            f"{metadata['position_bp']}:{metadata['project_variant_id']}:"
            f"{metadata['ref']}:{metadata['alt']}\n"
        ).encode("utf-8")
    ).hexdigest()
    return f"target_enriched_{digest[:12]}"


def _write_dataset_descriptor(
    path: Path,
    prefix: Path,
    base_descriptor: dict[str, Any],
    variant_set_id: str,
) -> dict[str, Any]:
    descriptor = {
        **base_descriptor,
        "dataset_id": "target_variant_dataset",
        "variant_set_id": variant_set_id,
        "variant_count": base_descriptor["variant_count"] + 1,
        "source_format": "ChAS_ACPA_PLUS_EXPLICIT_TARGET",
        "allele_orientation": "mixed_forward_and_reference",
        "files": {
            suffix: {
                "name": prefix.with_suffix(f".{suffix}").name,
                "sha256": sha256_file(prefix.with_suffix(f".{suffix}")),
            }
            for suffix in ("bed", "bim", "fam")
        },
    }
    validate_json_document(descriptor, "plink_dataset.schema.json")
    atomic_write_json(path, descriptor)
    return descriptor


def _dataset_artifacts(
    prefix: Path,
    descriptor_path: Path,
    published_dir: PurePosixPath,
    producer_stage: str,
    signature: str,
    assembly: str,
    sample_set_id: str,
    variant_set_id: str,
) -> list[dict[str, Any]]:
    artifacts = []
    media_types = {
        "bed": "application/octet-stream",
        "bim": "text/plain",
        "fam": "text/plain",
    }
    for suffix in ("bed", "bim", "fam"):
        physical_path = prefix.with_suffix(f".{suffix}")
        artifacts.append(
            build_file_artifact(
                physical_path=physical_path,
                published_path=str(published_dir / physical_path.name),
                artifact_id=f"target_variant_{suffix}",
                artifact_type=f"plink_{suffix}",
                media_type=media_types[suffix],
                producer_stage=producer_stage,
                producer_signature=signature,
                assembly=assembly,
                sample_set_id=sample_set_id,
                variant_set_id=variant_set_id,
                sensitivity="sensitive_genetic",
            )
        )
    artifacts.append(
        build_file_artifact(
            physical_path=descriptor_path,
            published_path=str(published_dir / descriptor_path.name),
            artifact_id="target_variant_dataset",
            artifact_type="plink_dataset_descriptor",
            media_type="application/json",
            producer_stage=producer_stage,
            producer_signature=signature,
            schema_name="plink_dataset.schema.json",
            schema_version="1.0.0",
            assembly=assembly,
            sample_set_id=sample_set_id,
            variant_set_id=variant_set_id,
            sensitivity="sensitive_genetic",
        )
    )
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Prépare le jeu cible enrichi sans modifier le jeu PLINK de l'étape 03."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    if config["policies"]["allow_clinical_genotype_inference"]:
        raise TargetVariantBlockError("clinical_genotype_inference_policy_forbidden")
    timeout_seconds = int(stage_inputs["parameters"].get("plink_timeout_seconds", 300))
    if timeout_seconds <= 0:
        raise TargetVariantInputError("invalid_plink_timeout")

    artifact_ids = (
        "config_input_target_variant_metadata",
        "config_input_target_variant_genotypes",
        "samples_master",
        "target_chromosome_base_bed",
        "target_chromosome_base_bim",
        "target_chromosome_base_fam",
        "target_chromosome_base_dataset",
    )
    artifacts = {
        artifact_id: _artifact_by_id(stage_inputs, artifact_id)
        for artifact_id in artifact_ids
    }
    paths = {
        artifact_id: _validated_artifact_path(artifact, run_dir)
        for artifact_id, artifact in artifacts.items()
    }
    metadata = _load_variant_metadata(paths["config_input_target_variant_metadata"])
    _validate_metadata_against_config(metadata, config)
    samples_table = validate_tsv_table(
        paths["samples_master"], "samples_master.schema.json"
    )
    genotype_rows = _validated_genotype_rows(
        paths["config_input_target_variant_genotypes"],
        samples_table.rows,
        metadata["ref"],
        metadata["alt"],
    )
    base_descriptor = read_json(paths["target_chromosome_base_dataset"])
    fam_rows, base_bim_rows = _validate_base_dataset(
        base_descriptor,
        paths["target_chromosome_base_dataset"],
        paths["target_chromosome_base_bed"],
        paths["target_chromosome_base_bim"],
        paths["target_chromosome_base_fam"],
        samples_table.rows,
        metadata,
    )
    if artifacts["samples_master"]["sample_set_id"] != base_descriptor["sample_set_id"]:
        raise TargetVariantBlockError("target_base_sample_set_id_mismatch")

    output_dir.mkdir(parents=True, exist_ok=True)
    genotype_audit_path = output_dir / "target_genotype_audit.tsv"
    _write_tsv(genotype_audit_path, GENOTYPE_AUDIT_COLUMNS, genotype_rows)
    validate_tsv_table(genotype_audit_path, "target_genotype_audit.schema.json")

    single_prefix = output_dir / "target_variant_single"
    enriched_prefix = output_dir / "target_variant"
    _write_single_variant_ped_map(single_prefix, fam_rows, genotype_rows, metadata)
    plink_executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(plink_executable)
    _run_plink(
        plink_executable,
        ["--file", str(single_prefix), "--make-bed", "--out", str(single_prefix)],
        timeout_seconds,
    )
    _complete_single_variant_alleles(single_prefix, metadata)
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(paths["target_chromosome_base_bed"].with_suffix("")),
            "--bmerge",
            str(single_prefix.with_suffix(".bed")),
            str(single_prefix.with_suffix(".bim")),
            str(single_prefix.with_suffix(".fam")),
            "--make-bed",
            "--out",
            str(enriched_prefix),
        ],
        timeout_seconds,
    )
    enriched_bim_rows = _validate_enriched_dataset(
        enriched_prefix,
        samples_table.row_count,
        len(base_bim_rows) + 1,
        metadata["project_variant_id"],
        fam_rows,
        metadata["ref"],
        metadata["alt"],
    )
    mendel_prefix = output_dir / "target_variant_mendel"
    _run_plink(
        plink_executable,
        [
            "--bfile",
            str(enriched_prefix),
            "--snp",
            metadata["project_variant_id"],
            "--mendel",
            "--out",
            str(mendel_prefix),
        ],
        timeout_seconds,
    )
    mendel_path = mendel_prefix.with_suffix(".mendel")
    evaluable_trio_count = _evaluable_trio_count(fam_rows)
    mendel_error_count = _mendel_error_count(
        mendel_path, metadata["project_variant_id"], evaluable_trio_count
    )
    mendel_summary_path = output_dir / "target_variant_mendel.tsv"
    if mendel_error_count:
        mendel_status = "FAIL"
    elif evaluable_trio_count:
        mendel_status = "PASS"
    else:
        mendel_status = "NOT_APPLICABLE"
    _write_tsv(
        mendel_summary_path,
        ("VARIANT_ID", "ERROR_COUNT", "STATUS"),
        [
            {
                "VARIANT_ID": metadata["project_variant_id"],
                "ERROR_COUNT": str(mendel_error_count),
                "STATUS": mendel_status,
            }
        ],
    )
    validate_tsv_table(mendel_summary_path, "target_variant_mendel.schema.json")

    before_checksums = {
        suffix: sha256_file(paths[f"target_chromosome_base_{suffix}"])
        for suffix in ("bed", "bim", "fam")
    }
    after_checksums = {
        suffix: sha256_file(enriched_prefix.with_suffix(f".{suffix}"))
        for suffix in ("bed", "bim", "fam")
    }
    report_path = output_dir / "target_variant_injection_report.json"
    report = {
        "schema_version": "1.0.0",
        "variant": metadata,
        "sample_count": samples_table.row_count,
        "variant_count_before": len(base_bim_rows),
        "variant_count_after": len(enriched_bim_rows),
        "mendel_error_count": mendel_error_count,
        "evaluable_trio_count": evaluable_trio_count,
        "mendel_status": mendel_status,
        "input_checksums": before_checksums,
        "output_checksums": after_checksums,
    }
    atomic_write_json(report_path, report)
    if mendel_error_count:
        raise TargetVariantBlockError("target_variant_mendel_errors")

    variant_set_id = _variant_set_id(base_descriptor["variant_set_id"], metadata)
    descriptor_path = output_dir / "target_variant.dataset.json"
    _write_dataset_descriptor(
        descriptor_path, enriched_prefix, base_descriptor, variant_set_id
    )
    temporary_suffixes = (
        "ped",
        "map",
        "log",
        "nosex",
        "hh",
        "lmendel",
        "fmendel",
        "imendel",
    )
    for prefix in (single_prefix, enriched_prefix, mendel_prefix):
        for suffix in temporary_suffixes:
            generated_path = prefix.with_suffix(f".{suffix}")
            if generated_path.exists():
                generated_path.unlink()
    for suffix in ("bed", "bim", "fam"):
        single_path = single_prefix.with_suffix(f".{suffix}")
        if single_path.exists():
            single_path.unlink()

    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    sample_set_id = base_descriptor["sample_set_id"]
    artifacts_out = _dataset_artifacts(
        enriched_prefix,
        descriptor_path,
        published_dir,
        producer_stage,
        stage_inputs["signature"],
        metadata["assembly"],
        sample_set_id,
        variant_set_id,
    )
    for artifact_id, artifact_type, path, schema_name in (
        (
            "target_genotype_audit",
            "target_genotype_audit",
            genotype_audit_path,
            "target_genotype_audit.schema.json",
        ),
        (
            "target_variant_mendel",
            "mendel_summary",
            mendel_summary_path,
            "target_variant_mendel.schema.json",
        ),
        (
            "target_variant_mendel_report",
            "plink_mendel_report",
            mendel_path,
            None,
        ),
        (
            "target_variant_injection_report",
            "target_variant_injection_report",
            report_path,
            None,
        ),
    ):
        artifacts_out.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type=(
                    "application/json"
                    if path.suffix == ".json"
                    else "text/tab-separated-values"
                ),
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=metadata["assembly"],
                sample_set_id=sample_set_id,
                variant_set_id=variant_set_id,
                sensitivity="sensitive_genetic",
            )
        )

    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": artifacts_out,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "explicit_target_variant_injection_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(artifacts.values()),
        "outputs": artifacts_out,
        "parameters": stage_inputs["parameters"],
        "tools": [
            {
                "tool": "plink",
                "configured": config["tools"]["plink"],
                "version": plink_version,
            }
        ],
        "counts": {
            "samples": samples_table.row_count,
            "variants_before": len(base_bim_rows),
            "variants_after": len(enriched_bim_rows),
            "mendel_errors": mendel_error_count,
            "evaluable_trios": evaluable_trio_count,
        },
        "metrics": {},
        "exclusions": [],
        "warnings": [],
        "checks": [
            {"check": "molecular_definition_confirmed", "status": "PASS"},
            {"check": "explicit_genotype_for_every_sample", "status": "PASS"},
            {"check": "clinical_inference_forbidden", "status": "PASS"},
            {"check": "base_dataset_unchanged", "status": "PASS"},
            {"check": "mendel_target_variant", "status": mendel_status},
        ],
        "known_limits": [
            "La confirmation REF/ALT repose sur la source externe explicitement "
            "auditée dans reference_validation.",
            "Le contrôle mendélien repose sur les liens PID/MID enregistrés dans "
            "le registre approuvé.",
        ],
        "expected_visualizations": [],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in artifacts_out
        ),
        encoding="utf-8",
    )
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except ExternalToolError as error:
        sys.stderr.write(f"{error}\n")
        return 3
    except TargetVariantBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        TargetVariantInputError,
        DocumentValidationError,
        TableValidationError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
