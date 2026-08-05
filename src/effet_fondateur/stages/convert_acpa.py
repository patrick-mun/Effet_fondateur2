"""Étape 03 : conversion ACPA en jeux PLINK genome-wide et chromosome cible."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from itertools import chain
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

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


REQUIRED_ACPA_COLUMNS = frozenset(
    {
        "Probe Set ID",
        "Forward Strand Base Calls",
        "dbSNP RS ID",
        "Chromosome",
        "Chromosomal Position",
    }
)
VALID_BASES = frozenset("ACGT")
AUTOSOMES = frozenset(str(chromosome) for chromosome in range(1, 23))
ALIGNMENT_COLUMNS = (
    "ROW_INDEX",
    "SAMPLE_ID",
    "FID",
    "IID",
    "SOURCE_FILE",
    "INVALID_GENOTYPE_COUNT",
    "SKIPPED_ROW_COUNT",
)
VARIANT_AUDIT_COLUMNS = (
    "VARIANT_ID",
    "PROBE_ID",
    "RSID",
    "RSID_STATUS",
    "CHROMOSOME",
    "POSITION_BP",
    "OBSERVED_ALLELES",
    "PRESENCE_COUNT",
    "STATUS",
    "EXCLUSION_CODE",
    "IN_GENOMEWIDE",
    "IN_TARGET_CHROMOSOME",
)


class AcpaConversionError(ValueError):
    """Signale une source ou un paramètre incompatible avec la conversion."""


class ExternalToolError(RuntimeError):
    """Signale l'échec de PLINK avec le code V2 réservé aux outils externes."""


class ScientificBlockError(RuntimeError):
    """Signale l'absence d'un jeu de variants exploitable après contrôles."""


@dataclass
class MarkerAggregate:
    """Métadonnées agrégées d'une sonde, sans conserver tous les génotypes."""

    probe_id: str
    chromosome: str
    position_bp: int
    presence_count: int = 0
    observed_alleles: set[str] = field(default_factory=set)
    rsids: set[str] = field(default_factory=set)
    coordinate_conflict: bool = False


@dataclass(frozen=True)
class SampleSpool:
    """Appels temporaires et métriques d'un export lu exactement une fois."""

    row: dict[str, Any]
    source_artifact: dict[str, Any]
    spool_path: Path
    invalid_genotype_count: int
    skipped_row_count: int
    chromosome_marker_counts: dict[str, int]
    chromosome_called_counts: dict[str, int]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="convert-acpa")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _normalize_chromosome(value: str) -> str:
    chromosome = value.strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    return chromosome.upper()


def _normalized_build(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = value.strip().lower().replace("grch", "hg")
    return normalized or None


def _parse_genotype(value: str) -> tuple[str, str]:
    genotype = "".join(value.strip().upper().split())
    if len(genotype) != 2 or any(base not in VALID_BASES for base in genotype):
        return "0", "0"
    return genotype[0], genotype[1]


def _read_source_to_spool(
    source_path: Path,
    spool_path: Path,
    expected_assembly: str,
    markers: dict[str, MarkerAggregate],
) -> tuple[int, int, dict[str, int], dict[str, int]]:
    """Lit une source une fois, écrit ses appels et agrège les sondes autosomiques."""
    invalid_genotype_count = 0
    skipped_row_count = 0
    metadata: dict[str, str] = {}
    seen_probes: set[str] = set()
    chromosome_marker_counts: dict[str, int] = {}
    chromosome_called_counts: dict[str, int] = {}
    with source_path.open("r", encoding="utf-8-sig", newline="") as source_file:
        header_line: str | None = None
        for raw_line in source_file:
            if raw_line.startswith("#"):
                key, separator, value = raw_line[1:].strip().partition(":")
                if separator:
                    metadata[key.strip()] = value.strip()
                continue
            if raw_line.strip():
                header_line = raw_line
                break
        if header_line is None:
            raise AcpaConversionError("missing_acpa_header")
        if _normalized_build(metadata.get("UCSC Genomic Version")) != _normalized_build(
            expected_assembly
        ):
            raise AcpaConversionError("incompatible_acpa_assembly")

        reader = csv.DictReader(chain([header_line], source_file), delimiter="\t")
        fieldnames = [field.strip() for field in (reader.fieldnames or [])]
        reader.fieldnames = fieldnames
        if REQUIRED_ACPA_COLUMNS - set(fieldnames):
            raise AcpaConversionError("missing_acpa_columns")
        with spool_path.open("w", encoding="utf-8", newline="") as spool_file:
            writer = csv.writer(spool_file, delimiter="\t", lineterminator="\n")
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise AcpaConversionError("malformed_acpa_row")
                probe_id = str(row["Probe Set ID"]).strip()
                if not probe_id:
                    skipped_row_count += 1
                    continue
                if any(character.isspace() for character in probe_id):
                    raise AcpaConversionError("probe_id_contains_whitespace")
                if probe_id in seen_probes:
                    raise AcpaConversionError("duplicate_probe_in_source")
                seen_probes.add(probe_id)
                chromosome = _normalize_chromosome(str(row["Chromosome"]))
                if chromosome not in AUTOSOMES:
                    continue
                try:
                    position_bp = int(float(str(row["Chromosomal Position"]).strip()))
                except ValueError:
                    skipped_row_count += 1
                    continue
                if position_bp <= 0:
                    skipped_row_count += 1
                    continue
                allele_1, allele_2 = _parse_genotype(
                    str(row["Forward Strand Base Calls"])
                )
                if allele_1 == "0":
                    invalid_genotype_count += 1
                chromosome_marker_counts[chromosome] = (
                    chromosome_marker_counts.get(chromosome, 0) + 1
                )
                if allele_1 != "0":
                    chromosome_called_counts[chromosome] = (
                        chromosome_called_counts.get(chromosome, 0) + 1
                    )
                rsid = str(row["dbSNP RS ID"]).strip()
                if rsid.casefold() in {"", ".", "na", "nan", "none"}:
                    rsid = ""
                writer.writerow((probe_id, allele_1, allele_2))

                marker = markers.get(probe_id)
                if marker is None:
                    marker = MarkerAggregate(probe_id, chromosome, position_bp)
                    markers[probe_id] = marker
                elif (
                    marker.chromosome != chromosome
                    or marker.position_bp != position_bp
                ):
                    marker.coordinate_conflict = True
                marker.presence_count += 1
                marker.observed_alleles.update(
                    allele for allele in (allele_1, allele_2) if allele != "0"
                )
                if rsid:
                    marker.rsids.add(rsid)
    return (
        invalid_genotype_count,
        skipped_row_count,
        chromosome_marker_counts,
        chromosome_called_counts,
    )


def _resolve_input_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    physical_path = _resolve_input_path(artifact, run_dir)
    if (
        not physical_path.is_file()
        or sha256_file(physical_path) != artifact["sha256"]
    ):
        raise AcpaConversionError("declared_input_modified")
    return physical_path


def _source_artifacts_by_relative_path(
    stage_inputs: dict[str, Any],
    source_root: Path,
    run_dir: Path,
) -> dict[str, tuple[dict[str, Any], Path]]:
    indexed: dict[str, tuple[dict[str, Any], Path]] = {}
    for artifact in stage_inputs["artifacts"]:
        if artifact["artifact_type"] != "acpa_source_file":
            continue
        physical_path = _validated_artifact_path(artifact, run_dir)
        try:
            relative_path = physical_path.resolve().relative_to(source_root.resolve())
        except ValueError as error:
            raise AcpaConversionError("source_outside_configured_root") from error
        relative_name = relative_path.as_posix()
        if relative_name in indexed:
            raise AcpaConversionError("duplicate_declared_source_path")
        indexed[relative_name] = (artifact, physical_path)
    return indexed


def _select_markers(
    markers: dict[str, MarkerAggregate],
    sample_count: int,
    marker_mode: str,
    target_chromosome: str,
) -> tuple[list[MarkerAggregate], list[MarkerAggregate], list[dict[str, str]]]:
    genomewide: list[MarkerAggregate] = []
    target: list[MarkerAggregate] = []
    audit_rows: list[dict[str, str]] = []
    for marker in sorted(
        markers.values(),
        key=lambda item: (int(item.chromosome), item.position_bp, item.probe_id),
    ):
        exclusion_code = ""
        if marker.coordinate_conflict:
            exclusion_code = "coordinate_conflict"
        elif len(marker.observed_alleles) > 2:
            exclusion_code = "more_than_two_alleles"
        elif marker_mode == "intersection" and marker.presence_count != sample_count:
            exclusion_code = "missing_in_samples"
        accepted = not exclusion_code
        if accepted:
            genomewide.append(marker)
            if marker.chromosome == target_chromosome:
                target.append(marker)
        rsid = next(iter(marker.rsids)) if len(marker.rsids) == 1 else ""
        if len(marker.rsids) == 1:
            rsid_status = "resolved"
        elif marker.rsids:
            rsid_status = "conflicting"
        else:
            rsid_status = "unresolved"
        audit_rows.append(
            {
                "VARIANT_ID": marker.probe_id,
                "PROBE_ID": marker.probe_id,
                "RSID": rsid,
                "RSID_STATUS": rsid_status,
                "CHROMOSOME": marker.chromosome,
                "POSITION_BP": str(marker.position_bp),
                "OBSERVED_ALLELES": ",".join(sorted(marker.observed_alleles)),
                "PRESENCE_COUNT": str(marker.presence_count),
                "STATUS": "ACCEPTED" if accepted else "EXCLUDED",
                "EXCLUSION_CODE": exclusion_code,
                "IN_GENOMEWIDE": str(accepted).lower(),
                "IN_TARGET_CHROMOSOME": str(
                    accepted and marker.chromosome == target_chromosome
                ).lower(),
            }
        )
    return genomewide, target, audit_rows


def _load_spooled_calls(spool_path: Path) -> dict[str, tuple[str, str]]:
    calls: dict[str, tuple[str, str]] = {}
    with spool_path.open("r", encoding="utf-8", newline="") as spool_file:
        for probe_id, allele_1, allele_2 in csv.reader(spool_file, delimiter="\t"):
            calls[probe_id] = (allele_1, allele_2)
    return calls


def _pedigree_columns(row: dict[str, Any]) -> list[str]:
    sex = {"MALE": "1", "FEMALE": "2"}.get(row["SEX"], "0")
    phenotype = {
        "UNAFFECTED": "1",
        "AFFECTED": "2",
    }.get(row["CLINICAL_STATUS"], "-9")
    return [row["FID"], row["IID"], row["PID"], row["MID"], sex, phenotype]


def _write_map(path: Path, markers: list[MarkerAggregate]) -> None:
    with path.open("w", encoding="utf-8", newline="") as map_file:
        writer = csv.writer(map_file, delimiter="\t", lineterminator="\n")
        for marker in markers:
            writer.writerow(
                (marker.chromosome, marker.probe_id, "0", marker.position_bp)
            )


def _write_ped_datasets(
    genomewide_path: Path,
    target_path: Path,
    samples: list[SampleSpool],
    genomewide_markers: list[MarkerAggregate],
    target_markers: list[MarkerAggregate],
) -> None:
    with genomewide_path.open("w", encoding="utf-8", newline="") as genome_file, (
        target_path.open("w", encoding="utf-8", newline="")
    ) as target_file:
        for sample in samples:
            calls = _load_spooled_calls(sample.spool_path)
            pedigree = _pedigree_columns(sample.row)
            genome_genotypes = [
                allele
                for marker in genomewide_markers
                for allele in calls.get(marker.probe_id, ("0", "0"))
            ]
            target_genotypes = [
                allele
                for marker in target_markers
                for allele in calls.get(marker.probe_id, ("0", "0"))
            ]
            genome_file.write("\t".join(pedigree + genome_genotypes) + "\n")
            target_file.write("\t".join(pedigree + target_genotypes) + "\n")


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


def _run_plink(executable: str, prefix: Path, timeout_seconds: int) -> None:
    try:
        completed = subprocess.run(
            [
                executable,
                "--file",
                str(prefix),
                "--make-bed",
                "--out",
                str(prefix),
            ],
            capture_output=True,
            check=False,
            text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExternalToolError("plink_execution_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError("plink_make_bed_failed")


def _validate_plink_files(prefix: Path, sample_count: int, variant_count: int) -> None:
    bed_path = prefix.with_suffix(".bed")
    bim_path = prefix.with_suffix(".bim")
    fam_path = prefix.with_suffix(".fam")
    try:
        bed_header = bed_path.read_bytes()[:3]
        bim_rows = [
            line for line in bim_path.read_text(encoding="utf-8").splitlines() if line
        ]
        fam_rows = [
            line for line in fam_path.read_text(encoding="utf-8").splitlines() if line
        ]
    except OSError as error:
        raise ExternalToolError("missing_plink_output") from error
    if bed_header != b"\x6c\x1b\x01":
        raise ExternalToolError("invalid_plink_bed_header")
    if len(bim_rows) != variant_count:
        raise ExternalToolError("unexpected_plink_variant_count")
    if len(fam_rows) != sample_count:
        raise ExternalToolError("unexpected_plink_sample_count")
    if any(len(row.split()) != 6 for row in bim_rows + fam_rows):
        raise ExternalToolError("malformed_plink_text_file")


def _variant_set_id(scope: str, markers: list[MarkerAggregate]) -> str:
    digest = hashlib.sha256()
    for marker in markers:
        digest.update(
            f"{marker.chromosome}:{marker.position_bp}:{marker.probe_id}\n".encode(
                "utf-8"
            )
        )
    return f"{scope}_{digest.hexdigest()[:12]}"


def _write_dataset_descriptor(
    path: Path,
    prefix: Path,
    dataset_id: str,
    assembly: str,
    scope: str,
    target_chromosome: int | None,
    sample_set_id: str,
    variant_set_id: str,
    sample_count: int,
    variant_count: int,
    marker_mode: str,
) -> dict[str, Any]:
    descriptor = {
        "schema_version": "1.0.0",
        "dataset_id": dataset_id,
        "assembly": assembly,
        "scope": scope,
        "target_chromosome": target_chromosome,
        "sample_set_id": sample_set_id,
        "variant_set_id": variant_set_id,
        "sample_count": sample_count,
        "variant_count": variant_count,
        "marker_mode": marker_mode,
        "source_format": "ChAS_ACPA",
        "allele_orientation": "forward_strand_calls",
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


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _dataset_artifacts(
    prefix: Path,
    descriptor_path: Path,
    artifact_prefix: str,
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
                artifact_id=f"{artifact_prefix}_{suffix}",
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
            artifact_id=f"{artifact_prefix}_dataset",
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
    """Convertit les sources approuvées en deux jeux PLINK cohérents et audités."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    target_chromosome_value = config["target"]["chromosome"]
    if target_chromosome_value is None:
        raise AcpaConversionError("missing_target_chromosome")
    target_chromosome = str(target_chromosome_value)
    configured_source_root = config["inputs"]["acpa_samples_dir"]
    if configured_source_root is None:
        raise AcpaConversionError("missing_acpa_samples_dir")
    source_root = Path(configured_source_root)
    if not source_root.is_absolute():
        source_root = Path.cwd() / source_root

    parameters = stage_inputs["parameters"]
    marker_mode = str(parameters.get("marker_mode", "intersection"))
    if marker_mode not in {"intersection", "union"}:
        raise AcpaConversionError("invalid_marker_mode")
    timeout_seconds = int(parameters.get("plink_timeout_seconds", 300))
    if timeout_seconds <= 0:
        raise AcpaConversionError("invalid_plink_timeout")

    samples_artifacts = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == "samples_master"
    ]
    if len(samples_artifacts) != 1:
        raise AcpaConversionError("samples_master_missing_or_ambiguous")
    samples_artifact = samples_artifacts[0]
    samples_path = _validated_artifact_path(samples_artifact, run_dir)
    samples_table = validate_tsv_table(samples_path, "samples_master.schema.json")
    source_artifacts = _source_artifacts_by_relative_path(
        stage_inputs,
        source_root,
        run_dir,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    spool_dir = output_dir / ".spool"
    spool_dir.mkdir()
    markers: dict[str, MarkerAggregate] = {}
    sample_spools: list[SampleSpool] = []
    selected_input_artifacts = [samples_artifact]
    for row_index, row in enumerate(samples_table.rows, start=1):
        source_entry = source_artifacts.get(row["SOURCE_FILE"])
        if source_entry is None:
            raise AcpaConversionError("sample_source_not_declared")
        source_artifact, source_path = source_entry
        spool_path = spool_dir / f"sample_{row_index:06d}.tsv"
        (
            invalid_count,
            skipped_count,
            chromosome_marker_counts,
            chromosome_called_counts,
        ) = _read_source_to_spool(
            source_path,
            spool_path,
            config["project"]["assembly"],
            markers,
        )
        sample_spools.append(
            SampleSpool(
                row=row,
                source_artifact=source_artifact,
                spool_path=spool_path,
                invalid_genotype_count=invalid_count,
                skipped_row_count=skipped_count,
                chromosome_marker_counts=chromosome_marker_counts,
                chromosome_called_counts=chromosome_called_counts,
            )
        )
        selected_input_artifacts.append(source_artifact)

    genomewide_markers, target_markers, variant_audit_rows = _select_markers(
        markers,
        samples_table.row_count,
        marker_mode,
        target_chromosome,
    )
    if not genomewide_markers:
        raise ScientificBlockError("no_valid_genomewide_markers")
    if not target_markers:
        raise ScientificBlockError("no_valid_target_chromosome_markers")

    genomewide_prefix = output_dir / "genomewide_base"
    target_prefix = output_dir / "target_chromosome_base"
    _write_map(genomewide_prefix.with_suffix(".map"), genomewide_markers)
    _write_map(target_prefix.with_suffix(".map"), target_markers)
    _write_ped_datasets(
        genomewide_prefix.with_suffix(".ped"),
        target_prefix.with_suffix(".ped"),
        sample_spools,
        genomewide_markers,
        target_markers,
    )

    plink_executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(plink_executable)
    _run_plink(plink_executable, genomewide_prefix, timeout_seconds)
    _run_plink(plink_executable, target_prefix, timeout_seconds)
    _validate_plink_files(
        genomewide_prefix,
        samples_table.row_count,
        len(genomewide_markers),
    )
    _validate_plink_files(
        target_prefix,
        samples_table.row_count,
        len(target_markers),
    )

    sample_set_id = samples_artifact["sample_set_id"]
    if sample_set_id is None:
        raise AcpaConversionError("samples_master_missing_sample_set_id")
    genomewide_variant_set_id = _variant_set_id("genomewide", genomewide_markers)
    target_variant_set_id = _variant_set_id("target", target_markers)
    genomewide_descriptor_path = output_dir / "genomewide_base.dataset.json"
    target_descriptor_path = output_dir / "target_chromosome_base.dataset.json"
    _write_dataset_descriptor(
        genomewide_descriptor_path,
        genomewide_prefix,
        "genomewide_base",
        config["project"]["assembly"],
        "autosomal_genomewide",
        None,
        sample_set_id,
        genomewide_variant_set_id,
        samples_table.row_count,
        len(genomewide_markers),
        marker_mode,
    )
    _write_dataset_descriptor(
        target_descriptor_path,
        target_prefix,
        "target_chromosome_base",
        config["project"]["assembly"],
        "target_chromosome",
        target_chromosome_value,
        sample_set_id,
        target_variant_set_id,
        samples_table.row_count,
        len(target_markers),
        marker_mode,
    )

    alignment_path = output_dir / "sample_alignment.tsv"
    alignment_rows = [
        {
            "ROW_INDEX": str(row_index),
            "SAMPLE_ID": sample.row["SAMPLE_ID"],
            "FID": sample.row["FID"],
            "IID": sample.row["IID"],
            "SOURCE_FILE": sample.row["SOURCE_FILE"],
            "INVALID_GENOTYPE_COUNT": str(sample.invalid_genotype_count),
            "SKIPPED_ROW_COUNT": str(sample.skipped_row_count),
        }
        for row_index, sample in enumerate(sample_spools, start=1)
    ]
    _write_tsv(alignment_path, ALIGNMENT_COLUMNS, alignment_rows)
    validate_tsv_table(alignment_path, "sample_alignment.schema.json")
    variant_audit_path = output_dir / "acpa_variant_audit.tsv"
    _write_tsv(variant_audit_path, VARIANT_AUDIT_COLUMNS, variant_audit_rows)
    validate_tsv_table(variant_audit_path, "acpa_variant_audit.schema.json")
    sample_chromosome_qc_path = output_dir / "acpa_sample_chromosome_qc.tsv"
    sample_chromosome_qc_rows = []
    for sample in sample_spools:
        for chromosome in sorted(
            sample.chromosome_marker_counts,
            key=int,
        ):
            marker_count = sample.chromosome_marker_counts[chromosome]
            called_count = sample.chromosome_called_counts.get(chromosome, 0)
            sample_chromosome_qc_rows.append(
                {
                    "SAMPLE_ID": sample.row["SAMPLE_ID"],
                    "CHROMOSOME": chromosome,
                    "MARKER_COUNT": str(marker_count),
                    "CALLED_GENOTYPE_COUNT": str(called_count),
                    "GENOTYPING_RATE": f"{called_count / marker_count:.6f}",
                }
            )
    _write_tsv(
        sample_chromosome_qc_path,
        (
            "SAMPLE_ID",
            "CHROMOSOME",
            "MARKER_COUNT",
            "CALLED_GENOTYPE_COUNT",
            "GENOTYPING_RATE",
        ),
        sample_chromosome_qc_rows,
    )
    validate_tsv_table(
        sample_chromosome_qc_path,
        "acpa_sample_chromosome_qc.schema.json",
    )

    report_path = output_dir / "acpa_conversion_report.json"
    report = {
        "schema_version": "1.0.0",
        "assembly": config["project"]["assembly"],
        "target_chromosome": target_chromosome_value,
        "marker_mode": marker_mode,
        "sample_count": samples_table.row_count,
        "source_read_counts": {
            sample.row["SOURCE_FILE"]: 1 for sample in sample_spools
        },
        "genomewide_variant_count": len(genomewide_markers),
        "target_chromosome_variant_count": len(target_markers),
        "accepted_marker_counts_by_chromosome": {
            chromosome: sum(
                marker.chromosome == chromosome for marker in genomewide_markers
            )
            for chromosome in sorted(AUTOSOMES, key=int)
        },
        "excluded_variant_count": sum(
            row["STATUS"] == "EXCLUDED" for row in variant_audit_rows
        ),
        "invalid_genotype_count": sum(
            sample.invalid_genotype_count for sample in sample_spools
        ),
        "skipped_row_count": sum(sample.skipped_row_count for sample in sample_spools),
    }
    atomic_write_json(report_path, report)

    for prefix in (genomewide_prefix, target_prefix):
        for suffix in ("ped", "map", "log", "nosex"):
            generated_path = prefix.with_suffix(f".{suffix}")
            if generated_path.exists():
                generated_path.unlink()
    shutil.rmtree(spool_dir)

    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    artifacts = _dataset_artifacts(
        genomewide_prefix,
        genomewide_descriptor_path,
        "genomewide_base",
        published_dir,
        producer_stage,
        stage_inputs["signature"],
        config["project"]["assembly"],
        sample_set_id,
        genomewide_variant_set_id,
    )
    artifacts.extend(
        _dataset_artifacts(
            target_prefix,
            target_descriptor_path,
            "target_chromosome_base",
            published_dir,
            producer_stage,
            stage_inputs["signature"],
            config["project"]["assembly"],
            sample_set_id,
            target_variant_set_id,
        )
    )
    for artifact_id, artifact_type, path, schema_name in (
        (
            "sample_alignment",
            "sample_alignment",
            alignment_path,
            "sample_alignment.schema.json",
        ),
        (
            "acpa_variant_audit",
            "variant_audit",
            variant_audit_path,
            "acpa_variant_audit.schema.json",
        ),
        (
            "acpa_sample_chromosome_qc",
            "sample_chromosome_qc",
            sample_chromosome_qc_path,
            "acpa_sample_chromosome_qc.schema.json",
        ),
        ("acpa_conversion_report", "conversion_report", report_path, None),
    ):
        artifacts.append(
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
                assembly=config["project"]["assembly"],
                sample_set_id=sample_set_id,
                sensitivity="sensitive_genetic",
            )
        )

    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "convert_acpa_dual_dataset_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": selected_input_artifacts,
        "outputs": artifacts,
        "parameters": parameters,
        "tools": [
            {
                "tool": "plink",
                "configured": config["tools"]["plink"],
                "version": plink_version,
            }
        ],
        "counts": {
            "samples": samples_table.row_count,
            "genomewide_variants": len(genomewide_markers),
            "target_chromosome_variants": len(target_markers),
            "excluded_variants": report["excluded_variant_count"],
            "invalid_genotypes": report["invalid_genotype_count"],
        },
        "metrics": {},
        "exclusions": [],
        "warnings": [],
        "checks": [
            {"check": "single_source_read", "status": "PASS"},
            {"check": "sample_order_alignment", "status": "PASS"},
            {"check": "plink_binary_contracts", "status": "PASS"},
        ],
        "known_limits": [
            "Les allèles ACPA sont conservés sur le brin forward et ne sont pas "
            "assimilés à REF/ALT.",
            "Les rsID absents ou conflictuels restent des annotations et n'excluent pas une sonde.",
            "Les visualisations QC de conversion seront ajoutées avec l'étape 18.",
        ],
        "expected_visualizations": [],
        "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in artifacts
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
    except ScientificBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        AcpaConversionError,
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
