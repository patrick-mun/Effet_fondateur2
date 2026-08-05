"""Étape 01 : inventaire et validation non destructive des exports ACPA."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import dataclass
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
DEFAULT_AUTOSOMES = tuple(str(chromosome) for chromosome in range(1, 23))
INVENTORY_COLUMNS = (
    "SOURCE_FILE",
    "SIZE_BYTES",
    "SHA256",
    "ENCODING",
    "ARRAY_TYPE",
    "GENOMIC_BUILD",
    "DATA_ROW_COUNT",
    "AUTOSOMES_PRESENT",
    "TARGET_CHROMOSOME_PRESENT",
    "STATUS",
)
QC_COLUMNS = ("SOURCE_FILE", "CHECK_ID", "STATUS", "DETAIL_CODE")


class SourceValidationError(ValueError):
    """Signale un contrat d'entrée invalide avant l'inspection des exports."""


@dataclass(frozen=True)
class SourceInspection:
    """Résumé non génotypique d'un export ACPA et de ses contrôles."""

    source_file: str
    size_bytes: int
    sha256: str
    array_type: str | None
    genomic_build: str | None
    data_row_count: int
    autosomes_present: tuple[str, ...]
    target_chromosome_present: bool
    checks: tuple[dict[str, str], ...]

    @property
    def passed(self) -> bool:
        return all(check["STATUS"] != "FAIL" for check in self.checks)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="validate-sources")
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


def _check(source_file: str, check_id: str, passed: bool, detail: str) -> dict[str, str]:
    return {
        "SOURCE_FILE": source_file,
        "CHECK_ID": check_id,
        "STATUS": "PASS" if passed else "FAIL",
        "DETAIL_CODE": detail if not passed else "valid",
    }


def _read_metadata_and_table(
    source_path: Path,
) -> tuple[dict[str, str], set[str], int, set[str]]:
    metadata: dict[str, str] = {}
    chromosomes: set[str] = set()
    data_row_count = 0
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
            raise SourceValidationError("missing_header")

        reader = csv.DictReader(chain([header_line], source_file), delimiter="\t")
        raw_fieldnames = reader.fieldnames or []
        fieldnames = [field.strip() for field in raw_fieldnames]
        if len(fieldnames) != len(set(fieldnames)):
            raise SourceValidationError("duplicate_columns")
        reader.fieldnames = fieldnames
        missing_columns = REQUIRED_ACPA_COLUMNS - set(fieldnames)
        if missing_columns:
            raise SourceValidationError("missing_required_columns")

        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise SourceValidationError("malformed_row")
            if not any(str(value).strip() for value in row.values()):
                continue
            data_row_count += 1
            chromosome = _normalize_chromosome(str(row["Chromosome"]))
            if chromosome:
                chromosomes.add(chromosome)
    return metadata, set(fieldnames), data_row_count, chromosomes


def inspect_source(
    source_path: Path,
    source_file: str,
    expected_sha256: str,
    assembly: str,
    target_chromosome: str,
    required_autosomes: tuple[str, ...],
    allowed_extensions: tuple[str, ...],
) -> SourceInspection:
    """Contrôle un export sans conserver ni interpréter ses génotypes."""
    checks = [
        _check(
            source_file,
            "declared_sha256",
            sha256_file(source_path) == expected_sha256,
            "source_modified_after_declaration",
        ),
        _check(
            source_file,
            "file_extension",
            source_path.suffix.lower() in allowed_extensions,
            "unsupported_extension",
        ),
    ]
    metadata: dict[str, str] = {}
    data_row_count = 0
    chromosomes: set[str] = set()
    parse_error: str | None = None
    try:
        metadata, _, data_row_count, chromosomes = _read_metadata_and_table(source_path)
    except UnicodeDecodeError:
        parse_error = "invalid_utf8"
    except SourceValidationError as error:
        parse_error = str(error)

    checks.append(
        _check(
            source_file,
            "acpa_structure",
            parse_error is None,
            parse_error or "valid",
        )
    )
    array_type = metadata.get("Array Type Name") or None
    genomic_build = metadata.get("UCSC Genomic Version") or None
    checks.extend(
        [
            _check(
                source_file,
                "array_type_metadata",
                array_type is not None,
                "missing_array_type",
            ),
            _check(
                source_file,
                "genomic_build",
                _normalized_build(genomic_build) == _normalized_build(assembly),
                "incompatible_or_missing_build",
            ),
            _check(
                source_file,
                "data_rows",
                data_row_count > 0,
                "empty_or_truncated_source",
            ),
            _check(
                source_file,
                "autosome_coverage",
                set(required_autosomes).issubset(chromosomes),
                "missing_required_autosomes",
            ),
            _check(
                source_file,
                "target_chromosome",
                target_chromosome in chromosomes,
                "missing_target_chromosome",
            ),
        ]
    )
    autosomes_present = tuple(
        str(chromosome)
        for chromosome in range(1, 23)
        if str(chromosome) in chromosomes
    )
    return SourceInspection(
        source_file=source_file,
        size_bytes=source_path.stat().st_size,
        sha256=sha256_file(source_path),
        array_type=array_type,
        genomic_build=genomic_build,
        data_row_count=data_row_count,
        autosomes_present=autosomes_present,
        target_chromosome_present=target_chromosome in chromosomes,
        checks=tuple(checks),
    )


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _source_root(configured_root: str) -> Path:
    path = Path(configured_root)
    return path if path.is_absolute() else Path.cwd() / path


def _declared_sources(
    stage_inputs: dict[str, Any],
    root_path: Path,
) -> list[tuple[dict[str, Any], Path, str]]:
    declared = []
    declared_paths: set[Path] = set()
    for artifact in stage_inputs["artifacts"]:
        if artifact["artifact_type"] != "acpa_source_file":
            raise SourceValidationError("unexpected_input_artifact")
        physical_path = Path(artifact["path"])
        if not physical_path.is_absolute():
            physical_path = Path.cwd() / physical_path
        try:
            relative_path = physical_path.resolve().relative_to(root_path.resolve())
        except ValueError as error:
            raise SourceValidationError("source_outside_configured_root") from error
        if physical_path.is_symlink() or not physical_path.is_file():
            raise SourceValidationError("unsafe_or_missing_source")
        declared_paths.add(physical_path.resolve())
        declared.append((artifact, physical_path, relative_path.as_posix()))
    current_entries = list(root_path.rglob("*"))
    if any(path.is_symlink() for path in current_entries):
        raise SourceValidationError("symlink_in_source_root")
    current_paths = {path.resolve() for path in current_entries if path.is_file()}
    if current_paths != declared_paths:
        raise SourceValidationError("declared_source_set_changed")
    return declared


def _global_checks(
    inspections: list[SourceInspection],
    expected_file_count: int | None,
) -> list[dict[str, str]]:
    source_file = "__SOURCE_SET__"
    basename_counts = Counter(Path(item.source_file).name for item in inspections)
    builds = {
        _normalized_build(item.genomic_build)
        for item in inspections
        if item.genomic_build is not None
    }
    return [
        _check(
            source_file,
            "expected_file_count",
            expected_file_count is None or len(inspections) == expected_file_count,
            "unexpected_file_count",
        ),
        _check(
            source_file,
            "unique_filenames",
            all(count == 1 for count in basename_counts.values()),
            "duplicate_basename",
        ),
        _check(
            source_file,
            "annotation_build_consistency",
            len(builds) == 1,
            "inconsistent_annotation_builds",
        ),
    ]


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Valide les sources déclarées et publie uniquement leur inventaire technique."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config_path = run_dir / "config.resolved.yaml"
    config = load_pipeline_config(config_path)
    configured_root = config["inputs"]["acpa_samples_dir"]
    if configured_root is None:
        raise SourceValidationError("missing_acpa_samples_dir")
    target_chromosome_value = config["target"]["chromosome"]
    if target_chromosome_value is None:
        raise SourceValidationError("missing_target_chromosome")

    parameters = stage_inputs["parameters"]
    allowed_extensions = tuple(
        str(extension).lower()
        for extension in parameters.get("allowed_extensions", [".txt"])
    )
    required_autosomes = tuple(
        str(chromosome)
        for chromosome in parameters.get("required_autosomes", DEFAULT_AUTOSOMES)
    )
    expected_file_count_value = parameters.get("expected_file_count")
    expected_file_count = (
        None if expected_file_count_value is None else int(expected_file_count_value)
    )
    root_path = _source_root(configured_root)
    declared_sources = _declared_sources(stage_inputs, root_path)
    if not declared_sources:
        raise SourceValidationError("no_declared_sources")

    inspections = [
        inspect_source(
            source_path=physical_path,
            source_file=relative_path,
            expected_sha256=artifact["sha256"],
            assembly=config["project"]["assembly"],
            target_chromosome=str(target_chromosome_value),
            required_autosomes=required_autosomes,
            allowed_extensions=allowed_extensions,
        )
        for artifact, physical_path, relative_path in declared_sources
    ]
    global_checks = _global_checks(inspections, expected_file_count)
    duplicate_basenames = {
        basename
        for basename, count in Counter(
            Path(item.source_file).name for item in inspections
        ).items()
        if count > 1
    }
    if duplicate_basenames:
        inspections = [
            SourceInspection(
                **{
                    **inspection.__dict__,
                    "checks": inspection.checks
                    + (
                        _check(
                            inspection.source_file,
                            "unique_filename",
                            Path(inspection.source_file).name not in duplicate_basenames,
                            "duplicate_basename",
                        ),
                    ),
                }
            )
            for inspection in inspections
        ]

    inventory_rows = [
        {
            "SOURCE_FILE": item.source_file,
            "SIZE_BYTES": str(item.size_bytes),
            "SHA256": item.sha256,
            "ENCODING": "utf-8-sig",
            "ARRAY_TYPE": item.array_type or "",
            "GENOMIC_BUILD": item.genomic_build or "",
            "DATA_ROW_COUNT": str(item.data_row_count),
            "AUTOSOMES_PRESENT": ",".join(item.autosomes_present),
            "TARGET_CHROMOSOME_PRESENT": str(item.target_chromosome_present).lower(),
            "STATUS": "PASS" if item.passed else "FAIL",
        }
        for item in inspections
    ]
    qc_rows = [check for item in inspections for check in item.checks] + global_checks
    output_dir.mkdir(parents=True, exist_ok=True)
    inventory_path = output_dir / "source_inventory.tsv"
    qc_path = output_dir / "source_qc.tsv"
    _write_tsv(inventory_path, INVENTORY_COLUMNS, inventory_rows)
    _write_tsv(qc_path, QC_COLUMNS, qc_rows)
    validate_tsv_table(inventory_path, "source_inventory.schema.json")
    validate_tsv_table(qc_path, "source_qc.schema.json")

    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    inventory_artifact = build_file_artifact(
        physical_path=inventory_path,
        published_path=str(published_dir / inventory_path.name),
        artifact_id="source_inventory",
        artifact_type="source_inventory",
        media_type="text/tab-separated-values",
        producer_stage=producer_stage,
        producer_signature=stage_inputs["signature"],
        schema_name="source_inventory.schema.json",
        schema_version="1.0.0",
        assembly=config["project"]["assembly"],
        sensitivity="sensitive_genetic",
    )
    qc_artifact = build_file_artifact(
        physical_path=qc_path,
        published_path=str(published_dir / qc_path.name),
        artifact_id="source_qc",
        artifact_type="source_qc",
        media_type="text/tab-separated-values",
        producer_stage=producer_stage,
        producer_signature=stage_inputs["signature"],
        schema_name="source_qc.schema.json",
        schema_version="1.0.0",
        assembly=config["project"]["assembly"],
        sensitivity="sensitive_genetic",
    )
    artifacts = [inventory_artifact, qc_artifact]
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

    all_checks = qc_rows
    failed_check_count = sum(check["STATUS"] == "FAIL" for check in all_checks)
    array_types = sorted(
        {item.array_type for item in inspections if item.array_type is not None}
    )
    warnings = (
        [{"code": "multiple_array_types"}] if len(array_types) > 1 else []
    )
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "validate_acpa_sources_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": stage_inputs["artifacts"],
        "outputs": artifacts,
        "parameters": parameters,
        "tools": [],
        "counts": {
            "source_files": len(inspections),
            "passed_files": sum(item.passed for item in inspections),
            "failed_checks": failed_check_count,
        },
        "metrics": {},
        "exclusions": [],
        "warnings": warnings,
        "checks": [
            {
                "check": "source_contracts",
                "status": "PASS" if failed_check_count == 0 else "FAIL",
            }
        ],
        "known_limits": [
            "Cette étape valide la structure des exports sans interpréter les génotypes."
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
    return 0 if failed_check_count == 0 else 4


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)
    try:
        return execute(parsed_arguments.stage_inputs, parsed_arguments.output_dir)
    except (
        DocumentValidationError,
        TableValidationError,
        SourceValidationError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
