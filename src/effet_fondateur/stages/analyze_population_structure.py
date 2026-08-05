"""Étape 08 : PCA populationnelle ajustée sur un ensemble indépendant."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

import numpy as np
from scipy.stats import chi2

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


MAX_COMPONENTS = 10
PC_COLUMNS = tuple(f"PC{index}" for index in range(1, MAX_COMPONENTS + 1))
SCORE_COLUMNS = (
    "SAMPLE_ID",
    "SAMPLE_SET_ID",
    "FID",
    "IID",
    "FAMILY_ID",
    "GROUP_LABEL",
    "ARRAY_BATCH",
    "REFERENCE_INCLUDED",
    "PROJECTED",
    *PC_COLUMNS,
    "OUTLIER_STATUS",
    "OUTLIER_REASON",
)
EIGENVALUE_COLUMNS = (
    "COMPONENT",
    "EIGENVALUE",
    "EXPLAINED_VARIANCE_RATIO",
    "CUMULATIVE_EXPLAINED_VARIANCE_RATIO",
    "SINGULAR_VALUE",
    "REFERENCE_SAMPLE_COUNT",
    "INFORMATIVE_VARIANT_COUNT",
)
OUTLIER_COLUMNS = (
    "SAMPLE_ID",
    "SAMPLE_SET_ID",
    "REFERENCE_INCLUDED",
    "PROJECTED",
    "COMPONENTS_USED",
    "DISTANCE_SQUARED",
    "P_VALUE",
    "OUTLIER_ALPHA",
    "OUTLIER_STATUS",
    "OUTLIER_REASON",
    "PROPOSED_EXCLUSION",
    "MANUAL_VALIDATION_REQUIRED",
)
LOADING_COLUMNS = (
    "VARIANT_ID",
    "CHROMOSOME",
    "POSITION_BP",
    "COUNTED_ALLELE",
    "REFERENCE_ALLELE_FREQUENCY",
    "REFERENCE_STANDARD_DEVIATION",
    "REFERENCE_CALL_RATE",
    "MODEL_STATUS",
    "EXCLUSION_CODE",
    *PC_COLUMNS,
)


class PopulationStructureInputError(ValueError):
    """Signale une entrée ou une configuration PCA invalide."""


class ExternalToolError(RuntimeError):
    """Signale un échec de PLINK ou une sortie de dosage invalide."""


class PopulationStructureBlockError(RuntimeError):
    """Signale une base de référence insuffisante pour ajuster une PCA."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="analyze-population-structure")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise PopulationStructureInputError(f"invalid_parameter:{name}")
    return value


def _probability(parameters: dict[str, Any], name: str, default: float) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PopulationStructureInputError(f"invalid_parameter:{name}")
    result = float(value)
    if not math.isfinite(result) or not 0 < result <= 1:
        raise PopulationStructureInputError(f"invalid_parameter:{name}")
    return result


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    requested_components = _positive_integer(parameters, "requested_components", 10)
    outlier_components = _positive_integer(parameters, "outlier_components", 5)
    if requested_components > MAX_COMPONENTS or outlier_components > requested_components:
        raise PopulationStructureInputError("invalid_component_configuration")
    return {
        "requested_components": requested_components,
        "outlier_components": outlier_components,
        "min_reference_samples": _positive_integer(
            parameters, "min_reference_samples", 3
        ),
        "min_informative_variants": _positive_integer(
            parameters, "min_informative_variants", 100
        ),
        "min_reference_call_rate": _probability(
            parameters, "min_reference_call_rate", 0.95
        ),
        "outlier_alpha": _probability(parameters, "outlier_alpha", 0.001),
        "plink_timeout_seconds": _positive_integer(
            parameters, "plink_timeout_seconds", 300
        ),
    }


def _resolve_input_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    artifact_path = Path(artifact["path"])
    if artifact_path.is_absolute():
        return artifact_path
    if PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        return run_dir / artifact_path
    return Path.cwd() / artifact_path


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [
        artifact
        for artifact in stage_inputs["artifacts"]
        if artifact["artifact_id"] == artifact_id
    ]
    if len(matches) != 1:
        raise PopulationStructureInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_artifact_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    physical_path = _resolve_input_path(artifact, run_dir)
    if not physical_path.is_file() or sha256_file(physical_path) != artifact["sha256"]:
        raise PopulationStructureInputError("declared_input_modified")
    return physical_path


def _read_tsv(path: Path, schema_name: str) -> list[dict[str, Any]]:
    return list(validate_tsv_table(path, schema_name).rows)


def _read_fam(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise PopulationStructureInputError("malformed_kinship_panel_fam")
    if len({(row[0], row[1]) for row in rows}) != len(rows):
        raise PopulationStructureInputError("duplicate_kinship_panel_sample")
    return rows


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise PopulationStructureInputError("malformed_kinship_panel_bim")
    if len({row[1] for row in rows}) != len(rows):
        raise PopulationStructureInputError("duplicate_kinship_panel_variant")
    return rows


def _validate_panel(
    descriptor: dict[str, Any],
    paths: dict[str, Path],
    artifacts: dict[str, dict[str, Any]],
) -> tuple[list[list[str]], list[list[str]]]:
    validate_json_document(descriptor, "plink_dataset.schema.json")
    fam_rows = _read_fam(paths["kinship_panel_fam"])
    bim_rows = _read_bim(paths["kinship_panel_bim"])
    if descriptor["source_format"] != "PLINK_KINSHIP_PANEL":
        raise PopulationStructureInputError("unexpected_panel_source_format")
    if descriptor["sample_count"] != len(fam_rows) or descriptor["variant_count"] != len(bim_rows):
        raise PopulationStructureInputError("panel_descriptor_count_mismatch")
    for suffix in ("bed", "bim", "fam"):
        artifact_id = f"kinship_panel_{suffix}"
        if descriptor["files"][suffix]["sha256"] != artifacts[artifact_id]["sha256"]:
            raise PopulationStructureInputError("panel_descriptor_checksum_mismatch")
    stems = {paths[f"kinship_panel_{suffix}"].with_suffix("") for suffix in ("bed", "bim", "fam")}
    if len(stems) != 1:
        raise PopulationStructureInputError("panel_files_do_not_share_prefix")
    return fam_rows, bim_rows


def _sample_context(
    samples_path: Path,
    proposal_path: Path,
    fam_rows: list[list[str]],
) -> tuple[list[str], dict[str, dict[str, Any]], dict[tuple[str, str], str], set[str]]:
    sample_rows = _read_tsv(samples_path, "samples_master.schema.json")
    proposal_rows = _read_tsv(proposal_path, "independent_set_proposal.schema.json")
    samples = {row["SAMPLE_ID"]: row for row in sample_rows}
    proposals = {row["SAMPLE_ID"]: row for row in proposal_rows}
    plink_to_sample = {(row["FID"], row["IID"]): row["SAMPLE_ID"] for row in sample_rows}
    ordered_sample_ids: list[str] = []
    for fam_row in fam_rows:
        sample_id = plink_to_sample.get((fam_row[0], fam_row[1]))
        if sample_id is None or sample_id not in proposals:
            raise PopulationStructureInputError("panel_sample_missing_from_context")
        proposal = proposals[sample_id]
        if proposal["FID"] != fam_row[0] or proposal["IID"] != fam_row[1]:
            raise PopulationStructureInputError("proposal_plink_identifier_mismatch")
        ordered_sample_ids.append(sample_id)
    if set(ordered_sample_ids) != set(proposals):
        raise PopulationStructureInputError("proposal_sample_set_mismatch")
    reference_ids = {
        sample_id for sample_id in ordered_sample_ids if proposals[sample_id]["PROPOSED_INCLUDED"]
    }
    return ordered_sample_ids, samples, plink_to_sample, reference_ids


def _resolve_executable(command: str | None) -> str:
    if command is None:
        raise ExternalToolError("plink_not_configured")
    executable = shutil.which(command)
    if executable is None:
        raise ExternalToolError("plink_not_available")
    return executable


def _plink_version(executable: str) -> str | None:
    try:
        completed = subprocess.run(
            [executable, "--version"], capture_output=True, check=False, text=True, timeout=10
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    output = (completed.stdout or completed.stderr).strip().splitlines()
    return output[0][:200] if output else None


def _export_dosages(
    executable: str,
    panel_prefix: Path,
    output_prefix: Path,
    timeout_seconds: int,
) -> Path:
    command = [
        executable,
        "--bfile",
        str(panel_prefix),
        "--recode",
        "A",
        "--out",
        str(output_prefix),
    ]
    try:
        completed = subprocess.run(
            command, capture_output=True, check=False, text=True, timeout=timeout_seconds
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExternalToolError("plink_dosage_export_failed") from error
    if completed.returncode != 0:
        raise ExternalToolError(f"plink_dosage_export_failed:{completed.returncode}")
    raw_path = output_prefix.with_suffix(".raw")
    if not raw_path.is_file():
        raise ExternalToolError("plink_dosage_output_missing")
    return raw_path


def _read_dosages(
    raw_path: Path,
    bim_rows: list[list[str]],
    fam_rows: list[list[str]],
) -> np.ndarray:
    with raw_path.open(encoding="utf-8", newline="") as input_file:
        reader = csv.reader(input_file, delimiter=" ", skipinitialspace=True)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ExternalToolError("empty_plink_raw") from error
        data_rows = [row for row in reader if row]
    expected_columns = [f"{row[1]}_{row[4]}" for row in bim_rows]
    if header[:6] != ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]:
        raise ExternalToolError("malformed_plink_raw_header")
    if header[6:] != expected_columns:
        raise ExternalToolError("plink_raw_variant_order_or_allele_mismatch")
    if len(data_rows) != len(fam_rows) or any(len(row) != len(header) for row in data_rows):
        raise ExternalToolError("malformed_plink_raw_rows")
    if [(row[0], row[1]) for row in data_rows] != [(row[0], row[1]) for row in fam_rows]:
        raise ExternalToolError("plink_raw_sample_order_mismatch")
    dosages = np.empty((len(data_rows), len(bim_rows)), dtype=np.float64)
    for sample_index, row in enumerate(data_rows):
        for variant_index, value in enumerate(row[6:]):
            if value == "NA":
                dosages[sample_index, variant_index] = np.nan
                continue
            try:
                dosage = float(value)
            except ValueError as error:
                raise ExternalToolError("invalid_plink_dosage") from error
            if not math.isfinite(dosage) or not 0 <= dosage <= 2:
                raise ExternalToolError("invalid_plink_dosage")
            dosages[sample_index, variant_index] = dosage
    return dosages


def _fit_and_project(
    dosages: np.ndarray,
    reference_indices: np.ndarray,
    parameters: dict[str, Any],
) -> dict[str, Any]:
    reference_dosages = dosages[reference_indices]
    nonmissing_counts = np.sum(~np.isnan(reference_dosages), axis=0)
    call_rates = nonmissing_counts / reference_dosages.shape[0]
    dosage_sums = np.nansum(reference_dosages, axis=0)
    means = np.divide(
        dosage_sums,
        nonmissing_counts,
        out=np.full(reference_dosages.shape[1], np.nan, dtype=np.float64),
        where=nonmissing_counts > 0,
    )
    frequencies = means / 2.0
    standard_deviations = np.sqrt(2.0 * frequencies * (1.0 - frequencies))
    call_rate_eligible = call_rates >= parameters["min_reference_call_rate"]
    polymorphic = np.isfinite(standard_deviations) & (standard_deviations > 0)
    informative = call_rate_eligible & polymorphic
    informative_count = int(np.sum(informative))
    if informative_count < parameters["min_informative_variants"]:
        raise PopulationStructureBlockError("insufficient_informative_variants")
    imputed = np.where(np.isnan(dosages[:, informative]), means[informative], dosages[:, informative])
    standardized = (imputed - means[informative]) / standard_deviations[informative]
    normalized = standardized / math.sqrt(informative_count)
    reference_matrix = normalized[reference_indices]
    left_vectors, singular_values_all, right_vectors_all = np.linalg.svd(
        reference_matrix, full_matrices=False
    )
    del left_vectors
    component_count = min(
        parameters["requested_components"],
        reference_matrix.shape[0] - 1,
        informative_count,
    )
    if component_count < 1:
        raise PopulationStructureBlockError("insufficient_pca_rank")
    singular_values = singular_values_all[:component_count]
    right_vectors = right_vectors_all[:component_count].copy()
    for component_index in range(component_count):
        anchor_index = int(np.argmax(np.abs(right_vectors[component_index])))
        if right_vectors[component_index, anchor_index] < 0:
            right_vectors[component_index] *= -1
    scores = normalized @ right_vectors.T
    eigenvalues = singular_values**2 / (reference_matrix.shape[0] - 1)
    total_variance = float(np.sum(reference_matrix**2) / (reference_matrix.shape[0] - 1))
    explained_ratios = eigenvalues / total_variance
    outlier_component_count = min(parameters["outlier_components"], component_count)
    positive_eigenvalues = eigenvalues[:outlier_component_count] > np.finfo(float).eps
    outlier_component_count = int(np.sum(positive_eigenvalues))
    if outlier_component_count < 1:
        raise PopulationStructureBlockError("zero_variance_pca")
    standardized_scores = scores[:, :outlier_component_count] / np.sqrt(
        eigenvalues[:outlier_component_count]
    )
    distances_squared = np.sum(standardized_scores**2, axis=1)
    p_values = chi2.sf(distances_squared, df=outlier_component_count)
    return {
        "means": means,
        "frequencies": frequencies,
        "standard_deviations": standard_deviations,
        "call_rates": call_rates,
        "call_rate_eligible": call_rate_eligible,
        "informative": informative,
        "scores": scores,
        "singular_values": singular_values,
        "eigenvalues": eigenvalues,
        "explained_ratios": explained_ratios,
        "right_vectors": right_vectors,
        "distances_squared": distances_squared,
        "p_values": p_values,
        "component_count": component_count,
        "outlier_component_count": outlier_component_count,
        "informative_count": informative_count,
    }


def _decimal(value: float) -> str:
    if abs(value) < 5e-16:
        value = 0.0
    return format(float(value), ".12g")


def _reference_set_id(reference_ids: set[str]) -> str:
    digest = hashlib.sha256("\n".join(sorted(reference_ids)).encode("utf-8")).hexdigest()
    return f"population_reference_{digest[:16]}"


def _result_rows(
    ordered_sample_ids: list[str],
    samples: dict[str, dict[str, Any]],
    reference_ids: set[str],
    sample_set_id: str,
    model: dict[str, Any],
    parameters: dict[str, Any],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    score_rows: list[dict[str, str]] = []
    outlier_rows: list[dict[str, str]] = []
    for sample_index, sample_id in enumerate(ordered_sample_ids):
        reference_included = sample_id in reference_ids
        outlier = bool(model["p_values"][sample_index] < parameters["outlier_alpha"])
        status = "OUTLIER" if outlier else "INLIER"
        reason = "MULTIVARIATE_PC_DISTANCE" if outlier else ""
        sample = samples[sample_id]
        score_row = {
            "SAMPLE_ID": sample_id,
            "SAMPLE_SET_ID": sample_set_id,
            "FID": sample["FID"],
            "IID": sample["IID"],
            "FAMILY_ID": sample["FID"],
            "GROUP_LABEL": sample["GROUP_LABEL"],
            "ARRAY_BATCH": sample["ARRAY_BATCH"] or "",
            "REFERENCE_INCLUDED": str(reference_included).lower(),
            "PROJECTED": str(not reference_included).lower(),
            "OUTLIER_STATUS": status,
            "OUTLIER_REASON": reason,
        }
        for component_index, column in enumerate(PC_COLUMNS):
            score_row[column] = (
                _decimal(model["scores"][sample_index, component_index])
                if component_index < model["component_count"]
                else ""
            )
        score_rows.append(score_row)
        outlier_rows.append(
            {
                "SAMPLE_ID": sample_id,
                "SAMPLE_SET_ID": sample_set_id,
                "REFERENCE_INCLUDED": str(reference_included).lower(),
                "PROJECTED": str(not reference_included).lower(),
                "COMPONENTS_USED": str(model["outlier_component_count"]),
                "DISTANCE_SQUARED": _decimal(model["distances_squared"][sample_index]),
                "P_VALUE": _decimal(model["p_values"][sample_index]),
                "OUTLIER_ALPHA": _decimal(parameters["outlier_alpha"]),
                "OUTLIER_STATUS": status,
                "OUTLIER_REASON": reason,
                "PROPOSED_EXCLUSION": str(outlier).lower(),
                "MANUAL_VALIDATION_REQUIRED": str(outlier).lower(),
            }
        )
    return score_rows, outlier_rows


def _eigenvalue_rows(model: dict[str, Any], reference_count: int) -> list[dict[str, str]]:
    rows = []
    cumulative = 0.0
    for component_index in range(model["component_count"]):
        cumulative += float(model["explained_ratios"][component_index])
        rows.append(
            {
                "COMPONENT": f"PC{component_index + 1}",
                "EIGENVALUE": _decimal(model["eigenvalues"][component_index]),
                "EXPLAINED_VARIANCE_RATIO": _decimal(model["explained_ratios"][component_index]),
                "CUMULATIVE_EXPLAINED_VARIANCE_RATIO": _decimal(cumulative),
                "SINGULAR_VALUE": _decimal(model["singular_values"][component_index]),
                "REFERENCE_SAMPLE_COUNT": str(reference_count),
                "INFORMATIVE_VARIANT_COUNT": str(model["informative_count"]),
            }
        )
    return rows


def _loading_rows(bim_rows: list[list[str]], model: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    informative_indices = np.flatnonzero(model["informative"])
    loading_index = {int(variant_index): position for position, variant_index in enumerate(informative_indices)}
    for variant_index, bim_row in enumerate(bim_rows):
        included = bool(model["informative"][variant_index])
        if not model["call_rate_eligible"][variant_index]:
            exclusion_code = "LOW_REFERENCE_CALL_RATE"
        elif not included:
            exclusion_code = "MONOMORPHIC_IN_REFERENCE"
        else:
            exclusion_code = ""
        row = {
            "VARIANT_ID": bim_row[1],
            "CHROMOSOME": bim_row[0],
            "POSITION_BP": bim_row[3],
            "COUNTED_ALLELE": bim_row[4],
            "REFERENCE_ALLELE_FREQUENCY": (
                _decimal(model["frequencies"][variant_index])
                if math.isfinite(float(model["frequencies"][variant_index]))
                else ""
            ),
            "REFERENCE_STANDARD_DEVIATION": (
                _decimal(model["standard_deviations"][variant_index])
                if math.isfinite(float(model["standard_deviations"][variant_index]))
                else ""
            ),
            "REFERENCE_CALL_RATE": _decimal(model["call_rates"][variant_index]),
            "MODEL_STATUS": "INCLUDED" if included else "EXCLUDED",
            "EXCLUSION_CODE": exclusion_code,
        }
        for component_index, column in enumerate(PC_COLUMNS):
            row[column] = (
                _decimal(model["right_vectors"][component_index, loading_index[variant_index]])
                if included and component_index < model["component_count"]
                else ""
            )
        rows.append(row)
    return rows


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row[column] for column in columns} for row in rows)


def _output_artifacts(
    stage_inputs: dict[str, Any],
    descriptor: dict[str, Any],
    sample_set_id: str,
    contracts: list[tuple[str, str, Path, str | None, str]],
) -> list[dict[str, Any]]:
    producer_stage = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    published_dir = PurePosixPath(stage_inputs["published_output_dir"])
    artifacts = []
    for artifact_id, artifact_type, path, schema_name, sensitivity in contracts:
        artifacts.append(
            build_file_artifact(
                physical_path=path,
                published_path=str(published_dir / path.name),
                artifact_id=artifact_id,
                artifact_type=artifact_type,
                media_type="application/json" if path.suffix == ".json" else "text/tab-separated-values",
                producer_stage=producer_stage,
                producer_signature=stage_inputs["signature"],
                schema_name=schema_name,
                schema_version="1.0.0" if schema_name else None,
                assembly=descriptor["assembly"],
                sample_set_id=sample_set_id,
                variant_set_id=descriptor["variant_set_id"],
                sensitivity=sensitivity,
            )
        )
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Ajuste la PCA indépendante et projette tous les individus du panel."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "samples_master",
        "kinship_panel_bed",
        "kinship_panel_bim",
        "kinship_panel_fam",
        "kinship_panel_dataset",
        "independent_set_proposal",
    )
    input_artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_artifact_path(artifact, run_dir) for artifact_id, artifact in input_artifacts.items()}
    descriptor = read_json(paths["kinship_panel_dataset"])
    fam_rows, bim_rows = _validate_panel(descriptor, paths, input_artifacts)
    ordered_sample_ids, samples, _, reference_ids = _sample_context(
        paths["samples_master"], paths["independent_set_proposal"], fam_rows
    )
    if len(reference_ids) < parameters["min_reference_samples"]:
        raise PopulationStructureBlockError("insufficient_independent_reference_samples")
    reference_indices = np.array(
        [index for index, sample_id in enumerate(ordered_sample_ids) if sample_id in reference_ids],
        dtype=np.int64,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    plink_executable = _resolve_executable(config["tools"]["plink"])
    plink_version = _plink_version(plink_executable)
    with tempfile.TemporaryDirectory(prefix="population_structure_", dir=output_dir) as work_dir:
        raw_path = _export_dosages(
            plink_executable,
            paths["kinship_panel_bed"].with_suffix(""),
            Path(work_dir) / "dosages",
            parameters["plink_timeout_seconds"],
        )
        dosages = _read_dosages(raw_path, bim_rows, fam_rows)
    model = _fit_and_project(dosages, reference_indices, parameters)
    reference_set_id = _reference_set_id(reference_ids)
    population_sample_set_id = descriptor["sample_set_id"]
    score_rows, outlier_rows = _result_rows(
        ordered_sample_ids,
        samples,
        reference_ids,
        population_sample_set_id,
        model,
        parameters,
    )
    eigenvalue_rows = _eigenvalue_rows(model, len(reference_ids))
    loading_rows = _loading_rows(bim_rows, model)

    score_path = output_dir / "population_scores.tsv"
    eigenvalue_path = output_dir / "population_eigenvalues.tsv"
    outlier_path = output_dir / "population_outliers.tsv"
    loading_path = output_dir / "population_variant_loadings.tsv"
    _write_tsv(score_path, SCORE_COLUMNS, score_rows)
    _write_tsv(eigenvalue_path, EIGENVALUE_COLUMNS, eigenvalue_rows)
    _write_tsv(outlier_path, OUTLIER_COLUMNS, outlier_rows)
    _write_tsv(loading_path, LOADING_COLUMNS, loading_rows)
    validate_tsv_table(score_path, "population_scores.schema.json")
    validate_tsv_table(eigenvalue_path, "population_eigenvalues.schema.json")
    validate_tsv_table(outlier_path, "population_outliers.schema.json")
    validate_tsv_table(loading_path, "population_variant_loadings.schema.json")

    outlier_count = sum(row["OUTLIER_STATUS"] == "OUTLIER" for row in outlier_rows)
    report_path = output_dir / "population_structure_report.json"
    report = {
        "schema_version": "1.0.0",
        "method_id": "independent_reference_pca_projection_v1",
        "reference_sample_set_id": reference_set_id,
        "population_sample_set_id": population_sample_set_id,
        "panel_sample_count": len(ordered_sample_ids),
        "reference_sample_count": len(reference_ids),
        "projected_sample_count": len(ordered_sample_ids) - len(reference_ids),
        "panel_variant_count": len(bim_rows),
        "informative_variant_count": model["informative_count"],
        "component_count": model["component_count"],
        "outlier_component_count": model["outlier_component_count"],
        "outlier_count": outlier_count,
        "outlier_method": "chi_square_distance_on_reference_scaled_pc_scores",
        "manual_validation_required": outlier_count > 0,
    }
    atomic_write_json(report_path, report)
    contracts = [
        ("population_scores", "population_scores", score_path, "population_scores.schema.json", "sensitive_genetic"),
        ("population_eigenvalues", "population_eigenvalues", eigenvalue_path, "population_eigenvalues.schema.json", "internal"),
        ("population_outliers", "population_outliers", outlier_path, "population_outliers.schema.json", "sensitive_genetic"),
        ("population_variant_loadings", "population_variant_loadings", loading_path, "population_variant_loadings.schema.json", "sensitive_genetic"),
        ("population_structure_report", "population_structure_report", report_path, None, "internal"),
    ]
    output_artifacts = _output_artifacts(
        stage_inputs, descriptor, population_sample_set_id, contracts
    )
    stage_outputs = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"],
        "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0",
        "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"],
        "method_id": "independent_reference_pca_projection_v1",
        "signature": stage_inputs["signature"],
        "started_at": started_at,
        "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(input_artifacts.values()),
        "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [
            {"tool": "plink", "configured": config["tools"]["plink"], "version": plink_version},
            {"tool": "numpy", "version": np.__version__},
        ],
        "counts": {
            "samples": len(ordered_sample_ids),
            "reference_samples": len(reference_ids),
            "projected_samples": len(ordered_sample_ids) - len(reference_ids),
            "panel_variants": len(bim_rows),
            "informative_variants": model["informative_count"],
            "outliers": outlier_count,
        },
        "metrics": {
            "component_count": model["component_count"],
            "outlier_component_count": model["outlier_component_count"],
            "first_component_explained_variance_ratio": float(model["explained_ratios"][0]),
        },
        "exclusions": [],
        "warnings": (
            [{"code": "manual_population_outlier_review_required", "outlier_count": outlier_count}]
            if outlier_count
            else []
        ),
        "checks": [
            {"check": "kinship_panel_integrity", "status": "PASS"},
            {"check": "independent_reference", "status": "PASS"},
            {"check": "pca_model", "status": "PASS"},
            {"check": "population_outliers", "status": "REVIEW_REQUIRED" if outlier_count else "PASS"},
        ],
        "known_limits": [
            "La référence indépendante provient de la proposition de l'étape 07 et peut rester en attente d'approbation humaine.",
            "Les p-values d'outlier sont exploratoires, reposent sur une approximation du chi-deux et ne constituent pas un critère clinique.",
            "Aucune référence populationnelle externe n'est utilisée ; les axes décrivent uniquement les individus de l'étude.",
            "Aucune exclusion proposée n'est appliquée automatiquement.",
        ],
        "expected_visualizations": [
            "population_scree_plot",
            "pseudonymized_population_pca",
            "population_pca_by_batch",
            "projected_related_samples",
        ],
        "manual_validation_required": outlier_count > 0,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {Path(str(artifact['path'])).name}\n"
            for artifact in output_artifacts
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
    except PopulationStructureBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except (
        PopulationStructureInputError,
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
