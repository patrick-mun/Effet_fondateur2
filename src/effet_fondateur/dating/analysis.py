"""Publication des entrées, estimations et sensibilités de datation V2."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json
from effet_fondateur.contracts import validate_json_document, validate_tsv_table
from effet_fondateur.dating.gamma import (
    GammaAgeError,
    GammaEstimate,
    chance_sharing_correction_cm,
    estimate_correlated_gamma,
    estimate_independent_gamma,
)


@dataclass(frozen=True)
class VariantAgePublication:
    """Sorties et décision synthétique de l'étape 14."""

    input_path: Path
    estimates_path: Path
    scenarios_path: Path
    summary_path: Path
    status: str
    included_unit_count: int
    scenario_count: int


INPUT_COLUMNS = (
    "INDEPENDENT_UNIT_ID", "FAMILY_ID", "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM",
    "INCLUDED", "EXCLUSION_CODE", "MAP_VERSION", "SEGMENT_METHOD",
)
ESTIMATE_COLUMNS = (
    "METHOD_ID", "MODEL", "ANALYSIS_STATUS", "N_UNITS", "EFFECTIVE_N", "RHO",
    "ESTIMATE_GENERATIONS", "CI_LOWER_GENERATIONS", "CI_UPPER_GENERATIONS",
    "CONFIDENCE_LEVEL", "PRIMARY", "EXCLUSION_CODE",
)
SCENARIO_COLUMNS = (
    "SCENARIO_ID", "SCENARIO_TYPE", "MODEL", "OMITTED_INDEPENDENT_UNIT_ID",
    "N_UNITS", "ESTIMATE_GENERATIONS", "CI_LOWER_GENERATIONS",
    "CI_UPPER_GENERATIONS", "GENERATION_YEARS", "ESTIMATE_YEARS",
    "CI_LOWER_YEARS", "CI_UPPER_YEARS", "STATUS", "EXCLUSION_CODE",
)


def _decimal(value: float) -> str:
    return f"{value:.12g}"


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _estimate_row(
    estimate: GammaEstimate | None, model: str, analysis_status: str,
    confidence_level: float, unit_count: int, exclusion_code: str,
) -> dict[str, str]:
    return {
        "METHOD_ID": "gamma_gandolfo_2014_v1", "MODEL": model,
        "ANALYSIS_STATUS": analysis_status, "N_UNITS": str(unit_count),
        "EFFECTIVE_N": "" if estimate is None else _decimal(estimate.effective_sample_size),
        "RHO": "" if estimate is None or estimate.correlation is None else _decimal(estimate.correlation),
        "ESTIMATE_GENERATIONS": "" if estimate is None else _decimal(estimate.generations),
        "CI_LOWER_GENERATIONS": "" if estimate is None else _decimal(estimate.confidence_lower),
        "CI_UPPER_GENERATIONS": "" if estimate is None else _decimal(estimate.confidence_upper),
        "CONFIDENCE_LEVEL": _decimal(confidence_level),
        "PRIMARY": "true" if analysis_status == "PRIMARY" and model == "CORRELATED" else "false",
        "EXCLUSION_CODE": exclusion_code,
    }


def _scenario_row(
    scenario_id: str, scenario_type: str, model: str,
    estimate: GammaEstimate | None, *, omitted_unit_id: str = "",
    generation_years: float | None = None, exclusion_code: str = "",
    unit_count: int = 0,
) -> dict[str, str]:
    multiplier = generation_years or 0.0
    return {
        "SCENARIO_ID": scenario_id, "SCENARIO_TYPE": scenario_type, "MODEL": model,
        "OMITTED_INDEPENDENT_UNIT_ID": omitted_unit_id,
        "N_UNITS": str(unit_count if estimate is None else estimate.sample_size),
        "ESTIMATE_GENERATIONS": "" if estimate is None else _decimal(estimate.generations),
        "CI_LOWER_GENERATIONS": "" if estimate is None else _decimal(estimate.confidence_lower),
        "CI_UPPER_GENERATIONS": "" if estimate is None else _decimal(estimate.confidence_upper),
        "GENERATION_YEARS": "" if generation_years is None else _decimal(generation_years),
        "ESTIMATE_YEARS": "" if estimate is None or generation_years is None else _decimal(estimate.generations * multiplier),
        "CI_LOWER_YEARS": "" if estimate is None or generation_years is None else _decimal(estimate.confidence_lower * multiplier),
        "CI_UPPER_YEARS": "" if estimate is None or generation_years is None else _decimal(estimate.confidence_upper * multiplier),
        "STATUS": "NOT_ESTIMATED" if estimate is None else "ESTIMATED",
        "EXCLUSION_CODE": exclusion_code,
    }


def publish_variant_age(
    *, founder_segments_path: Path, founder_consensus_path: Path,
    founder_summary_path: Path, output_dir: Path, confidence_level: float,
    minimum_units_for_estimate: int, minimum_units_for_primary: int,
    generation_years: Sequence[float], leave_one_family_out: bool,
    chance_sharing_parameters: dict[str, float | int] | None,
) -> VariantAgePublication:
    """Dérive les entrées Gamma et publie estimations et non-conclusions."""
    if minimum_units_for_estimate < 3:
        raise GammaAgeError("minimum_units_for_estimate_below_three")
    if minimum_units_for_primary < 5 or minimum_units_for_primary < minimum_units_for_estimate:
        raise GammaAgeError("minimum_units_for_primary_invalid")
    if (
        not generation_years or any(value <= 0 for value in generation_years)
        or len(set(generation_years)) != len(generation_years)
    ):
        raise GammaAgeError("generation_years_invalid")
    segments = validate_tsv_table(
        founder_segments_path, "founder_segments.schema.json"
    ).rows
    consensus_rows = validate_tsv_table(
        founder_consensus_path, "founder_consensus.schema.json"
    ).rows
    if len(consensus_rows) != 1:
        raise GammaAgeError("founder_consensus_missing_or_ambiguous")
    founder_summary = read_json(founder_summary_path)
    validate_json_document(founder_summary, "founder_analysis_summary.schema.json")
    map_version = consensus_rows[0]["MAP_VERSION"]
    founder_supported = founder_summary["status"] == "SUPPORTED_IBS_CANDIDATE"
    input_rows: list[dict[str, str]] = []
    included_source_rows: list[dict[str, Any]] = []
    for segment in segments:
        left_length = float(segment["LEFT_LENGTH_CM"])
        right_length = float(segment["RIGHT_LENGTH_CM"])
        exclusion_code = segment["EXCLUSION_CODE"] or ""
        included = segment["SEGMENT_STATUS"] == "INCLUDED" and founder_supported
        if included and segment["SEGMENT_METHOD"] != "target_centered_pairwise_max_ibs_v1":
            included = False
            exclusion_code = "NON_INDIVIDUAL_SEGMENT_METHOD"
        elif included and left_length <= 0:
            included = False
            exclusion_code = "ZERO_LEFT_LENGTH"
        elif included and right_length <= 0:
            included = False
            exclusion_code = "ZERO_RIGHT_LENGTH"
        elif not included and not exclusion_code:
            exclusion_code = founder_summary["status"]
        input_rows.append({
            "INDEPENDENT_UNIT_ID": segment["INDEPENDENT_UNIT_ID"],
            "FAMILY_ID": segment["FAMILY_ID"],
            "LEFT_LENGTH_CM": segment["LEFT_LENGTH_CM"],
            "RIGHT_LENGTH_CM": segment["RIGHT_LENGTH_CM"],
            "INCLUDED": "true" if included else "false",
            "EXCLUSION_CODE": "" if included else exclusion_code,
            "MAP_VERSION": map_version,
            "SEGMENT_METHOD": segment["SEGMENT_METHOD"],
        })
        if included:
            included_source_rows.append(segment)
    included_families = [row["FAMILY_ID"] for row in included_source_rows]
    if len(included_families) != len(set(included_families)):
        raise GammaAgeError("gamma_included_families_not_independent")
    left_lengths = [float(row["LEFT_LENGTH_CM"]) for row in included_source_rows]
    right_lengths = [float(row["RIGHT_LENGTH_CM"]) for row in included_source_rows]
    unit_count = len(included_source_rows)
    correction_cm = 0.0
    if chance_sharing_parameters is not None:
        correction_cm = chance_sharing_correction_cm(
            float(chance_sharing_parameters["median_allele_frequency"]),
            int(chance_sharing_parameters["markers_on_chromosome"]),
            float(chance_sharing_parameters["chromosome_length_cm"]),
        )
    correlated: GammaEstimate | None = None
    independent: GammaEstimate | None = None
    if unit_count >= minimum_units_for_estimate:
        correlated = estimate_correlated_gamma(
            left_lengths, right_lengths, confidence_level=confidence_level,
            chance_sharing_correction_cm=correction_cm,
        )
        independent = estimate_independent_gamma(
            left_lengths, right_lengths, confidence_level=confidence_level,
            chance_sharing_correction_cm=correction_cm,
        )
    status = "NOT_ESTIMATED" if correlated is None else (
        "PRIMARY_ESTIMATE" if unit_count >= minimum_units_for_primary else "EXPLORATORY_ONLY"
    )
    correlated_status = "NOT_ESTIMATED" if correlated is None else (
        "PRIMARY" if status == "PRIMARY_ESTIMATE" else "EXPLORATORY"
    )
    exclusion = "" if correlated is not None else (
        "FOUNDER_NOT_SUPPORTED" if not founder_supported else "INSUFFICIENT_INDEPENDENT_UNITS"
    )
    estimate_rows = [
        _estimate_row(correlated, "CORRELATED", correlated_status, confidence_level, unit_count, exclusion),
        _estimate_row(independent, "INDEPENDENT", "EXPLORATORY" if independent else "NOT_ESTIMATED", confidence_level, unit_count, exclusion),
    ]
    scenarios = [
        _scenario_row("gamma_correlated", "MODEL", "CORRELATED", correlated, exclusion_code=exclusion, unit_count=unit_count),
        _scenario_row("gamma_independent", "MODEL", "INDEPENDENT", independent, exclusion_code=exclusion, unit_count=unit_count),
    ]
    if correlated is not None:
        for years in generation_years:
            scenarios.append(_scenario_row(
                f"generation_{_decimal(years).replace('.', '_')}",
                "GENERATION_INTERVAL", "CORRELATED", correlated,
                generation_years=float(years),
            ))
    if leave_one_family_out and unit_count - 1 >= minimum_units_for_estimate:
        for omitted_index, omitted in enumerate(included_source_rows):
            retained_left = left_lengths[:omitted_index] + left_lengths[omitted_index + 1:]
            retained_right = right_lengths[:omitted_index] + right_lengths[omitted_index + 1:]
            estimate = estimate_correlated_gamma(
                retained_left, retained_right, confidence_level=confidence_level,
                chance_sharing_correction_cm=correction_cm,
            )
            scenarios.append(_scenario_row(
                f"lofo_{omitted['INDEPENDENT_UNIT_ID']}", "LEAVE_ONE_FAMILY_OUT",
                "CORRELATED", estimate,
                omitted_unit_id=omitted["INDEPENDENT_UNIT_ID"],
            ))
    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / "variant_age_input.tsv"
    estimates_path = output_dir / "variant_age_estimates.tsv"
    scenarios_path = output_dir / "variant_age_scenarios.tsv"
    summary_path = output_dir / "variant_age_summary.json"
    _write_tsv(input_path, INPUT_COLUMNS, input_rows)
    _write_tsv(estimates_path, ESTIMATE_COLUMNS, estimate_rows)
    _write_tsv(scenarios_path, SCENARIO_COLUMNS, scenarios)
    summary = {
        "schema_version": "1.0.0", "method_id": "gamma_gandolfo_2014_v1",
        "status": status, "primary_model": "CORRELATED",
        "included_unit_count": unit_count,
        "minimum_units_for_estimate": minimum_units_for_estimate,
        "minimum_units_for_primary": minimum_units_for_primary,
        "estimate_generations": None if correlated is None else correlated.generations,
        "confidence_lower_generations": None if correlated is None else correlated.confidence_lower,
        "confidence_upper_generations": None if correlated is None else correlated.confidence_upper,
        "confidence_level": confidence_level,
        "chance_sharing_correction_applied": chance_sharing_parameters is not None,
        "manual_validation_required": True,
    }
    atomic_write_json(summary_path, summary)
    validate_tsv_table(input_path, "variant_age_input.schema.json")
    validate_tsv_table(estimates_path, "variant_age_estimates.schema.json")
    validate_tsv_table(scenarios_path, "variant_age_scenarios.schema.json")
    validate_json_document(summary, "variant_age_summary.schema.json")
    return VariantAgePublication(
        input_path=input_path, estimates_path=estimates_path,
        scenarios_path=scenarios_path, summary_path=summary_path, status=status,
        included_unit_count=unit_count, scenario_count=len(scenarios),
    )
