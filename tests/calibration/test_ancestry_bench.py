"""Contrats ciblés de D27."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.ancestry_bench import (
    _align,
    _line_proportion,
    _local_sample_distances,
    load_configuration,
)
from effet_fondateur.calibration.report_ancestry import _render
from effet_fondateur.contracts import validate_json_document


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = REPOSITORY / "calibration/configurations/d27_ancestry_bench_v1.json"
SCHEMAS = REPOSITORY / "calibration/schemas"


def test_d27_configuration_preserves_unsealed_scientific_status() -> None:
    document = load_configuration(CONFIGURATION)
    assert document["reference_partition_role"] == "PHASING_REFERENCE_SHARED"
    assert document["reference_entity_count"] == 1403
    assert document["performance_thresholds_sealed"] is False
    assert document["production_promotion"] is False
    assert document["calibrated"] is False
    assert document["founder_effect_proven"] is False


def test_d27_line_projection_is_clipped_to_reference_segment() -> None:
    scores = np.asarray([[-1.0, 0.0], [0.5, 0.0], [2.0, 0.0]])
    proportions = _line_proportion(scores, np.asarray([1.0, 0.0]), np.asarray([0.0, 0.0]))
    assert proportions.tolist() == [0.0, 0.5, 1.0]


def test_d27_reference_only_procrustes_recovers_rotated_study() -> None:
    control_reference = np.asarray([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])
    rotation = np.asarray([[0.0, -1.0], [1.0, 0.0]])
    profile_reference = control_reference @ rotation
    control_study = np.asarray([[0.2, 0.7], [-0.3, 0.1]])
    profile_study = control_study @ rotation
    aligned = _align(
        profile_reference, profile_study, np.ones(2), control_reference
    )
    assert aligned == pytest.approx(control_study)


def test_d27_local_distance_ignores_copy_orientation_swap() -> None:
    ids = ["S1:H1", "S1:H2", "S2:H1", "S2:H2"]
    control = np.asarray([[0.0], [1.0], [2.0], [3.0]])
    profile = np.asarray([[1.0], [0.0], [2.0], [3.0]])
    distances, swaps = _local_sample_distances(ids, profile, control)
    assert distances.tolist() == [0.0, 0.0]
    assert swaps == 1


def test_d27_evaluation_schema_accepts_closed_document() -> None:
    document = {
        "schema_version": "1.0.0", "evaluation_state": "CLOSED",
        "report_id": "CAL-REP-16A", "method_id": "d27_ancestry_evaluation_v1",
        "configuration_sha256": "a" * 64, "protocol_sha256": "b" * 64,
        "calls_manifest_sha256": "c" * 64, "reference_role": "PHASING_REFERENCE_SHARED",
        "reference_entity_count": 1403, "analysis_scope": "REGIONAL",
        "stage_execution_scope": "CORE", "profiles": [str(index) for index in range(8)],
        "component_count": 10,
        "mixture": {"sample_count": 10, "mean_absolute_error": 0.1, "bootstrap_ci95_low": 0.05, "bootstrap_ci95_high": 0.15, "results_by_true_proportion": []},
        "stability_by_profile": [{
            "profile_id": str(index % 8), "scope": "GLOBAL", "entity_count": 10,
            "informative_variant_count": 100, "mean_standardized_distance": 0.1,
            "median_standardized_distance": 0.1, "p95_standardized_distance": 0.2,
            "maximum_standardized_distance": 0.3, "optimal_haplotype_swap_count": None,
        } for index in range(16)],
        "repeated_position_instability": 0.1,
        "local_founder_check": {"carrier_copy_count": 6, "independent_founder_haplotype_count": 1, "mean_founder_reference_axis_proportion": 0.5, "closer_to_founder_reference_count": 3},
        "criteria": [{"criterion_id": f"CAL-ANC-00{index}", "metric": "m", "value": 0.1, "threshold": None, "status": "REVIEW_REQUIRED"} for index in (1, 2)],
        "technical_status": "PASS_CORE_EVALUATOR", "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
        "production_promotion": False, "calibrated": False, "founder_effect_proven": False,
        "checks": ["PASS"], "limits": [f"limit {index}" for index in range(5)],
    }
    validate_json_document(document, "ancestry_evaluation.schema.json", schemas_dir=SCHEMAS)
    json.dumps(document)


def test_d27_protocol_refuses_geographic_assignment_and_retrofit() -> None:
    text = (REPOSITORY / "calibration/protocoles/D27_ANCESTRY_PRESPECIFICATION.md").read_text(encoding="utf-8")
    assert "proxy régional" in text
    assert "REVIEW_REQUIRED_UNSEALED" in text
    assert "pas une identité ethnique ou géographique" in text


def test_d27_report_is_detailed_and_keeps_nonqualification() -> None:
    evaluation = {
        "schema_version": "1.0.0", "evaluation_state": "CLOSED",
        "report_id": "CAL-REP-16A", "method_id": "d27_ancestry_evaluation_v1",
        "configuration_sha256": "a" * 64, "protocol_sha256": "b" * 64,
        "calls_manifest_sha256": "c" * 64, "reference_role": "PHASING_REFERENCE_SHARED",
        "reference_entity_count": 1403, "analysis_scope": "REGIONAL", "stage_execution_scope": "CORE",
        "profiles": ["control", "missing_01", "missing_05", "genotype_01", "phase_point_02", "phase_switch_02", "density_50", "combined"],
        "component_count": 10,
        "mixture": {"sample_count": 190, "mean_absolute_error": 0.16, "median_absolute_error": 0.15, "p95_absolute_error": 0.4, "bootstrap_ci95_low": 0.14, "bootstrap_ci95_high": 0.18, "results_by_true_proportion": [{"true_stratum_a_proportion": truth, "sample_count": 20, "estimated_mean": truth, "mean_absolute_error": 0.1} for truth in (0.0, 0.5, 1.0)]},
        "stability_by_profile": [{"profile_id": profile, "scope": scope, "entity_count": 908, "informative_variant_count": 244 if profile in {"density_50", "combined"} else 487, "mean_standardized_distance": 0.1, "median_standardized_distance": 0.08, "p95_standardized_distance": 0.2, "maximum_standardized_distance": 0.3, "optimal_haplotype_swap_count": 2 if scope == "LOCAL" else None} for scope in ("GLOBAL", "LOCAL") for profile in ("control", "missing_01", "missing_05", "genotype_01", "phase_point_02", "phase_switch_02", "density_50", "combined")],
        "repeated_position_instability": 0.6,
        "local_founder_check": {"carrier_copy_count": 6, "independent_founder_haplotype_count": 1, "mean_founder_reference_axis_proportion": 0.88, "closer_to_founder_reference_count": 6},
        "criteria": [{"criterion_id": f"CAL-ANC-00{index}", "metric": "metric", "value": 0.1, "threshold": None, "status": "REVIEW_REQUIRED"} for index in (1, 2)],
        "technical_status": "PASS_CORE_EVALUATOR", "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
        "production_promotion": False, "calibrated": False, "founder_effect_proven": False,
        "checks": ["PASS"], "limits": [f"limit {index}" for index in range(5)],
    }
    html = _render(evaluation, "f" * 64, "<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>")
    assert html.count('<section id="s') == 14
    assert html.count('<figure class="scientific-figure">') == 4
    assert "MAE" in html and "Procrustes" in html
    assert "calibrated=false" in html and "founder_effect_proven=false" in html
