"""D25b : censure de type II et isolation du banc Gamma."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.gamma_censoring_bench import (
    load_configuration,
    materialize_bench,
    pairwise_type2_observation,
)
from effet_fondateur.calibration.report_gamma_censoring import build_report
from effet_fondateur.dating import (
    GammaAgeError,
    estimate_independent_gamma_type2_censored,
)


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = (
    REPOSITORY / "calibration/configurations/d25b_gamma_censoring_v1.json"
)


def test_d25b_configuration_is_new_and_non_promoting() -> None:
    document = load_configuration(CONFIGURATION)
    assert document["simulation_seed"] == 25201
    assert document["bootstrap_seed"] == 25202
    assert document["replications_per_cell"] == 400
    assert document["production_promotion"] is False
    assert document["calibrated"] is False
    assert document["founder_effect_proven"] is False


def test_type2_analytic_fixture_matches_closed_form() -> None:
    estimate = estimate_independent_gamma_type2_censored(
        [1.0, 2.0, 2.0], [3.0, 4.0, 4.0]
    )
    assert estimate.generations == pytest.approx(18.75)
    assert estimate.model == "INDEPENDENT_TYPE2_CENSORED"
    assert estimate.sample_size == 3
    assert estimate.effective_sample_size == 2
    assert estimate.boundary_correction_applied is False
    assert estimate.confidence_lower == pytest.approx(6.81165858516453)
    assert estimate.confidence_upper == pytest.approx(54.795456685889505)


def test_type2_rejects_data_without_repeated_censoring_threshold() -> None:
    with pytest.raises(
        GammaAgeError, match="gamma_type2_censored_maximum_not_repeated"
    ):
        estimate_independent_gamma_type2_censored(
            [1.0, 2.0, 3.0], [2.0, 3.0, 4.0]
        )


def test_type2_rejects_unmodelled_corrections() -> None:
    with pytest.raises(
        GammaAgeError, match="gamma_type2_chance_sharing_not_modelled"
    ):
        estimate_independent_gamma_type2_censored(
            [1.0, 2.0, 2.0], [3.0, 4.0, 4.0],
            chance_sharing_correction_cm=0.01,
        )
    with pytest.raises(
        GammaAgeError, match="gamma_type2_requires_one_censored_arm_per_side"
    ):
        estimate_independent_gamma_type2_censored(
            [1.0, 2.0, 2.0], [3.0, 4.0, 4.0],
            right_censored_largest_arms_per_side=2,
        )


def test_pairwise_observation_changes_only_unique_maximum() -> None:
    source = [1.0, 5.0, 2.0, 3.0]
    assert pairwise_type2_observation(source) == [1.0, 3.0, 2.0, 3.0]
    assert source == [1.0, 5.0, 2.0, 3.0]


def test_d25b_materialization_keeps_age_out_of_blinded_input(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    root = materialize_bench(
        configuration_path=CONFIGURATION,
        workspace_root=workspace,
        run_id="d25b_materialization_test",
    )
    with gzip.open(
        root / "blinded_inputs/gamma_lengths.jsonl.gz", "rt", encoding="utf-8"
    ) as handle:
        blinded = json.loads(next(handle))
    with gzip.open(
        root / "private_truth/gamma_age_truth.jsonl.gz", "rt", encoding="utf-8"
    ) as handle:
        truth = json.loads(next(handle))
    assert "age_generations" not in blinded
    assert truth["age_generations"] == 5
    assert truth["case_id"] == blinded["case_id"]
    for side in ("left", "right"):
        observed = blinded[f"observed_{side}_lengths_cm"]
        assert observed.count(max(observed)) >= 2
    manifest = json.loads((root / "bench_manifest.json").read_text())
    assert manifest["case_count"] == 36_000
    assert manifest["truth_separated"] is True
    assert manifest["prior_panel_seed_reused"] is False


def test_d25b_report_is_detailed_and_immutable(tmp_path: Path) -> None:
    def metrics() -> dict[str, float | int]:
        return {
            "total_count": 10, "estimated_count": 10,
            "not_evaluable_count": 0, "bias": 0.1,
            "bias_ci95_low": -0.2, "bias_ci95_high": 0.4,
            "mae": 2.0, "mae_ci95_low": 1.0, "mae_ci95_high": 3.0,
            "rmse": 3.0, "rmse_ci95_low": 2.0, "rmse_ci95_high": 4.0,
            "coverage_numerator": 9, "coverage_denominator": 10,
            "coverage": 0.95, "coverage_ci95_low": 0.90,
            "coverage_ci95_high": 0.98,
            "mean_interval_length_generations": 20.0,
        }

    genealogies = ["INDEPENDENT", "CORRELATED_25", "CORRELATED_50"]
    methods = [
        "TYPE2_CENSORED_INDEPENDENT_V3", "R_FAITHFUL_INDEPENDENT_V2",
        "R_FAITHFUL_CORRELATED_V2", "FULL_INDEPENDENT_ORACLE_V1",
    ]
    method_results = [{
        "genealogy_profile": genealogy, "method": method, **metrics(),
    } for genealogy in genealogies for method in methods]
    age_results = [{
        "genealogy_profile": genealogy,
        "method": "TYPE2_CENSORED_INDEPENDENT_V3",
        "age_generations": age, **metrics(),
    } for genealogy in genealogies
      for age in [5, 10, 15, 20, 30, 40, 50, 75, 100, 150]]
    family_results = [{
        "genealogy_profile": genealogy,
        "method": "TYPE2_CENSORED_INDEPENDENT_V3",
        "family_count": family_count, **metrics(),
    } for genealogy in genealogies for family_count in [3, 5, 10]]
    criteria = [{
        "criterion_id": f"CAL-AGE-B0{index}", "metric": metric,
        "value": None if index == 7 else 0.1,
        "lower_threshold": None,
        "upper_threshold": None if index == 7 else 1.0,
        "status": "REVIEW_REQUIRED" if index == 7 else "PASS",
    } for index, metric in enumerate([
        "analytic_fixture_absolute_error", "absolute_primary_bias",
        "primary_coverage", "maximum_absolute_family_bias",
        "maximum_family_coverage_deviation_from_0_95",
        "primary_not_estimated_count", "correlated_stress_performance",
    ], 1)]
    document = {
        "schema_version": "1.0.0", "report_id": "CAL-REP-14B",
        "method_id": "d25b_gamma_type2_censoring_evaluation_v1",
        "evaluation_state": "CLOSED", "configuration_sha256": "a" * 64,
        "protocol_sha256": "b" * 64, "bench_manifest_sha256": "c" * 64,
        "run_manifest_sha256": "d" * 64, "technical_status": "PASS",
        "performance_qualification": "PASS_INDEPENDENT_ONLY",
        "independent_contract_qualification": "PASS",
        "correlated_stress_qualification": "REVIEW_REQUIRED",
        "production_promotion": False,
        "design": {
            "ages_generations": [5, 10, 15, 20, 30, 40, 50, 75, 100, 150],
            "family_counts": [3, 5, 10], "genealogy_profiles": genealogies,
            "replications_per_cell": 400, "scenario_count": 36000,
            "estimate_row_count": 144000, "confidence_level": 0.95,
            "observation_contract": "PAIRWISE_SHARED_TYPE2_CENSORED",
        },
        "analytic_fixture": {
            "event_count": 4, "total_exposure_morgans": 0.16,
            "expected_generations": 18.75, "observed_generations": 18.75,
            "absolute_error_generations": 0.0,
        },
        "primary_result": {
            "genealogy_profile": "INDEPENDENT",
            "method": "TYPE2_CENSORED_INDEPENDENT_V3", **metrics(),
            "lofo_evaluable_count": 5,
            "lofo_mean_max_absolute_delta_generations": 2.0,
            "lofo_p95_max_absolute_delta_generations": 5.0,
        },
        "method_results": method_results, "results_by_age": age_results,
        "results_by_family_count": family_results, "criteria": criteria,
        "limits": {
            "phase_error": "NOT_EVALUABLE",
            "marker_resolution": "NOT_EVALUABLE",
            "chance_sharing": "NOT_MODELLED",
            "d24_false_positives": "NOT_ADDRESSED",
        },
        "calibrated": False, "founder_effect_proven": False,
    }
    evaluation = tmp_path / "evaluation"
    evaluation.mkdir()
    source = evaluation / "CAL-REP-14B.json"
    source.write_text(json.dumps(document), encoding="utf-8")
    (evaluation / "checksums.sha256").write_text(
        f"{sha256_file(source)}  CAL-REP-14B.json\n", encoding="utf-8"
    )
    template = tmp_path / "template.html"
    template.write_text(
        "<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>",
        encoding="utf-8",
    )
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    report = build_report(
        evaluation_root=evaluation, template_path=template,
        output_root=workspace / "d25b_report",
    )
    html = report.read_text(encoding="utf-8")
    assert html.count('<section id="s') == 14
    assert html.count('<figure class="scientific-figure">') == 4
    assert "Équation 6" in html and "Tableau 5" in html
    assert "calibrated=false" in html
    assert "founder_effect_proven=false" in html
    assert "[To complete:" not in html
