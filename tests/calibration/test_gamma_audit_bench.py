"""D25a : fidélité R, contrats d'observation et isolation de vérité."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.gamma_audit_bench import (
    _pairwise_shared_observation,
    _rho_region,
    load_configuration,
    materialize_bench,
)
from effet_fondateur.calibration.report_gamma_audit import build_report
from effet_fondateur.dating import (
    GammaAgeError,
    estimate_correlated_gamma,
    estimate_correlated_gamma_r_faithful,
    estimate_independent_gamma_r_faithful,
)


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = REPOSITORY / "calibration/configurations/d25a_gamma_audit_v1.json"


def test_d25a_configuration_uses_new_seeds_and_full_factorial_grid() -> None:
    document = load_configuration(CONFIGURATION)
    assert document["simulation_seed"] == 25101
    assert document["bootstrap_seed"] == 25102
    assert document["ages_generations"] == [5, 10, 15, 20, 30, 40, 50, 75, 100, 150]
    assert document["family_counts"] == [3, 5, 10]
    assert [row["shared_path_fraction"] for row in document["genealogy_profiles"]] == [0, .25, .5]
    assert document["calibrated"] is False
    assert document["founder_effect_proven"] is False


def test_d25a_materialization_separates_age_and_truncates_only_longest(tmp_path: Path) -> None:
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    root = materialize_bench(
        configuration_path=CONFIGURATION,
        workspace_root=workspace,
        run_id="d25a_materialization_test",
    )
    with gzip.open(root / "blinded_inputs/gamma_lengths.jsonl.gz", "rt", encoding="utf-8") as handle:
        blinded = json.loads(next(handle))
    with gzip.open(root / "private_truth/gamma_age_truth.jsonl.gz", "rt", encoding="utf-8") as handle:
        truth = json.loads(next(handle))
    assert "age_generations" not in blinded
    assert truth["age_generations"] == 5
    assert truth["case_id"] == blinded["case_id"]
    for side in ("left", "right"):
        full = blinded[f"full_{side}_lengths_cm"]
        observed = blinded[f"truncated_{side}_lengths_cm"]
        changed = [index for index, pair in enumerate(zip(full, observed)) if pair[0] != pair[1]]
        assert len(changed) == 1
        assert observed[changed[0]] == sorted(full)[-2]
    manifest = json.loads((root / "bench_manifest.json").read_text())
    assert manifest["case_count"] == 18_000
    assert manifest["truth_separated"] is True
    assert manifest["d25_seed_reused"] is False


def test_pairwise_observation_preserves_all_nonmaximum_arms() -> None:
    values = [1.0, 5.0, 2.0, 3.0]
    assert _pairwise_shared_observation(values) == [1.0, 3.0, 2.0, 3.0]
    assert values == [1.0, 5.0, 2.0, 3.0]


def test_complete_breakpoints_do_not_receive_boundary_or_chance_correction() -> None:
    left = [1.0, 2.0, 3.0]
    right = [2.0, 3.0, 4.0]
    complete = estimate_independent_gamma_r_faithful(
        left, right, apply_boundary_correction=False
    )
    observed = estimate_independent_gamma_r_faithful(
        left, right, apply_boundary_correction=True
    )
    expected = (2 * len(left) - 1) / ((sum(left) + sum(right)) / 100)
    assert complete.generations == pytest.approx(expected)
    assert observed.generations < complete.generations
    assert complete.boundary_correction_applied is False
    with pytest.raises(
        GammaAgeError, match="chance_correction_requires_observed_segments"
    ):
        estimate_independent_gamma_r_faithful(
            left, right, apply_boundary_correction=False,
            chance_sharing_correction_cm=.01,
        )


def test_r_faithful_i1_uses_distinct_point_and_interval_effective_sizes() -> None:
    left = [7.91375903, 3.01864831, 2.19277389, 1.11229976, 0.22605547]
    right = [4.00602268, 2.27065984, 0.76940051, 0.96899146, 3.25233799]
    candidate = estimate_correlated_gamma_r_faithful(left, right)
    legacy = estimate_correlated_gamma(left, right)
    assert _rho_region(float(candidate.correlation), len(left)) == "I1"
    assert candidate.point_effective_sample_size == pytest.approx(-2.9483625259242037)
    assert candidate.effective_sample_size == pytest.approx(1.3528663483273233)
    assert candidate.generations == pytest.approx(37.87869299187613)
    assert candidate.confidence_lower == pytest.approx(5.93679665903694)
    assert candidate.confidence_upper == pytest.approx(80.82598956044556)
    assert legacy.generations == pytest.approx(20.416863984966565)


@pytest.mark.parametrize(
    ("fixture_id", "expected_region"),
    [
        ("NEGATIVE_RHO_I1", "I1"),
        ("NEGATIVE_RHO_I2", "I2"),
        ("NEGATIVE_RHO_I3", "I3"),
        ("NONNEGATIVE_RHO", "NONNEGATIVE"),
    ],
)
def test_r_reference_fixtures_cover_every_rho_branch(
    fixture_id: str, expected_region: str,
) -> None:
    config = load_configuration(CONFIGURATION)
    fixture = next(
        row for row in config["r_reference_fixtures"]
        if row["fixture_id"] == fixture_id
    )
    estimate = estimate_correlated_gamma_r_faithful(
        fixture["left_lengths_cm"], fixture["right_lengths_cm"]
    )
    assert _rho_region(
        float(estimate.correlation), len(fixture["left_lengths_cm"])
    ) == expected_region


def test_d25a_report_is_detailed_and_immutable(tmp_path: Path) -> None:
    def metrics() -> dict[str, float | int]:
        return {
            "total_count": 10, "estimated_count": 10, "not_evaluable_count": 0,
            "bias": .1, "bias_ci95_low": -.2, "bias_ci95_high": .4,
            "mae": 2.0, "mae_ci95_low": 1.0, "mae_ci95_high": 3.0,
            "rmse": 3.0, "rmse_ci95_low": 2.0, "rmse_ci95_high": 4.0,
            "coverage_numerator": 9, "coverage_denominator": 10,
            "coverage": .95, "coverage_ci95_low": .9,
            "coverage_ci95_high": .98, "mean_interval_length_generations": 20.0,
        }

    genealogies = ["INDEPENDENT", "CORRELATED_25", "CORRELATED_50"]
    observations = ["FULL_BREAKPOINTS", "PAIRWISE_SHARED_TRUNCATED"]
    def primary(genealogy: str, observation: str) -> str:
        if genealogy == "INDEPENDENT":
            return "EXACT_INDEPENDENT_ORACLE_V1" if observation == "FULL_BREAKPOINTS" else "R_FAITHFUL_INDEPENDENT_V2"
        return "R_FAITHFUL_CORRELATED_V2"

    contracts = [{
        "genealogy_profile": genealogy, "observation_profile": observation,
        "method": primary(genealogy, observation), **metrics(),
        "lofo_evaluable_count": 5,
        "lofo_mean_max_absolute_delta_generations": 2.0,
        "lofo_p95_max_absolute_delta_generations": 5.0,
    } for genealogy in genealogies for observation in observations]
    methods_by_observation = {
        "FULL_BREAKPOINTS": [
            "LEGACY_INDEPENDENT_V1", "LEGACY_CORRELATED_V1",
            "EXACT_INDEPENDENT_ORACLE_V1", "R_FAITHFUL_CORRELATED_V2",
        ],
        "PAIRWISE_SHARED_TRUNCATED": [
            "LEGACY_INDEPENDENT_V1", "LEGACY_CORRELATED_V1",
            "R_FAITHFUL_INDEPENDENT_V2", "R_FAITHFUL_CORRELATED_V2",
        ],
    }
    method_results = [{
        "genealogy_profile": genealogy, "observation_profile": observation,
        "method": method, **metrics(),
    } for genealogy in genealogies for observation in observations
      for method in methods_by_observation[observation]]
    age_results = [{
        "genealogy_profile": genealogy, "observation_profile": observation,
        "method": primary(genealogy, observation), "age_generations": age,
        **metrics(),
    } for genealogy in genealogies for observation in observations
      for age in [5, 10, 15, 20, 30, 40, 50, 75, 100, 150]]
    family_results = [{
        "genealogy_profile": genealogy, "observation_profile": observation,
        "method": primary(genealogy, observation), "family_count": family,
        **metrics(),
    } for genealogy in genealogies for observation in observations
      for family in [3, 5, 10]]
    parity_fixtures = [{
        "fixture_id": name, "expected_rho_region": region,
        "observed_rho_region": region, "rho": .1,
        "maximum_absolute_error": 1e-14, "region_matches": True,
    } for name, region in [
        ("OFFICIAL_SIX", "NONNEGATIVE"), ("I1", "I1"), ("I2", "I2"),
        ("I3", "I3"), ("POSITIVE", "NONNEGATIVE"),
    ]]
    document = {
        "schema_version": "1.0.0", "report_id": "CAL-REP-14A",
        "method_id": "d25a_gamma_domain_audit_v1", "evaluation_state": "CLOSED",
        "configuration_sha256": "a" * 64, "protocol_sha256": "b" * 64,
        "bench_manifest_sha256": "c" * 64, "run_manifest_sha256": "d" * 64,
        "r_parity_sha256": "e" * 64, "technical_status": "PASS",
        "performance_qualification": "FAIL",
        "design": {
            "ages_generations": [5, 10, 15, 20, 30, 40, 50, 75, 100, 150],
            "family_counts": [3, 5, 10], "genealogy_profiles": genealogies,
            "observation_profiles": observations, "replications_per_cell": 200,
            "independent_scenario_count": 18000, "estimate_row_count": 144000,
            "confidence_level": .95,
        },
        "r_parity": {
            "schema_version": "1.0.0", "state": "CLOSED",
            "configuration_sha256": "a" * 64,
            "official_r_script_sha256": "f" * 64,
            "r_version": "R 4.4.1", "fixture_count": 5,
            "maximum_absolute_error_generations": 1e-14,
            "all_regions_match": True, "fixtures": parity_fixtures,
        },
        "contract_results": contracts, "method_results": method_results,
        "results_by_age": age_results, "results_by_family_count": family_results,
        "rho_region_counts": [{
            "genealogy_profile": genealogy, "observation_profile": observation,
            "total_count": 10,
            "counts": {"I1": 1, "I2": 2, "I3": 3, "NONNEGATIVE": 4},
        } for genealogy in genealogies for observation in observations],
        "criteria": [{
            "criterion_id": f"CAL-AGE-A0{index}", "metric": metric,
            "value": None if index == 6 else .1, "lower_threshold": None,
            "upper_threshold": None if index == 6 else 1.0,
            "status": "REVIEW_REQUIRED" if index == 6 else "PASS",
        } for index, metric in enumerate([
            "maximum_r_python_error", "absolute_oracle_bias",
            "oracle_coverage", "absolute_truncated_independent_bias",
            "truncated_independent_coverage", "correlated_performance",
        ], 1)],
        "phase_error_evaluation": {
            "status": "NOT_EVALUABLE", "reason": "MATCHED_PHASE_REQUIRED",
        },
        "calibrated": False, "founder_effect_proven": False,
    }
    evaluation = tmp_path / "evaluation"
    evaluation.mkdir()
    source = evaluation / "CAL-REP-14A.json"
    source.write_text(json.dumps(document), encoding="utf-8")
    (evaluation / "checksums.sha256").write_text(
        f"{sha256_file(source)}  CAL-REP-14A.json\n", encoding="utf-8"
    )
    template = tmp_path / "template.html"
    template.write_text(
        '<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>',
        encoding="utf-8",
    )
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    report = build_report(
        evaluation_root=evaluation, template_path=template,
        output_root=workspace / "d25a_report",
    )
    html = report.read_text(encoding="utf-8")
    assert html.count('<section id="s') == 14
    assert html.count('<figure class="scientific-figure">') == 4
    assert "Équation 6" in html and "Tableau 5" in html
    assert "calibrated=false" in html and "founder_effect_proven=false" in html
    assert "[To complete:" not in html
