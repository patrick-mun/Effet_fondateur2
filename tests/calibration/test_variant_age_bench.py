"""D25 : grille d’âge, séparation de vérité et métriques Gamma."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.variant_age_bench import (
    _bootstrap_summary,
    _estimate,
    _profile_lengths,
    load_configuration,
    materialize_bench,
)
from effet_fondateur.calibration.report_variant_age import build_report


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = REPOSITORY / "calibration/configurations/d25_variant_age_bench_v1.json"


def test_d25_configuration_closes_full_age_family_grid() -> None:
    document = load_configuration(CONFIGURATION)
    assert document["ages_generations"] == [5, 10, 15, 20, 30, 40, 50, 75, 100, 150]
    assert document["family_counts"] == [3, 5, 10]
    assert document["replications_per_cell"] == 200
    assert document["calibrated"] is False
    assert document["founder_effect_proven"] is False


def test_d25_materialization_keeps_age_out_of_blinded_inputs(tmp_path: Path) -> None:
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    root = materialize_bench(
        configuration_path=CONFIGURATION,
        workspace_root=workspace,
        run_id="d25_materialization_test",
    )
    with gzip.open(root / "blinded_inputs/age_lengths.jsonl.gz", "rt", encoding="utf-8") as handle:
        blinded = json.loads(next(handle))
    with gzip.open(root / "private_truth/age_truth.jsonl.gz", "rt", encoding="utf-8") as handle:
        truth = json.loads(next(handle))
    assert "age_generations" not in blinded
    assert truth["age_generations"] == 5
    assert truth["case_id"] == blinded["case_id"]
    manifest = json.loads((root / "bench_manifest.json").read_text())
    assert manifest["case_count"] == 6000
    assert manifest["truth_separated"] is True


def test_d25_profiles_change_only_declared_measurement_component() -> None:
    left, right = [1.0, 4.0], [2.0, 5.0]
    assert _profile_lengths(left, right, {"kind": "EXACT", "value": 1}) == (left, right)
    assert _profile_lengths(left, right, {"kind": "ARM_CAP_CM", "value": 3}) == ([1.0, 3.0], [2.0, 3.0])
    scaled = _profile_lengths(left, right, {"kind": "MAP_SCALE", "value": .95})
    assert scaled[0] == pytest.approx([.95, 3.8])
    assert scaled[1] == pytest.approx([1.9, 4.75])


def test_d25_gamma_and_bootstrap_metrics_are_finite() -> None:
    result = _estimate("CORRELATED", [1.2, 3.4, 1.6], [2.0, 3.8, 2.7], .95)
    assert result["status"] == "ESTIMATED"
    assert result["ci_lower_generations"] < result["estimate_generations"] < result["ci_upper_generations"]
    summary = _bootstrap_summary([-2.0, 0.0, 2.0], 100, 2502)
    assert summary["bias"] == pytest.approx(0)
    assert summary["mae"] == pytest.approx(4 / 3)
    assert summary["rmse"] == pytest.approx((8 / 3) ** .5)


def test_d25_report_is_detailed_and_immutable(tmp_path: Path) -> None:
    def metrics() -> dict[str, float | int]:
        return {
            "total_count": 10, "estimated_count": 10, "not_evaluable_count": 0,
            "bias": -2.0, "bias_ci95_low": -3.0, "bias_ci95_high": -1.0,
            "mae": 3.0, "mae_ci95_low": 2.0, "mae_ci95_high": 4.0,
            "rmse": 4.0, "rmse_ci95_low": 3.0, "rmse_ci95_high": 5.0,
            "coverage_numerator": 9, "coverage_denominator": 10,
            "coverage": .9, "coverage_ci95_low": .6, "coverage_ci95_high": .98,
        }

    document = {
        "schema_version": "1.0.0", "report_id": "CAL-REP-14",
        "method_id": "gamma_age_truth_bench_v1", "evaluation_state": "CLOSED",
        "configuration_sha256": "a" * 64, "protocol_sha256": "b" * 64,
        "bench_manifest_sha256": "c" * 64, "run_manifest_sha256": "d" * 64,
        "technical_status": "PASS", "performance_qualification": "REVIEW_REQUIRED",
        "design": {
            "ages_generations": [5, 10, 15, 20, 30, 40, 50, 75, 100, 150],
            "family_counts": [3, 5, 10], "replications_per_cell": 1,
            "independent_scenario_count": 30, "primary_model": "CORRELATED",
            "primary_profile": "ORACLE_EXACT", "confidence_level": .95,
        },
        "primary_metrics": metrics(),
        "results_by_age": [{"age_generations": age, **metrics()} for age in [5, 10, 15, 20, 30, 40, 50, 75, 100, 150]],
        "results_by_family_count": [{"family_count": count, **metrics()} for count in [3, 5, 10]],
        "profile_results": [
            {"profile": profile, "model": model, **metrics()}
            for profile in ["ORACLE_EXACT", "WINDOW_3CM_CAP", "MAP_SCALE_0P95", "MAP_SCALE_1P05"]
            for model in ["CORRELATED", "INDEPENDENT"]
        ],
        "d24_propagation_results": [{"profile": profile, **metrics()} for profile in ["MAXIMUM", "BASELINE", "DENSITY_50", "MISSING_1PCT"]],
        "leave_one_family_out": {
            "evaluable_scenario_count": 20,
            "mean_max_absolute_delta_generations": 2.0,
            "p95_max_absolute_delta_generations": 5.0,
            "three_family_status": "NOT_EVALUABLE_AFTER_ONE_OMISSION",
        },
        "historical_reference": {
            "correlated_estimate_generations": 22.34,
            "independent_estimate_generations": 27.49,
            "maximum_absolute_reference_error_generations": .04,
        },
        "phase_error_evaluation": {"status": "NOT_EVALUABLE", "reason": "MATCHED_PHASE_REQUIRED"},
        "criteria": [
            {"criterion_id": f"CAL-AGE-00{index}", "metric": f"metric_{index}",
             "value": .04, "threshold": .1 if index == 6 else None,
             "status": "PASS" if index == 6 else "REVIEW_REQUIRED"}
            for index in range(1, 7)
        ],
        "calibrated": False, "founder_effect_proven": False,
    }
    evaluation = tmp_path / "evaluation"
    evaluation.mkdir()
    source = evaluation / "CAL-REP-14.json"
    source.write_text(json.dumps(document), encoding="utf-8")
    (evaluation / "checksums.sha256").write_text(
        f"{sha256_file(source)}  CAL-REP-14.json\n", encoding="utf-8"
    )
    template = tmp_path / "template.html"
    template.write_text(
        '<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>',
        encoding="utf-8",
    )
    report_workspace = tmp_path / "calibration_workspace"
    report_workspace.mkdir()
    report = build_report(
        evaluation_root=evaluation, template_path=template,
        output_root=report_workspace / "d25_report",
    )
    html = report.read_text(encoding="utf-8")
    assert html.count('<section id="s') == 14
    assert html.count('<figure class="scientific-figure">') == 4
    assert "Équation 5" in html and "Tableau 5" in html
    assert "calibrated=false" in html and "founder_effect_proven=false" in html
    assert "[To complete:" not in html
