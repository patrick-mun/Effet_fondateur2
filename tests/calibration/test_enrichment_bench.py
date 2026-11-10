"""Contrats ciblés du banc D28."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.enrichment_bench import (
    _profile_result,
    _wilson,
    load_configuration,
)
from effet_fondateur.calibration import report_enrichment
from effet_fondateur.contracts import validate_json_document


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = REPOSITORY / "calibration/configurations/d28_enrichment_bench_v1.json"
SCHEMAS = REPOSITORY / "calibration/schemas"


def test_d28_configuration_keeps_production_threshold_null() -> None:
    configuration = load_configuration(CONFIGURATION)
    assert configuration["benchmark_decision_alpha"] == 0.05
    assert configuration["empirical_classification_threshold"] is None
    assert configuration["performance_thresholds_sealed"] is False
    assert configuration["calibrated"] is False
    assert configuration["founder_effect_proven"] is False


def test_d28_configuration_rejects_implicit_production_classification(tmp_path: Path) -> None:
    document = json.loads(CONFIGURATION.read_text(encoding="utf-8"))
    document["empirical_classification_threshold"] = 0.05
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(Exception):
        load_configuration(path)


def test_profile_metrics_keep_non_evaluations_in_denominators() -> None:
    rows = [
        {"profile": "MAXIMUM", "truth_role": "POSITIVE", "evaluation_status": "EVALUATED", "benchmark_positive_call": True, "empirical_probability": 0.01},
        {"profile": "MAXIMUM", "truth_role": "POSITIVE", "evaluation_status": "NOT_EVALUATED", "benchmark_positive_call": False, "empirical_probability": None},
        {"profile": "MAXIMUM", "truth_role": "NEGATIVE", "evaluation_status": "EVALUATED", "benchmark_positive_call": False, "empirical_probability": 0.6},
        {"profile": "MAXIMUM", "truth_role": "NEGATIVE", "evaluation_status": "NOT_EVALUATED", "benchmark_positive_call": False, "empirical_probability": None},
    ]
    result = _profile_result(rows, "MAXIMUM")
    assert result["power"] == 0.5
    assert result["null_false_positive_rate"] == 0.0
    assert result["evaluable_positive"] == result["evaluable_null"] == 1
    assert result["non_evaluable_count"] == 2


def test_wilson_interval_exposes_small_sample_uncertainty() -> None:
    low, high = _wilson(0, 5)
    assert low == 0.0
    assert high is not None and high > 0.4
    assert _wilson(0, 0) == (None, None)


def test_enrichment_evaluation_schema_forbids_global_promotion() -> None:
    profile = {
        "profile": "MAXIMUM", "positive_total": 5, "null_total": 5,
        "evaluable_positive": 5, "evaluable_null": 5, "positive_calls": 4,
        "null_calls": 1, "power": 0.8, "power_ci95_low": 0.3,
        "power_ci95_high": 0.98, "null_false_positive_rate": 0.2,
        "null_fpr_ci95_low": 0.04, "null_fpr_ci95_high": 0.62,
        "median_positive_probability": 0.01, "median_null_probability": 0.5,
        "non_evaluable_count": 0,
    }
    document = {
        "schema_version": "1.0.0", "evaluation_state": "CLOSED",
        "report_id": "CAL-REP-16B", "method_id": "d28_enrichment_evaluation_v1",
        "configuration_sha256": "a" * 64, "protocol_sha256": "b" * 64,
        "calls_manifest_sha256": "c" * 64, "analysis_scope": "synthetic",
        "stage_execution_scope": "core", "primary_profile": "MAXIMUM",
        "scenario_counts": {"independent_total": 10, "positive": 5, "null": 5, "profiled_case_total": 40},
        "benchmark_decision_alpha": 0.05, "production_classification_threshold": None,
        "primary_metrics": {"null_false_positive_rate": 0.2, "positive_scenario_power": 0.8, "classification_emitted_without_prespecified_threshold_count": 0},
        "profile_results": [{**profile, "profile": name} for name in ("MAXIMUM", "BASELINE", "DENSITY_50", "MISSING_1PCT")],
        "case_results": [{
            "scenario_label": f"P{index % 5 + 1}", "truth_role": "POSITIVE",
            "age_generations": 100, "profile": ("MAXIMUM", "BASELINE", "DENSITY_50", "MISSING_1PCT")[index // 10],
            "evaluation_status": "EVALUATED", "non_evaluable_reason": None,
            "observed_total_shared_cm": 1.0, "attempted_null_draws": 10000,
            "evaluable_null_draws": 10000, "non_evaluable_null_draws": 0,
            "exceedance_count": 1, "empirical_probability": 0.0002,
            "monte_carlo_interval_low": 0.0, "monte_carlo_interval_high": 0.001,
            "benchmark_positive_call": True, "production_status": "NOT_CLASSIFIED",
        } for index in range(40)],
        "non_evaluations": [],
        "criteria": [
            {"criterion_id": "CAL-ENRICH-001", "metric": "m1", "value": 0, "denominator": 40, "threshold": 0, "lifecycle": "fixed_contract", "status": "PASS"},
            {"criterion_id": "CAL-ENRICH-002", "metric": "m2", "value": 0.2, "denominator": 5, "threshold": None, "lifecycle": "pilot_then_seal", "status": "REVIEW_REQUIRED"},
            {"criterion_id": "CAL-ENRICH-003", "metric": "m3", "value": 0.8, "denominator": 5, "threshold": None, "lifecycle": "pilot_then_seal", "status": "REVIEW_REQUIRED"},
        ],
        "technical_status": "PASS_CORE_EVALUATOR",
        "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
        "production_promotion": False, "calibrated": False,
        "founder_effect_proven": False, "checks": ["PASS"],
        "limits": [f"limit {index}" for index in range(5)],
    }
    validate_json_document(document, "enrichment_evaluation.schema.json", schemas_dir=SCHEMAS)
    document["calibrated"] = True
    with pytest.raises(Exception):
        validate_json_document(document, "enrichment_evaluation.schema.json", schemas_dir=SCHEMAS)


def test_d28_protocol_separates_benchmark_and_production_thresholds() -> None:
    text = (REPOSITORY / "calibration/protocoles/D28_ENRICHMENT_PRESPECIFICATION.md").read_text(encoding="utf-8")
    assert "SEALED_BEFORE_TRUTH" not in text  # le statut est produit par l'appel, pas promis comme résultat
    assert "empirical_classification_threshold" in text
    assert "0,05" in text
    assert "calibrated=false" in text


def test_d28_report_has_fourteen_anchored_sections(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        report_enrichment, "_sections",
        lambda document, source_sha: [f"<p>section {index}</p>" for index in range(14)],
    )
    template = "<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>"
    html = report_enrichment._render({}, "f" * 64, template)
    assert html.count('<h2 id="s') == 14
    assert html.count('class="toc-entry"') == 14
    assert html.count('<footer class="page-footer">') == 1
    assert "CAL-REP-16B" in html
