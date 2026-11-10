"""Contrats ciblés du banc D26 LD/ROH."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from effet_fondateur.calibration.ld_roh_bench import (
    _exact_r2,
    _match_segments,
    _roh_result,
    _target_genotypes,
    load_configuration,
)
from effet_fondateur.calibration.report_ld_roh import _render
from effet_fondateur.contracts import validate_json_document


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = REPOSITORY / "calibration/configurations/d26_ld_roh_bench_v1.json"
SCHEMAS = REPOSITORY / "calibration/schemas"


def test_d26_configuration_keeps_scientific_thresholds_unsealed() -> None:
    document = load_configuration(CONFIGURATION)
    assert document["profiles"] == [
        "control", "missing_01", "missing_05", "genotype_01",
        "density_50", "combined",
    ]
    assert document["nonphased_aliases"] == {
        "phase_point_02": "control", "phase_switch_02": "control",
    }
    assert document["performance_thresholds_sealed"] is False
    assert document["production_promotion"] is False
    assert document["calibrated"] is False
    assert document["founder_effect_proven"] is False


def test_d26_exact_r2_uses_jointly_called_dosages() -> None:
    first = np.asarray([0.0, 1.0, 2.0, np.nan, 1.0])
    second = np.asarray([0.0, 1.0, 2.0, 2.0, 1.0])
    assert _exact_r2(first, second) == pytest.approx(1.0)
    assert _exact_r2(np.ones(4), np.arange(4.0)) is None


def test_d26_target_resolution_uses_molecular_position_for_dense_truth() -> None:
    rows = [{
        "variant_id": "SYNTHETIC_TARGET_D14B_CHR5", "position_bp": 84991179,
        "ref": "T", "alt": "A", "genotypes": ["0|0", "0|1", "1|1", ".|."],
    }]
    observed = _target_genotypes(
        ["s1", "s2", "s3", "s4"], rows, "S-3IFSF", 84991179,
    )
    assert observed == {"s1": "T/T", "s2": "T/A", "s3": "A/A", "s4": "./."}


def test_d26_segment_matching_is_same_sample_one_to_one() -> None:
    reference = [
        {"sample_id": "s1", "start_bp": 100, "end_bp": 200},
        {"sample_id": "s1", "start_bp": 300, "end_bp": 400},
    ]
    observed = [
        {"sample_id": "s1", "start_bp": 90, "end_bp": 410},
        {"sample_id": "s2", "start_bp": 100, "end_bp": 200},
    ]
    matches = _match_segments(observed, reference, 0.5)
    assert len(matches) == 1
    assert matches[0][0] == 0


def test_d26_roh_denominators_follow_profile_evaluated_samples() -> None:
    reference = [
        {"sample_id": "called_positive", "start_bp": 100, "end_bp": 200},
        {"sample_id": "excluded_positive", "start_bp": 100, "end_bp": 200},
    ]
    result = _roh_result(
        "missing_01", [], [reference[0]], ["called_positive", "called_negative"],
        [{
            "is_positive": True, "sample_id": "excluded_positive",
            "start_bp": 100, "end_bp": 200,
        }], 0.5, np.asarray([1, 1000]), np.asarray([0.0, 1.0]),
    )
    assert result["dense_reference_segment_count"] == 1
    assert result["reference_negative_sample_count"] == 1
    assert result["injected_positive_count"] == 0


def _common(report_id: str, method_id: str) -> dict[str, object]:
    return {
        "schema_version": "1.0.0", "evaluation_state": "CLOSED",
        "configuration_sha256": "a" * 64, "protocol_sha256": "b" * 64,
        "run_manifest_sha256": "c" * 64,
        "source_d16_manifest_sha256": "d" * 64,
        "source_d17_manifest_sha256": "e" * 64,
        "profiles": ["control", "missing_01", "missing_05", "genotype_01", "density_50", "combined"],
        "paired_profiles": True, "performance_thresholds_sealed": False,
        "technical_status": "PASS", "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
        "production_promotion": False, "calibrated": False,
        "founder_effect_proven": False, "report_id": report_id,
        "method_id": method_id,
    }


def test_d26_reports_are_detailed_and_schema_valid() -> None:
    profiles = ["control", "missing_01", "missing_05", "genotype_01", "density_50", "combined"]
    ld = {
        **_common("CAL-REP-15", "d26_local_ld_evaluation_v1"),
        "r_squared_definition": "GENOTYPE_DOSAGE_PEARSON_CORRELATION_SQUARED",
        "results_by_profile": [{
            "profile_id": profile, "evaluated_pair_count": 100,
            "reference_pair_count": 100, "missing_reference_pair_count": 0,
            "r_squared_mean_absolute_error": 0.001,
            "r_squared_median_absolute_error": 0.0005,
            "r_squared_p95_absolute_error": 0.002,
            "r_squared_maximum_absolute_error": 0.01,
        } for profile in profiles],
        "criteria": [
            {"criterion_id": "CAL-LD-001", "metric": "r_squared_mean_absolute_error", "value": 0.001, "threshold": None, "status": "REVIEW_REQUIRED"},
            {"criterion_id": "CAL-LD-002", "metric": "target_chromosome_in_background_count", "value": 0, "threshold": 0, "status": "PASS"},
        ],
        "limits": ["limite scientifique une", "limite scientifique deux", "limite scientifique trois"],
    }
    roh = {
        **_common("CAL-REP-16", "d26_roh_evaluation_v1"),
        "reference_role": "DENSE_D16_PLINK_OPERATIONAL_REFERENCE",
        "dense_reference_segment_count": 22,
        "results_by_profile": [{
            "profile_id": profile, "observed_segment_count": 25,
            "dense_reference_segment_count": 22, "matched_segment_count": 20,
            "roh_detection_sensitivity": 20 / 22,
            "sensitivity_ci95_low": 0.72, "sensitivity_ci95_high": 0.98,
            "reference_negative_sample_count": 891,
            "false_positive_sample_count": 35, "roh_false_positive_rate": 35 / 891,
            "false_positive_rate_ci95_low": 0.028, "false_positive_rate_ci95_high": 0.054,
            "unmatched_observed_segment_fraction": 0.2,
            "roh_length_mean_absolute_error_cm": 0.2,
            "roh_length_median_absolute_error_cm": 0.1,
            "injected_positive_count": 3, "injected_positive_detected_count": 3,
            "injected_positive_sensitivity": 1.0,
        } for profile in profiles],
        "criteria": [
            {"criterion_id": "CAL-ROH-001", "metric": "roh_detection_sensitivity", "value": 20 / 22, "threshold": None, "status": "REVIEW_REQUIRED"},
            {"criterion_id": "CAL-ROH-002", "metric": "roh_false_positive_rate", "value": 35 / 891, "threshold": None, "status": "REVIEW_REQUIRED"},
            {"criterion_id": "CAL-ROH-003", "metric": "roh_length_mean_absolute_error", "value": 0.2, "threshold": None, "status": "REVIEW_REQUIRED"},
            {"criterion_id": "CAL-ROH-004", "metric": "roh_used_as_direct_gamma_length_count", "value": 0, "threshold": 0, "status": "PASS"},
        ],
        "limits": ["limite scientifique une", "limite scientifique deux", "limite scientifique trois"],
    }
    template = "<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>"
    for document, figures in ((ld, 3), (roh, 4)):
        validate_json_document(document, "ld_roh_evaluation.schema.json", schemas_dir=SCHEMAS)
        html = _render(document, "f" * 64, template)
        assert html.count('<section id="s') == 14
        assert html.count('<figure class="scientific-figure">') == figures
        assert "calibrated=false" in html
        assert "founder_effect_proven=false" in html


def test_d26_protocol_explicitly_refuses_retrofitted_pass() -> None:
    text = (REPOSITORY / "calibration/protocoles/D26_LD_ROH_PRESPECIFICATION.md").read_text(encoding="utf-8")
    assert "REVIEW_REQUIRED_UNSEALED" in text
    assert "référence **opérationnelle**" in text
    assert "calibrated" in text and "founder_effect_proven" in text
