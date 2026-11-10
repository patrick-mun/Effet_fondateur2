"""D24 : contrats, oracle indépendant, évaluateur clos et rendu CAL-REP-13."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.evaluate_founder_haplotype import (
    evaluate_founder_haplotype,
    load_configuration,
)
from effet_fondateur.calibration.founder_haplotype_oracle import (
    acceptance_results,
    aggregate_results,
    evaluate_case,
    wilson_interval,
)
from effet_fondateur.calibration.report_founder_haplotype import build_report
from effet_fondateur.contracts import validate_json_document


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIGURATION = REPOSITORY / "calibration/configurations/d24_founder_haplotype_bench_v1.json"
SCHEMAS = REPOSITORY / "calibration/schemas"


def _write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def _case(
    workspace: Path,
    scenario: dict[str, object],
    profile: str,
    *,
    detected: bool | None = None,
) -> Path:
    positive = scenario["truth_role"] == "POSITIVE"
    detected = positive if detected is None else detected
    root = workspace / "cases" / str(scenario["scenario_id"]) / profile
    truth_root, stage = root / "private_truth", root / "stage13"
    truth_root.mkdir(parents=True)
    stage.mkdir(parents=True)
    families = [f"family_{number}" for number in range(1, 4)]
    truth = {
        "schema_version": "1.0.0",
        "scenario_id": scenario["scenario_id"],
        "profile": profile,
        "truth_role": "POSITIVE" if positive else "NEGATIVE",
        "detectable": positive,
        "age_generations": scenario["age_generations"],
        "independent_unit_count": 3,
        "target": {
            "variant_id": "target",
            "chromosome": 5,
            "position_bp": 100,
            "position_cm": 1.0,
        },
        "source_manifest_sha256": f"{int(scenario['random_seed']):064x}",
        "segments": [
            {
                "family_id": family,
                "left_boundary_bp": 90,
                "right_boundary_bp": 110,
                "left_boundary_cm": 0.9,
                "right_boundary_cm": 1.1,
                "left_censored": False,
                "right_censored": False,
            }
            for family in families
        ] if positive else [],
    }
    (truth_root / "founder_haplotype_case.json").write_text(json.dumps(truth), encoding="utf-8")
    status = "SUPPORTED_IBS_CANDIDATE" if detected else "INSUFFICIENT_INFORMATIVE_MARKERS"
    interpretation = "IBS_SHARED_CANDIDATE" if detected else "NO_FOUNDER_CONCLUSION"
    summary = {
        "schema_version": "1.0.0",
        "method_id": "target_centered_exact_ibs_v1",
        "interpretation": interpretation,
        "status": status,
        "selected_carrier_count": 3,
        "excluded_carrier_count": 0,
        "background_haplotype_count": 6,
        "matching_background_haplotype_count": 0,
        "minimum_independent_carriers": 3,
        "minimum_flank_markers": 2,
        "ibd_claimed": False,
    }
    (stage / "founder_analysis_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    consensus_columns = [
        "ANALYSIS_ID", "TARGET_VARIANT_ID", "INTERPRETATION", "CONSENSUS_STATUS",
        "CARRIER_UNIT_COUNT", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP",
        "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "LEFT_MARKER_COUNT", "RIGHT_MARKER_COUNT",
        "SIGNATURE_MARKER_COUNT", "BACKGROUND_COHORT", "BACKGROUND_HAPLOTYPE_COUNT",
        "MATCHING_BACKGROUND_HAPLOTYPE_COUNT", "BACKGROUND_MATCH_FREQUENCY",
        "SEGMENT_METHOD", "MAP_VERSION",
    ]
    _write_tsv(stage / "founder_consensus.tsv", consensus_columns, [[
        "primary_exact_ibs", "target", interpretation, status, "3", "90", "100", "110",
        "0.1", "0.1", "2", "2", "4", "controls_unrelated_plus_noncarriers", "6", "0",
        "0", "target_centered_exact_ibs_v1", "a" * 64,
    ]])
    segment_columns = [
        "INDEPENDENT_UNIT_ID", "SAMPLE_ID", "FAMILY_ID", "CARRIER_HAPLOTYPE_ID",
        "TARGET_VARIANT_ID", "LEFT_BOUND_BP", "TARGET_BP", "RIGHT_BOUND_BP",
        "LEFT_LENGTH_CM", "RIGHT_LENGTH_CM", "PHASING_CONFIDENCE", "SEGMENT_METHOD",
        "SEGMENT_STATUS", "EXCLUSION_CODE",
    ]
    _write_tsv(stage / "founder_segments.tsv", segment_columns, [
        [
            f"unit_{number}", f"sample_{number}", family, "H1", "target", "90", "100", "110",
            "0.1", "0.1", "0.99", "target_centered_pairwise_max_ibs_v1",
            "INCLUDED" if detected else "INSUFFICIENT", "" if detected else "NO_SHARED_FLANK",
        ]
        for number, family in enumerate(families, 1)
    ])
    return root


def test_d24_configurations_and_schemas_are_closed() -> None:
    configuration = load_configuration(CONFIGURATION)
    assert configuration["primary_profile"] == "MAXIMUM"
    for name, schema in (
        ("d24_founder_chr5_pedigree_v1.json", "pedigree_config.schema.json"),
        ("d24_founder_chr5_age100_positive_v1.json", "founder_simulation_config.schema.json"),
    ):
        document = json.loads((REPOSITORY / "calibration/configurations" / name).read_text())
        validate_json_document(document, schema, schemas_dir=SCHEMAS)


def test_oracle_uses_scenarios_not_markers_and_keeps_ibs_distinct_from_ibd() -> None:
    truth = {
        "truth_role": "POSITIVE", "detectable": True, "independent_unit_count": 3,
        "scenario_id": "D24_POS_TEST", "profile": "MAXIMUM", "age_generations": 100,
        "source_manifest_sha256": "a" * 64,
        "target": {"variant_id": "target", "position_bp": 100, "position_cm": 1.0},
        "segments": [
            {"family_id": family, "left_boundary_bp": 90, "right_boundary_bp": 110,
             "left_boundary_cm": .9, "right_boundary_cm": 1.1,
             "left_censored": False, "right_censored": False}
            for family in ("F1", "F2", "F3")
        ],
    }
    summary = {
        "status": "SUPPORTED_IBS_CANDIDATE",
        "interpretation": "IBS_SHARED_CANDIDATE",
        "ibd_claimed": False,
        "minimum_independent_carriers": 3,
        "minimum_flank_markers": 2,
    }
    consensus = {"CONSENSUS_STATUS": "SUPPORTED_IBS_CANDIDATE", "INTERPRETATION": "IBS_SHARED_CANDIDATE",
                 "TARGET_VARIANT_ID": "target", "TARGET_BP": "100",
                 "LEFT_BOUND_BP": "88", "RIGHT_BOUND_BP": "112", "LEFT_LENGTH_CM": ".12", "RIGHT_LENGTH_CM": ".12"}
    segments = [
        {"FAMILY_ID": family, "SEGMENT_STATUS": "INCLUDED", "TARGET_VARIANT_ID": "target", "TARGET_BP": "100",
         "LEFT_BOUND_BP": "88", "RIGHT_BOUND_BP": "112",
         "LEFT_LENGTH_CM": ".12", "RIGHT_LENGTH_CM": ".12"}
        for family in ("F1", "F2", "F3")
    ]
    positive = evaluate_case(truth, summary, consensus, segments, minimum_independent_units=3)
    negative = {**positive, "truth_positive": False, "detectable": False, "expected_detection": False, "detected": False,
                "target_contained": None, "missing_family_count": 0, "errors": {key: [] for key in positive["errors"]}}
    aggregate = aggregate_results([positive, negative])
    assert aggregate["metrics"]["founder_segment_detection_sensitivity"]["denominator"] == 1
    assert aggregate["metrics"]["founder_segment_detection_specificity"]["estimate"] == 1
    assert aggregate["metrics"]["boundary_errors"]["boundary_mean_absolute_error"] == pytest.approx(.02)
    assert acceptance_results(aggregate["metrics"])[-1]["status"] == "PASS"
    assert 0 < wilson_interval(1, 1)[0] < 1


def test_evaluator_and_report_publish_only_aggregates(tmp_path: Path) -> None:
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    configuration = load_configuration(CONFIGURATION)
    profiles = [configuration["primary_profile"], *configuration["sensitivity_profiles"]]
    cases = [
        _case(
            workspace,
            scenario,
            profile,
            detected=(
                scenario["truth_role"] == "POSITIVE" and profile != "BASELINE"
            ),
        )
        for scenario in configuration["scenarios"]
        for profile in profiles
    ]
    document = evaluate_founder_haplotype(
        case_roots=cases,
        configuration_path=CONFIGURATION,
        workspace_root=workspace,
        run_id="d24_evaluation_test",
        created_at="2026-11-10T00:00:00+00:00",
    )
    assert document["technical_status"] == "PASS"
    assert document["performance_qualification"] == "REVIEW_REQUIRED"
    assert document["scenario_counts"] == {
        "total": 10,
        "positive": 5,
        "negative": 5,
        "detectable_positive": 5,
        "profiled_case_total": 40,
    }
    assert len(document["profile_results"]) == 4
    profiles_by_name = {row["profile"]: row for row in document["profile_results"]}
    assert document["metrics"]["founder_segment_detection_sensitivity"]["estimate"] == 1
    assert profiles_by_name["BASELINE"]["sensitivity"] == 0
    serialized = json.dumps(document)
    assert "family_0" not in serialized
    assert "sample_0" not in serialized
    assert "scenario_0" not in serialized
    template = tmp_path / "template.html"
    template.write_text('<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>', encoding="utf-8")
    report = build_report(
        evaluation_root=workspace / "d24_evaluation_test",
        template_path=template,
        output_root=workspace / "d24_reports_test",
    )
    html = report.read_text(encoding="utf-8")
    assert html.count('<section id="s') == 14
    assert html.count('class="toc-entry"') == 14
    assert html.count('<figure class="scientific-figure">') == 4
    assert "Équation 4" in html and "Tableau 3" in html
    assert "calibrated=false" in html and "founder_effect_proven=false" in html
    assert "IBS" in html and "preuve IBD" in html
    assert "phasage-vérité" in html and "aucun phasage SHAPEIT5" in html
    assert "quatre faux positifs parmi cinq négatifs" in html
    assert "[To complete:" not in html
    with pytest.raises(ValueError, match="run_exists"):
        build_report(
            evaluation_root=workspace / "d24_evaluation_test",
            template_path=template,
            output_root=workspace / "d24_reports_test",
        )


def test_current_two_family_panel_is_refused_by_d24_oracle() -> None:
    truth = {
        "truth_role": "POSITIVE", "detectable": True, "independent_unit_count": 2,
        "profile": "MAXIMUM", "target": {"position_bp": 100, "position_cm": 1.0},
        "segments": [
            {"family_id": family, "left_boundary_bp": 90, "right_boundary_bp": 110,
             "left_boundary_cm": .9, "right_boundary_cm": 1.1}
            for family in ("F1", "F2")
        ],
    }
    with pytest.raises(ValueError, match="independent_units_insufficient"):
        evaluate_case(truth, {}, {}, [], minimum_independent_units=3)
