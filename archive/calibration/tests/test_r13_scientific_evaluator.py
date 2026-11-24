from __future__ import annotations

import copy
from pathlib import Path

import pytest

from effet_fondateur.calibration.contracts import STAGE_IDS, load_yaml
from effet_fondateur.calibration.r13_scientific_evaluator import (
    R13ScientificEvaluatorError,
    assert_scientific_evaluator_ready,
    evaluate_closed_metric_document,
    evaluate_metric_observation,
    evaluate_stage_observations,
    load_metric_rule_bundle,
    validate_metric_rule_bundle,
    validate_metric_rule_manifest,
    validate_metric_rule_subset,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_RULE_DRAFT = (
    REPOSITORY_ROOT
    / "calibration/configurations/scientific_metric_rules_phasing_to_sensitivity_v1.draft.yaml"
)
SCIENTIFIC_RULE_DRAFT_STAGES = (
    "12",
    "12A",
    "13",
    "14",
    "15",
    "16",
    "16A",
    "16B",
    "16C",
    "17",
)
SCIENTIFIC_COMPLEMENT_DRAFT = (
    REPOSITORY_ROOT
    / "calibration/configurations/"
    "scientific_metric_rules_integrity_to_mapping_and_reporting_v1.draft.yaml"
)
SCIENTIFIC_COMPLEMENT_DRAFT_STAGES = (
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "10A",
    "11",
    "18",
    "19",
)


def _passing_context(condition: dict[str, object]) -> dict[str, object]:
    if "any_of" in condition:
        return _passing_context(condition["any_of"][0])  # type: ignore[index,arg-type]
    if "all_of" in condition:
        merged: dict[str, object] = {}
        for branch in condition["all_of"]:  # type: ignore[union-attr]
            merged.update(_passing_context(branch))
        return merged
    observed = {}
    for key, expected in condition.items():
        observed_key = key[:-4] if key.endswith(("_min", "_max")) else key
        observed[observed_key] = expected[0] if isinstance(expected, list) else expected
    return observed


def _passing_observation(
    rule: dict[str, object],
    context: dict[str, object] | None = None,
) -> dict[str, object]:
    selected = context or rule
    operator = selected["operator"]
    threshold = float(selected["value"])
    estimate = (
        threshold + max(abs(threshold) * 0.1, 0.01)
        if operator == "min"
        else threshold - max(abs(threshold) * 0.1, 0.01)
        if operator == "max"
        else 0.0
        if operator == "equivalence"
        else threshold
    )
    observation: dict[str, object] = {
        "metric_id": rule["metric_id"],
        "independent_n": rule["minimum_independent_n"],
        "integrity_passed": True,
        "estimate": estimate,
    }
    if context is not None:
        observation["context_id"] = context["context_id"]
        observation["context"] = _passing_context(context["when"])  # type: ignore[arg-type]
    decision_value = selected["decision_value"]
    if decision_value == "lcb95":
        observation["lcb95"] = threshold
    elif decision_value == "ucb95":
        observation["ucb95"] = threshold
    elif decision_value == "ci95":
        observation["ci95"] = [-threshold, threshold]
    return observation


def _passing_closed_metric_document() -> dict[str, object]:
    bundle = load_metric_rule_bundle(REPOSITORY_ROOT)
    stages = {}
    for stage_id in STAGE_IDS:
        observations = []
        for rule in bundle["stages"][stage_id]["rules"]:
            if rule["operator"] == "contextual":
                observations.extend(
                    _passing_observation(rule, context) for context in rule["contexts"]
                )
            else:
                observations.append(_passing_observation(rule))
        stages[stage_id] = {"observations": observations}
    return {
        "schema_version": "1.0.0",
        "release_id": "panel-master-generic-v1-r12c-scientific-validation",
        "scope": "R13_SCIENTIFIC_POST_CLOSURE_METRICS",
        "run_state": "CLOSED",
        "truth_access_state": "AUTHORIZED_POST_CLOSURE",
        "metric_producers_frozen": True,
        "closed_run_sha256": "a" * 64,
        "metric_producer_bundle_sha256": "b" * 64,
        "stages": stages,
    }


def _bundle() -> tuple[dict[str, object], dict[str, object]]:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    stages = {}
    for stage_id in STAGE_IDS:
        stages[stage_id] = {
            "rules": [
                {
                    "metric_id": metric,
                    "operator": "exact",
                    "value": 1,
                    "decision_value": "estimate",
                    "unit": "contract_test_only",
                    "independent_unit": "contract_test_unit",
                    "minimum_independent_n": 1,
                }
                for metric in criteria["stages"][stage_id]["primary_metrics"]
            ]
        }
    return (
        {
            "schema_version": "1.0.0",
            "release_id": "panel-master-generic-v1-r12c-scientific-validation",
            "scope": "R13_SCIENTIFIC_BLIND_VALIDATION",
            "criteria_version": "1.0.0",
            "stages": stages,
        },
        criteria,
    )


def test_scientific_evaluator_requires_complete_primary_metric_coverage() -> None:
    bundle, criteria = _bundle()
    audit = validate_metric_rule_bundle(bundle, criteria)

    assert audit["stage_count"] == 25
    assert audit["rule_count"] > 100
    assert audit["scientific_evaluation_authorized"] is False


def test_scientific_evaluator_rejects_missing_primary_metric() -> None:
    bundle, criteria = _bundle()
    broken = copy.deepcopy(bundle)
    broken["stages"]["12"]["rules"].pop()  # type: ignore[index]

    with pytest.raises(R13ScientificEvaluatorError, match="coverage_mismatch"):
        validate_metric_rule_bundle(broken, criteria)


def test_scientific_evaluator_is_ready_after_full_metric_binding() -> None:
    assert_scientific_evaluator_ready()


def test_phasing_to_sensitivity_rule_draft_matches_r08_exactly() -> None:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    draft = load_yaml(SCIENTIFIC_RULE_DRAFT)

    audit = validate_metric_rule_subset(
        draft,
        criteria,
        SCIENTIFIC_RULE_DRAFT_STAGES,
    )

    assert draft["draft_state"] == "DRAFT_SUBSET_ONLY"
    assert audit["stage_count"] == 10
    assert audit["evaluator_state"] == "DRAFT_SUBSET_ONLY"
    assert audit["scientific_evaluation_authorized"] is False

    stage_12_common_ser = draft["stages"]["12"]["rules"][1]
    assert [
        context["value"] for context in stage_12_common_ser["contexts"]
    ] == [0.01, 0.03]

    stage_14_bias = draft["stages"]["14"]["rules"][0]
    assert [context["value"] for context in stage_14_bias["contexts"]] == [
        0.10,
        0.15,
        0.25,
    ]

    stage_16c_sensitivity = draft["stages"]["16C"]["rules"][0]
    assert [
        context["value"] for context in stage_16c_sensitivity["contexts"]
    ] == [0.70, 0.80, 0.90, 0.95, 0.95]


def test_scientific_rule_draft_rejects_incomplete_context() -> None:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    draft = load_yaml(SCIENTIFIC_RULE_DRAFT)
    broken = copy.deepcopy(draft)
    del broken["stages"]["12"]["rules"][1]["contexts"][0]["when"]

    with pytest.raises(R13ScientificEvaluatorError, match="condition_missing"):
        validate_metric_rule_subset(
            broken,
            criteria,
            SCIENTIFIC_RULE_DRAFT_STAGES,
        )


def test_scientific_rule_draft_rejects_unknown_decision_value() -> None:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    draft = load_yaml(SCIENTIFIC_RULE_DRAFT)
    broken = copy.deepcopy(draft)
    broken["stages"]["16C"]["rules"][0]["contexts"][0][
        "decision_value"
    ] = "rounded_display"

    with pytest.raises(R13ScientificEvaluatorError, match="decision_value_invalid"):
        validate_metric_rule_subset(
            broken,
            criteria,
            SCIENTIFIC_RULE_DRAFT_STAGES,
        )


def test_integrity_to_reporting_rule_draft_matches_r08_exactly() -> None:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    draft = load_yaml(SCIENTIFIC_COMPLEMENT_DRAFT)

    audit = validate_metric_rule_subset(
        draft,
        criteria,
        SCIENTIFIC_COMPLEMENT_DRAFT_STAGES,
    )

    assert draft["draft_state"] == "DRAFT_SUBSET_ONLY"
    assert audit["stage_count"] == 15
    assert audit["evaluator_state"] == "DRAFT_SUBSET_ONLY"
    assert audit["scientific_evaluation_authorized"] is False

    stage_05 = {rule["metric_id"]: rule for rule in draft["stages"]["05"]["rules"]}
    assert stage_05["high_defect_sensitivity"]["decision_value"] == "lcb95"
    assert stage_05["missingness_error"]["decision_value"] == "ucb95"

    density_contexts = draft["stages"]["06"]["rules"][3]["contexts"]
    assert [context["value"] for context in density_contexts] == [5000, 3000, 1000]
    assert all(context["deterministic"] for context in density_contexts)

    stage_11 = {rule["metric_id"]: rule for rule in draft["stages"]["11"]["rules"]}
    assert stage_11["cm_error"]["value"] == 0.000001
    assert stage_11["cm_error"]["decision_value"] == "estimate"
    assert stage_11["cm_error"]["deterministic"] is True


def test_statistical_threshold_cannot_use_a_point_estimate() -> None:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    draft = load_yaml(SCIENTIFIC_COMPLEMENT_DRAFT)
    broken = copy.deepcopy(draft)
    broken["stages"]["05"]["rules"][0]["decision_value"] = "estimate"

    with pytest.raises(R13ScientificEvaluatorError, match="decision_value_invalid"):
        validate_metric_rule_subset(
            broken,
            criteria,
            SCIENTIFIC_COMPLEMENT_DRAFT_STAGES,
        )


def test_two_reviewed_drafts_cover_the_full_r08_stage_order() -> None:
    criteria = load_yaml(
        REPOSITORY_ROOT / "calibration/configurations/stage_criteria_v1.yaml"
    )
    complement = load_yaml(SCIENTIFIC_COMPLEMENT_DRAFT)
    phasing_to_sensitivity = load_yaml(SCIENTIFIC_RULE_DRAFT)
    fragments = (complement, phasing_to_sensitivity)
    merged_stages = {
        stage_id: next(
            fragment["stages"][stage_id]
            for fragment in fragments
            if stage_id in fragment["stages"]
        )
        for stage_id in STAGE_IDS
    }
    bundle = {
        "schema_version": "1.0.0",
        "release_id": "panel-master-generic-v1-r12c-scientific-validation",
        "scope": "R13_SCIENTIFIC_BLIND_VALIDATION",
        "criteria_version": "1.0.0",
        "stages": merged_stages,
    }

    audit = validate_metric_rule_bundle(bundle, criteria)

    assert tuple(merged_stages) == STAGE_IDS
    assert audit["stage_count"] == 25
    assert audit["scientific_evaluation_authorized"] is False


def test_frozen_metric_rule_manifest_covers_all_r08_stages() -> None:
    audit = validate_metric_rule_manifest(REPOSITORY_ROOT)

    assert audit["bundle_state"] == "FULL_RULE_SET_FROZEN"
    assert audit["fragment_count"] == 2
    assert audit["stage_count"] == 25
    assert audit["rule_count"] > 100
    assert audit["scientific_evaluation_authorized"] is False


def test_metric_decision_engine_distinguishes_point_failure_and_precision() -> None:
    rule = {
        "metric_id": "sensitivity",
        "operator": "min",
        "value": 0.90,
        "decision_value": "lcb95",
        "minimum_independent_n": 500,
    }

    passed = evaluate_metric_observation(
        rule,
        {
            "metric_id": "sensitivity",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.94,
            "lcb95": 0.91,
        },
    )
    imprecise = evaluate_metric_observation(
        rule,
        {
            "metric_id": "sensitivity",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.94,
            "lcb95": 0.89,
        },
    )
    failed = evaluate_metric_observation(
        rule,
        {
            "metric_id": "sensitivity",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.89,
            "lcb95": 0.86,
        },
    )

    assert passed["status"] == "PASS"
    assert imprecise["status"] == "NOT_EVALUABLE_PRECISION"
    assert failed["status"] == "FAIL_POINT_ESTIMATE"


def test_metric_decision_engine_applies_contextual_or_condition() -> None:
    draft = load_yaml(SCIENTIFIC_RULE_DRAFT)
    rule = draft["stages"]["12"]["rules"][1]

    decision = evaluate_metric_observation(
        rule,
        {
            "metric_id": "common_ser",
            "context_id": "common_ser_low_or_no_pedigree",
            "context": {"density": "D_MAX", "pedigree": "NO_PEDIGREE"},
            "independent_n": 100,
            "integrity_passed": True,
            "estimate": 0.02,
            "ucb95": 0.029,
        },
    )

    assert decision["status"] == "PASS"
    assert decision["threshold"] == 0.03

    with pytest.raises(R13ScientificEvaluatorError, match="context_mismatch"):
        evaluate_metric_observation(
            rule,
            {
                "metric_id": "common_ser",
                "context_id": "common_ser_low_or_no_pedigree",
                "context": {"density": "D_MAX", "pedigree": "COMPLETE"},
                "independent_n": 100,
                "integrity_passed": True,
                "estimate": 0.02,
                "ucb95": 0.029,
            },
        )


def test_metric_decision_engine_applies_equivalence_interval() -> None:
    rule = {
        "metric_id": "equivalence",
        "operator": "equivalence",
        "value": 0.20,
        "decision_value": "ci95",
        "minimum_independent_n": 500,
    }

    passed = evaluate_metric_observation(
        rule,
        {
            "metric_id": "equivalence",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.02,
            "ci95": [-0.12, 0.18],
        },
    )
    imprecise = evaluate_metric_observation(
        rule,
        {
            "metric_id": "equivalence",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.02,
            "ci95": [-0.12, 0.22],
        },
    )

    assert passed["status"] == "PASS"
    assert imprecise["status"] == "NOT_EVALUABLE_PRECISION"


def test_metric_decision_engine_enforces_independent_minimum() -> None:
    decision = evaluate_metric_observation(
        {
            "metric_id": "specificity",
            "operator": "min",
            "value": 0.99,
            "decision_value": "lcb95",
            "minimum_independent_n": 3000,
        },
        {
            "metric_id": "specificity",
            "independent_n": 2999,
            "integrity_passed": True,
            "estimate": 1.0,
            "lcb95": 0.999,
        },
    )

    assert decision["status"] == "NOT_EVALUABLE_STRUCTURAL"
    assert decision["reason"] == "INSUFFICIENT_INDEPENDENT_N"


def test_stage_decision_engine_requires_all_coprimary_metrics() -> None:
    rules = {
        "rules": [
            {
                "metric_id": "sensitivity",
                "operator": "min",
                "value": 0.90,
                "decision_value": "lcb95",
                "minimum_independent_n": 500,
            },
            {
                "metric_id": "error",
                "operator": "max",
                "value": 0.05,
                "decision_value": "ucb95",
                "minimum_independent_n": 500,
            },
        ]
    }
    observations = [
        {
            "metric_id": "sensitivity",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.95,
            "lcb95": 0.92,
        },
        {
            "metric_id": "error",
            "independent_n": 500,
            "integrity_passed": True,
            "estimate": 0.04,
            "ucb95": 0.06,
        },
    ]

    decision = evaluate_stage_observations("05", rules, observations)

    assert decision["status"] == "NOT_EVALUABLE_PRECISION"
    assert decision["metric_decision_count"] == 2

    with pytest.raises(R13ScientificEvaluatorError, match="coverage_mismatch"):
        evaluate_stage_observations("05", rules, observations[:1])


def test_unknown_structural_reason_requires_review() -> None:
    decision = evaluate_metric_observation(
        {
            "metric_id": "evaluability",
            "operator": "min",
            "value": 0.95,
            "decision_value": "lcb95",
            "minimum_independent_n": 500,
        },
        {
            "metric_id": "evaluability",
            "independent_n": 500,
            "integrity_passed": True,
            "structural_not_evaluable_reason": "UNPLANNED_DENOMINATOR_LOSS",
        },
    )

    assert decision["status"] == "REVIEW_REQUIRED"


def test_closed_metric_document_applies_all_frozen_rules_without_qualification() -> None:
    decision = evaluate_closed_metric_document(
        REPOSITORY_ROOT,
        _passing_closed_metric_document(),
    )

    assert decision["status"] == "PASS"
    assert decision["stage_count"] == 25
    assert decision["metric_decision_count"] > 135
    assert decision["decision_engine_state"] == "FULL_EVALUATION_READY"
    assert decision["scientific_qualification_authorized"] is False
    assert decision["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }


def test_closed_metric_document_propagates_coprimary_failure() -> None:
    document = _passing_closed_metric_document()
    first_observation = document["stages"]["00"]["observations"][0]  # type: ignore[index]
    first_observation["estimate"] = 0.0

    decision = evaluate_closed_metric_document(REPOSITORY_ROOT, document)

    assert decision["status"] == "FAIL_POINT_ESTIMATE"
    assert decision["stage_decisions"][0]["status"] == "FAIL_POINT_ESTIMATE"


def test_closed_metric_document_refuses_evaluation_before_run_closure() -> None:
    document = _passing_closed_metric_document()
    document["run_state"] = "RUNNING"

    with pytest.raises(R13ScientificEvaluatorError, match="run_not_closed"):
        evaluate_closed_metric_document(REPOSITORY_ROOT, document)
