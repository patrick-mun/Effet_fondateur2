"""Moteur générique de décision préspecifié à R08.

Le moteur applique l'ordre de décision R08.1. Il ne calcule pas les métriques
scientifiques : les évaluateurs d'étape lui remettent des estimations et leurs
bornes calculées sur les unités indépendantes prévues.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from effet_fondateur.calibration.contracts import STAGE_IDS


class EvaluationStatus(StrEnum):
    PASS = "PASS"
    FAIL = "FAIL"
    NOT_EVALUABLE = "NOT_EVALUABLE"
    REVIEW_REQUIRED = "REVIEW_REQUIRED"


STRUCTURAL_NOT_EVALUABLE_REASONS = frozenset(
    {
        "INFEASIBLE_TARGET_CATALOG",
        "INFEASIBLE_WINDOW",
        "INFEASIBLE_ROH",
        "INSUFFICIENT_INDEPENDENT_N",
        "INSUFFICIENT_DENSITY",
        "INSUFFICIENT_MINOR_COPIES",
        "INSUFFICIENT_INFORMATIVE_TRANSITIONS",
        "OPTIONAL_CELL_ABSENT_BY_DESIGN",
        "ANCESTRY_PROPORTION_METHOD_ABSENT",
        "CENSORING_MODEL_ABSENT",
        "RESOURCE_BLOCKED",
    }
)


class CalibrationEvaluationError(ValueError):
    """Contrat de métrique incohérent ou incomplet."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


@dataclass(frozen=True)
class MetricDecision:
    metric_id: str
    status: EvaluationStatus
    reason: str


@dataclass(frozen=True)
class StageDecision:
    stage_id: str
    status: EvaluationStatus
    metric_decisions: tuple[MetricDecision, ...]
    reasons: tuple[str, ...]


def _number(value: Any, field: str, metric_id: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CalibrationEvaluationError("metric_value_invalid", f"{metric_id}:{field}")
    numeric = float(value)
    if numeric != numeric or numeric in {float("inf"), float("-inf")}:
        raise CalibrationEvaluationError("metric_value_non_finite", f"{metric_id}:{field}")
    return numeric


def decide_metric(rule: Mapping[str, Any], observation: Mapping[str, Any]) -> MetricDecision:
    """Décide une métrique minimale, maximale, exacte ou d'équivalence."""
    metric_id = rule.get("metric_id")
    operator = rule.get("operator")
    if not isinstance(metric_id, str) or operator not in {"min", "max", "exact", "equivalence"}:
        raise CalibrationEvaluationError("metric_rule_invalid", str(metric_id))
    structural_reason = observation.get("structural_reason")
    if structural_reason is not None:
        if structural_reason not in STRUCTURAL_NOT_EVALUABLE_REASONS:
            return MetricDecision(metric_id, EvaluationStatus.REVIEW_REQUIRED, "unknown_missing_reason")
        return MetricDecision(metric_id, EvaluationStatus.NOT_EVALUABLE, str(structural_reason))
    estimate = _number(observation.get("estimate"), "estimate", metric_id)

    if operator == "exact":
        expected = _number(rule.get("value"), "value", metric_id)
        passed = estimate == expected
        return MetricDecision(
            metric_id,
            EvaluationStatus.PASS if passed else EvaluationStatus.FAIL,
            "exact_match" if passed else "point_estimate_failed",
        )

    threshold = _number(rule.get("value"), "value", metric_id)
    if operator == "min":
        if estimate < threshold:
            return MetricDecision(metric_id, EvaluationStatus.FAIL, "point_estimate_failed")
        lower = observation.get("lcb95")
        if lower is None or _number(lower, "lcb95", metric_id) < threshold:
            return MetricDecision(metric_id, EvaluationStatus.NOT_EVALUABLE, "precision_insufficient")
        return MetricDecision(metric_id, EvaluationStatus.PASS, "lower_bound_passed")
    if operator == "max":
        if estimate > threshold:
            return MetricDecision(metric_id, EvaluationStatus.FAIL, "point_estimate_failed")
        upper = observation.get("ucb95")
        if upper is None or _number(upper, "ucb95", metric_id) > threshold:
            return MetricDecision(metric_id, EvaluationStatus.NOT_EVALUABLE, "precision_insufficient")
        return MetricDecision(metric_id, EvaluationStatus.PASS, "upper_bound_passed")

    margin = threshold
    if abs(estimate) > margin:
        return MetricDecision(metric_id, EvaluationStatus.FAIL, "point_estimate_failed")
    lower = observation.get("ci95_low")
    upper = observation.get("ci95_high")
    if lower is None or upper is None:
        return MetricDecision(metric_id, EvaluationStatus.NOT_EVALUABLE, "precision_insufficient")
    if _number(lower, "ci95_low", metric_id) < -margin or _number(
        upper, "ci95_high", metric_id
    ) > margin:
        return MetricDecision(metric_id, EvaluationStatus.NOT_EVALUABLE, "precision_insufficient")
    return MetricDecision(metric_id, EvaluationStatus.PASS, "equivalence_interval_passed")


def decide_stage(
    *,
    stage_id: str,
    rules: Sequence[Mapping[str, Any]],
    observations: Mapping[str, Mapping[str, Any]],
    hard_fail_reasons: Sequence[str] = (),
    manual_review_pending: bool = False,
) -> StageDecision:
    """Applique l'intersection–union et l'ordre de statut R08.1."""
    if stage_id not in STAGE_IDS:
        raise CalibrationEvaluationError("unknown_stage", stage_id)
    if hard_fail_reasons:
        return StageDecision(
            stage_id,
            EvaluationStatus.FAIL,
            (),
            tuple(hard_fail_reasons),
        )
    decisions: list[MetricDecision] = []
    for rule in rules:
        metric_id = rule.get("metric_id")
        if not isinstance(metric_id, str) or metric_id not in observations:
            return StageDecision(
                stage_id,
                EvaluationStatus.FAIL,
                tuple(decisions),
                (f"required_metric_missing:{metric_id}",),
            )
        decisions.append(decide_metric(rule, observations[metric_id]))
    statuses = {decision.status for decision in decisions}
    if EvaluationStatus.FAIL in statuses:
        status = EvaluationStatus.FAIL
    elif EvaluationStatus.NOT_EVALUABLE in statuses:
        status = EvaluationStatus.NOT_EVALUABLE
    elif EvaluationStatus.REVIEW_REQUIRED in statuses or manual_review_pending:
        status = EvaluationStatus.REVIEW_REQUIRED
    else:
        status = EvaluationStatus.PASS
    reasons = tuple(
        f"{decision.metric_id}:{decision.reason}"
        for decision in decisions
        if decision.status is not EvaluationStatus.PASS
    )
    if manual_review_pending and status is EvaluationStatus.REVIEW_REQUIRED:
        reasons += ("manual_review_pending",)
    return StageDecision(stage_id, status, tuple(decisions), reasons)


def validate_criteria_document(document: Mapping[str, Any]) -> None:
    """Vérifie la couverture et le gel du contrat sérialisé des 25 étapes."""
    if document.get("criteria_version") != "1.0.0":
        raise CalibrationEvaluationError("criteria_version_mismatch")
    if document.get("decision_order") != [
        "FAIL_INTEGRITY",
        "NOT_EVALUABLE_STRUCTURAL",
        "FAIL_POINT_ESTIMATE",
        "NOT_EVALUABLE_PRECISION",
        "REVIEW_REQUIRED",
        "PASS",
    ]:
        raise CalibrationEvaluationError("decision_order_mismatch")
    stages = document.get("stages")
    if not isinstance(stages, dict) or tuple(stages) != STAGE_IDS:
        raise CalibrationEvaluationError("criteria_stage_set_mismatch")
    for stage_id, contract in stages.items():
        if not isinstance(contract, Mapping):
            raise CalibrationEvaluationError("stage_contract_invalid", stage_id)
        if not contract.get("plan_section") or not contract.get("primary_metrics"):
            raise CalibrationEvaluationError("stage_contract_incomplete", stage_id)
        if not isinstance(contract.get("pass_contract"), str) or not contract["pass_contract"]:
            raise CalibrationEvaluationError("stage_pass_contract_missing", stage_id)
    reasons = document.get("structural_not_evaluable_reasons")
    if set(reasons or ()) != STRUCTURAL_NOT_EVALUABLE_REASONS:
        raise CalibrationEvaluationError("not_evaluable_reason_set_mismatch")
