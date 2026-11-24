"""Contrat de l'évaluateur scientifique post-clôture R13.

L'évaluateur scientifique ne peut être déclaré prêt qu'après sérialisation de
toutes les règles R08 composites et validation des producteurs de métriques.
Ce module contrôle cette couverture sans inventer de seuil ni de résultat.
"""

from __future__ import annotations

from collections.abc import Mapping
import math
from pathlib import Path, PurePosixPath
from typing import Any

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    canonical_sha256,
    file_sha256,
    load_yaml,
)
from effet_fondateur.calibration.evaluation import validate_criteria_document
from effet_fondateur.calibration.release_guard import scan_blinded_document


RELEASE_ID = "panel-master-generic-v1-r12c-scientific-validation"
SCOPE = "R13_SCIENTIFIC_BLIND_VALIDATION"
POST_CLOSURE_SCOPE = "R13_SCIENTIFIC_POST_CLOSURE_METRICS"
SCIENTIFIC_EVALUATOR_STATE = "FULL_EVALUATION_READY"
METRIC_RULE_MANIFEST_PATH = Path(
    "calibration/configurations/scientific_metric_rules_v1.yaml"
)
ALLOWED_OPERATORS = {"min", "max", "exact", "equivalence", "contextual"}
DECISION_STATUS_ORDER = (
    "FAIL_INTEGRITY",
    "NOT_EVALUABLE_STRUCTURAL",
    "FAIL_POINT_ESTIMATE",
    "NOT_EVALUABLE_PRECISION",
    "REVIEW_REQUIRED",
    "PASS",
)
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
DECISION_VALUES_BY_OPERATOR = {
    "min": {"lcb95"},
    "max": {"ucb95"},
    "exact": {"estimate"},
    "equivalence": {"ci95"},
}


class R13ScientificEvaluatorError(RuntimeError):
    """Règles scientifiques incomplètes ou incohérentes avec R08."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13ScientificEvaluatorError(code, detail)


def _allowed_decision_values(operator: str, *, deterministic: bool) -> set[str]:
    """Retourne la statistique admise sans inventer d'IC pour un contrôle exact."""
    if deterministic and operator in {"min", "max"}:
        return {"estimate"}
    return DECISION_VALUES_BY_OPERATOR[operator]


def _validate_context_condition(condition: Mapping[str, Any], detail: str) -> None:
    """Valide le petit DSL contextuel sans interpréter de résultat."""
    if not condition:
        _fail("r13_context_condition_missing", detail)
    logical_keys = set(condition).intersection({"any_of", "all_of"})
    if logical_keys:
        if len(condition) != 1 or len(logical_keys) != 1:
            _fail("r13_context_condition_logic_invalid", detail)
        branches = condition[next(iter(logical_keys))]
        if not isinstance(branches, list) or not branches:
            _fail("r13_context_condition_logic_invalid", detail)
        for index, branch in enumerate(branches):
            if not isinstance(branch, Mapping):
                _fail("r13_context_condition_logic_invalid", f"{detail}[{index}]")
            _validate_context_condition(branch, f"{detail}[{index}]")
        return
    for key, expected in condition.items():
        if not isinstance(key, str) or not key:
            _fail("r13_context_condition_key_invalid", detail)
        if key.endswith(("_min", "_max")):
            if not isinstance(expected, (int, float)) or isinstance(expected, bool):
                _fail("r13_context_condition_bound_invalid", f"{detail}.{key}")
        elif isinstance(expected, list):
            if not expected or any(
                not isinstance(value, (str, int, float)) or isinstance(value, bool)
                for value in expected
            ):
                _fail("r13_context_condition_value_invalid", f"{detail}.{key}")
        elif not isinstance(expected, (str, int, float)) or isinstance(expected, bool):
            _fail("r13_context_condition_value_invalid", f"{detail}.{key}")


def _condition_matches(condition: Mapping[str, Any], observed: Mapping[str, Any]) -> bool:
    if "any_of" in condition:
        return any(_condition_matches(branch, observed) for branch in condition["any_of"])
    if "all_of" in condition:
        return all(_condition_matches(branch, observed) for branch in condition["all_of"])
    for key, expected in condition.items():
        if key.endswith("_min"):
            actual = observed.get(key[:-4])
            if not isinstance(actual, (int, float)) or actual < expected:
                return False
        elif key.endswith("_max"):
            actual = observed.get(key[:-4])
            if not isinstance(actual, (int, float)) or actual > expected:
                return False
        elif isinstance(expected, list):
            if observed.get(key) not in expected:
                return False
        elif observed.get(key) != expected:
            return False
    return True


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    numeric = float(value)
    return numeric if math.isfinite(numeric) else None


def _validate_metric_rule_stages(
    document: Mapping[str, Any],
    criteria: Mapping[str, Any],
    expected_stages: tuple[str, ...],
) -> dict[str, Any]:
    """Valide une couverture ordonnée de métriques primaires R08."""
    scan_blinded_document(document)
    validate_criteria_document(criteria)
    if document.get("schema_version") != "1.0.0":
        _fail("r13_metric_rules_schema_mismatch")
    if document.get("release_id") != RELEASE_ID or document.get("scope") != SCOPE:
        _fail("r13_metric_rules_identity_mismatch")
    if document.get("criteria_version") != criteria.get("criteria_version"):
        _fail("r13_metric_rules_criteria_version_mismatch")
    stages = document.get("stages")
    if not isinstance(stages, Mapping) or tuple(stages) != expected_stages:
        _fail("r13_metric_rules_stage_set_mismatch")

    rule_count = 0
    contextual_count = 0
    for stage_id in expected_stages:
        stage = stages[stage_id]
        if not isinstance(stage, Mapping):
            _fail("r13_metric_stage_invalid", stage_id)
        rules = stage.get("rules")
        if not isinstance(rules, list) or not rules:
            _fail("r13_metric_rules_missing", stage_id)
        expected = list(criteria["stages"][stage_id]["primary_metrics"])
        observed = [rule.get("metric_id") for rule in rules if isinstance(rule, Mapping)]
        if observed != expected:
            _fail("r13_primary_metric_coverage_mismatch", stage_id)
        for rule in rules:
            operator = rule.get("operator")
            if operator not in ALLOWED_OPERATORS:
                _fail("r13_metric_operator_invalid", stage_id)
            if not isinstance(rule.get("unit"), str) or not rule["unit"]:
                _fail("r13_metric_unit_missing", stage_id)
            if not isinstance(rule.get("independent_unit"), str) or not rule[
                "independent_unit"
            ]:
                _fail("r13_metric_independent_unit_missing", stage_id)
            minimum_n = rule.get("minimum_independent_n")
            if not isinstance(minimum_n, int) or minimum_n < 1:
                _fail("r13_metric_minimum_n_invalid", stage_id)
            if operator == "contextual":
                contexts = rule.get("contexts")
                if not isinstance(contexts, list) or not contexts:
                    _fail("r13_contextual_rule_missing", stage_id)
                context_ids: set[str] = set()
                for context in contexts:
                    if not isinstance(context, Mapping):
                        _fail("r13_context_invalid", stage_id)
                    context_id = context.get("context_id")
                    if (
                        not isinstance(context_id, str)
                        or not context_id
                        or context_id in context_ids
                    ):
                        _fail("r13_context_id_invalid", stage_id)
                    context_ids.add(context_id)
                    if not isinstance(context.get("when"), Mapping) or not context["when"]:
                        _fail("r13_context_condition_missing", context_id)
                    _validate_context_condition(context["when"], context_id)
                    if context.get("operator") not in ALLOWED_OPERATORS - {"contextual"}:
                        _fail("r13_context_operator_invalid", context_id)
                    if not isinstance(context.get("value"), (int, float)):
                        _fail("r13_context_threshold_missing", context_id)
                    context_operator = context["operator"]
                    deterministic = context.get("deterministic", False)
                    if not isinstance(deterministic, bool):
                        _fail("r13_context_deterministic_flag_invalid", context_id)
                    if context.get("decision_value") not in _allowed_decision_values(
                        context_operator,
                        deterministic=deterministic,
                    ):
                        _fail("r13_context_decision_value_invalid", context_id)
                contextual_count += 1
            else:
                if not isinstance(rule.get("value"), (int, float)):
                    _fail("r13_metric_threshold_missing", stage_id)
                deterministic = rule.get("deterministic", False)
                if not isinstance(deterministic, bool):
                    _fail("r13_metric_deterministic_flag_invalid", stage_id)
                if rule.get("decision_value") not in _allowed_decision_values(
                    operator,
                    deterministic=deterministic,
                ):
                    _fail("r13_metric_decision_value_invalid", stage_id)
            rule_count += 1
    return {
        "release_id": RELEASE_ID,
        "stage_count": len(expected_stages),
        "rule_count": rule_count,
        "contextual_rule_count": contextual_count,
        "evaluator_state": SCIENTIFIC_EVALUATOR_STATE,
        "scientific_evaluation_authorized": False,
    }


def validate_metric_rule_bundle(
    document: Mapping[str, Any],
    criteria: Mapping[str, Any],
) -> dict[str, Any]:
    """Exige une règle exécutable pour chaque métrique primaire des 25 étapes."""
    return _validate_metric_rule_stages(document, criteria, STAGE_IDS)


def validate_metric_rule_subset(
    document: Mapping[str, Any],
    criteria: Mapping[str, Any],
    expected_stages: tuple[str, ...],
) -> dict[str, Any]:
    """Valide un brouillon borné sans le rendre admissible à G5."""
    if not expected_stages or any(stage not in STAGE_IDS for stage in expected_stages):
        _fail("r13_metric_subset_invalid")
    audit = _validate_metric_rule_stages(document, criteria, expected_stages)
    return {
        **audit,
        "evaluator_state": "DRAFT_SUBSET_ONLY",
        "scientific_evaluation_authorized": False,
    }


def evaluate_metric_observation(
    rule: Mapping[str, Any],
    observation: Mapping[str, Any],
) -> dict[str, Any]:
    """Applique une règle R08 à une observation métrique déjà calculée.

    Cette fonction ne calcule ni intervalle ni vérité biologique. Elle exige
    les valeurs non arrondies produites après clôture et applique uniquement
    l'ordre de décision préspecifié.
    """
    metric_id = rule.get("metric_id")
    if not isinstance(metric_id, str) or observation.get("metric_id") != metric_id:
        _fail("r13_metric_observation_identity_mismatch", str(metric_id))
    operator = rule.get("operator")
    selected_rule: Mapping[str, Any] = rule
    context_id: str | None = None
    if operator == "contextual":
        context_id = observation.get("context_id")
        observed_context = observation.get("context")
        if not isinstance(context_id, str) or not isinstance(observed_context, Mapping):
            _fail("r13_metric_observation_context_missing", metric_id)
        contexts = rule.get("contexts")
        if not isinstance(contexts, list):
            _fail("r13_contextual_rule_missing", metric_id)
        matching_ids = [
            context
            for context in contexts
            if isinstance(context, Mapping) and context.get("context_id") == context_id
        ]
        if len(matching_ids) != 1:
            _fail("r13_metric_observation_context_unknown", context_id)
        selected_rule = matching_ids[0]
        condition = selected_rule.get("when")
        if not isinstance(condition, Mapping) or not _condition_matches(
            condition, observed_context
        ):
            _fail("r13_metric_observation_context_mismatch", context_id)
        operator = selected_rule.get("operator")

    independent_n = observation.get("independent_n")
    if not isinstance(independent_n, int) or isinstance(independent_n, bool) or independent_n < 0:
        _fail("r13_metric_observation_independent_n_invalid", metric_id)
    minimum_n = rule.get("minimum_independent_n")
    if not isinstance(minimum_n, int):
        _fail("r13_metric_minimum_n_invalid", metric_id)
    result = {
        "metric_id": metric_id,
        "context_id": context_id,
        "operator": operator,
        "decision_value": selected_rule.get("decision_value"),
        "threshold": selected_rule.get("value"),
        "independent_n": independent_n,
        "minimum_independent_n": minimum_n,
    }
    if observation.get("integrity_passed") is not True:
        return {**result, "status": "FAIL_INTEGRITY"}
    structural_reason = observation.get("structural_not_evaluable_reason")
    if structural_reason is not None:
        if not isinstance(structural_reason, str) or not structural_reason:
            _fail("r13_metric_observation_structural_reason_invalid", metric_id)
        if structural_reason not in STRUCTURAL_NOT_EVALUABLE_REASONS:
            return {
                **result,
                "status": "REVIEW_REQUIRED",
                "reason": structural_reason,
            }
        return {
            **result,
            "status": "NOT_EVALUABLE_STRUCTURAL",
            "reason": structural_reason,
        }
    if independent_n < minimum_n:
        return {
            **result,
            "status": "NOT_EVALUABLE_STRUCTURAL",
            "reason": "INSUFFICIENT_INDEPENDENT_N",
        }

    threshold = _finite_number(selected_rule.get("value"))
    estimate = _finite_number(observation.get("estimate"))
    if threshold is None or estimate is None:
        return {**result, "status": "FAIL_INTEGRITY", "reason": "NON_FINITE_VALUE"}
    if operator == "exact":
        status = "PASS" if estimate == threshold else "FAIL_POINT_ESTIMATE"
        return {**result, "estimate": estimate, "status": status}

    point_passes = (
        estimate >= threshold
        if operator == "min"
        else estimate <= threshold
        if operator == "max"
        else abs(estimate) <= threshold
        if operator == "equivalence"
        else False
    )
    if not point_passes:
        return {**result, "estimate": estimate, "status": "FAIL_POINT_ESTIMATE"}

    decision_value = selected_rule.get("decision_value")
    if decision_value == "estimate":
        return {**result, "estimate": estimate, "status": "PASS"}
    if decision_value in {"lcb95", "ucb95"}:
        bound = _finite_number(observation.get(decision_value))
        if bound is None:
            return {
                **result,
                "estimate": estimate,
                "status": "FAIL_INTEGRITY",
                "reason": "DECISION_BOUND_MISSING",
            }
        precision_passes = bound >= threshold if operator == "min" else bound <= threshold
        return {
            **result,
            "estimate": estimate,
            decision_value: bound,
            "status": "PASS" if precision_passes else "NOT_EVALUABLE_PRECISION",
        }
    if decision_value == "ci95":
        interval = observation.get("ci95")
        if (
            not isinstance(interval, list)
            or len(interval) != 2
            or _finite_number(interval[0]) is None
            or _finite_number(interval[1]) is None
        ):
            return {
                **result,
                "estimate": estimate,
                "status": "FAIL_INTEGRITY",
                "reason": "DECISION_INTERVAL_MISSING",
            }
        lower, upper = float(interval[0]), float(interval[1])
        if lower > upper:
            return {
                **result,
                "estimate": estimate,
                "status": "FAIL_INTEGRITY",
                "reason": "DECISION_INTERVAL_INVALID",
            }
        precision_passes = lower >= -threshold and upper <= threshold
        return {
            **result,
            "estimate": estimate,
            "ci95": [lower, upper],
            "status": "PASS" if precision_passes else "NOT_EVALUABLE_PRECISION",
        }
    return {
        **result,
        "estimate": estimate,
        "status": "FAIL_INTEGRITY",
        "reason": "DECISION_VALUE_UNSUPPORTED",
    }


def evaluate_stage_observations(
    stage_id: str,
    stage_rules: Mapping[str, Any],
    observations: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Évalue toutes les métriques co-primaires d'une étape R08."""
    if stage_id not in STAGE_IDS:
        _fail("r13_stage_observation_stage_invalid", stage_id)
    rules = stage_rules.get("rules")
    if not isinstance(rules, list) or not rules:
        _fail("r13_metric_rules_missing", stage_id)
    indexed: dict[tuple[str, str | None], Mapping[str, Any]] = {}
    for observation in observations:
        if not isinstance(observation, Mapping):
            _fail("r13_metric_observation_invalid", stage_id)
        metric_id = observation.get("metric_id")
        context_id = observation.get("context_id")
        if not isinstance(metric_id, str) or (
            context_id is not None and not isinstance(context_id, str)
        ):
            _fail("r13_metric_observation_identity_invalid", stage_id)
        key = (metric_id, context_id)
        if key in indexed:
            _fail("r13_metric_observation_duplicate", f"{stage_id}:{key}")
        indexed[key] = observation

    expected: list[tuple[Mapping[str, Any], tuple[str, str | None]]] = []
    for rule in rules:
        if not isinstance(rule, Mapping) or not isinstance(rule.get("metric_id"), str):
            _fail("r13_metric_rule_invalid", stage_id)
        metric_id = rule["metric_id"]
        if rule.get("operator") == "contextual":
            contexts = rule.get("contexts")
            if not isinstance(contexts, list) or not contexts:
                _fail("r13_contextual_rule_missing", metric_id)
            for context in contexts:
                if not isinstance(context, Mapping) or not isinstance(
                    context.get("context_id"), str
                ):
                    _fail("r13_context_id_invalid", metric_id)
                expected.append((rule, (metric_id, context["context_id"])))
        else:
            expected.append((rule, (metric_id, None)))
    expected_keys = {key for _, key in expected}
    if set(indexed) != expected_keys:
        missing = sorted(expected_keys - set(indexed), key=repr)
        unexpected = sorted(set(indexed) - expected_keys, key=repr)
        _fail(
            "r13_stage_observation_coverage_mismatch",
            f"{stage_id}:missing={missing}:unexpected={unexpected}",
        )

    decisions = [
        evaluate_metric_observation(rule, indexed[key]) for rule, key in expected
    ]
    stage_status = min(
        (decision["status"] for decision in decisions),
        key=DECISION_STATUS_ORDER.index,
    )
    return {
        "stage_id": stage_id,
        "status": stage_status,
        "metric_decision_count": len(decisions),
        "metric_decisions": decisions,
    }


def _load_metric_rule_manifest(
    repository_root: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Charge les fragments empreintés et leur couverture intégrale de R08."""
    repository_root = repository_root.resolve(strict=True)
    manifest_path = repository_root / METRIC_RULE_MANIFEST_PATH
    manifest = load_yaml(manifest_path)
    scan_blinded_document(manifest)
    if manifest.get("schema_version") != "1.0.0":
        _fail("r13_metric_manifest_schema_mismatch")
    if manifest.get("release_id") != RELEASE_ID or manifest.get("scope") != SCOPE:
        _fail("r13_metric_manifest_identity_mismatch")
    if manifest.get("bundle_state") != "FULL_RULE_SET_FROZEN":
        _fail("r13_metric_manifest_state_mismatch")

    criteria = load_yaml(
        repository_root / "calibration/configurations/stage_criteria_v1.yaml"
    )
    if manifest.get("criteria_version") != criteria.get("criteria_version"):
        _fail("r13_metric_manifest_criteria_version_mismatch")
    fragments = manifest.get("fragments")
    if not isinstance(fragments, list) or not fragments:
        _fail("r13_metric_manifest_fragments_missing")

    merged: dict[str, Any] = {}
    fragment_audits = []
    fragment_ids: set[str] = set()
    for fragment_entry in fragments:
        if not isinstance(fragment_entry, Mapping):
            _fail("r13_metric_manifest_fragment_invalid")
        fragment_id = fragment_entry.get("fragment_id")
        if (
            not isinstance(fragment_id, str)
            or not fragment_id
            or fragment_id in fragment_ids
        ):
            _fail("r13_metric_manifest_fragment_id_invalid")
        fragment_ids.add(fragment_id)
        relative_path = fragment_entry.get("relative_path")
        if not isinstance(relative_path, str) or not relative_path:
            _fail("r13_metric_manifest_fragment_path_invalid", fragment_id)
        relative = PurePosixPath(relative_path)
        if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
            _fail("r13_metric_manifest_fragment_path_invalid", fragment_id)
        fragment_path = repository_root.joinpath(*relative.parts)
        if fragment_path.is_symlink() or not fragment_path.is_file():
            _fail("r13_metric_manifest_fragment_missing", fragment_id)
        if file_sha256(fragment_path) != fragment_entry.get("sha256"):
            _fail("r13_metric_manifest_fragment_hash_mismatch", fragment_id)
        stages = fragment_entry.get("stages")
        if not isinstance(stages, list) or not stages or any(
            not isinstance(stage_id, str) for stage_id in stages
        ):
            _fail("r13_metric_manifest_fragment_stages_invalid", fragment_id)
        expected_stages = tuple(stages)
        fragment = load_yaml(fragment_path)
        fragment_audits.append(
            validate_metric_rule_subset(fragment, criteria, expected_stages)
        )
        for stage_id in expected_stages:
            if stage_id in merged:
                _fail("r13_metric_manifest_duplicate_stage", stage_id)
            merged[stage_id] = fragment["stages"][stage_id]

    if set(merged) != set(STAGE_IDS):
        _fail("r13_metric_manifest_stage_coverage_mismatch")
    ordered_bundle = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "scope": SCOPE,
        "criteria_version": criteria["criteria_version"],
        "stages": {stage_id: merged[stage_id] for stage_id in STAGE_IDS},
    }
    audit = validate_metric_rule_bundle(ordered_bundle, criteria)
    manifest_audit = {
        **audit,
        "bundle_state": "FULL_RULE_SET_FROZEN",
        "manifest_sha256": file_sha256(manifest_path),
        "expanded_bundle_sha256": canonical_sha256(ordered_bundle),
        "fragment_count": len(fragment_audits),
        "scientific_evaluation_authorized": False,
    }
    return manifest_audit, ordered_bundle


def validate_metric_rule_manifest(repository_root: Path) -> dict[str, Any]:
    """Valide le manifeste figé sans exposer son document développé."""
    audit, _ = _load_metric_rule_manifest(repository_root)
    return audit


def load_metric_rule_bundle(repository_root: Path) -> dict[str, Any]:
    """Retourne les 25 règles développées après validation de toutes les empreintes."""
    _, bundle = _load_metric_rule_manifest(repository_root)
    return bundle


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def evaluate_closed_metric_document(
    repository_root: Path,
    document: Mapping[str, Any],
) -> dict[str, Any]:
    """Applique les 135 règles à un document de métriques post-clôture.

    Le calcul scientifique des métriques et de leurs intervalles reste la
    responsabilité de producteurs figés distincts. Ce moteur vérifie leur
    couverture et applique seulement les décisions R08 non arrondies.
    """
    if document.get("schema_version") != "1.0.0":
        _fail("r13_metric_document_schema_mismatch")
    if document.get("release_id") != RELEASE_ID:
        _fail("r13_metric_document_release_mismatch")
    if document.get("scope") != POST_CLOSURE_SCOPE:
        _fail("r13_metric_document_scope_mismatch")
    if document.get("run_state") != "CLOSED":
        _fail("r13_metric_document_run_not_closed")
    if document.get("truth_access_state") != "AUTHORIZED_POST_CLOSURE":
        _fail("r13_metric_document_truth_access_invalid")
    if document.get("metric_producers_frozen") is not True:
        _fail("r13_metric_document_producers_not_frozen")
    for field in ("closed_run_sha256", "metric_producer_bundle_sha256"):
        if not _is_sha256(document.get(field)):
            _fail("r13_metric_document_evidence_hash_invalid", field)
    stages = document.get("stages")
    if not isinstance(stages, Mapping) or tuple(stages) != STAGE_IDS:
        _fail("r13_metric_document_stage_set_mismatch")

    bundle = load_metric_rule_bundle(repository_root)
    stage_decisions = []
    for stage_id in STAGE_IDS:
        stage = stages[stage_id]
        if not isinstance(stage, Mapping):
            _fail("r13_metric_document_stage_invalid", stage_id)
        observations = stage.get("observations")
        if not isinstance(observations, list):
            _fail("r13_metric_document_observations_invalid", stage_id)
        stage_decisions.append(
            evaluate_stage_observations(
                stage_id,
                bundle["stages"][stage_id],
                observations,
            )
        )
    overall_status = min(
        (decision["status"] for decision in stage_decisions),
        key=DECISION_STATUS_ORDER.index,
    )
    return {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "scope": POST_CLOSURE_SCOPE,
        "status": overall_status,
        "stage_count": len(stage_decisions),
        "metric_decision_count": sum(
            decision["metric_decision_count"] for decision in stage_decisions
        ),
        "stage_decisions": stage_decisions,
        "source_closed_run_sha256": document["closed_run_sha256"],
        "source_metric_producer_bundle_sha256": document[
            "metric_producer_bundle_sha256"
        ],
        "decision_engine_state": SCIENTIFIC_EVALUATOR_STATE,
        "scientific_qualification_authorized": False,
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
    }


def assert_scientific_evaluator_ready() -> None:
    """Refuse toute évaluation tant que les producteurs et règles ne sont pas gelés."""
    if SCIENTIFIC_EVALUATOR_STATE != "FULL_EVALUATION_READY":
        _fail("r13_scientific_evaluator_not_ready")
