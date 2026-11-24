"""Oracle indépendant D24 pour le segment fondateur candidat de l’étape 13."""

from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
from statistics import fmean
from typing import Any, Iterable


METHOD_ID = "founder_haplotype_truth_bench_v1"
PROFILES = ("MAXIMUM", "BASELINE", "DENSITY_50", "MISSING_1PCT")


def wilson_interval(successes: int, total: int, confidence: float = 0.95) -> tuple[float, float] | None:
    """Retourne l’intervalle binomial de Wilson, sans certitude artificielle."""
    if total == 0:
        return None
    if successes < 0 or successes > total or confidence != 0.95:
        raise ValueError("d24_binomial_inputs_invalid")
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def rate(successes: int, total: int) -> dict[str, int | float | None]:
    """Encode numérateur, dénominateur, estimation et IC descriptif à 95 %."""
    interval = wilson_interval(successes, total)
    return {
        "numerator": successes,
        "denominator": total,
        "estimate": None if total == 0 else successes / total,
        "ci95_low": None if interval is None else interval[0],
        "ci95_high": None if interval is None else interval[1],
    }


def _mean(values: list[float]) -> float | None:
    return None if not values else float(fmean(values))


def true_consensus(truth: dict[str, Any]) -> dict[str, float | int] | None:
    """Calcule l’intersection vraie des segments familiaux préspecifiés."""
    segments = truth["segments"]
    if not segments:
        return None
    target = truth["target"]
    result: dict[str, float | int | bool] = {
        "left_bp": max(row["left_boundary_bp"] for row in segments),
        "right_bp": min(row["right_boundary_bp"] for row in segments),
        "left_cm": max(row["left_boundary_cm"] for row in segments),
        "right_cm": min(row["right_boundary_cm"] for row in segments),
    }
    result["left_censored"] = all(
        row["left_censored"]
        for row in segments
        if row["left_boundary_cm"] == result["left_cm"]
    )
    result["right_censored"] = all(
        row["right_censored"]
        for row in segments
        if row["right_boundary_cm"] == result["right_cm"]
    )
    contains = (
        result["left_bp"] <= target["position_bp"] <= result["right_bp"]
        and result["left_cm"] <= target["position_cm"] <= result["right_cm"]
    )
    if truth["detectable"] and not contains:
        raise ValueError("d24_detectable_truth_has_empty_intersection")
    return result


def _predicted_consensus(
    truth: dict[str, Any], consensus: dict[str, Any]
) -> dict[str, float | int]:
    target_cm = float(truth["target"]["position_cm"])
    return {
        "left_bp": int(consensus["LEFT_BOUND_BP"]),
        "right_bp": int(consensus["RIGHT_BOUND_BP"]),
        "left_cm": target_cm - float(consensus["LEFT_LENGTH_CM"]),
        "right_cm": target_cm + float(consensus["RIGHT_LENGTH_CM"]),
    }


def evaluate_case(
    truth: dict[str, Any],
    summary: dict[str, Any],
    consensus: dict[str, Any],
    segments: list[dict[str, Any]],
    *,
    minimum_independent_units: int,
    minimum_flank_markers: int = 2,
) -> dict[str, Any]:
    """Compare un scénario clos à sa vérité sans exposer ses identifiants."""
    if truth["independent_unit_count"] < minimum_independent_units:
        raise ValueError("d24_independent_units_insufficient")
    truth_families = [row["family_id"] for row in truth["segments"]]
    if len(truth_families) != len(set(truth_families)):
        raise ValueError("d24_truth_family_duplicate")
    if truth["truth_role"] == "POSITIVE" and len(truth_families) < minimum_independent_units:
        raise ValueError("d24_positive_truth_segments_insufficient")
    if summary.get("minimum_independent_carriers") != minimum_independent_units:
        raise ValueError("d24_stage13_minimum_carriers_changed")
    if summary.get("minimum_flank_markers") != minimum_flank_markers:
        raise ValueError("d24_stage13_minimum_flanks_changed")
    target = truth["target"]
    if (
        consensus.get("TARGET_VARIANT_ID") != target["variant_id"]
        or int(consensus.get("TARGET_BP", 0)) != target["position_bp"]
        or any(
            row.get("TARGET_VARIANT_ID") != target["variant_id"]
            or int(row.get("TARGET_BP", 0)) != target["position_bp"]
            for row in segments
        )
    ):
        raise ValueError("d24_stage13_target_mismatch")
    expected_detection = truth["truth_role"] == "POSITIVE" and truth["detectable"]
    detected = (
        summary.get("status") == "SUPPORTED_IBS_CANDIDATE"
        and summary.get("interpretation") == "IBS_SHARED_CANDIDATE"
        and consensus.get("CONSENSUS_STATUS") == "SUPPORTED_IBS_CANDIDATE"
        and consensus.get("INTERPRETATION") == "IBS_SHARED_CANDIDATE"
    )
    ibd_claimed = bool(summary.get("ibd_claimed", False))
    target_contained: bool | None = None
    errors = {key: [] for key in (
        "left_bp", "right_bp", "left_cm", "right_cm",
        "family_left_bp", "family_right_bp", "family_left_cm", "family_right_cm",
    )}
    true_interval = true_consensus(truth)
    if detected and true_interval is not None:
        predicted = _predicted_consensus(truth, consensus)
        target_contained = (
            predicted["left_bp"] <= target["position_bp"] <= predicted["right_bp"]
            and predicted["left_cm"] <= target["position_cm"] <= predicted["right_cm"]
        )
        for side in ("left", "right"):
            if true_interval[f"{side}_censored"]:
                continue
            for unit in ("bp", "cm"):
                key = f"{side}_{unit}"
                errors[key].append(float(predicted[key]) - float(true_interval[key]))
        predicted_by_family = {
            row["FAMILY_ID"]: row
            for row in segments
            if row.get("SEGMENT_STATUS") == "INCLUDED"
        }
        if len(predicted_by_family) != sum(
            row.get("SEGMENT_STATUS") == "INCLUDED" for row in segments
        ):
            raise ValueError("d24_predicted_family_duplicate")
        target_cm = float(target["position_cm"])
        for row in truth["segments"]:
            predicted_family = predicted_by_family.get(row["family_id"])
            if predicted_family is None:
                continue
            values = {
                "left_bp": int(predicted_family["LEFT_BOUND_BP"]),
                "right_bp": int(predicted_family["RIGHT_BOUND_BP"]),
                "left_cm": target_cm - float(predicted_family["LEFT_LENGTH_CM"]),
                "right_cm": target_cm + float(predicted_family["RIGHT_LENGTH_CM"]),
            }
            for side in ("left", "right"):
                if row[f"{side}_censored"]:
                    continue
                for unit in ("bp", "cm"):
                    errors[f"family_{side}_{unit}"].append(
                        float(values[f"{side}_{unit}"])
                        - float(row[f"{side}_boundary_{unit}"])
                    )
    missing_family_count = 0
    if truth["truth_role"] == "POSITIVE":
        included = {
            row["FAMILY_ID"]
            for row in segments
            if row.get("SEGMENT_STATUS") == "INCLUDED"
        }
        missing_family_count = len(set(truth_families) - included)
    return {
        "scenario_id": truth["scenario_id"],
        "profile": truth["profile"],
        "truth_role": truth["truth_role"],
        "age_generations": truth["age_generations"],
        "source_manifest_sha256": truth["source_manifest_sha256"],
        "truth_signature": hashlib.sha256(json.dumps(
            {
                key: truth[key]
                for key in (
                    "scenario_id", "truth_role", "detectable", "age_generations",
                    "independent_unit_count", "target", "segments",
                )
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()).hexdigest(),
        "truth_positive": truth["truth_role"] == "POSITIVE",
        "detectable": bool(truth["detectable"]),
        "expected_detection": expected_detection,
        "detected": detected,
        "target_contained": target_contained,
        "ibd_claimed": ibd_claimed,
        "missing_family_count": missing_family_count,
        "errors": errors,
    }


def _boundary_metrics(results: Iterable[dict[str, Any]]) -> dict[str, Any]:
    combined: defaultdict[str, list[float]] = defaultdict(list)
    evaluated = 0
    for result in results:
        if result["errors"]["left_cm"] or result["errors"]["right_cm"]:
            evaluated += 1
        for key, values in result["errors"].items():
            combined[key].extend(values)
    return {
        "evaluated_segment_count": evaluated,
        "evaluated_boundary_count": len(combined["left_cm"]) + len(combined["right_cm"]),
        "left_mean_signed_error_bp": _mean(combined["left_bp"]),
        "right_mean_signed_error_bp": _mean(combined["right_bp"]),
        "left_mean_absolute_error_bp": _mean([abs(value) for value in combined["left_bp"]]),
        "right_mean_absolute_error_bp": _mean([abs(value) for value in combined["right_bp"]]),
        "left_mean_signed_error_cm": _mean(combined["left_cm"]),
        "right_mean_signed_error_cm": _mean(combined["right_cm"]),
        "left_mean_absolute_error_cm": _mean([abs(value) for value in combined["left_cm"]]),
        "right_mean_absolute_error_cm": _mean([abs(value) for value in combined["right_cm"]]),
        "boundary_mean_absolute_error": _mean(
            [abs(value) for value in combined["left_cm"] + combined["right_cm"]]
        ),
        "family_boundary_pair_count": min(
            len(combined["family_left_cm"]), len(combined["family_right_cm"])
        ),
    }


def _profile_aggregate(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Calcule les métriques d’un profil sans mélanger les profils appariés."""
    positives = [row for row in results if row["truth_positive"]]
    detectable = [row for row in positives if row["detectable"]]
    negatives = [row for row in results if not row["truth_positive"]]
    true_positive = sum(row["detected"] for row in detectable)
    false_negative = len(detectable) - true_positive
    false_positive = sum(row["detected"] for row in negatives)
    true_negative = len(negatives) - false_positive
    detected_rows = [row for row in detectable if row["detected"]]
    contained = sum(row["target_contained"] is True for row in detectable)
    ibd_count = sum(row["ibd_claimed"] for row in results)
    boundary = _boundary_metrics(detected_rows)
    metrics = {
        "founder_segment_detection_sensitivity": rate(true_positive, len(detectable)),
        "founder_segment_detection_specificity": rate(true_negative, len(negatives)),
        "target_containment_rate_when_detectable": rate(contained, len(detectable)),
        "boundary_errors": boundary,
        "ibs_reclassified_as_ibd_without_explicit_ibd_count": ibd_count,
    }
    return {
        "scenario_counts": {
            "total": len(results),
            "positive": len(positives),
            "negative": len(negatives),
            "detectable_positive": len(detectable),
        },
        "confusion": {
            "true_positive": true_positive,
            "false_negative": false_negative,
            "true_negative": true_negative,
            "false_positive": false_positive,
        },
        "metrics": metrics,
    }


def aggregate_results(
    results: list[dict[str, Any]], *, primary_profile: str = "MAXIMUM"
) -> dict[str, Any]:
    """Agrège le primaire seul et décrit séparément les profils appariés."""
    primary = [row for row in results if row["profile"] == primary_profile]
    aggregate = _profile_aggregate(primary)
    aggregate["scenario_counts"]["profiled_case_total"] = len(results)
    profiles = []
    for profile in PROFILES:
        rows = [row for row in results if row["profile"] == profile]
        if rows:
            profile_aggregate = _profile_aggregate(rows)
            profile_metrics = profile_aggregate["metrics"]
            profiles.append({
                "profile": profile,
                **profile_aggregate["scenario_counts"],
                "sensitivity": profile_metrics["founder_segment_detection_sensitivity"]["estimate"],
                "specificity": profile_metrics["founder_segment_detection_specificity"]["estimate"],
                "target_containment": profile_metrics["target_containment_rate_when_detectable"]["estimate"],
                "boundary_mean_absolute_error_cm": profile_metrics["boundary_errors"]["boundary_mean_absolute_error"],
            })
    positives = [row for row in results if row["truth_positive"]]
    anomalies = []
    non_detectable = sum(not row["detectable"] for row in positives)
    missing_families = sum(row["missing_family_count"] for row in positives)
    ibd_count = sum(row["ibd_claimed"] for row in results)
    if non_detectable:
        anomalies.append(f"positive_non_detectable_count={non_detectable}")
    if missing_families:
        anomalies.append(f"missing_or_excluded_family_count={missing_families}")
    if ibd_count:
        anomalies.append(f"ibs_reclassified_as_ibd_count={ibd_count}")
    return {
        **aggregate,
        "profile_results": profiles,
        "anomalies": anomalies,
    }


def acceptance_results(metrics: dict[str, Any]) -> list[dict[str, Any]]:
    """Raccorde les cinq métriques D04 sans sceller les seuils du pilote."""
    mapping = (
        ("CAL-FOUND-001", "founder_segment_detection_sensitivity"),
        ("CAL-FOUND-002", "founder_segment_detection_specificity"),
        ("CAL-FOUND-003", "target_containment_rate_when_detectable"),
        ("CAL-FOUND-004", "boundary_mean_absolute_error"),
    )
    decisions = []
    for criterion_id, name in mapping:
        measure = metrics["boundary_errors"] if criterion_id == "CAL-FOUND-004" else metrics[name]
        value = measure[name] if criterion_id == "CAL-FOUND-004" else measure["estimate"]
        denominator = measure["evaluated_boundary_count"] if criterion_id == "CAL-FOUND-004" else measure["denominator"]
        decisions.append({
            "criterion_id": criterion_id,
            "metric": name,
            "value": value,
            "denominator": denominator,
            "threshold": None,
            "lifecycle": "pilot_then_seal",
            "status": "NOT_EVALUABLE" if value is None else "REVIEW_REQUIRED",
        })
    count = metrics["ibs_reclassified_as_ibd_without_explicit_ibd_count"]
    decisions.append({
        "criterion_id": "CAL-FOUND-005",
        "metric": "ibs_reclassified_as_ibd_without_explicit_ibd_count",
        "value": count,
        "denominator": None,
        "threshold": 0,
        "lifecycle": "fixed_contract",
        "status": "PASS" if count == 0 else "FAIL",
    })
    return decisions


def validate_case_balance(results: list[dict[str, Any]], config: dict[str, Any]) -> None:
    """Exige la matrice scénario × profil sans pseudo-réplication."""
    definitions = {row["scenario_id"]: row for row in config["scenarios"]}
    if len(definitions) != len(config["scenarios"]):
        raise ValueError("d24_scenario_definition_duplicate")
    seeds = [row["random_seed"] for row in config["scenarios"]]
    if len(seeds) != len(set(seeds)):
        raise ValueError("d24_scenario_seed_duplicate")
    expected_profiles = {config["primary_profile"], *config["sensitivity_profiles"]}
    observed_pairs = [(row["scenario_id"], row["profile"]) for row in results]
    if len(observed_pairs) != len(set(observed_pairs)):
        raise ValueError("d24_scenario_profile_duplicate")
    by_scenario: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in results:
        definition = definitions.get(row["scenario_id"])
        if definition is None:
            raise ValueError("d24_scenario_not_prespecified")
        if (
            row["truth_role"] != definition["truth_role"]
            or row["age_generations"] != definition["age_generations"]
        ):
            raise ValueError("d24_scenario_definition_mismatch")
        by_scenario[row["scenario_id"]].append(row)
    if set(by_scenario) != set(definitions):
        raise ValueError("d24_scenario_matrix_incomplete")
    for rows in by_scenario.values():
        if {row["profile"] for row in rows} != expected_profiles:
            raise ValueError("d24_scenario_profile_matrix_incomplete")
        if len({row["truth_signature"] for row in rows}) != 1:
            raise ValueError("d24_paired_truth_changed_between_profiles")
        if len({row["source_manifest_sha256"] for row in rows}) != 1:
            raise ValueError("d24_truth_source_changed_between_profiles")
    counts = Counter(row["truth_role"] for row in definitions.values())
    if counts["POSITIVE"] < config["minimum_positive_scenarios"]:
        raise ValueError("d24_positive_scenario_count_insufficient")
    if counts["NEGATIVE"] < config["minimum_negative_scenarios"]:
        raise ValueError("d24_negative_scenario_count_insufficient")
