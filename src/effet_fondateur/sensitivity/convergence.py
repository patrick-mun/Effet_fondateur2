"""Convergence exploratoire des familles porteuses dans les PCA 1000G."""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from typing import Any

from effet_fondateur.audit import atomic_write_json
from effet_fondateur.contracts import validate_json_document, validate_tsv_table


def _mean_vector(vectors: list[tuple[float, ...]]) -> tuple[float, ...]:
    return tuple(sum(vector[index] for vector in vectors) / len(vectors) for index in range(10))


def _mean_pairwise_distance(vectors: list[tuple[float, ...]]) -> float:
    distances = [math.dist(left, right) for left, right in itertools.combinations(vectors, 2)]
    return sum(distances) / len(distances)


def publish_population_convergence(
    *, scores_path: Path, cohorts_path: Path, output_path: Path,
    global_matching_tolerance: float = 0.10,
) -> dict[str, Any]:
    """Compare trois unités porteuses à tous les triplets témoins indépendants."""
    scores = validate_tsv_table(scores_path, "ancestry_scores.schema.json")
    cohorts = validate_tsv_table(cohorts_path, "cohorts_frozen.schema.json")
    carriers = {
        row["SAMPLE_ID"] for row in cohorts.rows
        if row["COHORT_ID"] == "target_carriers_independent" and row["INCLUDED"]
    }
    controls = {
        row["SAMPLE_ID"] for row in cohorts.rows
        if row["COHORT_ID"] == "controls_unrelated" and row["INCLUDED"]
    }
    result: dict[str, Any] = {
        "schema_version": "1.0.0", "method_id": "global_local_triplet_convergence_v1",
        "analysis_role": "EXPLORATORY_POST_HOC", "component_count": 10,
        "distance": "MEAN_PAIRWISE_EUCLIDEAN", "null_mode": "EXHAUSTIVE_CONTROL_TRIPLETS",
        "random_seed": None, "global_matching_relative_tolerance": global_matching_tolerance,
        "carrier_unit_count": len(carriers), "control_unit_count": len(controls),
        "local_entity_definition": "MEAN_OF_PROJECTED_HAPLOTYPE_SCORES_PER_SAMPLE",
        "status": "NOT_EVALUATED", "global": None, "local": None, "conditional_local": None,
        "negative_control_windows": {"status": "NOT_EVALUATED", "reason": "COMPARABLE_WINDOWS_NOT_PRESPECIFIED"},
        "interpretation": "Analyse secondaire exploratoire post hoc ; aucune attribution ethnique, généalogique ou causale.",
    }
    if len(carriers) != 3 or len(controls) < 3:
        atomic_write_json(output_path, result); validate_json_document(result, "population_convergence.schema.json")
        return result
    vectors: dict[str, dict[str, list[tuple[float, ...]]]] = {"GLOBAL": {}, "LOCAL": {}}
    for row in scores.rows:
        sample_id = row["SAMPLE_ID"]
        scope = row["ANALYSIS_SCOPE"]
        if not row["PROJECTED"] or sample_id not in carriers | controls or scope not in vectors:
            continue
        vector_values = tuple(row[f"PC{index}"] for index in range(1, 11))
        if any(value is None for value in vector_values):
            continue
        vectors[scope].setdefault(sample_id, []).append(tuple(float(value) for value in vector_values))
    if any(sample not in vectors[scope] for scope in vectors for sample in carriers | controls):
        atomic_write_json(output_path, result); validate_json_document(result, "population_convergence.schema.json")
        return result
    points = {
        scope: {sample: _mean_vector(sample_vectors) for sample, sample_vectors in scope_vectors.items()}
        for scope, scope_vectors in vectors.items()
    }
    ordered_carriers, ordered_controls = sorted(carriers), sorted(controls)
    observed = {
        scope: _mean_pairwise_distance([points[scope][sample] for sample in ordered_carriers])
        for scope in points
    }
    null_distances: list[tuple[float, float]] = []
    for triplet in itertools.combinations(ordered_controls, 3):
        null_distances.append(tuple(
            _mean_pairwise_distance([points[scope][sample] for sample in triplet])
            for scope in ("GLOBAL", "LOCAL")
        ))
    def empirical(index: int, candidates: list[tuple[float, float]]) -> dict[str, Any]:
        exceedance = sum(row[index] <= observed[("GLOBAL", "LOCAL")[index]] for row in candidates)
        return {"observed_distance": observed[("GLOBAL", "LOCAL")[index]], "triplet_count": len(candidates),
                "as_or_more_compact_count": exceedance, "empirical_probability": (exceedance + 1) / (len(candidates) + 1)}
    matched = [row for row in null_distances if abs(row[0] - observed["GLOBAL"]) / observed["GLOBAL"] <= global_matching_tolerance]
    result.update({"status": "EVALUATED", "global": empirical(0, null_distances),
                   "local": empirical(1, null_distances),
                   "conditional_local": empirical(1, matched) if matched else None})
    atomic_write_json(output_path, result)
    validate_json_document(result, "population_convergence.schema.json")
    return result
