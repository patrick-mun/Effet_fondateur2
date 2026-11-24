"""Agrégats graphiques D22 issus des sorties closes, sans recalcul d'analyse.

Aucun sujet, score individuel ou génotype n'est exporté. Les barycentres et
écarts-types sont des descriptions pour figures, pas de nouvelles inférences.
"""

from __future__ import annotations

from collections import Counter
import math
from pathlib import Path
from typing import Any

import numpy as np

from effet_fondateur.calibration.evaluate_structure import _closed_upstream, _component_paths
from effet_fondateur.calibration.evaluate_technical import raw_table
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.structure_bench import CASES, COHORTS, checked_document
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError


def aggregate_population(scores: list[dict[str, str]], outliers: list[dict[str, str]], truth: dict[str, Any],
                         alpha: float) -> dict[str, Any]:
    """Barycentres ± écart-type et histogrammes communs ; pas de scores publiés."""
    if [row["SAMPLE_ID"] for row in scores] != [row["SAMPLE_ID"] for row in outliers] or len(scores) != len(truth["groups"]):
        raise TechnicalEvaluationError("d22_visual_sample_axis_mismatch")
    points = np.array([[float(row["PC1"]), float(row["PC2"])] for row in scores])
    groups = np.array(truth["groups"])
    reference = np.array([row["REFERENCE_INCLUDED"] == "true" for row in scores])
    if int(reference.sum()) != truth["reference_count"]:
        raise TechnicalEvaluationError("d22_visual_reference_count_mismatch")
    centers = []
    for group in ("A", "B", "MIXED", "RARE"):
        for role, selected in (("REFERENCE", reference), ("PROJECTED", ~reference)):
            values = points[selected & (groups == group)]
            if len(values):
                centers.append(dict(group=group, role=role, count=len(values), mean=values.mean(axis=0).tolist(),
                                    standard_deviation=values.std(axis=0, ddof=1).tolist() if len(values) > 1 else [0., 0.]))
    distances = np.array([float(row["DISTANCE_SQUARED"]) for row in outliers])
    threshold = -2 * math.log(alpha)
    edges = np.linspace(0, max(float(distances.max()), threshold) * 1.05, 13)
    histogram = dict(edges=edges.tolist(), threshold=threshold, degrees_of_freedom=2,
                     reference_counts=np.histogram(distances[reference], bins=edges)[0].tolist(),
                     projected_counts=np.histogram(distances[~reference], bins=edges)[0].tolist())
    return dict(centroids=centers, distance_histogram=histogram, sample_count=len(scores))


def aggregate_cohorts(rows: list[dict[str, str]], sample_count: int) -> dict[str, Any]:
    """Compte les décisions parmi tous les sujets, pas seulement les candidats."""
    result = {}
    for cohort in COHORTS:
        members = [row for row in rows if row["COHORT_ID"] == cohort]
        if len(members) != sample_count or len({row["SAMPLE_ID"] for row in members}) != sample_count:
            raise TechnicalEvaluationError("d22_visual_cohort_axis_mismatch")
        included = sum(row["INCLUDED"] == "true" for row in members)
        reasons = Counter(row["EXCLUSION_CODE"] for row in members if row["INCLUDED"] == "false")
        if "" in reasons or included + sum(reasons.values()) != sample_count:
            raise TechnicalEvaluationError("d22_visual_exclusion_partition_invalid")
        result[cohort] = dict(included=included, not_included=sample_count - included, exclusion_counts=dict(sorted(reasons.items())))
    return result


def collect_visual_evidence(*, panel_root: Path, pipeline_root: Path, evaluation_root: Path) -> dict[str, Any]:
    """Extrait après clôture et vérifie les liaisons aux verdicts déjà publiés."""
    panel = checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    receipt = checked_document(pipeline_root, "run_index.json", "closure_state", "CLOSED")
    evaluation = checked_document(evaluation_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
    fingerprints = dict(panel_manifest_sha256=_sha(panel_root / "panel_manifest.json"), run_index_sha256=_sha(pipeline_root / "run_index.json"),
                        evaluation_manifest_sha256=_sha(evaluation_root / "evaluation_manifest.json"))
    if evaluation["panel_manifest_sha256"] != fingerprints["panel_manifest_sha256"] or receipt["panel_manifest_sha256"] != fingerprints["panel_manifest_sha256"] or evaluation["run_index_sha256"] != fingerprints["run_index_sha256"]:
        raise TechnicalEvaluationError("d22_visual_sources_not_bound")
    truth = _json(panel_root / "private_truth/structure_truth.json")
    reports = {stage: _json(evaluation_root / ("CAL-REP-" + stage + ".json")) for stage in ("08", "09")}
    cases, sources = [], []
    for index, entry in enumerate(receipt["runs"]):
        run = _closed_upstream(pipeline_root, entry)
        paths = _component_paths(run, "stages/08_analyze_population_structure", "08", "analyze_population_structure")
        cohort_paths = _component_paths(run, entry["freeze_directories"][0], "09", "freeze_cohorts")
        population = aggregate_population(raw_table(paths["population_scores"]), raw_table(paths["population_outliers"]), truth,
                                          evaluation["configuration"]["population_parameters"]["outlier_alpha"])
        cohorts = aggregate_cohorts(raw_table(cohort_paths["cohorts_frozen"]), population["sample_count"])
        if any(cohorts[cohort]["included"] != reports["09"]["evaluations"][index]["metrics"]["included_counts"][cohort] for cohort in COHORTS):
            raise TechnicalEvaluationError("d22_visual_cohort_counts_differ_from_evaluation")
        eigenvalues = [dict(component=row["COMPONENT"], eigenvalue=float(row["EIGENVALUE"]), explained_variance_ratio=float(row["EXPLAINED_VARIANCE_RATIO"]))
                       for row in raw_table(paths["population_eigenvalues"])]
        cases.append(dict(case=CASES[index], **population, cohorts=cohorts, eigenvalues=eigenvalues))
        for name, path in [(name, paths[name]) for name in ("population_scores", "population_outliers", "population_eigenvalues")] + [("cohorts_frozen", cohort_paths["cohorts_frozen"])]:
            sources.append(dict(bundle_id=entry["bundle_id"], artifact_id=name, sha256=_sha(path)))
    for root, filename, state_key, state in ((panel_root, "panel_manifest.json", "manifest_state", "FINALIZED"),
                                           (pipeline_root, "run_index.json", "closure_state", "CLOSED"), (evaluation_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")):
        checked_document(root, filename, state_key, state)
    if fingerprints != dict(panel_manifest_sha256=_sha(panel_root / "panel_manifest.json"), run_index_sha256=_sha(pipeline_root / "run_index.json"), evaluation_manifest_sha256=_sha(evaluation_root / "evaluation_manifest.json")):
        raise TechnicalEvaluationError("d22_visual_sources_changed")
    return dict(schema_version="1.0.0", presentation_version="2.0", method_id=evaluation["method_id"], source_fingerprints=fingerprints,
                source_artifacts=sources, cases=cases, aggregation_only=True, individual_identifiers_exported=False, analyses_rerun=False)
