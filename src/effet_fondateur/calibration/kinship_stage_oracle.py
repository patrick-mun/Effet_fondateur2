"""D21 : conformité 07, vérité des degrés et proposition indépendante.

Les estimations sont contrôlées depuis les génotypes ; la précision de
classification contre le pedigree reste un résultat distinct. Les sorties
ne conservent aucune clé individuelle, seulement des agrégats par unité.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from effet_fondateur.calibration.evaluate_qc import compare_metrics
from effet_fondateur.calibration.evaluate_technical import raw_table
from effet_fondateur.calibration.kinship_oracle import classification_metrics, degree_and_relation, independent_proposal, king_metrics
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.simulation.kinship_bench import DEGREES
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix


def check_pairs(matrix: GenotypeMatrix, rows: list[dict[str, str]], registry: list[dict[str, str]], parameters: dict[str, Any], tolerance: float) -> tuple[dict[str, Any], dict[tuple[str, str], dict[str, str]], dict[tuple[str, str], set[str]]]:
    """Vérifie les lignes KING rapportées ; aucune paire absente n'est imputée."""
    indices = {sample: index for index, (_, sample) in enumerate(matrix.samples)}
    samples = {row["SAMPLE_ID"]: row for row in registry if row["SAMPLE_ID"] in indices}
    by_plink = {(row["FID"], row["IID"]): row["SAMPLE_ID"] for row in samples.values()}
    declared = set()
    for sample, row in samples.items():
        for field in ("PID", "MID"):
            parent = by_plink.get((row["FID"], row[field])) if row[field] != "0" else None
            if parent:
                declared.add(tuple(sorted((sample, parent))))
    expected, numeric_errors, malformed = [], 0, 0
    pairs, edges = {}, {pair: {"DECLARED_PARENT_CHILD"} for pair in declared}
    for row in rows:
        pair = (row["SAMPLE_ID_1"], row["SAMPLE_ID_2"])
        if pair[0] >= pair[1] or pair in pairs or any(sample not in indices for sample in pair):
            malformed += 1
            continue
        pairs[pair] = row
        metrics = king_metrics(matrix, indices[pair[0]], indices[pair[1]])
        for field, value in metrics.items():
            if value is None:
                numeric_errors += 1
            elif field == "N_SNP":
                numeric_errors += row[field] != str(value)
            else:
                numeric_errors += not np.isclose(float(row[field]), value, rtol=0, atol=tolerance)
        degree, relation = degree_and_relation(float(row["KINSHIP"]), float(row["IBS0"]), parameters)
        same_family = samples[pair[0]]["FID"] == samples[pair[1]]["FID"]
        declared_relation = "PARENT_CHILD" if pair in declared else "SAME_FAMILY_UNSPECIFIED" if same_family else "NONE"
        if declared_relation == "PARENT_CHILD":
            concordance = "DISCORDANT_DEGREE" if degree != "FIRST" else "CONCORDANT" if relation == "PARENT_CHILD_COMPATIBLE" else "DISCORDANT_RELATION"
        else:
            concordance = "FAMILY_RELATION_UNRESOLVED" if same_family and degree != "UNRELATED" else "CRYPTIC_RELATEDNESS" if degree != "UNRELATED" else "NO_DECLARED_RELATION"
        expected.append(dict(PAIR_ID=f"kinship_pair_{len(expected)+1:06d}", SAMPLE_ID_1=pair[0], SAMPLE_ID_2=pair[1],
                             FID_1=samples[pair[0]]["FID"], IID_1=pair[0], FID_2=samples[pair[1]]["FID"], IID_2=pair[1],
                             INFERRED_DEGREE=degree, INFERRED_RELATION=relation, DECLARED_RELATION=declared_relation,
                             PAIR_ORIGIN="WITHIN_FAMILY" if same_family else "BETWEEN_FAMILY", CONCORDANCE_STATUS=concordance))
        if degree != "UNRELATED" and DEGREES.index(degree) <= parameters["relatedness_max_degree_for_independence"]:
            edges.setdefault(pair, set()).add("KING_RELATEDNESS")
    row_comparison = compare_metrics(expected, rows, "PAIR_ID")
    row_comparison.update(numeric=numeric_errors, malformed=malformed)
    if list(pairs) != sorted(pairs):
        row_comparison["order"] += 1
    return row_comparison, pairs, edges


def _component_labels(samples: list[str], neighbors: dict[str, set[str]]) -> dict[str, str]:
    """Étiquette les composantes dans l'ordre des sujets, indépendamment de 07."""
    remaining, result, index = set(samples), {}, 0
    for sample in samples:
        if sample not in remaining:
            continue
        index += 1
        frontier = [sample]
        while frontier:
            current = frontier.pop()
            if current not in remaining:
                continue
            remaining.remove(current)
            result[current] = f"kinship_component_{index:06d}"
            frontier.extend(neighbors[current] & remaining)
    return result


def evaluate_kinship_stage(matrix: GenotypeMatrix, registry: list[dict[str, str]], qc_rows: list[dict[str, str]],
                           artifacts: dict[str, Path], units: list[dict[str, Any]], original_sample_ids: tuple[str, ...],
                           config: dict[str, Any], *, origins: np.ndarray | None = None,
                           before_pruning: GenotypeMatrix | None = None) -> dict[str, Any]:
    """Mesure 07 contre génotypes et pedigree, sans approuver les exclusions."""
    parameters = config["kinship_parameters"]
    actual_pairs = raw_table(artifacts["kinship_pairs"])
    comparisons, pairs, edges = check_pairs(matrix, actual_pairs, registry, parameters, config["king_absolute_tolerance"])
    reasons = ["pairs_" + key for key, count in comparisons.items() if key != "compared" and count]
    samples = [sample for _, sample in matrix.samples]
    source_metadata = {row["SAMPLE_ID"]: row for row in registry}
    by_plink = {(row["FID"], row["IID"]): row["SAMPLE_ID"] for row in registry if row["SAMPLE_ID"] in set(samples)}
    proposal, neighbors = independent_proposal(samples, edges, qc_rows)
    labels = _component_labels(samples, neighbors)
    qc = {row["SAMPLE_ID"]: row for row in qc_rows}
    expected_proposal = []
    for sample in samples:
        selected = sample in proposal["selected"]
        representative = proposal["representatives"].get(sample, "")
        sources = edges[tuple(sorted((sample, representative)))] if representative else set()
        basis = "KING_AND_DECLARED" if len(sources) == 2 else next(iter(sources)) if sources else ""
        expected_proposal.append(dict(SAMPLE_ID=sample, FID=source_metadata[sample]["FID"], IID=sample,
                                      COMPONENT_ID=labels[sample], PROPOSED_INCLUDED=str(selected).lower(),
                                      PROPOSAL_STATUS="RETAIN" if selected else "EXCLUDE_RELATED", REPRESENTATIVE_SAMPLE_ID=representative,
                                      RELATION_BASIS=basis, F_MISS=float(qc[sample]["F_MISS"]),
                                      ABS_HETEROZYGOSITY_Z=abs(float(qc[sample]["HETEROZYGOSITY_Z"])) if qc[sample]["HETEROZYGOSITY_Z"] else None,
                                      DECISION_SOURCE="stage_07_quality_aware_greedy_maximal", MANUAL_VALIDATION_REQUIRED=str(not selected).lower()))
    actual_proposal = raw_table(artifacts["independent_set_proposal"])
    proposal_comparison = compare_metrics(expected_proposal, actual_proposal, "SAMPLE_ID")
    reasons += ["proposal_" + key for key, count in proposal_comparison.items() if key != "compared" and count]
    selected = {row["SAMPLE_ID"] for row in actual_proposal if row["PROPOSED_INCLUDED"] == "true"}
    graph_violations = sum(left in selected and right in selected for left, right in edges)
    nonmaximal = sum(sample not in selected and not neighbors[sample] & selected for sample in samples)
    if graph_violations or nonmaximal:
        reasons.append("observed_graph_independence_or_maximality_failed")
    expected_pedigree = []
    for sample in samples:
        child = source_metadata[sample]
        for role, field in (("FATHER", "PID"), ("MOTHER", "MID")):
            if child[field] == "0":
                continue
            parent = by_plink.get((child["FID"], child[field]), "")
            pair = pairs.get(tuple(sorted((sample, parent)))) if parent else None
            status = pair["CONCORDANCE_STATUS"] if pair else "NOT_EVALUATED_KING_PAIR_MISSING" if parent else "NOT_EVALUATED_PARENT_ABSENT"
            expected_pedigree.append(dict(RELATION_ID=f"pedigree_relation_{len(expected_pedigree)+1:06d}", CHILD_SAMPLE_ID=sample,
                                         PARENT_ROLE=role, DECLARED_PARENT_IID=child[field], PARENT_SAMPLE_ID=parent,
                                         INFERRED_DEGREE=pair["INFERRED_DEGREE"] if pair else "", INFERRED_RELATION=pair["INFERRED_RELATION"] if pair else "",
                                         KINSHIP=float(pair["KINSHIP"]) if pair else None, IBS0=float(pair["IBS0"]) if pair else None, CONCORDANCE_STATUS=status))
    pedigree = compare_metrics(expected_pedigree, raw_table(artifacts["pedigree_concordance"]), "RELATION_ID")
    reasons += ["pedigree_" + key for key, count in pedigree.items() if key != "compared" and count]
    summary = []
    for degree in DEGREES:
        members = [row for row in actual_pairs if row["INFERRED_DEGREE"] == degree]
        summary.append(dict(INFERRED_DEGREE=degree, PAIR_COUNT=len(members), DECLARED_PAIR_COUNT=sum(row["DECLARED_RELATION"] != "NONE" for row in members),
                            CRYPTIC_PAIR_COUNT=sum(row["CONCORDANCE_STATUS"] == "CRYPTIC_RELATEDNESS" for row in members),
                            CONCORDANT_PAIR_COUNT=sum(row["CONCORDANCE_STATUS"] == "CONCORDANT" for row in members),
                            DISCORDANT_PAIR_COUNT=sum(row["CONCORDANCE_STATUS"].startswith("DISCORDANT") for row in members)))
    degree_comparison = compare_metrics(summary, raw_table(artifacts["kinship_degree_summary"]), "INFERRED_DEGREE")
    reasons += ["degree_summary_" + key for key, count in degree_comparison.items() if key != "compared" and count]
    functional_reasons = list(reasons)
    records, before_records, diagnostics, lost_units, forbidden, cryptic_positive = [], [], [], 0, 0, 0
    indices = {sample: index for index, (_, sample) in enumerate(matrix.samples)}
    pre_indices = {sample: index for index, (_, sample) in enumerate(before_pruning.samples)} if before_pruning is not None else {}
    truth_degree_by_pair = {}
    for unit in units:
        members = [original_sample_ids[index] for index in unit["sample_indices"]]
        if not set(members) <= set(samples):
            lost_units += 1
            continue
        pair = pairs.get(tuple(sorted(members)))
        truth_degree_by_pair[tuple(sorted(members))] = unit["degree"]
        observed = pair["INFERRED_DEGREE"] if pair else "MISSING"
        records.append(dict(true_degree=unit["degree"], observed_degree=observed, stratum_id=unit["stratum_id"]))
        forbidden += unit["degree"] != "UNRELATED" and set(members) <= selected
        cryptic_positive += bool(unit["cryptic"] and pair is not None)
        after_metrics = king_metrics(matrix, indices[members[0]], indices[members[1]])
        realized = None
        if origins is not None:
            left, right = indices[members[0]], indices[members[1]]
            called = (matrix.alleles[:, left, 0] >= 0) & (matrix.alleles[:, right, 0] >= 0)
            realized = float(np.mean(np.sum(origins[called, left, :, None] == origins[called, right, None, :], axis=(1, 2)) / 4))
        before_phi = None
        if before_pruning is not None and set(members) <= set(pre_indices):
            before_metrics = king_metrics(before_pruning, pre_indices[members[0]], pre_indices[members[1]])
            before_phi = before_metrics["KINSHIP"]
            before_degree = degree_and_relation(before_phi, before_metrics["IBS0"], parameters)[0] if before_phi is not None else "MISSING"
            before_records.append(dict(true_degree=unit["degree"], observed_degree=before_degree, stratum_id=unit["stratum_id"]))
        diagnostics.append(dict(degree=unit["degree"], stratum_id=unit["stratum_id"], theoretical=unit["theoretical_phi"],
                                realized=realized, estimate=after_metrics["KINSHIP"], before=before_phi))
    if lost_units:
        reasons.append("target_units_lost_after_qc")
        functional_reasons.append("target_units_lost_after_qc")
    if forbidden:
        reasons.append("true_related_pair_in_independent_proposal")
    performance = classification_metrics(records, config)
    stratified = {stratum["stratum_id"]: classification_metrics([record for record in records if record["stratum_id"] == stratum["stratum_id"]], config) for stratum in config["strata"]}
    if performance["missing_pair_count"]:
        reasons.append("mandatory_king_pair_missing")
        functional_reasons.append("mandatory_king_pair_missing")
    report = _json(artifacts["kinship_report"])
    if report["proposed_inclusion_count"] != len(selected) or report["reported_pair_count"] != len(actual_pairs) or not report["manual_validation_required"]:
        reasons.append("kinship_report_count_or_review_mismatch")
        functional_reasons.append("kinship_report_count_or_review_mismatch")
    estimates = {}
    for degree in DEGREES:
        group = [row for row in diagnostics if row["degree"] == degree]
        estimates[degree] = {field + "_mean": float(np.mean([row[field] for row in group if row[field] is not None]))
                            if any(row[field] is not None for row in group) else None for field in ("theoretical", "realized", "estimate", "before")}
    false_related = sum(row["INFERRED_DEGREE"] != "UNRELATED" and truth_degree_by_pair.get(pair, "UNRELATED") == "UNRELATED" for pair, row in pairs.items())
    return dict(status="FAIL" if reasons else "PASS", reasons=sorted(set(reasons)),
                functional_concordance="FAIL" if functional_reasons else "PASS",
                comparisons=dict(pairs=comparisons, proposal=proposal_comparison, pedigree=pedigree, degree_summary=degree_comparison),
                metrics=dict(reported_pairs=len(actual_pairs), target_units=len(records), lost_units=lost_units, cryptic_positive_pairs_detected=cryptic_positive,
                             independent_proposal_count=len(selected), forbidden_true_pair_count=int(forbidden), observed_graph_violation_count=graph_violations,
                             maximality_error_count=nonmaximal, declared_parent_rows=len(expected_pedigree),
                             reported_false_related_pairs=false_related),
                performance=performance, stratified_performance=stratified,
                estimation_diagnostics=estimates,
                pre_pruning_oracle_performance=classification_metrics(before_records, config) if before_pruning is not None else None,
                performance_decision="REVIEW_REQUIRED", frozen_cohort_09_evaluated=False)
