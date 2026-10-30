"""Évaluateurs D20 05/10 après clôture, sans lancement de processus.

Concordance au contrat et détection du clone sont séparées. Tous les taux
restent descriptifs : cas appariés IID, aucun intervalle de qualification.
Les rapports publiés ne contiennent que des agrégats, jamais une clé de sujet.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import itertools
import json
import math
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np

from effet_fondateur.calibration.evaluate_technical import _dataset, contained_file, raw_table
from effet_fondateur.calibration.materialize_ibd import _json, _sha, _text
from effet_fondateur.calibration.qc_bench import METHOD_ID, PRODUCER, checked_inventory, generate_cases, file_inventory, with_target, write_json
from effet_fondateur.calibration.qc_oracle import COHORTS, QcExpectation, confusion, final_expectation, preliminary_expectation, subset
from effet_fondateur.calibration.run_qc_bench import checked_stage
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, TechnicalEvaluationError, compare_genotypes, observation_matrix
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document, validate_tsv_table


RELATIVE_TOLERANCE = 5e-4
ABSOLUTE_TOLERANCE = 5e-6


def verify_final_input_binding(inputs: dict[str, Any], audit: dict[str, Any], parameters: dict[str, Any]) -> None:
    """Relie les entrées consommées ; distingue les annexes inutilisées du v1.

    Le premier runner a aussi déclaré trois rapports 04 non consommés par 10.
    Leur présence n'est pas une différence de génotypes ou de paramètres.
    Ils restent empreintés dans le reçu, et ne sont pas pris comme preuves QC.
    Aucune autre entrée inattendue et aucune divergence consommée n'est admise.
    """
    declared = {item["artifact_id"]: item for item in inputs["artifacts"]}
    consumed = {item["artifact_id"]: item for item in audit["inputs"]}
    unused = {"target_variant_injection_report", "target_variant_mendel", "target_variant_mendel_report"}
    if len(declared) != len(inputs["artifacts"]) or len(consumed) != len(audit["inputs"]) or not consumed.keys() <= declared.keys() or declared.keys() - consumed.keys() - unused or any(
        declared[key] != item for key, item in consumed.items()) or inputs["parameters"] != parameters or inputs["signature"] != audit["signature"]:
        raise TechnicalEvaluationError("d20_final_input_audit_binding_invalid")


def compare_metrics(expected: list[dict[str, Any]], actual: list[dict[str, str]], key: str) -> dict[str, int]:
    """Vérifie exhaustivité, unicité, ordre, décisions exactes et mesures arrondies.

    Tolérances numériques : 5e-4 relative, 5e-6 absolue (zéro absolu pour
    HWE_P, afin de ne pas assimiler une très petite p-value à zéro). Les codes,
    booléens, effectifs entiers et statuts ne bénéficient d'aucune tolérance.
    """
    expected_keys = [row[key] for row in expected]
    actual_keys = [row[key] for row in actual]
    expected_by_key = {row[key]: row for row in expected}
    counts = Counter(actual_keys)
    output = dict(lost=len(set(expected_keys) - set(actual_keys)), extra=len(set(actual_keys) - set(expected_keys)),
                  duplicate=sum(count - 1 for count in counts.values()),
                  order=sum(left != right for left, right in zip(expected_keys, actual_keys)) + abs(len(expected_keys) - len(actual_keys)),
                  decision=0, numeric=0, compared=len(actual))
    for row in actual:
        if row[key] not in expected_by_key:
            continue
        for field, expected_value in expected_by_key[row[key]].items():
            observed = row.get(field)
            if isinstance(expected_value, float) or expected_value is None:
                try:
                    matched = observed == "" if expected_value is None or not math.isfinite(expected_value) else math.isfinite(float(observed)) and math.isclose(float(observed), expected_value, rel_tol=RELATIVE_TOLERANCE, abs_tol=0.0 if field == "HWE_P" else ABSOLUTE_TOLERANCE)
                except (ValueError, TypeError):
                    matched = False
                output["numeric"] += not matched
            else:
                output["decision"] += observed != str(expected_value)
    return output


def compare_preliminary_alerts(expected: QcExpectation, actual: list[dict[str, str]], parameters: dict[str, Any]) -> dict[str, int]:
    """Contrôle tous les motifs non-PI_HAT, leur portée, statut, valeur et seuil."""
    oracle = []
    definitions = ((expected.samples, "IID", "SAMPLE", "EXCLUSION_CODE", "sample_missingness", "F_MISS", "F_MISS", "sample_missingness_max", "EXCLUDED"),
                   (expected.samples, "IID", "SAMPLE", "QC_STATUS", "ALERT", "HETEROZYGOSITY_Z", "HETEROZYGOSITY_Z", "heterozygosity_z_alert", "ALERT"),
                   (expected.variants, "VARIANT_ID", "VARIANT", "EXCLUSION_CODE", "variant_missingness", "F_MISS", "F_MISS", "variant_missingness_max", "EXCLUDED"),
                   (expected.variants, "VARIANT_ID", "VARIANT", "ALERT_CODES", "low_maf_for_kinship", "MAF", "MAF", "kinship_maf_alert_threshold", "ALERT"),
                   (expected.variants, "VARIANT_ID", "VARIANT", "ALERT_CODES", "batch_missingness_delta", "BATCH_MISSINGNESS_DELTA", "BATCH_MISSINGNESS_DELTA", "batch_missingness_delta_alert", "ALERT"))
    for rows, key, scope, field, code, metric, value_field, threshold, status in definitions:
        for row in rows:
            if code not in row[field].split(","):
                continue
            reason = "heterozygosity_outlier" if metric == "HETEROZYGOSITY_Z" else code
            oracle.append(dict(SCOPE=scope, ENTITY_ID=row[key], METRIC=metric, OBSERVED_VALUE=row[value_field],
                               THRESHOLD=float(parameters[threshold]), STATUS=status, REASON_CODE=reason))
    oracle.extend((dict(SCOPE="GLOBAL", ENTITY_ID="genomewide_pre_qc", METRIC="SEX_DISCORDANCE", OBSERVED_VALUE="", THRESHOLD="",
                        STATUS="NOT_EVALUATED", REASON_CODE="autosomal_dataset_without_sex_chromosomes"),
                   dict(SCOPE="GLOBAL", ENTITY_ID="genomewide_pre_qc", METRIC="HWE_FILTER", OBSERVED_VALUE="not_applied", THRESHOLD="",
                        STATUS="INFO", REASON_CODE="cohorts_not_frozen")))
    actual = [dict(row) for row in actual if row["REASON_CODE"] != "potential_technical_duplicate"]
    for rows in (oracle, actual):
        for row in rows:
            row["RULE_KEY"] = "|".join(str(row[key]) for key in ("SCOPE", "ENTITY_ID", "REASON_CODE"))
        rows.sort(key=lambda row: row["RULE_KEY"])
    return compare_metrics(oracle, actual, "RULE_KEY")


def _paths(root: Path, directory: str) -> dict[str, Path]:
    """Résout les seuls artefacts publiés et valides, sans ouvrir une source réelle."""
    artifacts = checked_stage(root, directory)
    paths = {}
    for artifact_id, artifact in artifacts.items():
        path = contained_file(root, artifact["path"])
        if artifact["schema_name"]:
            schema_name = artifact["schema_name"]
            if Path(schema_name).name != schema_name:
                raise TechnicalEvaluationError("d20_schema_path_invalid")
            schema = _json(Path(__file__).resolve().parents[3] / "schemas" / schema_name)
            if "x-tsv" in schema:
                validate_tsv_table(path, schema_name)
            else:
                validate_json_document(_json(path), schema_name)
        paths[artifact_id] = path
    return paths


def _rules(expected: QcExpectation, actual_samples: list[dict[str, str]], actual_variants: list[dict[str, str]],
           actual_batches: list[dict[str, str]], *, final: bool, cohort: str = "") -> dict[str, dict[str, Any]]:
    """Produit les matrices de conformité sur les unités propres à chaque règle."""
    definitions = [("sample_missingness", expected.samples, actual_samples, "IID", "EXCLUSION_CODE", "sample_missingness", "sample")]
    if final:
        definitions += [("variant_missingness", expected.variants, actual_variants, "VARIANT_ID", "FILTER_CODES", "variant_missingness", "variant")]
        if cohort == COHORTS[0]:
            definitions += [(code, expected.variants, actual_variants, "VARIANT_ID", "FILTER_CODES", code, "variant") for code in ("low_maf_controls", "hwe_controls")]
        if cohort == COHORTS[2]:
            definitions += [("low_maf_target_region", expected.variants, actual_variants, "VARIANT_ID", "FILTER_CODES", "low_maf_target_region", "variant"),
                            ("target_exception", expected.variants, actual_variants, "VARIANT_ID", "FILTER_EXCEPTION", "true", "variant")]
        definitions += [("batch_missingness_delta", expected.batches, actual_batches, "ARRAY_BATCH", "ALERT_CODE", "batch_missingness_delta", "batch")]
    else:
        definitions += [("heterozygosity_outlier", expected.samples, actual_samples, "IID", "QC_STATUS", "ALERT", "sample"),
                        ("variant_missingness", expected.variants, actual_variants, "VARIANT_ID", "EXCLUSION_CODE", "variant_missingness", "variant")]
        definitions += [(code, expected.variants, actual_variants, "VARIANT_ID", "ALERT_CODES", code, "variant")
                        for code in ("low_maf_for_kinship", "batch_missingness_delta")]
    output = {}
    for name, oracle_rows, observed_rows, key, field, code, unit in definitions:
        universe = {row[key] for row in oracle_rows}
        expected_positive = {row[key] for row in oracle_rows if code in str(row[field]).split(",")}
        actual_positive = {row[key] for row in observed_rows if code in row[field].split(",")}
        output[name] = dict(unit=unit, **confusion(expected_positive, actual_positive, universe))
    return output


def _metric_checks(expected: QcExpectation, actual_samples: list[dict[str, str]], actual_variants: list[dict[str, str]],
                   actual_batches: list[dict[str, str]]) -> tuple[dict[str, Any], list[str]]:
    """Agrège les comparaisons de tables et leurs motifs, sans données individuelles."""
    comparisons = dict(samples=compare_metrics(expected.samples, actual_samples, "IID"),
                       variants=compare_metrics(expected.variants, actual_variants, "VARIANT_ID"))
    if expected.batches:
        comparisons["batches"] = compare_metrics(expected.batches, actual_batches, "ARRAY_BATCH")
    reasons = [f"{table}_{metric}_mismatch" for table, counts in comparisons.items() for metric, count in counts.items() if metric != "compared" and count]
    return comparisons, reasons


def _evaluate_05(matrix: GenotypeMatrix, expected: QcExpectation, artifacts: dict[str, Path], pairs: list[list[str]], batches: list[str], parameters: dict[str, Any]) -> dict[str, Any]:
    """Évalue 05 : tables, alertes, clone connu, politiques et génotypes retenus."""
    samples = raw_table(artifacts["qc_individual_metrics"])
    variants = raw_table(artifacts["qc_variant_metrics"])
    comparisons, reasons = _metric_checks(expected, samples, variants, [])
    genotypes = compare_genotypes(expected.retained, _dataset(artifacts, "genomewide_pre_qc"), stage="05")
    reasons += genotypes["reasons"]
    rules = _rules(expected, samples, variants, [], final=False)
    alerts = raw_table(artifacts["qc_alerts"])
    comparisons["alerts"] = compare_preliminary_alerts(expected, alerts, parameters)
    reasons += ["alerts_" + key for key, value in comparisons["alerts"].items() if key != "compared" and value]
    known_reasons = {"sample_missingness", "variant_missingness", "heterozygosity_outlier", "low_maf_for_kinship", "batch_missingness_delta",
                     "potential_technical_duplicate", "autosomal_dataset_without_sex_chromosomes", "cohorts_not_frozen"}
    if any(row["REASON_CODE"] not in known_reasons for row in alerts):
        reasons.append("unknown_alert_reason")
    for name, field, rows, key, value, scope in (
        ("sample_missingness", "EXCLUSION_CODE", expected.samples, "IID", "sample_missingness", "SAMPLE"),
        ("heterozygosity_outlier", "QC_STATUS", expected.samples, "IID", "ALERT", "SAMPLE"),
        ("variant_missingness", "EXCLUSION_CODE", expected.variants, "VARIANT_ID", "variant_missingness", "VARIANT"),
        ("low_maf_for_kinship", "ALERT_CODES", expected.variants, "VARIANT_ID", "low_maf_for_kinship", "VARIANT"),
        ("batch_missingness_delta", "ALERT_CODES", expected.variants, "VARIANT_ID", "batch_missingness_delta", "VARIANT")):
        expected_entities = {row[key] for row in rows if value in row[field].split(",")}
        observed_entities = [row["ENTITY_ID"] for row in alerts if row["REASON_CODE"] == name and row["SCOPE"] == scope]
        if set(observed_entities) != expected_entities or len(observed_entities) != len(set(observed_entities)):
            reasons.append(name + "_alert_reason_incomplete")
    universe_pairs = {"|".join(sorted((f"{left}:{left}", f"{right}:{right}"))) for (_, left), (_, right) in itertools.combinations(matrix.samples, 2)}
    expected_pairs = {"|".join(sorted(f"{sample}:{sample}" for sample in pair)) for pair in pairs}
    observed_pairs = set()
    for row in alerts:
        if row["REASON_CODE"] == "potential_technical_duplicate":
            pair = "|".join(sorted(row["ENTITY_ID"].split("|")))
            if pair in observed_pairs:
                reasons.append("duplicate_pair_alert_repeated")
            observed_pairs.add(pair)
            if row["SCOPE"] != "PAIR" or row["STATUS"] != "ALERT" or row["METRIC"] != "PI_HAT" or not math.isclose(float(row["THRESHOLD"]), parameters["duplicate_pi_hat_alert"]) or not parameters["duplicate_pi_hat_alert"] <= float(row["OBSERVED_VALUE"]) <= 1:
                reasons.append("duplicate_alert_metric_invalid")
    rules["potential_technical_duplicate"] = dict(unit="pair", basis="injected_clone_not_estimator_concordance", **confusion(expected_pairs, observed_pairs, universe_pairs))
    if expected_pairs - observed_pairs:
        reasons.append("injected_clone_not_detected")
    report = _json(artifacts["qc_preliminary_report"])
    if report["duplicate_scan_status"] != "EVALUATED" or report["hwe_filter_applied"] or report["sex_check_status"] != "NOT_EVALUATED":
        reasons.append("preliminary_policy_or_evaluability_invalid")
    # Les résumés de lot sont contrôlés indépendamment, y compris les maxima.
    batch_rows = []
    missing = matrix.alleles[:, :, 0] < 0
    for batch in sorted(set(batches)):
        selected = np.array(batches) == batch
        sample_fractions = missing[:, selected].mean(axis=0)
        variant_fractions = missing[:, selected].mean(axis=1)
        batch_rows.append(dict(ARRAY_BATCH=batch, SAMPLE_COUNT=int(selected.sum()), VARIANT_COUNT=len(matrix.variants),
                               MEAN_SAMPLE_F_MISS=float(sample_fractions.mean()), MAX_SAMPLE_F_MISS=float(sample_fractions.max()),
                               MEAN_VARIANT_F_MISS=float(variant_fractions.mean()), MAX_VARIANT_F_MISS=float(variant_fractions.max())))
    comparisons["batches"] = compare_metrics(batch_rows, raw_table(artifacts["qc_batch_summary"]), "ARRAY_BATCH")
    reasons += ["batch_summary_" + key for key, count in comparisons["batches"].items() if key != "compared" and count]
    return dict(status="FAIL" if reasons else "PASS", reasons=sorted(set(reasons)), comparisons=comparisons,
                genotype_metrics=genotypes["metrics"], rules=rules, sex_check="NOT_EVALUATED_BY_AUTOSOMAL_CONTRACT")


def _evaluate_10(local: GenotypeMatrix, preliminary: QcExpectation, memberships: dict[str, list[str]], batches_by_sample: dict[str, str],
                 parameters: dict[str, Any], target_id: str, artifacts: dict[str, Path]) -> dict[str, Any]:
    """Évalue chaque cohorte 10 contre ses entrées préspecifiées, sans prétendre à 09."""
    all_samples = raw_table(artifacts["sample_qc_final"])
    all_variants = raw_table(artifacts["variant_qc_final"])
    all_batches = raw_table(artifacts["batch_qc_final"])
    reasons = []
    results = {}
    for cohort in COHORTS:
        source = local if cohort == COHORTS[2] else preliminary.retained
        membership = set(memberships[cohort])
        selected = np.array([sample in membership for _, sample in source.samples])
        if {sample for (_, sample), keep in zip(source.samples, selected, strict=True) if keep} != membership:
            raise TechnicalEvaluationError("d20_prespecified_cohort_not_in_oracle_source")
        subset_matrix = subset(source, selected, np.ones(len(source.variants), dtype=bool))
        expected = final_expectation(subset_matrix, [batches_by_sample[sample] for _, sample in subset_matrix.samples], cohort, parameters, target_id)
        dataset_id = cohort if cohort == COHORTS[2] else cohort + "_qc"
        sample_rows, variant_rows, batch_rows = [[row for row in rows if row["DATASET_ID"] == dataset_id] for rows in (all_samples, all_variants, all_batches)]
        comparisons, errors = _metric_checks(expected, sample_rows, variant_rows, batch_rows)
        genotype_result = compare_genotypes(expected.retained, _dataset(artifacts, dataset_id), stage="10")
        errors += genotype_result["reasons"]
        reasons += errors
        results[cohort] = dict(status="FAIL" if errors else "PASS", comparisons=comparisons, genotype_metrics=genotype_result["metrics"],
                               rules=_rules(expected, sample_rows, variant_rows, batch_rows, final=True, cohort=cohort))
    valid_datasets = {cohort if cohort == COHORTS[2] else cohort + "_qc" for cohort in COHORTS}
    if any(row["DATASET_ID"] not in valid_datasets for row in all_samples + all_variants + all_batches):
        reasons.append("unexpected_final_dataset")
    if any(row["COHORT_ID"] != next(cohort for cohort in COHORTS if row["DATASET_ID"] == (cohort if cohort == COHORTS[2] else cohort + "_qc")) for row in all_samples + all_variants + all_batches if row["DATASET_ID"] in valid_datasets):
        reasons.append("final_cohort_identity_mismatch")
    return dict(status="FAIL" if reasons else "PASS", reasons=sorted(set(reasons)), cohorts=results,
                execution_scope="ISOLATED_FINAL_QC_NOT_STAGE_09")


def aggregate_rules(evaluations: list[dict[str, Any]], stage: str) -> dict[str, Any]:
    """Additionne seulement des unités identiques ; ne moyenne pas des taux."""
    counts: dict[str, Counter] = {}
    units = {}
    for evaluation in evaluations:
        groups = [("", evaluation["rules"])] if stage == "05" else [(cohort + "/", value["rules"]) for cohort, value in evaluation["cohorts"].items()]
        for prefix, rules in groups:
            for name, values in rules.items():
                key = prefix + name
                counts.setdefault(key, Counter()).update({metric: values[metric] for metric in ("tp", "fp", "fn", "tn")})
                units[key] = values["unit"]
    return {key: dict(unit=units[key], **values, sensitivity=values["tp"] / (values["tp"] + values["fn"]) if values["tp"] + values["fn"] else None,
                      specificity=values["tn"] / (values["tn"] + values["fp"]) if values["tn"] + values["fp"] else None,
                      coverage="POSITIVE_AND_NEGATIVE" if values["tp"] + values["fn"] and values["tn"] + values["fp"] else "MISSING_CLASS")
            for key, values in counts.items()}


def evaluate_qc(*, panel_root: Path, run_root: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Évalue les quatorze cas clos ; aucune absence changée en PASS artificiel."""
    destination = new_workspace_run(workspace_root, run_id)
    panel_root, run_root = panel_root.absolute(), run_root.absolute()
    if panel_root.parent != destination.parent or run_root.parent != destination.parent:
        raise TechnicalEvaluationError("d20_evaluation_roots_invalid")
    panel = checked_inventory(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    receipt = checked_inventory(run_root, "run_index.json", "closure_state", "CLOSED")
    panel_sha, receipt_sha = _sha(panel_root / "panel_manifest.json"), _sha(run_root / "run_index.json")
    if receipt["panel_manifest_sha256"] != panel_sha or [item["bundle_id"] for item in receipt["runs"]] != [item["bundle_id"] for item in panel["bundles"]]:
        raise TechnicalEvaluationError("d20_receipt_panel_binding_invalid")
    config = _json(panel_root / "private_truth/configuration.json")
    if _sha(panel_root / "private_truth/configuration.json") != panel["configuration_sha256"]:
        raise TechnicalEvaluationError("d20_configuration_binding_invalid")
    target = _json(panel_root / "private_truth/target_configuration.json")["target_metadata"]
    expected_cases = generate_cases(config)
    if len(expected_cases) != len(receipt["runs"]):
        raise TechnicalEvaluationError("d20_case_coverage_incomplete")
    evaluations: dict[str, list[dict[str, Any]]] = {"05": [], "10": []}
    tool_versions = set()
    for (case, expected_panel, pairs), entry in zip(expected_cases, receipt["runs"], strict=True):
        bundle_id = entry["bundle_id"]
        source = panel_root / "blinded_inputs/bundles" / bundle_id
        truth = _json(panel_root / "private_truth" / f"{bundle_id}.json")
        metadata = raw_table(source / "metadata/samples.tsv")
        with _text(source / "observed.vcf.gz") as handle:
            observed = read_panel(handle, phased=False)
        if truth["case"] != case or observed.samples != expected_panel.samples or observed.variants != expected_panel.variants or not np.array_equal(observed.alleles, expected_panel.alleles) or truth["duplicate_pairs"] != [list(pair) for pair in pairs]:
            raise TechnicalEvaluationError("d20_injection_replay_mismatch")
        matrix = observation_matrix(observed, metadata)
        batches = [row["ARRAY_BATCH"] for row in metadata]
        if batches != ["LOT_A", "LOT_B"] * 96:
            raise TechnicalEvaluationError("d20_batch_truth_mismatch")
        expected_preliminary = preliminary_expectation(matrix, batches, config["preliminary_parameters"])
        local = with_target(matrix, target, truth["target_genotypes"])
        pipeline = contained_file(run_root, entry["run_path"] + "/manifest.json").parent
        if _sha(pipeline / "manifest.json") != entry["manifest_sha256"] or _sha(pipeline / "stages/10_qc_final/stage_inputs.json") != entry["final_inputs_sha256"] or entry["preliminary_returncode"] != 0 or entry["final_returncode"] != 0 or entry["cohort_producer"] != PRODUCER:
            raise TechnicalEvaluationError("d20_run_receipt_invalid")
        final_inputs = _json(pipeline / "stages/10_qc_final/stage_inputs.json")
        final_audit = _json(pipeline / "stages/10_qc_final/audit.json")
        tool_versions.update(tool["version"] for tool in final_audit["tools"] if tool["version"] is not None)
        verify_final_input_binding(final_inputs, final_audit, config["final_parameters"])
        for item in final_inputs["artifacts"]:
            if item["artifact_id"].startswith("cohort") and item["producer_stage"] != PRODUCER:
                raise TechnicalEvaluationError("d20_fake_cohort_producer")
        for cohort in COHORTS:
            expected_keep = "".join(f"{sample}\t{sample}\n" for sample in truth["cohorts"][cohort])
            if (pipeline / "stages/D20_cohort_contract_fixture" / f"{cohort}.keep").read_text() != expected_keep:
                raise TechnicalEvaluationError("d20_cohort_fixture_changed")
        artifacts_05 = _paths(pipeline, "stages/05_qc_preliminary")
        result_05 = _evaluate_05(matrix, expected_preliminary, artifacts_05, truth["duplicate_pairs"], batches, config["preliminary_parameters"])
        # La conversion et l'ajout cible sont des prérequis, pas des sources de vérité.
        base_result = compare_genotypes(matrix, _dataset(_paths(pipeline, "stages/03_convert_acpa"), "genomewide_base"))
        target_result = compare_genotypes(local, _dataset(_paths(pipeline, "stages/04_prepare_target_variant_dataset"), "target_variant"))
        if base_result["status"] != "PASS" or target_result["status"] != "PASS":
            result_05["status"] = "FAIL"
            result_05["reasons"].append("upstream_genotype_concordance_failed")
        result_10 = _evaluate_10(local, expected_preliminary, truth["cohorts"], dict(zip((sample for _, sample in matrix.samples), batches, strict=True)),
                                 config["final_parameters"], target["project_variant_id"], _paths(pipeline, "stages/10_qc_final"))
        for stage, result in (("05", result_05), ("10", result_10)):
            evaluations[stage].append(dict(bundle_id=bundle_id, case=case, **result))
    documents = {}
    for stage, cases in evaluations.items():
        rules = aggregate_rules(cases, stage)
        missing_classes = [name for name, values in rules.items() if values["coverage"] != "POSITIVE_AND_NEGATIVE"]
        documents[stage] = dict(schema_version="1.0.0", method_id=METHOD_ID, stage_id=stage, report_id="CAL-REP-" + stage,
                                status="FAIL" if any(item["status"] == "FAIL" for item in cases) or missing_classes else "PASS",
                                evaluable=True, evaluations=cases, rules=rules, missing_classes=missing_classes,
                                performance_decision="REVIEW_REQUIRED", uncertainty="NOT_ESTIMATED_PAIRED_COMPONENT_FIXTURES",
                                calibrated=False, founder_effect_proven=False)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        for stage, document in documents.items():
            write_json(temporary / f"CAL-REP-{stage}.json", document)
        manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, evaluation_state="CLOSED", created_at=datetime.now(timezone.utc).isoformat(),
                        panel_manifest_sha256=panel_sha, run_index_sha256=receipt_sha,
                        configuration=config, plink_versions=sorted(tool_versions),
                        stages=[dict(stage_id=stage, status=document["status"], case_count=len(document["evaluations"]),
                                     fail_count=sum(row["status"] == "FAIL" for row in document["evaluations"])) for stage, document in documents.items()],
                        artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        checked_inventory(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        checked_inventory(run_root, "run_index.json", "closure_state", "CLOSED")
        if _sha(panel_root / "panel_manifest.json") != panel_sha or _sha(run_root / "run_index.json") != receipt_sha:
            raise TechnicalEvaluationError("d20_sources_changed_during_evaluation")
        write_json(temporary / "evaluation_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def main() -> int:
    """CLI d'évaluation en lecture seule des sources, avec publication d'agrégats."""
    parser = argparse.ArgumentParser(description="Évaluation QC D20 après clôture, sans outil externe")
    for name in ("panel-root", "run-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    manifest = evaluate_qc(**vars(parser.parse_args()))
    print(json.dumps(manifest["stages"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
