"""D22 : évaluation après clôture de 08/09, sans processus externe."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np

from effet_fondateur.calibration.cohort_oracle import cohort_expectation
from effet_fondateur.calibration.evaluate_qc import _paths
from effet_fondateur.calibration.evaluate_technical import _dataset, contained_file, raw_table
from effet_fondateur.calibration.materialize_ibd import _json, _sha, _text
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.structure_bench import CASES, COHORTS, METHOD_ID, PRODUCER, checked_document, fixture_tables, generate_cases, registry_rows, upstream_configuration
from effet_fondateur.calibration.structure_oracle import compare_tables, counted_dosages, fit_reference_pca, structure_metrics
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, TechnicalEvaluationError, compare_genotypes, observation_matrix
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document


def _closed_upstream(run_root: Path, entry: dict[str, Any]) -> Path:
    """Ancre chaque audit et sortie 00–03 à son manifest, sans étape fabriquée."""
    path = contained_file(run_root, entry["run_path"] + "/manifest.json")
    manifest = _json(path)
    validate_json_document(manifest, "run_manifest.schema.json")
    if _sha(path) != entry["manifest_sha256"] or [row["stage_id"] for row in manifest["stages"]] != ["00", "01", "02", "03"] or manifest["manual_decisions_required"]:
        raise TechnicalEvaluationError("d22_upstream_manifest_changed")
    for stage in manifest["stages"]:
        directory = "stages/" + stage["stage_id"] + "_" + stage["stage_name"]
        audit = contained_file(path.parent, directory + "/audit.json")
        outputs = contained_file(path.parent, directory + "/stage_outputs.json")
        if stage["state"] != "SUCCEEDED" or stage["audit_sha256"] != _sha(audit) or stage["stage_outputs_sha256"] != _sha(outputs) or stage["signature"] != _json(audit)["signature"]:
            raise TechnicalEvaluationError("d22_upstream_stage_binding_invalid")
        _paths(path.parent, directory)
    return path.parent


def _component_paths(run_root: Path, directory: str, stage: str, name: str) -> dict[str, Path]:
    """Vérifie l'enveloppe du composant, les vraies entrées et leur signature."""
    inputs = _json(contained_file(run_root, directory + "/stage_inputs.json"))
    validate_json_document(inputs, "stage_inputs.schema.json")
    signature = hashlib.sha256(json.dumps(dict(parameters=inputs["parameters"], artifacts=inputs["artifacts"]), sort_keys=True).encode()).hexdigest()
    if inputs["stage_id"] != stage or inputs["stage_name"] != name or inputs["published_output_dir"] != directory or inputs["signature"] != signature:
        raise TechnicalEvaluationError("d22_component_identity_invalid")
    for artifact in inputs["artifacts"]:
        if artifact["producer_stage"] not in {PRODUCER, "02_build_sample_registry", "08_analyze_population_structure"} or _sha(contained_file(run_root, artifact["path"])) != artifact["sha256"]:
            raise TechnicalEvaluationError("d22_component_input_changed")
    outputs = _json(contained_file(run_root, directory + "/stage_outputs.json"))
    if any(outputs[key] != inputs[key] for key in ("run_id", "stage_id", "stage_name", "signature")):
        raise TechnicalEvaluationError("d22_component_output_identity_changed")
    return _paths(run_root, directory)


def _population_expectation(matrix: GenotypeMatrix, bim_path: Path, registry: list[dict[str, str]], proposals: list[dict[str, str]],
                            sample_set_id: str, config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    """Tables de référence construites depuis les génotypes et les entrées 08."""
    parameters = config["population_parameters"]
    reference_ids = {row["SAMPLE_ID"] for row in proposals if row["PROPOSED_INCLUDED"] == "true"}
    reference = np.array([sample in reference_ids for _, sample in matrix.samples])
    model = fit_reference_pca(counted_dosages(matrix, bim_path), reference, parameters)
    metadata = {row["SAMPLE_ID"]: row for row in registry}
    scores, outliers, loadings, eigenvalues = [], [], [], []
    for index, (family, sample) in enumerate(matrix.samples):
        row = metadata[sample]
        is_reference, flagged = bool(reference[index]), bool(model["outliers"][index])
        context = dict(SAMPLE_ID=sample, SAMPLE_SET_ID=sample_set_id, REFERENCE_INCLUDED=str(is_reference).lower(), PROJECTED=str(not is_reference).lower())
        status, reason = ("OUTLIER", "MULTIVARIATE_PC_DISTANCE") if flagged else ("INLIER", "")
        coordinates = {f"PC{axis+1}": float(model["scores"][index, axis]) if axis < model["components"] else "" for axis in range(10)}
        scores.append(dict(**context, FID=family, IID=sample, FAMILY_ID=family, GROUP_LABEL=row["GROUP_LABEL"], ARRAY_BATCH=row["ARRAY_BATCH"], **coordinates,
                           OUTLIER_STATUS=status, OUTLIER_REASON=reason))
        outliers.append(dict(**context, COMPONENTS_USED=2, DISTANCE_SQUARED=float(model["distances"][index]), P_VALUE=float(model["p_values"][index]),
                             OUTLIER_ALPHA=parameters["outlier_alpha"], OUTLIER_STATUS=status, OUTLIER_REASON=reason,
                             PROPOSED_EXCLUSION=str(flagged).lower(), MANUAL_VALIDATION_REQUIRED=str(flagged).lower()))
    informative_position = {int(index): rank for rank, index in enumerate(np.flatnonzero(model["included"]))}
    bim = [line.split() for line in bim_path.read_text().splitlines() if line.strip()]
    for index, (chromosome, position, variant) in enumerate(matrix.variants):
        included = bool(model["included"][index])
        code = "" if included else "LOW_REFERENCE_CALL_RATE" if not model["call_rate_eligible"][index] else "MONOMORPHIC_IN_REFERENCE"
        weights = {f"PC{axis+1}": float(model["weights"][informative_position[index], axis]) if included and axis < model["components"] else "" for axis in range(10)}
        loadings.append(dict(VARIANT_ID=variant, CHROMOSOME=chromosome, POSITION_BP=position, COUNTED_ALLELE=bim[index][4],
                             REFERENCE_ALLELE_FREQUENCY=float(model["frequencies"][index]), REFERENCE_STANDARD_DEVIATION=float(model["scales"][index]),
                             REFERENCE_CALL_RATE=float(model["call_rates"][index]), MODEL_STATUS="INCLUDED" if included else "EXCLUDED", EXCLUSION_CODE=code, **weights))
    for axis in range(model["components"]):
        eigenvalues.append(dict(COMPONENT=f"PC{axis+1}", EIGENVALUE=float(model["eigenvalues"][axis]),
                                EXPLAINED_VARIANCE_RATIO=float(model["explained_ratios"][axis]), CUMULATIVE_EXPLAINED_VARIANCE_RATIO=float(model["explained_ratios"][:axis+1].sum()),
                                SINGULAR_VALUE=float(model["singular_values"][axis]), REFERENCE_SAMPLE_COUNT=int(reference.sum()), INFORMATIVE_VARIANT_COUNT=model["informative_count"]))
    return model, dict(population_scores=scores, population_outliers=outliers, population_variant_loadings=loadings, population_eigenvalues=eigenvalues)


def evaluate_population(matrix: GenotypeMatrix, artifacts: dict[str, Path], registry: list[dict[str, str]], proposals: list[dict[str, str]],
                        descriptor: dict[str, Any], bim_path: Path, truth: dict[str, Any], config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Conformité exhaustive 08 et performance distincte contre les injections."""
    model, expected = _population_expectation(matrix, bim_path, registry, proposals, descriptor["sample_set_id"], config)
    keys = dict(population_scores="SAMPLE_ID", population_outliers="SAMPLE_ID", population_variant_loadings="VARIANT_ID", population_eigenvalues="COMPONENT")
    comparisons = {name: compare_tables(rows, raw_table(artifacts[name]), keys[name], config) for name, rows in expected.items()}
    reasons = [name + "_" + key for name, comparison in comparisons.items() for key, count in comparison.items() if key != "compared" and count]
    report = _json(artifacts["population_structure_report"])
    counts = dict(panel_sample_count=len(matrix.samples), reference_sample_count=truth["reference_count"], projected_sample_count=len(matrix.samples) - truth["reference_count"],
                  panel_variant_count=len(matrix.variants), informative_variant_count=model["informative_count"], component_count=model["components"], outlier_component_count=2,
                  outlier_count=int(model["outliers"].sum()), manual_validation_required=bool(model["outliers"].any()))
    if any(report.get(key) != value for key, value in counts.items()):
        reasons.append("population_report_count_mismatch")
    actual = raw_table(artifacts["population_outliers"])
    flags = np.array([row["PROPOSED_EXCLUSION"] == "true" for row in actual])
    if [row["SAMPLE_ID"] for row in actual] != [sample for _, sample in matrix.samples]:
        raise TechnicalEvaluationError("d22_population_output_axis_mismatch")
    published_scores = np.array([[float(row[f"PC{axis+1}"]) for axis in range(model["components"])] for row in raw_table(artifacts["population_scores"])])
    measured = structure_metrics(published_scores, truth, flags, config["confidence_level"])
    published_model = dict(reference=model["reference"], scores=published_scores,
                           eigenvalues=np.array([float(row["EIGENVALUE"]) for row in raw_table(artifacts["population_eigenvalues"])]))
    return dict(status="FAIL" if reasons else "PASS", reasons=reasons, comparisons=comparisons, metrics={**counts, **measured}, performance_decision="REVIEW_REQUIRED"), published_model


def evaluate_freeze(run_root: Path, directories: list[str], artifacts_08: dict[str, Path], registry: list[dict[str, str]],
                    tables: dict[str, list[dict[str, str]]], alt: str, config: dict[str, Any]) -> dict[str, Any]:
    """Oracle 09, comparaison des gels et de tous les fichiers keep."""
    first_inputs = _json(run_root / directories[0] / "stage_inputs.json")
    parameters = first_inputs["parameters"]
    if parameters["control_group_labels"] != ["CONTROL"] or parameters["required_nonempty_cohorts"] != list(COHORTS):
        raise TechnicalEvaluationError("d22_freeze_policy_changed")
    by_id = {item["artifact_id"]: item for item in first_inputs["artifacts"]}
    expected_reviews = {}
    for decision, artifact_id, rows, field, value in (
        ("kinship_exclusion_approval", "independent_set_proposal", tables["independent_set_proposal"], "PROPOSED_INCLUDED", "false"),
        ("population_structure_exclusion_approval", "population_outliers", raw_table(artifacts_08["population_outliers"]), "PROPOSED_EXCLUSION", "true")):
        excluded = sorted(row["SAMPLE_ID"] for row in rows if row[field] == value)
        if excluded:
            expected_reviews[decision] = dict(review_id="D22_fixture_" + decision, reviewer_role="calibration_contract_fixture", reviewed_at=config["fixture_review_policy"]["reviewed_at"],
                                              reviewed_artifact_sha256=by_id[artifact_id]["sha256"], approved_exclusion_sample_ids=excluded)
    if parameters["decision_reviews"] != expected_reviews:
        raise TechnicalEvaluationError("d22_fixture_review_binding_changed")
    expected = cohort_expectation(registry, tables["target_genotype_audit"], tables["qc_individual_metrics"], tables["independent_set_proposal"],
                                  raw_table(artifacts_08["population_scores"]), raw_table(artifacts_08["population_outliers"]), parameters, alt)
    results, repeat_files, keep_mismatch = [], [], 0
    for directory in directories:
        paths = _component_paths(run_root, directory, "09", "freeze_cohorts")
        inputs = _json(run_root / directory / "stage_inputs.json")
        if inputs["parameters"] != parameters or inputs["artifacts"] != first_inputs["artifacts"]:
            raise TechnicalEvaluationError("d22_repeat_scientific_inputs_differ")
        actual_members = raw_table(paths["cohorts_frozen"])
        actual_decisions = raw_table(paths["cohort_decisions"])
        wanted_members = [dict(ROW_KEY=row["COHORT_ID"] + "/" + row["SAMPLE_ID"], **row) for row in expected["rows"]]
        observed_members = [dict(ROW_KEY=row["COHORT_ID"] + "/" + row["SAMPLE_ID"], **row) for row in actual_members]
        wanted_decisions = [dict(ROW_KEY="/".join(row[key] for key in ("MANUAL_DECISION_ID", "SAMPLE_ID", "COHORT_ID")), **row) for row in expected["decisions"]]
        observed_decisions = [dict(ROW_KEY="/".join(row[key] for key in ("MANUAL_DECISION_ID", "SAMPLE_ID", "COHORT_ID")), **row) for row in actual_decisions]
        comparisons = dict(members=compare_tables(wanted_members, observed_members, "ROW_KEY", config),
                           summary=compare_tables(expected["summary"], raw_table(paths["cohort_summary"]), "COHORT_ID", config),
                           decisions=compare_tables(wanted_decisions, observed_decisions, "ROW_KEY", config))
        for cohort, text in expected["keeps"].items():
            keep_mismatch += paths["cohort_keep_" + cohort].read_text() != text
        files = {key: path.read_bytes() for key, path in paths.items()}
        repeat_files.append(files)
        results.append(comparisons)
    mismatches = sum(repeat_files[0].get(key) != repeat_files[1].get(key) for key in repeat_files[0].keys() | repeat_files[1].keys())
    reasons = ["repeat_" + str(index) + "_" + table + "_" + key for index, result in enumerate(results) for table, comparison in result.items() for key, count in comparison.items() if key != "compared" and count]
    if keep_mismatch:
        reasons.append("keep_files_mismatch")
    if mismatches:
        reasons.append("cohort_repeat_mismatch")
    report = _json(run_root / directories[0] / "cohort_freeze_report.json")
    expected_counts = {row["COHORT_ID"]: row["INCLUDED_COUNT"] for row in expected["summary"]}
    if report.get("cohort_included_counts") != expected_counts or report.get("manual_validation_required") is not False or report.get("sample_count") != len(registry):
        reasons.append("freeze_report_counts_or_review_mismatch")
    return dict(status="FAIL" if reasons else "PASS", reasons=reasons, comparisons=results,
                metrics=dict(cohort_repeat_mismatch_count=mismatches, keep_file_mismatch_count=keep_mismatch, membership_count=len(expected["rows"]),
                             decision_row_count=len(expected["decisions"]), included_counts=expected_counts, repeats=len(directories)),
                performance_decision="REVIEW_REQUIRED", real_kinship_qualification=False)


def _invariance(baseline: dict[str, Any], observed: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    """Compare les scores de référence après les seules ambiguïtés de signe."""
    reference = baseline["reference"]
    first, second = baseline["scores"][reference], observed["scores"][reference].copy()
    sign = np.where((first * second).sum(axis=0) < 0, -1, 1)
    second *= sign
    maximum = float(np.max(np.abs(first - second)))
    matched = np.allclose(first, second, rtol=config["relative_tolerance"], atol=config["absolute_tolerance"]) and np.allclose(baseline["eigenvalues"], observed["eigenvalues"], rtol=config["relative_tolerance"], atol=config["absolute_tolerance"])
    return dict(status="PASS" if matched else "FAIL", reference_score_max_absolute_difference_after_sign_alignment=maximum,
                axis_sign_alignment_only=True, arbitrary_rotation_allowed=False)


def evaluate_structure(*, panel_root: Path, run_root: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Publie les résultats agrégés après contrôle des sources et clôtures."""
    destination = new_workspace_run(workspace_root, run_id)
    panel_root, run_root = panel_root.absolute(), run_root.absolute()
    if panel_root.parent != destination.parent or run_root.parent != destination.parent:
        raise TechnicalEvaluationError("d22_evaluation_roots_invalid")
    panel = checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    receipt = checked_document(run_root, "run_index.json", "closure_state", "CLOSED")
    panel_sha, receipt_sha = _sha(panel_root / "panel_manifest.json"), _sha(run_root / "run_index.json")
    if receipt["panel_manifest_sha256"] != panel_sha or [entry["bundle_id"] for entry in receipt["runs"]] != [f"B{index:02d}" for index in range(len(CASES))]:
        raise TechnicalEvaluationError("d22_receipt_binding_invalid")
    config = _json(panel_root / "private_truth/configuration.json")
    if _sha(panel_root / "private_truth/configuration.json") != panel["configuration_sha256"]:
        raise TechnicalEvaluationError("d22_configuration_binding_invalid")
    panels, truth = generate_cases(config)
    if _json(panel_root / "private_truth/structure_truth.json") != truth:
        raise TechnicalEvaluationError("d22_private_truth_replay_failed")
    target_path = Path(__file__).resolve().parents[3] / "calibration/configurations" / config["target_configuration"]
    if _sha(target_path) != config["target_configuration_sha256"]:
        raise TechnicalEvaluationError("d22_public_annotation_source_changed")
    target = _json(target_path)
    results = {"08": [], "09": []}
    models, versions = [], set()
    for index, entry in enumerate(receipt["runs"]):
        source = panel_root / "blinded_inputs/bundles" / entry["bundle_id"]
        if load_pipeline_config(source / "pipeline.yaml") != upstream_configuration(config, target):
            raise TechnicalEvaluationError("d22_public_profile_changed")
        registry = raw_table(source / "metadata/samples.tsv")
        if registry != registry_rows(panels[index], CASES[index]):
            raise TechnicalEvaluationError("d22_public_registry_changed")
        with _text(source / "observed.vcf.gz") as handle:
            observed = read_panel(handle, phased=False)
        if observed.samples != panels[index].samples or observed.variants != panels[index].variants or not np.array_equal(observed.alleles, panels[index].alleles):
            raise TechnicalEvaluationError("d22_observation_replay_failed")
        pipeline = _closed_upstream(run_root, entry)
        conversion = _paths(pipeline, "stages/03_convert_acpa")
        matrix = _dataset(conversion, "genomewide_base")
        if compare_genotypes(observation_matrix(observed, registry), matrix)["status"] != "PASS":
            raise TechnicalEvaluationError("d22_conversion_not_concordant")
        tables = {name: raw_table(source / "contract_fixtures" / (name + ".tsv")) for name in ("target_genotype_audit", "qc_individual_metrics", "independent_set_proposal")}
        if tables != fixture_tables(panels[index], registry, config, target["target_metadata"]):
            raise TechnicalEvaluationError("d22_fixture_truth_mismatch")
        fixture = pipeline / "stages/D22_upstream_contract_fixture"
        for name in ("target_genotype_audit.tsv", "qc_individual_metrics.tsv", "independent_set_proposal.tsv", "target_metadata.yaml"):
            if _sha(fixture / name) != _sha(source / "contract_fixtures" / name):
                raise TechnicalEvaluationError("d22_copied_fixture_changed")
        for suffix in ("bed", "bim", "fam"):
            if _sha(fixture / ("kinship_panel." + suffix)) != _sha(conversion["genomewide_base_" + suffix]):
                raise TechnicalEvaluationError("d22_fixture_panel_not_exact_copy")
        descriptor = _json(fixture / "kinship_panel.dataset.json")
        population = _component_paths(pipeline, "stages/08_analyze_population_structure", "08", "analyze_population_structure")
        if _json(pipeline / "stages/08_analyze_population_structure/stage_inputs.json")["parameters"] != config["population_parameters"]:
            raise TechnicalEvaluationError("d22_population_parameters_changed")
        evaluation_08, model = evaluate_population(matrix, population, registry, tables["independent_set_proposal"], descriptor, fixture / "kinship_panel.bim", truth, config)
        models.append(model)
        evaluation_09 = evaluate_freeze(pipeline, entry["freeze_directories"], population, registry, tables, target["target_metadata"]["alt"], config)
        rejections = entry["review_rejections"]
        if index == 0 and [row["case"] for row in rejections] != ["MISSING_REVIEWS", "WRONG_SHA", "INCOMPLETE_LIST", "NO_TIMEZONE"]:
            raise TechnicalEvaluationError("d22_negative_case_coverage_missing")
        for rejection in rejections:
            directory = contained_file(pipeline, rejection["directory"] + "/runner.log").parent
            if rejection["reason"] not in (directory / "runner.log").read_text() or any((directory / name).exists() for name in ("stage_outputs.json", "cohorts.frozen.tsv", "keep")):
                raise TechnicalEvaluationError("d22_negative_review_evidence_invalid")
        evaluation_09["review_rejections"] = [dict(case=row["case"], returncode=row["returncode"], reason=row["reason"]) for row in rejections]
        for stage, evaluation in (("08", evaluation_08), ("09", evaluation_09)):
            results[stage].append(dict(bundle_id=entry["bundle_id"], case=CASES[index], **evaluation))
        audit = _json(pipeline / "stages/08_analyze_population_structure/audit.json")
        versions.update(f'{tool["tool"]}: {tool["version"]}' for tool in audit["tools"])
    for index, label in ((1, "projected_genotype_invariance"), (2, "descriptive_label_invariance")):
        check = _invariance(models[0], models[index], config)
        results["08"][index][label] = check
        if check["status"] != "PASS":
            results["08"][index]["reasons"].append(label + "_failed")
            results["08"][index]["status"] = "FAIL"
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        stages = []
        for stage, evaluations in results.items():
            failed = sum(row["status"] == "FAIL" for row in evaluations)
            document = dict(schema_version="1.0.0", method_id=METHOD_ID, stage_id=stage, status="FAIL" if failed else "PASS", evaluable=True, evaluations=evaluations,
                            performance_decision="REVIEW_REQUIRED", calibrated=False, founder_effect_proven=False)
            write_json(temporary / ("CAL-REP-" + stage + ".json"), document)
            stages.append(dict(stage_id=stage, status=document["status"], case_count=len(evaluations), fail_count=failed))
        manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, evaluation_state="CLOSED", created_at=datetime.now(timezone.utc).isoformat(),
                        configuration=config, tool_versions=sorted(versions), stages=stages, panel_manifest_sha256=panel_sha, run_index_sha256=receipt_sha,
                        scope="ISOLATED_COMPONENTS_08_09_NOT_04_07", kinship_07_issue="OPEN_DEFERRED_UNTIL_GLOBAL_REVIEW_BY_USER", artifacts=file_inventory(temporary),
                        calibrated=False, founder_effect_proven=False)
        checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        checked_document(run_root, "run_index.json", "closure_state", "CLOSED")
        if panel_sha != _sha(panel_root / "panel_manifest.json") or receipt_sha != _sha(run_root / "run_index.json"):
            raise TechnicalEvaluationError("d22_sources_changed_during_evaluation")
        write_json(temporary / "evaluation_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def main() -> int:
    """CLI d'évaluation uniquement, sans processus PLINK ou pipeline."""
    parser = argparse.ArgumentParser(description="Évaluateurs D22 après clôture")
    for name in ("panel-root", "run-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    print(evaluate_structure(**vars(parser.parse_args()))["stages"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
