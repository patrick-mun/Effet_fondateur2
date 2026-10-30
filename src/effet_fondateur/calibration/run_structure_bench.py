"""D22 : 00–03 réelles puis 08/09 isolées, amont et revues de fixture déclarés."""

from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

from effet_fondateur.calibration.evaluate_technical import raw_table
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.run_qc_bench import checked_stage
from effet_fondateur.calibration.structure_bench import CASES, METHOD_ID, PRODUCER, UPSTREAM_STAGES, checked_document
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import build_file_artifact, load_pipeline_config, validate_json_document


def fixture_artifacts(run_root: Path, source: Path) -> list[dict[str, Any]]:
    """Copie le triplet 03 et les contrats publics ; aucun 04–07 n'est prétendu."""
    fixture = run_root / "stages/D22_upstream_contract_fixture"
    shutil.copytree(source / "contract_fixtures", fixture)
    converted = checked_stage(run_root, "stages/03_convert_acpa")
    descriptor = deepcopy(_json(run_root / converted["genomewide_base_dataset"]["path"]))
    for suffix in ("bed", "bim", "fam"):
        path = fixture / ("kinship_panel." + suffix)
        shutil.copyfile(run_root / converted["genomewide_base_" + suffix]["path"], path)
        descriptor["files"][suffix] = dict(name=path.name, sha256=_sha(path))
    descriptor.update(dataset_id="d22_fixture_kinship_panel", source_format="PLINK_KINSHIP_PANEL")
    write_json(fixture / "kinship_panel.dataset.json", descriptor)
    signature = hashlib.sha256(json.dumps(file_inventory(fixture), sort_keys=True).encode()).hexdigest()
    definitions = [("kinship_panel_" + suffix, "kinship_panel." + suffix, None) for suffix in ("bed", "bim", "fam")]
    definitions += [("kinship_panel_dataset", "kinship_panel.dataset.json", "plink_dataset.schema.json"),
                    ("config_input_target_variant_metadata", "target_metadata.yaml", None)]
    definitions += [(name, name + ".tsv", name + ".schema.json") for name in ("qc_individual_metrics", "target_genotype_audit", "independent_set_proposal")]
    result = [checked_stage(run_root, "stages/02_build_sample_registry")["samples_master"]]
    for artifact_id, name, schema in definitions:
        path = fixture / name
        result.append(build_file_artifact(path, path.relative_to(run_root).as_posix(), artifact_id, artifact_id,
                                         "application/json" if path.suffix == ".json" else "text/plain", PRODUCER, signature,
                                         schema_name=schema, schema_version="1.0.0" if schema else None,
                                         assembly="GRCh38", sample_set_id=descriptor["sample_set_id"], variant_set_id=descriptor["variant_set_id"]))
    return result


def component_inputs(run_root: Path, stage: str, name: str, directory: str, parameters: dict[str, Any], artifacts: list[dict[str, Any]]) -> Path:
    """Enveloppe de composant avec signature des vraies entrées et paramètres."""
    output = run_root / "stages" / directory
    output.mkdir()
    signature = hashlib.sha256(json.dumps(dict(parameters=parameters, artifacts=artifacts), sort_keys=True).encode()).hexdigest()
    inputs = dict(schema_version="1.0.0", run_id=_json(run_root / "manifest.json")["run_id"], stage_id=stage, stage_name=name,
                  signature=signature, attempt_number=1, published_output_dir=output.relative_to(run_root).as_posix(), parameters=parameters, artifacts=artifacts)
    validate_json_document(inputs, "stage_inputs.schema.json")
    write_json(output / "stage_inputs.json", inputs)
    return output


def fixture_reviews(run_root: Path, upstream: list[dict[str, Any]], population: dict[str, dict[str, Any]], policy: dict[str, Any]) -> dict[str, Any]:
    """Lie les revues techniques autorisées, sans ouvrir une vérité biologique."""
    if policy["scope"] != "SYNTHETIC_D22_ONLY_NOT_SCIENTIFIC_APPROVAL" or policy["reviewer_role"] != "calibration_contract_fixture":
        raise TechnicalEvaluationError("d22_fixture_review_scope_invalid")
    proposal = next(item for item in upstream if item["artifact_id"] == "independent_set_proposal")
    reviews = {}
    for decision, item, field in (("kinship_exclusion_approval", proposal, "PROPOSED_INCLUDED"),
                                  ("population_structure_exclusion_approval", population["population_outliers"], "PROPOSED_EXCLUSION")):
        rows = raw_table(run_root / item["path"])
        excluded = sorted(row["SAMPLE_ID"] for row in rows if row[field] == ("false" if field == "PROPOSED_INCLUDED" else "true"))
        if excluded:
            reviews[decision] = dict(review_id="D22_fixture_" + decision, reviewer_role=policy["reviewer_role"], reviewed_at=policy["reviewed_at"],
                                     reviewed_artifact_sha256=item["sha256"], approved_exclusion_sample_ids=excluded)
    return reviews


def _execute_component(source: Path, output: Path, module: str, timeout_seconds: int) -> tuple[int, str]:
    """Exécute le point d'entrée réel ; le log reste dans le reçu sensible local."""
    with (output / "runner.log").open("wb") as log:
        completed = subprocess.run([sys.executable, "-m", "effet_fondateur.stages." + module, "--stage-inputs", str(output / "stage_inputs.json"), "--output-dir", str(output)],
                                   cwd=source, stdout=log, stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
    return completed.returncode, (output / "runner.log").read_text(encoding="utf-8")


def run_structure_bench(*, panel_root: Path, workspace_root: Path, run_id: str, timeout_seconds: int = 600) -> dict[str, Any]:
    """Clôture quatre cas, deux gels par cas et quatre refus de revue en base."""
    destination = new_workspace_run(workspace_root, run_id)
    panel_root = panel_root.absolute()
    if panel_root.parent != destination.parent or timeout_seconds < 1 or destination.with_name(run_id + ".failed").exists():
        raise TechnicalEvaluationError("d22_runner_arguments_invalid")
    panel = checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    panel_sha = _sha(panel_root / "panel_manifest.json")
    execution = _json(panel_root / "blinded_inputs/execution.json")
    policy = execution["fixture_review_policy"]
    if execution["producer"] != PRODUCER or execution["scope"] != "ISOLATED_COMPONENTS_08_09_NOT_04_07" or policy["scope"] != "SYNTHETIC_D22_ONLY_NOT_SCIENTIFIC_APPROVAL" or policy["reviewer_role"] != "calibration_contract_fixture" or [row["bundle_id"] for row in panel["bundles"]] != [f"B{index:02d}" for index in range(len(CASES))]:
        raise TechnicalEvaluationError("d22_execution_scope_invalid")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        receipts = []
        for index, entry in enumerate(panel["bundles"]):
            source = panel_root / "blinded_inputs/bundles" / entry["bundle_id"]
            config = load_pipeline_config(source / "pipeline.yaml")
            approval = config["stages"]["build_sample_registry"]["parameters"]["manual_approval"]
            if set(config["stages"]) != set(UPSTREAM_STAGES) or not all(value["enabled"] for value in config["stages"].values()) or any(value is not None for key, value in config["tools"].items() if key != "plink") or approval["decision_id"] != "CAL-DEC-006" or not approval["approved"]:
                raise TechnicalEvaluationError("d22_upstream_scope_forbidden")
            output = temporary / "pipeline_runs" / entry["bundle_id"]
            output.mkdir(parents=True)
            with (output / "runner.log").open("wb") as log:
                process = subprocess.run([str(Path(sys.executable).parent / "effet-fondateur"), "run", "--config", "pipeline.yaml", "--runs-dir", str(output)],
                                         cwd=source, stdout=log, stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
            manifests = list(output.glob("*/manifest.json"))
            if process.returncode or len(manifests) != 1:
                raise TechnicalEvaluationError("d22_upstream_process_failed")
            manifest_path = manifests[0]
            manifest = _json(manifest_path)
            if [row["stage_id"] for row in manifest["stages"]] != ["00", "01", "02", "03"] or any(row["state"] != "SUCCEEDED" for row in manifest["stages"]) or manifest["manual_decisions_required"]:
                raise TechnicalEvaluationError("d22_upstream_closure_invalid")
            run_root = manifest_path.parent
            upstream = fixture_artifacts(run_root, source)
            population_artifacts = [item for item in upstream if item["artifact_id"] in {"samples_master", "independent_set_proposal"} or item["artifact_id"].startswith("kinship_panel_")]
            directory = component_inputs(run_root, "08", "analyze_population_structure", "08_analyze_population_structure", execution["population_parameters"], population_artifacts)
            code, _ = _execute_component(source, directory, "analyze_population_structure", timeout_seconds)
            if code:
                raise TechnicalEvaluationError("d22_population_process_failed")
            population = checked_stage(run_root, "stages/08_analyze_population_structure")
            artifacts = [item for item in upstream if not item["artifact_id"].startswith("kinship_panel_")]
            artifacts += [population[key] for key in ("population_scores", "population_outliers")]
            parameters = deepcopy(execution["freeze_parameters"])
            parameters["decision_reviews"] = fixture_reviews(run_root, upstream, population, execution["fixture_review_policy"])
            repeated = []
            for repeat in ("a", "b"):
                directory = component_inputs(run_root, "09", "freeze_cohorts", "09_freeze_cohorts_" + repeat, parameters, artifacts)
                code, _ = _execute_component(source, directory, "freeze_cohorts", timeout_seconds)
                if code:
                    raise TechnicalEvaluationError("d22_freeze_process_failed")
                checked_stage(run_root, directory.relative_to(run_root).as_posix())
                repeated.append(directory.relative_to(run_root).as_posix())
            rejections = []
            if index == 0:
                definitions = (("MISSING_REVIEWS", 4, "manual_decision_unresolved:kinship_exclusion_approval"),
                               ("WRONG_SHA", 4, "reviewed_artifact_mismatch:kinship_exclusion_approval"),
                               ("INCOMPLETE_LIST", 4, "incomplete_approved_exclusions:kinship_exclusion_approval"),
                               ("NO_TIMEZONE", 2, "reviewed_at_requires_timezone:kinship_exclusion_approval"))
                for name, expected_code, expected_reason in definitions:
                    negative = deepcopy(parameters)
                    review = negative["decision_reviews"]["kinship_exclusion_approval"]
                    if name == "MISSING_REVIEWS":
                        negative["decision_reviews"] = {}
                    elif name == "WRONG_SHA":
                        review["reviewed_artifact_sha256"] = "0" * 64
                    elif name == "INCOMPLETE_LIST":
                        review["approved_exclusion_sample_ids"] = review["approved_exclusion_sample_ids"][:-1]
                    else:
                        review["reviewed_at"] = "2026-10-30T07:00:00"
                    directory = component_inputs(run_root, "09", "freeze_cohorts", "09_reject_" + name.lower(), negative, artifacts)
                    code, log_text = _execute_component(source, directory, "freeze_cohorts", timeout_seconds)
                    if code != expected_code or expected_reason not in log_text or any((directory / name).exists() for name in ("stage_outputs.json", "cohorts.frozen.tsv", "cohort_decisions.tsv", "keep")):
                        raise TechnicalEvaluationError("d22_expected_review_rejection_failed")
                    rejections.append(dict(case=name, directory=directory.relative_to(run_root).as_posix(), returncode=code, reason=expected_reason))
            receipts.append(dict(bundle_id=entry["bundle_id"], run_path=run_root.relative_to(temporary).as_posix(), manifest_sha256=_sha(manifest_path),
                                 upstream_returncode=process.returncode, population_returncode=0, freeze_directories=repeated, review_rejections=rejections,
                                 scope=execution["scope"], fixture_producer=PRODUCER))
        checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        if panel_sha != _sha(panel_root / "panel_manifest.json"):
            raise TechnicalEvaluationError("d22_sources_changed_during_run")
        receipt = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", panel_manifest_sha256=panel_sha,
                       runs=receipts, artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        write_json(temporary / "run_index.json", receipt)
        temporary.rename(destination)
    except Exception:
        temporary.rename(destination.with_name(run_id + ".failed"))
        raise
    return receipt


def main() -> int:
    """CLI de l'intégration bornée D22, jamais du pipeline complet."""
    parser = argparse.ArgumentParser(description="D22 : 00–03 puis composants 08/09 isolés")
    for name in ("panel-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    receipt = run_structure_bench(**vars(parser.parse_args()))
    print(dict(closed_runs=len(receipt["runs"]), state=receipt["closure_state"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
