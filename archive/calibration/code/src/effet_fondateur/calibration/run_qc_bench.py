"""Exécute seulement 00–05 puis le composant 10 avec cohortes de fixture.

Le processus scientifique ne reçoit aucun chemin privé. Aucun faux succès
09/10 n'est ajouté au manifest d'orchestration 00–05. Le reçu spécifique D20
empreinte l'appel séparé de 10, sa vraie provenance et ses sorties terminales.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import METHOD_ID, PRODUCER, COHORTS, checked_inventory, file_inventory, write_json
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import build_file_artifact, load_pipeline_config, validate_json_document


STAGES = ("initialize_run", "validate_sources", "build_sample_registry", "convert_acpa", "prepare_target_variant_dataset", "qc_preliminary")


def checked_stage(run_root: Path, directory: str) -> dict[str, dict[str, Any]]:
    """Relie audit, sorties, schémas et empreintes publiés par le vrai module."""
    from effet_fondateur.calibration.evaluate_technical import contained_file
    outputs = _json(contained_file(run_root, directory + "/stage_outputs.json"))
    audit = _json(contained_file(run_root, directory + "/audit.json"))
    validate_json_document(outputs, "stage_outputs.schema.json")
    validate_json_document(audit, "stage_audit.schema.json")
    if audit["outputs"] != outputs["artifacts"] or any(audit[key] != outputs[key] for key in ("run_id", "stage_id", "stage_name", "signature")):
        raise TechnicalEvaluationError("d20_stage_audit_binding_invalid")
    artifacts = {}
    for item in outputs["artifacts"]:
        path = contained_file(run_root, item["path"])
        if item["artifact_id"] in artifacts or not item["path"].startswith(directory + "/") or _sha(path) != item["sha256"] or item["producer_signature"] != outputs["signature"] or item["producer_stage"] != f'{outputs["stage_id"]}_{outputs["stage_name"]}':
            raise TechnicalEvaluationError("d20_stage_artifact_invalid")
        artifacts[item["artifact_id"]] = item
    return artifacts


def _final_inputs(run_root: Path, fixture: Path, parameters: dict[str, Any]) -> Path:
    """Assemble le contrat 10 avec producteurs honnêtes, sans exécuter 09."""
    artifacts = []
    for directory, prefixes in (("02_build_sample_registry", ("samples_master",)),
                                ("04_prepare_target_variant_dataset", ("target_variant_",)),
                                ("05_qc_preliminary", ("genomewide_pre_qc_",))):
        for item in checked_stage(run_root, "stages/" + directory).values():
            if any(item["artifact_id"] == prefix or item["artifact_id"] in {prefix + suffix for suffix in ("bed", "bim", "fam", "dataset")} for prefix in prefixes):
                artifacts.append(item)
    original = _json(run_root / "stages/04_prepare_target_variant_dataset/stage_inputs.json")
    artifacts.extend(item for item in original["artifacts"] if item["artifact_id"] == "config_input_target_variant_metadata")
    fixture_dir = run_root / "stages/D20_cohort_contract_fixture"
    shutil.copytree(fixture, fixture_dir)
    signature = hashlib.sha256(json.dumps(file_inventory(fixture_dir), sort_keys=True).encode()).hexdigest()
    for artifact_id, name in [("cohorts_frozen", "cohorts.tsv"), *[("cohort_keep_" + cohort, cohort + ".keep") for cohort in COHORTS]]:
        path = fixture_dir / name
        artifacts.append(build_file_artifact(path, path.relative_to(run_root).as_posix(), artifact_id, artifact_id,
                                             "text/tab-separated-values", PRODUCER, signature, assembly="GRCh38"))
    stage_signature = hashlib.sha256(json.dumps(dict(artifacts=artifacts, parameters=parameters), sort_keys=True).encode()).hexdigest()
    inputs = dict(schema_version="1.0.0", run_id=_json(run_root / "manifest.json")["run_id"], stage_id="10", stage_name="qc_final",
                  signature=stage_signature, attempt_number=1, published_output_dir="stages/10_qc_final", parameters=parameters, artifacts=artifacts)
    validate_json_document(inputs, "stage_inputs.schema.json")
    output_dir = run_root / "stages/10_qc_final"
    output_dir.mkdir()
    write_json(output_dir / "stage_inputs.json", inputs)
    return output_dir


def run_qc_bench(*, panel_root: Path, workspace_root: Path, run_id: str, timeout_seconds: int = 600) -> dict[str, Any]:
    """Exécute le banc fermé ; garde les diagnostics d'échec sans reçu CLOSED."""
    destination = new_workspace_run(workspace_root, run_id)
    panel_root = panel_root.absolute()
    if panel_root.parent != destination.parent or type(timeout_seconds) is not int or timeout_seconds < 1 or destination.with_name(run_id + ".failed").exists():
        raise TechnicalEvaluationError("d20_runner_arguments_invalid")
    panel = checked_inventory(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    panel_sha = _sha(panel_root / "panel_manifest.json")
    execution = _json(panel_root / "blinded_inputs/qc_execution.json")
    if execution["cohort_producer"] != PRODUCER or execution["scope"] != "ISOLATED_FINAL_QC_NOT_STAGE_09":
        raise TechnicalEvaluationError("d20_fixture_scope_invalid")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        receipts = []
        for bundle in panel["bundles"]:
            bundle_id = bundle["bundle_id"]
            source = panel_root / "blinded_inputs/bundles" / bundle_id
            config = load_pipeline_config(source / "pipeline.yaml")
            if set(config["stages"]) != set(STAGES) or not all(stage["enabled"] for stage in config["stages"].values()) or any(
                value is not None for key, value in config["tools"].items() if key != "plink") or not config["stages"]["build_sample_registry"]["parameters"]["manual_approval"]["approved"] or config["stages"]["build_sample_registry"]["parameters"]["manual_approval"]["decision_id"] != "CAL-DEC-004":
                raise TechnicalEvaluationError("d20_runner_scope_forbidden")
            output = temporary / "pipeline_runs" / bundle_id
            output.mkdir(parents=True)
            with (output / "runner.log").open("wb") as log:
                completed = subprocess.run([str(Path(sys.executable).parent / "effet-fondateur"), "run", "--config", "pipeline.yaml", "--runs-dir", str(output)],
                                           cwd=source, stdout=log, stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
            manifests = list(output.glob("*/manifest.json"))
            if completed.returncode != 0 or len(manifests) != 1:
                raise TechnicalEvaluationError("d20_preliminary_process_failed")
            manifest_path = manifests[0]
            manifest = _json(manifest_path)
            validate_json_document(manifest, "run_manifest.schema.json")
            if [entry["stage_id"] for entry in manifest["stages"]] != [f"{index:02d}" for index in range(6)] or any(entry["state"] != "SUCCEEDED" for entry in manifest["stages"]) or manifest["manual_decisions_required"]:
                raise TechnicalEvaluationError("d20_preliminary_closure_invalid")
            run_root = manifest_path.parent
            output_dir = _final_inputs(run_root, panel_root / "blinded_inputs/cohort_fixtures" / bundle_id, execution["final_parameters"])
            with (output / "final_runner.log").open("wb") as log:
                final_process = subprocess.run([sys.executable, "-m", "effet_fondateur.stages.qc_final", "--stage-inputs", str(output_dir / "stage_inputs.json"), "--output-dir", str(output_dir)],
                                               cwd=source, stdout=log, stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
            if final_process.returncode != 0:
                raise TechnicalEvaluationError("d20_final_process_failed")
            checked_stage(run_root, "stages/10_qc_final")
            receipts.append(dict(bundle_id=bundle_id, run_path=run_root.relative_to(temporary).as_posix(),
                                 manifest_sha256=_sha(manifest_path), final_inputs_sha256=_sha(output_dir / "stage_inputs.json"),
                                 preliminary_returncode=completed.returncode, final_returncode=final_process.returncode,
                                 final_scope="ISOLATED_FINAL_QC_NOT_STAGE_09", cohort_producer=PRODUCER))
        checked_inventory(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        if _sha(panel_root / "panel_manifest.json") != panel_sha:
            raise TechnicalEvaluationError("d20_panel_changed")
        receipt = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", panel_manifest_sha256=panel_sha,
                       runs=receipts, artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        write_json(temporary / "run_index.json", receipt)
        temporary.rename(destination)
    except Exception:
        temporary.rename(destination.with_name(run_id + ".failed"))
        raise
    return receipt


def main() -> int:
    """CLI du rejeu borné avec reçu CLOSED, jamais du pipeline complet."""
    parser = argparse.ArgumentParser(description="D20 : 00–05 puis QC 10 isolé, pas de pipeline complet")
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    receipt = run_qc_bench(**vars(parser.parse_args()))
    print(json.dumps(dict(closed_runs=len(receipt["runs"]), state=receipt["closure_state"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
