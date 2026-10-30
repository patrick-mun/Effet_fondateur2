"""Rejeu diagnostique D19 borné 00–04, avec cible technique explicite.

Le runner ne crée aucune approbation et n'interprète aucune vérité cible.
Le processus pipeline reçoit seulement l'arbre aveugle. Le reçu CLOSED est
publié après terminaison ; une erreur attendue reste FAILED dans le pipeline
et sera évaluée séparément comme un contrôle de rejet, pas un succès d'ajout.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any

from effet_fondateur.calibration.evaluate_technical import validate_schema
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_scenarios import MANIFEST_PATH, SCENARIOS, load_checked_panel, new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document


def _closed_receipt(output: Path, completed: subprocess.CompletedProcess, *, rejection: bool, bundle_id: str, temporary: Path, run_id: str, panel_root: Path, workspace_root: Path) -> dict[str, Any]:
    """Ancre les artefacts terminaux sans juger la réussite du contrôle négatif."""
    paths = list(output.glob("*/manifest.json"))
    # La CLI transforme StageExecutionError en parser.error (code shell 2).
    # Le code métier V2=4 reste dans le manifest et est vérifié séparément.
    if len(paths) != 1 or completed.returncode != (2 if rejection else 0):
        raise TechnicalEvaluationError("d19_0304_process_or_manifest_invalid")
    path = paths[0]
    manifest = _json(path)
    validate_json_document(manifest, "run_manifest.schema.json")
    if [item["stage_id"] for item in manifest["stages"]] != ["00", "01", "02", "03", "04"] or any(
        item["state"] != ("FAILED" if rejection and item["stage_id"] == "04" else "SUCCEEDED") for item in manifest["stages"]
    ) or manifest["manual_decisions_required"]:
        raise TechnicalEvaluationError("d19_0304_closure_invalid")
    receipt = dict(bundle_id=bundle_id, run_path=(Path(run_id) / path.parent.relative_to(temporary)).as_posix(),
                   manifest_sha256=_sha(path), environment_sha256=_sha(path.parent / "environment.json"),
                   target_truth_path=(panel_root.relative_to(workspace_root) / "private_truth/d19_inputs" / f"{bundle_id}.target.json").as_posix(),
                   target_truth_sha256=_sha(panel_root / "private_truth/d19_inputs" / f"{bundle_id}.target.json"))
    if rejection:
        evidence = list((path.parent / "stages/attempts").glob("04_prepare_target_variant_dataset.*.failed/stderr.log"))
        if len(evidence) != 1:
            raise TechnicalEvaluationError("d19_rejection_evidence_missing")
        receipt["rejection_evidence"] = dict(path=evidence[0].relative_to(path.parent).as_posix(), sha256=_sha(evidence[0]))
    return receipt


def run_technical_scenarios(*, panel_root: Path, workspace_root: Path, run_id: str, timeout_seconds: int = 600) -> dict[str, Any]:
    """Lance séquentiellement les lots préparés, hors de l'arbre aveugle immuable."""
    target = new_workspace_run(workspace_root, run_id)
    workspace_root, panel_root = workspace_root.absolute(), panel_root.absolute()
    if (workspace_root / f"{run_id}.failed").exists():
        raise TechnicalEvaluationError("d19_run_failure_diagnostics_exist")
    if panel_root.parent != workspace_root or type(timeout_seconds) is not int or timeout_seconds < 1:
        raise TechnicalEvaluationError("d19_0304_runner_arguments_invalid")
    panel = load_checked_panel(panel_root)
    panel_sha = _sha(panel_root / MANIFEST_PATH)
    rejection = panel["scenario"] == SCENARIOS[1]
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        receipts = []
        for bundle in panel["inventory"]["bundles"]:
            bundle_id = bundle["bundle_id"]
            source = panel_root / "blinded_inputs/bundles" / bundle_id
            config = load_pipeline_config(source / "pipeline.yaml")
            # Défense en profondeur : une préparation future ne doit pas activer
            # silencieusement QC, parenté, phasage ou analyse scientifique.
            if set(config["stages"]) != {"initialize_run", "validate_sources", "build_sample_registry", "convert_acpa", "prepare_target_variant_dataset"} or not all(stage["enabled"] for stage in config["stages"].values()) or any(
                value is not None for name, value in config["tools"].items() if name != "plink"
            ) or not config["tools"]["plink"] or not config["stages"]["build_sample_registry"]["parameters"].get("manual_approval", {}).get("approved"):
                raise TechnicalEvaluationError("d19_0304_runner_scope_forbidden")
            output = temporary / "pipeline_runs" / bundle_id
            output.mkdir(parents=True)
            command = [str(Path(sys.executable).parent / "effet-fondateur"), "run", "--config", "pipeline.yaml", "--runs-dir", str(output)]
            with (output / "runner.log").open("wb") as log:
                completed = subprocess.run(command, cwd=source, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_seconds, check=False)
            receipts.append(_closed_receipt(output, completed, rejection=rejection, bundle_id=bundle_id,
                                             temporary=temporary, run_id=run_id, panel_root=panel_root, workspace_root=workspace_root))
        load_checked_panel(panel_root)
        if _sha(panel_root / MANIFEST_PATH) != panel_sha:
            raise TechnicalEvaluationError("d19_0304_panel_changed")
        index = dict(schema_version="1.0.0", closure_state="CLOSED", d18_manifest_sha256=panel["source_d18_manifest_sha256"],
                     technical_panel_manifest_sha256=panel_sha, runs=receipts)
        validate_schema(index, "technical_run_index")
        (temporary / "run_index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(target)
    except Exception:
        # Garder les diagnostics d'un échec inattendu sans publier un reçu CLOSED.
        failure = workspace_root / f"{run_id}.failed"
        if not failure.exists():
            temporary.rename(failure)
        raise
    return index


def main() -> int:
    parser = argparse.ArgumentParser(description="Rejeu D19 technique 00–04 uniquement")
    for name in ("panel-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=600)
    index = run_technical_scenarios(**vars(parser.parse_args()))
    print(json.dumps(dict(closed_runs=len(index["runs"]), stages=["00", "01", "02", "03", "04"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
