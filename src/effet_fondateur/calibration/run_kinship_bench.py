"""Rejeu D21 réel et borné 00–03/05–07 ; propositions jamais approuvées."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any

from effet_fondateur.calibration.kinship_bench import METHOD_ID, STAGES, STAGE_IDS, checked_document
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.run_qc_bench import checked_stage
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document


def run_kinship_bench(*, panel_root: Path, workspace_root: Path, run_id: str, timeout_seconds: int = 1800) -> dict[str, Any]:
    """Publie CLOSED après sept vraies étapes, en gardant la revue 07 en attente.

    Le runner n'ouvre aucune vérité de pedigree. Il transmet uniquement la
    configuration de l'arbre aveugle et conserve les diagnostics d'échec.
    """
    destination = new_workspace_run(workspace_root, run_id)
    panel_root = panel_root.absolute()
    if panel_root.parent != destination.parent or timeout_seconds < 1 or destination.with_name(run_id + ".failed").exists():
        raise TechnicalEvaluationError("d21_runner_arguments_invalid")
    panel = checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    source_sha = _sha(panel_root / "panel_manifest.json")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        receipts = []
        for entry in panel["bundles"]:
            bundle_id = entry["bundle_id"]
            if bundle_id not in {"B00", "B01", "B02"}:
                raise TechnicalEvaluationError("d21_bundle_invalid")
            bundle = panel_root / "blinded_inputs/bundles" / bundle_id
            config = load_pipeline_config(bundle / "pipeline.yaml")
            approval = config["stages"]["build_sample_registry"]["parameters"].get("manual_approval", {})
            if set(config["stages"]) != set(STAGES) or not all(value["enabled"] for value in config["stages"].values()) or any(
                value is not None for key, value in config["tools"].items() if key not in {"plink", "king"}) or approval.get("decision_id") != "CAL-DEC-005" or not approval.get("approved"):
                raise TechnicalEvaluationError("d21_runner_scope_forbidden")
            output = temporary / "pipeline_runs" / bundle_id
            output.mkdir(parents=True)
            with (output / "runner.log").open("wb") as log:
                completed = subprocess.run([str(Path(sys.executable).parent / "effet-fondateur"), "run", "--config", "pipeline.yaml", "--runs-dir", str(output)],
                                           cwd=bundle, stdout=log, stderr=subprocess.STDOUT, check=False, timeout=timeout_seconds)
            manifests = list(output.glob("*/manifest.json"))
            if completed.returncode or len(manifests) != 1:
                raise TechnicalEvaluationError("d21_pipeline_process_failed")
            path = manifests[0]
            manifest = _json(path)
            validate_json_document(manifest, "run_manifest.schema.json")
            if [stage["stage_id"] for stage in manifest["stages"]] != list(STAGE_IDS) or any(stage["state"] != "SUCCEEDED" for stage in manifest["stages"]) or manifest["manual_decisions_required"] != ["kinship_exclusion_approval"]:
                raise TechnicalEvaluationError("d21_closure_or_manual_review_invalid")
            for stage_id, stage_name in zip(STAGE_IDS, STAGES, strict=True):
                checked_stage(path.parent, f"stages/{stage_id}_{stage_name}")
            receipts.append(dict(bundle_id=bundle_id, run_path=path.parent.relative_to(temporary).as_posix(),
                                 manifest_sha256=_sha(path), returncode=completed.returncode,
                                 pending_manual_decisions=manifest["manual_decisions_required"]))
        checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        if _sha(panel_root / "panel_manifest.json") != source_sha:
            raise TechnicalEvaluationError("d21_source_changed")
        receipt = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", panel_manifest_sha256=source_sha,
                       runs=receipts, artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        write_json(temporary / "run_index.json", receipt)
        temporary.rename(destination)
    except Exception:
        temporary.rename(destination.with_name(run_id + ".failed"))
        raise
    return receipt


def main() -> int:
    """CLI du seul diagnostic D21 ; affiche des comptes, jamais des sujets."""
    parser = argparse.ArgumentParser(description="D21 : 00–03/05–07 réels, sans approbation d'exclusion")
    for name in ("panel-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=1800)
    receipt = run_kinship_bench(**vars(parser.parse_args()))
    print({"closed_runs": len(receipt["runs"]), "closure": receipt["closure_state"]})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
