"""Rejeu diagnostique D19 strictement limité aux configurations D18 00–02.

Lancement séparé de l'évaluateur : aucune vérité ni correspondance de profils
n'est passée au processus pipeline. Le reçu est scellé après la fin de tous
les processus ; l'approbation humaine reste absente. Ce n'est pas un pilote
de qualification et aucun outil génétique externe n'est activé.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any

from effet_fondateur.calibration.evaluate_technical import validate_schema
from effet_fondateur.calibration.export_blinded_inputs import load_checked_export
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.contracts import load_pipeline_config


def run_technical_inputs(*, d18_root: Path, workspace_root: Path, run_id: str, timeout_seconds: int = 300) -> dict[str, Any]:
    """Exécute uniquement 00–02 dans un dossier neuf, avec timeout par lot."""
    workspace_root = workspace_root.absolute()
    parts = tuple(part.casefold() for part in workspace_root.parts)
    if workspace_root.name != "calibration_workspace" or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(path.is_symlink() for path in (workspace_root, *workspace_root.parents)):
        raise TechnicalEvaluationError("d19_workspace_invalid")
    if d18_root.absolute().parent != workspace_root or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id) or type(timeout_seconds) is not int or timeout_seconds < 1:
        raise TechnicalEvaluationError("d19_runner_arguments_invalid")
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise TechnicalEvaluationError("d19_run_exists")
    load_checked_export(d18_root)
    inventory = _json(d18_root / "blinded_inputs/inventory.json")
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        entries = []
        for bundle in inventory["bundles"]:
            bundle_id = bundle["bundle_id"]
            source_root = (d18_root / "blinded_inputs/bundles" / bundle_id).absolute()
            config = load_pipeline_config(source_root / "pipeline.yaml")
            if set(config["stages"]) != {"initialize_run", "validate_sources", "build_sample_registry"} or not all(stage["enabled"] for stage in config["stages"].values()) or config["stages"]["build_sample_registry"]["parameters"] or any(value is not None for value in config["tools"].values()):
                raise TechnicalEvaluationError("d19_runner_scope_forbidden")
            output = temporary / "pipeline_runs" / bundle_id
            output.mkdir(parents=True)
            command = [str(Path(sys.executable).parent / "effet-fondateur"), "run", "--config", "pipeline.yaml", "--runs-dir", str(output)]
            with (output / "runner.log").open("wb") as log:
                completed = subprocess.run(command, cwd=source_root, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_seconds, check=False)
            if completed.returncode != 0:
                raise TechnicalEvaluationError("d19_partial_pipeline_process_failed")
            manifests = list(output.glob("*/manifest.json"))
            if len(manifests) != 1:
                raise TechnicalEvaluationError("d19_partial_pipeline_manifest_invalid")
            path = manifests[0]
            manifest = _json(path)
            if {stage["stage_id"] for stage in manifest["stages"]} != {"00", "01", "02"} or any(stage["state"] != "SUCCEEDED" for stage in manifest["stages"]) or "sample_registry_approval" not in manifest["manual_decisions_required"]:
                raise TechnicalEvaluationError("d19_partial_pipeline_closure_invalid")
            relative_run = path.parent.relative_to(temporary)
            entries.append(dict(bundle_id=bundle_id, run_path=(Path(run_id) / relative_run).as_posix(),
                                manifest_sha256=_sha(path), environment_sha256=_sha(path.parent / "environment.json"),
                                target_truth_path=None, target_truth_sha256=None))
        index = dict(schema_version="1.0.0", closure_state="CLOSED",
                     d18_manifest_sha256=_sha(d18_root / "private_truth/d18/blinded_export_manifest.json"), runs=entries)
        validate_schema(index, "technical_run_index")
        load_checked_export(d18_root)
        (temporary / "run_index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(target_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return index


def main() -> int:
    """CLI diagnostique ; aucun identifiant individuel en sortie."""
    parser = argparse.ArgumentParser(description="Rejeu régional 00–02 D19 uniquement")
    for name in ("d18-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    index = run_technical_inputs(**vars(parser.parse_args()))
    print(json.dumps(dict(closed_runs=len(index["runs"]), stages=["00", "01", "02"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
