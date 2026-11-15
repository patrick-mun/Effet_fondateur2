"""Exécuteur figé de la validation aveugle technique R13.

Le runner copie uniquement la liste blanche ``VALIDATION_BLIND`` dans un
workspace transactionnel dépourvu de ``private_truth`` et de ``governance``.
Il ferme ensuite les sorties des 25 étapes. Il n'ouvre jamais la vérité et ne
produit aucune décision scientifique.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
from typing import Any

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    canonical_sha256,
    file_sha256,
    load_json,
    validate_document,
)
from effet_fondateur.calibration.micro_stage_worker import execute_stage
from effet_fondateur.calibration.r12b_validation import (
    PARTITION_ID,
    RELEASE_ID,
    validate_frozen_candidate,
)
from effet_fondateur.calibration.release_guard import validate_blinded_inventory


DEFAULT_RUN_ID = "RUN-R13-VALIDATION-BLIND-001"
EXECUTION_MODE = "R13_TECHNICAL_BLIND_VALIDATION"


class R13RunnerError(RuntimeError):
    """Erreur bloquante avant la clôture d'un run R13."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13RunnerError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _copy_closed_blind_inputs(candidate_root: Path, workspace: Path) -> dict[str, Any]:
    source_root = candidate_root / f"blinded_inputs/{PARTITION_ID}"
    inventory = load_json(source_root / "blinded_input_inventory.json")
    destination_root = workspace / f"blinded_inputs/{PARTITION_ID}"
    destination_root.mkdir(parents=True)
    for artifact in inventory["artifacts"]:
        source = source_root / artifact["relative_path"]
        destination = destination_root / artifact["relative_path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    _write_json(destination_root / "blinded_input_inventory.json", inventory)
    validate_blinded_inventory(
        inventory,
        view_root=destination_root,
        repository_root=workspace,
        verify_files=True,
    )
    return inventory


def _close_run(
    repository_root: Path,
    candidate_root: Path,
    workspace: Path,
    run_id: str,
    inventory: dict[str, Any],
) -> dict[str, Any]:
    run_root = workspace / "runs" / run_id
    outputs: list[dict[str, Any]] = []
    for stage_id in STAGE_IDS:
        stage_manifest = load_json(run_root / stage_id / "stage_outputs.json")
        for artifact in stage_manifest["artifacts"]:
            path = run_root / artifact["relative_path"]
            if file_sha256(path) != artifact["sha256"]:
                _fail("stage_output_hash_mismatch", stage_id)
            outputs.append(
                {
                    "stage_id": stage_id,
                    "artifact_id": artifact["artifact_id"],
                    "relative_path": path.relative_to(workspace).as_posix(),
                    "sha256": artifact["sha256"],
                }
            )
    frozen_manifest = candidate_root / "governance/r12b_frozen_manifest.json"
    manifest: dict[str, Any] = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "release_id": RELEASE_ID,
        "partition_id": PARTITION_ID,
        "run_state": "CLOSED",
        "closed_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "manifest_sha256": "0" * 64,
        "consumed_inventory_sha256": canonical_sha256(inventory),
        "stage_outputs": outputs,
        "dependencies": {
            "scope": "TECHNICAL_INTEGRITY_ONLY",
            "r12b_frozen_manifest_sha256": file_sha256(frozen_manifest),
            "runner_sha256": file_sha256(
                repository_root / "src/effet_fondateur/calibration/r13_runner.py"
            ),
            "truth_opened": False,
            "scientific_tools_executed": False,
            "real_data_accessed": False,
        },
    }
    manifest["manifest_sha256"] = canonical_sha256(
        {**manifest, "manifest_sha256": None}
    )
    _write_json(run_root / "pipeline_run_manifest.json", manifest)
    validate_document(
        manifest,
        repository_root / "calibration/schemas/pipeline_run_manifest.schema.json",
    )
    _write_json(
        run_root / "run_closure_receipt.json",
        {
            "schema_version": "1.0.0",
            "run_id": run_id,
            "release_id": RELEASE_ID,
            "partition_id": PARTITION_ID,
            "run_state": "CLOSED",
            "manifest_sha256": manifest["manifest_sha256"],
            "stage_count": len(STAGE_IDS),
            "truth_opened": False,
            "decision_scope": "TECHNICAL_INTEGRITY",
            "scientific_performance_evaluated": False,
        },
    )
    return manifest


def run_blind_validation(
    repository_root: Path,
    candidate_root: Path,
    workspace: Path,
    run_id: str = DEFAULT_RUN_ID,
) -> dict[str, Any]:
    """Exécute les étapes techniques sur une copie aveugle puis clôt le run."""
    repository_root = repository_root.resolve(strict=True)
    candidate_root = candidate_root.resolve(strict=True)
    workspace = workspace.resolve(strict=False)
    validate_frozen_candidate(repository_root, candidate_root)
    if workspace.exists():
        _fail("workspace_already_exists", str(workspace))
    if not run_id or "/" in run_id or "\\" in run_id or run_id in {".", ".."}:
        _fail("run_id_invalid", run_id)
    staging = workspace.with_name(f".{workspace.name}.building")
    if staging.exists():
        _fail("staging_already_exists", str(staging))
    try:
        staging.mkdir(parents=True)
        inventory = _copy_closed_blind_inputs(candidate_root, staging)
        run_root = staging / "runs" / run_id
        run_root.mkdir(parents=True)
        input_relative = Path(f"blinded_inputs/{PARTITION_ID}")
        run_relative = Path("runs") / run_id
        for stage_id in STAGE_IDS:
            execute_stage(
                release_root=staging,
                input_root=input_relative,
                run_root=run_relative,
                stage_id=stage_id,
                release_id=RELEASE_ID,
                execution_mode=EXECUTION_MODE,
            )
        manifest = _close_run(
            repository_root,
            candidate_root,
            staging,
            run_id,
            inventory,
        )
        staging.rename(workspace)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    validate_closed_run(repository_root, candidate_root, workspace, run_id)
    return manifest


def validate_closed_run(
    repository_root: Path,
    candidate_root: Path,
    workspace: Path,
    run_id: str = DEFAULT_RUN_ID,
) -> dict[str, Any]:
    """Vérifie la clôture sans ouvrir la vérité privée."""
    repository_root = repository_root.resolve(strict=True)
    candidate_root = candidate_root.resolve(strict=True)
    workspace = workspace.resolve(strict=True)
    validate_frozen_candidate(repository_root, candidate_root)
    if any((workspace / name).exists() for name in ("private_truth", "governance")):
        _fail("private_tree_present_in_run_workspace")
    inventory = load_json(
        workspace / f"blinded_inputs/{PARTITION_ID}/blinded_input_inventory.json"
    )
    validate_blinded_inventory(
        inventory,
        view_root=workspace / f"blinded_inputs/{PARTITION_ID}",
        repository_root=workspace,
        verify_files=True,
    )
    manifest_path = workspace / "runs" / run_id / "pipeline_run_manifest.json"
    manifest = load_json(manifest_path)
    validate_document(
        manifest,
        repository_root / "calibration/schemas/pipeline_run_manifest.schema.json",
    )
    if manifest.get("release_id") != RELEASE_ID or manifest.get("partition_id") != PARTITION_ID:
        _fail("closed_run_identity_mismatch")
    if manifest.get("consumed_inventory_sha256") != canonical_sha256(inventory):
        _fail("closed_run_inventory_mismatch")
    expected_manifest_hash = canonical_sha256({**manifest, "manifest_sha256": None})
    if manifest.get("manifest_sha256") != expected_manifest_hash:
        _fail("closed_run_manifest_hash_mismatch")
    stages = {output["stage_id"] for output in manifest["stage_outputs"]}
    if stages != set(STAGE_IDS):
        _fail("closed_run_stage_set_mismatch")
    for output in manifest["stage_outputs"]:
        path = workspace / output["relative_path"]
        if not path.is_file() or file_sha256(path) != output["sha256"]:
            _fail("closed_run_output_hash_mismatch", output["artifact_id"])
    for path in workspace.rglob("*"):
        if path.is_symlink():
            _fail("run_workspace_symlink_forbidden", str(path))
        if path.is_file() and path.stat(follow_symlinks=False).st_nlink != 1:
            _fail("run_workspace_hardlink_forbidden", str(path))
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument("--validate-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.validate_only:
        manifest = validate_closed_run(
            arguments.repository_root,
            arguments.candidate_root,
            arguments.workspace,
            arguments.run_id,
        )
    else:
        manifest = run_blind_validation(
            arguments.repository_root,
            arguments.candidate_root,
            arguments.workspace,
            arguments.run_id,
        )
    print(json.dumps(manifest, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
