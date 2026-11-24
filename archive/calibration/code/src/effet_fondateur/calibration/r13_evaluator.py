"""Évaluateur figé de la validation aveugle technique R13.

La vérité R12b n'est routée qu'après validation du manifeste de run clos. Les
verdicts produits portent exclusivement sur l'intégrité technique ; aucune
métrique de performance scientifique n'est calculée.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys
from typing import Any

from effet_fondateur.calibration.contracts import (
    STAGE_IDS,
    file_sha256,
    load_json,
    validate_document,
)
from effet_fondateur.calibration.r12b_validation import RELEASE_ID
from effet_fondateur.calibration.r13_runner import (
    DEFAULT_RUN_ID,
    validate_closed_run,
)
from effet_fondateur.calibration.truth_router import (
    allowed_truth_for_stage,
    authorize_truth_access,
)


class R13EvaluatorError(RuntimeError):
    """Erreur bloquante de l'évaluation post-clôture R13."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13EvaluatorError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_metrics(path: Path, stage_id: str, output_count: int) -> None:
    path.write_text(
        "stage_id\tmetric_id\tvalue\tunit\tdecision_scope\t"
        "scientific_performance_evaluated\n"
        f"{stage_id}\tclosed_outputs_verified\t{output_count}\tcount\t"
        "TECHNICAL_INTEGRITY\tfalse\n",
        encoding="utf-8",
    )


def evaluate_closed_run(
    repository_root: Path,
    candidate_root: Path,
    workspace: Path,
    run_id: str = DEFAULT_RUN_ID,
) -> dict[str, Any]:
    """Ouvre la vérité après clôture et produit 25 décisions techniques."""
    repository_root = repository_root.resolve(strict=True)
    candidate_root = candidate_root.resolve(strict=True)
    workspace = workspace.resolve(strict=True)
    run_manifest = validate_closed_run(
        repository_root,
        candidate_root,
        workspace,
        run_id,
    )
    evaluation_root = workspace / "evaluations"
    if evaluation_root.exists():
        _fail("evaluation_already_exists", str(evaluation_root))
    staging = workspace / ".evaluations.building"
    if staging.exists():
        _fail("evaluation_staging_exists", str(staging))
    truth_manifest = load_json(candidate_root / "private_truth/master_truth_manifest.json")
    if truth_manifest.get("disclosure_state") != "SEALED_PRIVATE":
        _fail("truth_manifest_not_sealed")
    output_index: dict[str, list[dict[str, Any]]] = {stage: [] for stage in STAGE_IDS}
    for output in run_manifest["stage_outputs"]:
        output_index[output["stage_id"]].append(output)
    try:
        staging.mkdir()
        evaluations: list[dict[str, Any]] = []
        for stage_id in STAGE_IDS:
            truth_ids = sorted(allowed_truth_for_stage(stage_id))
            grant = authorize_truth_access(
                evaluation_id=f"EVAL-{stage_id}",
                stage_id=stage_id,
                requested_truth_ids=truth_ids,
                truth_manifest=truth_manifest,
                run_manifest=run_manifest,
            )
            for truth_id, relative, digest in grant.artifacts:
                truth_path = candidate_root / relative
                if not truth_path.is_file() or file_sha256(truth_path) != digest:
                    _fail("granted_truth_hash_mismatch", truth_id)
            stage_root = staging / stage_id
            stage_root.mkdir()
            metrics = stage_root / f"ETAPE_{stage_id}.metrics.tsv"
            _write_metrics(metrics, stage_id, len(output_index[stage_id]))
            evaluation = {
                "schema_version": "1.0.0",
                "evaluation_id": f"EVAL-{stage_id}",
                "stage_id": stage_id,
                "release_id": RELEASE_ID,
                "closed_run_id": run_id,
                "truth_grants": truth_ids,
                "metrics_path": f"evaluations/{stage_id}/{metrics.name}",
                "metrics_sha256": file_sha256(metrics),
                "status": "PASS",
                "decision_scope": "TECHNICAL_INTEGRITY",
                "scientific_performance_evaluated": False,
                "technical_checks": [
                    "RUN_CLOSED_BEFORE_TRUTH_OPENING",
                    "OUTPUTS_DECLARED_AND_HASHED",
                    "TRUTH_SCOPE_LEAST_PRIVILEGE",
                    "VALIDATION_BLIND_PARTITION",
                    "NO_REAL_DATA",
                ],
                "counts": {
                    "N_EXPECTED": len(output_index[stage_id]),
                    "N_ELIGIBLE": len(output_index[stage_id]),
                    "N_EVALUABLE": len(output_index[stage_id]),
                    "N_EXCLUDED": 0,
                },
                "criteria_version": "1.0.0",
                "closed": True,
            }
            validate_document(
                evaluation,
                repository_root / "calibration/schemas/stage_evaluation_manifest.schema.json",
            )
            _write_json(stage_root / f"ETAPE_{stage_id}.evaluation.json", evaluation)
            evaluations.append(evaluation)
        decision = {
            "schema_version": "1.0.0",
            "decision_id": "R13-TECHNICAL-VALIDATION-001",
            "release_id": RELEASE_ID,
            "closed_run_id": run_id,
            "scope": "TECHNICAL_INTEGRITY_ONLY",
            "status": "PASS",
            "stage_count": len(evaluations),
            "stage_statuses": {item["stage_id"]: item["status"] for item in evaluations},
            "scientific_performance_evaluated": False,
            "r07_status": "REVISION_REQUIRED",
            "qualification_flags": {
                "calibrated": False,
                "founder_effect_proven": False,
                "production_promotion": False,
            },
            "restrictions": [
                "NO_SCIENTIFIC_QUALIFICATION_CLAIM",
                "NO_FOUNDER_EFFECT_CLAIM",
                "NO_PRODUCTION_PROMOTION",
                "R07_REVISION_REQUIRED_BEFORE_FULL_QUALIFICATION",
            ],
        }
        _write_json(staging / "R13_TECHNICAL_DECISION.json", decision)
        staging.rename(evaluation_root)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    validate_evaluations(repository_root, workspace, run_id)
    return decision


def validate_evaluations(
    repository_root: Path,
    workspace: Path,
    run_id: str = DEFAULT_RUN_ID,
) -> dict[str, Any]:
    """Valide les 25 manifestes techniques et leur synthèse restrictive."""
    repository_root = repository_root.resolve(strict=True)
    workspace = workspace.resolve(strict=True)
    root = workspace / "evaluations"
    decision = load_json(root / "R13_TECHNICAL_DECISION.json")
    if decision.get("scope") != "TECHNICAL_INTEGRITY_ONLY":
        _fail("evaluation_scope_mismatch")
    if decision.get("scientific_performance_evaluated") is not False:
        _fail("scientific_performance_flag_forbidden")
    flags = decision.get("qualification_flags")
    if flags != {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }:
        _fail("qualification_flags_mismatch")
    for stage_id in STAGE_IDS:
        stage_root = root / stage_id
        evaluation = load_json(stage_root / f"ETAPE_{stage_id}.evaluation.json")
        validate_document(
            evaluation,
            repository_root / "calibration/schemas/stage_evaluation_manifest.schema.json",
        )
        metrics = stage_root / f"ETAPE_{stage_id}.metrics.tsv"
        if evaluation.get("closed_run_id") != run_id:
            _fail("evaluation_run_mismatch", stage_id)
        if evaluation.get("decision_scope") != "TECHNICAL_INTEGRITY":
            _fail("evaluation_stage_scope_mismatch", stage_id)
        if not metrics.is_file() or file_sha256(metrics) != evaluation.get("metrics_sha256"):
            _fail("evaluation_metrics_hash_mismatch", stage_id)
    return decision


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
        decision = validate_evaluations(
            arguments.repository_root,
            arguments.workspace,
            arguments.run_id,
        )
    else:
        decision = evaluate_closed_run(
            arguments.repository_root,
            arguments.candidate_root,
            arguments.workspace,
            arguments.run_id,
        )
    print(json.dumps(decision, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
