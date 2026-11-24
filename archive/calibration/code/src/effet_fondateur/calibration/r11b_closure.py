"""Clôture suivie de la confirmation de ressources R11b.

Le workspace R11b est volontairement ignoré par Git. Ce module vérifie ses
reçus physiques, lie leur empreinte à R07 v2 et produit le reçu suivi qui
autorise seulement la construction d'une nouvelle release scientifique.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from effet_fondateur.calibration.contracts import (
    canonical_sha256,
    file_sha256,
    load_json,
)
from effet_fondateur.calibration.r07_revision import (
    CONFIGURATION_PATH as R07_CONFIGURATION_PATH,
    RECEIPT_PATH as R07_RECEIPT_PATH,
    validate_tracked_receipt,
)


WORKSPACE_RELATIVE = Path(
    "calibration_workspace/resource_probes/r11b_scientific_resources_v1"
)
DECISION_RELATIVE = WORKSPACE_RELATIVE / "evaluations/R11B_RESOURCE_DECISION.json"
RECEIPT_RELATIVE = Path("calibration/readiness/R11B_RESOURCE_CONFIRMATION.json")
EXPECTED_FLAGS = {
    "calibrated": False,
    "founder_effect_proven": False,
    "production_promotion": False,
}
SHAPEIT_RECEIPTS = {
    "dlow": WORKSPACE_RELATIVE
    / "probes/shapeit_dlow_attempt01/governance/measurement.json",
    "bridge": WORKSPACE_RELATIVE
    / "probes/shapeit_bridge_attempt02/governance/measurement.json",
    "target": WORKSPACE_RELATIVE
    / "probes/shapeit_target_attempt01/governance/measurement.json",
}


class R11bClosureError(RuntimeError):
    """Violation bloquante de la clôture R11b."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R11bClosureError(code, detail)


def _load_completed_receipt(path: Path, expected_scope: str) -> dict[str, Any]:
    if not path.is_file():
        _fail("r11b_measurement_missing", str(path))
    receipt = load_json(path)
    if receipt.get("state") != "COMPLETED":
        _fail("r11b_measurement_not_completed", str(path))
    if receipt.get("scope") != expected_scope:
        _fail("r11b_measurement_scope_mismatch", str(path))
    if receipt.get("qualification_flags") != EXPECTED_FLAGS:
        _fail("r11b_measurement_flags_mismatch", str(path))
    return receipt


def receipt_document(repository_root: Path) -> dict[str, Any]:
    """Valide les preuves R11b et retourne leur reçu de clôture déterministe."""
    repository_root = repository_root.resolve(strict=True)
    r07 = validate_tracked_receipt(repository_root)
    decision_path = repository_root / DECISION_RELATIVE
    if not decision_path.is_file():
        _fail("r11b_decision_missing")
    decision = load_json(decision_path)
    if decision.get("decision_id") != "R11B-SCIENTIFIC-RESOURCE-CONFIRMATION-V1":
        _fail("r11b_decision_identity_mismatch")
    if decision.get("status") != "RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE":
        _fail("r11b_resource_gate_not_ready")
    checks = decision.get("checks")
    if not isinstance(checks, list) or len(checks) != 5:
        _fail("r11b_check_count_mismatch")
    if any(check.get("status") != "PASS" for check in checks):
        _fail("r11b_resource_check_failed")
    if decision.get("qualification_flags") != EXPECTED_FLAGS:
        _fail("r11b_decision_flags_mismatch")
    if decision.get("scientific_g5_authorized") is not False:
        _fail("r11b_scientific_g5_flag_invalid")
    if decision.get("next_phase") != "BUILD_NEW_SCIENTIFIC_RELEASE":
        _fail("r11b_next_phase_mismatch")

    shapeit = {
        name: _load_completed_receipt(repository_root / relative, "RESOURCE_CONFIRMATION_ONLY")
        for name, relative in SHAPEIT_RECEIPTS.items()
    }
    expected_counts = {
        "dlow": (250, 250),
        "bridge": (2_000, 1_500),
        "target": (4_000, 3_000),
    }
    for name, receipt in shapeit.items():
        study, reference = expected_counts[name]
        counts = receipt.get("counts", {})
        if (
            counts.get("study_samples") != study
            or counts.get("reference_samples") != reference
        ):
            _fail("r11b_shapeit_count_mismatch", name)

    ibd_paths = sorted(
        (repository_root / WORKSPACE_RELATIVE / "probes").glob(
            "r11b_ibd_*/governance/measurement.json"
        )
    )
    if len(ibd_paths) != 18:
        _fail("r11b_ibd_receipt_count_mismatch", str(len(ibd_paths)))
    ibd_receipts = [
        _load_completed_receipt(path, "RESOURCE_CONFIRMATION_ONLY")
        for path in ibd_paths
    ]
    measurement_bindings = [
        {
            "relative_path": relative.as_posix(),
            "sha256": file_sha256(repository_root / relative),
        }
        for relative in SHAPEIT_RECEIPTS.values()
    ] + [
        {
            "relative_path": path.relative_to(repository_root).as_posix(),
            "sha256": file_sha256(path),
        }
        for path in ibd_paths
    ]
    if len(ibd_receipts) != 18:
        _fail("r11b_ibd_receipt_count_mismatch")

    return {
        "schema_version": "1.0.0",
        "receipt_id": "R11B-RESOURCE-CONFIRMATION-READY",
        "issued_on": str(decision["created_at"])[:10],
        "r07_revision_id": r07["revision_id"],
        "r07_status": r07["status"],
        "decision_id": decision["decision_id"],
        "status": decision["status"],
        "checks": checks,
        "shapeit_probe_count": len(shapeit),
        "ibd_probe_count": len(ibd_receipts),
        "measurement_binding_sha256": canonical_sha256(measurement_bindings),
        "evidence_sha256": {
            R07_CONFIGURATION_PATH.as_posix(): file_sha256(
                repository_root / R07_CONFIGURATION_PATH
            ),
            R07_RECEIPT_PATH.as_posix(): file_sha256(
                repository_root / R07_RECEIPT_PATH
            ),
            DECISION_RELATIVE.as_posix(): file_sha256(decision_path),
        },
        "new_scientific_release_construction_authorized": True,
        "scientific_g5_passed": False,
        "scientific_validation_executed": False,
        "qualification_flags": EXPECTED_FLAGS,
        "restrictions": [
            "NEW_INDEPENDENT_RELEASE_REQUIRED",
            "VALIDATION_BLIND_MUST_BE_MATERIALIZED_AND_SEALED",
            "SCIENTIFIC_RUNNER_AND_EVALUATOR_MUST_BE_FROZEN_BEFORE_G5",
            "NO_REAL_DATA_AUTHORIZED",
            "NO_SCIENTIFIC_QUALIFICATION_CLAIM",
        ],
        "next_phase": "R12C_BUILD_AND_FREEZE_NEW_SCIENTIFIC_RELEASE",
    }


def write_tracked_receipt(repository_root: Path) -> dict[str, Any]:
    """Écrit le reçu suivi après validation physique complète des preuves."""
    repository_root = repository_root.resolve(strict=True)
    document = receipt_document(repository_root)
    path = repository_root / RECEIPT_RELATIVE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return document


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--write", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    document = (
        write_tracked_receipt(arguments.repository_root)
        if arguments.write
        else receipt_document(arguments.repository_root)
    )
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
