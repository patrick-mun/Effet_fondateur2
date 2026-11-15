from __future__ import annotations

from pathlib import Path

import pytest

from effet_fondateur.calibration.contracts import load_json
from effet_fondateur.calibration.r12c_release import (
    FROZEN_COMPONENTS,
    G5_RECEIPT_PATH,
    R12cReleaseError,
    build_release,
    g5_receipt_document,
    validate_release,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r12c_freezes_granular_checkpoint_ledger() -> None:
    assert Path(
        "src/effet_fondateur/calibration/r13_checkpoint_ledger.py"
    ) in FROZEN_COMPONENTS
    assert Path(
        "src/effet_fondateur/calibration/r13_scientific_materializer.py"
    ) in FROZEN_COMPONENTS
    assert Path(
        "src/effet_fondateur/calibration/r13_scientific_call_worker.py"
    ) in FROZEN_COMPONENTS
    assert Path(
        "calibration/configurations/scientific_materialization_v1.yaml"
    ) in FROZEN_COMPONENTS


def test_r12c_builds_sealed_scientific_release_without_r13_results(
    tmp_path: Path,
) -> None:
    release_root = tmp_path / "r12c"

    audit = build_release(REPOSITORY_ROOT, release_root)
    receipt = g5_receipt_document(REPOSITORY_ROOT, release_root)

    assert audit["identity_count"] == 4_000
    assert audit["stage_count"] == 25
    assert audit["r13_result_count"] == 0
    assert audit["execution_workers"] == 8
    assert audit["materializable_call_count"] == 1_180
    assert receipt["gate_g5_passed"] is True
    assert receipt["receipt_revision"] == 2
    assert receipt["freeze_revision"] == "R12C-EXECUTABLE-REFREEZE-002"
    assert receipt["scientific_validation_authorized"] is True
    assert receipt["scientific_validation_executed"] is False
    assert receipt["scientific_volume_reduced"] is False
    assert receipt["materializer_state"] == "FULL_CALL_MATERIALIZATION_READY"
    assert receipt["single_call_worker_state"] == "FULL_SINGLE_CALL_EXECUTION_READY"
    execution_contract = load_json(
        release_root / "governance/r13_execution_contract.json"
    )
    assert execution_contract["call_count"] == 1_180
    assert execution_contract["truth_available_to_worker"] is False
    assert not (release_root / "runs").exists()
    assert not (release_root / "evaluations").exists()


def test_r12c_release_rejects_registry_tampering(tmp_path: Path) -> None:
    release_root = tmp_path / "r12c"
    build_release(REPOSITORY_ROOT, release_root)
    registry = release_root / "blinded_inputs/VALIDATION_BLIND/study_registry.json"
    registry.write_text("{}\n", encoding="utf-8")

    with pytest.raises(R12cReleaseError, match="checksum_mismatch"):
        validate_release(REPOSITORY_ROOT, release_root)


def test_tracked_r12c_g5_receipt_matches_frozen_release() -> None:
    release_root = (
        REPOSITORY_ROOT
        / "calibration_workspace/releases/panel_master_generic_v1_r12c_scientific_validation"
    )

    assert load_json(REPOSITORY_ROOT / G5_RECEIPT_PATH) == g5_receipt_document(
        REPOSITORY_ROOT,
        release_root,
    )
