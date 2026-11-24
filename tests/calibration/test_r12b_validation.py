from __future__ import annotations

import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.contracts import STAGE_IDS, file_sha256, load_json
from effet_fondateur.calibration.r12b_validation import (
    FROZEN_COMPONENTS,
    PARTITION_ID,
    RELEASE_ID,
    R12bValidationError,
    build_frozen_candidate,
    g5_receipt_document,
    validate_frozen_candidate,
)
from effet_fondateur.calibration.release_guard import ReleaseGuardError
from effet_fondateur.calibration.r13_evaluator import evaluate_closed_run
from effet_fondateur.calibration.r13_runner import run_blind_validation


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r12b_builds_sealed_validation_candidate_without_r13_outputs(
    tmp_path: Path,
) -> None:
    release_root = tmp_path / "r12b"

    audit = build_frozen_candidate(REPOSITORY_ROOT, release_root)
    receipt = g5_receipt_document(REPOSITORY_ROOT, release_root)

    assert audit.release_id == RELEASE_ID
    assert audit.partition_id == PARTITION_ID
    assert audit.truth_family_count == 14
    assert audit.view_count == 11
    assert audit.blind_artifact_count == 11
    assert audit.run_count == 0
    assert audit.evaluation_count == 0
    assert audit.r07_status == "REVISION_REQUIRED"
    assert receipt["gate_g5_passed"] is True
    assert receipt["scope"] == "R13_TECHNICAL_BLIND_VALIDATION"
    assert receipt["scientific_validation_authorized"] is False
    assert receipt["qualification_flags"] == {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }
    assert not (release_root / "runs").exists()
    assert not (release_root / "evaluations").exists()


def test_tracked_g5_receipt_matches_rebuilt_candidate(tmp_path: Path) -> None:
    release_root = tmp_path / "r12b"
    build_frozen_candidate(REPOSITORY_ROOT, release_root)

    assert load_json(
        REPOSITORY_ROOT / "calibration/readiness/R12B_G5_READY.json"
    ) == g5_receipt_document(REPOSITORY_ROOT, release_root)


def test_r12b_blind_documents_contain_only_pseudonyms(tmp_path: Path) -> None:
    release_root = tmp_path / "r12b"
    build_frozen_candidate(REPOSITORY_ROOT, release_root)

    inventory = load_json(
        release_root / "blinded_inputs/VALIDATION_BLIND/blinded_input_inventory.json"
    )
    for artifact in inventory["artifacts"]:
        document = load_json(
            release_root / "blinded_inputs/VALIDATION_BLIND" / artifact["relative_path"]
        )
        serialized = json.dumps(document, sort_keys=True)
        assert "VAL-R12B-" not in serialized
        assert "master_id" not in serialized
        assert all(subject["sample_id"].startswith("SYN-") for subject in document["subjects"])


def test_r12b_freezes_exact_r13_components(tmp_path: Path) -> None:
    release_root = tmp_path / "r12b"
    build_frozen_candidate(REPOSITORY_ROOT, release_root)

    frozen = load_json(release_root / "governance/r12b_frozen_manifest.json")
    bindings = {item["relative_path"]: item["sha256"] for item in frozen["frozen_components"]}

    assert set(bindings) == {path.as_posix() for path in FROZEN_COMPONENTS}
    assert bindings["src/effet_fondateur/calibration/r13_runner.py"] == file_sha256(
        REPOSITORY_ROOT / "src/effet_fondateur/calibration/r13_runner.py"
    )
    assert bindings["src/effet_fondateur/calibration/r13_evaluator.py"] == file_sha256(
        REPOSITORY_ROOT / "src/effet_fondateur/calibration/r13_evaluator.py"
    )
    assert len(STAGE_IDS) == 25


def test_frozen_r13_components_complete_ephemeral_technical_dry_run(
    tmp_path: Path,
) -> None:
    release_root = tmp_path / "r12b"
    workspace = tmp_path / "r13-dry-run"
    build_frozen_candidate(REPOSITORY_ROOT, release_root)

    manifest = run_blind_validation(REPOSITORY_ROOT, release_root, workspace)
    decision = evaluate_closed_run(REPOSITORY_ROOT, release_root, workspace)

    assert manifest["run_state"] == "CLOSED"
    assert manifest["partition_id"] == "VALIDATION_BLIND"
    assert manifest["dependencies"]["truth_opened"] is False
    assert {item["stage_id"] for item in manifest["stage_outputs"]} == set(STAGE_IDS)
    assert decision["status"] == "PASS"
    assert decision["scope"] == "TECHNICAL_INTEGRITY_ONLY"
    assert decision["scientific_performance_evaluated"] is False
    assert decision["r07_status"] == "REVISION_REQUIRED"
    assert not (release_root / "runs").exists()


def test_r12b_rejects_blind_artifact_tampering(tmp_path: Path) -> None:
    release_root = tmp_path / "r12b"
    build_frozen_candidate(REPOSITORY_ROOT, release_root)
    artifact = next((release_root / "blinded_inputs/VALIDATION_BLIND/views").glob("*.json"))
    artifact.write_text("{}\n", encoding="utf-8")

    with pytest.raises(
        (R12bValidationError, ReleaseGuardError),
        match="release_checksum_mismatch|artifact_",
    ):
        validate_frozen_candidate(REPOSITORY_ROOT, release_root)
