from __future__ import annotations

import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.r13_checkpoint_ledger import (
    LEDGER_RELATIVE,
    MARKDOWN_RELATIVE,
    TSV_RELATIVE,
    R13CheckpointError,
    checkpoint_status,
    claim_next_batch,
    complete_call,
    fail_call,
    initialize_checkpoint_ledger,
    interrupt_running_calls,
    load_checkpoint_ledger,
    main,
)


def test_checkpoint_ledger_materializes_1180_call_table(tmp_path: Path) -> None:
    run_root = tmp_path / "run"

    ledger = initialize_checkpoint_ledger(run_root)
    status = checkpoint_status(run_root)

    assert len(ledger["calls"]) == 1_180
    assert status == {
        "ledger_state": "OPEN",
        "total": 1_180,
        "pending": 1_180,
        "running": 0,
        "completed": 0,
        "failed": 0,
        "remaining": 1_180,
        "percent_complete": 0.0,
        "invalid_completed_outputs": 0,
        "resume_safe": True,
        "output_hashes_verified": False,
    }
    tsv_lines = (run_root / TSV_RELATIVE).read_text(encoding="utf-8").splitlines()
    markdown = (run_root / MARKDOWN_RELATIVE).read_text(encoding="utf-8")
    assert len(tsv_lines) == 1_181
    assert tsv_lines[0].endswith("\tdone")
    assert tsv_lines[1].endswith("\tNON")
    assert "| # | Fait | Appel |" in markdown
    assert "| 1 | NON | `SHAPEIT-DMAX-CHR01-SEED01` |" in markdown


def test_checkpoint_batches_are_claimed_completed_and_retried_safely(
    tmp_path: Path,
) -> None:
    run_root = tmp_path / "run"
    initialize_checkpoint_ledger(run_root)

    claimed = claim_next_batch(
        run_root,
        batch_size=8,
        worker_id="LOCAL-8CPU",
        timestamp="2026-11-14T10:00:00Z",
    )
    first_id = claimed[0]["call_id"]
    second_id = claimed[1]["call_id"]
    output = run_root / "outputs" / f"{first_id}.json"
    output.parent.mkdir(parents=True)
    output.write_text('{"closed": true}\n', encoding="utf-8")

    completed = complete_call(
        run_root,
        call_id=first_id,
        output_paths=[output],
        timestamp="2026-11-14T10:10:00Z",
    )
    failed = fail_call(
        run_root,
        call_id=second_id,
        error_code="TOOL_EXIT_1",
        timestamp="2026-11-14T10:11:00Z",
    )
    interrupt_running_calls(
        run_root,
        worker_id="LOCAL-8CPU",
        timestamp="2026-11-14T10:12:00Z",
    )
    retried = claim_next_batch(
        run_root,
        batch_size=1,
        worker_id="LOCAL-8CPU",
        timestamp="2026-11-14T10:13:00Z",
    )

    assert completed["status"] == "COMPLETED"
    assert completed["outputs"][0]["relative_path"] == f"outputs/{first_id}.json"
    assert failed["status"] == "FAILED"
    assert retried[0]["call_id"] == second_id
    assert retried[0]["attempts"] == 2
    assert checkpoint_status(run_root)["completed"] == 1


def test_missing_completed_output_is_reclaimed(tmp_path: Path) -> None:
    run_root = tmp_path / "run"
    initialize_checkpoint_ledger(run_root)
    call = claim_next_batch(run_root, batch_size=1, worker_id="WORKER")[0]
    output = run_root / "outputs/result.txt"
    output.parent.mkdir(parents=True)
    output.write_text("closed\n", encoding="utf-8")
    complete_call(run_root, call_id=call["call_id"], output_paths=[output])
    output.unlink()

    assert checkpoint_status(run_root)["resume_safe"] is False
    reclaimed = claim_next_batch(run_root, batch_size=1, worker_id="WORKER")

    assert reclaimed[0]["call_id"] == call["call_id"]
    assert reclaimed[0]["attempts"] == 2


def test_checkpoint_ledger_rejects_schedule_tampering_and_large_batches(
    tmp_path: Path,
) -> None:
    run_root = tmp_path / "run"
    initialize_checkpoint_ledger(run_root)

    with pytest.raises(R13CheckpointError, match="batch_size_invalid"):
        claim_next_batch(run_root, batch_size=9, worker_id="WORKER")

    ledger_path = run_root / LEDGER_RELATIVE
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["calls"][0]["density"] = "D_LOW"
    ledger_path.write_text(json.dumps(ledger), encoding="utf-8")

    with pytest.raises(R13CheckpointError, match="schedule_drift"):
        load_checkpoint_ledger(run_root)


def test_status_cli_can_verify_output_hashes(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    run_root = tmp_path / "run"
    initialize_checkpoint_ledger(run_root)

    exit_code = main(
        [
            "status",
            "--run-root",
            str(run_root),
            "--verify-output-hashes",
        ]
    )
    result = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert result["output_hashes_verified"] is True
