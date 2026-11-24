"""Registre granulaire et reprenable des 1 180 appels scientifiques R13."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import sys
from typing import Any, Iterator, Mapping, Sequence

from effet_fondateur.calibration.contracts import canonical_sha256, file_sha256, load_json
from effet_fondateur.calibration.r13_scientific_runner import (
    RELEASE_ID,
    build_external_call_schedule,
)


LEDGER_RELATIVE = Path("governance/external_call_checkpoints.json")
TSV_RELATIVE = Path("governance/external_call_checkpoints.tsv")
MARKDOWN_RELATIVE = Path("governance/external_call_checkpoints.md")
LOCK_RELATIVE = Path("governance/external_call_checkpoints.lock")
AMENDMENT_RELATIVE = Path("governance/r13_execution_amendment.json")
ALLOWED_STATES = ("PENDING", "RUNNING", "COMPLETED", "FAILED")
MAX_BATCH_SIZE = 8
AMENDMENT_ID = "R13-LIGHTWEIGHT-INTEGRATION-001"


class R13CheckpointError(RuntimeError):
    """Registre absent, incohérent ou transition de checkpoint interdite."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13CheckpointError(code, detail)


def _timestamp() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _atomic_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.tmp-{os.getpid()}")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


@contextmanager
def _ledger_lock(run_root: Path) -> Iterator[None]:
    lock_path = run_root / LOCK_RELATIVE
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _schedule_rows() -> list[dict[str, Any]]:
    schedule = build_external_call_schedule()
    rows: list[dict[str, Any]] = []
    ordinal = 0
    for family in ("shapeit", "ibd"):
        for call in schedule[family]:
            ordinal += 1
            rows.append(
                {
                    "ordinal": ordinal,
                    "call_id": call["call_id"],
                    "family": family.upper(),
                    "density": call["density"],
                    "autosome_or_world": str(
                        call.get("autosome", call.get("world_id", ""))
                    ),
                    "tool": str(call.get("tool", "SHAPEIT5")),
                    "status": "PENDING",
                    "attempts": 0,
                    "claimed_by": None,
                    "started_at": None,
                    "finished_at": None,
                    "last_error": None,
                    "outputs": [],
                    "call_contract": dict(call),
                    "call_contract_sha256": canonical_sha256(call),
                }
            )
    return rows


def initialize_checkpoint_ledger(run_root: Path) -> dict[str, Any]:
    """Crée une seule fois le registre exact des 640 + 540 appels."""
    run_root = run_root.resolve(strict=False)
    with _ledger_lock(run_root):
        ledger_path = run_root / LEDGER_RELATIVE
        if ledger_path.exists():
            _fail("checkpoint_ledger_already_exists", str(ledger_path))
        schedule = build_external_call_schedule()
        ledger = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "ledger_state": "OPEN",
            "schedule_sha256": canonical_sha256(schedule),
            "call_count": 1_180,
            "calls": _schedule_rows(),
        }
        _atomic_json(ledger_path, ledger)
        _write_tables(run_root, ledger)
    return ledger


def validate_checkpoint_ledger(ledger: Mapping[str, Any]) -> dict[str, int]:
    if ledger.get("schema_version") != "1.0.0" or ledger.get("release_id") != RELEASE_ID:
        _fail("checkpoint_ledger_identity_mismatch")
    expected_schedule = build_external_call_schedule()
    if ledger.get("schedule_sha256") != canonical_sha256(expected_schedule):
        _fail("checkpoint_ledger_schedule_mismatch")
    calls = ledger.get("calls")
    if (
        ledger.get("call_count") != 1_180
        or not isinstance(calls, list)
        or len(calls) != 1_180
    ):
        _fail("checkpoint_ledger_call_count_mismatch")
    expected_rows = _schedule_rows()
    identifiers: set[str] = set()
    counts = {state: 0 for state in ALLOWED_STATES}
    for expected, row in zip(expected_rows, calls, strict=True):
        if not isinstance(row, Mapping):
            _fail("checkpoint_row_invalid")
        for field in (
            "ordinal",
            "call_id",
            "family",
            "density",
            "autosome_or_world",
            "tool",
            "call_contract",
            "call_contract_sha256",
        ):
            if row.get(field) != expected[field]:
                _fail("checkpoint_row_schedule_drift", str(row.get("call_id")))
        call_id = row["call_id"]
        if call_id in identifiers:
            _fail("checkpoint_call_id_duplicate", call_id)
        identifiers.add(call_id)
        status = row.get("status")
        if status not in ALLOWED_STATES:
            _fail("checkpoint_status_invalid", call_id)
        attempts = row.get("attempts")
        if isinstance(attempts, bool) or not isinstance(attempts, int) or attempts < 0:
            _fail("checkpoint_attempts_invalid", call_id)
        outputs = row.get("outputs")
        if not isinstance(outputs, list):
            _fail("checkpoint_outputs_invalid", call_id)
        if status == "COMPLETED" and not outputs:
            _fail("checkpoint_completed_outputs_missing", call_id)
        if status == "RUNNING" and (
            not row.get("claimed_by") or not row.get("started_at")
        ):
            _fail("checkpoint_running_provenance_missing", call_id)
        if status == "COMPLETED" and (
            not row.get("claimed_by")
            or not row.get("started_at")
            or not row.get("finished_at")
        ):
            _fail("checkpoint_completed_provenance_missing", call_id)
        if status == "FAILED" and (
            not row.get("claimed_by")
            or not row.get("started_at")
            or not row.get("finished_at")
            or not row.get("last_error")
        ):
            _fail("checkpoint_failed_provenance_missing", call_id)
        counts[status] += 1
    if sum(counts.values()) != 1_180:
        _fail("checkpoint_summary_count_mismatch")
    return counts


def load_checkpoint_ledger(run_root: Path) -> dict[str, Any]:
    ledger = load_json(run_root.resolve(strict=True) / LEDGER_RELATIVE)
    validate_checkpoint_ledger(ledger)
    return ledger


def _summary(ledger: Mapping[str, Any]) -> dict[str, int | float]:
    counts = validate_checkpoint_ledger(ledger)
    return {
        "total": 1_180,
        "pending": counts["PENDING"],
        "running": counts["RUNNING"],
        "completed": counts["COMPLETED"],
        "failed": counts["FAILED"],
        "remaining": counts["PENDING"] + counts["RUNNING"] + counts["FAILED"],
        "percent_complete": round(100 * counts["COMPLETED"] / 1_180, 2),
    }


def _write_tables(run_root: Path, ledger: Mapping[str, Any]) -> None:
    summary = _summary(ledger)
    tsv_path = run_root / TSV_RELATIVE
    markdown_path = run_root / MARKDOWN_RELATIVE
    columns = (
        "ordinal",
        "call_id",
        "family",
        "density",
        "autosome_or_world",
        "tool",
        "status",
        "attempts",
        "claimed_by",
        "started_at",
        "finished_at",
        "last_error",
        "done",
    )
    tsv_lines = ["\t".join(columns)]
    for row in ledger["calls"]:
        tsv_lines.append(
            "\t".join(
                (
                    ("OUI" if row["status"] == "COMPLETED" else "NON")
                    if column == "done"
                    else "" if row.get(column) is None
                    else str(row[column])
                )
                for column in columns
            )
        )
    _atomic_text(tsv_path, "\n".join(tsv_lines) + "\n")
    markdown_lines = [
        "# Suivi des checkpoints R13",
        "",
        f"- Total : {summary['total']}",
        f"- Terminés : {summary['completed']}",
        f"- En cours : {summary['running']}",
        f"- En attente : {summary['pending']}",
        f"- Échecs à reprendre : {summary['failed']}",
        f"- Restants : {summary['remaining']}",
        f"- Avancement : {summary['percent_complete']:.2f} %",
        "",
        "| # | Fait | Appel | Famille | Densité | Chromosome/monde | Outil | État | Tentatives |",
        "|---:|---|---|---|---|---|---|---|---:|",
    ]
    markdown_lines.extend(
        "| {ordinal} | {done} | `{call_id}` | {family} | {density} | {autosome_or_world} | {tool} | {status} | {attempts} |".format(
            done="OUI" if row["status"] == "COMPLETED" else "NON", **row
        )
        for row in ledger["calls"]
    )
    _atomic_text(markdown_path, "\n".join(markdown_lines) + "\n")


def _completed_outputs_valid(
    run_root: Path,
    row: Mapping[str, Any],
    *,
    verify_hashes: bool,
) -> bool:
    if row.get("status") != "COMPLETED":
        return True
    for output in row.get("outputs", []):
        if not isinstance(output, Mapping):
            return False
        relative = output.get("relative_path")
        if not isinstance(relative, str):
            return False
        path = run_root / relative
        try:
            resolved = path.resolve(strict=True)
            resolved.relative_to(run_root)
        except (FileNotFoundError, ValueError):
            return False
        if not resolved.is_file():
            return False
        if resolved.stat().st_size != output.get("size_bytes"):
            return False
        if verify_hashes and file_sha256(resolved) != output.get("sha256"):
            return False
    return True


def _mutate_ledger(run_root: Path, mutation: Any) -> dict[str, Any]:
    run_root = run_root.resolve(strict=True)
    with _ledger_lock(run_root):
        ledger_path = run_root / LEDGER_RELATIVE
        ledger = load_json(ledger_path)
        validate_checkpoint_ledger(ledger)
        mutation(ledger)
        validate_checkpoint_ledger(ledger)
        ledger["ledger_state"] = (
            "COMPLETE"
            if all(row["status"] == "COMPLETED" for row in ledger["calls"])
            else "OPEN"
        )
        _atomic_json(ledger_path, ledger)
        _write_tables(run_root, ledger)
    return ledger


def _assert_execution_open(run_root: Path) -> None:
    """Refuse toute nouvelle réservation après un amendement de clôture.

    La simple présence d'un overlay illisible ou incohérent bloque aussi
    l'exécution : une gouvernance endommagée ne doit jamais rouvrir des appels.
    """
    amendment_path = run_root.resolve(strict=True) / AMENDMENT_RELATIVE
    if not amendment_path.exists():
        return
    try:
        amendment = load_json(amendment_path)
    except Exception as error:
        raise R13CheckpointError(
            "checkpoint_execution_amendment_invalid", str(amendment_path)
        ) from error
    required = {
        "schema_version": "1.0.0",
        "amendment_id": AMENDMENT_ID,
        "release_id": RELEASE_ID,
        "run_id": "RUN-R13-SCIENTIFIC-001",
        "state": "EXECUTION_CLOSED_BY_AMENDMENT",
        "historical_call_count": 1_180,
        "retained_call_ids": ["SHAPEIT-DMAX-CHR01-SEED01"],
        "retained_completed_calls": 1,
        "withdrawn_unexecuted_calls": 1_179,
        "authorized_call_ids": [],
    }
    if any(amendment.get(key) != expected for key, expected in required.items()):
        _fail("checkpoint_execution_amendment_invalid", str(amendment_path))
    if amendment.get("original_ledger_sha256") != file_sha256(
        run_root.resolve(strict=True) / LEDGER_RELATIVE
    ):
        _fail("checkpoint_execution_amendment_invalid", "original_ledger_sha256")
    _fail("checkpoint_execution_closed_by_amendment", AMENDMENT_ID)


def claim_next_batch(
    run_root: Path,
    *,
    batch_size: int,
    worker_id: str,
    retry_failed: bool = True,
    timestamp: str | None = None,
) -> list[dict[str, Any]]:
    """Réserve atomiquement les prochains appels non terminés."""
    if (
        isinstance(batch_size, bool)
        or not isinstance(batch_size, int)
        or not 1 <= batch_size <= MAX_BATCH_SIZE
    ):
        _fail("checkpoint_batch_size_invalid")
    if not isinstance(worker_id, str) or not worker_id:
        _fail("checkpoint_worker_id_invalid")
    _assert_execution_open(run_root)
    claimed: list[dict[str, Any]] = []
    allowed = {"PENDING", "FAILED"} if retry_failed else {"PENDING"}

    def mutation(ledger: dict[str, Any]) -> None:
        for row in ledger["calls"]:
            if row["status"] == "COMPLETED" and not _completed_outputs_valid(
                run_root.resolve(strict=True), row, verify_hashes=False
            ):
                row["status"] = "FAILED"
                row["finished_at"] = timestamp or _timestamp()
                row["last_error"] = "COMPLETED_OUTPUT_INTEGRITY_MISMATCH"
        for row in ledger["calls"]:
            if row["status"] not in allowed:
                continue
            row["status"] = "RUNNING"
            row["attempts"] += 1
            row["claimed_by"] = worker_id
            row["started_at"] = timestamp or _timestamp()
            row["finished_at"] = None
            row["last_error"] = None
            row["outputs"] = []
            claimed.append(dict(row))
            if len(claimed) == batch_size:
                break

    _mutate_ledger(run_root, mutation)
    return claimed


def interrupt_running_calls(
    run_root: Path,
    *,
    worker_id: str | None = None,
    timestamp: str | None = None,
) -> list[dict[str, Any]]:
    """Rend reprenables les appels restés RUNNING après un arrêt propre."""
    interrupted: list[dict[str, Any]] = []

    def mutation(ledger: dict[str, Any]) -> None:
        for row in ledger["calls"]:
            if row["status"] != "RUNNING":
                continue
            if worker_id is not None and row["claimed_by"] != worker_id:
                continue
            row["status"] = "FAILED"
            row["finished_at"] = timestamp or _timestamp()
            row["last_error"] = "INTERRUPTED_BY_USER"
            row["outputs"] = []
            interrupted.append(dict(row))

    _mutate_ledger(run_root, mutation)
    return interrupted


def _find_row(ledger: Mapping[str, Any], call_id: str) -> dict[str, Any]:
    matches = [row for row in ledger["calls"] if row["call_id"] == call_id]
    if len(matches) != 1:
        _fail("checkpoint_call_unknown", call_id)
    return matches[0]


def complete_call(
    run_root: Path,
    *,
    call_id: str,
    output_paths: Sequence[Path],
    timestamp: str | None = None,
) -> dict[str, Any]:
    """Marque un appel terminé après empreinte de toutes ses sorties."""
    run_root = run_root.resolve(strict=True)
    if not output_paths:
        _fail("checkpoint_completion_outputs_missing", call_id)
    outputs = []
    for output_path in output_paths:
        resolved = output_path.resolve(strict=True)
        try:
            relative = resolved.relative_to(run_root)
        except ValueError as error:
            raise R13CheckpointError("checkpoint_output_outside_run", call_id) from error
        outputs.append(
            {
                "relative_path": relative.as_posix(),
                "sha256": file_sha256(resolved),
                "size_bytes": resolved.stat().st_size,
            }
        )
    completed: dict[str, Any] = {}

    def mutation(ledger: dict[str, Any]) -> None:
        row = _find_row(ledger, call_id)
        if row["status"] != "RUNNING":
            _fail("checkpoint_completion_state_invalid", call_id)
        row["status"] = "COMPLETED"
        row["finished_at"] = timestamp or _timestamp()
        row["outputs"] = outputs
        completed.update(row)

    _mutate_ledger(run_root, mutation)
    return completed


def fail_call(
    run_root: Path,
    *,
    call_id: str,
    error_code: str,
    timestamp: str | None = None,
) -> dict[str, Any]:
    """Conserve un échec explicite afin qu'il soit repris au lot suivant."""
    if not isinstance(error_code, str) or not error_code:
        _fail("checkpoint_failure_code_invalid", call_id)
    failed: dict[str, Any] = {}

    def mutation(ledger: dict[str, Any]) -> None:
        row = _find_row(ledger, call_id)
        if row["status"] != "RUNNING":
            _fail("checkpoint_failure_state_invalid", call_id)
        row["status"] = "FAILED"
        row["finished_at"] = timestamp or _timestamp()
        row["last_error"] = error_code
        row["outputs"] = []
        failed.update(row)

    _mutate_ledger(run_root, mutation)
    return failed


def checkpoint_status(
    run_root: Path,
    *,
    verify_output_hashes: bool = False,
) -> dict[str, Any]:
    resolved_root = run_root.resolve(strict=True)
    ledger = load_checkpoint_ledger(resolved_root)
    invalid_completed = sum(
        row["status"] == "COMPLETED"
        and not _completed_outputs_valid(
            resolved_root,
            row,
            verify_hashes=verify_output_hashes,
        )
        for row in ledger["calls"]
    )
    return {
        "ledger_state": ledger["ledger_state"],
        **_summary(ledger),
        "invalid_completed_outputs": invalid_completed,
        "resume_safe": invalid_completed == 0,
        "output_hashes_verified": verify_output_hashes,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=("init", "status", "claim", "complete", "fail", "interrupt"),
    )
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--worker-id", default="LOCAL-8CPU")
    parser.add_argument("--call-id")
    parser.add_argument("--output", type=Path, action="append", default=[])
    parser.add_argument("--error-code")
    parser.add_argument("--only-worker-id")
    parser.add_argument("--verify-output-hashes", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.action == "init":
        initialize_checkpoint_ledger(arguments.run_root)
        document = checkpoint_status(arguments.run_root)
    elif arguments.action == "status":
        document = checkpoint_status(
            arguments.run_root,
            verify_output_hashes=arguments.verify_output_hashes,
        )
    elif arguments.action == "claim":
        document = {
            "claimed": claim_next_batch(
                arguments.run_root,
                batch_size=arguments.batch_size,
                worker_id=arguments.worker_id,
            )
        }
    elif arguments.action == "complete":
        if not arguments.call_id:
            _fail("checkpoint_call_id_required")
        document = complete_call(
            arguments.run_root,
            call_id=arguments.call_id,
            output_paths=arguments.output,
        )
    elif arguments.action == "fail":
        if not arguments.call_id or not arguments.error_code:
            _fail("checkpoint_failure_arguments_required")
        document = fail_call(
            arguments.run_root,
            call_id=arguments.call_id,
            error_code=arguments.error_code,
        )
    else:
        document = {
            "interrupted": interrupt_running_calls(
                arguments.run_root,
                worker_id=arguments.only_worker_id,
            )
        }
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
