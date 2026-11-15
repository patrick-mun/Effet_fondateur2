from __future__ import annotations

from pathlib import Path
import shutil

from effet_fondateur.calibration.r13_checkpoint_ledger import (
    checkpoint_status,
    claim_next_batch,
    initialize_checkpoint_ledger,
)
from effet_fondateur.calibration.r13_scientific_call_worker import (
    execute_call,
    preflight_call,
    worker_contract_audit,
)
from effet_fondateur.calibration.r13_scientific_materializer import (
    SmokeDimensions,
    materialize_shapeit_shard,
)
from effet_fondateur.calibration.r13_scientific_runner import (
    build_external_call_schedule,
)

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = (
    REPOSITORY_ROOT
    / "calibration_workspace/releases/panel_master_generic_v1_r12c_scientific_validation"
)


def _fake_compress_and_index(source: Path, output: Path, _bcftools: Path) -> None:
    shutil.copyfile(source, output)
    Path(f"{output}.tbi").write_text("synthetic-index\n", encoding="utf-8")


def test_worker_contract_binds_tools_and_all_calls() -> None:
    audit = worker_contract_audit(REPOSITORY_ROOT)

    assert audit["state"] == "FULL_SINGLE_CALL_EXECUTION_READY"
    assert audit["shapeit_version"] == "5.1.1"
    assert audit["supported_call_count"] == 1_180
    assert audit["shell_used"] is False
    assert audit["truth_available_to_worker"] is False


def test_single_worker_preflights_and_closes_one_claimed_checkpoint(
    tmp_path: Path,
    monkeypatch,
) -> None:
    run_root = tmp_path / "run"
    truth_root = tmp_path / "truth"
    initialize_checkpoint_ledger(run_root)
    call = build_external_call_schedule()["shapeit"][0]
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_materializer._compress_and_index_vcf",
        _fake_compress_and_index,
    )
    materialize_shapeit_shard(
        REPOSITORY_ROOT,
        RELEASE_ROOT,
        run_root,
        truth_root,
        call,
        test_dimensions=SmokeDimensions(12, 9, 2, 120),
    )
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_call_worker._validate_g5",
        lambda *_args: {},
    )
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_call_worker._probe_shapeit_version",
        lambda _path: "5.1.1",
    )
    fake_tools = {
        "bcftools": tmp_path / "bcftools",
        "shapeit": tmp_path / "shapeit",
        "java": tmp_path / "java",
        "hap_ibd": tmp_path / "hap.jar",
        "refined_ibd": tmp_path / "refined.jar",
    }
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_call_worker._tool_paths",
        lambda _root: fake_tools,
    )

    def fake_run(command, *, stdout_path, stderr_path, timeout_seconds):
        output = Path(command[command.index("--output") + 1])
        log = Path(command[command.index("--log") + 1])
        output.write_text("synthetic-bcf\n", encoding="utf-8")
        log.write_text("synthetic-shapeit-log\n", encoding="utf-8")
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text("", encoding="utf-8")
        return {
            "returncode": 0,
            "wall_seconds": 1.0,
            "user_cpu_seconds": 0.8,
            "system_cpu_seconds": 0.1,
            "maximum_rss_native": 1024,
        }

    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_call_worker._run_command",
        fake_run,
    )
    monkeypatch.setattr(
        "effet_fondateur.calibration.r13_scientific_call_worker._index_shapeit_output",
        lambda _bcftools, output: Path(f"{output}.csi").write_text(
            "synthetic-index\n", encoding="utf-8"
        ),
    )

    before = preflight_call(
        REPOSITORY_ROOT, RELEASE_ROOT, run_root, call["call_id"]
    )
    claimed = claim_next_batch(run_root, batch_size=1, worker_id="TEST-WORKER")
    receipt = execute_call(
        REPOSITORY_ROOT, RELEASE_ROOT, run_root, call["call_id"]
    )

    assert before["status"] == "READY_FOR_SINGLE_CALL_EXECUTION"
    assert claimed[0]["call_id"] == call["call_id"]
    assert receipt["state"] == "COMPLETED"
    assert receipt["truth_available_to_worker"] is False
    assert checkpoint_status(run_root)["completed"] == 1
