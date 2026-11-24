from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.contracts import STAGE_IDS
from effet_fondateur.calibration.r13_scientific_runner import (
    EXPECTED_SCHEDULE,
    R13ScientificRunnerError,
    assert_full_campaign_capacity,
    build_external_call_schedule,
    execute_frozen_plan,
    validate_execution_plan,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _plan() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "release_id": "panel-master-generic-v1-r12c-scientific-validation",
        "partition_id": "VALIDATION_BLIND",
        "scope": "R13_SCIENTIFIC_BLIND_VALIDATION",
        "stage_ids": list(STAGE_IDS),
        "required_workers": 8,
        "schedule": copy.deepcopy(EXPECTED_SCHEDULE),
        "external_calls": build_external_call_schedule(),
        "truth_available_to_runner": False,
        "work_items": [
            {
                "work_item_id": f"WORK-{stage}",
                "stage_id": stage,
                "command": [".venv/bin/python", "-m", "frozen.stage", stage],
                "dependencies": [],
                "expected_outputs": [f"runs/RUN-R13/{stage}/metrics.json"],
            }
            for stage in STAGE_IDS
        ],
    }


def test_scientific_runner_accepts_only_closed_blind_full_schedule() -> None:
    audit = validate_execution_plan(_plan())

    assert audit["stage_count"] == 25
    assert audit["work_item_count"] == 25
    assert audit["shapeit_call_count"] == 640
    assert audit["ibd_call_count"] == 540
    assert audit["required_workers"] == 8
    assert audit["runner_implementation_state"] == "FULL_EXECUTION_READY"
    assert audit["execution_authorized"] is False


def test_scientific_runner_rejects_truth_access() -> None:
    plan = _plan()
    plan["truth_available_to_runner"] = True

    with pytest.raises(
        (R13ScientificRunnerError, ValueError),
        match="truth|forbidden",
    ):
        validate_execution_plan(plan)


def test_scientific_runner_rejects_reduced_schedule() -> None:
    plan = _plan()
    plan["schedule"]["shapeit"]["total"] = 639  # type: ignore[index]

    with pytest.raises(R13ScientificRunnerError, match="schedule_mismatch"):
        validate_execution_plan(plan)


def test_scientific_runner_builds_exact_r07_external_call_inventory() -> None:
    schedule = build_external_call_schedule()

    assert len(schedule["shapeit"]) == 640
    assert len(schedule["ibd"]) == 540
    assert sum(call["density"] == "D_MAX" for call in schedule["shapeit"]) == 340
    assert sum(call["density"] == "D_MID" for call in schedule["shapeit"]) == 150
    assert sum(call["density"] == "D_LOW" for call in schedule["shapeit"]) == 150
    assert len({call["world_id"] for call in schedule["ibd"]}) == 90
    assert len({call["call_id"] for call in schedule["ibd"]}) == 540


def test_scientific_runner_rejects_external_call_inventory_drift() -> None:
    plan = _plan()
    plan["external_calls"]["shapeit"].pop()  # type: ignore[index]

    with pytest.raises(R13ScientificRunnerError, match="schedule_mismatch"):
        validate_execution_plan(plan)


def test_scientific_runner_refuses_host_below_eight_cpus() -> None:
    with pytest.raises(R13ScientificRunnerError, match="underprovisioned"):
        assert_full_campaign_capacity(7)


def test_scientific_runner_accepts_local_eight_cpu_capacity() -> None:
    assert_full_campaign_capacity(8)


def test_scientific_runner_executes_and_resumes_eight_worker_batches(
    tmp_path: Path,
) -> None:
    plan = _plan()
    for item in plan["work_items"]:  # type: ignore[index]
        output = item["expected_outputs"][0]
        marker = f"launches/{item['work_item_id']}"
        item["command"] = [
            str(REPOSITORY_ROOT / ".venv/bin/python"),
            "-c",
            (
                "import os; from pathlib import Path; "
                "r=Path(os.environ['EFFET_FONDATEUR_R13_RUN_ROOT']); "
                f"m=r/{marker!r}; m.parent.mkdir(parents=True, exist_ok=True); m.mkdir(); "
                f"p=r/{output!r}; "
                "p.parent.mkdir(parents=True, exist_ok=True); p.write_text('closed\\n')"
            ),
        ]
    receipt_path = tmp_path / "g5.json"
    receipt_path.write_text(
        json.dumps(
            {
                "receipt_id": "R12C-SCIENTIFIC-G5-READY",
                "release_id": "panel-master-generic-v1-r12c-scientific-validation",
                "gate_g5_passed": True,
                "scientific_validation_authorized": True,
            }
        ),
        encoding="utf-8",
    )
    run_root = tmp_path / "run"

    first = execute_frozen_plan(REPOSITORY_ROOT, plan, run_root, receipt_path)
    resumed = execute_frozen_plan(REPOSITORY_ROOT, plan, run_root, receipt_path)

    assert first["state"] == "CLOSED_EXECUTION_COMPLETE"
    assert resumed["completed_work_item_count"] == 25
    assert len(list((run_root / "launches").iterdir())) == 25
