import json
from pathlib import Path

import pytest
import yaml

from effet_fondateur.orchestrator import IntegrityError, StageExecutionError
from effet_fondateur.orchestrator.pipeline import resume_pipeline, run_pipeline


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def write_test_config(
    path: Path,
    *,
    target_complete: bool = True,
    synthetic_parameters: dict[str, object] | None = None,
) -> None:
    config = yaml.safe_load(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
            encoding="utf-8"
        )
    )
    if not target_complete:
        for field in ("position_bp", "ref", "alt", "project_variant_id"):
            config["target"][field] = None
    if synthetic_parameters is not None:
        config["stages"]["synthetic_stage"] = {
            "enabled": True,
            "parameters": synthetic_parameters,
        }
    path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def read_manifest(run_dir: Path) -> dict[str, object]:
    return json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))


def test_run_initializes_immutable_audited_structure(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path)

    run_dir = run_pipeline(config_path, tmp_path / "runs")

    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "TECHNICALLY_VALID"
    assert manifest["stages"][0]["state"] == "SUCCEEDED"
    assert (run_dir / "config.resolved.yaml").is_file()
    assert (run_dir / "environment.json").is_file()
    assert (run_dir / "events.jsonl").is_file()
    assert (run_dir / "stages" / "00_initialize_run" / "audit.json").is_file()


def test_incomplete_target_keeps_run_incomplete(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path, target_complete=False)

    run_dir = run_pipeline(config_path, tmp_path / "runs")

    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "INCOMPLETE"
    assert manifest["manual_decisions_required"] == ["target_variant_definition"]


def test_synthetic_stage_runs_in_subprocess_and_is_reused(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path, synthetic_parameters={"message": "validation-ok"})
    run_dir = run_pipeline(config_path, tmp_path / "runs")
    result_path = run_dir / "stages" / "T00_synthetic_stage" / "synthetic_result.json"
    modification_time = result_path.stat().st_mtime_ns

    resume_pipeline(run_dir)

    assert json.loads(result_path.read_text(encoding="utf-8"))["message"] == "validation-ok"
    assert result_path.stat().st_mtime_ns == modification_time
    command = json.loads(
        (run_dir / "stages" / "T00_synthetic_stage" / "command.json").read_text(
            encoding="utf-8"
        )
    )
    assert command["shell"] is False
    assert command["module"] == "effet_fondateur.stages.synthetic_stage"
    events = (run_dir / "events.jsonl").read_text(encoding="utf-8")
    assert "stage_reused" in events


def test_failed_stage_is_audited_and_can_be_resumed(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path, synthetic_parameters={"fail_attempts": 1})
    runs_dir = tmp_path / "runs"

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4
    run_dir = next(path for path in runs_dir.iterdir() if not path.name.startswith("."))
    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "BLOCKED"
    assert manifest["stages"][-1]["state"] == "FAILED"
    assert manifest["stages"][-1]["attempt_count"] == 1
    assert any((run_dir / "stages" / "attempts").iterdir())

    resume_pipeline(run_dir)

    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "TECHNICALLY_VALID"
    assert manifest["stages"][-1]["state"] == "SUCCEEDED"
    assert manifest["stages"][-1]["attempt_count"] == 2


def test_resume_blocks_modified_resolved_configuration(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path, synthetic_parameters={"message": "original"})
    run_dir = run_pipeline(config_path, tmp_path / "runs")
    resolved_config_path = run_dir / "config.resolved.yaml"
    resolved_config_path.write_text(
        resolved_config_path.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )

    with pytest.raises(IntegrityError, match="configuration résolue"):
        resume_pipeline(run_dir)

    assert read_manifest(run_dir)["global_status"] == "BLOCKED"


def test_resume_blocks_modified_published_artifact(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path, synthetic_parameters={"message": "original"})
    run_dir = run_pipeline(config_path, tmp_path / "runs")
    result_path = run_dir / "stages" / "T00_synthetic_stage" / "synthetic_result.json"
    result_path.write_text('{"message": "modified"}\n', encoding="utf-8")

    with pytest.raises(IntegrityError, match="modifié"):
        resume_pipeline(run_dir)

    manifest = read_manifest(run_dir)
    assert manifest["global_status"] == "BLOCKED"
    assert manifest["stages"][-1]["state"] == "FAILED"


def test_resume_blocks_modified_stage_audit(tmp_path: Path) -> None:
    config_path = tmp_path / "config.yaml"
    write_test_config(config_path, synthetic_parameters={"message": "original"})
    run_dir = run_pipeline(config_path, tmp_path / "runs")
    audit_path = run_dir / "stages" / "T00_synthetic_stage" / "audit.json"
    audit_path.write_text(audit_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(IntegrityError, match="Audit d'étape modifié"):
        resume_pipeline(run_dir)
