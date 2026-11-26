import subprocess

import pytest

from effet_fondateur.external_tools import ExternalToolNotFoundError
from effet_fondateur.orchestrator import environment


def test_missing_tool_remains_a_non_blocking_environment_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*_args, **_kwargs):
        raise ExternalToolNotFoundError("bcftools_not_found")

    monkeypatch.setattr(environment, "resolve_executable", fail)

    assert environment.capture_tool_environment("bcftools", "bcftools") == {
        "configured": "bcftools",
        "available": False,
        "version": None,
    }


def test_failed_probe_keeps_the_tool_available_without_a_version(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        environment,
        "resolve_executable",
        lambda command, _tool: f"/tools/{command}",
    )
    monkeypatch.setattr(environment, "run_probe", lambda *_args, **_kwargs: None)

    assert environment.capture_tool_environment("bcftools", "bcftools") == {
        "configured": "bcftools",
        "available": True,
        "version": None,
    }


def test_shapeit_probe_keeps_the_version_line(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        environment,
        "resolve_executable",
        lambda command, _tool: f"/tools/{command}",
    )
    completed = subprocess.CompletedProcess(
        ["shapeit5", "--help"],
        0,
        "header\n  * Version       : 5.1.1 / commit = test\n",
        "",
    )
    monkeypatch.setattr(
        environment,
        "run_probe",
        lambda *_args, **_kwargs: completed,
    )

    result = environment.capture_tool_environment(
        "shapeit5_phase_common", "shapeit5_phase_common"
    )

    assert result["available"] is True
    assert result["version"] == "* Version       : 5.1.1 / commit = test"
