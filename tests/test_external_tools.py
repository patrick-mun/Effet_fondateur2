import subprocess

import pytest

from effet_fondateur.external_tools import (
    ExternalCommandError,
    ExternalToolConfigurationError,
    ExternalToolNotFoundError,
    resolve_executable,
    run_checked,
    run_checked_to_text_file,
    run_probe,
)


def test_resolve_executable_rejects_empty_configuration() -> None:
    with pytest.raises(ExternalToolConfigurationError, match="plink_not_configured"):
        resolve_executable("", "plink")


def test_resolve_executable_reports_missing_tool(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("effet_fondateur.external_tools.shutil.which", lambda _: None)

    with pytest.raises(ExternalToolNotFoundError, match="bcftools_not_found"):
        resolve_executable("bcftools", "bcftools")


def test_run_checked_preserves_nonzero_return_code() -> None:
    def runner(
        command: list[str], timeout_seconds: float
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(command, 7, "", "failure")

    with pytest.raises(ExternalCommandError) as error:
        run_checked(runner, ["tool", "--check"], 10, "tool_check")

    assert str(error.value) == "tool_check_failed:7"
    assert error.value.returncode == 7
    assert error.value.stdout == ""
    assert error.value.stderr == "failure"


def test_run_checked_normalizes_timeout() -> None:
    def runner(command: list[str], timeout_seconds: float):
        raise subprocess.TimeoutExpired(command, timeout_seconds)

    with pytest.raises(ExternalCommandError, match="tool_check_failed") as error:
        run_checked(runner, ["tool", "--check"], 10, "tool_check")

    assert error.value.returncode is None


def test_run_probe_returns_none_on_startup_error() -> None:
    def runner(command: list[str], timeout_seconds: float):
        raise OSError("unavailable")

    assert run_probe(runner, ["tool", "--version"], 10) is None


def test_run_probe_preserves_nonzero_result_for_caller_interpretation() -> None:
    result = subprocess.CompletedProcess(["tool", "--version"], 2, "version", "")

    assert run_probe(lambda command, timeout: result, ["tool", "--version"], 10) is result


def test_run_checked_to_text_file_writes_stdout(tmp_path) -> None:
    output_path = tmp_path / "query.tsv"

    run_checked_to_text_file(
        ["printf", "variant\\n"], 10, "query", output_path
    )

    assert output_path.read_text(encoding="utf-8") == "variant\n"
