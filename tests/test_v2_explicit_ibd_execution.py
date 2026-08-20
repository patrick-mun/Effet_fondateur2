import subprocess
from pathlib import Path

import pytest

from effet_fondateur.audit import sha256_file
from effet_fondateur.explicit_ibd.execution import ExplicitIbdExternalError, run_tool, validate_adapters
from effet_fondateur.orchestrator.environment import build_environment


def _adapters(tmp_path: Path):
    hap, refined = tmp_path / "hap.jar", tmp_path / "refined.jar"
    hap.write_bytes(b"hap")
    refined.write_bytes(b"refined")
    return {"java_command": "java", "expected_java_major": 17, "hap_ibd_jar": str(hap), "hap_ibd_version": "1.0", "hap_ibd_sha256": sha256_file(hap), "refined_ibd_jar": str(refined), "refined_ibd_version": "1.0", "refined_ibd_sha256": sha256_file(refined)}


def test_adapter_refuses_bad_jar_hash(tmp_path: Path):
    adapters = _adapters(tmp_path)
    adapters["hap_ibd_sha256"] = "0" * 64
    with pytest.raises(ExplicitIbdExternalError, match="sha256"):
        validate_adapters(adapters)


def test_tool_timeout_is_controlled(tmp_path: Path):
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 1)
    with pytest.raises(ExplicitIbdExternalError, match="timeout"):
        run_tool(tool="HAP_IBD", adapters=_adapters(tmp_path), vcf_path=tmp_path / "x.vcf.gz", map_path=tmp_path / "map", output_prefix=tmp_path / "out", minimum_cm=2.0, minimum_markers=100, threads=1, memory_mb=256, timeout_seconds=1, runner=timeout)


def test_nonzero_tool_exit_is_controlled(tmp_path: Path):
    def failed(*args, **kwargs):
        return subprocess.CompletedProcess(args[0], 7, "", "failure")
    with pytest.raises(ExplicitIbdExternalError, match="failed:7"):
        run_tool(tool="REFINED_IBD", adapters=_adapters(tmp_path), vcf_path=tmp_path / "x.vcf.gz", map_path=tmp_path / "map", output_prefix=tmp_path / "out", minimum_cm=2.0, minimum_markers=100, threads=1, memory_mb=256, timeout_seconds=1, runner=failed)


def test_environment_handles_structured_explicit_ibd_adapter(tmp_path: Path):
    adapters = _adapters(tmp_path)
    config = {
        "tools": {
            "explicit_ibd_adapters": adapters,
        }
    }

    environment = build_environment(config)["tools"]["explicit_ibd_adapters"]

    assert environment["expected_java_major"] == 17
    assert environment["hap_ibd"]["sha256"] == adapters["hap_ibd_sha256"]
    assert environment["refined_ibd"]["version"] == "1.0"
