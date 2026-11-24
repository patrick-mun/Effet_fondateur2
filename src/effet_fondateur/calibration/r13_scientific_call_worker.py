"""Worker unitaire figé pour un appel externe scientifique R13."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path, PurePosixPath
import resource
import shutil
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence

from effet_fondateur.calibration.contracts import file_sha256, load_json, load_yaml
from effet_fondateur.calibration.r13_checkpoint_ledger import (
    complete_call,
    fail_call,
    load_checkpoint_ledger,
)
from effet_fondateur.calibration.r13_scientific_materializer import (
    CONFIGURATION_PATH,
    load_materializer_configuration,
    materializer_contract_audit,
    resolve_shapeit_chromosome,
)
from effet_fondateur.calibration.r13_scientific_runner import RELEASE_ID
from effet_fondateur.phasing.shapeit5 import SHAPEIT5_VERSION_PATTERN


WORKER_STATE = "FULL_SINGLE_CALL_EXECUTION_READY"
SHAPEIT_RELATIVE = Path("calibration_workspace/toolchains/r11/bin/SHAPEIT5_phase_common")
JAVA_RELATIVE = Path(
    "calibration_workspace/toolchains/r11/java/"
    "jdk-17.0.21.0.2+1/Contents/Home/bin/java"
)
HAP_IBD_RELATIVE = Path("calibration_workspace/toolchains/r11/jars/hap-ibd.jar")
REFINED_IBD_RELATIVE = Path(
    "calibration_workspace/toolchains/r11/jars/refined-ibd.17Jan20.102.jar"
)
EXPECTED_TOOL_HASHES = {
    "shapeit": "88c8c19060e9d5fdff6538aa7bd4ab7d165a179777565fb035e11d53a72f7e6e",
    "hap_ibd": "846ddf08db5590339c58072fa9944fff8466f9c2825f9c02612cb9d38492752c",
    "refined_ibd": "a2c2d7ee6c1dff5c06831a5987ba3a3ce831cb5798a13d413ab0c10a577db919",
}


class R13CallWorkerError(RuntimeError):
    """Appel inconnu, environnement invalide ou commande externe en échec."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13CallWorkerError(code, detail)


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _tool_paths(repository_root: Path) -> dict[str, Path]:
    repository_root = repository_root.resolve(strict=True)
    bcftools_name = shutil.which("bcftools")
    if bcftools_name is None:
        _fail("r13_worker_bcftools_missing")
    paths = {
        "bcftools": Path(bcftools_name).resolve(strict=True),
        "shapeit": (repository_root / SHAPEIT_RELATIVE).resolve(strict=True),
        "java": (repository_root / JAVA_RELATIVE).resolve(strict=True),
        "hap_ibd": (repository_root / HAP_IBD_RELATIVE).resolve(strict=True),
        "refined_ibd": (repository_root / REFINED_IBD_RELATIVE).resolve(strict=True),
    }
    for name, expected in EXPECTED_TOOL_HASHES.items():
        if file_sha256(paths[name]) != expected:
            _fail("r13_worker_tool_hash_mismatch", name)
    return paths


def _probe_shapeit_version(path: Path) -> str:
    completed = subprocess.run(
        [str(path), "--help"],
        capture_output=True,
        text=True,
        check=False,
        timeout=20,
    )
    output = "\n".join((completed.stdout, completed.stderr))
    match = SHAPEIT5_VERSION_PATTERN.search(output)
    if completed.returncode != 0 or match is None:
        _fail("r13_worker_shapeit_version_unreadable")
    if match.group(1) != "5.1.1":
        _fail("r13_worker_shapeit_version_mismatch", match.group(1))
    return match.group(1)


def _probe_version(command: Sequence[str], code: str) -> str:
    completed = subprocess.run(
        list(command), capture_output=True, text=True, check=False, timeout=20
    )
    output = "\n".join((completed.stdout, completed.stderr)).strip().splitlines()
    if completed.returncode != 0 or not output:
        _fail(code)
    return output[0]


def worker_contract_audit(repository_root: Path) -> dict[str, Any]:
    tools = _tool_paths(repository_root)
    version = _probe_shapeit_version(tools["shapeit"])
    materializer = materializer_contract_audit(repository_root)
    return {
        "state": WORKER_STATE,
        "release_id": RELEASE_ID,
        "shapeit_version": version,
        "bcftools_version": _probe_version(
            [str(tools["bcftools"]), "--version"],
            "r13_worker_bcftools_version_unreadable",
        ),
        "java_version": _probe_version(
            [str(tools["java"]), "-version"],
            "r13_worker_java_version_unreadable",
        ),
        "tool_sha256": {
            name: file_sha256(tools[name])
            for name in ("shapeit", "hap_ibd", "refined_ibd")
        },
        "supported_call_count": materializer["supported_call_count"],
        "shell_used": False,
        "truth_available_to_worker": False,
        "checkpoint_transitions": ["RUNNING_TO_COMPLETED", "RUNNING_TO_FAILED"],
    }


def _ledger_row(run_root: Path, call_id: str) -> dict[str, Any]:
    ledger = load_checkpoint_ledger(run_root)
    matches = [row for row in ledger["calls"] if row["call_id"] == call_id]
    if len(matches) != 1:
        _fail("r13_worker_call_unknown", call_id)
    return matches[0]


def _validate_g5(repository_root: Path, release_root: Path) -> dict[str, Any]:
    receipt = load_json(repository_root / "calibration/readiness/R12C_G5_READY.json")
    if (
        receipt.get("receipt_id") != "R12C-SCIENTIFIC-G5-READY"
        or receipt.get("release_id") != RELEASE_ID
        or receipt.get("gate_g5_passed") is not True
        or receipt.get("scientific_validation_authorized") is not True
        or receipt.get("materializer_state") != "FULL_CALL_MATERIALIZATION_READY"
        or receipt.get("single_call_worker_state") != WORKER_STATE
    ):
        _fail("r13_worker_g5_invalid")
    frozen = load_json(release_root / "governance/r12c_frozen_manifest.json")
    frozen_path = release_root / "governance/r12c_frozen_manifest.json"
    checksums_path = release_root / "release_checksums.sha256"
    evidence = receipt.get("evidence_sha256")
    if (
        frozen.get("materializer_state") != "FULL_CALL_MATERIALIZATION_READY"
        or frozen.get("single_call_worker_state") != WORKER_STATE
        or not isinstance(evidence, Mapping)
        or evidence.get("r12c_frozen_manifest") != file_sha256(frozen_path)
        or evidence.get("release_checksums") != file_sha256(checksums_path)
    ):
        _fail("r13_worker_release_materializer_not_frozen")
    for binding in frozen.get("frozen_components", []):
        relative = binding.get("relative_path") if isinstance(binding, Mapping) else None
        if (
            not isinstance(relative, str)
            or not (repository_root / relative).is_file()
            or file_sha256(repository_root / relative) != binding.get("sha256")
        ):
            _fail("r13_worker_frozen_component_drift", str(relative))
    return receipt


def _artifact_map(run_root: Path, manifest: Mapping[str, Any]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        _fail("r13_worker_materialization_artifacts_missing")
    for artifact in artifacts:
        if not isinstance(artifact, Mapping):
            _fail("r13_worker_materialization_artifact_invalid")
        relative = artifact.get("relative_path")
        if not isinstance(relative, str):
            _fail("r13_worker_materialization_path_invalid")
        pure = PurePosixPath(relative)
        if pure.is_absolute() or ".." in pure.parts or "private_truth" in {
            part.casefold() for part in pure.parts
        }:
            _fail("r13_worker_private_truth_path_forbidden")
        path = run_root / relative
        if (
            not path.is_file()
            or path.stat().st_size != artifact.get("size_bytes")
            or file_sha256(path) != artifact.get("sha256")
        ):
            _fail("r13_worker_materialization_integrity_mismatch", relative)
        result[path.name] = path
    return result


def _materialization_manifest_path(
    repository_root: Path,
    run_root: Path,
    call_contract: Mapping[str, Any],
) -> Path:
    configuration = load_materializer_configuration(repository_root)
    call_id = str(call_contract["call_id"])
    density = str(call_contract["density"])
    if call_id.startswith("SHAPEIT-"):
        chromosome = resolve_shapeit_chromosome(call_contract, configuration)
        return (
            run_root
            / "blinded_materialized/shapeit"
            / f"chr{chromosome:02d}"
            / density
            / "manifest.json"
        )
    return (
        run_root
        / "blinded_materialized/ibd"
        / str(call_contract["world_id"])
        / density
        / "manifest.json"
    )


def preflight_call(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    call_id: str,
) -> dict[str, Any]:
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    run_root = run_root.resolve(strict=True)
    _validate_g5(repository_root, release_root)
    tools = _tool_paths(repository_root)
    row = _ledger_row(run_root, call_id)
    if row["status"] not in {"PENDING", "RUNNING", "FAILED"}:
        _fail("r13_worker_call_not_preflightable", call_id)
    manifest_path = _materialization_manifest_path(
        repository_root, run_root, row["call_contract"]
    )
    if not manifest_path.is_file():
        _fail("r13_worker_materialization_missing", str(manifest_path))
    manifest = load_json(manifest_path)
    if (
        manifest.get("release_id") != RELEASE_ID
        or manifest.get("density") != row["density"]
        or manifest.get("truth_available_to_worker") is not False
    ):
        _fail("r13_worker_materialization_contract_mismatch", call_id)
    artifacts = _artifact_map(run_root, manifest)
    required = (
        {"study.vcf.gz", "study.vcf.gz.tbi", "reference.vcf.gz", "reference.vcf.gz.tbi", "pedigree.tsv", "genetic_map.tsv.gz"}
        if row["family"] == "SHAPEIT"
        else {"phased.vcf.gz", "phased.vcf.gz.tbi", "ibd.map"}
    )
    if not required.issubset(artifacts):
        _fail("r13_worker_required_materialization_missing", call_id)
    return {
        "call_id": call_id,
        "status": "READY_FOR_SINGLE_CALL_EXECUTION",
        "checkpoint_state": row["status"],
        "attempts": row["attempts"],
        "materialization_manifest": manifest_path.relative_to(run_root).as_posix(),
        "materialization_artifact_count": len(artifacts),
        "shapeit_version": _probe_shapeit_version(tools["shapeit"]),
        "truth_available_to_worker": False,
    }


def _build_command(
    repository_root: Path,
    run_root: Path,
    row: Mapping[str, Any],
    output_root: Path,
) -> tuple[list[str], list[Path]]:
    tools = _tool_paths(repository_root)
    configuration = load_yaml(repository_root / CONFIGURATION_PATH)
    manifest = load_json(
        _materialization_manifest_path(repository_root, run_root, row["call_contract"])
    )
    artifacts = _artifact_map(run_root, manifest)
    call_contract = row["call_contract"]
    if row["family"] == "SHAPEIT":
        output = output_root / "phased.bcf"
        log = output_root / "shapeit.log"
        seed = int(configuration["shapeit"]["seed_base"]) + int(
            call_contract["seed_index"]
        )
        command = [
            str(tools["shapeit"]),
            "--input", str(artifacts["study.vcf.gz"]),
            "--reference", str(artifacts["reference.vcf.gz"]),
            "--map", str(artifacts["genetic_map.tsv.gz"]),
            "--region", str(manifest["chromosome"]),
            "--output", str(output),
            "--output-format", "bcf",
            "--log", str(log),
            "--thread", str(configuration["shapeit"]["threads_per_call"]),
            "--seed", str(seed),
            "--pedigree", str(artifacts["pedigree.tsv"]),
        ]
        return command, [output, log]
    prefix = output_root / "segments"
    tool_key = "hap_ibd" if call_contract["tool"] == "HAP_IBD" else "refined_ibd"
    command = [
        str(tools["java"]),
        f"-Xmx{int(configuration['ibd']['heap_mib'])}m",
        "-jar", str(tools[tool_key]),
        f"gt={artifacts['phased.vcf.gz']}",
        f"map={artifacts['ibd.map']}",
        f"out={prefix}",
    ]
    if tool_key == "hap_ibd":
        command.extend(("min-seed=1.0", "min-output=1.0", "min-markers=20", "nthreads=4"))
    else:
        command.extend(("length=1.0", "nthreads=4"))
    return command, [Path(f"{prefix}.ibd.gz")]


def _run_command(
    command: Sequence[str],
    *,
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    try:
        with (
            stdout_path.open("w", encoding="utf-8") as stdout_handle,
            stderr_path.open("w", encoding="utf-8") as stderr_handle,
        ):
            completed = subprocess.run(
                list(command),
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
                timeout=timeout_seconds,
            )
    except subprocess.TimeoutExpired as error:
        raise R13CallWorkerError("r13_worker_call_timeout", command[0]) from error
    finished = time.monotonic()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "returncode": completed.returncode,
        "wall_seconds": finished - started,
        "user_cpu_seconds": after.ru_utime - before.ru_utime,
        "system_cpu_seconds": after.ru_stime - before.ru_stime,
        "maximum_rss_native": after.ru_maxrss,
    }


def _index_shapeit_output(bcftools: Path, output: Path) -> None:
    if Path(f"{output}.csi").is_file():
        return
    indexed = subprocess.run(
        [str(bcftools), "index", "--csi", str(output)],
        capture_output=True,
        text=True,
        check=False,
    )
    if indexed.returncode != 0:
        _fail("r13_worker_output_index_failed", indexed.stderr[-1000:])


def execute_call(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    call_id: str,
) -> dict[str, Any]:
    """Exécute exactement un checkpoint déjà réservé et clôt son état."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    run_root = run_root.resolve(strict=True)
    preflight_call(repository_root, release_root, run_root, call_id)
    row = _ledger_row(run_root, call_id)
    if row["status"] != "RUNNING" or row["attempts"] < 1:
        _fail("r13_worker_call_must_be_claimed", call_id)
    attempt_root = (
        run_root
        / "outputs/external_calls"
        / call_id
        / f"attempt-{int(row['attempts']):02d}"
    )
    if attempt_root.exists():
        _fail("r13_worker_attempt_collision", str(attempt_root))
    attempt_root.mkdir(parents=True)
    stdout_path = attempt_root / "stdout.log"
    stderr_path = attempt_root / "stderr.log"
    receipt_path = attempt_root / "receipt.json"
    command, expected_outputs = _build_command(
        repository_root, run_root, row, attempt_root
    )
    configuration = load_materializer_configuration(repository_root)
    timeout = int(
        configuration["shapeit" if row["family"] == "SHAPEIT" else "ibd"][
            "timeout_seconds"
        ]
    )
    try:
        measurement = _run_command(
            command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            timeout_seconds=timeout,
        )
        if measurement["returncode"] != 0:
            _fail("r13_worker_external_nonzero", f"returncode={measurement['returncode']}")
        tools = _tool_paths(repository_root)
        if row["family"] == "SHAPEIT":
            _index_shapeit_output(tools["bcftools"], expected_outputs[0])
            expected_outputs.append(Path(f"{expected_outputs[0]}.csi"))
        if any(not path.is_file() for path in expected_outputs):
            _fail("r13_worker_expected_output_missing", call_id)
        receipt = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "call_id": call_id,
            "attempt": row["attempts"],
            "state": "COMPLETED",
            "completed_at": _timestamp(),
            "command": command,
            "shell_used": False,
            "truth_available_to_worker": False,
            "measurement": measurement,
            "outputs": [
                {
                    "relative_path": path.relative_to(run_root).as_posix(),
                    "sha256": file_sha256(path),
                    "size_bytes": path.stat().st_size,
                }
                for path in expected_outputs
            ],
        }
        _write_json(receipt_path, receipt)
        complete_call(
            run_root,
            call_id=call_id,
            output_paths=[*expected_outputs, stdout_path, stderr_path, receipt_path],
        )
        return receipt
    except BaseException as error:
        failure = {
            "schema_version": "1.0.0",
            "release_id": RELEASE_ID,
            "call_id": call_id,
            "attempt": row["attempts"],
            "state": "FAILED",
            "failed_at": _timestamp(),
            "error_type": type(error).__name__,
            "error": str(error),
        }
        _write_json(receipt_path, failure)
        current = _ledger_row(run_root, call_id)
        if current["status"] == "RUNNING":
            fail_call(
                run_root,
                call_id=call_id,
                error_code=(error.code if isinstance(error, R13CallWorkerError) else "UNEXPECTED_WORKER_ERROR"),
            )
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("audit", "preflight", "execute"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--call-id")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.action == "audit":
        document = worker_contract_audit(arguments.repository_root)
    else:
        if not all((arguments.release_root, arguments.run_root, arguments.call_id)):
            _fail("r13_worker_cli_arguments_missing")
        function = preflight_call if arguments.action == "preflight" else execute_call
        document = function(
            arguments.repository_root,
            arguments.release_root,
            arguments.run_root,
            arguments.call_id,
        )
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
