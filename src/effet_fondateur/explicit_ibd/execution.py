"""Adaptateurs reproductibles Hap-IBD et Refined IBD."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from effet_fondateur.audit import sha256_file


class ExplicitIbdExternalError(RuntimeError):
    """Signale une indisponibilité, un timeout ou une sortie externe invalide."""


@dataclass(frozen=True)
class ToolRun:
    tool: str
    command: tuple[str, ...]
    output_path: Path
    log_path: Path


def validate_adapters(adapters: dict[str, Any] | None) -> dict[str, Any]:
    """Contrôle Java, les JAR épinglés et leurs empreintes avant tout appel."""
    if adapters is None:
        raise ExplicitIbdExternalError("explicit_ibd_adapters_not_configured")
    for name in ("hap_ibd", "refined_ibd"):
        jar = Path(adapters[f"{name}_jar"])
        if jar.is_symlink() or not jar.is_file():
            raise ExplicitIbdExternalError(f"{name}_jar_missing")
        if sha256_file(jar) != adapters[f"{name}_sha256"]:
            raise ExplicitIbdExternalError(f"{name}_jar_sha256_mismatch")
    try:
        probe = subprocess.run([adapters["java_command"], "-version"], capture_output=True, text=True, check=False, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExplicitIbdExternalError("java_unavailable_or_timeout") from error
    version_text = probe.stderr + probe.stdout
    match = re.search(r'version "(\d+)', version_text)
    if probe.returncode != 0 or match is None or int(match.group(1)) != adapters["expected_java_major"]:
        raise ExplicitIbdExternalError("java_version_mismatch")
    return adapters


def run_tool(
    *, tool: str, adapters: dict[str, Any], vcf_path: Path, map_path: Path,
    output_prefix: Path, minimum_cm: float, minimum_markers: int,
    threads: int, memory_mb: int, timeout_seconds: int,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> ToolRun:
    """Exécute un JAR borné et exige la sortie `.ibd.gz` attendue."""
    if tool not in {"HAP_IBD", "REFINED_IBD"}:
        raise ExplicitIbdExternalError("explicit_ibd_tool_invalid")
    key = "hap_ibd" if tool == "HAP_IBD" else "refined_ibd"
    command = [
        adapters["java_command"], f"-Xmx{memory_mb}m", "-jar", adapters[f"{key}_jar"],
        f"gt={vcf_path}", f"map={map_path}", f"out={output_prefix}",
    ]
    if tool == "HAP_IBD":
        command.extend((f"min-seed={minimum_cm}", f"min-output={minimum_cm}", f"min-markers={minimum_markers}", f"nthreads={threads}"))
    else:
        command.extend((f"length={minimum_cm}", f"nthreads={threads}"))
    log_path = output_prefix.with_suffix(".log")
    try:
        completed = runner(command, capture_output=True, text=True, check=False, timeout=timeout_seconds)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExplicitIbdExternalError(f"{key}_unavailable_or_timeout") from error
    log_path.write_text((completed.stdout or "") + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode != 0:
        raise ExplicitIbdExternalError(f"{key}_failed:{completed.returncode}")
    output_path = Path(f"{output_prefix}.ibd.gz")
    if not output_path.is_file():
        raise ExplicitIbdExternalError(f"{key}_output_missing")
    return ToolRun(tool, tuple(command), output_path, log_path)
