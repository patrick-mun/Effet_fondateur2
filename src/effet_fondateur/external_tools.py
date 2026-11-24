"""Primitives communes pour invoquer les outils externes du pipeline."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Callable, Sequence


CommandRunner = Callable[[Sequence[str], float], subprocess.CompletedProcess[str]]


class ExternalToolConfigurationError(ValueError):
    """Signale qu'une commande externe n'est pas configurée."""


class ExternalToolNotFoundError(RuntimeError):
    """Signale qu'un exécutable configuré est absent du ``PATH``."""


class ExternalCommandError(RuntimeError):
    """Signale l'échec d'exécution ou le code retour non nul d'une commande."""

    def __init__(
        self,
        operation: str,
        returncode: int | None = None,
        *,
        stdout: str = "",
        stderr: str = "",
        timed_out: bool = False,
    ) -> None:
        self.operation = operation
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.timed_out = timed_out
        suffix = "" if returncode is None else f":{returncode}"
        super().__init__(f"{operation}_failed{suffix}")


def default_command_runner(
    command: Sequence[str], timeout_seconds: float
) -> subprocess.CompletedProcess[str]:
    """Exécute une commande sans interpréteur de shell et capture ses sorties."""
    return subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )


def resolve_executable(command: str, tool: str) -> str:
    """Résout une commande configurée en chemin exécutable explicite."""
    if not isinstance(command, str) or not command.strip():
        raise ExternalToolConfigurationError(f"{tool}_not_configured")
    executable = shutil.which(command)
    if executable is None:
        raise ExternalToolNotFoundError(f"{tool}_not_found")
    return executable


def run_checked(
    runner: CommandRunner,
    command: Sequence[str],
    timeout_seconds: float,
    operation: str,
) -> subprocess.CompletedProcess[str]:
    """Exécute une commande et normalise ses erreurs techniques."""
    try:
        completed = runner(command, timeout_seconds)
    except subprocess.TimeoutExpired as error:
        raise ExternalCommandError(operation, timed_out=True) from error
    except OSError as error:
        raise ExternalCommandError(operation) from error
    if completed.returncode != 0:
        raise ExternalCommandError(
            operation,
            completed.returncode,
            stdout=completed.stdout or "",
            stderr=completed.stderr or "",
        )
    return completed


def run_checked_to_text_file(
    command: Sequence[str],
    timeout_seconds: float,
    operation: str,
    output_path: Path,
) -> None:
    """Exécute une commande en écrivant sa sortie standard dans un nouveau fichier."""
    try:
        with output_path.open("x", encoding="utf-8") as output_handle:
            completed = subprocess.run(
                list(command),
                check=False,
                stdout=output_handle,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout_seconds,
            )
    except subprocess.TimeoutExpired as error:
        raise ExternalCommandError(operation, timed_out=True) from error
    except OSError as error:
        raise ExternalCommandError(operation) from error
    if completed.returncode != 0:
        raise ExternalCommandError(
            operation,
            completed.returncode,
            stderr=completed.stderr or "",
        )


def run_probe(
    runner: CommandRunner,
    command: Sequence[str],
    timeout_seconds: float,
) -> subprocess.CompletedProcess[str] | None:
    """Exécute une sonde non bloquante et retourne ``None`` si elle ne démarre pas."""
    try:
        return runner(command, timeout_seconds)
    except (OSError, subprocess.TimeoutExpired):
        return None
