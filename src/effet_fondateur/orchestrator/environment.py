"""Inventaire technique non sensible de l'environnement d'un run."""

from __future__ import annotations

import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from effet_fondateur import __version__


TOOL_VERSION_ARGUMENTS = {
    "plink": ("--version",),
    "king": ("--version",),
    "bcftools": ("--version",),
    "rscript": ("--version",),
}


def capture_tool_environment(
    command_name: str | None,
    tool_key: str,
) -> dict[str, Any]:
    """Détecte un outil sans rendre son absence bloquante au bootstrap."""
    if command_name is None:
        return {"configured": None, "available": False, "version": None}
    resolved_path = shutil.which(command_name)
    if resolved_path is None:
        return {"configured": command_name, "available": False, "version": None}

    version = None
    arguments = TOOL_VERSION_ARGUMENTS.get(tool_key, ("--version",))
    try:
        completed_process = subprocess.run(
            [resolved_path, *arguments],
            capture_output=True,
            check=False,
            text=True,
            timeout=5,
        )
        version_output = completed_process.stdout.strip() or completed_process.stderr.strip()
        if version_output:
            version = version_output.splitlines()[0][:500]
    except (OSError, subprocess.TimeoutExpired):
        version = None
    # L'étape qui utilise réellement l'outil décidera si une absence ou une
    # version illisible est bloquante. Le bootstrap reste générique.
    return {"configured": command_name, "available": True, "version": version}


def build_environment(config: dict[str, Any]) -> dict[str, Any]:
    """Construit l'inventaire Python, plateforme et outils déclaré dans le run."""
    return {
        "schema_version": "1.0.0",
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": Path(sys.executable).name,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "pipeline_version": __version__,
        "tools": {
            tool_key: capture_tool_environment(command_name, tool_key)
            for tool_key, command_name in config["tools"].items()
        },
    }
