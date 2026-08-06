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
    "shapeit5_phase_common": ("--help",),
    "shapeit5_phase_rare": ("--help",),
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
            version_lines = version_output.splitlines()
            if tool_key.startswith("shapeit5_"):
                version_line = next(
                    (line.strip() for line in version_lines if "* Version" in line),
                    version_lines[0],
                )
                version = version_line[:500]
            else:
                version = version_lines[0][:500]
    except (OSError, subprocess.TimeoutExpired):
        version = None
    # L'étape qui utilise réellement l'outil décidera si une absence ou une
    # version illisible est bloquante. Le bootstrap reste générique.
    return {"configured": command_name, "available": True, "version": version}


def build_environment(config: dict[str, Any]) -> dict[str, Any]:
    """Construit l'inventaire Python, plateforme et outils déclaré dans le run."""
    tools_environment: dict[str, Any] = {}
    for tool_key, tool_config in config["tools"].items():
        if tool_key == "phasing_adapter" and isinstance(tool_config, dict):
            tools_environment[tool_key] = {
                "adapter_id": tool_config["adapter_id"],
                "expected_version": tool_config["expected_version"],
                "phase_common": capture_tool_environment(
                    tool_config["phase_common_command"], "shapeit5_phase_common"
                ),
                "phase_rare": capture_tool_environment(
                    tool_config["phase_rare_command"], "shapeit5_phase_rare"
                ),
            }
        else:
            tools_environment[tool_key] = capture_tool_environment(tool_config, tool_key)
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
        "tools": tools_environment,
    }
