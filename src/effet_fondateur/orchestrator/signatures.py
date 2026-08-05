"""Construction des signatures déterministes des étapes V2."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.orchestrator.errors import StageExecutionError
from effet_fondateur.orchestrator.models import StageDefinition


def module_sha256(module_name: str) -> str:
    """Empreinte le fichier Python réellement importé pour une étape."""
    module_specification = importlib.util.find_spec(module_name)
    if module_specification is None or module_specification.origin is None:
        raise StageExecutionError(f"Module d'étape introuvable : {module_name}", 2)
    return sha256_file(Path(module_specification.origin))


def build_stage_signature(
    definition: StageDefinition,
    parameters: dict[str, Any],
    config_sha256: str,
    input_artifacts: list[dict[str, Any]] | None = None,
) -> str:
    """Lie le code, la configuration, les paramètres et les contrats d'E/S."""
    # Les versions de schéma font partie de la signature : un contrat modifié
    # invalide volontairement le cache même si le script Python n'a pas changé.
    signature_payload = {
        "stage_id": definition.stage_id,
        "stage_name": definition.stage_name,
        "module_sha256": module_sha256(definition.module),
        "parameters": parameters,
        "config_sha256": config_sha256,
        "input_artifacts": input_artifacts or [],
        "input_schema_version": "1.0.0",
        "output_schema_version": "1.0.0",
        "audit_schema_version": "1.0.0",
    }
    serialized_payload = json.dumps(
        signature_payload,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(serialized_payload.encode("utf-8")).hexdigest()
