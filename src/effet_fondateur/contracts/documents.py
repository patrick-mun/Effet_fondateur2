"""Validation des documents JSON versionnés du pipeline V2."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema


DEFAULT_SCHEMAS_DIR = Path(__file__).resolve().parents[3] / "schemas"


class DocumentValidationError(ValueError):
    """Signale qu'un document ne respecte pas son schéma versionné."""


def validate_json_document(
    document: dict[str, Any],
    schema_name: str,
    schemas_dir: Path = DEFAULT_SCHEMAS_DIR,
) -> None:
    """Valide un objet Python avec un schéma JSON local nommé explicitement."""
    schema_path = schemas_dir / schema_name
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise DocumentValidationError(f"Schéma introuvable : {schema_path}") from error
    except json.JSONDecodeError as error:
        raise DocumentValidationError(f"Schéma JSON invalide : {schema_path}") from error

    try:
        jsonschema.Draft202012Validator.check_schema(schema)
    except jsonschema.SchemaError as error:
        raise DocumentValidationError(f"Schéma non conforme : {schema_path}") from error

    validator = jsonschema.Draft202012Validator(
        schema,
        format_checker=jsonschema.Draft202012Validator.FORMAT_CHECKER,
    )
    errors = sorted(validator.iter_errors(document), key=lambda item: list(item.absolute_path))
    if errors:
        first_error = errors[0]
        validation_path = ".".join(str(part) for part in first_error.absolute_path)
        validation_path = validation_path or "<racine>"
        raise DocumentValidationError(
            f"Document invalide pour {schema_name} à {validation_path} : "
            f"{first_error.message}"
        )
