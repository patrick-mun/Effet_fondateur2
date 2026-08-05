"""Chargement strict de la configuration du pipeline V2."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
import yaml


DEFAULT_SCHEMA_PATH = Path(__file__).resolve().parents[3] / "schemas" / "pipeline_config.schema.json"


class ConfigurationError(ValueError):
    """Signale une configuration absente, illisible ou contraire au schéma."""


def _format_validation_path(error: jsonschema.ValidationError) -> str:
    path_parts = [str(part) for part in error.absolute_path]
    return ".".join(path_parts) if path_parts else "<racine>"


def load_pipeline_config(
    config_path: Path,
    schema_path: Path = DEFAULT_SCHEMA_PATH,
) -> dict[str, Any]:
    """Charge un YAML et le valide sans résoudre ni modifier ses chemins."""
    try:
        raw_config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ConfigurationError(f"Configuration introuvable : {config_path}") from error
    except yaml.YAMLError as error:
        raise ConfigurationError(f"YAML invalide dans {config_path} : {error}") from error

    if not isinstance(raw_config, dict):
        raise ConfigurationError(
            f"La configuration {config_path} doit contenir un objet YAML à la racine."
        )

    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ConfigurationError(f"Schéma introuvable : {schema_path}") from error
    except json.JSONDecodeError as error:
        raise ConfigurationError(f"Schéma JSON invalide : {schema_path}") from error

    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(raw_config), key=lambda item: list(item.absolute_path))
    if errors:
        first_error = errors[0]
        validation_path = _format_validation_path(first_error)
        raise ConfigurationError(
            f"Configuration invalide à {validation_path} : {first_error.message}"
        )

    return raw_config
