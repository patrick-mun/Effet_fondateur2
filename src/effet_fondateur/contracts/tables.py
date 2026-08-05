"""Lecture et validation génériques des tables TSV versionnées."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jsonschema

from effet_fondateur.audit import sha256_file
from effet_fondateur.contracts.documents import DEFAULT_SCHEMAS_DIR


class TableValidationError(ValueError):
    """Signale une structure TSV ou une contrainte de ligne invalide."""


@dataclass(frozen=True)
class ValidatedTable:
    """Table TSV validée, convertie selon son schéma et associée à son empreinte."""

    path: Path
    schema_name: str
    schema_version: str
    sha256: str
    rows: tuple[dict[str, Any], ...]

    @property
    def row_count(self) -> int:
        """Retourne le nombre de lignes de données validées."""
        return len(self.rows)


def _load_table_schema(schema_name: str, schemas_dir: Path) -> dict[str, Any]:
    """Charge un schéma possédant les extensions tabulaires V2 attendues."""
    schema_path = schemas_dir / schema_name
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise TableValidationError(f"Schéma de table introuvable : {schema_path}") from error
    except json.JSONDecodeError as error:
        raise TableValidationError(f"Schéma de table JSON invalide : {schema_path}") from error
    if "x-tsv" not in schema or "row" not in schema.get("$defs", {}):
        raise TableValidationError(f"Schéma sans contrat TSV : {schema_path}")
    return schema


def _resolve_property_schema(
    property_schema: dict[str, Any], schema: dict[str, Any]
) -> dict[str, Any]:
    """Résout les références locales nécessaires au décodage des cellules."""
    reference = property_schema.get("$ref")
    if reference is None:
        return property_schema
    reference_prefix = "#/$defs/"
    if not reference.startswith(reference_prefix):
        raise TableValidationError(f"Référence de schéma TSV non locale : {reference}")
    definition_name = reference.removeprefix(reference_prefix)
    try:
        return schema["$defs"][definition_name]
    except KeyError as error:
        raise TableValidationError(f"Définition de schéma TSV absente : {reference}") from error


def _allows_null(property_schema: dict[str, Any]) -> bool:
    property_type = property_schema.get("type")
    if property_type == "null":
        return True
    if isinstance(property_type, list) and "null" in property_type:
        return True
    return any(option.get("type") == "null" for option in property_schema.get("oneOf", []))


def _expects_boolean(property_schema: dict[str, Any]) -> bool:
    property_type = property_schema.get("type")
    return property_type == "boolean" or (
        isinstance(property_type, list) and "boolean" in property_type
    )


def _decode_value(
    raw_value: str,
    column_schema: dict[str, Any],
    null_value: str,
    boolean_values: dict[str, bool],
    schema: dict[str, Any],
) -> Any:
    """Convertit uniquement les valeurs dont le type est explicite dans le schéma."""
    column_schema = _resolve_property_schema(column_schema, schema)
    if raw_value == null_value and _allows_null(column_schema):
        return None
    if _expects_boolean(column_schema) and raw_value in boolean_values:
        return boolean_values[raw_value]
    return raw_value


def _validation_column(error: jsonschema.ValidationError) -> str:
    path_parts = list(error.absolute_path)
    if path_parts:
        return str(path_parts[0])
    return "<ligne>"


def validate_tsv_table(
    table_path: Path,
    schema_name: str,
    schemas_dir: Path = DEFAULT_SCHEMAS_DIR,
) -> ValidatedTable:
    """Valide l'en-tête et chaque ligne TSV sans journaliser leur contenu."""
    schema = _load_table_schema(schema_name, schemas_dir)
    table_contract = schema["x-tsv"]
    expected_columns = table_contract["columns"]
    row_schema = {
        "$schema": schema["$schema"],
        "$defs": schema["$defs"],
        **schema["$defs"]["row"],
    }
    try:
        jsonschema.Draft202012Validator.check_schema(row_schema)
    except jsonschema.SchemaError as error:
        raise TableValidationError(f"Schéma de ligne non conforme : {schema_name}") from error
    validator = jsonschema.Draft202012Validator(row_schema)
    observed_primary_keys: set[tuple[Any, ...]] = set()

    try:
        input_file = table_path.open("r", encoding="utf-8", newline="")
    except FileNotFoundError as error:
        raise TableValidationError(f"Table TSV introuvable : {table_path}") from error

    rows: list[dict[str, Any]] = []
    with input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if reader.fieldnames != expected_columns:
            raise TableValidationError(
                f"En-tête invalide pour {table_path}; ordre ou colonnes non conformes."
            )
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row or any(raw_row[column] is None for column in expected_columns):
                raise TableValidationError(
                    f"Ligne {line_number} : nombre de champs différent de l'en-tête."
                )
            decoded_row = {
                column: _decode_value(
                    raw_row[column],
                    schema["$defs"]["row"]["properties"][column],
                    table_contract["null_value"],
                    table_contract["boolean_values"],
                    schema,
                )
                for column in expected_columns
            }
            errors = sorted(
                validator.iter_errors(decoded_row),
                key=lambda item: list(item.absolute_path),
            )
            if errors:
                first_error = errors[0]
                # La valeur individuelle n'est volontairement pas recopiée dans
                # l'erreur : ligne, colonne et contrainte suffisent au diagnostic.
                raise TableValidationError(
                    f"Ligne {line_number}, colonne {_validation_column(first_error)} : "
                    f"contrainte {first_error.validator} non satisfaite."
                )
            primary_key = tuple(
                decoded_row[column] for column in table_contract["primary_key"]
            )
            if primary_key in observed_primary_keys:
                raise TableValidationError(
                    f"Ligne {line_number} : clé primaire TSV dupliquée."
                )
            observed_primary_keys.add(primary_key)
            rows.append(decoded_row)

    if not rows and not table_contract.get("allow_empty", False):
        raise TableValidationError(f"Table TSV sans ligne de données : {table_path}")
    return ValidatedTable(
        path=table_path,
        schema_name=schema_name,
        schema_version=table_contract["schema_version"],
        sha256=sha256_file(table_path),
        rows=tuple(rows),
    )
