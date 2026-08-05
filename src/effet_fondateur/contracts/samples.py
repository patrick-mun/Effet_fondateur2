"""Contraintes métier de la table maître des échantillons."""

from __future__ import annotations

from pathlib import Path, PurePosixPath

from effet_fondateur.contracts.tables import (
    TableValidationError,
    ValidatedTable,
    validate_tsv_table,
)


SAMPLES_SCHEMA_NAME = "samples_master.schema.json"
FORBIDDEN_GENOTYPE_SOURCES = {"clinical_status", "group_label"}


def _validate_source_file(source_file: str, source_root: Path) -> None:
    """Confine une source déclarée sous la racine autorisée, liens inclus."""
    declared_path = PurePosixPath(source_file)
    if (
        declared_path.is_absolute()
        or ".." in declared_path.parts
        or "\\" in source_file
    ):
        raise TableValidationError("SOURCE_FILE doit être un chemin relatif sûr.")
    resolved_root = source_root.resolve()
    resolved_source = source_root.joinpath(*declared_path.parts).resolve()
    if not resolved_source.is_relative_to(resolved_root):
        raise TableValidationError("SOURCE_FILE sort du répertoire source autorisé.")
    if not resolved_source.is_file():
        raise TableValidationError("Un SOURCE_FILE déclaré est absent ou n'est pas un fichier.")


def validate_samples_master(
    table_path: Path,
    source_root: Path,
) -> ValidatedTable:
    """Valide le registre maître, ses clés, filiations et fichiers sources."""
    validated_table = validate_tsv_table(table_path, SAMPLES_SCHEMA_NAME)
    sample_ids: set[str] = set()
    plink_ids: set[tuple[str, str]] = set()
    source_files: set[str] = set()
    rows_by_plink_id = {
        (row["FID"], row["IID"]): row for row in validated_table.rows
    }

    for row in validated_table.rows:
        if row["SAMPLE_ID"] in sample_ids:
            raise TableValidationError("SAMPLE_ID dupliqué dans la table maître.")
        sample_ids.add(row["SAMPLE_ID"])

        plink_id = (row["FID"], row["IID"])
        if plink_id in plink_ids:
            raise TableValidationError("Couple FID + IID dupliqué dans la table maître.")
        plink_ids.add(plink_id)

        if row["SOURCE_FILE"] in source_files:
            raise TableValidationError(
                "SOURCE_FILE est assigné à plusieurs individus dans la table maître."
            )
        source_files.add(row["SOURCE_FILE"])
        _validate_source_file(row["SOURCE_FILE"], source_root)
        genotype_source = row["TARGET_GENOTYPE_SOURCE"]
        # Cette défense reste dans le code en plus du schéma : une future
        # extension du vocabulaire ne doit jamais réintroduire cette inférence.
        if (
            genotype_source is not None
            and genotype_source.casefold() in FORBIDDEN_GENOTYPE_SOURCES
        ):
            raise TableValidationError(
                "TARGET_GENOTYPE_SOURCE ne peut pas provenir du statut clinique ou du groupe."
            )
        if row["PID"] == row["IID"] or row["MID"] == row["IID"]:
            raise TableValidationError("Un individu ne peut pas être son propre parent.")
        if row["PID"] != "0" and row["PID"] == row["MID"]:
            raise TableValidationError("PID et MID ne peuvent pas désigner le même parent.")

    for row in validated_table.rows:
        father = rows_by_plink_id.get((row["FID"], row["PID"]))
        mother = rows_by_plink_id.get((row["FID"], row["MID"]))
        if father is not None and father["SEX"] == "FEMALE":
            raise TableValidationError("Un PID présent dans la table est déclaré FEMALE.")
        if mother is not None and mother["SEX"] == "MALE":
            raise TableValidationError("Un MID présent dans la table est déclaré MALE.")
    return validated_table
