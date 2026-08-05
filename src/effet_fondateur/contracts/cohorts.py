"""Contraintes métier des cohortes analytiques figées."""

from __future__ import annotations

from pathlib import Path

from effet_fondateur.contracts.tables import (
    TableValidationError,
    ValidatedTable,
    validate_tsv_table,
)


COHORTS_SCHEMA_NAME = "cohorts_frozen.schema.json"


def validate_cohorts_frozen(
    table_path: Path,
    known_sample_ids: set[str] | None = None,
) -> ValidatedTable:
    """Valide les décisions de cohorte sans déduire de rôle depuis les groupes."""
    validated_table = validate_tsv_table(table_path, COHORTS_SCHEMA_NAME)
    cohort_memberships: set[tuple[str, str]] = set()
    for row in validated_table.rows:
        membership = (row["COHORT_ID"], row["SAMPLE_ID"])
        if membership in cohort_memberships:
            raise TableValidationError(
                "Couple COHORT_ID + SAMPLE_ID dupliqué dans les cohortes."
            )
        cohort_memberships.add(membership)
        if known_sample_ids is not None and row["SAMPLE_ID"] not in known_sample_ids:
            raise TableValidationError(
                "Une cohorte référence un SAMPLE_ID absent de la table maître."
            )
        if row["FAMILY_REPRESENTATIVE"] and row["INDEPENDENT_UNIT_ID"] is None:
            # Le représentant familial n'est utile analytiquement que si toutes
            # les étapes peuvent retrouver l'unité indépendante correspondante.
            raise TableValidationError(
                "FAMILY_REPRESENTATIVE exige un INDEPENDENT_UNIT_ID."
            )
    return validated_table
