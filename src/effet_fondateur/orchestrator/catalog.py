"""Catalogue explicite des étapes que l'orchestrateur est autorisé à lancer."""

from __future__ import annotations

from collections.abc import Iterable

from effet_fondateur.orchestrator.errors import PipelineError
from effet_fondateur.orchestrator.models import StageDefinition


BOOTSTRAP_STAGE_NAMES = {"initialize_run"}


def build_stage_catalog(
    definitions: Iterable[StageDefinition],
) -> dict[str, StageDefinition]:
    """Valide l'unicité et les dépendances, puis indexe les étapes par nom."""
    definitions_tuple = tuple(definitions)
    stage_names = [definition.stage_name for definition in definitions_tuple]
    stage_ids = [definition.stage_id for definition in definitions_tuple]
    if len(stage_names) != len(set(stage_names)):
        raise PipelineError("Le catalogue contient un nom d'étape dupliqué.")
    if len(stage_ids) != len(set(stage_ids)):
        raise PipelineError("Le catalogue contient un identifiant d'étape dupliqué.")

    known_dependencies = set(stage_names) | BOOTSTRAP_STAGE_NAMES
    for definition in definitions_tuple:
        unknown_dependencies = set(definition.dependencies) - known_dependencies
        if unknown_dependencies:
            raise PipelineError(
                f"Dépendance inconnue dans le catalogue pour {definition.stage_name}."
            )
        if definition.stage_name in definition.dependencies:
            raise PipelineError(
                f"Une étape ne peut pas dépendre d'elle-même : {definition.stage_name}."
            )
        if (
            definition.manual_decision_id is not None
            and definition.manual_decision_id
            in definition.resolves_manual_decision_ids
        ):
            raise PipelineError(
                f"Une étape ne peut pas produire et résoudre la même décision : "
                f"{definition.stage_name}."
            )

    # Le catalogue est une liste blanche : une valeur YAML ne peut jamais choisir
    # librement un module Python à importer et exécuter.
    return {definition.stage_name: definition for definition in definitions_tuple}
