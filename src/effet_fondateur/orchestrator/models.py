"""Modèles immuables utilisés par l'orchestrateur."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StageDefinition:
    """Décrit un script d'étape exécutable sans logique scientifique embarquée."""

    stage_id: str
    stage_name: str
    module: str
    critical: bool = True
    dependencies: tuple[str, ...] = ()

    @property
    def directory_name(self) -> str:
        """Retourne le nom stable du dossier publié de l'étape."""
        return f"{self.stage_id}_{self.stage_name}"
