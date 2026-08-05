"""Entrées-sorties durables utilisées par l'orchestrateur V2."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any
from uuid import uuid4


def atomic_write_json(path: Path, document: dict[str, Any]) -> None:
    """Écrit un document JSON par remplacement atomique dans son dossier cible."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with temporary_path.open("w", encoding="utf-8") as output_file:
            json.dump(document, output_file, ensure_ascii=False, indent=2, sort_keys=True)
            output_file.write("\n")
            output_file.flush()
            os.fsync(output_file.fileno())
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def append_json_event(path: Path, event: dict[str, Any]) -> None:
    """Ajoute un événement JSON sur une ligne et force son écriture sur disque."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as output_file:
        output_file.write(json.dumps(event, ensure_ascii=False, sort_keys=True))
        output_file.write("\n")
        output_file.flush()
        os.fsync(output_file.fileno())


def read_json(path: Path) -> dict[str, Any]:
    """Lit un objet JSON et refuse toute racine qui n'est pas un objet."""
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ValueError(f"Document JSON introuvable : {path}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"Document JSON invalide : {path}") from error
    if not isinstance(document, dict):
        raise ValueError(f"Le document JSON doit être un objet : {path}")
    return document
