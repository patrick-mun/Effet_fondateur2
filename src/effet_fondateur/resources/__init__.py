"""Accès aux ressources immuables distribuées avec le paquet."""

from __future__ import annotations

from importlib.resources import files
from importlib.resources.abc import Traversable


def package_resource(*relative_parts: str) -> Traversable:
    """Retourne une ressource embarquée sans dépendre de la racine du dépôt."""
    return files("effet_fondateur").joinpath(*relative_parts)


SCHEMAS_DIR = package_resource("resources", "schemas")
REFERENCES_DIR = package_resource("resources", "references")
REPORTING_ASSETS_DIR = package_resource("reporting", "assets")
