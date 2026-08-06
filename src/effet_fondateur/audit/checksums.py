"""Calcul d'empreintes reproductibles pour les artefacts du pipeline."""

from __future__ import annotations

import hashlib
from pathlib import Path


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Calcule le MD5 fournisseur d'un fichier public sans usage cryptographique."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as input_file:
        while chunk := input_file.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Calcule l'empreinte SHA-256 d'un fichier sans le charger en mémoire."""
    digest = hashlib.sha256()
    with path.open("rb") as input_file:
        while chunk := input_file.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()
