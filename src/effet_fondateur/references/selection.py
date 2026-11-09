"""Sélection déclarée de référence, sans accès à la vérité de simulation.

La préparation D23a vérifie D10 et ses groupes d'indépendance. Le pipeline
consomme seulement sa liste de références autorisées, liée à la signature
d'entrée ; il ne lit ni donneurs, ni transmissions, ni phase vraie.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any, Sequence

from effet_fondateur.audit import read_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document


@dataclass(frozen=True)
class ReferenceSampleSelection:
    """Liste ordonnée privée ; seules les empreintes sont publiables en audit."""

    sample_ids: tuple[str, ...]
    source_panel_id: str
    selection_sha256: str
    membership_sha256: str
    partition_manifest_sha256: str

    def audit_summary(self) -> dict[str, Any]:
        """Retourne la provenance agrégée, sans identifiant de référence."""
        return dict(selection_sha256=self.selection_sha256,
                    membership_sha256=self.membership_sha256,
                    partition_manifest_sha256=self.partition_manifest_sha256,
                    sample_count=len(self.sample_ids))

    def checked_samples(self, available: Sequence[str]) -> list[str]:
        """Refuse perte, doublon ou nom inconnu au lieu de réduire la liste."""
        if (not self.sample_ids or len(self.sample_ids) != len({sample.casefold() for sample in self.sample_ids})
                or len(available) != len(set(available)) or not set(self.sample_ids) <= set(available)):
            raise ValueError("phasing_reference_selection_not_in_source")
        return list(self.sample_ids)


def load_reference_selection(path: Path, *, expected_sha256: str, panel_id: str) -> ReferenceSampleSelection:
    """Valide fichier, schéma, effectif, ordre canonique et engagement SHA-256."""
    if any(item.is_symlink() for item in (path, *path.parents)) or not path.is_file() or sha256_file(path) != expected_sha256:
        raise ValueError("phasing_reference_selection_integrity_invalid")
    document = read_json(path)
    try:
        validate_json_document(document, "phasing_reference_selection.schema.json")
    except DocumentValidationError:
        # uniqueItems peut inclure toute la liste dans le message JSONSchema.
        raise ValueError("phasing_reference_selection_schema_invalid") from None
    samples = document["sample_ids"]
    payload = "".join(sample + "\n" for sample in samples).encode()
    if (document["source_panel_id"] != panel_id or document["sample_count"] != len(samples)
            or samples != sorted(samples, key=str.casefold)
            or len({sample.casefold() for sample in samples}) != len(samples)
            or hashlib.sha256(payload).hexdigest() != document["membership_sha256"]):
        raise ValueError("phasing_reference_selection_content_invalid")
    return ReferenceSampleSelection(tuple(samples), panel_id, expected_sha256,
                                    document["membership_sha256"], document["partition_manifest_sha256"])


def stage_reference_selection(stage_inputs: dict[str, Any], config: dict[str, Any], run_dir: Path,
                              panel_id: str) -> ReferenceSampleSelection | None:
    """Exige une entrée signée si la sélection est configurée ou obligatoire.

    L'absence de sélection conserve le comportement historique seulement si
    aucun chemin n'est configuré et si le garde-fou explicite est désactivé.
    """
    required = stage_inputs["parameters"].get("require_reference_sample_selection", False)
    if not isinstance(required, bool):
        raise ValueError("require_reference_sample_selection_invalid")
    configured = config["inputs"].get("phasing_reference_selection")
    matches = [item for item in stage_inputs["artifacts"] if item["artifact_id"] == "config_input_phasing_reference_selection"]
    if not matches and configured is None and not required:
        return None
    if len(matches) != 1 or configured is None:
        raise ValueError("phasing_reference_selection_missing_or_ambiguous")
    artifact = matches[0]
    path = Path(artifact["path"])
    if not path.is_absolute():
        path = run_dir / path if path.parts[:1] == ("stages",) else Path.cwd() / path
    configured_path = Path(configured)
    if not configured_path.is_absolute():
        configured_path = Path.cwd() / configured_path
    if path.resolve() != configured_path.resolve():
        raise ValueError("phasing_reference_selection_config_binding_invalid")
    return load_reference_selection(path, expected_sha256=artifact["sha256"], panel_id=panel_id)
