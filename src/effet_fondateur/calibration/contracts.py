"""Contrats statiques de la calibration unifiée R00–R14.

Ce module ne génère aucune donnée. Il charge les contrats de R09, refuse les
clés dupliquées et vérifie que les 25 étapes, 14 vérités et 11 vues arrêtées
dans le plan unifié sont toutes présentes, sans extension implicite.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import jsonschema
import yaml
from yaml.constructor import ConstructorError
from yaml.resolver import BaseResolver


STAGE_IDS = (
    "00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
    "10", "10A", "11", "12", "12A", "13", "14", "15", "16",
    "16A", "16B", "16C", "17", "18", "19",
)

TRUTH_CONSUMERS: dict[str, tuple[str, ...]] = {
    "T-SOURCE": ("01",),
    "T-SUBJECT": ("02", "07", "09"),
    "T-GENOTYPE": ("03", "05", "10", "10A"),
    "T-TARGET": ("04", "09", "10", "12"),
    "T-MEASURE": ("03", "05", "10", "17"),
    "T-KINSHIP": ("06", "07", "08", "09"),
    "T-POPULATION": ("08", "16A"),
    "T-PHASE": ("11", "12", "12A", "13"),
    "T-FOUNDER": ("13", "14", "16B"),
    "T-LD": ("15",),
    "T-ROH": ("16",),
    "T-IBD": ("16C",),
    "T-SCENARIO": ("17",),
    "T-RENDER": ("18", "19"),
}

VIEW_CONSUMERS: dict[str, tuple[str, ...]] = {
    "V-SOURCE": ("01", "02", "03"),
    "V-TARGET": ("04", "09", "10"),
    "V-DENSITY": ("03", "05", "06", "08", "10", "15", "16", "17"),
    "V-MEASUREMENT": ("03", "05", "10", "17"),
    "V-PEDIGREE": ("02", "07", "09", "12", "17"),
    "V-REGIONAL": ("11", "12", "13", "14", "15", "16", "16A", "16B", "16C"),
    "V-AUTOSOMAL": ("05", "06", "07", "08", "10A", "12A", "16C"),
    "V-PHASING-REFERENCE": ("12", "12A", "16A", "16B"),
    "V-FOUNDER-NULL": ("13", "16B"),
    "V-ROH-IBD": ("16", "16C"),
    "V-SENSITIVITY": ("17",),
}

SCHEMA_FILES = (
    "panel_release_manifest.schema.json",
    "subject_registry.schema.json",
    "master_truth_manifest.schema.json",
    "panel_view_manifest.schema.json",
    "blinded_input_inventory.schema.json",
    "pipeline_run_manifest.schema.json",
    "stage_evaluation_manifest.schema.json",
    "readiness_receipt.schema.json",
    "calibration_decision.schema.json",
)


class CalibrationContractError(ValueError):
    """Erreur de contrat avec un code stable utilisable dans un reçu."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


class _UniqueKeyLoader(yaml.SafeLoader):
    """Chargeur YAML sûr qui refuse les clés dupliquées."""


def _construct_unique_mapping(
    loader: _UniqueKeyLoader,
    node: yaml.MappingNode,
    deep: bool = False,
) -> dict[Any, Any]:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise ConstructorError(
                "while constructing a mapping",
                node.start_mark,
                f"duplicate key: {key!r}",
                key_node.start_mark,
            )
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_UniqueKeyLoader.add_constructor(
    BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    document: dict[str, Any] = {}
    for key, value in pairs:
        if key in document:
            raise CalibrationContractError("duplicate_json_key", key)
        document[key] = value
    return document


def load_json(path: Path) -> dict[str, Any]:
    """Charge un objet JSON UTF-8 en refusant les clés répétées."""
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except (OSError, json.JSONDecodeError) as error:
        raise CalibrationContractError("json_contract_invalid", str(path)) from error
    if not isinstance(value, dict):
        raise CalibrationContractError("contract_root_not_object", str(path))
    return value


def load_yaml(path: Path) -> dict[str, Any]:
    """Charge un objet YAML sûr en refusant les clés répétées."""
    try:
        value = yaml.load(path.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader)
    except (OSError, yaml.YAMLError) as error:
        raise CalibrationContractError("yaml_contract_invalid", str(path)) from error
    if not isinstance(value, dict):
        raise CalibrationContractError("contract_root_not_object", str(path))
    return value


def file_sha256(path: Path) -> str:
    """Retourne l'empreinte SHA-256 du contenu exact d'un contrat."""
    return sha256(path.read_bytes()).hexdigest()


def canonical_sha256(document: Any) -> str:
    """Empreinte un objet selon un JSON canonique indépendant de l'indentation."""
    payload = json.dumps(
        document,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def r08_section_sha256(plan_path: Path) -> str:
    """Empreinte la préspecification R08, bornée par ses titres normatifs."""
    text = plan_path.read_text(encoding="utf-8")
    try:
        section = text[text.index("## 8. Préspecification"):text.index("## 9. Contrats")]
    except ValueError as error:
        raise CalibrationContractError("r08_section_missing", str(plan_path)) from error
    return sha256(section.encode("utf-8")).hexdigest()


def validate_schema_set(schema_root: Path) -> dict[str, str]:
    """Valide le socle fermé des neuf schémas R09 et retourne leurs empreintes."""
    actual = {path.name for path in schema_root.glob("*.schema.json")}
    expected = set(SCHEMA_FILES)
    if actual != expected:
        detail = f"missing={sorted(expected - actual)} extra={sorted(actual - expected)}"
        raise CalibrationContractError("schema_set_mismatch", detail)
    hashes: dict[str, str] = {}
    for name in SCHEMA_FILES:
        path = schema_root / name
        schema = load_json(path)
        try:
            jsonschema.Draft202012Validator.check_schema(schema)
        except jsonschema.SchemaError as error:
            raise CalibrationContractError("json_schema_invalid", name) from error
        hashes[name] = file_sha256(path)
    return hashes


def validate_document(document: dict[str, Any], schema_path: Path) -> None:
    """Valide un document contre un schéma R09 autonome."""
    schema = load_json(schema_path)
    validator = jsonschema.Draft202012Validator(
        schema,
        format_checker=jsonschema.Draft202012Validator.FORMAT_CHECKER,
    )
    errors = sorted(
        validator.iter_errors(document),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        first = errors[0]
        location = "/".join(str(part) for part in first.absolute_path) or "$"
        raise CalibrationContractError(
            "document_schema_mismatch",
            f"{schema_path.name}:{location}:{first.message}",
        )


def _normalized_mapping(document: Any, field: str) -> dict[str, tuple[str, ...]]:
    if not isinstance(document, dict):
        raise CalibrationContractError("bundle_mapping_invalid", field)
    normalized: dict[str, tuple[str, ...]] = {}
    for key, values in document.items():
        if not isinstance(key, str) or not isinstance(values, list) or any(
            not isinstance(value, str) for value in values
        ):
            raise CalibrationContractError("bundle_mapping_invalid", field)
        normalized[key] = tuple(values)
    return normalized


def validate_contract_bundle(bundle: dict[str, Any], plan_path: Path) -> None:
    """Vérifie la fermeture du bundle R09 contre les décisions R03–R08."""
    if bundle.get("contract_version") != "1.0.0":
        raise CalibrationContractError("contract_version_mismatch")
    if tuple(bundle.get("stage_ids", ())) != STAGE_IDS:
        raise CalibrationContractError("stage_set_mismatch")
    if _normalized_mapping(bundle.get("truth_consumers"), "truth_consumers") != TRUTH_CONSUMERS:
        raise CalibrationContractError("truth_policy_mismatch")
    if _normalized_mapping(bundle.get("view_consumers"), "view_consumers") != VIEW_CONSUMERS:
        raise CalibrationContractError("view_policy_mismatch")
    if tuple(bundle.get("schema_files", ())) != SCHEMA_FILES:
        raise CalibrationContractError("schema_index_mismatch")
    if bundle.get("r08_prespecification_sha256") != r08_section_sha256(plan_path):
        raise CalibrationContractError("r08_prespecification_hash_mismatch")
    roots = bundle.get("governance_roots")
    if not isinstance(roots, list) or len(roots) != 8 or len(set(roots)) != 8:
        raise CalibrationContractError("governance_roots_mismatch")
    gaps = bundle.get("resolved_gaps")
    expected_gaps = {f"GAP-R04-{index:02d}" for index in range(1, 9)} | {
        "GAP-R06-01", "GAP-R06-02", "GAP-R06-03"
    }
    if not isinstance(gaps, dict) or set(gaps) != expected_gaps:
        raise CalibrationContractError("resolved_gap_set_mismatch")
    if any(not isinstance(value, str) or not value for value in gaps.values()):
        raise CalibrationContractError("resolved_gap_empty")
