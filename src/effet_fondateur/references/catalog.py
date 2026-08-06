"""Résolution hors ligne d'un panel phasé depuis un catalogue versionné."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urljoin, urlparse

from effet_fondateur.audit import sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document


DEFAULT_REFERENCE_CATALOG_PATH = (
    Path(__file__).resolve().parents[3]
    / "config"
    / "references"
    / "1000g_high_coverage_grch38_phased.json"
)
OFFICIAL_REFERENCE_HOST = "ftp.1000genomes.ebi.ac.uk"
REQUIRED_AUTOSOMES = frozenset(range(1, 23))


class ReferenceCatalogError(ValueError):
    """Signale un catalogue ambigu, incomplet ou non approuvé."""


@dataclass(frozen=True)
class ReferenceAsset:
    """Décrit un fichier public attendu et son MD5 publié par le fournisseur."""

    filename: str
    url: str
    official_md5: str


@dataclass(frozen=True)
class ResolvedReferencePanel:
    """Résultat reproductible de la résolution d'un chromosome de référence."""

    catalog_id: str
    catalog_sha256: str
    panel_id: str
    provider: str
    release_id: str
    assembly: str
    chromosome: int
    sample_count: int
    sample_scope: str
    population_scope: str
    readme_url: str
    readme_sha256: str
    manifest_url: str
    manifest_sha256: str
    vcf: ReferenceAsset
    index: ReferenceAsset


def _read_catalog(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ReferenceCatalogError(f"reference_catalog_not_found:{path}") from error
    except json.JSONDecodeError as error:
        raise ReferenceCatalogError(f"reference_catalog_invalid_json:{path}") from error
    if not isinstance(document, dict):
        raise ReferenceCatalogError("reference_catalog_root_not_object")
    try:
        validate_json_document(document, "reference_panel_catalog.schema.json")
    except DocumentValidationError as error:
        raise ReferenceCatalogError("reference_catalog_schema_invalid") from error
    return document


def _validate_official_url(url: str, base_url: str, field: str) -> None:
    parsed_url = urlparse(url)
    parsed_base = urlparse(base_url)
    if (
        parsed_url.scheme != "https"
        or parsed_url.hostname != OFFICIAL_REFERENCE_HOST
        or parsed_base.scheme != "https"
        or parsed_base.hostname != OFFICIAL_REFERENCE_HOST
        or not parsed_url.path.startswith(parsed_base.path.rstrip("/") + "/")
        or parsed_url.query
        or parsed_url.fragment
    ):
        raise ReferenceCatalogError(f"reference_catalog_unapproved_url:{field}")


def _validate_panel(panel: dict[str, Any]) -> dict[int, dict[str, Any]]:
    _validate_official_url(panel["readme_url"], panel["base_url"], "readme_url")
    _validate_official_url(
        panel["manifest_url"], panel["base_url"], "manifest_url"
    )
    chromosome_entries: dict[int, dict[str, Any]] = {}
    for entry in panel["chromosomes"]:
        chromosome = entry["chromosome"]
        if chromosome in chromosome_entries:
            raise ReferenceCatalogError("reference_catalog_duplicate_chromosome")
        chromosome_entries[chromosome] = entry
    if set(chromosome_entries) != REQUIRED_AUTOSOMES:
        raise ReferenceCatalogError("reference_catalog_incomplete_autosomes")
    return chromosome_entries


def load_reference_catalog(
    catalog_path: Path = DEFAULT_REFERENCE_CATALOG_PATH,
) -> dict[str, Any]:
    """Charge le catalogue et applique les contrôles sémantiques de confiance."""
    catalog = _read_catalog(catalog_path)
    panel_ids: set[str] = set()
    for panel in catalog["panels"]:
        if panel["panel_id"] in panel_ids:
            raise ReferenceCatalogError("reference_catalog_duplicate_panel_id")
        panel_ids.add(panel["panel_id"])
        _validate_panel(panel)
    return catalog


def resolve_reference_panel(
    *,
    panel_id: str,
    assembly: str,
    chromosome: int,
    catalog_path: Path = DEFAULT_REFERENCE_CATALOG_PATH,
) -> ResolvedReferencePanel:
    """Résout un autosome sans accès réseau et sans sélection de population."""
    if (
        isinstance(chromosome, bool)
        or not isinstance(chromosome, int)
        or chromosome not in REQUIRED_AUTOSOMES
    ):
        raise ReferenceCatalogError("unsupported_reference_chromosome")
    catalog = load_reference_catalog(catalog_path)
    matches = [panel for panel in catalog["panels"] if panel["panel_id"] == panel_id]
    if len(matches) != 1:
        raise ReferenceCatalogError("reference_panel_not_found_or_ambiguous")
    panel = matches[0]
    if panel["assembly"] != assembly:
        raise ReferenceCatalogError("reference_panel_assembly_mismatch")
    entries = _validate_panel(panel)
    entry = entries[chromosome]
    vcf_filename = panel["vcf_filename_template"].format(chromosome=chromosome)
    index_filename = panel["index_filename_template"].format(chromosome=chromosome)
    base_url = panel["base_url"].rstrip("/") + "/"
    vcf_url = urljoin(base_url, vcf_filename)
    index_url = urljoin(base_url, index_filename)
    _validate_official_url(vcf_url, panel["base_url"], "vcf_url")
    _validate_official_url(index_url, panel["base_url"], "index_url")
    return ResolvedReferencePanel(
        catalog_id=catalog["catalog_id"],
        catalog_sha256=sha256_file(catalog_path),
        panel_id=panel["panel_id"],
        provider=panel["provider"],
        release_id=panel["release_id"],
        assembly=panel["assembly"],
        chromosome=chromosome,
        sample_count=panel["sample_count"],
        sample_scope=panel["sample_scope"],
        population_scope=panel["population_scope"],
        readme_url=panel["readme_url"],
        readme_sha256=panel["readme_sha256"],
        manifest_url=panel["manifest_url"],
        manifest_sha256=panel["manifest_sha256"],
        vcf=ReferenceAsset(vcf_filename, vcf_url, entry["vcf_md5"]),
        index=ReferenceAsset(index_filename, index_url, entry["index_md5"]),
    )
