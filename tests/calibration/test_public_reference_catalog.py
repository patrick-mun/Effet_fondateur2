"""Tests du catalogue D09 des références publiques de calibration."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "references"
    / "public_reference_catalog.json"
)
SCHEMA_PATH = (
    REPOSITORY_ROOT
    / "calibration"
    / "schemas"
    / "public_reference_catalog.schema.json"
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Refuse les clés JSON dupliquées dans les contrats D09."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Clé JSON dupliquée dans le contrat D09 : {key}")
        result[key] = value
    return result


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_keys,
    )


def _catalog() -> dict[str, Any]:
    return _load_json(CATALOG_PATH)


def _validator() -> jsonschema.Draft202012Validator:
    schema = _load_json(SCHEMA_PATH)
    jsonschema.Draft202012Validator.check_schema(schema)
    return jsonschema.Draft202012Validator(
        schema,
        format_checker=jsonschema.FormatChecker(),
    )


def _errors(catalog: dict[str, Any]) -> list[jsonschema.ValidationError]:
    return sorted(
        _validator().iter_errors(catalog),
        key=lambda error: list(error.absolute_path),
    )


def _source_by_type(catalog: dict[str, Any], source_type: str) -> dict[str, Any]:
    return next(
        source for source in catalog["sources"] if source["source_type"] == source_type
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_schema_and_frozen_catalog_are_valid() -> None:
    assert _errors(_catalog()) == []


def test_catalog_has_unique_source_identifiers_and_types() -> None:
    catalog = _catalog()
    source_ids = [source["source_id"] for source in catalog["sources"]]
    source_types = [source["source_type"] for source in catalog["sources"]]

    assert len(source_ids) == len(set(source_ids))
    assert len(source_types) == len(set(source_types)) == 5


@pytest.mark.parametrize("missing_field", ["source_id", "version", "origin", "fingerprint"])
def test_identity_version_origin_and_sha256_are_required(missing_field: str) -> None:
    catalog = _catalog()
    del catalog["sources"][0][missing_field]

    assert _errors(catalog)


@pytest.mark.parametrize(
    "invalid_sha256",
    ["", "abc", "A" * 64, "0" * 63, "g" * 64],
)
def test_invalid_source_sha256_is_rejected(invalid_sha256: str) -> None:
    catalog = _catalog()
    catalog["sources"][0]["fingerprint"]["sha256"] = invalid_sha256

    assert _errors(catalog)


@pytest.mark.parametrize(
    "invalid_url",
    [
        "http://example.org/reference.vcf.gz",
        "https://example.org/reference.vcf.gz",
        "file:///tmp/reference.vcf.gz",
        "../reference.vcf.gz",
        "/tmp/reference.vcf.gz",
    ],
)
def test_non_https_or_local_origins_are_rejected(invalid_url: str) -> None:
    catalog = _catalog()
    catalog["sources"][0]["origin"]["canonical_url"] = invalid_url

    assert _errors(catalog)


def test_every_source_is_limited_to_local_calibration_without_redistribution() -> None:
    for source in _catalog()["sources"]:
        assert source["rights"]["use_decision"] == "ALLOWED_LOCAL_CALIBRATION"
        assert source["rights"]["redistribution_allowed_by_project"] is False
        assert source["acquisition"]["download_during_d09"] is False
        assert source["acquisition"]["local_sha256_required_before_use"] is True


def test_unreviewed_array_design_source_is_rejected() -> None:
    catalog = _catalog()
    source = _source_by_type(catalog, "ARRAY_DESIGN_MANIFEST")
    source["source_id"] = "unreviewed_source"

    assert _errors(catalog)


def test_reviewed_cytoscan_source_is_bound_to_the_acquired_package() -> None:
    source = _source_by_type(_catalog(), "ARRAY_DESIGN_MANIFEST")

    assert source["source_id"] == "thermofisher_cytoscan750k_accel_na36_r1"
    assert source["version"] == "NA36.r1"
    assert source["assembly"]["name"] == "GRCh38"
    assert source["fingerprint"]["sha256"] == (
        "6a696ef3d2e61e9db74ead91c6410e15efc90da6a8bfc814f5da034e364607a3"
    )
    assert source["rights"]["redistribution_allowed_by_project"] is False


def test_missing_provider_checksum_cannot_claim_an_algorithm() -> None:
    catalog = _catalog()
    source = _source_by_type(catalog, "GENETIC_MAP")
    source["acquisition"]["provider_checksum_algorithm"] = "SHA-256"
    source["acquisition"]["provider_checksum_required"] = False

    assert _errors(catalog)


def test_assembly_is_explicit_and_verified_for_every_source() -> None:
    for source in _catalog()["sources"]:
        assembly = source["assembly"]
        assert assembly["verified"] is True
        assert assembly["name"] in {"GRCh38", "ASSEMBLY_INDEPENDENT"}
        assert assembly["associated_panel_assembly"] == "GRCh38"


def test_repository_catalog_bindings_match_current_files() -> None:
    for source in _catalog()["sources"]:
        binding = source["catalog_binding"]
        bound_path = REPOSITORY_ROOT / binding["repository_path"]

        assert bound_path.is_file()
        assert _sha256(bound_path) == binding["sha256"]


def test_phased_panel_matches_the_existing_pipeline_catalog() -> None:
    source = _source_by_type(_catalog(), "PHASED_HAPLOTYPE_PANEL")
    pipeline_catalog = _load_json(
        REPOSITORY_ROOT
        / "config"
        / "references"
        / "1000g_high_coverage_grch38_phased.json"
    )
    panel = pipeline_catalog["panels"][0]

    assert source["version"] == panel["release_id"]
    assert source["assembly"]["name"] == panel["assembly"]
    assert source["origin"]["canonical_url"] == panel["base_url"]
    assert source["fingerprint"]["asset_url"] == panel["manifest_url"]
    assert source["fingerprint"]["sha256"] == panel["manifest_sha256"]
    assert source["fingerprint"]["scope"] == "PROVIDER_MANIFEST"
    assert source["acquisition"]["provider_checksum_algorithm"] == "MD5"
    assert len(panel["chromosomes"]) == 22


def test_population_metadata_and_unrelated_index_match_pipeline_catalog() -> None:
    catalog = _catalog()
    pipeline_catalog = _load_json(
        REPOSITORY_ROOT
        / "config"
        / "references"
        / "1000g_ancestry_grch38.json"
    )
    expected_by_type = {
        "POPULATION_METADATA": pipeline_catalog["population_metadata"],
        "UNRELATED_SAMPLE_INDEX": pipeline_catalog["unrelated_sample_index"],
    }

    for source_type, expected in expected_by_type.items():
        source = _source_by_type(catalog, source_type)
        assert source["origin"]["canonical_url"] == expected["url"]
        assert source["fingerprint"]["asset_url"] == expected["url"]
        assert source["fingerprint"]["sha256"] == expected["sha256"]


def test_genetic_map_matches_the_existing_pipeline_catalog() -> None:
    source = _source_by_type(_catalog(), "GENETIC_MAP")
    pipeline_catalog = _load_json(
        REPOSITORY_ROOT
        / "config"
        / "references"
        / "shapeit4_hapmap_grch38.json"
    )
    genetic_map = pipeline_catalog["maps"][0]

    assert source["version"] == genetic_map["release_id"]
    assert source["assembly"]["name"] == genetic_map["assembly"]
    assert source["origin"]["canonical_url"] == genetic_map["archive_url"]
    assert source["fingerprint"]["asset_url"] == genetic_map["archive_url"]
    assert source["fingerprint"]["sha256"] == genetic_map["archive_sha256"]
    assert source["rights"]["review_status"] == (
        "REPOSITORY_TERMS_VERIFIED_SOURCE_RIGHTS_PARTIAL"
    )
    assert "HapMap" in source["rights"]["notes"]
    assert genetic_map["chromosomes"] == list(range(1, 23))


def test_catalog_rejects_an_attempt_to_enable_real_study_data() -> None:
    catalog = deepcopy(_catalog())
    catalog["scope"]["real_study_data_allowed"] = True

    assert _errors(catalog)


def test_duplicate_json_keys_are_rejected() -> None:
    with pytest.raises(ValueError, match="Clé JSON dupliquée"):
        json.loads(
            '{"schema_version": "1.0.0", "schema_version": "2.0.0"}',
            object_pairs_hook=_reject_duplicate_keys,
        )
