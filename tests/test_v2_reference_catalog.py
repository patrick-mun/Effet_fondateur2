import json
from pathlib import Path

import pytest

from effet_fondateur.references import (
    DEFAULT_REFERENCE_CATALOG_PATH,
    ReferenceCatalogError,
    load_reference_catalog,
    resolve_reference_panel,
)


PANEL_ID = "1kg_3202_high_coverage_20220422"


def write_modified_catalog(tmp_path: Path, mutation) -> Path:
    catalog = json.loads(DEFAULT_REFERENCE_CATALOG_PATH.read_text(encoding="utf-8"))
    mutation(catalog)
    catalog_path = tmp_path / "catalog.json"
    catalog_path.write_text(
        json.dumps(catalog, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return catalog_path


def test_catalog_covers_all_autosomes_with_all_release_samples() -> None:
    catalog = load_reference_catalog()

    assert catalog["schema_version"] == "1.0.0"
    assert len(catalog["panels"]) == 1
    panel = catalog["panels"][0]
    assert panel["sample_count"] == 3202
    assert panel["sample_scope"] == "ALL_RELEASE_SAMPLES"
    assert panel["population_scope"] == "ALL_AVAILABLE_1KGP_POPULATIONS"
    assert {entry["chromosome"] for entry in panel["chromosomes"]} == set(
        range(1, 23)
    )


def test_resolver_returns_official_chr19_assets_without_network() -> None:
    resolved = resolve_reference_panel(
        panel_id=PANEL_ID,
        assembly="GRCh38",
        chromosome=19,
    )

    assert resolved.release_id == "20220422_3202_phased_SNV_INDEL_SV"
    assert resolved.provider == "IGSR_1000_GENOMES"
    assert resolved.sample_count == 3202
    assert resolved.population_scope == "ALL_AVAILABLE_1KGP_POPULATIONS"
    assert resolved.vcf.filename.startswith("1kGP_high_coverage_Illumina.chr19.")
    assert resolved.vcf.url == (
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000G_2504_high_coverage/working/"
        "20220422_3202_phased_SNV_INDEL_SV/"
        "1kGP_high_coverage_Illumina.chr19.filtered."
        "SNV_INDEL_SV_phased_panel.vcf.gz"
    )
    assert resolved.vcf.official_md5 == "cf5c060fbcc0fcd2fa089f27fd6ac78e"
    assert resolved.index.official_md5 == "f70b288343a7d5cc40124b672fd2727f"
    assert resolved.readme_url.endswith("README_1kGP_phased_panel_110722.pdf")
    assert resolved.readme_sha256 == (
        "bbc6b3fe5d30b9b9a999ec600d88d15e9135078ca56c2e1dcc6d12d1f44dda13"
    )
    assert resolved.manifest_url.endswith("20220804_manifest.txt")
    assert resolved.manifest_sha256 == (
        "71366add1f73a5a3283b8d89d5253f0b5b6e963032eb4a7e1862c19228c4579d"
    )
    assert len(resolved.catalog_sha256) == 64


def test_resolver_rejects_assembly_mismatch() -> None:
    with pytest.raises(ReferenceCatalogError, match="assembly_mismatch"):
        resolve_reference_panel(
            panel_id=PANEL_ID,
            assembly="GRCh37",
            chromosome=19,
        )


def test_resolver_rejects_sex_chromosomes() -> None:
    with pytest.raises(ReferenceCatalogError, match="unsupported_reference_chromosome"):
        resolve_reference_panel(
            panel_id=PANEL_ID,
            assembly="GRCh38",
            chromosome=23,
        )


def test_resolver_rejects_noninteger_chromosome() -> None:
    with pytest.raises(ReferenceCatalogError, match="unsupported_reference_chromosome"):
        resolve_reference_panel(
            panel_id=PANEL_ID,
            assembly="GRCh38",
            chromosome=19.0,
        )


def test_catalog_rejects_missing_autosome(tmp_path: Path) -> None:
    catalog_path = write_modified_catalog(
        tmp_path,
        lambda catalog: catalog["panels"][0]["chromosomes"].pop(),
    )

    with pytest.raises(ReferenceCatalogError, match="schema_invalid"):
        load_reference_catalog(catalog_path)


def test_catalog_rejects_duplicate_autosome(tmp_path: Path) -> None:
    def duplicate_chromosome(catalog: dict) -> None:
        catalog["panels"][0]["chromosomes"][-1]["chromosome"] = 21

    catalog_path = write_modified_catalog(tmp_path, duplicate_chromosome)

    with pytest.raises(ReferenceCatalogError, match="duplicate_chromosome"):
        load_reference_catalog(catalog_path)


def test_catalog_rejects_nonofficial_host(tmp_path: Path) -> None:
    def replace_host(catalog: dict) -> None:
        panel = catalog["panels"][0]
        panel["base_url"] = panel["base_url"].replace(
            "ftp.1000genomes.ebi.ac.uk", "example.org"
        )

    catalog_path = write_modified_catalog(tmp_path, replace_host)

    with pytest.raises(ReferenceCatalogError, match="unapproved_url"):
        load_reference_catalog(catalog_path)


def test_catalog_rejects_embedded_url_credentials(tmp_path: Path) -> None:
    def add_credentials(catalog: dict) -> None:
        panel = catalog["panels"][0]
        panel["base_url"] = panel["base_url"].replace(
            "https://", "https://user:password@"
        )

    catalog_path = write_modified_catalog(tmp_path, add_credentials)

    with pytest.raises(ReferenceCatalogError, match="unapproved_url"):
        load_reference_catalog(catalog_path)
