"""Tests du contrat privé D11 des haplotypes de base."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import jsonschema
import pytest


ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "calibration/schemas/base_haplotype_manifest.schema.json"
SHA = "a" * 64


def _schema() -> dict[str, object]:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def _manifest() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "base_haplotype_manifest_id": "d11.calibration.chr19",
        "manifest_state": "FINALIZED",
        "created_at": "2026-10-28T10:00:00Z",
        "panel": {
            "panel_id": "calibration_base_chr19",
            "role": "CALIBRATION",
            "synthetic_identifier_prefix": "SYN_CAL_",
            "individual_count": 4,
            "haplotype_count": 8,
        },
        "source": {
            "source_id": "igsr_1kg_3202_phased_grch38_20220422",
            "release_id": "20220422_3202_phased_SNV_INDEL_SV",
            "assembly": "GRCh38",
            "chromosome": 19,
            "vcf_filename": "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
            "vcf_sha256": SHA,
            "index_sha256": SHA,
            "genetic_map_source_id": "shapeit4_hapmap_grch38_git_7a0cab7",
            "genetic_map_sha256": SHA,
        },
        "partition": {
            "partition_manifest_id": "d10.partition.manifest",
            "partition_manifest_sha256": SHA,
            "partition_role": "DONOR_CALIBRATION",
            "membership_sha256": SHA,
            "donor_count": 4,
        },
        "method": {
            "method_id": "population_stratified_whole_haplotype_repairing_v1",
            "method_version": "1.0.0",
            "random_seed": 1002,
            "pairing_unit": "PUBLIC_SOURCE_POPULATION",
            "haplotype_usage": "EACH_DONOR_HAPLOTYPE_EXACTLY_ONCE",
            "same_source_pair_policy": "FORBIDDEN",
            "variant_policy": "BIALLELIC_SNV_COMPLETE_PHASED",
            "genotype_derivation": "UNORDERED_PAIR_OF_TRUE_HAPLOTYPES",
        },
        "domain": {
            "chromosome": 19,
            "start_bp": 100,
            "end_bp": 300,
            "base_position_unit": "BP_1_BASED_CLOSED",
            "genetic_distance_unit": "cM",
            "genetic_map_interpolation": "PIECEWISE_LINEAR_NO_EXTRAPOLATION",
        },
        "confidentiality": {
            "raw_donor_identifiers_in_manifest": False,
            "local_source_paths_in_manifest": False,
            "donor_provenance_private": True,
            "artifacts_version_control_allowed": False,
            "aggregate_counts_publishable": True,
        },
        "artifacts": [
            {
                "artifact_id": "d11.phased",
                "role": "PHASED_HAPLOTYPES",
                "path": "private_truth/d11/base_haplotypes.phased.vcf.gz",
                "sha256": SHA,
                "size_bytes": 100,
                "record_count": 3,
            },
            {
                "artifact_id": "d11.unphased",
                "role": "UNPHASED_GENOTYPES",
                "path": "private_truth/d11/base_genotypes.unphased.vcf.gz",
                "sha256": SHA,
                "size_bytes": 100,
                "record_count": 3,
            },
            {
                "artifact_id": "d11.provenance",
                "role": "DONOR_PROVENANCE",
                "path": "private_truth/d11/donor_provenance.tsv.gz",
                "sha256": SHA,
                "size_bytes": 100,
                "record_count": 8,
            },
            {
                "artifact_id": "d11.map",
                "role": "GENETIC_MAP",
                "path": "private_truth/d11/genetic_map.tsv",
                "sha256": SHA,
                "size_bytes": 100,
                "record_count": 3,
            },
        ],
        "integrity": {
            "source_donor_count": 4,
            "synthetic_individual_count": 4,
            "source_haplotype_count": 8,
            "synthetic_haplotype_count": 8,
            "variant_count": 3,
            "population_stratum_count": 2,
            "haplotype_copy_error_count": 0,
            "same_source_pair_count": 0,
            "allele_count_difference": 0,
            "genotype_haplotype_discordance_count": 0,
            "genetic_map_coordinate_error_count": 0,
            "genetic_map_domain_covered": True,
            "mendelian_evaluation": "DEFERRED_TO_D13_NO_PEDIGREE",
            "checks": [
                "SOURCE_HASHES_MATCHED",
                "PARTITION_BINDING_MATCHED",
                "HAPLOTYPE_COPY_BIJECTION",
                "SOURCE_PAIRS_DISJOINT",
                "ALLELE_FREQUENCIES_EXACT",
                "HAPLOTYPE_LD_EXACT_BY_PERMUTATION",
                "GENOTYPE_HAPLOTYPE_CONCORDANCE_100_PERCENT",
                "GENETIC_MAP_COORDINATES_VALID",
                "RAW_DONOR_IDENTIFIERS_NOT_IN_MANIFEST",
            ],
        },
    }


def _validate(manifest: dict[str, object]) -> None:
    jsonschema.Draft202012Validator(_schema()).validate(manifest)


def test_schema_is_valid_draft_2020_12() -> None:
    jsonschema.Draft202012Validator.check_schema(_schema())


def test_valid_manifest_is_accepted() -> None:
    _validate(_manifest())


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("manifest_state",), "DRAFT"),
        (("panel", "synthetic_identifier_prefix"), "DONOR_001"),
        (("source", "assembly"), "GRCh37"),
        (("partition", "partition_role"), "PHASING_REFERENCE_SHARED"),
        (("method", "haplotype_usage"), "WITH_REPLACEMENT"),
        (("method", "same_source_pair_policy"), "ALLOWED"),
        (("domain", "genetic_distance_unit"), "Morgans"),
        (("confidentiality", "raw_donor_identifiers_in_manifest"), True),
        (("integrity", "same_source_pair_count"), 1),
        (("integrity", "allele_count_difference"), 1),
        (("integrity", "mendelian_evaluation"), "PASS"),
    ],
)
def test_fixed_contract_violations_are_rejected(
    path: tuple[str, ...], value: object
) -> None:
    manifest = deepcopy(_manifest())
    target: dict[str, object] = manifest
    for key in path[:-1]:
        target = target[key]  # type: ignore[assignment]
    target[path[-1]] = value
    with pytest.raises(jsonschema.ValidationError):
        _validate(manifest)


def test_absolute_local_artifact_path_is_rejected() -> None:
    manifest = _manifest()
    manifest["artifacts"][0]["path"] = "/Users/example/Desktop/reference.vcf.gz"  # type: ignore[index]
    with pytest.raises(jsonschema.ValidationError):
        _validate(manifest)


def test_artifact_outside_private_d11_is_rejected() -> None:
    manifest = _manifest()
    manifest["artifacts"][0]["path"] = "blinded_inputs/base.vcf.gz"  # type: ignore[index]
    with pytest.raises(jsonschema.ValidationError):
        _validate(manifest)


def test_each_artifact_role_is_required_exactly_once() -> None:
    manifest = _manifest()
    manifest["artifacts"][1]["role"] = "PHASED_HAPLOTYPES"  # type: ignore[index]
    with pytest.raises(jsonschema.ValidationError):
        _validate(manifest)


def test_unknown_property_is_rejected() -> None:
    manifest = _manifest()
    manifest["source"]["local_path"] = "/tmp/source.vcf.gz"  # type: ignore[index]
    with pytest.raises(jsonschema.ValidationError):
        _validate(manifest)


def test_validation_prefix_must_match_role() -> None:
    manifest = _manifest()
    manifest["panel"]["role"] = "VALIDATION"  # type: ignore[index]
    with pytest.raises(jsonschema.ValidationError):
        _validate(manifest)
