"""Tests des contrats JSON D13."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema


ROOT = Path(__file__).resolve().parents[2]
CONFIG_SCHEMA = ROOT / "calibration/schemas/pedigree_config.schema.json"
MANIFEST_SCHEMA = ROOT / "calibration/schemas/pedigree_manifest.schema.json"
SHA = "a" * 64


def _config() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "pedigree_config_id": "d13.test",
        "config_state": "FINALIZED",
        "method": {
            "method_id": "first_cousin_whole_haplotype_transmission_v1",
            "method_version": "1.0.0",
            "random_seed": 1003,
            "family_template": "FIRST_COUSIN_PEDIGREE",
            "source_selection": "SAME_STRATUM_WITHOUT_REPLACEMENT",
            "recombination_policy": "FORBIDDEN_DEFERRED_TO_D14",
        },
        "family_groups": [
            {
                "family_group_id": "PED_A",
                "stratum_id": "POP_A",
                "family_count": 1,
            }
        ],
        "unrelated_pair_count": 1,
    }


def _artifact(role: str, suffix: str, count: int) -> dict[str, object]:
    return {
        "artifact_id": f"d13.{suffix}",
        "role": role,
        "path": f"private_truth/d13/{suffix}",
        "sha256": SHA,
        "size_bytes": 10,
        "record_count": count,
    }


def _manifest() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "pedigree_manifest_id": "d13.manifest",
        "manifest_state": "FINALIZED",
        "created_at": "2026-10-29T10:00:00Z",
        "input": {
            "manifest_id": "d12.manifest",
            "manifest_sha256": SHA,
            "panel_role": "CALIBRATION",
            "dense_artifact_id": "d12.dense",
            "dense_artifact_sha256": SHA,
            "observed_artifact_id": "d12.observed",
            "observed_artifact_sha256": SHA,
            "ancestry_artifact_id": "d12.ancestry",
            "ancestry_artifact_sha256": SHA,
            "provenance_artifact_id": "d12.provenance",
            "provenance_artifact_sha256": SHA,
        },
        "configuration": {
            "configuration_id": "d13.test",
            "configuration_sha256": SHA,
            "random_seed": 1003,
        },
        "panel": {
            "panel_id": "d13.calibration",
            "role": "CALIBRATION",
            "synthetic_identifier_prefix": "SYN_CAL_D13_",
            "source_individual_count": 8,
            "output_individual_count": 12,
            "descendant_count": 4,
            "family_count": 1,
            "trio_count": 4,
            "unrelated_pair_count": 1,
        },
        "domain": {
            "chromosome": 19,
            "start_bp": 100,
            "end_bp": 200,
            "position_unit": "BP_1_BASED_CLOSED",
        },
        "method": {
            "method_id": "first_cousin_whole_haplotype_transmission_v1",
            "method_version": "1.0.0",
            "family_template": "FIRST_COUSIN_PEDIGREE",
            "source_selection": "SAME_STRATUM_WITHOUT_REPLACEMENT",
            "transmission_unit": "WHOLE_PARENTAL_HAPLOTYPE",
            "kinship_definition": "PEDIGREE_COANCESTRY_COEFFICIENT",
            "recombination_policy": "FORBIDDEN_DEFERRED_TO_D14",
        },
        "confidentiality": {
            "real_sample_files_read": False,
            "d12_identifiers_in_blinded_panel": False,
            "pedigree_labels_in_blinded_panel": False,
            "ancestry_labels_in_blinded_panel": False,
            "local_source_paths_in_manifest": False,
            "artifacts_version_control_allowed": False,
        },
        "artifacts": [
            _artifact("DENSE_PHASED_PEDIGREE_PANEL", "dense.vcf.gz", 2),
            _artifact("BLINDED_OBSERVED_PEDIGREE_PANEL", "observed.vcf.gz", 1),
            _artifact("PEDIGREE_TRUTH", "pedigree.tsv.gz", 12),
            _artifact("KINSHIP_TRUTH", "kinship.tsv.gz", 29),
            _artifact("MENDELIAN_TRANSMISSION_TRUTH", "transmission.tsv.gz", 8),
            _artifact("D12_HAPLOTYPE_PROVENANCE", "provenance.tsv.gz", 24),
            _artifact("ANCESTRY_TRUTH", "ancestry.tsv.gz", 12),
        ],
        "integrity": {
            "dense_variant_count": 2,
            "observed_variant_count": 1,
            "transmission_record_count": 8,
            "kinship_record_count": 29,
            "relationship_class_counts": {
                "UNRELATED": 13,
                "FIRST_DEGREE": 9,
                "SECOND_DEGREE": 6,
                "THIRD_DEGREE": 1,
            },
            "founder_haplotype_copy_error_count": 0,
            "founder_allele_count_difference": 0,
            "mendelian_error_count": 0,
            "family_founder_overlap_count": 0,
            "invalid_parent_reference_count": 0,
            "checks": [
                "D12_BINDING_MATCHED",
                "SOURCE_SAMPLE_UNIVERSE_MATCHED",
                "FOUNDER_HAPLOTYPE_BIJECTION",
                "FAMILY_FOUNDERS_DISJOINT",
                "MENDELIAN_TRANSMISSIONS_EXACT",
                "TRIO_COUNT_MATCHED",
                "KINSHIP_CLASSES_COVERED",
                "ANCESTRY_PROPAGATED",
                "BLINDED_PANEL_HAS_NO_TRUTH_LABEL",
            ],
        },
    }


def test_d13_config_and_manifest_schemas_accept_valid_documents() -> None:
    for path, instance in ((CONFIG_SCHEMA, _config()), (MANIFEST_SCHEMA, _manifest())):
        schema = json.loads(path.read_text(encoding="utf-8"))
        jsonschema.Draft202012Validator.check_schema(schema)
        jsonschema.Draft202012Validator(
            schema, format_checker=jsonschema.FormatChecker()
        ).validate(instance)
