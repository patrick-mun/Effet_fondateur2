"""Validation des contrats JSON D12."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import jsonschema
import pytest


ROOT = Path(__file__).resolve().parents[2]
CONFIG_SCHEMA = ROOT / "calibration/schemas/population_structure_config.schema.json"
MANIFEST_SCHEMA = ROOT / "calibration/schemas/population_structure_manifest.schema.json"
SHA = "a" * 64


def _config() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "population_structure_config_id": "d12.micro.v1",
        "config_state": "FINALIZED",
        "method": {
            "method_id": "reciprocal_whole_haplotype_mixture_v1",
            "method_version": "1.0.0",
            "random_seed": 1002,
            "pure_policy": "RETAIN_UNUSED_D11_DIPLOIDS",
            "mixture_unit": "RECIPROCAL_WHOLE_HAPLOTYPE_PAIR",
            "supported_proportions": "PURE_OR_TWO_WAY_EQUAL",
        },
        "mixtures": [
            {
                "mixture_id": "ADMIX_A_B",
                "individual_count": 2,
                "components": [
                    {"stratum_id": "POP_A", "proportion": 0.5},
                    {"stratum_id": "POP_B", "proportion": 0.5},
                ],
            }
        ],
    }


def _artifact(role: str, name: str, count: int) -> dict[str, object]:
    return {
        "artifact_id": f"d12.{name}",
        "role": role,
        "path": f"private_truth/d12/{name}",
        "sha256": SHA,
        "size_bytes": 10,
        "record_count": count,
    }


def _manifest() -> dict[str, object]:
    return {
        "schema_version": "1.0.0",
        "population_structure_manifest_id": "d12.manifest",
        "manifest_state": "FINALIZED",
        "created_at": "2026-10-29T08:00:00Z",
        "inputs": {
            "d11a": {"manifest_id": "d11a.manifest", "manifest_sha256": SHA},
            "d11b": {"manifest_id": "d11b.manifest", "manifest_sha256": SHA},
            "panel_role": "CALIBRATION",
        },
        "configuration": {
            "configuration_id": "d12.micro.v1",
            "configuration_sha256": SHA,
            "random_seed": 1002,
        },
        "panel": {
            "panel_id": "d12.calibration",
            "role": "CALIBRATION",
            "synthetic_identifier_prefix": "SYN_CAL_D12_",
            "individual_count": 4,
            "haplotype_count": 8,
            "pure_individual_count": 2,
            "admixed_individual_count": 2,
            "source_stratum_count": 2,
            "output_stratum_count": 3,
            "mixture_profile_count": 1,
        },
        "domain": {
            "chromosome": 19,
            "start_bp": 100,
            "end_bp": 200,
            "position_unit": "BP_1_BASED_CLOSED",
        },
        "method": {
            "method_id": "reciprocal_whole_haplotype_mixture_v1",
            "method_version": "1.0.0",
            "mixture_unit": "RECIPROCAL_WHOLE_HAPLOTYPE_PAIR",
            "supported_proportions": "PURE_OR_TWO_WAY_EQUAL",
            "unequal_mixtures": "DEFERRED_TO_D14_RECOMBINATION_MOSAICS",
        },
        "confidentiality": {
            "real_sample_files_read": False,
            "raw_donor_identifiers_in_d12_artifacts": False,
            "ancestry_labels_in_blinded_panel": False,
            "local_source_paths_in_manifest": False,
            "artifacts_version_control_allowed": False,
        },
        "artifacts": [
            _artifact("DENSE_PHASED_POPULATION_PANEL", "panel.vcf.gz", 2),
            _artifact("ANCESTRY_TRUTH", "truth.tsv.gz", 4),
            _artifact("D11_HAPLOTYPE_PROVENANCE", "provenance.tsv.gz", 8),
            _artifact("BLINDED_OBSERVED_POPULATION_PANEL", "observed.vcf.gz", 1),
        ],
        "integrity": {
            "dense_variant_count": 2,
            "observed_variant_count": 1,
            "individual_count": 4,
            "source_haplotype_count": 8,
            "output_haplotype_count": 8,
            "haplotype_copy_error_count": 0,
            "allele_count_difference": 0,
            "genotype_haplotype_discordance_count": 0,
            "invalid_proportion_count": 0,
            "proportion_sum_error_count": 0,
            "panel_role_collision_count": 0,
            "checks": [
                "D11A_D11B_BINDING_MATCHED",
                "PANEL_AND_DOMAIN_MATCHED",
                "HAPLOTYPE_COPY_BIJECTION",
                "ANCESTRY_PROPORTIONS_VALID",
                "ANCESTRY_PROPORTIONS_SUM_TO_ONE",
                "CALIBRATION_VALIDATION_ROLE_SEPARATED",
                "OBSERVED_MARKER_UNIVERSE_PRESERVED",
                "BLINDED_PANEL_HAS_NO_ANCESTRY_LABEL",
            ],
        },
    }


@pytest.mark.parametrize(
    ("schema_path", "instance"),
    [(CONFIG_SCHEMA, _config()), (MANIFEST_SCHEMA, _manifest())],
)
def test_d12_contracts_are_valid(
    schema_path: Path, instance: dict[str, object]
) -> None:
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.Draft202012Validator(
        schema, format_checker=jsonschema.FormatChecker()
    ).validate(instance)


def test_config_requires_even_admixed_count() -> None:
    config = deepcopy(_config())
    config["mixtures"][0]["individual_count"] = 3
    schema = json.loads(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(config)


def test_manifest_forbids_ancestry_label_in_blinded_panel() -> None:
    manifest = deepcopy(_manifest())
    manifest["confidentiality"]["ancestry_labels_in_blinded_panel"] = True
    schema = json.loads(MANIFEST_SCHEMA.read_text(encoding="utf-8"))
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(manifest)
