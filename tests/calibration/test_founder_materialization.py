"""Test d'intégration locale de la matérialisation D14."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import jsonschema

from effet_fondateur.calibration.materialize_founder_recombination import (
    materialize_founder_recombination,
)


ROOT = Path(__file__).resolve().parents[2]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _artifact(identifier: str, role: str, relative: str, path: Path, count: int):
    return {
        "artifact_id": identifier,
        "role": role,
        "path": relative,
        "sha256": _sha(path),
        "record_count": count,
    }


def _fixture(tmp_path: Path) -> dict[str, object]:
    workspace = tmp_path / "calibration_workspace"
    d11 = workspace / "d11"
    d12 = workspace / "d12"
    d13 = workspace / "d13"
    d11_private = d11 / "private_truth/d11"
    d12_private = d12 / "private_truth/d12"
    d13_private = d13 / "private_truth/d13"
    for path in (d11_private, d12_private, d13_private):
        path.mkdir(parents=True)

    map_path = d11_private / "genetic_map.tsv"
    map_path.write_text(
        "CHROMOSOME\tPOSITION_BP\tPOSITION_CM\n"
        "19\t100\t0\n19\t200\t50\n19\t300\t100\n",
        encoding="utf-8",
    )
    d11_manifest = {
        "base_haplotype_manifest_id": "d11.manifest",
        "manifest_state": "FINALIZED",
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 300},
        "artifacts": [
            _artifact(
                "d11.map",
                "GENETIC_MAP",
                "private_truth/d11/genetic_map.tsv",
                map_path,
                3,
            )
        ],
    }
    d11_manifest_path = d11_private / "base_haplotype_manifest.json"
    d11_manifest_path.write_text(json.dumps(d11_manifest), encoding="utf-8")

    d12_provenance = d12_private / "d11_haplotype_provenance.tsv.gz"
    _gzip(
        d12_provenance,
        "sample_id\toutput_haplotype_index\td11_sample_id\t"
        "d11_haplotype_index\tsource_stratum_id\n"
        + "".join(
            f"D{sample}\t{hap}\tX{sample}\t{hap}\tPOP_A\n"
            for sample in range(1, 6)
            for hap in (1, 2)
        ),
    )
    d12_manifest = {
        "population_structure_manifest_id": "d12.manifest",
        "manifest_state": "FINALIZED",
        "inputs": {
            "d11a": {
                "manifest_id": "d11.manifest",
                "manifest_sha256": _sha(d11_manifest_path),
            }
        },
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 300},
        "artifacts": [
            _artifact(
                "d12.provenance",
                "D11_HAPLOTYPE_PROVENANCE",
                "private_truth/d12/d11_haplotype_provenance.tsv.gz",
                d12_provenance,
                10,
            )
        ],
    }
    d12_manifest_path = d12_private / "population_structure_manifest.json"
    d12_manifest_path.write_text(json.dumps(d12_manifest), encoding="utf-8")

    samples = ("S1", "S2", "S3", "S4", "S5")
    dense = d13_private / "pedigree_panel.phased.vcf.gz"
    _gzip(
        dense,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\t1|0\t1|1\t1|0\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\n"
        "chr19\t300\tv3\tG\tA\t.\tPASS\t.\tGT\t0|1\t0|0\t1|0\t1|1\t1|0\n",
    )
    observed = d13_private / "blinded_observed_panel.vcf.gz"
    _gzip(
        observed,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        "chr19\t100\tCHIP-1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\t0/1\n",
    )
    pedigree = d13_private / "pedigree_truth.tsv.gz"
    _gzip(
        pedigree,
        "sample_id\tfamily_unit_id\tfather_id\tmother_id\tsex\t"
        "generation\tis_founder\n"
        "S1\tFAM_A\t0\t0\tMALE\t0\ttrue\n"
        "S2\tFAM_A\tS1\tS5\tFEMALE\t2\tfalse\n"
        "S3\tFAM_B\t0\t0\tMALE\t0\ttrue\n"
        "S4\tFAM_B\tS3\tS5\tFEMALE\t2\tfalse\n"
        "S5\tUNRELATED_S5\t0\t0\tFEMALE\t0\ttrue\n",
    )
    provenance = d13_private / "d12_haplotype_provenance.tsv.gz"
    _gzip(
        provenance,
        "sample_id\td13_haplotype_index\td12_sample_id\td12_haplotype_index\n"
        "S1\t1\tD1\t1\nS1\t2\tD1\t2\n"
        "S2\t1\tD1\t1\nS2\t2\tD2\t2\n"
        "S3\t1\tD3\t1\nS3\t2\tD3\t2\n"
        "S4\t1\tD3\t1\nS4\t2\tD4\t2\n"
        "S5\t1\tD5\t1\nS5\t2\tD5\t2\n",
    )
    transmissions = d13_private / "mendelian_transmission_truth.tsv.gz"
    _gzip(
        transmissions,
        "family_unit_id\tchild_id\tparent_id\tparental_haplotype_index\t"
        "child_haplotype_index\trecombination_count\n"
        "FAM_A\tS2\tS1\t1\t1\t0\n"
        "FAM_B\tS4\tS3\t1\t1\t0\n",
    )
    d13_manifest = {
        "pedigree_manifest_id": "d13.manifest",
        "manifest_state": "FINALIZED",
        "input": {
            "manifest_id": "d12.manifest",
            "manifest_sha256": _sha(d12_manifest_path),
        },
        "panel": {"role": "CALIBRATION", "output_individual_count": 5},
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 300},
        "artifacts": [
            _artifact(
                "d13.dense",
                "DENSE_PHASED_PEDIGREE_PANEL",
                "private_truth/d13/pedigree_panel.phased.vcf.gz",
                dense,
                3,
            ),
            _artifact(
                "d13.observed",
                "BLINDED_OBSERVED_PEDIGREE_PANEL",
                "private_truth/d13/blinded_observed_panel.vcf.gz",
                observed,
                1,
            ),
            _artifact(
                "d13.pedigree",
                "PEDIGREE_TRUTH",
                "private_truth/d13/pedigree_truth.tsv.gz",
                pedigree,
                5,
            ),
            _artifact(
                "d13.transmissions",
                "MENDELIAN_TRANSMISSION_TRUTH",
                "private_truth/d13/mendelian_transmission_truth.tsv.gz",
                transmissions,
                2,
            ),
            _artifact(
                "d13.provenance",
                "D12_HAPLOTYPE_PROVENANCE",
                "private_truth/d13/d12_haplotype_provenance.tsv.gz",
                provenance,
                10,
            ),
        ],
    }
    (d13_private / "pedigree_manifest.json").write_text(
        json.dumps(d13_manifest), encoding="utf-8"
    )
    config = tmp_path / "d14.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "founder_simulation_config_id": "d14_test",
                "config_state": "FINALIZED",
                "replication_id": "D14_TEST_R001",
                "method": {
                    "method_id": "map_first_breakpoint_founder_mosaic_v1",
                    "method_version": "1.0.0",
                    "random_seed": 1004,
                    "age_generations": 150,
                    "family_selection": "ALL_D13_PEDIGREE_FAMILIES",
                    "carrier_lineage_selection": (
                        "RANDOM_TERMINAL_DESCENDANT_HAPLOTYPE"
                    ),
                    "founder_haplotype_selection": (
                        "RANDOM_UNRELATED_D13_HAPLOTYPE"
                    ),
                    "breakpoint_process": (
                        "INDEPENDENT_EXPONENTIAL_FIRST_EVENT_EACH_SIDE"
                    ),
                    "edge_policy": "RIGHT_CENSOR_AT_MAP_DOMAIN",
                },
                "target_variant": {
                    "variant_id": "SYNTHETIC_TARGET_D14",
                    "chromosome": 19,
                    "position_bp": 200,
                    "reference_allele": "C",
                    "alternate_allele": "T",
                    "position_policy": (
                        "EXISTING_BIALLELIC_SNV_RESET_TO_SYNTHETIC_MUTATION"
                    ),
                },
            }
        ),
        encoding="utf-8",
    )
    return {
        "d13_root": d13,
        "d12_root": d12,
        "d11a_root": d11,
        "founder_config_path": config,
        "workspace_root": workspace,
        "run_id": "d14_test",
        "created_at": "2026-10-29T10:00:00Z",
    }


def test_materialization_preserves_mendel_and_validates_manifest(
    tmp_path: Path,
) -> None:
    arguments = _fixture(tmp_path)
    manifest = materialize_founder_recombination(**arguments)
    schema = json.loads(
        (
            ROOT / "calibration/schemas/founder_simulation_manifest.schema.json"
        ).read_text()
    )
    jsonschema.Draft202012Validator(
        schema, format_checker=jsonschema.FormatChecker()
    ).validate(manifest)
    assert manifest["integrity"]["mendelian_error_count"] == 0
    assert manifest["integrity"]["target_variant_count"] == 1
    assert manifest["panel"]["family_count"] == 2
    private = arguments["workspace_root"] / "d14_test/private_truth/d14"
    with gzip.open(
        private / "blinded_observed_founder_panel.vcf.gz", "rt", encoding="utf-8"
    ) as handle:
        blinded = handle.read()
    assert "SYNTHETIC_TARGET_D14" not in blinded
    assert "FAM_A" not in blinded
    assert all(path.stat().st_mode & 0o077 == 0 for path in private.iterdir())
