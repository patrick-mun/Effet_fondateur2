"""Tests d'intégration locale de la matérialisation D12."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_population_structure import (
    materialize_population_structure,
)
from effet_fondateur.calibration.simulation.population_structure import (
    PopulationStructureError,
)


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "calibration/schemas/population_structure_manifest.schema.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_gzip(path: Path, content: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(content)


def _fixture(tmp_path: Path) -> dict[str, object]:
    workspace = tmp_path / "calibration_workspace"
    d11a_root = workspace / "d11a_test"
    d11a_private = d11a_root / "private_truth/d11"
    d11a_private.mkdir(parents=True)
    samples = [f"SYN_CAL_{index:06d}" for index in range(1, 5)]
    phased_path = d11a_private / "base_haplotypes.phased.vcf.gz"
    _write_gzip(
        phased_path,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|0\t1|1\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\n",
    )
    provenance_path = d11a_private / "donor_provenance.tsv.gz"
    provenance_rows = [
        "synthetic_identifier\tsynthetic_haplotype_index\tsource_identifier\t"
        "source_haplotype_index\tsource_population"
    ]
    for index, sample in enumerate(samples, start=1):
        population = "POP_A" if index <= 2 else "POP_B"
        provenance_rows.extend(
            (
                f"{sample}\t1\tDONOR_{index}A\t1\t{population}",
                f"{sample}\t2\tDONOR_{index}B\t2\t{population}",
            )
        )
    _write_gzip(provenance_path, "\n".join(provenance_rows) + "\n")
    d11a_manifest = {
        "base_haplotype_manifest_id": "d11a_test.d11",
        "manifest_state": "FINALIZED",
        "panel": {
            "role": "CALIBRATION",
            "individual_count": 4,
            "haplotype_count": 8,
        },
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 200},
        "artifacts": [
            {
                "artifact_id": "d11a_test.phased",
                "role": "PHASED_HAPLOTYPES",
                "path": "private_truth/d11/base_haplotypes.phased.vcf.gz",
                "sha256": _sha256(phased_path),
                "record_count": 2,
            },
            {
                "artifact_id": "d11a_test.provenance",
                "role": "DONOR_PROVENANCE",
                "path": "private_truth/d11/donor_provenance.tsv.gz",
                "sha256": _sha256(provenance_path),
                "record_count": 8,
            },
        ],
        "integrity": {"variant_count": 2},
    }
    d11a_manifest_path = d11a_private / "base_haplotype_manifest.json"
    d11a_manifest_path.write_text(
        json.dumps(d11a_manifest, indent=2) + "\n", encoding="utf-8"
    )

    d11b_root = workspace / "d11b_test"
    d11b_private = d11b_root / "private_truth/d11b"
    d11b_private.mkdir(parents=True)
    observed_path = d11b_private / "observed_chip_genotypes.vcf.gz"
    _write_gzip(
        observed_path,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        "chr19\t200\tCHIP-2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\n",
    )
    d11b_manifest = {
        "chip_harmonization_manifest_id": "d11b_test.d11b",
        "manifest_state": "FINALIZED",
        "dense_truth": {
            "base_haplotype_manifest_id": "d11a_test.d11",
            "base_haplotype_manifest_sha256": _sha256(d11a_manifest_path),
            "panel_role": "CALIBRATION",
        },
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 200},
        "artifacts": [
            {
                "artifact_id": "d11b_test.observed",
                "role": "CHIP_OBSERVED_GENOTYPES",
                "path": "private_truth/d11b/observed_chip_genotypes.vcf.gz",
                "sha256": _sha256(observed_path),
                "record_count": 1,
            }
        ],
        "integrity": {"observed_vcf_variant_count": 1},
    }
    (d11b_private / "chip_harmonization_manifest.json").write_text(
        json.dumps(d11b_manifest, indent=2) + "\n", encoding="utf-8"
    )

    config_path = tmp_path / "d12_config.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "population_structure_config_id": "d12_test.v1",
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
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "d11a_root": d11a_root,
        "d11b_root": d11b_root,
        "population_config_path": config_path,
        "workspace_root": workspace,
        "run_id": "d12_test",
        "created_at": "2026-10-29T08:00:00Z",
    }


def test_materialization_links_d11_and_separates_truth_from_blinded_panel(
    tmp_path: Path,
) -> None:
    arguments = _fixture(tmp_path)
    manifest = materialize_population_structure(**arguments)
    jsonschema.Draft202012Validator(
        json.loads(SCHEMA.read_text(encoding="utf-8")),
        format_checker=jsonschema.FormatChecker(),
    ).validate(manifest)

    private_root = arguments["workspace_root"] / "d12_test/private_truth/d12"
    assert manifest["panel"]["admixed_individual_count"] == 2
    assert manifest["panel"]["pure_individual_count"] == 2
    assert manifest["integrity"]["allele_count_difference"] == 0
    with gzip.open(
        private_root / "ancestry_truth.tsv.gz", "rt", encoding="utf-8"
    ) as handle:
        truth_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(truth_rows) == 4
    assert sum(row["stratum_id"] == "ADMIX_A_B" for row in truth_rows) == 2
    with gzip.open(
        private_root / "blinded_observed_panel.vcf.gz", "rt", encoding="utf-8"
    ) as handle:
        blinded = handle.read()
    assert "POP_A" not in blinded
    assert "ADMIX_A_B" not in blinded
    assert "DONOR_" not in blinded
    assert "\t100\t" not in blinded
    assert "\t200\t" in blinded
    assert "\tCHIP-2\t" in blinded
    assert all(path.stat().st_mode & 0o077 == 0 for path in private_root.iterdir())
    assert str(tmp_path) not in json.dumps(manifest)


def test_materialization_rejects_d11b_bound_to_another_d11a(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    manifest_path = (
        arguments["d11b_root"]
        / "private_truth/d11b/chip_harmonization_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["dense_truth"]["base_haplotype_manifest_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")
    with pytest.raises(
        PopulationStructureError, match="population_structure_d11_binding_mismatch"
    ):
        materialize_population_structure(**arguments)


def test_materialization_rejects_traversing_run_identifier(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["run_id"] = "../outside"
    with pytest.raises(
        PopulationStructureError,
        match="population_structure_run_identifier_invalid",
    ):
        materialize_population_structure(**arguments)
    assert not (tmp_path / "outside").exists()
