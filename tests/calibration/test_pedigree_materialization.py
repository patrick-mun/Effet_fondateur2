"""Tests d'intégration locale de la matérialisation D13."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_pedigrees import (
    materialize_pedigrees,
)
from effet_fondateur.calibration.simulation.pedigrees import (
    PedigreeSimulationError,
)


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "calibration/schemas/pedigree_manifest.schema.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gzip(path: Path, payload: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(payload)


def _fixture(tmp_path: Path) -> dict[str, object]:
    workspace = tmp_path / "calibration_workspace"
    d12_root = workspace / "d12_test"
    private = d12_root / "private_truth/d12"
    private.mkdir(parents=True)
    samples = [f"SYN_CAL_D12_{index:06d}" for index in range(1, 9)]

    dense_path = private / "population_panel.phased.vcf.gz"
    _gzip(
        dense_path,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        "chr19\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t"
        "0|0\t0|1\t1|0\t1|1\t0|1\t1|0\t0|0\t1|1\n"
        "chr19\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t"
        "0|1\t1|0\t0|0\t1|1\t0|1\t1|0\t0|0\t1|1\n",
    )
    observed_path = private / "blinded_observed_panel.vcf.gz"
    _gzip(
        observed_path,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        "chr19\t200\tCHIP-2\tC\tT\t.\tPASS\t.\tGT\t"
        "0/1\t0/1\t0/0\t1/1\t0/1\t0/1\t0/0\t1/1\n",
    )
    ancestry_path = private / "ancestry_truth.tsv.gz"
    _gzip(
        ancestry_path,
        "sample_id\tstratum_id\tmixture_proportions\n"
        + "".join(f'{sample}\tPOP_A\t"{{""POP_A"":1.0}}"\n' for sample in samples),
    )
    provenance_path = private / "d11_haplotype_provenance.tsv.gz"
    _gzip(
        provenance_path,
        "sample_id\toutput_haplotype_index\td11_sample_id\t"
        "d11_haplotype_index\tsource_stratum_id\n"
        + "".join(
            f"{sample}\t{haplotype_index}\t"
            f"SYN_CAL_D11_{sample_index:06d}\t{haplotype_index}\tPOP_A\n"
            for sample_index, sample in enumerate(samples, start=1)
            for haplotype_index in (1, 2)
        ),
    )
    artifacts = [
        {
            "artifact_id": "d12.dense",
            "role": "DENSE_PHASED_POPULATION_PANEL",
            "path": "private_truth/d12/population_panel.phased.vcf.gz",
            "sha256": _sha256(dense_path),
            "record_count": 2,
        },
        {
            "artifact_id": "d12.observed",
            "role": "BLINDED_OBSERVED_POPULATION_PANEL",
            "path": "private_truth/d12/blinded_observed_panel.vcf.gz",
            "sha256": _sha256(observed_path),
            "record_count": 1,
        },
        {
            "artifact_id": "d12.ancestry",
            "role": "ANCESTRY_TRUTH",
            "path": "private_truth/d12/ancestry_truth.tsv.gz",
            "sha256": _sha256(ancestry_path),
            "record_count": 8,
        },
        {
            "artifact_id": "d12.provenance",
            "role": "D11_HAPLOTYPE_PROVENANCE",
            "path": "private_truth/d12/d11_haplotype_provenance.tsv.gz",
            "sha256": _sha256(provenance_path),
            "record_count": 16,
        },
    ]
    manifest = {
        "population_structure_manifest_id": "d12_test.d12",
        "manifest_state": "FINALIZED",
        "panel": {"role": "CALIBRATION", "individual_count": 8},
        "domain": {"chromosome": 19, "start_bp": 100, "end_bp": 200},
        "artifacts": artifacts,
        "integrity": {"dense_variant_count": 2, "observed_variant_count": 1},
    }
    (private / "population_structure_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    config_path = tmp_path / "pedigree_config.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "pedigree_config_id": "d13_test",
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
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "d12_root": d12_root,
        "pedigree_config_path": config_path,
        "workspace_root": workspace,
        "run_id": "d13_test",
        "created_at": "2026-10-29T10:00:00Z",
    }


def test_materialization_produces_exact_mendelian_pedigrees(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    manifest = materialize_pedigrees(**arguments)
    jsonschema.Draft202012Validator(
        json.loads(SCHEMA.read_text(encoding="utf-8")),
        format_checker=jsonschema.FormatChecker(),
    ).validate(manifest)
    private = arguments["workspace_root"] / "d13_test/private_truth/d13"
    assert manifest["panel"]["output_individual_count"] == 12
    assert manifest["panel"]["trio_count"] == 4
    assert manifest["integrity"]["mendelian_error_count"] == 0
    assert manifest["integrity"]["relationship_class_counts"] == {
        "UNRELATED": 13,
        "FIRST_DEGREE": 9,
        "SECOND_DEGREE": 6,
        "THIRD_DEGREE": 1,
    }
    with gzip.open(private / "pedigree_truth.tsv.gz", "rt", encoding="utf-8") as handle:
        pedigree = list(csv.DictReader(handle, delimiter="\t"))
    assert len(pedigree) == 12
    assert sum(row["is_founder"] == "false" for row in pedigree) == 4
    with gzip.open(
        private / "blinded_observed_panel.vcf.gz", "rt", encoding="utf-8"
    ) as handle:
        blinded = handle.read()
    assert "SYN_CAL_D12_" not in blinded
    assert "PED_A" not in blinded
    assert "POP_A" not in blinded
    assert "\tCHIP-2\t" in blinded
    assert all(path.stat().st_mode & 0o077 == 0 for path in private.iterdir())
    assert str(tmp_path) not in json.dumps(manifest)


def test_materialization_rejects_corrupted_d12_artifact(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    dense = (
        arguments["d12_root"]
        / "private_truth/d12/population_panel.phased.vcf.gz"
    )
    dense.write_bytes(dense.read_bytes() + b"corrupt")
    with pytest.raises(
        PedigreeSimulationError, match="pedigree_d12_artifact_sha256_mismatch"
    ):
        materialize_pedigrees(**arguments)


def test_materialization_rejects_traversing_run_id(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    arguments["run_id"] = "../outside"
    with pytest.raises(
        PedigreeSimulationError, match="pedigree_run_identifier_invalid"
    ):
        materialize_pedigrees(**arguments)
    assert not (tmp_path / "outside").exists()
