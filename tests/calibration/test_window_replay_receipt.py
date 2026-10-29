"""Tests du reçu D14b liant les manifests D11a à D14."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from effet_fondateur.calibration.materialize_window_replay_receipt import (
    materialize_window_replay_receipt,
)
from effet_fondateur.calibration.simulation.windows import SimulationWindowError


ROOT = Path(__file__).resolve().parents[2]
DOMAIN = {"chromosome": 5, "start_bp": 300, "end_bp": 900}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(root: Path, relative: str, document: dict[str, object]) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def _fixture(tmp_path: Path) -> dict[str, object]:
    workspace = tmp_path / "calibration_workspace"
    roots = {stage: workspace / stage.casefold() for stage in (
        "D14B", "D11A", "D11B", "D12", "D13", "D14"
    )}
    source_sha = "1" * 64
    index_sha = "2" * 64
    map_sha = "3" * 64
    _write(
        roots["D14B"],
        "private_truth/d14b/simulation_window_manifest.json",
        {
            "simulation_window_manifest_id": "d14b.design",
            "manifest_state": "FINALIZED",
            "configuration": {"configuration_id": "d14b.config"},
            "panel": {"role": "CALIBRATION"},
            "inputs": {
                "source_vcf_sha256": source_sha,
                "source_index_sha256": index_sha,
                "genetic_map_sha256": map_sha,
            },
            "windows": [
                {
                    "window_id": "primary",
                    "purpose": "PRIMARY",
                    "chromosome": 5,
                    "domain": {"start_bp": 290, "end_bp": 910},
                    "target": {"position_bp": 600},
                },
                {"window_id": "challenge", "purpose": "CENSORING_CHALLENGE"},
            ],
        },
    )
    d11_path = _write(
        roots["D11A"],
        "private_truth/d11/base_haplotype_manifest.json",
        {
            "base_haplotype_manifest_id": "d11.manifest",
            "manifest_state": "FINALIZED",
            "panel": {"role": "CALIBRATION"},
            "source": {
                "vcf_sha256": source_sha,
                "index_sha256": index_sha,
                "genetic_map_sha256": map_sha,
            },
            "domain": DOMAIN,
        },
    )
    d11_sha = _sha(d11_path)
    d11b_path = _write(
        roots["D11B"],
        "private_truth/d11b/chip_harmonization_manifest.json",
        {
            "chip_harmonization_manifest_id": "d11b.manifest",
            "manifest_state": "FINALIZED",
            "dense_truth": {
                "base_haplotype_manifest_id": "d11.manifest",
                "base_haplotype_manifest_sha256": d11_sha,
                "panel_role": "CALIBRATION",
            },
            "domain": DOMAIN,
        },
    )
    d11b_sha = _sha(d11b_path)
    d12_path = _write(
        roots["D12"],
        "private_truth/d12/population_structure_manifest.json",
        {
            "population_structure_manifest_id": "d12.manifest",
            "manifest_state": "FINALIZED",
            "inputs": {
                "d11a": {
                    "manifest_id": "d11.manifest",
                    "manifest_sha256": d11_sha,
                },
                "d11b": {
                    "manifest_id": "d11b.manifest",
                    "manifest_sha256": d11b_sha,
                },
            },
            "panel": {"role": "CALIBRATION"},
            "domain": DOMAIN,
        },
    )
    d12_sha = _sha(d12_path)
    d13_path = _write(
        roots["D13"],
        "private_truth/d13/pedigree_manifest.json",
        {
            "pedigree_manifest_id": "d13.manifest",
            "manifest_state": "FINALIZED",
            "input": {
                "manifest_id": "d12.manifest",
                "manifest_sha256": d12_sha,
            },
            "panel": {"role": "CALIBRATION"},
            "domain": DOMAIN,
        },
    )
    d13_sha = _sha(d13_path)
    _write(
        roots["D14"],
        "private_truth/d14/founder_simulation_manifest.json",
        {
            "founder_simulation_manifest_id": "d14.manifest",
            "manifest_state": "FINALIZED",
            "inputs": {
                "d11a_manifest_id": "d11.manifest",
                "d11a_manifest_sha256": d11_sha,
                "d12_manifest_id": "d12.manifest",
                "d12_manifest_sha256": d12_sha,
                "d13_manifest_id": "d13.manifest",
                "d13_manifest_sha256": d13_sha,
            },
            "panel": {"role": "CALIBRATION"},
            "domain": {**DOMAIN, "target_position_bp": 600},
            "integrity": {
                "dense_variant_count": 100,
                "observed_variant_count": 20,
                "family_count": 2,
                "recombination_event_count": 4,
                "left_censored_count": 0,
                "right_censored_count": 0,
                "mendelian_error_count": 0,
            },
        },
    )
    return {
        "d14b_root": roots["D14B"],
        "d11a_root": roots["D11A"],
        "d11b_root": roots["D11B"],
        "d12_root": roots["D12"],
        "d13_root": roots["D13"],
        "d14_root": roots["D14"],
        "workspace_root": workspace,
        "run_id": "d14b_replay_test",
        "created_at": "2026-10-29T10:00:00Z",
    }


def test_replay_receipt_validates_complete_hash_chain(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    receipt = materialize_window_replay_receipt(**arguments)
    schema = json.loads(
        (ROOT / "calibration/schemas/simulation_window_replay.schema.json").read_text(
            encoding="utf-8"
        )
    )
    jsonschema.Draft202012Validator(
        schema, format_checker=jsonschema.FormatChecker()
    ).validate(receipt)
    assert receipt["result"]["recombination_event_count"] == 4
    assert receipt["result"]["mendelian_error_count"] == 0


def test_replay_receipt_refuses_tampered_d12_binding(tmp_path: Path) -> None:
    arguments = _fixture(tmp_path)
    manifest_path = (
        arguments["d12_root"]
        / "private_truth/d12/population_structure_manifest.json"
    )
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    document["inputs"]["d11a"]["manifest_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(
        SimulationWindowError, match="simulation_window_replay_d12_binding_mismatch"
    ):
        materialize_window_replay_receipt(**arguments)
