"""Scellement D14b du rejeu D11a–D14 sur la fenêtre principale."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.calibration.simulation.windows import SimulationWindowError


_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_CHECKS = (
    "D14B_DESIGN_HASH_RECORDED",
    "D11_PUBLIC_SOURCE_BINDING_MATCHED",
    "D11B_D11_BINDING_MATCHED",
    "D12_D11A_D11B_BINDINGS_MATCHED",
    "D13_D12_BINDING_MATCHED",
    "D14_D11_D12_D13_BINDINGS_MATCHED",
    "PRIMARY_DOMAIN_CONTAINED",
    "TARGET_POSITION_MATCHED",
    "INTERNAL_RECOMBINATION_OBSERVED",
    "MENDELIAN_ERRORS_ZERO",
)
_MANIFESTS = {
    "D14B": "private_truth/d14b/simulation_window_manifest.json",
    "D11A": "private_truth/d11/base_haplotype_manifest.json",
    "D11B": "private_truth/d11b/chip_harmonization_manifest.json",
    "D12": "private_truth/d12/population_structure_manifest.json",
    "D13": "private_truth/d13/pedigree_manifest.json",
    "D14": "private_truth/d14/founder_simulation_manifest.json",
}
_ID_FIELDS = {
    "D11A": "base_haplotype_manifest_id",
    "D11B": "chip_harmonization_manifest_id",
    "D12": "population_structure_manifest_id",
    "D13": "pedigree_manifest_id",
    "D14": "founder_simulation_manifest_id",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_manifest(root: Path, stage: str) -> tuple[dict[str, Any], Path, str]:
    if (
        not root.is_dir()
        or root.is_symlink()
        or root.parent.name != "calibration_workspace"
    ):
        raise SimulationWindowError("simulation_window_replay_root_invalid", stage)
    path = root / _MANIFESTS[stage]
    if not path.is_file() or path.is_symlink():
        raise SimulationWindowError("simulation_window_replay_manifest_missing", stage)
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SimulationWindowError(
            "simulation_window_replay_manifest_invalid", stage
        ) from error
    if not isinstance(document, dict) or document.get("manifest_state") != "FINALIZED":
        raise SimulationWindowError("simulation_window_replay_manifest_invalid", stage)
    return document, path, _sha256(path)


def _domain(document: dict[str, Any], stage: str) -> dict[str, Any]:
    domain = document.get("domain")
    if not isinstance(domain, dict):
        raise SimulationWindowError("simulation_window_replay_domain_invalid", stage)
    required = ("chromosome", "start_bp", "end_bp")
    if any(type(domain.get(key)) is not int for key in required):
        raise SimulationWindowError("simulation_window_replay_domain_invalid", stage)
    return {key: domain[key] for key in required}


def _require_equal(actual: object, expected: object, code: str) -> None:
    if actual != expected:
        raise SimulationWindowError(code)


def materialize_window_replay_receipt(
    *,
    d14b_root: Path,
    d11a_root: Path,
    d11b_root: Path,
    d12_root: Path,
    d13_root: Path,
    d14_root: Path,
    workspace_root: Path,
    run_id: str,
    created_at: str,
) -> dict[str, Any]:
    """Vérifie les liens D14b→D11a→D11b→D12→D13→D14 et écrit un reçu."""
    if _IDENTIFIER.fullmatch(run_id) is None:
        raise SimulationWindowError("simulation_window_replay_id_invalid")
    if workspace_root.name != "calibration_workspace" or workspace_root.is_symlink():
        raise SimulationWindowError("simulation_window_workspace_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise SimulationWindowError("simulation_window_run_exists")
    design, _, design_sha = _load_manifest(d14b_root, "D14B")
    loaded = {
        stage: _load_manifest(root, stage)
        for stage, root in (
            ("D11A", d11a_root),
            ("D11B", d11b_root),
            ("D12", d12_root),
            ("D13", d13_root),
            ("D14", d14_root),
        )
    }
    primary_windows = [
        window
        for window in design.get("windows", [])
        if isinstance(window, dict) and window.get("purpose") == "PRIMARY"
    ]
    if len(primary_windows) != 1:
        raise SimulationWindowError("simulation_window_replay_primary_invalid")
    primary = primary_windows[0]
    design_domain = primary.get("domain")
    target = primary.get("target")
    if not isinstance(design_domain, dict) or not isinstance(target, dict):
        raise SimulationWindowError("simulation_window_replay_primary_invalid")
    d11, d11_path, d11_sha = loaded["D11A"]
    d11b, d11b_path, d11b_sha = loaded["D11B"]
    d12, d12_path, d12_sha = loaded["D12"]
    d13, d13_path, d13_sha = loaded["D13"]
    d14, d14_path, d14_sha = loaded["D14"]
    design_inputs = design.get("inputs", {})
    d11_source = d11.get("source", {})
    for d11_key, design_key in (
        ("vcf_sha256", "source_vcf_sha256"),
        ("index_sha256", "source_index_sha256"),
        ("genetic_map_sha256", "genetic_map_sha256"),
    ):
        _require_equal(
            d11_source.get(d11_key),
            design_inputs.get(design_key),
            "simulation_window_replay_public_binding_mismatch",
        )
    d11b_dense = d11b.get("dense_truth", {})
    _require_equal(
        d11b_dense.get("base_haplotype_manifest_id"),
        d11.get("base_haplotype_manifest_id"),
        "simulation_window_replay_d11b_binding_mismatch",
    )
    _require_equal(
        d11b_dense.get("base_haplotype_manifest_sha256"),
        d11_sha,
        "simulation_window_replay_d11b_binding_mismatch",
    )
    d12_inputs = d12.get("inputs", {})
    for key, expected_id, expected_sha in (
        ("d11a", d11.get("base_haplotype_manifest_id"), d11_sha),
        ("d11b", d11b.get("chip_harmonization_manifest_id"), d11b_sha),
    ):
        binding = d12_inputs.get(key, {})
        _require_equal(
            (binding.get("manifest_id"), binding.get("manifest_sha256")),
            (expected_id, expected_sha),
            "simulation_window_replay_d12_binding_mismatch",
        )
    d13_input = d13.get("input", {})
    _require_equal(
        (d13_input.get("manifest_id"), d13_input.get("manifest_sha256")),
        (d12.get("population_structure_manifest_id"), d12_sha),
        "simulation_window_replay_d13_binding_mismatch",
    )
    d14_inputs = d14.get("inputs", {})
    for id_key, sha_key, expected_id, expected_sha in (
        (
            "d11a_manifest_id",
            "d11a_manifest_sha256",
            d11.get("base_haplotype_manifest_id"),
            d11_sha,
        ),
        (
            "d12_manifest_id",
            "d12_manifest_sha256",
            d12.get("population_structure_manifest_id"),
            d12_sha,
        ),
        (
            "d13_manifest_id",
            "d13_manifest_sha256",
            d13.get("pedigree_manifest_id"),
            d13_sha,
        ),
    ):
        _require_equal(
            (d14_inputs.get(id_key), d14_inputs.get(sha_key)),
            (expected_id, expected_sha),
            "simulation_window_replay_d14_binding_mismatch",
        )
    stage_documents = {
        "D11A": d11,
        "D11B": d11b,
        "D12": d12,
        "D13": d13,
        "D14": d14,
    }
    stage_paths = {
        "D11A": (d11_path, d11_sha),
        "D11B": (d11b_path, d11b_sha),
        "D12": (d12_path, d12_sha),
        "D13": (d13_path, d13_sha),
        "D14": (d14_path, d14_sha),
    }
    domains = {stage: _domain(document, stage) for stage, document in stage_documents.items()}
    canonical_domain = domains["D11A"]
    if any(domain != canonical_domain for domain in domains.values()):
        raise SimulationWindowError("simulation_window_replay_stage_domain_mismatch")
    if (
        canonical_domain["chromosome"] != primary.get("chromosome")
        or canonical_domain["start_bp"] < design_domain.get("start_bp", 0)
        or canonical_domain["end_bp"] > design_domain.get("end_bp", 0)
    ):
        raise SimulationWindowError("simulation_window_replay_domain_not_contained")
    d14_domain = d14.get("domain", {})
    _require_equal(
        d14_domain.get("target_position_bp"),
        target.get("position_bp"),
        "simulation_window_replay_target_mismatch",
    )
    d14_integrity = d14.get("integrity", {})
    if d14_integrity.get("recombination_event_count", 0) < 1:
        raise SimulationWindowError("simulation_window_replay_internal_event_missing")
    if d14_integrity.get("mendelian_error_count") != 0:
        raise SimulationWindowError("simulation_window_replay_mendelian_error")
    panel_roles = {
        d11.get("panel", {}).get("role"),
        d11b.get("dense_truth", {}).get("panel_role"),
        d12.get("panel", {}).get("role"),
        d13.get("panel", {}).get("role"),
        d14.get("panel", {}).get("role"),
    }
    if panel_roles != {design.get("panel", {}).get("role")}:
        raise SimulationWindowError("simulation_window_replay_panel_role_mismatch")
    receipt = {
        "schema_version": "1.0.0",
        "simulation_window_replay_id": f"{run_id}.replay",
        "receipt_state": "FINALIZED",
        "created_at": created_at,
        "design": {
            "manifest_id": design["simulation_window_manifest_id"],
            "manifest_sha256": design_sha,
            "configuration_id": design["configuration"]["configuration_id"],
            "primary_window_id": primary["window_id"],
            "chromosome": primary["chromosome"],
            "start_bp": design_domain["start_bp"],
            "end_bp": design_domain["end_bp"],
            "target_position_bp": target["position_bp"],
        },
        "stages": [
            {
                "stage": stage,
                "manifest_id": stage_documents[stage][_ID_FIELDS[stage]],
                "manifest_sha256": stage_paths[stage][1],
                "chromosome": domains[stage]["chromosome"],
                "start_bp": domains[stage]["start_bp"],
                "end_bp": domains[stage]["end_bp"],
            }
            for stage in ("D11A", "D11B", "D12", "D13", "D14")
        ],
        "result": {
            "panel_role": d14["panel"]["role"],
            "dense_variant_count": d14_integrity["dense_variant_count"],
            "observed_variant_count": d14_integrity["observed_variant_count"],
            "family_count": d14_integrity["family_count"],
            "recombination_event_count": d14_integrity["recombination_event_count"],
            "left_censored_count": d14_integrity["left_censored_count"],
            "right_censored_count": d14_integrity["right_censored_count"],
            "mendelian_error_count": d14_integrity["mendelian_error_count"],
        },
        "confidentiality": {
            "real_sample_files_read": False,
            "sample_identifiers_recorded": False,
            "local_source_paths_recorded": False,
            "version_control_allowed": False,
        },
        "checks": list(_CHECKS),
    }
    temporary_root = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private_root = temporary_root / "private_truth/d14b_replay"
        private_root.mkdir(parents=True, mode=0o700)
        receipt_path = private_root / "simulation_window_replay.json"
        receipt_path.write_text(
            json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        receipt_path.chmod(0o600)
        temporary_root.rename(target_root)
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    return receipt


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Scelle le rejeu D14b de D11a à D14.")
    parser.add_argument("--d14b-root", type=Path, required=True)
    parser.add_argument("--d11a-root", type=Path, required=True)
    parser.add_argument("--d11b-root", type=Path, required=True)
    parser.add_argument("--d12-root", type=Path, required=True)
    parser.add_argument("--d13-root", type=Path, required=True)
    parser.add_argument("--d14-root", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    return parser


def main() -> int:
    """Point d'entrée du reçu D14b ; n'affiche que les résultats agrégés."""
    arguments = _parser().parse_args()
    receipt = materialize_window_replay_receipt(
        d14b_root=arguments.d14b_root,
        d11a_root=arguments.d11a_root,
        d11b_root=arguments.d11b_root,
        d12_root=arguments.d12_root,
        d13_root=arguments.d13_root,
        d14_root=arguments.d14_root,
        workspace_root=arguments.workspace_root,
        run_id=arguments.run_id,
        created_at=datetime.now(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z"),
    )
    print(json.dumps(receipt["result"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
