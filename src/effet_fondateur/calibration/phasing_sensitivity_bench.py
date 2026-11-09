"""Banc D23b : densité emboîtée et usage du pedigree, sans données réelles."""
from __future__ import annotations

import argparse
import gzip
import json
import logging
from pathlib import Path
import shutil
from typing import Any

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration import phasing_bench
from effet_fondateur.calibration.phasing_preparation import verify_preparation
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.run_qc_bench import checked_stage
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel, read_panel
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document

LOGGER = logging.getLogger(__name__)
METHOD_ID = "paired_density_pedigree_phasing_bench_v1"


def _schemas_dir() -> Path:
    return Path(__file__).resolve().parents[3] / "calibration/schemas"


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge le contrat fermé D23b et vérifie la séparation des graines."""
    config = json.loads(path.read_text())
    validate_json_document(config, "phasing_sensitivity_config.schema.json", schemas_dir=_schemas_dir())
    if config["phasing_seed"] == config["bootstrap_seed"]:
        raise ValueError("d23b_seed_collision")
    return config


def density_profile(panel: MeasurementPanel, profile: dict[str, Any], target: dict[str, Any]) -> MeasurementPanel:
    """Applique le masque de rang préspecifié, séparément sur chaque autosome."""
    chromosomes = np.array([int(row[0].removeprefix("chr")) for row in panel.variants])
    protected = np.array([
        int(row[0].removeprefix("chr")) == target["chromosome"]
        and int(row[1]) == target["position_bp"]
        for row in panel.variants
    ])
    retained = np.zeros(len(panel.variants), dtype=bool)
    residues = set(profile["rank_residues"])
    modulus = profile["rank_modulus"]
    for chromosome in range(1, 23):
        indices = np.flatnonzero((chromosomes == chromosome) & ~protected)
        retained[indices] = [rank % modulus in residues for rank in range(len(indices))]
    retained |= protected
    return MeasurementPanel(
        panel.samples,
        tuple(row for row, keep in zip(panel.variants, retained, strict=True) if keep),
        panel.alleles[retained].copy(),
    )


def _resolved_config(root: Path, base: dict[str, Any]) -> None:
    (root / "config.resolved.yaml").write_text(yaml.safe_dump(base, sort_keys=False))
    load_pipeline_config(root / "config.resolved.yaml")


def _fixture_root(root: Path, scope: str, profile_id: str, pedigree_mode: str | None) -> None:
    root.mkdir(parents=True)
    write_json(root / "manifest.json", {
        "run_id": "d23b_" + "_".join(filter(None, (scope.lower(), profile_id.lower(), (pedigree_mode or "").lower()))),
        "scope": "COMPONENT_BENCH_NOT_FULL_PIPELINE",
    })


def run_bench(*, preparation_root: Path, configuration_path: Path, workspace_root: Path, run_id: str,
              reference_source_dir: Path, shapeit_bin_dir: Path) -> dict[str, Any]:
    """Exécute les parcours 10A/12A et 11/12 préspecifiés pour D23b."""
    destination = new_workspace_run(workspace_root, run_id).absolute()
    preparation = preparation_root.absolute()
    config = load_configuration(configuration_path)
    if preparation.parent != destination.parent or reference_source_dir.is_symlink() or not reference_source_dir.is_dir():
        raise ValueError("d23b_source_scope_invalid")
    verify_preparation(preparation)
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23b_preparation_anchor_mismatch")
    for tool in ("plink", "bcftools"):
        if shutil.which(tool) is None:
            raise ValueError("d23b_required_tool_missing")
    destination.mkdir(mode=0o700)
    write_json(destination / "configuration.json", config)
    target = yaml.safe_load((preparation / "blinded_inputs/target_metadata.yaml").read_text())
    base = yaml.safe_load((preparation / "blinded_inputs/pipeline.yaml").read_text())
    base["tools"].update(
        plink="plink",
        bcftools="bcftools",
        phasing_adapter={
            "adapter_id": "shapeit5_phase_common_rare_v1",
            "phase_common_command": str(shapeit_bin_dir / "SHAPEIT5_phase_common"),
            "phase_rare_command": str(shapeit_bin_dir / "SHAPEIT5_phase_rare"),
            "expected_version": "5.1.1",
        },
    )
    maps = destination / "genetic_maps"
    shutil.copytree(preparation / "reference_assets/genetic_maps", maps)
    cache = destination / "reference_cache"
    phasing_bench.prepare_reference_cache(
        reference_source_dir.absolute(), cache, Path(base["inputs"]["reference_panel_catalog"]),
        target["chromosome"], config["command_timeout_seconds"],
    )
    with gzip.open(preparation / "blinded_inputs/observed.vcf.gz", "rt") as handle:
        observed = read_panel(handle, phased=False)
    timeout = config["command_timeout_seconds"]
    common = {
        "reference_panel_id": "1kg_3202_high_coverage_20220422",
        "reference_cache_dir": str(cache),
        "require_reference_sample_selection": True,
        "shapeit5_threads": 1,
        "shapeit5_seed": config["phasing_seed"],
    }
    autosomal_receipts, regional_receipts = [], []
    for density in config["density_profiles"]:
        profile_id = density["profile_id"]
        profile = density_profile(observed, density, target)
        LOGGER.info("D23b autosomal %s", profile_id)
        autosomal_root = destination / "AUTOSOMAL" / profile_id
        _fixture_root(autosomal_root, "AUTOSOMAL", profile_id, None)
        _resolved_config(autosomal_root, base)
        inputs = phasing_bench.fixture_inputs(
            autosomal_root, preparation, profile, config, target, include_pedigree=False,
        )
        stage10a = phasing_bench.execute_component(
            autosomal_root, "10A", "prepare_autosomal_phasing_panel",
            {"minimum_total_variants": config["minimum_total_variants"],
             "minimum_variants_per_autosome": config["minimum_variants_per_autosome"]},
            inputs, timeout,
        )
        stages = [stage10a]
        if stage10a["returncode"] == 0:
            params12a = dict(
                **common,
                genetic_map_cache_dir=str(maps),
                genetic_map_id="shapeit4_hapmap_grch38_7a0cab7",
                reference_source_dir=str(reference_source_dir.absolute()),
                reference_extract_workers=2,
                tool_timeout_seconds=timeout,
                reference_extract_timeout_seconds=timeout,
                minimum_matched_variants_per_autosome=config["minimum_variants_per_autosome"],
            )
            stages.append(phasing_bench.execute_component(
                autosomal_root, "12A", "phase_autosomal_panel", params12a,
                inputs + stage10a["outputs"], timeout * 22,
            ))
        autosomal_receipts.append({
            "profile": profile_id,
            "observed_markers": len(profile.variants),
            "pedigree_mode": "WITHOUT_PEDIGREE_BY_CONTRACT",
            "stages": stages,
        })
        write_json(destination / "progress.json", {
            "autosomal": autosomal_receipts, "regional": regional_receipts, "state": "RUNNING",
        })

        for pedigree_mode in config["regional_pedigree_modes"]:
            LOGGER.info("D23b régional %s %s", profile_id, pedigree_mode)
            regional_root = destination / "REGIONAL" / profile_id / pedigree_mode
            _fixture_root(regional_root, "REGIONAL", profile_id, pedigree_mode)
            _resolved_config(regional_root, base)
            regional_inputs = phasing_bench.fixture_inputs(
                regional_root, preparation, profile, config, target,
                include_pedigree=pedigree_mode == "WITH_PEDIGREE",
            )
            params11 = {
                "window_left_bp": config["target_flank_bp"],
                "window_right_bp": config["target_flank_bp"],
                "genetic_map_id": "shapeit4_hapmap_grch38_7a0cab7",
                "genetic_map_cache_dir": str(maps),
                "genetic_map_cache_offline": True,
            }
            stage11 = phasing_bench.execute_component(
                regional_root, "11", "prepare_target_region", params11,
                [row for row in regional_inputs if row["artifact_id"] != "config_input_phasing_reference_selection"],
                timeout,
            )
            regional_stages = [stage11]
            if stage11["returncode"] == 0:
                params12 = dict(
                    **common,
                    reference_cache_offline=True,
                    minimum_common_variants=config["minimum_variants_per_autosome"],
                    shapeit5_input_timeout_seconds=timeout,
                    shapeit5_execution_timeout_seconds=timeout,
                )
                regional_stages.append(phasing_bench.execute_component(
                    regional_root, "12", "phase_target_region", params12,
                    regional_inputs + stage11["outputs"], timeout,
                ))
            regional_receipts.append({
                "profile": profile_id,
                "observed_markers": len(profile.variants),
                "pedigree_mode": pedigree_mode,
                "stages": regional_stages,
            })
            write_json(destination / "progress.json", {
                "autosomal": autosomal_receipts, "regional": regional_receipts, "state": "RUNNING",
            })
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23b_preparation_changed_during_run")
    receipt = {
        "schema_version": "1.0.0",
        "method_id": METHOD_ID,
        "closure_state": "CLOSED",
        "preparation_manifest_sha256": config["preparation_manifest_sha256"],
        "configuration_sha256": sha256_file(configuration_path),
        "scope": "ISOLATED_COMPONENTS_10A_11_12_12A_NOT_04_09",
        "fixture_producer": phasing_bench.PRODUCER,
        "implementation_sha256": sha256_file(Path(__file__)),
        "autosomal": autosomal_receipts,
        "regional": regional_receipts,
        "artifacts": file_inventory(destination, ("run_index.json",)),
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(receipt, "phasing_sensitivity_run.schema.json", schemas_dir=_schemas_dir())
    write_json(destination / "run_index.json", receipt)
    for path in destination.rglob("*"):
        if path.is_file() and path.stat().st_mode & 0o222:
            path.chmod(0o600)
    return receipt


def finalize_bench(*, run_root: Path, preparation_root: Path, configuration_path: Path) -> dict[str, Any]:
    """Ferme un run complet dont seule la publication finale a échoué.

    Chaque composant et sa liste de sorties sont revérifiés. Cette reprise ne
    lance aucun outil génétique et refuse un reçu déjà présent ou incomplet.
    """
    root, preparation = run_root.absolute(), preparation_root.absolute()
    if not root.is_dir() or (root / "run_index.json").exists():
        raise ValueError("d23b_finalize_run_state_invalid")
    config = load_configuration(configuration_path)
    verify_preparation(preparation)
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23b_finalize_preparation_anchor_mismatch")
    if json.loads((root / "configuration.json").read_text()) != config:
        raise ValueError("d23b_finalize_configuration_mismatch")
    progress = json.loads((root / "progress.json").read_text())
    autosomal, regional = progress.get("autosomal"), progress.get("regional")
    if not isinstance(autosomal, list) or len(autosomal) != 4 or not isinstance(regional, list) or len(regional) != 8:
        raise ValueError("d23b_finalize_progress_incomplete")
    expected_autosomal = [(row["profile_id"], "WITHOUT_PEDIGREE_BY_CONTRACT") for row in config["density_profiles"]]
    expected_regional = [(row["profile_id"], mode) for row in config["density_profiles"] for mode in config["regional_pedigree_modes"]]
    if [(row["profile"], row["pedigree_mode"]) for row in autosomal] != expected_autosomal:
        raise ValueError("d23b_finalize_autosomal_cells_invalid")
    if [(row["profile"], row["pedigree_mode"]) for row in regional] != expected_regional:
        raise ValueError("d23b_finalize_regional_cells_invalid")
    cells = [
        *((root / "AUTOSOMAL" / row["profile"], row, ("10A", "12A")) for row in autosomal),
        *((root / "REGIONAL" / row["profile"] / row["pedigree_mode"], row, ("11", "12")) for row in regional),
    ]
    for cell_root, entry, expected_stages in cells:
        if tuple(stage["stage_id"] for stage in entry["stages"]) != expected_stages:
            raise ValueError("d23b_finalize_stage_set_invalid")
        for stage in entry["stages"]:
            if stage["returncode"] != 0:
                raise ValueError("d23b_finalize_component_failed")
            outputs = checked_stage(cell_root, stage["directory"])
            if list(outputs.values()) != stage["outputs"]:
                raise ValueError("d23b_finalize_component_outputs_changed")
    receipt = {
        "schema_version": "1.0.0",
        "method_id": METHOD_ID,
        "closure_state": "CLOSED",
        "preparation_manifest_sha256": config["preparation_manifest_sha256"],
        "configuration_sha256": sha256_file(configuration_path),
        "scope": "ISOLATED_COMPONENTS_10A_11_12_12A_NOT_04_09",
        "fixture_producer": phasing_bench.PRODUCER,
        "implementation_sha256": sha256_file(Path(__file__)),
        "autosomal": autosomal,
        "regional": regional,
        "artifacts": file_inventory(root, ("run_index.json",)),
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(receipt, "phasing_sensitivity_run.schema.json", schemas_dir=_schemas_dir())
    write_json(root / "run_index.json", receipt)
    for path in root.rglob("*"):
        if path.is_file() and path.stat().st_mode & 0o222:
            path.chmod(0o600)
    return receipt


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    parser = argparse.ArgumentParser(description="Banc D23b densité et pedigree")
    for name in ("preparation-root", "configuration-path", "workspace-root", "reference-source-dir", "shapeit-bin-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    result = run_bench(**vars(parser.parse_args()))
    print({"state": result["closure_state"], "autosomal": len(result["autosomal"]), "regional": len(result["regional"])})
    return 0


def finalize_main() -> int:
    """CLI de reprise bornée après échec de la seule publication finale."""
    parser = argparse.ArgumentParser(description="Fermeture D23b sans relance génétique")
    for name in ("run-root", "preparation-root", "configuration-path"):
        parser.add_argument("--" + name, type=Path, required=True)
    result = finalize_bench(**vars(parser.parse_args()))
    print({"state": result["closure_state"], "autosomal": len(result["autosomal"]), "regional": len(result["regional"])})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
