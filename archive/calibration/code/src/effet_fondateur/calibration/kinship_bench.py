"""Préparation D21 : entrées publiques neutres et vérité privée séparées."""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np
import yaml

from effet_fondateur.calibration.evaluate_technical import contained_file, validate_schema
from effet_fondateur.calibration.export_blinded_inputs import SAMPLE_COLUMNS, _write_bundle, pipeline_configuration
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import file_inventory, write_json, write_table
from effet_fondateur.calibration.simulation.kinship_bench import KinshipSimulation, generate_simulation
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel
from effet_fondateur.calibration.technical_evaluation import TechnicalEvaluationError
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config
from effet_fondateur.contracts.samples import validate_samples_master


METHOD_ID = "iid_pedigree_kinship_bench_v1"
STAGES = ("initialize_run", "validate_sources", "build_sample_registry", "convert_acpa", "qc_preliminary", "build_kinship_panel", "infer_kinship")
STAGE_IDS = ("00", "01", "02", "03", "05", "06", "07")


def validate_configuration(config: dict[str, Any]) -> None:
    """Contrôle schéma strict, graines distinctes et géométrie du banc borné."""
    validate_schema(config, "kinship_bench_config")
    if len({config[key] for key in ("simulation_seed", "missingness_seed", "bootstrap_seed")}) != 3 or len({row["stratum_id"] for row in config["strata"]}) != 2 or len({row["profile_id"] for row in config["profiles"]}) != 3:
        raise TechnicalEvaluationError("d21_seed_or_profile_collision")
    if max(config["ld_copy_offsets"] + [config["monomorphic_offset"]]) >= config["markers_per_autosome"] or config["monomorphic_offset"] in config["ld_copy_offsets"]:
        raise TechnicalEvaluationError("d21_marker_interventions_overlap")
    if config["panel_parameters"]["required_autosomes"] != list(range(1, 23)):
        raise TechnicalEvaluationError("d21_autosome_coverage_not_prespecified")
    profiles = [(row["profile_id"], row["missingness_fraction"], row["marker_stride"]) for row in config["profiles"]]
    if profiles != [("COMPLETE", 0, 1), ("MISSING_5PCT", .05, 1), ("DENSITY_25PCT", 0, 4)]:
        raise TechnicalEvaluationError("d21_profile_semantics_changed")


def checked_document(root: Path, name: str, state_key: str, expected_state: str) -> dict[str, Any]:
    """Refuse clôture manquante, chemins non confinés et inventaire incomplet."""
    path = contained_file(root, name)
    if not path.is_file():
        raise TechnicalEvaluationError("d21_closure_missing")
    document = _json(path)
    validate_schema(document, "kinship_bench_document")
    if document.get("method_id") != METHOD_ID or document.get(state_key) != expected_state:
        raise TechnicalEvaluationError("d21_state_invalid")
    artifacts = document.get("artifacts", [])
    if not artifacts or len({row["path"] for row in artifacts}) != len(artifacts):
        raise TechnicalEvaluationError("d21_inventory_invalid")
    for artifact in artifacts:
        source = contained_file(root, artifact["path"])
        if not source.is_file() or _sha(source) != artifact["sha256"]:
            raise TechnicalEvaluationError("d21_artifact_changed")
    if file_inventory(root, (name,)) != artifacts:
        raise TechnicalEvaluationError("d21_uninventoried_file")
    return document


def profile_panel(simulation: KinshipSimulation, config: dict[str, Any], profile_index: int) -> tuple[MeasurementPanel, np.ndarray]:
    """Applique densité et manque sans modifier la vérité de transmission."""
    profile = config["profiles"][profile_index]
    marker_count, sample_count, _ = simulation.alleles.shape
    selected = np.arange(0, marker_count, profile["marker_stride"])
    alleles = np.sort(simulation.alleles[selected], axis=2)
    missing = np.random.default_rng(config["missingness_seed"] + profile_index).random(alleles.shape[:2]) < profile["missingness_fraction"]
    alleles[missing] = -1
    samples = tuple("S" + hashlib.sha256(f'{METHOD_ID}:{config["simulation_seed"]}:{profile_index}:{index}'.encode()).hexdigest()[:16] for index in range(sample_count))
    per_chromosome = config["markers_per_autosome"]
    variants = tuple((str(int(index) // per_chromosome + 1), str(config["marker_start_bp"] + int(index) % per_chromosome * config["marker_spacing_bp"]),
                      f"kin_chr{int(index)//per_chromosome+1}_{int(index)%per_chromosome:04d}", "A", "G", ".", "PASS", ".", "GT") for index in selected)
    return MeasurementPanel(samples, variants, alleles), selected


def public_registry(panel: MeasurementPanel, simulation: KinshipSimulation, config: dict[str, Any], profile_index: int) -> list[dict[str, str]]:
    """Publie seulement les liens PC observables, jamais un degré ou une strate."""
    rows = [dict(zip(SAMPLE_COLUMNS, (sample, f"{sample}.txt", sample, sample, "0", "0", "UNKNOWN", "UNKNOWN", "UNSPECIFIED", "", "", "", "true", "true", ""), strict=True)) for sample in panel.samples]
    for unit in simulation.units:
        family = "F" + hashlib.sha256(f'{METHOD_ID}:family:{config["simulation_seed"]}:{profile_index}:{unit["unit_index"]}'.encode()).hexdigest()[:16]
        for member, rank in enumerate(unit["sample_indices"]):
            rows[rank].update(FID=family if not unit["cryptic"] else family + str(member), INCLUDE_GENOMEWIDE="true", ARRAY_BATCH="LOT_A" if rank % 2 == 0 else "LOT_B")
        if unit["relation"] == "PARENT_CHILD":
            parent, child = unit["sample_indices"]
            rows[parent]["SEX"] = "MALE"
            rows[child]["PID"] = panel.samples[parent]
    return rows


def bench_configuration(config: dict[str, Any]) -> dict[str, Any]:
    """Profil réel 00–03/05–07, sans mutation ni cohorte gelée fictive."""
    result = pipeline_configuration(dict(chromosome=5, position_bp=config["marker_start_bp"], ref="A", alt="G"))
    result["project"]["project_id"] = "d21_kinship_component_bench"
    result["tools"].update(config["tools"])
    result["stages"]["validate_sources"]["parameters"]["required_autosomes"] = list(range(1, 23))
    result["stages"]["build_sample_registry"]["parameters"]["manual_approval"] = deepcopy(config["manual_approval"])
    for stage, parameters in (("convert_acpa", dict(marker_mode="intersection", plink_timeout_seconds=600)),
                              ("qc_preliminary", config["preliminary_parameters"]), ("build_kinship_panel", config["panel_parameters"]),
                              ("infer_kinship", config["kinship_parameters"])):
        result["stages"][stage] = dict(enabled=True, parameters=deepcopy(parameters))
    return result


def prepare_kinship_bench(*, configuration_path: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Matérialise un banc neuf fermé ; aucune exécution scientifique n'est lancée."""
    destination = new_workspace_run(workspace_root, run_id)
    config = _json(configuration_path)
    validate_configuration(config)
    simulation = generate_simulation(config)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        private = temporary / "private_truth"
        private.mkdir(mode=0o700)
        shutil.copyfile(configuration_path, private / "configuration.json")
        np.savez_compressed(private / "gene_drop.npz", alleles=simulation.alleles, origins=simulation.origins)
        write_json(private / "pedigrees.json", dict(units=simulation.units, mendelian_error_count=simulation.mendelian_error_count))
        bundles = []
        for index, profile in enumerate(config["profiles"]):
            bundle_id = f"B{index:02d}"
            panel, selected = profile_panel(simulation, config, index)
            bundle = temporary / "blinded_inputs/bundles" / bundle_id
            _write_bundle(bundle, panel, dict(target=dict(chromosome=5, position_bp=config["marker_start_bp"], ref="A", alt="G")))
            rows = public_registry(panel, simulation, config, index)
            write_table(bundle / "metadata/samples.tsv", tuple(rows[0]), rows)
            (bundle / "pipeline.yaml").write_text(yaml.safe_dump(bench_configuration(config), sort_keys=False), encoding="utf-8")
            validate_samples_master(bundle / "metadata/samples.tsv", bundle / "acpa_samples")
            load_pipeline_config(bundle / "pipeline.yaml")
            write_json(private / f"{bundle_id}.json", dict(profile=profile, source_marker_indices=selected.tolist()))
            bundles.append(dict(bundle_id=bundle_id, sample_count=len(panel.samples), marker_count=len(panel.variants)))
        document = dict(schema_version="1.0.0", method_id=METHOD_ID, manifest_state="FINALIZED", created_at=datetime.now(timezone.utc).isoformat(),
                        configuration_sha256=_sha(configuration_path), bundles=bundles, artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        write_json(temporary / "panel_manifest.json", document)
        checked_document(temporary, "panel_manifest.json", "manifest_state", "FINALIZED")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return document


def main() -> int:
    """CLI de préparation bornée ; aucune clé individuelle sur stdout."""
    parser = argparse.ArgumentParser(description="Préparation du banc D21, sans outil externe")
    for name in ("configuration-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    document = prepare_kinship_bench(**vars(parser.parse_args()))
    print({"bundles": len(document["bundles"]), "state": document["manifest_state"]})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
