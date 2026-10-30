"""Reprise régionale D23 : correction du seul adaptateur de contrats amont.

Les essais d'origine restent immuables. Les composants 11/12 sont exécutés
réellement et séparément ; aucun résultat 10A/12A n'est réinventé ou réécrit.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import logging
from pathlib import Path
import shutil
from typing import Any

import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.phasing_bench import CASES, METHOD_ID, execute_component, load_configuration
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.qc_bench import write_table
from effet_fondateur.calibration.evaluate_technical import raw_table
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document, validate_tsv_table


def bind_dataset_artifacts(artifacts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Renseigne les ensembles depuis les descripteurs, sans modifier les données."""
    result = deepcopy(artifacts)
    for prefix in ("genomewide_pre_qc", "target_chromosome_all_qc"):
        descriptor_artifact = next(item for item in result if item["artifact_id"] == prefix + "_dataset")
        path = Path(descriptor_artifact["path"])
        if sha256_file(path) != descriptor_artifact["sha256"]:
            raise ValueError("d23_fixture_descriptor_changed")
        descriptor = json.loads(path.read_text())
        validate_json_document(descriptor, "plink_dataset.schema.json")
        for item in result:
            if item["artifact_id"] in {prefix + "_" + suffix for suffix in ("bed", "bim", "fam", "dataset")}:
                for key in ("sample_set_id", "variant_set_id"):
                    item[key] = descriptor[key]
                if not item["artifact_id"].endswith("_dataset"):
                    suffix = item["artifact_id"].rsplit("_", 1)[1]
                    if descriptor["files"][suffix] != dict(name=Path(item["path"]).name, sha256=item["sha256"]):
                        raise ValueError("d23_fixture_descriptor_artifact_mismatch")
    return result


def canonical_genotype_audit(artifacts: list[dict[str, Any]], destination: Path) -> list[dict[str, Any]]:
    """Réordonne uniquement les colonnes ; refuse tout champ manquant ou supplémentaire."""
    result = deepcopy(artifacts)
    item = next(row for row in result if row["artifact_id"] == "target_genotype_audit")
    rows = raw_table(Path(item["path"]))
    columns = ("SAMPLE_ID", "FID", "IID", "GENOTYPE", "GENOTYPE_SOURCE", "MEASUREMENT_METHOD", "STATUS")
    if not rows or any(set(row) != set(columns) for row in rows):
        raise ValueError("d23_genotype_audit_columns_invalid")
    destination.mkdir()
    output = destination / "target_genotype_audit.tsv"
    write_table(output, columns, rows)
    validate_tsv_table(output, "target_genotype_audit.schema.json")
    if raw_table(output) != rows:
        raise ValueError("d23_genotype_audit_values_changed")
    item.update(path=str(output), sha256=sha256_file(output))
    return result


def run_profile(*, source_root: Path, configuration_path: Path, workspace_root: Path, run_id: str, profile: str) -> dict[str, Any]:
    """Nouveau reçu régional ancré à l'essai initial, sans écrasement ni QC fictif."""
    if profile not in CASES:
        raise ValueError("d23_regional_profile_unknown")
    destination = new_workspace_run(workspace_root, run_id).absolute()
    source_root = source_root.absolute()
    if source_root.parent != destination.parent:
        raise ValueError("d23_regional_source_scope_invalid")
    source = source_root / profile
    anchor = source / "stages/10A_prepare_autosomal_phasing_panel/stage_inputs.json"
    inputs = json.loads(anchor.read_text())
    validate_json_document(inputs, "stage_inputs.schema.json")
    for artifact in inputs["artifacts"]:
        path = Path(artifact["path"])
        if path.is_symlink() or not path.is_file() or sha256_file(path) != artifact["sha256"]:
            raise ValueError("d23_regional_input_changed")
    artifacts = bind_dataset_artifacts(inputs["artifacts"])
    config = load_configuration(configuration_path)
    destination.mkdir(mode=0o700)
    (destination / "stages").mkdir()
    artifacts = canonical_genotype_audit(artifacts, destination / "stages/D23_CANONICAL_AUDIT_FIXTURE")
    for item in artifacts:
        if item["artifact_id"] in {"target_genotype_audit", "samples_master", "variant_qc_final"}:
            validate_tsv_table(Path(item["path"]), item["artifact_id"] + ".schema.json")
    shutil.copyfile(source / "config.resolved.yaml", destination / "config.resolved.yaml")
    load_pipeline_config(destination / "config.resolved.yaml")
    write_json(destination / "manifest.json", dict(run_id="d23_regional_" + profile.lower(), scope="COMPONENT_BENCH_11_12"))
    write_json(destination / "configuration.json", config)
    timeout = config["command_timeout_seconds"]
    params11 = dict(window_left_bp=config["target_flank_bp"], window_right_bp=config["target_flank_bp"], genetic_map_id="shapeit4_hapmap_grch38_7a0cab7", genetic_map_cache_dir=str(source_root / "genetic_maps"), genetic_map_cache_offline=True)
    steps = [execute_component(destination, "11", "prepare_target_region", params11, [item for item in artifacts if item["artifact_id"] != "config_input_phasing_reference_selection"], timeout)]
    if steps[0]["returncode"] == 0:
        params12 = dict(reference_panel_id="1kg_3202_high_coverage_20220422", reference_cache_dir=str(source_root / "reference_cache"),
                        reference_cache_offline=True, require_reference_sample_selection=True, shapeit5_threads=1, shapeit5_seed=config["phasing_seed"],
                        minimum_common_variants=config["minimum_variants_per_autosome"], shapeit5_input_timeout_seconds=timeout, shapeit5_execution_timeout_seconds=timeout)
        # Les chemins des sorties régionales sont relatifs à ce nouveau reçu.
        steps.append(execute_component(destination, "12", "phase_target_region", params12, artifacts + steps[0]["outputs"], timeout))
    receipt = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", profile=profile, source_profile_input_sha256=sha256_file(anchor),
                   preparation_manifest_sha256=config["preparation_manifest_sha256"], configuration_sha256=sha256_file(configuration_path),
                   implementation_sha256=sha256_file(Path(__file__)), correction="FIXTURE_SET_IDENTITIES_AND_CANONICAL_TSV_ORDER_NO_GENOTYPE_CHANGE", stages=steps,
                   artifacts=file_inventory(destination, ("run_index.json",)), calibrated=False, founder_effect_proven=False)
    write_json(destination / "run_index.json", receipt)
    return receipt


def main() -> int:
    """CLI d'une reprise régionale privée, à un profil explicitement nommé."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    parser = argparse.ArgumentParser(description="Reprise régionale D23 du contrat de fixture")
    for name in ("source-root", "configuration-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--profile", choices=CASES, required=True)
    receipt = run_profile(**vars(parser.parse_args()))
    print(dict(state=receipt["closure_state"], profile=receipt["profile"], codes=[row["returncode"] for row in receipt["stages"]]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
