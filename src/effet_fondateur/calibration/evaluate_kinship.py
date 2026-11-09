"""Évaluation D21 après clôture : oracle 06, parenté et proposition 07.

Les vérités sont ouvertes seulement après contrôle de la clôture et des
empreintes. Aucune commande externe n'est exécutée par cet évaluateur.
Les sorties ne contiennent que des agrégats, sans sujets, génotypes ou paires.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile
from typing import Any

import numpy as np

from effet_fondateur.calibration.evaluate_qc import _paths
from effet_fondateur.calibration.evaluate_technical import _dataset, contained_file, raw_table
from effet_fondateur.calibration.kinship_bench import METHOD_ID, STAGES, STAGE_IDS, bench_configuration, checked_document, profile_panel, public_registry, validate_configuration
from effet_fondateur.calibration.kinship_panel_oracle import evaluate_panel
from effet_fondateur.calibration.kinship_stage_oracle import evaluate_kinship_stage
from effet_fondateur.calibration.materialize_ibd import _json, _sha, _text
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.qc_oracle import subset
from effet_fondateur.calibration.simulation.kinship_bench import KinshipSimulation, generate_simulation
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.technical_evaluation import GenotypeMatrix, TechnicalEvaluationError, compare_genotypes, observation_matrix
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import load_pipeline_config, validate_json_document


def _closed_run(run_root: Path, entry: dict[str, Any]) -> Path:
    """Ancre manifest, audits et sorties des sept étapes, sans approuver 07."""
    path = contained_file(run_root, entry["run_path"] + "/manifest.json")
    manifest = _json(path)
    validate_json_document(manifest, "run_manifest.schema.json")
    if _sha(path) != entry["manifest_sha256"] or entry["returncode"] != 0 or manifest["manual_decisions_required"] != ["kinship_exclusion_approval"] or entry["pending_manual_decisions"] != manifest["manual_decisions_required"]:
        raise TechnicalEvaluationError("d21_receipt_or_manual_decision_changed")
    if [stage["stage_id"] for stage in manifest["stages"]] != list(STAGE_IDS):
        raise TechnicalEvaluationError("d21_stage_scope_changed")
    for stage, stage_id, stage_name in zip(manifest["stages"], STAGE_IDS, STAGES, strict=True):
        directory = f"stages/{stage_id}_{stage_name}"
        audit = contained_file(path.parent, directory + "/audit.json")
        outputs = contained_file(path.parent, directory + "/stage_outputs.json")
        if stage["state"] != "SUCCEEDED" or stage["stage_name"] != stage_name or stage["audit_sha256"] != _sha(audit) or stage["stage_outputs_sha256"] != _sha(outputs) or _json(audit)["signature"] != stage["signature"]:
            raise TechnicalEvaluationError("d21_stage_closure_changed")
    return path.parent


def _check_sources(panel_root: Path, config: dict[str, Any]) -> KinshipSimulation:
    """Rejoue la simulation et vérifie la vérité privée, sans ouvrir le pipeline."""
    validate_configuration(config)
    simulation = generate_simulation(config)
    with np.load(panel_root / "private_truth/gene_drop.npz", allow_pickle=False) as saved:
        if set(saved.files) != {"alleles", "origins"} or not np.array_equal(saved["alleles"], simulation.alleles) or not np.array_equal(saved["origins"], simulation.origins):
            raise TechnicalEvaluationError("d21_gene_drop_replay_mismatch")
    if _json(panel_root / "private_truth/pedigrees.json") != dict(units=json.loads(json.dumps(simulation.units)), mendelian_error_count=0):
        raise TechnicalEvaluationError("d21_pedigree_truth_changed")
    return simulation


def _expected_preliminary(matrix: GenotypeMatrix, parameters: dict[str, Any]) -> GenotypeMatrix:
    """Exclusions de manque simultanées, calculées depuis les observations seules."""
    missing = matrix.alleles[:, :, 0] < 0
    return subset(matrix, missing.mean(axis=0) <= parameters["sample_missingness_max"], missing.mean(axis=1) <= parameters["variant_missingness_max"])


def evaluate_kinship(*, panel_root: Path, run_root: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Publie les deux analyses D21 mesurées, avec incertitude et limites explicites."""
    destination = new_workspace_run(workspace_root, run_id)
    panel_root, run_root = panel_root.absolute(), run_root.absolute()
    if panel_root.parent != destination.parent or run_root.parent != destination.parent:
        raise TechnicalEvaluationError("d21_evaluation_roots_invalid")
    panel = checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
    receipt = checked_document(run_root, "run_index.json", "closure_state", "CLOSED")
    panel_sha, run_sha = _sha(panel_root / "panel_manifest.json"), _sha(run_root / "run_index.json")
    if receipt["panel_manifest_sha256"] != panel_sha or [row["bundle_id"] for row in receipt["runs"]] != [f"B{index:02d}" for index in range(3)] or [row["bundle_id"] for row in panel["bundles"]] != [row["bundle_id"] for row in receipt["runs"]]:
        raise TechnicalEvaluationError("d21_case_or_panel_binding_invalid")
    config = _json(panel_root / "private_truth/configuration.json")
    if _sha(panel_root / "private_truth/configuration.json") != panel["configuration_sha256"]:
        raise TechnicalEvaluationError("d21_configuration_binding_invalid")
    simulation = _check_sources(panel_root, config)
    evaluations = {"06": [], "07": []}
    tool_versions = set()
    for index, entry in enumerate(receipt["runs"]):
        bundle_id = entry["bundle_id"]
        source = panel_root / "blinded_inputs/bundles" / bundle_id
        if load_pipeline_config(source / "pipeline.yaml") != bench_configuration(config):
            raise TechnicalEvaluationError("d21_public_configuration_changed")
        expected_panel, source_indices = profile_panel(simulation, config, index)
        if panel["bundles"][index]["sample_count"] != len(expected_panel.samples) or panel["bundles"][index]["marker_count"] != len(expected_panel.variants):
            raise TechnicalEvaluationError("d21_panel_count_mismatch")
        truth = _json(panel_root / "private_truth" / f"{bundle_id}.json")
        if truth != dict(profile=config["profiles"][index], source_marker_indices=source_indices.tolist()):
            raise TechnicalEvaluationError("d21_profile_truth_changed")
        with _text(source / "observed.vcf.gz") as handle:
            observed = read_panel(handle, phased=False)
        if observed.samples != expected_panel.samples or observed.variants != expected_panel.variants or not np.array_equal(observed.alleles, expected_panel.alleles):
            raise TechnicalEvaluationError("d21_observation_replay_mismatch")
        registry = raw_table(source / "metadata/samples.tsv")
        if registry != public_registry(expected_panel, simulation, config, index):
            raise TechnicalEvaluationError("d21_non_neutral_public_registry")
        matrix = observation_matrix(observed, registry)
        pipeline = _closed_run(run_root, entry)
        conversion = _dataset(_paths(pipeline, "stages/03_convert_acpa"), "genomewide_base")
        if compare_genotypes(matrix, conversion)["status"] != "PASS":
            raise TechnicalEvaluationError("d21_upstream_conversion_failed")
        pre_artifacts = _paths(pipeline, "stages/05_qc_preliminary")
        preliminary = _dataset(pre_artifacts, "genomewide_pre_qc")
        if compare_genotypes(_expected_preliminary(matrix, config["preliminary_parameters"]), preliminary)["status"] != "PASS":
            raise TechnicalEvaluationError("d21_upstream_filtering_failed")
        panel_artifacts = _paths(pipeline, "stages/06_build_kinship_panel")
        evaluation_06 = evaluate_panel(preliminary, registry, raw_table(pre_artifacts["qc_variant_metrics"]), panel_artifacts, config["panel_parameters"], config)
        panel_matrix = _dataset(panel_artifacts, "kinship_panel")
        source_markers = {variant[2]: source_indices[position] for position, variant in enumerate(matrix.variants)}
        source_samples = {sample: position for position, sample in enumerate(expected_panel.samples)}
        panel_origins = simulation.origins[
            np.array([source_markers[variant[2]] for variant in panel_matrix.variants])[:, None],
            np.array([source_samples[sample] for _, sample in panel_matrix.samples])[None, :],
        ]
        evaluation_07 = evaluate_kinship_stage(panel_matrix, registry, raw_table(pre_artifacts["qc_individual_metrics"]),
                                              _paths(pipeline, "stages/07_infer_kinship"), simulation.units, expected_panel.samples, config,
                                              origins=panel_origins, before_pruning=preliminary)
        context = dict(input_markers=len(matrix.variants), post_qc_markers=len(preliminary.variants),
                       post_ld_markers=len(panel_matrix.variants), observed_missingness=float((matrix.alleles[:, :, 0] < 0).mean()))
        for stage, evaluation in (("06", evaluation_06), ("07", evaluation_07)):
            evaluations[stage].append(dict(bundle_id=bundle_id, profile=config["profiles"][index], context=context, **evaluation))
        for stage in ("06_build_kinship_panel", "07_infer_kinship"):
            audit = _json(pipeline / "stages" / stage / "audit.json")
            tool_versions.update(f'{tool["tool"]}: {tool["version"]}' for tool in audit["tools"])
    documents = {}
    for stage, results in evaluations.items():
        documents[stage] = dict(schema_version="1.0.0", method_id=METHOD_ID, stage_id=stage, report_id="CAL-REP-" + stage,
                                status="FAIL" if any(row["status"] == "FAIL" for row in results) else "PASS", evaluable=True,
                                evaluations=results, performance_decision="REVIEW_REQUIRED", calibrated=False, founder_effect_proven=False)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=destination.parent))
    try:
        for stage, document in documents.items():
            write_json(temporary / f"CAL-REP-{stage}.json", document)
        manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, evaluation_state="CLOSED", created_at=datetime.now(timezone.utc).isoformat(),
                        configuration=config, tool_versions=sorted(tool_versions), panel_manifest_sha256=panel_sha, run_index_sha256=run_sha,
                        stages=[dict(stage_id=stage, status=document["status"], case_count=len(document["evaluations"]),
                                     fail_count=sum(row["status"] == "FAIL" for row in document["evaluations"])) for stage, document in documents.items()],
                        artifacts=file_inventory(temporary), calibrated=False, founder_effect_proven=False)
        checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        checked_document(run_root, "run_index.json", "closure_state", "CLOSED")
        if _sha(panel_root / "panel_manifest.json") != panel_sha or _sha(run_root / "run_index.json") != run_sha:
            raise TechnicalEvaluationError("d21_sources_changed_during_evaluation")
        write_json(temporary / "evaluation_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return manifest


def main() -> int:
    """CLI de l'évaluateur clos, sans appel d'outil génétique."""
    parser = argparse.ArgumentParser(description="Évaluation D21 06/07 après clôture")
    for name in ("panel-root", "run-root", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    manifest = evaluate_kinship(**vars(parser.parse_args()))
    print(manifest["stages"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
