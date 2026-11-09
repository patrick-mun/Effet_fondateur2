"""Évaluation fermée D23b, distincte des processus ayant exécuté SHAPEIT5."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.evaluate_phasing import (
    component_paths,
    evaluate_autosomal_input,
    evaluate_phase_outputs,
    evaluate_region,
    phase_subset,
)
from effet_fondateur.calibration.evaluate_technical import contained_file, raw_table
from effet_fondateur.calibration.phasing_bench import PRODUCER
from effet_fondateur.calibration.phasing_oracle import parse_vcf
from effet_fondateur.calibration.phasing_preparation import verify_preparation
from effet_fondateur.calibration.phasing_sensitivity_bench import (
    METHOD_ID,
    density_profile,
    load_configuration,
)
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.technical_evaluation import compare_genotypes, observation_matrix, read_plink
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document

METRICS = (
    "haplotype_hamming_rate",
    "isolated_point_motif_rate",
    "rare_variant_phase_error_rate",
    "switch_error_rate",
    "target_mutation_orientation_accuracy",
)


def checked_closure(root: Path) -> dict[str, Any]:
    """Vérifie le reçu fermé et l'inventaire exhaustif du banc D23b."""
    document = json.loads(contained_file(root, "run_index.json").read_text())
    if (document.get("schema_version"), document.get("method_id"), document.get("closure_state")) != ("1.0.0", METHOD_ID, "CLOSED"):
        raise ValueError("d23b_closure_invalid")
    if document.get("calibrated") is not False or document.get("founder_effect_proven") is not False:
        raise ValueError("d23b_verdict_invalid")
    artifacts = document.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts or file_inventory(root, ("run_index.json",)) != artifacts:
        raise ValueError("d23b_inventory_invalid")
    for item in artifacts:
        if sha256_file(contained_file(root, item["path"])) != item["sha256"]:
            raise ValueError("d23b_closed_artifact_changed")
    return document


def paired_rate_contrast(first: list[dict[str, int]], second: list[dict[str, int]], *, seed: int, replicates: int) -> dict[str, Any]:
    """Différence de ratios avec rééchantillonnage apparié des unités."""
    if len(first) != len(second) or len(first) < 1 or replicates < 1:
        raise ValueError("d23b_paired_counts_invalid")
    first_values = np.array([[row["numerator"], row["denominator"]] for row in first], dtype=np.int64)
    second_values = np.array([[row["numerator"], row["denominator"]] for row in second], dtype=np.int64)
    if np.any(first_values < 0) or np.any(second_values < 0) or np.any(first_values[:, 0] > first_values[:, 1]) or np.any(second_values[:, 0] > second_values[:, 1]):
        raise ValueError("d23b_paired_counts_invalid")
    denominators = (int(first_values[:, 1].sum()), int(second_values[:, 1].sum()))
    result = {
        "first_numerator": int(first_values[:, 0].sum()),
        "first_denominator": denominators[0],
        "second_numerator": int(second_values[:, 0].sum()),
        "second_denominator": denominators[1],
        "estimate": None,
        "interval_95": None,
        "bootstrap_valid_replicates": 0,
        "bootstrap_excluded_zero_denominator": replicates,
        "interval_method": "PAIRED_FAMILY_CLUSTER_PERCENTILE_DESCRIPTIVE",
    }
    if 0 in denominators:
        result["status"] = "NOT_EVALUABLE"
        return result
    result["estimate"] = first_values[:, 0].sum() / denominators[0] - second_values[:, 0].sum() / denominators[1]
    indices = np.random.default_rng(seed).integers(len(first), size=(replicates, len(first)))
    first_num = first_values[indices, 0].sum(axis=1)
    first_den = first_values[indices, 1].sum(axis=1)
    second_num = second_values[indices, 0].sum(axis=1)
    second_den = second_values[indices, 1].sum(axis=1)
    valid = (first_den > 0) & (second_den > 0)
    values = first_num[valid] / first_den[valid] - second_num[valid] / second_den[valid]
    result.update(
        status="EVALUABLE",
        interval_95=np.quantile(values, [.025, .975]).tolist() if len(values) else None,
        bootstrap_valid_replicates=int(valid.sum()),
        bootstrap_excluded_zero_denominator=int((~valid).sum()),
    )
    return result


def _performance(root: Path, entry: dict[str, Any], stage: str, observation: Any, truth: Any,
                 registry: list[dict[str, str]], selection: dict[str, Any], target: dict[str, Any],
                 config: dict[str, Any], preparation: Path, workspace: Path,
                 *, expected_pedigree_records: int | None = None) -> dict[str, Any]:
    steps = {row["stage_id"]: row for row in entry["stages"]}
    if stage not in steps:
        return {"status": "NOT_EVALUABLE", "reason": "UPSTREAM_COMPONENT_FAILED"}
    step = steps[stage]
    if step["returncode"]:
        return {"status": "FAIL", "returncode": step["returncode"], "metrics_status": "NOT_EVALUABLE"}
    paths = component_paths(root, step, workspace)
    source = root / "stages" / PRODUCER
    if stage == "10A":
        return evaluate_autosomal_input(source / "genomewide_pre_qc", paths["autosomal_phasing_panel_bed"].with_suffix(""))
    if stage == "11":
        original_map = next((preparation / "reference_assets/genetic_maps").rglob(f"chr{target['chromosome']}.grch38.tsv"))
        return evaluate_region(source / "target_chromosome_all_qc", paths["target_region_bed"].with_suffix(""), paths["target_genetic_map"], original_map, target, config)
    upstream_stage = "11" if stage == "12" else "10A"
    upstream = steps[upstream_stage]
    if upstream["returncode"]:
        return {"status": "NOT_EVALUABLE", "reason": "UPSTREAM_COMPONENT_FAILED"}
    upstream_paths = component_paths(root, upstream, workspace)
    bed = upstream_paths["target_region_bed" if stage == "12" else "autosomal_phasing_panel_bed"]
    candidates = read_plink(bed.with_suffix(""))
    keys = {row[:2] for row in candidates.variants}
    if stage == "12A":
        keys.discard((target["chromosome"], target["position_bp"]))
    expected_phase = phase_subset(truth, keys)
    result = evaluate_phase_outputs(
        paths, stage, expected_phase, registry, config, target,
        [sample for _, sample in candidates.samples], selection, observation, preparation,
        include_unit_counts=True,
    )
    if stage == "12" and expected_pedigree_records is not None:
        inputs_manifest = json.loads(paths["shapeit5_inputs_manifest"].read_text())
        observed_records = inputs_manifest["pedigree_record_count"]
        result["pedigree_record_count"] = observed_records
        result["pedigree_record_count_expected"] = expected_pedigree_records
        if observed_records != expected_pedigree_records:
            result["status"] = "FAIL"
            result["pedigree_record_count_error"] = True
    return result


def _check_fixture_conversion(root: Path, observation: Any, registry: list[dict[str, str]]) -> None:
    matrix = observation_matrix(observation, registry)
    converted = read_plink(root / "stages" / PRODUCER / "genomewide_pre_qc")
    if compare_genotypes(matrix, converted)["status"] != "PASS":
        raise ValueError("d23b_fixture_conversion_not_concordant")


def _contrasts(performance: dict[tuple[str, str, str], dict[str, Any]], config: dict[str, Any]) -> list[dict[str, Any]]:
    results = []
    def append(kind: str, scope: str, profile: str, mode: str, first_key: tuple[str, str, str], second_key: tuple[str, str, str]) -> None:
        for metric in METRICS:
            common = {
                "contrast_type": kind, "scope": scope, "profile": profile,
                "pedigree_mode": mode, "metric": metric, "stratum": "ALL",
            }
            if "unit_counts" not in performance[first_key] or "unit_counts" not in performance[second_key]:
                results.append({
                    **common,
                    "status": "NOT_EVALUABLE",
                    "reason": "SOURCE_CELL_NOT_EVALUABLE",
                    "estimate": None,
                    "interval_95": None,
                    "bootstrap_valid_replicates": 0,
                    "bootstrap_excluded_zero_denominator": config["bootstrap_replicates"],
                })
                continue
            first = performance[first_key]["unit_counts"][metric]["ALL"]
            second = performance[second_key]["unit_counts"][metric]["ALL"]
            results.append({
                **common,
                **paired_rate_contrast(first, second, seed=config["bootstrap_seed"], replicates=config["bootstrap_replicates"]),
            })
    for profile in ("DENSITY_25", "DENSITY_50", "DENSITY_75"):
        append("DENSITY_MINUS_DENSITY_100", "AUTOSOMAL", profile, "WITHOUT_PEDIGREE_BY_CONTRACT",
               ("AUTOSOMAL", profile, "WITHOUT_PEDIGREE_BY_CONTRACT"), ("AUTOSOMAL", "DENSITY_100", "WITHOUT_PEDIGREE_BY_CONTRACT"))
        for mode in ("WITH_PEDIGREE", "WITHOUT_PEDIGREE"):
            append("DENSITY_MINUS_DENSITY_100", "REGIONAL", profile, mode,
                   ("REGIONAL", profile, mode), ("REGIONAL", "DENSITY_100", mode))
    for profile in ("DENSITY_25", "DENSITY_50", "DENSITY_75", "DENSITY_100"):
        append("WITHOUT_MINUS_WITH_PEDIGREE", "REGIONAL", profile, "PAIRED",
               ("REGIONAL", profile, "WITHOUT_PEDIGREE"), ("REGIONAL", profile, "WITH_PEDIGREE"))
    return results


def evaluate_bench(*, preparation_root: Path, run_root: Path, configuration_path: Path,
                   workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Évalue les 12 cellules D23b et publie uniquement des agrégats."""
    destination = new_workspace_run(workspace_root, run_id).absolute()
    workspace = workspace_root.absolute()
    preparation, run_root = preparation_root.absolute(), run_root.absolute()
    config = load_configuration(configuration_path)
    preparation_manifest = verify_preparation(preparation)
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23b_evaluation_preparation_changed")
    run = checked_closure(run_root)
    if run["configuration_sha256"] != sha256_file(configuration_path):
        raise ValueError("d23b_evaluation_run_configuration_mismatch")
    with gzip.open(preparation / "private_truth/truth.phased.vcf.gz", "rt") as handle:
        truth = parse_vcf(handle)
    with gzip.open(preparation / "blinded_inputs/observed.vcf.gz", "rt") as handle:
        observed = read_panel(handle, phased=False)
    target = yaml.safe_load((preparation / "blinded_inputs/target_metadata.yaml").read_text())
    registry = raw_table(preparation / "blinded_inputs/metadata/samples.tsv")
    selection = json.loads((preparation / "blinded_inputs/reference_selection.json").read_text())
    definitions = {row["profile_id"]: row for row in config["density_profiles"]}
    performance: dict[tuple[str, str, str], dict[str, Any]] = {}
    technical: list[dict[str, Any]] = []
    for entry in run["autosomal"]:
        profile = entry["profile"]
        root = run_root / "AUTOSOMAL" / profile
        observation = density_profile(observed, definitions[profile], target)
        _check_fixture_conversion(root, observation, registry)
        stage10a = _performance(root, entry, "10A", observation, truth, registry, selection, target, config, preparation, workspace)
        stage12a = _performance(root, entry, "12A", observation, truth, registry, selection, target, config, preparation, workspace)
        technical.append({"scope": "AUTOSOMAL", "profile": profile, "pedigree_mode": entry["pedigree_mode"], "stage_id": "10A", **stage10a})
        performance[("AUTOSOMAL", profile, entry["pedigree_mode"])] = stage12a
    for entry in run["regional"]:
        profile, mode = entry["profile"], entry["pedigree_mode"]
        root = run_root / "REGIONAL" / profile / mode
        observation = density_profile(observed, definitions[profile], target)
        _check_fixture_conversion(root, observation, registry)
        stage11 = _performance(root, entry, "11", observation, truth, registry, selection, target, config, preparation, workspace)
        stage12 = _performance(
            root, entry, "12", observation, truth, registry, selection, target,
            config, preparation, workspace,
            expected_pedigree_records=8 if mode == "WITH_PEDIGREE" else 0,
        )
        technical.append({"scope": "REGIONAL", "profile": profile, "pedigree_mode": mode, "stage_id": "11", **stage11})
        performance[("REGIONAL", profile, mode)] = stage12
    destination.mkdir(mode=0o700)
    public_performance = []
    for (scope, profile, mode), result in performance.items():
        counts = result.pop("unit_counts", None)
        if counts is None and result.get("status") == "PASS":
            raise ValueError("d23b_unit_counts_missing")
        public_performance.append({"scope": scope, "profile": profile, "pedigree_mode": mode, **result})
        if counts is not None:
            result["unit_counts"] = counts
    contrasts = _contrasts(performance, config)
    write_json(destination / "technical_evaluations.json", technical)
    write_json(destination / "performance_evaluations.json", public_performance)
    write_json(destination / "paired_contrasts.json", contrasts)
    receipt = {
        "schema_version": "1.0.0",
        "method_id": "paired_density_pedigree_phasing_evaluation_v1",
        "closure_state": "CLOSED",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "preparation_manifest_sha256": config["preparation_manifest_sha256"],
        "source_run_sha256": sha256_file(run_root / "run_index.json"),
        "configuration_sha256": sha256_file(configuration_path),
        "implementation_sha256": sha256_file(Path(__file__)),
        "oracle_implementation_sha256": sha256_file(Path(__file__).with_name("phasing_oracle.py")),
        "prespecification_sha256": sha256_file(Path(__file__).resolve().parents[3] / "calibration/protocoles/D23B_DENSITY_PEDIGREE_PRESPECIFICATION.md"),
        "independent_units": preparation_manifest["counts"]["independent_simulation_units"],
        "performance_qualification": "REVIEW_REQUIRED",
        "limitations": [
            "SYNTHETIC_TECHNICAL_TARGET_NOT_DOCK6",
            "SINGLE_PANEL_SINGLE_PHASING_SOFTWARE",
            "DESCRIPTIVE_PAIRED_INTERVALS_NOT_CAUSAL",
            "PERFORMANCE_THRESHOLDS_NOT_SEALED",
            "D21_SCIENTIFIC_FAIL_07_UNCHANGED",
        ],
        "artifacts": file_inventory(destination, ("run_index.json",)),
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(
        receipt, "phasing_sensitivity_evaluation_manifest.schema.json",
        schemas_dir=Path(__file__).resolve().parents[3] / "calibration/schemas",
    )
    write_json(destination / "run_index.json", receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description="Évaluation indépendante D23b")
    for name in ("preparation-root", "run-root", "configuration-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    result = evaluate_bench(**vars(parser.parse_args()))
    print({"state": result["closure_state"], "qualification": result["performance_qualification"]})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
