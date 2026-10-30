"""D23 : concordance indépendante et performances après clôture des outils.

Seul ce processus ouvre la vérité. Les rapports ne reçoivent que des agrégats.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.evaluate_technical import contained_file, raw_table
from effet_fondateur.calibration.phasing_bench import CASES, METHOD_ID, PRODUCER, load_configuration, observed_profile
from effet_fondateur.calibration.phasing_oracle import PhasePanel, align_prediction, parse_vcf, score_phase
from effet_fondateur.calibration.phasing_preparation import verify_preparation
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.qc_oracle import subset
from effet_fondateur.calibration.run_qc_bench import checked_stage
from effet_fondateur.calibration.simulation.measurement_errors import read_panel
from effet_fondateur.calibration.technical_evaluation import compare_genotypes, observation_matrix, read_plink
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import DocumentValidationError, validate_json_document


STAGES = ("10A", "11", "12", "12A")


def validate_evaluation(document: dict[str, Any], *, manifest: bool = False) -> None:
    """Contrats fermés : interdit notamment l'ajout de données individuelles au rendu."""
    schemas = Path(__file__).resolve().parents[3] / "calibration/schemas"
    name = "phasing_evaluation_manifest.schema.json" if manifest else "phasing_evaluation.schema.json"
    try:
        validate_json_document(document, name, schemas_dir=schemas)
    except DocumentValidationError:
        raise ValueError("d23_aggregate_evaluation_schema_invalid") from None


def checked_closure(root: Path) -> dict[str, Any]:
    """Refuse fichiers ajoutés, manquants, liens et reçus non clos avant analyse."""
    path = contained_file(root, "run_index.json")
    document = json.loads(path.read_text())
    if document.get("schema_version") != "1.0.0" or document.get("method_id") != METHOD_ID or document.get("closure_state") != "CLOSED" or document.get("calibrated") is not False or document.get("founder_effect_proven") is not False:
        raise ValueError("d23_closure_invalid")
    artifacts = document.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts or len({row["path"] for row in artifacts}) != len(artifacts):
        raise ValueError("d23_inventory_invalid")
    for item in artifacts:
        if set(item) != {"path", "sha256"} or sha256_file(contained_file(root, item["path"])) != item["sha256"]:
            raise ValueError("d23_closed_artifact_changed")
    if file_inventory(root, ("run_index.json",)) != artifacts:
        raise ValueError("d23_uninventoried_file")
    return document


def component_paths(root: Path, entry: dict[str, Any], workspace: Path) -> dict[str, Path]:
    """Relie enveloppe signée, entrées physiques, audit et sorties du composant."""
    directory = entry["directory"]
    inputs = json.loads(contained_file(root, directory + "/stage_inputs.json").read_text())
    validate_json_document(inputs, "stage_inputs.schema.json")
    signature = hashlib.sha256(json.dumps(dict(parameters=inputs["parameters"], artifacts=inputs["artifacts"]), sort_keys=True).encode()).hexdigest()
    if inputs["signature"] != signature or inputs["published_output_dir"] != directory or inputs["stage_id"] != entry["stage_id"]:
        raise ValueError("d23_component_binding_invalid")
    for artifact in inputs["artifacts"]:
        path = Path(artifact["path"])
        path = root / path if not path.is_absolute() else path
        if not path.absolute().is_relative_to(workspace.absolute()):
            # Seuls les trois catalogues publics versionnés sont hors workspace.
            allowed = Path(__file__).resolve().parents[3] / "config/references"
            if path.parent != allowed or path.suffix != ".json":
                raise ValueError("d23_component_input_outside_scope")
            if path.is_symlink():
                raise ValueError("d23_component_input_symlink")
        else:
            contained_file(workspace, path.relative_to(workspace).as_posix())
        if sha256_file(path) != artifact["sha256"]:
            raise ValueError("d23_component_input_changed")
    if entry["returncode"]:
        return {}
    outputs = checked_stage(root, directory)
    if list(outputs.values()) != entry["outputs"]:
        raise ValueError("d23_receipt_output_mismatch")
    for name in ("stage_outputs.json", "audit.json"):
        value = json.loads(contained_file(root, directory + "/" + name).read_text())
        if any(value[key] != inputs[key] for key in ("run_id", "stage_id", "stage_name", "signature")):
            raise ValueError("d23_component_output_binding_invalid")
    return {key: contained_file(root, item["path"]) for key, item in outputs.items()}


def evaluate_autosomal_input(source_prefix: Path, actual_prefix: Path) -> dict[str, Any]:
    """Attente 10A par comptage des allèles, indépendante du filtre PLINK."""
    source, actual = read_plink(source_prefix), read_plink(actual_prefix)
    complete = np.all(source.alleles >= 0, axis=(1, 2))
    fam = [line.split() for line in source_prefix.with_suffix(".fam").read_text().splitlines()]
    founders = np.array([row[2] == row[3] == "0" for row in fam])
    if not np.any(founders):
        raise ValueError("d23_10a_founders_missing")
    # PLINK exclut les non-fondateurs des fréquences, sauf --nonfounders.
    founder_calls = source.alleles[:, founders]
    maximum = founder_calls.max(axis=(1, 2))
    minimum = founder_calls.min(axis=(1, 2))
    minor_count = np.minimum((founder_calls == maximum[:, None, None]).sum(axis=(1, 2)), (founder_calls == minimum[:, None, None]).sum(axis=(1, 2)))
    keep = complete & (maximum != minimum) & (minor_count / (2 * int(founders.sum())) >= .000001)
    expected = subset(source, np.ones(len(source.samples), dtype=bool), keep)
    comparison = compare_genotypes(expected, actual, stage="10A")
    counts = [dict(chromosome=chromosome, expected=sum(row[0] == chromosome for row in expected.variants), observed=sum(row[0] == chromosome for row in actual.variants)) for chromosome in range(1, 23)]
    return dict(status=comparison["status"], concordance=comparison, source_markers=len(source.variants), expected_markers=int(keep.sum()),
                observed_markers=len(actual.variants), excluded_incomplete=int((~complete).sum()), excluded_monomorphic=int((complete & (maximum == minimum)).sum()), chromosome_counts=counts)


def evaluate_region(source_prefix: Path, actual_prefix: Path, map_path: Path, original_map: Path, target: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    """Attente 11 : fenêtre physique, dosages et interpolation NumPy indépendante."""
    source, actual = read_plink(source_prefix), read_plink(actual_prefix)
    keep = np.array([row[0] == target["chromosome"] and abs(row[1]-target["position_bp"]) <= config["target_flank_bp"] for row in source.variants])
    expected = subset(source, np.ones(len(source.samples), dtype=bool), keep)
    comparison = compare_genotypes(expected, actual, stage="11")
    reference = raw_table(original_map)
    bp = np.array([int(row["POSITION_BP"]) for row in reference])
    cm = np.array([float(row["POSITION_CM"]) for row in reference])
    if np.any(np.diff(bp) <= 0) or np.any(np.diff(cm) < 0):
        raise ValueError("d23_oracle_map_invalid")
    map_rows = raw_table(map_path)
    actual_positions = np.array([int(row["POSITION_BP"]) for row in map_rows])
    wanted = np.array([row[1] for row in expected.variants])
    observed_cm = np.array([float(row["POSITION_CM"]) for row in map_rows])
    expected_cm = np.interp(actual_positions, bp, cm)
    residuals = np.abs(observed_cm - expected_cm)
    coordinate_errors = int(np.count_nonzero(residuals > config["map_absolute_tolerance_cm"])) + int(np.any(actual_positions < bp[0]) or np.any(actual_positions > bp[-1]))
    order_ok = actual_positions.tolist() == wanted.tolist()
    target_retained = sum(row[0] == target["chromosome"] and row[1] == target["position_bp"] for row in actual.variants) == 1
    return dict(status="PASS" if comparison["status"] == "PASS" and coordinate_errors == 0 and order_ok and target_retained else "FAIL",
                concordance=comparison, expected_markers=len(expected.variants), observed_markers=len(actual.variants), map_coordinate_error_count=coordinate_errors,
                map_order_correct=order_ok, map_max_absolute_error_cm=float(residuals.max()) if len(residuals) else None, target_retained=target_retained,
                map_plot=dict(position_bp=actual_positions.tolist(), expected_cm=expected_cm.tolist(), observed_cm=observed_cm.tolist()))


def _read_phased(path: Path, timeout: int) -> PhasePanel:
    result = subprocess.run(["bcftools", "view", "--no-version", "-Ov", str(path)], capture_output=True, text=True, timeout=timeout, check=False)
    if result.returncode:
        raise ValueError("d23_evaluation_bcf_decode_failed")
    return parse_vcf(result.stdout.splitlines())


def phase_subset(panel: PhasePanel, keys: set[tuple[int, int]]) -> PhasePanel:
    indices = [index for index, locus in enumerate(panel.loci) if locus[:2] in keys]
    if len(indices) != len(keys):
        raise ValueError("d23_truth_locus_universe_mismatch")
    return replace(panel, loci=tuple(panel.loci[index] for index in indices), alleles=panel.alleles[indices], phased=panel.phased[indices], phase_sets=panel.phase_sets[indices],
                   variant_ids=tuple(panel.variant_ids[index] for index in indices) if panel.variant_ids else ())


def acceptance_results(stage: str, result: dict[str, Any]) -> list[dict[str, Any]]:
    """Raccorde les noms et décisions aux six critères D04, sans sceller de seuil."""
    path = Path(__file__).resolve().parents[3] / "calibration/protocoles/CRITERES_ACCEPTATION.yaml"
    catalogue = yaml.safe_load(path.read_text())
    if catalogue["document"]["version"] != "0.2.0" or catalogue["document"]["status"] != "APPROVED_PRE_PILOT":
        raise ValueError("d23_acceptance_version_requires_protocol_revision")
    decisions = []
    for rule in catalogue["criteria"]:
        if not rule["id"].startswith("CAL-PHASE-") or "CAL-REP-" + stage not in rule["reports"]:
            continue
        name = rule["metric"]
        value, denominator = None, None
        if name == "autosomal_sample_or_variant_universe_mismatch_count":
            if "concordance" in result:
                value = sum(count for key, count in result["concordance"]["metrics"].items() if key.startswith(("sample_", "variant_")))
            elif "sample_order_correct" in result:
                value = int(not result["sample_order_correct"]) + int(not result["variant_order_correct"])
        elif name == "genetic_map_or_coordinate_consistency_error_count" and "map_coordinate_error_count" in result:
            value = result["map_coordinate_error_count"] + int(not result["map_order_correct"])
        elif name == "target_variant_retained":
            value = result.get("target_retained")
        elif name in result.get("metrics", {}):
            measure = result["metrics"][name]["ALL"]
            value, denominator = measure["estimate"], measure["denominator"]
        acceptance = rule["acceptance"]
        if acceptance["lifecycle"] == "pilot_then_seal" and acceptance["value"] is not None:
            raise ValueError("d23_threshold_sealing_requires_new_protocol")
        status = "NOT_EVALUABLE" if value is None else "REVIEW_REQUIRED" if acceptance["lifecycle"] == "pilot_then_seal" else "PASS" if value == acceptance["value"] else "FAIL"
        decisions.append(dict(criterion_id=rule["id"], metric=name, value=value, denominator=denominator, threshold=acceptance["value"], lifecycle=acceptance["lifecycle"], status=status))
    return decisions


def evaluate_phase_outputs(paths: dict[str, Path], stage: str, truth: PhasePanel, registry: list[dict[str, str]], config: dict[str, Any],
                           target: dict[str, Any], expected_samples: list[str], reference_selection: dict[str, Any],
                           observation: Any, preparation: Path) -> dict[str, Any]:
    """Métriques et conformité indépendantes, sans prendre les PASS publiés pour oracle."""
    files = [paths["shapeit5_final_bcf"]] if stage == "12" else [paths[f"autosomal_phased_chr{chromosome}_bcf"] for chromosome in range(1, 23)]
    panels = [_read_phased(path, config["command_timeout_seconds"]) for path in files]
    sample_order_ok = all(list(panel.samples) == expected_samples for panel in panels)
    if any(panel.samples != panels[0].samples for panel in panels):
        raise ValueError("d23_phased_chromosome_sample_order_mismatch")
    predicted = PhasePanel(panels[0].samples, tuple(locus for panel in panels for locus in panel.loci), np.concatenate([panel.alleles for panel in panels]),
                           np.concatenate([panel.phased for panel in panels]), np.concatenate([panel.phase_sets for panel in panels]), tuple(marker for panel in panels for marker in panel.variant_ids))
    units = {family: index for index, family in enumerate(sorted({row["FID"] for row in registry}))}
    rows = {row["SAMPLE_ID"]: row for row in registry}
    founders = [index for index, sample in enumerate(truth.samples) if rows[sample]["PID"] == rows[sample]["MID"] == "0"]
    unit_indices = [units[rows[sample]["FID"]] for sample in truth.samples]
    performance = score_phase(truth, predicted, founder_indices=founders, unit_indices=unit_indices, target=(target["chromosome"], target["position_bp"]),
                              rare_maf_upper_bound=config["rare_maf_upper_bound"], density_gap_bp=config["density_gap_bp"], bootstrap_seed=config["bootstrap_seed"], bootstrap_replicates=config["bootstrap_replicates"])
    observed_by_locus = {(int(row[0].removeprefix("chr")), int(row[1])): index for index, row in enumerate(observation.variants)}
    observed_calls = observation.alleles[[observed_by_locus[locus[:2]] for locus in truth.loci]]
    aligned, _, _, _ = align_prediction(truth, predicted)
    observed_missing = np.any(observed_calls < 0, axis=2)
    predicted_missing = np.any(aligned < 0, axis=2)
    preserved_errors = int(np.count_nonzero(observed_missing != predicted_missing))
    preserved_errors += int(np.count_nonzero(~observed_missing & ~predicted_missing & (observed_calls.sum(axis=2) != aligned.sum(axis=2))))
    preserved_errors += int(np.count_nonzero(predicted_missing & np.any(aligned >= 0, axis=2)))
    variant_order_ok = [locus[:2] for locus in predicted.loci] == [locus[:2] for locus in truth.loci] and predicted.variant_ids == truth.variant_ids
    target_retained = sum(locus[:2] == (target["chromosome"], target["position_bp"]) for locus in predicted.loci) == 1
    reference_ok = True
    manual_review = False
    reliable_carriers = unreliable_carriers = None
    map_errors = 0
    if stage == "12A":
        manifest = json.loads(paths["autosomal_phasing_manifest"].read_text())
        reference_ok = manifest.get("sample_selection", {}).get("membership_sha256") == reference_selection["membership_sha256"]
        for chromosome in range(1, 23):
            oracle_map = {int(row["POSITION_BP"]): float(row["POSITION_CM"]) for row in raw_table(preparation / f"reference_assets/chr{chromosome}.study_map.tsv")}
            rows_map = [line.split() for line in paths[f"autosomal_phased_chr{chromosome}_ibd_map"].read_text().splitlines()]
            actual_keys = [(int(row[0].removeprefix("chr")), int(row[3])) for row in rows_map]
            predicted_keys = [locus[:2] for locus in panels[chromosome-1].loci]
            map_errors += int(actual_keys != predicted_keys)
            map_errors += sum(int(row[3]) not in oracle_map or abs(float(row[2])-oracle_map.get(int(row[3]), 0)) > config["map_absolute_tolerance_cm"] for row in rows_map)
    else:
        manifest = json.loads(paths["reference_window_manifest"].read_text())
        reference_ok = manifest.get("sample_selection", {}).get("membership_sha256") == reference_selection["membership_sha256"]
        reference_header = subprocess.run(["bcftools", "query", "--list-samples", str(paths["reference_window_vcf"])], capture_output=True, text=True, timeout=config["command_timeout_seconds"], check=False)
        reference_ok &= reference_header.returncode == 0 and reference_header.stdout.splitlines() == reference_selection["sample_ids"]
        summary = json.loads(paths["phase_target_region_summary"].read_text())
        manual_review = summary["manual_validation_required"]
        reliable_carriers, unreliable_carriers = summary["reliable_carrier_count"], summary["unreliable_carrier_count"]
    performance.update(sample_order_correct=sample_order_ok, variant_order_correct=variant_order_ok, reference_membership_correct=reference_ok, target_retained=target_retained,
                       observed_genotype_preservation_errors=preserved_errors,
                       map_coordinate_error_count=map_errors, published_manual_review_required=manual_review,
                       reliable_carrier_count=reliable_carriers, unreliable_carrier_count=unreliable_carriers,
                       status="PASS" if sample_order_ok and variant_order_ok and reference_ok and preserved_errors == 0 and map_errors == 0 and performance["genotype_discordances"] == 0 and performance["unphased_heterozygotes"] == 0 and (target_retained or stage == "12A") else "FAIL",
                       chromosome_counts=[dict(chromosome=panel.loci[0][0], observed=len(panel.loci)) for panel in panels])
    return performance


def evaluate_bench(*, preparation_root: Path, run_root: Path, regional_roots: list[Path], configuration_path: Path, workspace_root: Path, run_id: str) -> dict[str, Any]:
    """Clôture quatre résultats indépendants et leurs agrégations de graphiques."""
    destination = new_workspace_run(workspace_root, run_id).absolute()
    workspace = workspace_root.absolute()
    preparation, original = preparation_root.absolute(), run_root.absolute()
    config = load_configuration(configuration_path)
    manifest = verify_preparation(preparation)
    if sha256_file(preparation / "preparation_manifest.json") != config["preparation_manifest_sha256"]:
        raise ValueError("d23_evaluation_preparation_changed")
    run = checked_closure(original)
    if [row["profile"] for row in run["profiles"]] != list(CASES) or run["configuration_sha256"] != sha256_file(configuration_path):
        raise ValueError("d23_evaluation_run_configuration_mismatch")
    regional = {}
    for root in regional_roots:
        root = root.absolute()
        receipt = checked_closure(root)
        if receipt["profile"] in regional or receipt["configuration_sha256"] != run["configuration_sha256"]:
            raise ValueError("d23_regional_duplicate_or_config_changed")
        regional[receipt["profile"]] = root, receipt
    if set(regional) != set(CASES):
        raise ValueError("d23_regional_profile_missing")
    with gzip.open(preparation / "private_truth/truth.phased.vcf.gz", "rt") as handle:
        truth = parse_vcf(handle)
    with gzip.open(preparation / "blinded_inputs/observed.vcf.gz", "rt") as handle:
        base_observed = read_panel(handle, phased=False)
    target = yaml.safe_load((preparation / "blinded_inputs/target_metadata.yaml").read_text())
    registry = raw_table(preparation / "blinded_inputs/metadata/samples.tsv")
    selection = json.loads((preparation / "blinded_inputs/reference_selection.json").read_text())
    results: dict[str, list[dict[str, Any]]] = {stage: [] for stage in STAGES}
    destination.mkdir(mode=0o700)
    for entry in run["profiles"]:
        case = entry["profile"]
        root = original / case
        regional_root, regional_receipt = regional[case]
        anchor = root / "stages/10A_prepare_autosomal_phasing_panel/stage_inputs.json"
        if sha256_file(anchor) != regional_receipt["source_profile_input_sha256"]:
            raise ValueError("d23_regional_parent_binding_invalid")
        source = root / "stages" / PRODUCER
        observation = observed_profile(base_observed, case, config, target)
        matrix = observation_matrix(observation, registry)
        converted = read_plink(source / "genomewide_pre_qc")
        conversion = compare_genotypes(matrix, converted)
        if conversion["status"] != "PASS":
            raise ValueError("d23_fixture_conversion_not_concordant")
        steps = {item["stage_id"]: (root, item) for item in entry["stages"] if item["stage_id"] in {"10A", "12A"}}
        steps.update({item["stage_id"]: (regional_root, item) for item in regional_receipt["stages"]})
        for stage in STAGES:
            if stage not in steps:
                results[stage].append(dict(profile=case, status="NOT_EVALUABLE", reason="UPSTREAM_COMPONENT_FAILED"))
                continue
            stage_root, step = steps[stage]
            paths = component_paths(stage_root, step, workspace)
            if step["returncode"]:
                results[stage].append(dict(profile=case, status="FAIL", returncode=step["returncode"], metrics_status="NOT_EVALUABLE"))
                continue
            if stage == "10A":
                result = evaluate_autosomal_input(source / "genomewide_pre_qc", paths["autosomal_phasing_panel_bed"].with_suffix(""))
            elif stage == "11":
                original_map = next((preparation / "reference_assets/genetic_maps").rglob(f"chr{target['chromosome']}.grch38.tsv"))
                result = evaluate_region(source / "target_chromosome_all_qc", paths["target_region_bed"].with_suffix(""), paths["target_genetic_map"], original_map, target, config)
            else:
                source_stage = "11" if stage == "12" else "10A"
                upstream_root, upstream_step = steps[source_stage]
                upstream_paths = component_paths(upstream_root, upstream_step, workspace)
                bed = upstream_paths["target_region_bed" if stage == "12" else "autosomal_phasing_panel_bed"]
                candidates = read_plink(bed.with_suffix(""))
                keys = {row[:2] for row in candidates.variants}
                if stage == "12A":
                    keys.discard((target["chromosome"], target["position_bp"]))
                expected_phase = phase_subset(truth, keys)
                result = evaluate_phase_outputs(paths, stage, expected_phase, registry, config, target, [sample for _, sample in candidates.samples], selection, observation, preparation)
            results[stage].append(dict(profile=case, **result))
    documents = []
    for stage in STAGES:
        rows = results[stage]
        for row in rows:
            row["criteria"] = acceptance_results(stage, row)
        document = dict(schema_version="1.0.0", method_id=METHOD_ID, stage_id=stage, evaluations=rows,
                        status="FAIL" if any(row["status"] == "FAIL" for row in rows) else "NOT_EVALUABLE" if any(row["status"] == "NOT_EVALUABLE" for row in rows) else "PASS",
                        performance_qualification="REVIEW_REQUIRED" if stage in {"12", "12A"} else "NOT_APPLICABLE", calibrated=False, founder_effect_proven=False)
        validate_evaluation(document)
        write_json(destination / (stage + ".evaluation.json"), document)
        documents.append(dict(stage_id=stage, status=document["status"], performance_qualification=document["performance_qualification"]))
    receipt = dict(schema_version="1.0.0", method_id=METHOD_ID, closure_state="CLOSED", created_at=datetime.now(timezone.utc).isoformat(),
                   preparation_manifest_sha256=config["preparation_manifest_sha256"], source_run_sha256=sha256_file(original / "run_index.json"),
                   regional_source_sha256=[sha256_file(root / "run_index.json") for root, _ in regional.values()], configuration_sha256=sha256_file(configuration_path),
                   implementation_sha256=sha256_file(Path(__file__)), oracle_implementation_sha256=sha256_file(Path(__file__).with_name("phasing_oracle.py")),
                   config=config, stages=documents, independent_units=manifest["counts"]["independent_simulation_units"],
                   acceptance_criteria_sha256=sha256_file(Path(__file__).resolve().parents[3] / "calibration/protocoles/CRITERES_ACCEPTATION.yaml"),
                   prespecification_sha256=sha256_file(Path(__file__).resolve().parents[3] / "calibration/protocoles/D23_PHASING_PRESPECIFICATION.md"),
                   limitations=["SYNTHETIC_TECHNICAL_TARGET_NOT_DOCK6", "SUBSAMPLED_CHIP_NOT_FULL_DENSITY", "PERFORMANCE_THRESHOLDS_NOT_SEALED", "D21_SCIENTIFIC_FAIL_07_UNCHANGED", "FIXTURE_CONTRACTS_NOT_UPSTREAM_QUALIFICATION"],
                   artifacts=file_inventory(destination, ("run_index.json",)), calibrated=False, founder_effect_proven=False)
    validate_evaluation(receipt, manifest=True)
    write_json(destination / "run_index.json", receipt)
    return receipt


def main() -> int:
    """CLI d'évaluation distincte du processus de phasage."""
    parser = argparse.ArgumentParser(description="Évaluation D23 indépendante")
    for name in ("preparation-root", "run-root", "configuration-path", "workspace-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--regional-roots", type=Path, nargs=3, required=True)
    parser.add_argument("--run-id", required=True)
    receipt = evaluate_bench(**vars(parser.parse_args()))
    print(receipt["stages"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
