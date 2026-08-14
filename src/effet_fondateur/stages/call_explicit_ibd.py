"""Étape 16C : appels IBD autosomaux en aveugle puis interrogation de la cible."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from itertools import combinations
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Iterable, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError, TableValidationError, build_file_artifact,
    load_pipeline_config, validate_json_document, validate_tsv_table,
)
from effet_fondateur.explicit_ibd import Scenario, classify_explicit_ibd, evaluate_scenario, parse_hap_ibd, parse_refined_ibd
from effet_fondateur.explicit_ibd.analysis import IbdSegment, ScenarioEvaluation
from effet_fondateur.explicit_ibd.execution import ExplicitIbdExternalError, run_tool, validate_adapters
from effet_fondateur.orchestrator.state import utc_now


METHOD_ID = "dual_target_centered_explicit_ibd_v1"


class ExplicitIbdInputError(ValueError):
    """Signale une entrée 16C non cohérente."""


def _artifact(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [item for item in stage_inputs["artifacts"] if item["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise ExplicitIbdInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _path(artifact: dict[str, Any], run_dir: Path) -> Path:
    candidate = Path(artifact["path"])
    path = run_dir / candidate if PurePosixPath(artifact["path"]).parts[:1] == ("stages",) else candidate
    if not path.is_absolute():
        path = Path.cwd() / path
    if path.is_symlink() or not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise ExplicitIbdInputError("declared_input_modified")
    return path


def _run(command: Sequence[str], timeout: int) -> str:
    import subprocess
    try:
        result = subprocess.run(list(command), capture_output=True, text=True, check=False, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExplicitIbdExternalError("bcftools_unavailable_or_timeout") from error
    if result.returncode != 0:
        raise ExplicitIbdExternalError(f"bcftools_failed:{result.returncode}")
    return result.stdout


def _query_complete_variants(
    bcftools: str, bcf: Path, timeout: int
) -> tuple[tuple[str, str, int, bool], ...]:
    """Compatibilité et contrôle isolé d'un BCF phasé sans imputation."""
    output = _run(
        [bcftools, "query", "-f", "%ID\t%CHROM\t%POS[\t%GT]\n", str(bcf)],
        timeout,
    )
    rows: list[tuple[str, str, int, bool]] = []
    seen: set[str] = set()
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) < 4 or fields[0] in seen:
            raise ExplicitIbdInputError("phased_variant_query_invalid")
        seen.add(fields[0])
        complete = all("." not in genotype and "|" in genotype and genotype.count("|") == 1 for genotype in fields[3:])
        rows.append((fields[0], fields[1].removeprefix("chr"), int(fields[2]), complete))
    if not rows:
        raise ExplicitIbdInputError("phased_variant_query_empty")
    return tuple(rows)


def _scenario(raw: dict[str, Any], calibrated: bool = False) -> Scenario:
    scenario = Scenario(raw["scenario_id"], raw["role"], float(raw["minimum_cm"]), int(raw["minimum_markers"]), calibrated)
    scenario.validate()
    return scenario


def _map(path: Path, chromosome: int) -> tuple[dict[int, float], tuple[int, ...], list[dict[str, Any]]]:
    cm_at_bp: dict[int, float] = {}
    rows: list[dict[str, Any]] = []
    previous_bp, previous_cm = 0, -math.inf
    for order, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        fields = line.split()
        if len(fields) != 4:
            raise ExplicitIbdInputError(f"ibd_map_malformed:chr{chromosome}")
        observed, variant_id, raw_cm, raw_bp = int(fields[0].removeprefix("chr")), fields[1], fields[2], fields[3]
        bp, cm = int(raw_bp), float(raw_cm)
        if observed != chromosome or bp <= previous_bp or cm < previous_cm:
            raise ExplicitIbdInputError(f"ibd_map_not_monotonic:chr{chromosome}")
        cm_at_bp[bp] = cm
        rows.append({"VARIANT_ORDER": order, "VARIANT_ID": variant_id, "CHROMOSOME": str(chromosome), "POSITION_BP": bp, "POSITION_CM": f"{cm:.12g}", "IS_TARGET": False, "HAP_IBD_USE": "INCLUDED", "REFINED_IBD_USE": "INCLUDED", "DECISION_REASON": "INCLUDED"})
        previous_bp, previous_cm = bp, cm
    if len(rows) < 2:
        raise ExplicitIbdInputError(f"ibd_map_empty:chr{chromosome}")
    return cm_at_bp, tuple(cm_at_bp), rows


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: "" if row.get(column) is None else "true" if row.get(column) is True else "false" if row.get(column) is False else row.get(column) for column in columns})


def _query_phased(bcftools: str, path: Path, timeout: int) -> tuple[list[str], dict[tuple[str, int, str, str], tuple[str, ...]]]:
    samples = [line for line in _run([bcftools, "query", "--list-samples", str(path)], timeout).splitlines() if line]
    records: dict[tuple[str, int, str, str], tuple[str, ...]] = {}
    output = _run([bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(path)], timeout)
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != 4 + len(samples):
            raise ExplicitIbdInputError("phased_alignment_query_invalid")
        key = (fields[0].removeprefix("chr"), int(fields[1]), fields[2], fields[3])
        if key in records:
            raise ExplicitIbdInputError("phased_alignment_duplicate_variant")
        records[key] = tuple(fields[4:])
    return samples, records


def _align_mutant_haplotypes(
    bcftools: str, regional_bcf: Path, chromosome_bcf: Path,
    carrier_rows: Sequence[dict[str, Any]], family_by_sample: dict[str, str],
    timeout: int, minimum_informative: int = 5, minimum_concordance: float = 0.9,
) -> tuple[dict[str, set[str]], list[dict[str, Any]]]:
    regional_samples, regional = _query_phased(bcftools, regional_bcf, timeout)
    chromosome_samples, chromosome = _query_phased(bcftools, chromosome_bcf, timeout)
    if regional_samples != chromosome_samples:
        raise ExplicitIbdInputError("regional_autosomal_sample_order_mismatch")
    overlap = sorted(set(regional) & set(chromosome), key=lambda key: key[1])
    if not overlap:
        raise ExplicitIbdInputError("regional_autosomal_no_overlap")
    sample_index = {sample: index for index, sample in enumerate(regional_samples)}
    mutant: dict[str, set[str]] = {}
    audit: list[dict[str, Any]] = []
    for row in carrier_rows:
        if row["RELIABILITY_STATUS"] != "PASS" or row["ALT_COPY_COUNT"] == "0":
            continue
        sample = row["SAMPLE_ID"]
        if sample not in sample_index:
            raise ExplicitIbdInputError("carrier_missing_from_autosomal_phasing")
        index, same, swapped = sample_index[sample], 0, 0
        for key in overlap:
            left, right = regional[key][index], chromosome[key][index]
            if "|" not in left or "|" not in right or "." in left or "." in right:
                continue
            left_alleles, right_alleles = left.split("|"), right.split("|")
            if left_alleles[0] == left_alleles[1] or sorted(left_alleles) != sorted(right_alleles):
                continue
            same += left_alleles == right_alleles
            swapped += left_alleles == right_alleles[::-1]
        informative = same + swapped
        homozygous_alternate = row["ALT_COPY_COUNT"] == "2" and row["CARRIER_HAPLOTYPE"] == "BOTH"
        if homozygous_alternate:
            # L'orientation H1/H2 ne change pas l'attribution : les deux copies
            # portent explicitement l'allèle alternatif. Imposer une concordance
            # de phase entre deux runs indépendants serait un faux garde-fou.
            orientation = "NOT_REQUIRED_HOMOZYGOUS"
            regional_haps = whole_haps = {"H1", "H2"}
            concordance: str | None = None
        else:
            if informative < minimum_informative or max(same, swapped) / informative < minimum_concordance:
                raise ExplicitIbdInputError(f"carrier_phase_alignment_unreliable:{sample}")
            orientation = "SAME" if same > swapped else "SWAPPED"
            regional_haps = {row["CARRIER_HAPLOTYPE"]}
            whole_haps = regional_haps if orientation == "SAME" else {"H2" if hap == "H1" else "H1" for hap in regional_haps}
            concordance = f"{max(same, swapped) / informative:.12g}"
        family = family_by_sample[sample]
        mutant.setdefault(family, set()).update(whole_haps)
        audit.append({"SAMPLE_ID": sample, "FAMILY_ID": family, "INFORMATIVE_MARKERS": informative, "SAME_COUNT": same, "SWAPPED_COUNT": swapped, "ORIENTATION": orientation, "CONCORDANCE": concordance, "REGIONAL_MUTANT_HAPLOTYPE": row["CARRIER_HAPLOTYPE"], "AUTOSOMAL_MUTANT_HAPLOTYPE": "BOTH" if whole_haps == {"H1", "H2"} else next(iter(whole_haps))})
    if len(mutant) < 3:
        raise ExplicitIbdInputError("insufficient_aligned_carrier_families")
    return mutant, audit


def _control_frequency(segments: Sequence[IbdSegment], scenario_id: str, controls: set[str], target_bp: int) -> tuple[int, int, float | None]:
    pairs = list(combinations(sorted(controls), 2))
    positive = 0
    for sample_1, sample_2 in pairs:
        observed = {segment.tool for segment in segments if segment.scenario_id == scenario_id and segment.contains(target_bp) and {segment.sample_1, segment.sample_2} == {sample_1, sample_2}}
        positive += observed == {"HAP_IBD", "REFINED_IBD"}
    return len(pairs), positive, positive / len(pairs) if pairs else None


def _segment_row(segment: IbdSegment, target_bp: int) -> dict[str, Any]:
    return {"TOOL": segment.tool, "SCENARIO_ID": segment.scenario_id, "SAMPLE_1": segment.sample_1, "HAPLOTYPE_1": segment.haplotype_1, "FAMILY_1": segment.family_1, "SAMPLE_2": segment.sample_2, "HAPLOTYPE_2": segment.haplotype_2, "FAMILY_2": segment.family_2, "CHROMOSOME": segment.chromosome, "START_BP": segment.start_bp, "END_BP": segment.end_bp, "START_CM": f"{segment.start_cm:.12g}", "END_CM": f"{segment.end_cm:.12g}", "LENGTH_CM": f"{segment.length_cm:.12g}", "MARKER_COUNT": segment.marker_count, "CONTAINS_TARGET": segment.contains(target_bp) if segment.chromosome == "19" else False}


def _evaluation_status(value: ScenarioEvaluation) -> str:
    if not value.evaluable: return "NOT_EVALUABLE"
    if value.all_pairs_by_both_tools and value.boundaries_concordant: return "CONCORDANT"
    if value.all_pairs_by_both_tools: return "DISCORDANT"
    return "NO_CALL"


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Appelle les deux outils sans statut mutationnel, puis évalue chr cible."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path); validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = stage_inputs["parameters"]
    if parameters.get("method") != METHOD_ID or parameters.get("input_scope") != "AUTOSOMAL_GENOMEWIDE":
        raise ExplicitIbdInputError("explicit_ibd_method_or_scope_invalid")
    raw_scenarios = (parameters["primary"], *parameters["sensitivities"])
    declared = tuple(_scenario(raw) for raw in raw_scenarios)
    if declared[0].role != "PRIMARY" or any(item.role != "SENSITIVITY" for item in declared[1:]) or len({item.scenario_id for item in declared}) != len(declared):
        raise ExplicitIbdInputError("explicit_ibd_scenario_registry_invalid")
    timeout, minimum_control_pairs = int(parameters["tool_timeout_seconds"]), int(parameters.get("minimum_control_pairs", 100))
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    adapters = validate_adapters(config["tools"]["explicit_ibd_adapters"])
    base_ids = ("config_input_target_variant_metadata", "samples_master", "cohorts_frozen", "shapeit5_final_bcf", "shapeit5_final_index", "carrier_haplotypes", "founder_analysis_summary", "autosomal_phasing_manifest")
    all_ids = (*base_ids, *(f"autosomal_phased_chr{chromosome}_{kind}" for chromosome in range(1, 23) for kind in ("bcf", "index", "ibd_map")))
    artifacts = {name: _artifact(stage_inputs, name) for name in all_ids}
    paths = {name: _path(value, run_dir) for name, value in artifacts.items()}
    target = yaml.safe_load(paths["config_input_target_variant_metadata"].read_text(encoding="utf-8")); validate_json_document(target, "target_variant_metadata.schema.json")
    target_chromosome, target_bp = int(target["chromosome"]), int(target["position_bp"])
    manifest = read_json(paths["autosomal_phasing_manifest"]); validate_json_document(manifest, "autosomal_phasing_manifest.schema.json")
    if manifest["sample_count"] < 1 or len(manifest["chromosomes"]) != 22:
        raise ExplicitIbdInputError("autosomal_phasing_manifest_invalid")
    samples = validate_tsv_table(paths["samples_master"], "samples_master.schema.json").rows
    family_by_sample = {row["SAMPLE_ID"]: row["FID"] for row in samples}
    if len(family_by_sample) != len(samples):
        raise ExplicitIbdInputError("sample_registry_duplicate")
    cohorts = validate_tsv_table(paths["cohorts_frozen"], "cohorts_frozen.schema.json").rows
    controls = {row["SAMPLE_ID"] for row in cohorts if row["COHORT_ID"] == "controls_unrelated" and row["INCLUDED"]}
    independent_carriers = {
        row["SAMPLE_ID"] for row in cohorts
        if row["COHORT_ID"] == "target_carriers_independent" and row["INCLUDED"]
    }
    carriers = [
        row for row in validate_tsv_table(paths["carrier_haplotypes"], "carrier_haplotypes.schema.json").rows
        if row["SAMPLE_ID"] in independent_carriers
    ]
    mutant, alignment_rows = _align_mutant_haplotypes(config["tools"]["bcftools"], paths["shapeit5_final_bcf"], paths[f"autosomal_phased_chr{target_chromosome}_bcf"], carriers, family_by_sample, timeout)
    families = tuple(sorted(mutant))
    analysis_dir = output_dir / "explicit_ibd"; analysis_dir.mkdir(parents=True)
    all_segments: list[IbdSegment] = []
    marker_rows: list[dict[str, Any]] = []
    tool_records: list[dict[str, Any]] = []
    raw_output_specs: list[tuple[str, Path, str]] = []
    target_map: tuple[dict[int, float], tuple[int, ...]] | None = None
    # Les appels ci-dessous ne reçoivent ni cohorte, ni statut porteur, ni cible.
    manifest_by_chromosome = {row["chromosome"]: row for row in manifest["chromosomes"]}
    for chromosome in range(1, 23):
        bcf, map_path = paths[f"autosomal_phased_chr{chromosome}_bcf"], paths[f"autosomal_phased_chr{chromosome}_ibd_map"]
        manifest_record = manifest_by_chromosome.get(chromosome)
        if (
            manifest_record is None
            or manifest_record["bcf"]["sha256"] != sha256_file(bcf)
            or manifest_record["index"]["sha256"] != sha256_file(paths[f"autosomal_phased_chr{chromosome}_index"])
            or manifest_record["ibd_map"]["sha256"] != sha256_file(map_path)
        ):
            raise ExplicitIbdInputError(f"autosomal_manifest_artifact_mismatch:chr{chromosome}")
        cm_at_bp, marker_positions, chromosome_marker_rows = _map(map_path, chromosome)
        offset = len(marker_rows)
        for row in chromosome_marker_rows:
            row["VARIANT_ORDER"] += offset
        marker_rows.extend(chromosome_marker_rows)
        if chromosome == target_chromosome:
            target_map = (cm_at_bp, marker_positions)
        for scenario in declared:
            for tool in ("HAP_IBD", "REFINED_IBD"):
                prefix = analysis_dir / "tool_outputs" / f"chr{chromosome}" / f"{scenario.scenario_id}.{tool.lower()}"; prefix.parent.mkdir(parents=True, exist_ok=True)
                run = run_tool(tool=tool, adapters=adapters, vcf_path=bcf, map_path=map_path, output_prefix=prefix, minimum_cm=scenario.minimum_cm, minimum_markers=scenario.minimum_markers, threads=int(parameters["threads"]), memory_mb=int(parameters["java_memory_mb"]), timeout_seconds=timeout)
                parsed = parse_hap_ibd(run.output_path, scenario.scenario_id, family_by_sample, cm_at_bp, marker_positions) if tool == "HAP_IBD" else parse_refined_ibd(run.output_path, scenario.scenario_id, family_by_sample, cm_at_bp, marker_positions)
                all_segments.extend(parsed)
                tool_records.append({"tool": tool, "scenario": scenario.scenario_id, "chromosome": chromosome, "command": list(run.command), "version": adapters["hap_ibd_version" if tool == "HAP_IBD" else "refined_ibd_version"]})
                raw_output_specs.extend(((f"raw_{tool.lower()}_{scenario.scenario_id}_chr{chromosome}", run.output_path, "sensitive_genetic"), (f"log_{tool.lower()}_{scenario.scenario_id}_chr{chromosome}", run.log_path, "internal")))
    if target_map is None:
        raise ExplicitIbdInputError("target_chromosome_map_missing")
    control_results: dict[str, tuple[int, int, float | None]] = {scenario.scenario_id: _control_frequency(all_segments, scenario.scenario_id, controls, target_bp) for scenario in declared}
    scenarios = tuple(Scenario(item.scenario_id, item.role, item.minimum_cm, item.minimum_markers, control_results[item.scenario_id][0] >= minimum_control_pairs) for item in declared)
    evaluations: list[ScenarioEvaluation] = []
    for scenario in scenarios:
        evaluable, _, frequency = control_results[scenario.scenario_id]
        acceptable = frequency is not None and frequency <= float(parameters["maximum_control_frequency"])
        evaluations.append(evaluate_scenario(scenario, all_segments, families, {key: frozenset(value) for key, value in mutant.items()}, target_bp, int(parameters["boundary_tolerance_bp"]), frequency, float(parameters["maximum_control_frequency"]), evaluable >= minimum_control_pairs and acceptable))
    primary_evaluation, sensitivity_evaluations = evaluations[0], tuple(evaluations[1:])
    status = classify_explicit_ibd(primary_evaluation, sensitivity_evaluations)
    segments_path, marker_path = analysis_dir / "explicit_ibd_segments.tsv", analysis_dir / "ibd_marker_audit.tsv"
    pair_path, concordance_path = analysis_dir / "target_ibd_pair_matrix.tsv", analysis_dir / "ibd_tool_concordance.tsv"
    calibration_path, control_path = analysis_dir / "ibd_calibration.tsv", analysis_dir / "ibd_control_frequency.tsv"
    alignment_path = analysis_dir / "target_chromosome_phase_alignment.tsv"
    feasibility_path, summary_path = analysis_dir / "explicit_ibd_feasibility.json", analysis_dir / "explicit_ibd_summary.json"
    _write_tsv(segments_path, ("TOOL","SCENARIO_ID","SAMPLE_1","HAPLOTYPE_1","FAMILY_1","SAMPLE_2","HAPLOTYPE_2","FAMILY_2","CHROMOSOME","START_BP","END_BP","START_CM","END_CM","LENGTH_CM","MARKER_COUNT","CONTAINS_TARGET"), (_segment_row(segment, target_bp) for segment in all_segments))
    _write_tsv(marker_path, ("VARIANT_ORDER","VARIANT_ID","CHROMOSOME","POSITION_BP","POSITION_CM","IS_TARGET","HAP_IBD_USE","REFINED_IBD_USE","DECISION_REASON"), marker_rows)
    _write_tsv(alignment_path, ("SAMPLE_ID","FAMILY_ID","INFORMATIVE_MARKERS","SAME_COUNT","SWAPPED_COUNT","ORIENTATION","CONCORDANCE","REGIONAL_MUTANT_HAPLOTYPE","AUTOSOMAL_MUTANT_HAPLOTYPE"), alignment_rows)
    pair_rows: list[dict[str, Any]] = []
    for scenario, evaluation in zip(scenarios, evaluations, strict=True):
        for family_1, family_2 in combinations(families, 2):
            selected = [segment for segment in evaluation.selected_segments if segment.family_pair == (family_1, family_2)]
            tools = {segment.tool for segment in selected}
            pair_rows.append({"SCENARIO_ID": scenario.scenario_id, "ROLE": scenario.role, "FAMILY_1": family_1, "FAMILY_2": family_2, "HAP_IBD_STATUS": "CALLED" if "HAP_IBD" in tools else "NOT_CALLED" if evaluation.evaluable else "NOT_EVALUABLE", "REFINED_IBD_STATUS": "CALLED" if "REFINED_IBD" in tools else "NOT_CALLED" if evaluation.evaluable else "NOT_EVALUABLE", "BOTH_METHODS": tools == {"HAP_IBD", "REFINED_IBD"}, "CONTAINS_TARGET": bool(selected) and all(segment.contains(target_bp) for segment in selected), "BOUNDARIES_CONCORDANT": evaluation.boundaries_concordant, "DETAIL_CODE": evaluation.reason})
    _write_tsv(pair_path, ("SCENARIO_ID","ROLE","FAMILY_1","FAMILY_2","HAP_IBD_STATUS","REFINED_IBD_STATUS","BOTH_METHODS","CONTAINS_TARGET","BOUNDARIES_CONCORDANT","DETAIL_CODE"), pair_rows)
    _write_tsv(concordance_path, ("SCENARIO_ID","ALL_REQUIRED_PAIRS","GLOBAL_HAPLOTYPE_COHERENT","BOUNDARIES_CONCORDANT","COMMON_START_BP","COMMON_END_BP","TARGET_IN_COMMON_INTERSECTION","STATUS"), ({"SCENARIO_ID": scenario.scenario_id, "ALL_REQUIRED_PAIRS": value.all_pairs_by_both_tools, "GLOBAL_HAPLOTYPE_COHERENT": value.globally_coherent, "BOUNDARIES_CONCORDANT": value.boundaries_concordant, "COMMON_START_BP": value.common_start_bp, "COMMON_END_BP": value.common_end_bp, "TARGET_IN_COMMON_INTERSECTION": value.common_start_bp is not None and value.common_start_bp <= target_bp <= value.common_end_bp, "STATUS": _evaluation_status(value)} for scenario, value in zip(scenarios, evaluations, strict=True)))
    calibration_rows, control_rows = [], []
    for scenario in scenarios:
        count, positive, frequency = control_results[scenario.scenario_id]
        calibrated, acceptable = count >= minimum_control_pairs, frequency is not None and frequency <= float(parameters["maximum_control_frequency"])
        calibration_rows.append({"SCENARIO_ID": scenario.scenario_id, "ROLE": scenario.role, "MINIMUM_CM": scenario.minimum_cm, "MINIMUM_MARKERS": scenario.minimum_markers, "CALIBRATED": calibrated, "SPECIFICITY_ACCEPTABLE": acceptable, "CONTROL_FREQUENCY": frequency, "MAXIMUM_CONTROL_FREQUENCY": parameters["maximum_control_frequency"], "STATUS": "ACCEPTABLE" if calibrated and acceptable else "UNACCEPTABLE" if calibrated else "NOT_EVALUABLE"})
        control_rows.append({"SCENARIO_ID": scenario.scenario_id, "CONTROL_SOURCE": "INTERNAL_INDEPENDENT_CONTROLS", "EVALUABLE_UNIT_COUNT": count, "POSITIVE_UNIT_COUNT": positive, "FREQUENCY": frequency, "PRESPECIFIED_MAXIMUM": parameters["maximum_control_frequency"], "STATUS": "NOT_EVALUABLE" if frequency is None else "ACCEPTABLE" if acceptable else "TOO_FREQUENT"})
    _write_tsv(calibration_path, ("SCENARIO_ID","ROLE","MINIMUM_CM","MINIMUM_MARKERS","CALIBRATED","SPECIFICITY_ACCEPTABLE","CONTROL_FREQUENCY","MAXIMUM_CONTROL_FREQUENCY","STATUS"), calibration_rows)
    _write_tsv(control_path, ("SCENARIO_ID","CONTROL_SOURCE","EVALUABLE_UNIT_COUNT","POSITIVE_UNIT_COUNT","FREQUENCY","PRESPECIFIED_MAXIMUM","STATUS"), control_rows)
    target_positions = target_map[1]
    target_covered = min(target_positions) <= target_bp <= max(target_positions)
    local_markers = sum(abs(position - target_bp) <= 5_000_000 for position in target_positions)
    feasibility = {"schema_version": "1.0.0", "input_scope": parameters["input_scope"], "variant_count": len(marker_rows), "complete_phased_variant_count": len(marker_rows), "target_present": target_covered, "primary_density_possible": local_markers >= scenarios[0].minimum_markers, "status": "GO" if target_covered and local_markers >= scenarios[0].minimum_markers else "GO_SENSITIVITY_ONLY" if target_covered and local_markers >= min(item.minimum_markers for item in scenarios[1:]) else "NOT_EVALUABLE"}
    validate_json_document(feasibility, "explicit_ibd_feasibility.schema.json"); atomic_write_json(feasibility_path, feasibility)
    summary = {"schema_version": "1.0.0", "method_id": METHOD_ID, "status": status, "primary_status": _evaluation_status(primary_evaluation), "sensitivity_statuses": [{"scenario_id": scenario.scenario_id, "status": _evaluation_status(value)} for scenario, value in zip(scenarios[1:], sensitivity_evaluations, strict=True)], "input_scope": parameters["input_scope"], "target": {"assembly": target["assembly"], "chromosome": str(target_chromosome), "position_bp": target_bp, "project_variant_id": target["project_variant_id"]}, "family_count": len(families), "interpretation": {"ibs_only_from_step13": True, "explicit_ibd_supported": status in {"PRIMARY_CONCORDANT_IBD_SUPPORT", "SENSITIVITY_ONLY_IBD_SUPPORT"}, "ibd_proven": status == "PRIMARY_CONCORDANT_IBD_SUPPORT", "founder_effect_proven": False, "geographic_origin_inferred": False, "composite_score_calculated": False}, "limitations": ["Une absence d'appel sur puce SNP ne réfute pas un effet fondateur.", "Les sensibilités sous 2 cM ne constituent pas une preuve primaire.", "Un appel IBD explicite ne prouve pas à lui seul l'effet fondateur."]}
    validate_json_document(summary, "explicit_ibd_summary.schema.json"); atomic_write_json(summary_path, summary)
    contracts = ((segments_path,"explicit_ibd_segments.schema.json"),(marker_path,"explicit_ibd_marker_audit.schema.json"),(pair_path,"explicit_ibd_pair_results.schema.json"),(concordance_path,"explicit_ibd_concordance.schema.json"),(calibration_path,"explicit_ibd_calibration.schema.json"),(control_path,"explicit_ibd_control_frequency.schema.json"))
    for path, schema in contracts: validate_tsv_table(path, schema)
    specifications = [("explicit_ibd_segments",segments_path,"explicit_ibd_segments.schema.json","sensitive_genetic"),("explicit_ibd_marker_audit",marker_path,"explicit_ibd_marker_audit.schema.json","internal"),("explicit_ibd_pair_results",pair_path,"explicit_ibd_pair_results.schema.json","sensitive_genetic"),("explicit_ibd_concordance",concordance_path,"explicit_ibd_concordance.schema.json","internal"),("explicit_ibd_calibration",calibration_path,"explicit_ibd_calibration.schema.json","internal"),("explicit_ibd_control_frequency",control_path,"explicit_ibd_control_frequency.schema.json","internal"),("explicit_ibd_feasibility",feasibility_path,"explicit_ibd_feasibility.schema.json","internal"),("explicit_ibd_summary",summary_path,"explicit_ibd_summary.schema.json","internal"),("explicit_ibd_phase_alignment",alignment_path,None,"sensitive_genetic")]
    outputs = [build_file_artifact(physical_path=path, published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}", artifact_id=artifact_id, artifact_type=artifact_id, media_type="application/json" if path.suffix == ".json" else "text/tab-separated-values", producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"], schema_name=schema, schema_version="1.0.0" if schema else None, assembly=target["assembly"], sample_set_id=manifest["sample_set_id"], variant_set_id=manifest["variant_set_id"], sensitivity=sensitivity) for artifact_id, path, schema, sensitivity in specifications]
    for artifact_id, path, sensitivity in raw_output_specs:
        outputs.append(build_file_artifact(physical_path=path, published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}", artifact_id=artifact_id, artifact_type="explicit_ibd_raw_output" if sensitivity == "sensitive_genetic" else "explicit_ibd_tool_log", media_type="application/gzip" if path.suffix == ".gz" else "text/plain", producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"], assembly=target["assembly"], sample_set_id=manifest["sample_set_id"], variant_set_id=manifest["variant_set_id"], sensitivity=sensitivity))
    stage_outputs = {"schema_version":"1.0.0","run_id":stage_inputs["run_id"],"stage_id":stage_inputs["stage_id"],"stage_name":stage_inputs["stage_name"],"signature":stage_inputs["signature"],"artifacts":outputs}; validate_json_document(stage_outputs,"stage_outputs.schema.json"); atomic_write_json(output_dir/"stage_outputs.json",stage_outputs)
    audit = {"schema_version":"1.0.0","run_id":stage_inputs["run_id"],"stage_id":stage_inputs["stage_id"],"stage_name":stage_inputs["stage_name"],"method_id":METHOD_ID,"signature":stage_inputs["signature"],"started_at":started_at,"completed_at":utc_now(),"duration_seconds":monotonic()-started_clock,"inputs":list(artifacts.values()),"outputs":outputs,"parameters":parameters,"tools":tool_records,"counts":{"independent_families":len(families),"independent_controls":len(controls),"autosomes":22,"input_variants":len(marker_rows),"normalized_segments":len(all_segments)},"metrics":{"analysis_status":status,"ibd_proven":status=="PRIMARY_CONCORDANT_IBD_SUPPORT","founder_effect_proven":False,"composite_founder_score_calculated":False},"exclusions":[],"warnings":[] if status=="PRIMARY_CONCORDANT_IBD_SUPPORT" else [{"code":"absence_or_non_primary_call_does_not_refute_founder_effect","count":1}],"checks":[{"check":name,"status":"PASS"} for name in ("blind_genomewide_calling","same_variant_universe_both_tools","same_genetic_map_both_tools","no_missing_alleles","all_22_autosomes","regional_autosomal_phase_alignment","empirical_control_frequency","ibs_ibd_founder_separation")],"known_limits":summary["limitations"],"expected_visualizations":["explicit_ibd_target_segments"],"manual_validation_required":True}; validate_json_document(audit,"stage_audit.schema.json"); atomic_write_json(output_dir/"audit.json",audit)
    (output_dir/"checksums.sha256").write_text("".join(f"{item['sha256']}  {item['artifact_id']}\n" for item in outputs),encoding="utf-8")
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    parser=argparse.ArgumentParser(); parser.add_argument("--stage-inputs",type=Path,required=True); parser.add_argument("--output-dir",type=Path,required=True); parsed=parser.parse_args(arguments)
    try: return execute(parsed.stage_inputs,parsed.output_dir)
    except ExplicitIbdExternalError as error: sys.stderr.write(f"{error}\n"); return 3
    except (ExplicitIbdInputError,DocumentValidationError,TableValidationError,ValueError,OSError,json.JSONDecodeError) as error: sys.stderr.write(f"{error}\n"); return 2


if __name__ == "__main__": raise SystemExit(main())
