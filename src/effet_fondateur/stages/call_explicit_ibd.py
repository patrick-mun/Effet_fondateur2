"""Étape 16C : appel IBD explicite Hap-IBD/Refined IBD centré cible."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
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
    """Signale une entrée ou provenance 16C invalide."""


def _artifact(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [item for item in stage_inputs["artifacts"] if item["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise ExplicitIbdInputError(f"artifact_missing_or_ambiguous:{artifact_id}")
    return matches[0]


def _path(item: dict[str, Any], run_dir: Path) -> Path:
    raw = Path(item["path"])
    path = run_dir / raw if PurePosixPath(item["path"]).parts[:1] == ("stages",) else (raw if raw.is_absolute() else Path.cwd() / raw)
    if path.is_symlink() or not path.is_file() or sha256_file(path) != item["sha256"]:
        raise ExplicitIbdInputError(f"artifact_integrity_failure:{item['artifact_id']}")
    return path


def _run(command: list[str], timeout: int) -> str:
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ExplicitIbdExternalError("bcftools_unavailable_or_timeout") from error
    if completed.returncode != 0:
        raise ExplicitIbdExternalError(f"bcftools_failed:{completed.returncode}")
    return completed.stdout


def _query_complete_variants(bcftools: str, bcf: Path, timeout: int) -> tuple[tuple[str, str, int, bool], ...]:
    """Inventorie les variants et exclut tout marqueur ayant un GT absent/non phasé."""
    output = _run([bcftools, "query", "-f", "%ID\t%CHROM\t%POS[\t%GT]\n", str(bcf)], timeout)
    rows: list[tuple[str, str, int, bool]] = []
    seen: set[str] = set()
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) < 4 or fields[0] in seen:
            raise ExplicitIbdInputError("phased_variant_query_invalid")
        seen.add(fields[0])
        complete = all("." not in gt and "|" in gt and gt.count("|") == 1 for gt in fields[3:])
        rows.append((fields[0], fields[1].removeprefix("chr"), int(fields[2]), complete))
    if not rows:
        raise ExplicitIbdInputError("phased_variant_query_empty")
    return tuple(rows)


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: "" if row.get(name) is None else ("true" if row.get(name) is True else "false" if row.get(name) is False else row.get(name)) for name in columns})


def _scenario(raw: dict[str, Any]) -> Scenario:
    result = Scenario(raw["scenario_id"], raw["role"], float(raw["minimum_cm"]), int(raw["minimum_markers"]), bool(raw["calibrated"]))
    result.validate()
    return result


def _segment_row(segment: IbdSegment, target_bp: int) -> dict[str, Any]:
    return {"TOOL": segment.tool, "SCENARIO_ID": segment.scenario_id, "SAMPLE_1": segment.sample_1, "HAPLOTYPE_1": segment.haplotype_1, "FAMILY_1": segment.family_1, "SAMPLE_2": segment.sample_2, "HAPLOTYPE_2": segment.haplotype_2, "FAMILY_2": segment.family_2, "CHROMOSOME": segment.chromosome, "START_BP": segment.start_bp, "END_BP": segment.end_bp, "START_CM": f"{segment.start_cm:.12g}", "END_CM": f"{segment.end_cm:.12g}", "LENGTH_CM": f"{segment.length_cm:.12g}", "MARKER_COUNT": segment.marker_count, "CONTAINS_TARGET": segment.contains(target_bp)}


def _evaluation_status(value: ScenarioEvaluation) -> str:
    if not value.evaluable:
        return "NOT_EVALUABLE"
    if value.all_pairs_by_both_tools and value.boundaries_concordant:
        return "CONCORDANT"
    if value.all_pairs_by_both_tools:
        return "DISCORDANT"
    return "NO_CALL"


def _prepare_tool_inputs(bcftools: str, bcf: Path, output_dir: Path, map_rows: Sequence[dict[str, Any]], query_rows: Sequence[tuple[str, str, int, bool]], timeout: int) -> tuple[Path, Path, list[dict[str, Any]]]:
    by_id = {row["VARIANT_ID"]: row for row in map_rows}
    audit: list[dict[str, Any]] = []
    included: list[dict[str, Any]] = []
    for order, (variant_id, chrom, position, complete) in enumerate(query_rows, start=1):
        map_row = by_id.get(variant_id)
        if map_row is None or int(map_row["POSITION_BP"]) != position or str(map_row["CHROMOSOME"]) != chrom:
            raise ExplicitIbdInputError("ibd_map_bcf_universe_mismatch")
        reason = "INCLUDED" if complete else "TRUE_MISSING_CALL"
        if complete:
            included.append(map_row)
        audit.append({"VARIANT_ORDER": order, "VARIANT_ID": variant_id, "CHROMOSOME": chrom, "POSITION_BP": position, "POSITION_CM": map_row["POSITION_CM"], "IS_TARGET": map_row["IS_TARGET_VARIANT"], "HAP_IBD_USE": "INCLUDED" if complete else "EXCLUDED", "REFINED_IBD_USE": "INCLUDED" if complete else "EXCLUDED", "DECISION_REASON": reason})
    if not included or not any(row["IS_TARGET_VARIANT"] for row in included):
        raise ExplicitIbdInputError("target_missing_or_incomplete_for_explicit_ibd")
    positions = output_dir / "common_variants.positions.tsv"
    positions.write_text("".join(f"{row['CHROMOSOME']}\t{row['POSITION_BP']}\n" for row in included), encoding="utf-8")
    vcf = output_dir / "explicit_ibd_common.vcf.gz"
    _run([bcftools, "view", "-T", str(positions), "-Oz", "-o", str(vcf), str(bcf)], timeout)
    _run([bcftools, "index", "-f", "-t", str(vcf)], timeout)
    tool_map = output_dir / "explicit_ibd_common.map"
    # Hap-IBD et Refined IBD lisent la carte PLINK quatre colonnes sans en-tête.
    tool_map.write_text("".join(f"{row['CHROMOSOME']}\t{row['VARIANT_ID']}\t{row['POSITION_CM']}\t{row['POSITION_BP']}\n" for row in included), encoding="utf-8")
    return vcf, tool_map, audit


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Prépare un univers commun, appelle deux méthodes et publie les contrats."""
    started_at, clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path); validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = stage_inputs["parameters"]
    if parameters.get("method") != METHOD_ID:
        raise ExplicitIbdInputError("explicit_ibd_method_invalid")
    primary = _scenario(parameters["primary"])
    sensitivities = tuple(_scenario(raw) for raw in parameters["sensitivities"])
    if primary.role != "PRIMARY" or any(item.role != "SENSITIVITY" for item in sensitivities) or len({primary.scenario_id, *(item.scenario_id for item in sensitivities)}) != 1 + len(sensitivities):
        raise ExplicitIbdInputError("explicit_ibd_scenario_registry_invalid")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    adapters = validate_adapters(config["tools"]["explicit_ibd_adapters"])
    ids = ("config_input_target_variant_metadata", "samples_master", "cohorts_frozen", "target_genetic_map", "shapeit5_final_bcf", "shapeit5_final_index", "carrier_haplotypes", "founder_analysis_summary")
    artifacts = {name: _artifact(stage_inputs, name) for name in ids}
    paths = {name: _path(item, run_dir) for name, item in artifacts.items()}
    target = yaml.safe_load(paths["config_input_target_variant_metadata"].read_text(encoding="utf-8")); validate_json_document(target, "target_variant_metadata.schema.json")
    samples = validate_tsv_table(paths["samples_master"], "samples_master.schema.json").rows
    family_by_sample = {row["SAMPLE_ID"]: row["FID"] for row in samples}
    carriers = validate_tsv_table(paths["carrier_haplotypes"], "carrier_haplotypes.schema.json").rows
    mutant: dict[str, set[str]] = {}
    for row in carriers:
        if row["RELIABILITY_STATUS"] != "PASS" or row["ALT_COPY_COUNT"] == "0":
            continue
        family = family_by_sample[row["SAMPLE_ID"]]
        mutant.setdefault(family, set()).update(("H1", "H2") if row["CARRIER_HAPLOTYPE"] == "BOTH" else (row["CARRIER_HAPLOTYPE"],))
    families = tuple(sorted(mutant))
    map_rows = validate_tsv_table(paths["target_genetic_map"], "target_genetic_map.schema.json").rows
    query_rows = _query_complete_variants(config["tools"]["bcftools"], paths["shapeit5_final_bcf"], parameters["tool_timeout_seconds"])
    analysis_dir = output_dir / "explicit_ibd"; analysis_dir.mkdir(parents=True)
    vcf, tool_map, marker_audit = _prepare_tool_inputs(config["tools"]["bcftools"], paths["shapeit5_final_bcf"], analysis_dir, map_rows, query_rows, parameters["tool_timeout_seconds"])
    complete_rows = [row for row in map_rows if any(item[0] == row["VARIANT_ID"] and item[3] for item in query_rows)]
    cm_at_bp = {int(row["POSITION_BP"]): float(row["POSITION_CM"]) for row in complete_rows}
    marker_positions = tuple(sorted(cm_at_bp))
    target_bp = int(target["position_bp"])
    all_segments: list[IbdSegment] = []
    tool_records: list[dict[str, Any]] = []
    raw_by_scenario = {item["scenario_id"]: item for item in (parameters["primary"], *parameters["sensitivities"])}
    for scenario in (primary, *sensitivities):
        for tool in ("HAP_IBD", "REFINED_IBD"):
            prefix = analysis_dir / "tool_outputs" / f"{scenario.scenario_id}.{tool.lower()}"; prefix.parent.mkdir(parents=True, exist_ok=True)
            run = run_tool(tool=tool, adapters=adapters, vcf_path=vcf, map_path=tool_map, output_prefix=prefix, minimum_cm=scenario.minimum_cm, minimum_markers=scenario.minimum_markers, threads=parameters["threads"], memory_mb=parameters["java_memory_mb"], timeout_seconds=parameters["tool_timeout_seconds"])
            parsed = parse_hap_ibd(run.output_path, scenario.scenario_id, family_by_sample, cm_at_bp, marker_positions) if tool == "HAP_IBD" else parse_refined_ibd(run.output_path, scenario.scenario_id, family_by_sample, cm_at_bp, marker_positions)
            all_segments.extend(parsed); tool_records.append({"tool": tool, "scenario": scenario.scenario_id, "command": list(run.command), "version": adapters["hap_ibd_version" if tool == "HAP_IBD" else "refined_ibd_version"]})
    evaluations = []
    for scenario in (primary, *sensitivities):
        raw = raw_by_scenario[scenario.scenario_id]
        evaluations.append(evaluate_scenario(scenario, all_segments, families, {key: frozenset(value) for key, value in mutant.items()}, target_bp, parameters["boundary_tolerance_bp"], raw["control_frequency"], parameters["maximum_control_frequency"], raw["calibration_specificity_acceptable"]))
    primary_evaluation, sensitivity_evaluations = evaluations[0], tuple(evaluations[1:])
    status = classify_explicit_ibd(primary_evaluation, sensitivity_evaluations)
    segments_path = analysis_dir / "explicit_ibd_segments.tsv"
    marker_path = analysis_dir / "ibd_marker_audit.tsv"
    pair_path = analysis_dir / "target_ibd_pair_matrix.tsv"
    concordance_path = analysis_dir / "ibd_tool_concordance.tsv"
    calibration_path = analysis_dir / "ibd_calibration.tsv"
    control_path = analysis_dir / "ibd_control_frequency.tsv"
    feasibility_path = analysis_dir / "explicit_ibd_feasibility.json"
    summary_path = analysis_dir / "explicit_ibd_summary.json"
    segment_columns = ("TOOL","SCENARIO_ID","SAMPLE_1","HAPLOTYPE_1","FAMILY_1","SAMPLE_2","HAPLOTYPE_2","FAMILY_2","CHROMOSOME","START_BP","END_BP","START_CM","END_CM","LENGTH_CM","MARKER_COUNT","CONTAINS_TARGET")
    _write_tsv(segments_path, segment_columns, (_segment_row(item, target_bp) for item in all_segments))
    _write_tsv(marker_path, ("VARIANT_ORDER","VARIANT_ID","CHROMOSOME","POSITION_BP","POSITION_CM","IS_TARGET","HAP_IBD_USE","REFINED_IBD_USE","DECISION_REASON"), marker_audit)
    pair_rows = []
    for scenario, evaluation in zip((primary, *sensitivities), evaluations, strict=True):
        for family_1, family_2 in combinations(families, 2):
            selected = [item for item in evaluation.selected_segments if item.family_pair == (family_1, family_2)]
            tools = {item.tool for item in selected}
            pair_rows.append({"SCENARIO_ID": scenario.scenario_id, "ROLE": scenario.role, "FAMILY_1": family_1, "FAMILY_2": family_2, "HAP_IBD_STATUS": "CALLED" if "HAP_IBD" in tools else ("NOT_CALLED" if evaluation.evaluable else "NOT_EVALUABLE"), "REFINED_IBD_STATUS": "CALLED" if "REFINED_IBD" in tools else ("NOT_CALLED" if evaluation.evaluable else "NOT_EVALUABLE"), "BOTH_METHODS": tools == {"HAP_IBD", "REFINED_IBD"}, "CONTAINS_TARGET": bool(selected) and all(item.contains(target_bp) for item in selected), "BOUNDARIES_CONCORDANT": evaluation.boundaries_concordant, "DETAIL_CODE": evaluation.reason})
    _write_tsv(pair_path, ("SCENARIO_ID","ROLE","FAMILY_1","FAMILY_2","HAP_IBD_STATUS","REFINED_IBD_STATUS","BOTH_METHODS","CONTAINS_TARGET","BOUNDARIES_CONCORDANT","DETAIL_CODE"), pair_rows)
    _write_tsv(concordance_path, ("SCENARIO_ID","ALL_REQUIRED_PAIRS","GLOBAL_HAPLOTYPE_COHERENT","BOUNDARIES_CONCORDANT","COMMON_START_BP","COMMON_END_BP","TARGET_IN_COMMON_INTERSECTION","STATUS"), ({"SCENARIO_ID": scenario.scenario_id, "ALL_REQUIRED_PAIRS": value.all_pairs_by_both_tools, "GLOBAL_HAPLOTYPE_COHERENT": value.globally_coherent, "BOUNDARIES_CONCORDANT": value.boundaries_concordant, "COMMON_START_BP": value.common_start_bp, "COMMON_END_BP": value.common_end_bp, "TARGET_IN_COMMON_INTERSECTION": value.common_start_bp is not None and value.common_start_bp <= target_bp <= value.common_end_bp, "STATUS": _evaluation_status(value)} for scenario, value in zip((primary, *sensitivities), evaluations, strict=True)))
    calibration_rows, control_rows = [], []
    for scenario in (primary, *sensitivities):
        raw = raw_by_scenario[scenario.scenario_id]; acceptable = raw["calibrated"] and raw["calibration_specificity_acceptable"]
        calibration_rows.append({"SCENARIO_ID": scenario.scenario_id, "ROLE": scenario.role, "MINIMUM_CM": scenario.minimum_cm, "MINIMUM_MARKERS": scenario.minimum_markers, "CALIBRATED": raw["calibrated"], "SPECIFICITY_ACCEPTABLE": raw["calibration_specificity_acceptable"], "CONTROL_FREQUENCY": raw["control_frequency"], "MAXIMUM_CONTROL_FREQUENCY": parameters["maximum_control_frequency"], "STATUS": "ACCEPTABLE" if acceptable and raw["control_frequency"] is not None else ("UNACCEPTABLE" if raw["calibrated"] else "NOT_EVALUABLE")})
        control_rows.append({"SCENARIO_ID": scenario.scenario_id, "CONTROL_SOURCE": "CALIBRATION_NULL", "EVALUABLE_UNIT_COUNT": 0, "POSITIVE_UNIT_COUNT": 0, "FREQUENCY": raw["control_frequency"], "PRESPECIFIED_MAXIMUM": parameters["maximum_control_frequency"], "STATUS": "NOT_EVALUABLE" if raw["control_frequency"] is None else ("ACCEPTABLE" if raw["control_frequency"] <= parameters["maximum_control_frequency"] else "TOO_FREQUENT")})
    _write_tsv(calibration_path, ("SCENARIO_ID","ROLE","MINIMUM_CM","MINIMUM_MARKERS","CALIBRATED","SPECIFICITY_ACCEPTABLE","CONTROL_FREQUENCY","MAXIMUM_CONTROL_FREQUENCY","STATUS"), calibration_rows)
    _write_tsv(control_path, ("SCENARIO_ID","CONTROL_SOURCE","EVALUABLE_UNIT_COUNT","POSITIVE_UNIT_COUNT","FREQUENCY","PRESPECIFIED_MAXIMUM","STATUS"), control_rows)
    feasibility = {"schema_version":"1.0.0","input_scope":parameters["input_scope"],"variant_count":len(query_rows),"complete_phased_variant_count":len(complete_rows),"target_present":any(row["IS_TARGET_VARIANT"] for row in complete_rows),"primary_density_possible":len(complete_rows)>=primary.minimum_markers,"status":"GO" if len(complete_rows)>=primary.minimum_markers else ("GO_SENSITIVITY_ONLY" if len(complete_rows)>=min(item.minimum_markers for item in sensitivities) else "NOT_EVALUABLE")}
    validate_json_document(feasibility, "explicit_ibd_feasibility.schema.json"); atomic_write_json(feasibility_path, feasibility)
    summary = {"schema_version":"1.0.0","method_id":METHOD_ID,"status":status,"primary_status":_evaluation_status(primary_evaluation),"sensitivity_statuses":[{"scenario_id":scenario.scenario_id,"status":_evaluation_status(value)} for scenario,value in zip(sensitivities,sensitivity_evaluations,strict=True)],"input_scope":parameters["input_scope"],"target":{"assembly":target["assembly"],"chromosome":str(target["chromosome"]),"position_bp":target_bp,"project_variant_id":target["project_variant_id"]},"family_count":len(families),"interpretation":{"ibs_only_from_step13":True,"explicit_ibd_supported":status in {"PRIMARY_CONCORDANT_IBD_SUPPORT","SENSITIVITY_ONLY_IBD_SUPPORT"},"ibd_proven":status=="PRIMARY_CONCORDANT_IBD_SUPPORT","founder_effect_proven":False,"geographic_origin_inferred":False,"composite_score_calculated":False},"limitations":["Une absence d'appel sur puce SNP ne réfute pas un effet fondateur.","Les sensibilités sous 2 cM ne constituent pas une preuve primaire.","Un appel IBD explicite ne prouve pas à lui seul l'effet fondateur."]}
    validate_json_document(summary, "explicit_ibd_summary.schema.json"); atomic_write_json(summary_path, summary)
    contracts = ((segments_path,"explicit_ibd_segments.schema.json"),(marker_path,"explicit_ibd_marker_audit.schema.json"),(pair_path,"explicit_ibd_pair_results.schema.json"),(concordance_path,"explicit_ibd_concordance.schema.json"),(calibration_path,"explicit_ibd_calibration.schema.json"),(control_path,"explicit_ibd_control_frequency.schema.json"))
    for path, schema in contracts: validate_tsv_table(path, schema)
    specs = [("explicit_ibd_segments",segments_path,"explicit_ibd_segments.schema.json","sensitive_genetic"),("explicit_ibd_marker_audit",marker_path,"explicit_ibd_marker_audit.schema.json","internal"),("explicit_ibd_pair_results",pair_path,"explicit_ibd_pair_results.schema.json","sensitive_genetic"),("explicit_ibd_concordance",concordance_path,"explicit_ibd_concordance.schema.json","internal"),("explicit_ibd_calibration",calibration_path,"explicit_ibd_calibration.schema.json","internal"),("explicit_ibd_control_frequency",control_path,"explicit_ibd_control_frequency.schema.json","internal"),("explicit_ibd_feasibility",feasibility_path,"explicit_ibd_feasibility.schema.json","internal"),("explicit_ibd_summary",summary_path,"explicit_ibd_summary.schema.json","internal")]
    outputs = [build_file_artifact(physical_path=path,published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",artifact_id=artifact_id,artifact_type=artifact_id,media_type="application/json" if path.suffix==".json" else "text/tab-separated-values",producer_stage=stage_inputs["stage_name"],producer_signature=stage_inputs["signature"],schema_name=schema,schema_version="1.0.0",assembly=target["assembly"],sample_set_id=None,variant_set_id=target["project_variant_id"],sensitivity=sensitivity) for artifact_id,path,schema,sensitivity in specs]
    stage_outputs={"schema_version":"1.0.0","run_id":stage_inputs["run_id"],"stage_id":stage_inputs["stage_id"],"stage_name":stage_inputs["stage_name"],"signature":stage_inputs["signature"],"artifacts":outputs}; validate_json_document(stage_outputs,"stage_outputs.schema.json"); atomic_write_json(output_dir/"stage_outputs.json",stage_outputs)
    audit={"schema_version":"1.0.0","run_id":stage_inputs["run_id"],"stage_id":stage_inputs["stage_id"],"stage_name":stage_inputs["stage_name"],"method_id":METHOD_ID,"signature":stage_inputs["signature"],"started_at":started_at,"completed_at":utc_now(),"duration_seconds":monotonic()-clock,"inputs":list(artifacts.values()),"outputs":outputs,"parameters":parameters,"tools":tool_records,"counts":{"independent_families":len(families),"input_variants":len(query_rows),"complete_variants":len(complete_rows),"normalized_segments":len(all_segments)},"metrics":{"analysis_status":status,"ibd_proven":status=="PRIMARY_CONCORDANT_IBD_SUPPORT","founder_effect_proven":False,"composite_founder_score_calculated":False},"exclusions":[{"code":"true_missing_call","count":sum(not item[3] for item in query_rows)}],"warnings":[] if status=="PRIMARY_CONCORDANT_IBD_SUPPORT" else [{"code":"absence_or_non_primary_call_does_not_refute_founder_effect","count":1}],"checks":[{"check":"same_variant_universe_both_tools","status":"PASS"},{"check":"same_genetic_map_both_tools","status":"PASS"},{"check":"no_missing_alleles_hap_ibd","status":"PASS"},{"check":"global_mutant_haplotype_assignment","status":"PASS"},{"check":"ibs_ibd_founder_separation","status":"PASS"}],"known_limits":summary["limitations"],"expected_visualizations":["explicit_ibd_target_segments"],"manual_validation_required":True}; validate_json_document(audit,"stage_audit.schema.json"); atomic_write_json(output_dir/"audit.json",audit)
    (output_dir/"checksums.sha256").write_text("".join(f"{item['sha256']}  {item['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n" for item in outputs),encoding="utf-8")
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    parser=argparse.ArgumentParser(); parser.add_argument("--stage-inputs",type=Path,required=True); parser.add_argument("--output-dir",type=Path,required=True); parsed=parser.parse_args(arguments)
    try: return execute(parsed.stage_inputs,parsed.output_dir)
    except ExplicitIbdExternalError as error: sys.stderr.write(f"{error}\n"); return 3
    except (ExplicitIbdInputError,DocumentValidationError,TableValidationError,ValueError,OSError,json.JSONDecodeError) as error: sys.stderr.write(f"{error}\n"); return 2
    except Exception as error: sys.stderr.write(f"{error}\n"); return 5


if __name__ == "__main__": raise SystemExit(main())
