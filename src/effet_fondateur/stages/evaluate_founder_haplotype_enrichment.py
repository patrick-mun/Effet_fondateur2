"""Étape 16B : rareté empirique d'un partage IBS exact centré cible."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from itertools import chain
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

import yaml

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError, TableValidationError, build_file_artifact,
    load_pipeline_config, validate_json_document, validate_tsv_table,
)
from effet_fondateur.founder_enrichment import (
    HaplotypeProfile, Marker, enumerate_internal_null, evaluate_exact_sharing,
    sample_external_null, summarize_null,
)
from effet_fondateur.founder_enrichment.observed import (
    EnrichmentAnalysisError, distinct_background_count,
    validate_independent_family_profiles,
)
from effet_fondateur.founder_enrichment.publication import (
    NULL_DRAW_COLUMNS, null_draw_rows, write_tsv,
)
from effet_fondateur.orchestrator.state import utc_now


METHOD_ID = "target_centered_empirical_haplotype_sharing_v1"


class FounderEnrichmentInputError(ValueError):
    """Signale une configuration, provenance ou entrée 16B invalide."""


class FounderEnrichmentExternalError(RuntimeError):
    """Signale un échec contrôlé de bcftools."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="evaluate-founder-haplotype-enrichment")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _parameters(raw: dict[str, Any]) -> dict[str, Any]:
    expected = {
        "method", "primary_statistic", "minimum_independent_units",
        "minimum_flank_markers", "internal_null_mode", "external_null_draws",
        "external_null_max_attempts", "random_seed",
        "minimum_evaluable_null_draws", "empirical_classification_threshold",
        "run_superpopulation_sensitivities", "bcftools_timeout_seconds",
    }
    if set(raw) != expected:
        raise FounderEnrichmentInputError("founder_enrichment_parameter_set_invalid")
    if raw["method"] != METHOD_ID or raw["primary_statistic"] != "total_shared_cm" or raw["internal_null_mode"] != "exhaustive":
        raise FounderEnrichmentInputError("founder_enrichment_method_invalid")
    integer_names = (
        "minimum_independent_units", "minimum_flank_markers", "external_null_draws",
        "external_null_max_attempts", "minimum_evaluable_null_draws",
        "bcftools_timeout_seconds",
    )
    if any(isinstance(raw[name], bool) or not isinstance(raw[name], int) or raw[name] < 1 for name in integer_names):
        raise FounderEnrichmentInputError("founder_enrichment_integer_parameter_invalid")
    if raw["minimum_independent_units"] < 3 or raw["random_seed"] < 0:
        raise FounderEnrichmentInputError("founder_enrichment_units_or_seed_invalid")
    if raw["external_null_max_attempts"] < raw["external_null_draws"] or raw["minimum_evaluable_null_draws"] > raw["external_null_draws"]:
        raise FounderEnrichmentInputError("founder_enrichment_draw_counts_incoherent")
    threshold = raw["empirical_classification_threshold"]
    if threshold is not None and (isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not 0 < threshold <= 1):
        raise FounderEnrichmentInputError("founder_enrichment_threshold_invalid")
    if not isinstance(raw["run_superpopulation_sensitivities"], bool):
        raise FounderEnrichmentInputError("founder_enrichment_sensitivity_flag_invalid")
    return dict(raw)


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [item for item in stage_inputs["artifacts"] if item["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise FounderEnrichmentInputError(f"artifact_missing_or_ambiguous:{artifact_id}")
    return matches[0]


def _validated_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    raw_path = Path(artifact["path"])
    path = run_dir / raw_path if PurePosixPath(artifact["path"]).parts[:1] == ("stages",) else (raw_path if raw_path.is_absolute() else Path.cwd() / raw_path)
    if path.is_symlink() or not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise FounderEnrichmentInputError(f"artifact_integrity_failure:{artifact['artifact_id']}")
    return path


def _run_bcftools(command: str, arguments: list[str], timeout: int) -> str:
    try:
        completed = subprocess.run([command, *arguments], capture_output=True, text=True, check=False, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise FounderEnrichmentExternalError("bcftools_unavailable_or_timeout") from error
    if completed.returncode != 0:
        raise FounderEnrichmentExternalError(f"bcftools_failed:{completed.returncode}")
    return completed.stdout


def _sample_ids(command: str, path: Path, timeout: int) -> tuple[str, ...]:
    samples = tuple(line for line in _run_bcftools(command, ["query", "--list-samples", str(path)], timeout).splitlines() if line)
    if not samples or len(samples) != len(set(samples)):
        raise FounderEnrichmentInputError("vcf_sample_set_invalid")
    return samples


def _query_records(command: str, path: Path, sample_ids: tuple[str, ...], timeout: int) -> list[tuple[str, int, str, str, str, tuple[str, ...]]]:
    arguments = ["query", "--samples", ",".join(sample_ids), "--force-samples", "--format", "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n", str(path)]
    records: list[tuple[str, int, str, str, str, tuple[str, ...]]] = []
    for line in _run_bcftools(command, arguments, timeout).splitlines():
        fields = line.split("\t")
        if len(fields) != 5 + len(sample_ids) or "," in fields[4]:
            raise FounderEnrichmentInputError("vcf_query_record_invalid")
        records.append((fields[0].removeprefix("chr"), int(fields[1]), fields[2], fields[3], fields[4], tuple(fields[5:])))
    return records


def _split_gt(genotype: str, reverse: bool = False) -> tuple[str | None, str | None]:
    if "|" not in genotype:
        return None, None
    values = genotype.split("|")
    if len(values) != 2 or any(value not in {"0", "1", "."} for value in values):
        return None, None
    converted = tuple(None if value == "." else (str(1 - int(value)) if reverse else value) for value in values)
    return converted[0], converted[1]


def _study_panel(command: str, path: Path, map_path: Path, timeout: int) -> tuple[tuple[Marker, ...], dict[str, tuple[tuple[str | None, ...], tuple[str | None, ...]]], dict[str, tuple[str, str]]]:
    samples = _sample_ids(command, path, timeout)
    records = _query_records(command, path, samples, timeout)
    by_id = {record[2]: record for record in records}
    if len(by_id) != len(records):
        raise FounderEnrichmentInputError("study_variant_id_duplicate")
    map_rows = validate_tsv_table(map_path, "target_genetic_map.schema.json").rows
    selected_rows = [row for row in map_rows if row["VARIANT_ID"] in by_id]
    markers = tuple(Marker(row["VARIANT_ID"], int(row["POSITION_BP"]), float(row["POSITION_CM"]), bool(row["IS_TARGET_VARIANT"])) for row in selected_rows)
    if len(markers) != len(records) or sum(marker.is_target for marker in markers) != 1:
        raise FounderEnrichmentInputError("study_map_variant_set_mismatch")
    record_by_marker = [by_id[marker.variant_id] for marker in markers]
    alleles_by_sample: dict[str, tuple[tuple[str | None, ...], tuple[str | None, ...]]] = {}
    for sample_index, sample_id in enumerate(samples):
        split = [_split_gt(record[5][sample_index]) for record in record_by_marker]
        alleles_by_sample[sample_id] = (tuple(value[0] for value in split), tuple(value[1] for value in split))
    allele_pairs = {marker.variant_id: (record[3], record[4]) for marker, record in zip(markers, record_by_marker, strict=True)}
    return markers, alleles_by_sample, allele_pairs


def _reference_profiles(command: str, path: Path, sample_ids: tuple[str, ...], markers: tuple[Marker, ...], study_alleles: dict[str, tuple[str, str]], timeout: int) -> tuple[HaplotypeProfile, ...]:
    records = _query_records(command, path, sample_ids, timeout)
    records_by_locus: dict[tuple[str, int], list[tuple[str, int, str, str, str, tuple[str, ...]]]] = defaultdict(list)
    for record in records:
        records_by_locus[(record[0], record[1])].append(record)
    aligned: list[tuple[tuple[str, int, str, str, str, tuple[str, ...]] | None, bool]] = []
    chromosome = records[0][0] if records else ""
    for marker in markers:
        if marker.is_target:
            aligned.append((None, False))
            continue
        expected_ref, expected_alt = study_alleles[marker.variant_id]
        candidates = [record for record in records_by_locus[(chromosome, marker.position_bp)] if (record[3], record[4]) in {(expected_ref, expected_alt), (expected_alt, expected_ref)}]
        if len(candidates) != 1:
            raise FounderEnrichmentInputError("reference_marker_missing_or_ambiguous")
        record = candidates[0]
        aligned.append((record, (record[3], record[4]) == (expected_alt, expected_ref)))
    profiles: list[HaplotypeProfile] = []
    for sample_index, sample_id in enumerate(sample_ids):
        pairs = [((None, None) if record is None else _split_gt(record[5][sample_index], reverse)) for record, reverse in aligned]
        profiles.extend((
            HaplotypeProfile(sample_id, "H1", tuple(pair[0] for pair in pairs)),
            HaplotypeProfile(sample_id, "H2", tuple(pair[1] for pair in pairs)),
        ))
    return tuple(profiles)


def _profile(alleles_by_sample: dict[str, tuple[tuple[str | None, ...], tuple[str | None, ...]]], sample_id: str, haplotype_id: str, family_id: str | None = None) -> HaplotypeProfile:
    try:
        first, second = alleles_by_sample[sample_id]
    except KeyError as error:
        raise FounderEnrichmentInputError("selected_study_sample_missing") from error
    if haplotype_id == "H1":
        alleles = first
    elif haplotype_id == "H2":
        alleles = second
    elif haplotype_id == "BOTH":
        alleles = tuple(left if left == right else None for left, right in zip(first, second, strict=True))
    else:
        raise FounderEnrichmentInputError("selected_haplotype_invalid")
    return HaplotypeProfile(sample_id, haplotype_id, alleles, family_id)


def _decimal(value: float | None) -> str | None:
    return None if value is None else f"{value:.12g}"


def _classification(probability: float | None, threshold: float | None, evaluable: bool, multiple_backgrounds: bool) -> str:
    if not evaluable:
        return "NOT_EVALUATED"
    if multiple_backgrounds:
        return "MULTIPLE_CARRIER_BACKGROUNDS"
    if threshold is None:
        return "NOT_CLASSIFIED"
    return "UNUSUAL_TARGET_CENTERED_SHARING" if probability is not None and probability <= threshold else "NO_UNUSUAL_SHARING_DETECTED"


def _load_target_metadata(path: Path) -> dict[str, Any]:
    """Charge les métadonnées cible YAML/JSON et valide leur contrat."""
    try:
        target = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise FounderEnrichmentInputError("invalid_target_variant_metadata") from error
    if not isinstance(target, dict):
        raise FounderEnrichmentInputError("invalid_target_variant_metadata")
    validate_json_document(target, "target_variant_metadata.schema.json")
    return target


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Valide les producteurs, calcule 16B et publie sans modifier aucun run antérieur."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path); validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = _parameters(stage_inputs["parameters"])
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    artifact_ids = (
        "config_input_target_variant_metadata", "samples_master", "cohorts_frozen", "target_genetic_map",
        "shapeit5_final_bcf", "shapeit5_final_index", "carrier_haplotypes", "harmonized_reference_vcf",
        "harmonized_reference_index", "reference_harmonization_manifest", "founder_segments", "founder_consensus",
        "founder_sharing_matrix", "founder_analysis_summary", "ancestry_scores", "ancestry_variant_audit", "reference_ancestry_summary",
    )
    artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_path(artifact, run_dir) for artifact_id, artifact in artifacts.items()}
    target = _load_target_metadata(paths["config_input_target_variant_metadata"])
    founder_summary = read_json(paths["founder_analysis_summary"]); validate_json_document(founder_summary, "founder_analysis_summary.schema.json")
    ancestry_summary = read_json(paths["reference_ancestry_summary"]); validate_json_document(ancestry_summary, "reference_ancestry_summary.schema.json")
    harmonization = read_json(paths["reference_harmonization_manifest"]); validate_json_document(harmonization, "reference_harmonization_manifest.schema.json")
    if founder_summary["method_id"] != "target_centered_exact_ibs_v1" or founder_summary["status"] != "SUPPORTED_IBS_CANDIDATE":
        raise FounderEnrichmentInputError("step13_primary_method_or_status_invalid")
    if any(target[key] != ancestry_summary["target"][key] for key in ("chromosome", "position_bp", "ref", "alt")) or target["project_variant_id"] != ancestry_summary["target"]["variant_id"]:
        raise FounderEnrichmentInputError("target_identity_mismatch")
    if harmonization["assembly"] != target["assembly"] or harmonization["target_variant_id"] != target["project_variant_id"]:
        raise FounderEnrichmentInputError("reference_provenance_mismatch")
    bcftools = config["tools"]["bcftools"]
    if not isinstance(bcftools, str) or not bcftools:
        raise FounderEnrichmentInputError("bcftools_not_configured")
    version = _run_bcftools(bcftools, ["--version"], parameters["bcftools_timeout_seconds"]).splitlines()[0]

    markers, study_alleles, allele_pairs = _study_panel(bcftools, paths["shapeit5_final_bcf"], paths["target_genetic_map"], parameters["bcftools_timeout_seconds"])
    target_marker = next(marker for marker in markers if marker.is_target)
    if target_marker.variant_id != target["project_variant_id"] or target_marker.position_bp != target["position_bp"] or allele_pairs[target_marker.variant_id] != (target["ref"], target["alt"]):
        raise FounderEnrichmentInputError("target_bcf_identity_mismatch")
    sample_rows = validate_tsv_table(paths["samples_master"], "samples_master.schema.json").rows
    family_by_sample = {row["SAMPLE_ID"]: row["FID"] for row in sample_rows}
    cohort_rows = validate_tsv_table(paths["cohorts_frozen"], "cohorts_frozen.schema.json").rows
    carrier_rows = validate_tsv_table(paths["carrier_haplotypes"], "carrier_haplotypes.schema.json").rows
    segment_rows = validate_tsv_table(paths["founder_segments"], "founder_segments.schema.json").rows
    consensus = validate_tsv_table(paths["founder_consensus"], "founder_consensus.schema.json").rows[0]
    validate_tsv_table(paths["founder_sharing_matrix"], "founder_sharing_matrix.schema.json")
    validate_tsv_table(paths["ancestry_variant_audit"], "ancestry_variant_audit.schema.json")
    selected_rows = [row for row in segment_rows if row["SEGMENT_STATUS"] == "INCLUDED"]
    if len(selected_rows) < parameters["minimum_independent_units"] or len({row["FAMILY_ID"] for row in selected_rows}) != len(selected_rows):
        raise FounderEnrichmentInputError("independent_family_units_invalid")
    representatives = tuple(_profile(study_alleles, row["SAMPLE_ID"], row["CARRIER_HAPLOTYPE_ID"], row["FAMILY_ID"]) for row in selected_rows)
    validate_independent_family_profiles(representatives, parameters["minimum_independent_units"])
    observed = evaluate_exact_sharing(markers, representatives, minimum_flank_markers=parameters["minimum_flank_markers"])
    if observed.evaluation_status != "EVALUATED" or any(abs(actual - float(expected)) > 1e-9 for actual, expected in ((observed.left_shared_cm, consensus["LEFT_LENGTH_CM"]), (observed.right_shared_cm, consensus["RIGHT_LENGTH_CM"]))) or observed.left_marker_count != int(consensus["LEFT_MARKER_COUNT"]) or observed.right_marker_count != int(consensus["RIGHT_MARKER_COUNT"]):
        raise FounderEnrichmentInputError("step13_observed_statistic_not_reproduced")

    qc_samples = {row["SAMPLE_ID"] for row in cohort_rows if row["COHORT_ID"] == "target_chromosome_all_qc" and row["INCLUDED"]}
    reliable_noncarriers = [row for row in carrier_rows if row["SAMPLE_ID"] in qc_samples and row["RELIABILITY_STATUS"] == "PASS" and row["ALT_COPY_COUNT"] == "0"]
    internal_profiles = tuple(_profile(study_alleles, row["SAMPLE_ID"], haplotype_id, family_by_sample[row["SAMPLE_ID"]]) for row in reliable_noncarriers for haplotype_id in ("H1", "H2"))
    internal_draws = enumerate_internal_null(markers, internal_profiles, unit_count=len(representatives), minimum_flank_markers=parameters["minimum_flank_markers"])

    ancestry_rows = validate_tsv_table(paths["ancestry_scores"], "ancestry_scores.schema.json").rows
    reference_rows = [row for row in ancestry_rows if row["ANALYSIS_SCOPE"] == "LOCAL" and row["ENTITY_TYPE"] == "REFERENCE_HAPLOTYPE" and row["REFERENCE_INCLUDED"]]
    reference_ids = tuple(sorted({row["SAMPLE_ID"] for row in reference_rows}))
    if len(reference_ids) != 2504 or len(reference_rows) != 5008 or ancestry_summary["local"]["reference_entity_count"] != 5008:
        raise FounderEnrichmentInputError("reference_2504_individuals_5008_haplotypes_required")
    reference_profiles = _reference_profiles(bcftools, paths["harmonized_reference_vcf"], reference_ids, markers, allele_pairs, parameters["bcftools_timeout_seconds"])
    external_draws_by_stratum = {"ALL": sample_external_null(markers, reference_profiles, unit_count=len(representatives), evaluable_draws=parameters["external_null_draws"], max_attempts=parameters["external_null_max_attempts"], random_seed=parameters["random_seed"], minimum_flank_markers=parameters["minimum_flank_markers"])}
    if parameters["run_superpopulation_sensitivities"]:
        superpopulation_by_sample = {row["SAMPLE_ID"]: row["SUPERPOPULATION"] for row in reference_rows}
        for stratum_index, superpopulation in enumerate(sorted({value for value in superpopulation_by_sample.values() if value is not None}), start=1):
            stratum_profiles = tuple(profile for profile in reference_profiles if superpopulation_by_sample[profile.individual_id] == superpopulation)
            if len({profile.individual_id for profile in stratum_profiles}) < len(representatives):
                continue
            external_draws_by_stratum[superpopulation] = sample_external_null(markers, stratum_profiles, unit_count=len(representatives), evaluable_draws=parameters["external_null_draws"], max_attempts=parameters["external_null_max_attempts"], random_seed=parameters["random_seed"] + stratum_index, minimum_flank_markers=parameters["minimum_flank_markers"], stratum=superpopulation)
    external_draws = external_draws_by_stratum["ALL"]

    left, right = observed.left_bound_index, observed.right_bound_index
    assert left is not None and right is not None and observed.total_shared_cm is not None
    signature_indexes = tuple(index for index in range(left, right + 1) if not markers[index].is_target)
    representative_signature = tuple(representatives[0].alleles[index] for index in signature_indexes)
    consistency_rows = []
    multiple_backgrounds = False
    for family_id in sorted({family_by_sample[row["SAMPLE_ID"]] for row in carrier_rows if row["ALT_COPY_COUNT"] != "0" and row["RELIABILITY_STATUS"] == "PASS"}):
        mutant_profiles = []
        for row in carrier_rows:
            if family_by_sample[row["SAMPLE_ID"]] != family_id or row["ALT_COPY_COUNT"] == "0" or row["RELIABILITY_STATUS"] != "PASS":
                continue
            copies = ("H1", "H2") if row["CARRIER_HAPLOTYPE"] == "BOTH" else (row["CARRIER_HAPLOTYPE"],)
            mutant_profiles.extend(_profile(study_alleles, row["SAMPLE_ID"], copy, family_id) for copy in copies)
        compatible = sum(tuple(profile.alleles[index] for index in signature_indexes) == representative_signature for profile in mutant_profiles)
        backgrounds = distinct_background_count(mutant_profiles, signature_indexes)
        family_multiple = backgrounds > 1
        multiple_backgrounds |= family_multiple
        consistency_rows.append({"FAMILY_ID": family_id, "EXPLICIT_MUTANT_COPY_COUNT": len(mutant_profiles), "EVALUABLE_MUTANT_COPY_COUNT": len(mutant_profiles), "COMPATIBLE_COPY_COUNT": compatible, "BACKGROUND_COUNT": backgrounds, "STATUS": "MULTIPLE_CARRIER_BACKGROUNDS" if family_multiple else "CONSISTENT", "DETAIL_CODE": "SIGNATURE_DISCORDANCE" if family_multiple else None})

    internal_summary = summarize_null(internal_draws, observed.total_shared_cm, exhaustive=True)
    external_summary = summarize_null(external_draws, observed.total_shared_cm, exhaustive=False)
    minimum_null = parameters["minimum_evaluable_null_draws"]
    evaluable = (
        internal_summary.evaluable_draws >= minimum_null
        and external_summary.evaluable_draws >= minimum_null
    )
    status = _classification(external_summary.empirical_probability, parameters["empirical_classification_threshold"], evaluable, multiple_backgrounds)
    summaries = [("INTERNAL", "ALL", None, internal_summary)] + [("EXTERNAL", stratum, parameters["external_null_draws"], summarize_null(draws, observed.total_shared_cm, exhaustive=False)) for stratum, draws in external_draws_by_stratum.items()]

    analysis_dir = output_dir / "founder_enrichment"; analysis_dir.mkdir(parents=True)
    units_path = analysis_dir / "founder_haplotype_units.tsv"
    boundaries_path = analysis_dir / "founder_haplotype_boundaries.tsv"
    consistency_path = analysis_dir / "founder_haplotype_family_consistency.tsv"
    draws_path = analysis_dir / "founder_haplotype_null_draws.tsv.gz"
    variant_path = analysis_dir / "founder_haplotype_variant_audit.tsv"
    summary_table_path = analysis_dir / "founder_haplotype_enrichment_summary.tsv"
    summary_json_path = analysis_dir / "founder_haplotype_enrichment_summary.json"
    write_tsv(units_path, ("UNIT_ID", "FAMILY_ID", "SAMPLE_ID", "HAPLOTYPE_ID", "ROLE", "SELECTION_SOURCE", "STATUS", "EXCLUSION_CODE"), ({"UNIT_ID": row["INDEPENDENT_UNIT_ID"], "FAMILY_ID": row["FAMILY_ID"], "SAMPLE_ID": row["SAMPLE_ID"], "HAPLOTYPE_ID": row["CARRIER_HAPLOTYPE_ID"], "ROLE": "INDEPENDENT_CARRIER_FAMILY", "SELECTION_SOURCE": "STEP13_PRESELECTED_REPRESENTATIVE", "STATUS": "INCLUDED", "EXCLUSION_CODE": None} for row in selected_rows))
    boundary_rows = []
    for side, bound_index, shared_cm, shared_bp, marker_count in (("LEFT", left, observed.left_shared_cm, observed.left_shared_bp, observed.left_marker_count), ("RIGHT", right, observed.right_shared_cm, observed.right_shared_bp, observed.right_marker_count)):
        boundary_rows.append({"ANALYSIS_ID": "primary_exact_ibs", "SIDE": side, "BOUND_VARIANT_ID": markers[bound_index].variant_id, "BOUND_BP": markers[bound_index].position_bp, "TARGET_BP": target_marker.position_bp, "SHARED_CM": _decimal(shared_cm), "SHARED_BP": shared_bp, "MARKER_COUNT": marker_count, "STATUS": "EVALUATED", "STOP_REASON": "FIRST_MISSING_OR_DISCORDANT"})
    write_tsv(boundaries_path, ("ANALYSIS_ID", "SIDE", "BOUND_VARIANT_ID", "BOUND_BP", "TARGET_BP", "SHARED_CM", "SHARED_BP", "MARKER_COUNT", "STATUS", "STOP_REASON"), boundary_rows)
    write_tsv(consistency_path, ("FAMILY_ID", "EXPLICIT_MUTANT_COPY_COUNT", "EVALUABLE_MUTANT_COPY_COUNT", "COMPATIBLE_COPY_COUNT", "BACKGROUND_COUNT", "STATUS", "DETAIL_CODE"), consistency_rows)
    write_tsv(
        draws_path,
        NULL_DRAW_COLUMNS,
        null_draw_rows(chain(internal_draws, *external_draws_by_stratum.values())),
    )
    variant_rows = [{"VARIANT_ORDER": index + 1, "VARIANT_ID": marker.variant_id, "CHROMOSOME": target["chromosome"], "POSITION_BP": marker.position_bp, "POSITION_CM": _decimal(marker.position_cm), "IS_TARGET": marker.is_target, "OBSERVED_USE": "ANCHOR_ONLY" if marker.is_target else "TESTED", "INTERNAL_NULL_USE": "ANCHOR_ONLY" if marker.is_target else "TESTED", "EXTERNAL_NULL_USE": "ANCHOR_ONLY" if marker.is_target else "TESTED", "DECISION_REASON": "TARGET_EXCLUDED_FROM_SIGNATURE" if marker.is_target else "COMMON_EXACT_MARKER"} for index, marker in enumerate(markers)]
    write_tsv(variant_path, ("VARIANT_ORDER", "VARIANT_ID", "CHROMOSOME", "POSITION_BP", "POSITION_CM", "IS_TARGET", "OBSERVED_USE", "INTERNAL_NULL_USE", "EXTERNAL_NULL_USE", "DECISION_REASON"), variant_rows)
    summary_rows = [{"NULL_SOURCE": source, "STRATUM": stratum, "ANALYSIS_STATUS": status, "INDEPENDENT_FAMILY_COUNT": len(representatives), "OBSERVED_LEFT_CM": _decimal(observed.left_shared_cm), "OBSERVED_RIGHT_CM": _decimal(observed.right_shared_cm), "OBSERVED_TOTAL_CM": _decimal(observed.total_shared_cm), "REQUESTED_DRAWS": requested, "ATTEMPTED_DRAWS": summary.attempted_draws, "EVALUABLE_DRAWS": summary.evaluable_draws, "NON_EVALUABLE_DRAWS": summary.non_evaluable_draws, "EXCEEDANCE_COUNT": summary.exceedance_count, "EMPIRICAL_PROBABILITY": _decimal(summary.empirical_probability), "EXACT_PROBABILITY": _decimal(summary.exact_probability), "INTERVAL_LOW": _decimal(summary.interval_low), "INTERVAL_HIGH": _decimal(summary.interval_high), "CLASSIFICATION_THRESHOLD": _decimal(parameters["empirical_classification_threshold"])} for source, stratum, requested, summary in summaries]
    write_tsv(summary_table_path, ("NULL_SOURCE", "STRATUM", "ANALYSIS_STATUS", "INDEPENDENT_FAMILY_COUNT", "OBSERVED_LEFT_CM", "OBSERVED_RIGHT_CM", "OBSERVED_TOTAL_CM", "REQUESTED_DRAWS", "ATTEMPTED_DRAWS", "EVALUABLE_DRAWS", "NON_EVALUABLE_DRAWS", "EXCEEDANCE_COUNT", "EMPIRICAL_PROBABILITY", "EXACT_PROBABILITY", "INTERVAL_LOW", "INTERVAL_HIGH", "CLASSIFICATION_THRESHOLD"), summary_rows)
    summary_json = {"schema_version": "1.0.0", "method_id": METHOD_ID, "primary_statistic": "total_shared_cm", "status": status, "independent_family_count": len(representatives), "observed": {"evaluation_status": observed.evaluation_status, "left_shared_cm": observed.left_shared_cm, "right_shared_cm": observed.right_shared_cm, "total_shared_cm": observed.total_shared_cm, "left_marker_count": observed.left_marker_count, "right_marker_count": observed.right_marker_count}, "null_results": [{"source": source, "stratum": stratum, "requested_draws": requested, "attempted_draws": summary.attempted_draws, "evaluable_draws": summary.evaluable_draws, "non_evaluable_draws": summary.non_evaluable_draws, "exceedance_count": summary.exceedance_count, "empirical_probability": summary.empirical_probability, "exact_probability": summary.exact_probability, "interval_low": summary.interval_low, "interval_high": summary.interval_high} for source, stratum, requested, summary in summaries], "classification_threshold": parameters["empirical_classification_threshold"], "random_seed": parameters["random_seed"], "provenance": {"assembly": target["assembly"], "target_variant_id": target["project_variant_id"], "target_ref": target["ref"], "target_alt": target["alt"], "map_sha256": artifacts["target_genetic_map"]["sha256"], "study_bcf_sha256": artifacts["shapeit5_final_bcf"]["sha256"], "reference_vcf_sha256": artifacts["harmonized_reference_vcf"]["sha256"], "step13_summary_sha256": artifacts["founder_analysis_summary"]["sha256"], "step16a_summary_sha256": artifacts["reference_ancestry_summary"]["sha256"]}, "interpretation": {"ibs_only": True, "ibd_proven": False, "founder_effect_proven": False, "geographic_origin_inferred": False, "composite_score_calculated": False, "statement": "Partage IBS centré cible, pas preuve IBD."}}
    validate_json_document(summary_json, "founder_haplotype_enrichment_summary.schema.json"); atomic_write_json(summary_json_path, summary_json)
    # Le fichier des tirages contient plusieurs millions de lignes : valider
    # d'abord toutes les petites sorties évite une longue passe inutile si leur
    # sérialisation ou leur contrat présente une erreur.
    contracts = ((units_path, "founder_haplotype_units.schema.json"), (boundaries_path, "founder_haplotype_boundaries.schema.json"), (consistency_path, "founder_haplotype_family_consistency.schema.json"), (variant_path, "founder_haplotype_variant_audit.schema.json"), (summary_table_path, "founder_haplotype_enrichment_summary_table.schema.json"), (draws_path, "founder_haplotype_null_draws.schema.json"))
    for path, schema in contracts: validate_tsv_table(path, schema)

    specs = (("founder_haplotype_units", units_path, "founder_haplotype_units.schema.json", "sensitive_genetic"), ("founder_haplotype_boundaries", boundaries_path, "founder_haplotype_boundaries.schema.json", "sensitive_genetic"), ("founder_haplotype_family_consistency", consistency_path, "founder_haplotype_family_consistency.schema.json", "sensitive_genetic"), ("founder_haplotype_null_draws", draws_path, "founder_haplotype_null_draws.schema.json", "internal"), ("founder_haplotype_variant_audit", variant_path, "founder_haplotype_variant_audit.schema.json", "internal"), ("founder_haplotype_enrichment_summary", summary_table_path, "founder_haplotype_enrichment_summary_table.schema.json", "internal"), ("founder_haplotype_enrichment_summary_json", summary_json_path, "founder_haplotype_enrichment_summary.schema.json", "internal"))
    output_artifacts = [build_file_artifact(physical_path=path, published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}", artifact_id=artifact_id, artifact_type=artifact_id, media_type="application/gzip" if path.suffix == ".gz" else ("application/json" if path.suffix == ".json" else "text/tab-separated-values"), producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"], schema_name=schema, schema_version="1.0.0", assembly=target["assembly"], sample_set_id=None, variant_set_id=target["project_variant_id"], sensitivity=sensitivity) for artifact_id, path, schema, sensitivity in specs]
    stage_outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "signature": stage_inputs["signature"], "artifacts": output_artifacts}; validate_json_document(stage_outputs, "stage_outputs.schema.json"); atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "method_id": METHOD_ID, "signature": stage_inputs["signature"], "started_at": started_at, "completed_at": utc_now(), "duration_seconds": monotonic() - started_clock, "inputs": list(artifacts.values()), "outputs": output_artifacts, "parameters": parameters, "tools": [{"tool": "bcftools", "configured": bcftools, "version": version}], "counts": {"independent_families": len(representatives), "internal_evaluable_draws": internal_summary.evaluable_draws, "internal_non_evaluable_draws": internal_summary.non_evaluable_draws, "external_evaluable_draws": external_summary.evaluable_draws, "external_non_evaluable_draws": external_summary.non_evaluable_draws, "reference_individuals": len(reference_ids), "reference_haplotypes": len(reference_profiles)}, "metrics": {"analysis_status": status, "observed_total_cm": observed.total_shared_cm, "external_empirical_probability": external_summary.empirical_probability, "ibd_claimed": False, "founder_effect_proven": False, "composite_founder_score_calculated": False}, "exclusions": [{"code": "null_draw_not_evaluable", "count": internal_summary.non_evaluable_draws + external_summary.non_evaluable_draws}], "warnings": [] if status not in {"NOT_EVALUATED", "MULTIPLE_CARRIER_BACKGROUNDS"} else [{"code": status.lower(), "count": 1}], "checks": [{"check": "producer_artifact_integrity", "status": "PASS"}, {"check": "step13_method_and_representatives_reused", "status": "PASS"}, {"check": "target_excluded_from_signature", "status": "PASS"}, {"check": "independent_family_unit", "status": "PASS"}, {"check": "distinct_individuals_per_null_draw", "status": "PASS"}, {"check": "reference_2504_individuals_5008_haplotypes", "status": "PASS"}, {"check": "no_network_transfer", "status": "PASS"}, {"check": "ibs_ibd_separation", "status": "PASS"}], "known_limits": ["Trois familles donnent une puissance limitée et un résultat exploratoire.", "La densité de puce, le phasage et la carte conditionnent les limites.", "1000 Genomes représente imparfaitement la démographie réunionnaise.", "Un partage inhabituel reste IBS et ne prouve ni IBD ni effet fondateur."], "expected_visualizations": ["founder_haplotype_enrichment_null_distribution"], "manual_validation_required": True}; validate_json_document(audit, "stage_audit.schema.json"); atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text("".join(f"{artifact['sha256']}  {artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n" for artifact in output_artifacts), encoding="utf-8")
    return 0


def _classify_error(error: Exception) -> int:
    if isinstance(error, FounderEnrichmentExternalError): return 3
    if isinstance(error, EnrichmentAnalysisError): return 4
    if isinstance(error, (FounderEnrichmentInputError, DocumentValidationError, TableValidationError, OSError, ValueError, json.JSONDecodeError)): return 2
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    parsed = _build_parser().parse_args(arguments)
    try: return execute(parsed.stage_inputs, parsed.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return _classify_error(error)


if __name__ == "__main__":
    raise SystemExit(main())
