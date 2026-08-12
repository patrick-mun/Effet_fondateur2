"""Étape 16A : PCA globale et haplotypique locale sur références externes."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Iterable, Sequence

import numpy as np
import yaml

from effet_fondateur.ancestry import (
    AncestryAnalysisError,
    AncestryExtractCacheError,
    AncestryPcaError,
    AncestryReferenceError,
    GenotypePanel,
    HarmonizedPca,
    ReferenceSample,
    Variant,
    cache_ancestry_metadata,
    cache_reference_extract,
    harmonize_alt_dosages,
    load_reference_samples,
    parse_vcf_query_panel,
    population_centroids,
    read_bim_variants,
    read_plink_raw_panel,
)
from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.orchestrator.state import utc_now


MAX_COMPONENTS = 10
PC_COLUMNS = tuple(f"PC{index}" for index in range(1, MAX_COMPONENTS + 1))
SCORE_COLUMNS = (
    "ANALYSIS_SCOPE", "ENTITY_TYPE", "ENTITY_ID", "SAMPLE_ID", "HAPLOTYPE",
    "POPULATION", "SUPERPOPULATION", "TARGET_COPY_STATUS",
    "REFERENCE_INCLUDED", "PROJECTED", *PC_COLUMNS,
)
EIGENVALUE_COLUMNS = (
    "ANALYSIS_SCOPE", "COMPONENT", "EIGENVALUE", "EXPLAINED_VARIANCE_RATIO",
    "REFERENCE_ENTITY_COUNT", "INFORMATIVE_VARIANT_COUNT",
)
LOADING_COLUMNS = (
    "ANALYSIS_SCOPE", "CHROMOSOME", "POSITION_BP", "VARIANT_ID", "REF", "ALT",
    "REFERENCE_ALT_FREQUENCY", "REFERENCE_SCALE", *PC_COLUMNS,
)
AUDIT_COLUMNS = (
    "ANALYSIS_SCOPE", "CHROMOSOME", "POSITION_BP", "STUDY_VARIANT_ID",
    "STUDY_REF", "STUDY_ALT", "REFERENCE_VARIANT_ID", "REFERENCE_REF",
    "REFERENCE_ALT", "STATUS", "EXCLUSION_CODE",
)
CENTROID_COLUMNS = (
    "ANALYSIS_SCOPE", "GROUP_LEVEL", "GROUP_ID", "REFERENCE_ENTITY_COUNT", *PC_COLUMNS,
)


class AnalyzeReferenceAncestryInputError(ValueError):
    """Signale une configuration ou une entrée de l'étape 16A invalide."""


class AnalyzeReferenceAncestryExternalError(RuntimeError):
    """Signale un échec de PLINK ou bcftools."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="analyze-reference-ancestry")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise AnalyzeReferenceAncestryInputError(f"invalid_parameter:{name}")
    return value


def _probability(parameters: dict[str, Any], name: str, default: float) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AnalyzeReferenceAncestryInputError(f"invalid_parameter:{name}")
    result = float(value)
    if not math.isfinite(result) or not 0 < result <= 1:
        raise AnalyzeReferenceAncestryInputError(f"invalid_parameter:{name}")
    return result


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    if parameters.get("method", "reference_only_global_local_pca_v1") != "reference_only_global_local_pca_v1":
        raise AnalyzeReferenceAncestryInputError("invalid_parameter:method")
    result = {
        "method": "reference_only_global_local_pca_v1",
        "reference_panel_id": parameters.get("reference_panel_id", "1kg_3202_high_coverage_20220422"),
        "ancestry_cache_dir": parameters.get("ancestry_cache_dir", "data/cache/references"),
        "ancestry_cache_offline": parameters.get("ancestry_cache_offline", False),
        "global_requested_components": _positive_integer(parameters, "global_requested_components", 10),
        "local_requested_components": _positive_integer(parameters, "local_requested_components", 10),
        "minimum_global_variants": _positive_integer(parameters, "minimum_global_variants", 100),
        "minimum_local_variants": _positive_integer(parameters, "minimum_local_variants", 10),
        "minimum_reference_call_rate": _probability(parameters, "minimum_reference_call_rate", 0.95),
        "download_timeout_seconds": _positive_integer(parameters, "download_timeout_seconds", 7200),
        "bcftools_timeout_seconds": _positive_integer(parameters, "bcftools_timeout_seconds", 7200),
        "plink_timeout_seconds": _positive_integer(parameters, "plink_timeout_seconds", 300),
    }
    if (
        not isinstance(result["reference_panel_id"], str)
        or not result["reference_panel_id"]
        or not isinstance(result["ancestry_cache_dir"], str)
        or not result["ancestry_cache_dir"]
        or not isinstance(result["ancestry_cache_offline"], bool)
        or result["global_requested_components"] > MAX_COMPONENTS
        or result["local_requested_components"] > MAX_COMPONENTS
    ):
        raise AnalyzeReferenceAncestryInputError("invalid_ancestry_configuration")
    return result


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [artifact for artifact in stage_inputs["artifacts"] if artifact["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise AnalyzeReferenceAncestryInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = Path(artifact["path"])
    if not path.is_absolute() and PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        path = run_dir / path
    elif not path.is_absolute():
        path = Path.cwd() / path
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise AnalyzeReferenceAncestryInputError(f"artifact_missing_or_modified:{artifact['artifact_id']}")
    return path


def _resolve_tool(command: str | None, name: str) -> str:
    executable = shutil.which(command) if command else None
    if executable is None:
        raise AnalyzeReferenceAncestryExternalError(f"{name}_not_available")
    return executable


def _tool_version(executable: str) -> str | None:
    completed = subprocess.run([executable, "--version"], capture_output=True, text=True, check=False, timeout=10)
    lines = (completed.stdout or completed.stderr).strip().splitlines()
    return lines[0][:200] if lines else None


def _run(command: list[str], timeout: int, code: str, stdout_path: Path | None = None) -> subprocess.CompletedProcess[str] | None:
    try:
        if stdout_path is None:
            completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout)
        else:
            with stdout_path.open("x", encoding="utf-8") as handle:
                binary = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=False, timeout=timeout)
            completed = binary
    except (OSError, subprocess.TimeoutExpired) as error:
        raise AnalyzeReferenceAncestryExternalError(code) from error
    if completed.returncode != 0:
        raise AnalyzeReferenceAncestryExternalError(f"{code}:{(completed.stderr or '')[:300]}")
    return completed if stdout_path is None else None


def _samples(bcftools: str, path: Path, timeout: int) -> tuple[str, ...]:
    completed = _run([bcftools, "query", "--list-samples", str(path)], timeout, "bcftools_ancestry_samples_failed")
    assert completed is not None
    samples = tuple(line for line in completed.stdout.splitlines() if line)
    if not samples or len(samples) != len(set(samples)):
        raise AnalyzeReferenceAncestryInputError("ancestry_vcf_samples_invalid")
    return samples


def _query_panel(
    bcftools: str,
    path: Path,
    sample_ids: tuple[str, ...],
    output_path: Path,
    timeout: int,
    *,
    haplotypes: bool,
    samples_file: Path | None = None,
) -> GenotypePanel:
    command = [bcftools, "query"]
    if samples_file is not None:
        command.extend(["--samples-file", str(samples_file), "--force-samples"])
    command.extend(["--format", "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n", str(path)])
    _run(command, timeout, "bcftools_ancestry_query_failed", output_path)
    with output_path.open(encoding="utf-8") as handle:
        return parse_vcf_query_panel(
            handle,
            sample_ids,
            haplotypes=haplotypes,
            phased_required=haplotypes,
        )


def _combine_panels(panels: list[GenotypePanel]) -> GenotypePanel:
    nonempty = [panel for panel in panels if panel.variants]
    if not nonempty:
        raise AnalyzeReferenceAncestryInputError("global_reference_extract_empty")
    samples = nonempty[0].sample_ids
    if any(panel.sample_ids != samples or panel.ploidy != 2 for panel in nonempty):
        raise AnalyzeReferenceAncestryInputError("global_reference_sample_order_mismatch")
    return GenotypePanel(
        samples,
        tuple(variant for panel in nonempty for variant in panel.variants),
        np.concatenate([panel.alt_dosages for panel in nonempty], axis=1),
        2,
    )


def _load_target(path: Path) -> dict[str, Any]:
    try:
        target = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise AnalyzeReferenceAncestryInputError("target_metadata_invalid") from error
    validate_json_document(target, "target_variant_metadata.schema.json")
    return target


def _reference_panel(path: Path, panel_id: str) -> dict[str, Any]:
    document = read_json(path)
    validate_json_document(document, "reference_panel_catalog.schema.json")
    panels = [panel for panel in document["panels"] if panel["panel_id"] == panel_id]
    if len(panels) != 1:
        raise AnalyzeReferenceAncestryInputError("ancestry_reference_panel_missing_or_ambiguous")
    return panels[0]


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if value is None else str(value).lower() if isinstance(value, bool) else value for key, value in row.items()})


def _pc_values(values: np.ndarray) -> dict[str, float | None]:
    return {column: float(values[index]) if index < len(values) else None for index, column in enumerate(PC_COLUMNS)}


def _score_rows(
    scope: str,
    result: HarmonizedPca,
    metadata: dict[str, ReferenceSample],
    carrier_status: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    local = scope == "LOCAL"
    rows: list[dict[str, Any]] = []
    for entity_id, score in zip(result.reference_sample_ids, result.model.reference_scores):
        sample_id, _, haplotype = entity_id.partition(":")
        reference = metadata[sample_id]
        rows.append({
            "ANALYSIS_SCOPE": scope,
            "ENTITY_TYPE": "REFERENCE_HAPLOTYPE" if local else "REFERENCE_INDIVIDUAL",
            "ENTITY_ID": f"REF:{entity_id}", "SAMPLE_ID": sample_id,
            "HAPLOTYPE": haplotype or None, "POPULATION": reference.population,
            "SUPERPOPULATION": reference.superpopulation,
            "TARGET_COPY_STATUS": "REFERENCE_UNKNOWN" if local else "NOT_APPLICABLE",
            "REFERENCE_INCLUDED": True, "PROJECTED": False, **_pc_values(score),
        })
    for entity_id, score in zip(result.study_sample_ids, result.study_scores):
        sample_id, _, haplotype = entity_id.partition(":")
        rows.append({
            "ANALYSIS_SCOPE": scope,
            "ENTITY_TYPE": "STUDY_HAPLOTYPE" if local else "STUDY_INDIVIDUAL",
            "ENTITY_ID": f"STUDY:{entity_id}", "SAMPLE_ID": sample_id,
            "HAPLOTYPE": haplotype or None, "POPULATION": None,
            "SUPERPOPULATION": None,
            "TARGET_COPY_STATUS": carrier_status[entity_id] if carrier_status else "NOT_APPLICABLE",
            "REFERENCE_INCLUDED": False, "PROJECTED": True, **_pc_values(score),
        })
    return rows


def _eigenvalue_rows(scope: str, result: HarmonizedPca) -> list[dict[str, Any]]:
    return [{
        "ANALYSIS_SCOPE": scope, "COMPONENT": f"PC{index}", "EIGENVALUE": float(value),
        "EXPLAINED_VARIANCE_RATIO": float(value / result.model.total_variance),
        "REFERENCE_ENTITY_COUNT": len(result.reference_sample_ids),
        "INFORMATIVE_VARIANT_COUNT": len(result.variants),
    } for index, value in enumerate(result.model.eigenvalues, start=1)]


def _loading_rows(scope: str, result: HarmonizedPca) -> list[dict[str, Any]]:
    frequencies = result.model.allele_frequencies[result.model.informative_variant_mask]
    return [{
        "ANALYSIS_SCOPE": scope, "CHROMOSOME": variant.chromosome.removeprefix("chr"),
        "POSITION_BP": variant.position_bp, "VARIANT_ID": variant.variant_id,
        "REF": variant.ref, "ALT": variant.alt,
        "REFERENCE_ALT_FREQUENCY": float(frequencies[index]),
        "REFERENCE_SCALE": float(result.model.scales[index]),
        **_pc_values(result.model.loadings[index]),
    } for index, variant in enumerate(result.variants)]


def _variant_audit_rows(
    scope: str,
    reference: GenotypePanel,
    study: GenotypePanel,
    result: HarmonizedPca,
    minimum_call_rate: float,
) -> list[dict[str, Any]]:
    reference_by_locus = {variant.locus: (index, variant) for index, variant in enumerate(reference.variants)}
    informative = {variant.locus for variant in result.variants}
    rows: list[dict[str, Any]] = []
    for variant in study.variants:
        match = reference_by_locus.get(variant.locus)
        reference_variant = match[1] if match else None
        if match is None:
            status, exclusion = "ABSENT_FROM_REFERENCE", "absent_from_reference"
        elif {variant.ref, variant.alt} != {reference_variant.ref, reference_variant.alt}:
            status, exclusion = "ALLELE_MISMATCH", "allele_mismatch"
        elif variant.locus in informative:
            status, exclusion = "INFORMATIVE", None
        else:
            values = reference.alt_dosages[:, match[0]]
            if float(np.mean(~np.isnan(values))) < minimum_call_rate:
                status, exclusion = "LOW_REFERENCE_CALL_RATE", "low_reference_call_rate"
            else:
                status, exclusion = "MONOMORPHIC_REFERENCE", "monomorphic_reference"
        rows.append({
            "ANALYSIS_SCOPE": scope, "CHROMOSOME": variant.chromosome.removeprefix("chr"),
            "POSITION_BP": variant.position_bp, "STUDY_VARIANT_ID": variant.variant_id,
            "STUDY_REF": variant.ref, "STUDY_ALT": variant.alt,
            "REFERENCE_VARIANT_ID": reference_variant.variant_id if reference_variant else None,
            "REFERENCE_REF": reference_variant.ref if reference_variant else None,
            "REFERENCE_ALT": reference_variant.alt if reference_variant else None,
            "STATUS": status, "EXCLUSION_CODE": exclusion,
        })
    return rows


def _centroid_rows(scope: str, result: HarmonizedPca, metadata: dict[str, ReferenceSample]) -> list[dict[str, Any]]:
    base_samples = [entity.split(":", 1)[0] for entity in result.reference_sample_ids]
    rows: list[dict[str, Any]] = []
    for level, labels in (
        ("POPULATION", [metadata[sample].population for sample in base_samples]),
        ("SUPERPOPULATION", [metadata[sample].superpopulation for sample in base_samples]),
    ):
        for group, (count, values) in population_centroids(result.model.reference_scores, labels).items():
            rows.append({"ANALYSIS_SCOPE": scope, "GROUP_LEVEL": level, "GROUP_ID": group, "REFERENCE_ENTITY_COUNT": count, **_pc_values(values)})
    return rows


def _carrier_copy_status(path: Path, study_haplotype_ids: tuple[str, ...]) -> dict[str, str]:
    rows = validate_tsv_table(path, "carrier_haplotypes.schema.json").rows
    by_sample = {row["SAMPLE_ID"]: row for row in rows}
    if {entity.split(":", 1)[0] for entity in study_haplotype_ids} != set(by_sample):
        raise AnalyzeReferenceAncestryInputError("carrier_haplotype_sample_set_mismatch")
    result: dict[str, str] = {}
    for entity in study_haplotype_ids:
        sample, haplotype = entity.split(":", 1)
        row = by_sample[sample]
        is_carrier = row["CARRIER_HAPLOTYPE"] in {haplotype, "BOTH"}
        if is_carrier and row["RELIABILITY_STATUS"] != "PASS":
            result[entity] = "CARRIER_COPY_UNRELIABLE"
        else:
            result[entity] = "CARRIER_COPY" if is_carrier else "NON_CARRIER_COPY"
    return result


def _require_complete_local_target(
    panel: GenotypePanel, target: dict[str, Any]
) -> Variant:
    """Exige une cible unique et complète avant toute projection locale."""
    indexes = [
        index
        for index, variant in enumerate(panel.variants)
        if (
            variant.variant_id == target["project_variant_id"]
            and variant.locus == (str(target["chromosome"]), target["position_bp"])
            and (variant.ref, variant.alt) == (target["ref"], target["alt"])
        )
    ]
    if len(indexes) != 1:
        raise AnalyzeReferenceAncestryInputError("local_target_missing_or_mismatched")
    if np.isnan(np.asarray(panel.alt_dosages, dtype=float)[:, indexes[0]]).any():
        raise AnalyzeReferenceAncestryInputError("local_target_genotype_missing")
    return panel.variants[indexes[0]]


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Publie deux PCA séparées sans jamais ajuster les axes sur l'étude."""

    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "config_input_target_variant_metadata", "config_input_reference_panel_catalog",
        "config_input_ancestry_reference_catalog", "samples_master", "kinship_panel_bed",
        "kinship_panel_bim", "kinship_panel_fam", "kinship_panel_dataset",
        "shapeit5_final_bcf", "shapeit5_final_index", "carrier_haplotypes",
        "harmonized_reference_vcf", "harmonized_reference_index",
        "reference_harmonization_manifest", "roh_analysis_summary",
    )
    artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_path(artifact, run_dir) for artifact_id, artifact in artifacts.items()}
    target = _load_target(paths["config_input_target_variant_metadata"])
    if target["assembly"] != config["project"]["assembly"] or target["chromosome"] not in range(1, 23):
        raise AnalyzeReferenceAncestryInputError("ancestry_target_not_supported")
    configured_target = config["target"]
    for key in ("chromosome", "position_bp", "ref", "alt", "project_variant_id"):
        if configured_target[key] != target[key]:
            raise AnalyzeReferenceAncestryInputError("ancestry_target_config_mismatch")

    plink = _resolve_tool(config["tools"]["plink"], "plink")
    bcftools = _resolve_tool(config["tools"]["bcftools"], "bcftools")
    cache_root = Path(parameters["ancestry_cache_dir"])
    if not cache_root.is_absolute():
        cache_root = Path.cwd() / cache_root
    cached_metadata = cache_ancestry_metadata(
        catalog_path=paths["config_input_ancestry_reference_catalog"],
        cache_root=cache_root,
        offline=parameters["ancestry_cache_offline"],
        timeout_seconds=parameters["download_timeout_seconds"],
    )
    reference_samples = load_reference_samples(cached_metadata)
    if len(reference_samples) != 2504:
        raise AnalyzeReferenceAncestryInputError("unrelated_reference_count_mismatch")
    metadata = {sample.sample_id: sample for sample in reference_samples}
    reference_sample_ids = tuple(sorted(metadata))
    reference_panel = _reference_panel(paths["config_input_reference_panel_catalog"], parameters["reference_panel_id"])
    if reference_panel["assembly"] != target["assembly"]:
        raise AnalyzeReferenceAncestryInputError("reference_target_assembly_mismatch")

    descriptor = read_json(paths["kinship_panel_dataset"])
    validate_json_document(descriptor, "plink_dataset.schema.json")
    study_variants = read_bim_variants(paths["kinship_panel_bim"])
    fam_rows = [line.split() for line in paths["kinship_panel_fam"].read_text(encoding="utf-8").splitlines() if line]
    if (
        descriptor["source_format"] != "PLINK_KINSHIP_PANEL"
        or descriptor["sample_count"] != len(fam_rows)
        or descriptor["variant_count"] != len(study_variants)
        or paths["kinship_panel_bed"].read_bytes()[:3] != b"\x6c\x1b\x01"
    ):
        raise AnalyzeReferenceAncestryInputError("kinship_panel_integrity_mismatch")
    for suffix in ("bed", "bim", "fam"):
        if descriptor["files"][suffix]["sha256"] != artifacts[f"kinship_panel_{suffix}"]["sha256"]:
            raise AnalyzeReferenceAncestryInputError("kinship_panel_descriptor_checksum_mismatch")
    if len(fam_rows) == 0 or any(len(row) != 6 for row in fam_rows):
        raise AnalyzeReferenceAncestryInputError("kinship_panel_fam_invalid")
    sample_rows = validate_tsv_table(paths["samples_master"], "samples_master.schema.json").rows
    sample_by_plink = {(row["FID"], row["IID"]): row["SAMPLE_ID"] for row in sample_rows}
    if len(sample_by_plink) != len(sample_rows) or not {(row[0], row[1]) for row in fam_rows} <= set(sample_by_plink):
        raise AnalyzeReferenceAncestryInputError("kinship_panel_sample_registry_mismatch")
    harmonization_manifest = read_json(paths["reference_harmonization_manifest"])
    validate_json_document(harmonization_manifest, "reference_harmonization_manifest.schema.json")
    if (
        harmonization_manifest["assembly"] != target["assembly"]
        or harmonization_manifest["chromosome"] != target["chromosome"]
        or harmonization_manifest["target_variant_id"] != target["project_variant_id"]
    ):
        raise AnalyzeReferenceAncestryInputError("local_harmonization_target_mismatch")
    validate_json_document(read_json(paths["roh_analysis_summary"]), "roh_analysis_summary.schema.json")

    analysis_dir = output_dir / "ancestry"
    analysis_dir.mkdir()
    cache_statuses: list[str] = []
    with tempfile.TemporaryDirectory(prefix=".reference_ancestry_", dir=output_dir) as temporary_name:
        temporary = Path(temporary_name)
        samples_file = temporary / "unrelated.samples.txt"
        samples_file.write_text("".join(f"{sample}\n" for sample in reference_sample_ids), encoding="utf-8")

        raw_prefix = temporary / "study_global"
        _run([
            plink, "--bfile", str(paths["kinship_panel_bed"].with_suffix("")),
            "--recode", "A", "--out", str(raw_prefix),
        ], parameters["plink_timeout_seconds"], "plink_global_ancestry_export_failed")
        global_study = read_plink_raw_panel(raw_prefix.with_suffix(".raw"), study_variants, sample_by_plink)

        variants_by_chromosome: dict[int, list[Variant]] = defaultdict(list)
        for variant in study_variants:
            chromosome = int(variant.chromosome.removeprefix("chr"))
            variants_by_chromosome[chromosome].append(variant)
        source_by_chromosome = {row["chromosome"]: row for row in reference_panel["chromosomes"]}
        global_reference_parts: list[GenotypePanel] = []

        def extractor(source_url: str, positions: Path, samples: Path, vcf: Path, index: Path, timeout: int) -> None:
            _run([
                bcftools, "view", "--samples-file", str(samples), "--force-samples",
                "--regions-file", str(positions), "--types", "snps", "--min-alleles", "2",
                "--max-alleles", "2", "--output-type", "z", "--output", str(vcf), source_url,
            ], timeout, "bcftools_global_reference_extract_failed")
            _run([bcftools, "index", "--tbi", str(vcf)], timeout, "bcftools_global_reference_index_failed")
            generated = Path(f"{vcf}.tbi")
            if generated != index:
                generated.replace(index)

        for chromosome in sorted(variants_by_chromosome):
            positions = temporary / f"chr{chromosome}.positions.tsv"
            positions.write_text("".join(
                f"chr{chromosome}\t{variant.position_bp}\t{variant.position_bp}\n"
                for variant in sorted(variants_by_chromosome[chromosome], key=lambda value: value.position_bp)
            ), encoding="utf-8")
            source = source_by_chromosome[chromosome]
            filename = reference_panel["vcf_filename_template"].format(chromosome=chromosome)
            cached = cache_reference_extract(
                cache_root=cache_root, panel_id=reference_panel["panel_id"], assembly=reference_panel["assembly"],
                chromosome=chromosome, source_url=f"{reference_panel['base_url']}/{filename}",
                source_vcf_md5=source["vcf_md5"], source_index_md5=source["index_md5"],
                positions_path=positions, samples_path=samples_file,
                offline=parameters["ancestry_cache_offline"], timeout_seconds=parameters["bcftools_timeout_seconds"],
                extractor=extractor,
            )
            cache_statuses.append(cached.status)
            if _samples(bcftools, cached.vcf_path, parameters["bcftools_timeout_seconds"]) != reference_sample_ids:
                raise AnalyzeReferenceAncestryInputError("global_reference_extract_sample_mismatch")
            query_path = temporary / f"chr{chromosome}.query.tsv"
            global_reference_parts.append(_query_panel(
                bcftools, cached.vcf_path, reference_sample_ids, query_path,
                parameters["bcftools_timeout_seconds"], haplotypes=False,
            ))
        global_reference = _combine_panels(global_reference_parts)
        global_result = harmonize_alt_dosages(
            global_reference, global_study,
            requested_components=parameters["global_requested_components"],
            minimum_variants=parameters["minimum_global_variants"],
            minimum_reference_call_rate=parameters["minimum_reference_call_rate"],
        )

        study_local_samples = _samples(bcftools, paths["shapeit5_final_bcf"], parameters["bcftools_timeout_seconds"])
        local_study = _query_panel(
            bcftools, paths["shapeit5_final_bcf"], study_local_samples,
            temporary / "local.study.query.tsv", parameters["bcftools_timeout_seconds"], haplotypes=True,
        )
        _require_complete_local_target(local_study, target)
        local_reference_samples = _samples(bcftools, paths["harmonized_reference_vcf"], parameters["bcftools_timeout_seconds"])
        if not set(reference_sample_ids) <= set(local_reference_samples):
            raise AnalyzeReferenceAncestryInputError("local_unrelated_reference_missing")
        local_reference = _query_panel(
            bcftools, paths["harmonized_reference_vcf"], reference_sample_ids,
            temporary / "local.reference.query.tsv", parameters["bcftools_timeout_seconds"],
            haplotypes=True, samples_file=samples_file,
        )
        local_result = harmonize_alt_dosages(
            local_reference, local_study,
            requested_components=parameters["local_requested_components"],
            minimum_variants=parameters["minimum_local_variants"],
            minimum_reference_call_rate=parameters["minimum_reference_call_rate"],
        )
        carrier_status = _carrier_copy_status(paths["carrier_haplotypes"], local_result.study_sample_ids)

    score_rows = _score_rows("GLOBAL", global_result, metadata) + _score_rows("LOCAL", local_result, metadata, carrier_status)
    eigenvalue_rows = _eigenvalue_rows("GLOBAL", global_result) + _eigenvalue_rows("LOCAL", local_result)
    loading_rows = _loading_rows("GLOBAL", global_result) + _loading_rows("LOCAL", local_result)
    audit_rows = _variant_audit_rows("GLOBAL", global_reference, global_study, global_result, parameters["minimum_reference_call_rate"]) + _variant_audit_rows("LOCAL", local_reference, local_study, local_result, parameters["minimum_reference_call_rate"])
    centroid_rows = _centroid_rows("GLOBAL", global_result, metadata) + _centroid_rows("LOCAL", local_result, metadata)
    outputs = {
        "ancestry_scores": (analysis_dir / "ancestry_scores.tsv", SCORE_COLUMNS, score_rows, "ancestry_scores.schema.json", "sensitive_genetic"),
        "ancestry_eigenvalues": (analysis_dir / "ancestry_eigenvalues.tsv", EIGENVALUE_COLUMNS, eigenvalue_rows, "ancestry_eigenvalues.schema.json", "internal"),
        "ancestry_variant_loadings": (analysis_dir / "ancestry_variant_loadings.tsv", LOADING_COLUMNS, loading_rows, "ancestry_variant_loadings.schema.json", "internal"),
        "ancestry_variant_audit": (analysis_dir / "ancestry_variant_audit.tsv", AUDIT_COLUMNS, audit_rows, "ancestry_variant_audit.schema.json", "sensitive_genetic"),
        "ancestry_population_centroids": (analysis_dir / "ancestry_population_centroids.tsv", CENTROID_COLUMNS, centroid_rows, "ancestry_population_centroids.schema.json", "public"),
    }
    for path, columns, rows, schema, _ in outputs.values():
        _write_tsv(path, columns, rows)
        validate_tsv_table(path, schema)
    local_positions = [variant.position_bp for variant in local_study.variants]
    summary = {
        "schema_version": "1.0.0", "method_id": "reference_only_global_local_pca_v1",
        "assembly": target["assembly"],
        "target": {"variant_id": target["project_variant_id"], "chromosome": target["chromosome"], "position_bp": target["position_bp"], "ref": target["ref"], "alt": target["alt"]},
        "global": {"reference_entity_count": len(global_result.reference_sample_ids), "study_entity_count": len(global_result.study_sample_ids), "candidate_variant_count": global_result.candidate_variant_count, "informative_variant_count": len(global_result.variants), "component_count": len(global_result.model.eigenvalues)},
        "local": {"reference_entity_count": len(local_result.reference_sample_ids), "study_entity_count": len(local_result.study_sample_ids), "candidate_variant_count": local_result.candidate_variant_count, "informative_variant_count": len(local_result.variants), "component_count": len(local_result.model.eigenvalues), "region_start_bp": min(local_positions), "region_end_bp": max(local_positions)},
        "cache": {"metadata_status": cached_metadata.status, "extract_hits": cache_statuses.count("HIT"), "extract_populated": cache_statuses.count("POPULATED"), "offline": parameters["ancestry_cache_offline"]},
        "interpretation": {"policy": "RELATIVE_REFERENCE_POSITIONING_ONLY", "ethnic_identity_assigned": False, "genealogical_ancestor_identified": False, "local_ancestry_proven": False, "ibd_proven": False},
        "checks": {"reference_unrelated_only": "PASS", "study_not_used_for_axes": "PASS", "global_local_separated": "PASS", "target_and_region_resolved": "PASS", "cache_integrity": "PASS", "variant_harmonization": "PASS"},
    }
    validate_json_document(summary, "reference_ancestry_summary.schema.json")
    summary_path = analysis_dir / "reference_ancestry_summary.json"
    atomic_write_json(summary_path, summary)

    source = artifacts["shapeit5_final_bcf"]
    specifications = [(artifact_id, values[0], "text/tab-separated-values", values[3], values[4]) for artifact_id, values in outputs.items()]
    specifications.append(("reference_ancestry_summary", summary_path, "application/json", "reference_ancestry_summary.schema.json", "internal"))
    output_artifacts = [build_file_artifact(
        physical_path=path,
        published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
        artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
        producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
        schema_name=schema, schema_version="1.0.0", assembly=target["assembly"],
        sample_set_id=source.get("sample_set_id"), variant_set_id=source.get("variant_set_id"),
        sensitivity=sensitivity,
    ) for artifact_id, path, media_type, schema, sensitivity in specifications]
    stage_outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "signature": stage_inputs["signature"], "artifacts": output_artifacts}
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": "reference_only_global_local_pca_v1", "signature": stage_inputs["signature"], "started_at": started_at, "completed_at": utc_now(), "duration_seconds": monotonic() - started_clock,
        "inputs": list(artifacts.values()), "outputs": output_artifacts, "parameters": parameters,
        "tools": [{"tool": "python", "configured": Path(sys.executable).name, "version": sys.version.split()[0]}, {"tool": "numpy", "configured": "Python package", "version": importlib.metadata.version("numpy")}, {"tool": "plink", "configured": config["tools"]["plink"], "version": _tool_version(plink)}, {"tool": "bcftools", "configured": config["tools"]["bcftools"], "version": _tool_version(bcftools)}],
        "counts": {"unrelated_reference_individuals": len(reference_sample_ids), "global_study_individuals": len(global_result.study_sample_ids), "global_informative_variants": len(global_result.variants), "local_reference_haplotypes": len(local_result.reference_sample_ids), "local_study_haplotypes": len(local_result.study_sample_ids), "local_informative_variants": len(local_result.variants)},
        "metrics": {"analysis_status": "EXPLORATORY_REFERENCE_POSITIONING", "global_components": len(global_result.model.eigenvalues), "local_components": len(local_result.model.eigenvalues)},
        "exclusions": [{"scope": scope, "reason": reason, "count": count} for (scope, reason), count in sorted(Counter((row["ANALYSIS_SCOPE"], row["EXCLUSION_CODE"]) for row in audit_rows if row["EXCLUSION_CODE"]).items())],
        "warnings": [{"code": "reference_populations_are_broad_proxies", "count": 1}],
        "checks": [{"check": key, "status": value} for key, value in summary["checks"].items()],
        "known_limits": ["Les populations 1000 Genomes sont des repères larges et non une attribution ethnique.", "La PCA locale décrit une proximité haplotypique dans la région configurée ; elle ne prouve ni ascendance locale, ni ancêtre généalogique, ni IBD.", "Les individus de l'étude sont projetés et ne déterminent jamais les axes."],
        "expected_visualizations": ["plot_reference_ancestry_global", "plot_reference_ancestry_local"],
        "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text("".join(f"{artifact['sha256']}  {artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n" for artifact in output_artifacts), encoding="utf-8")
    return 0


def _classify_error(error: Exception) -> int:
    if isinstance(error, (AncestryAnalysisError, AncestryPcaError, AncestryReferenceError, AncestryExtractCacheError)):
        return 4
    if isinstance(error, (AnalyzeReferenceAncestryInputError, DocumentValidationError, TableValidationError, OSError, ValueError, json.JSONDecodeError)):
        return 2
    if isinstance(error, AnalyzeReferenceAncestryExternalError):
        return 3
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    parsed = _build_parser().parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return _classify_error(error)


if __name__ == "__main__":
    raise SystemExit(main())
