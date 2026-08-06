"""Exécution contrôlée de SHAPEIT5 et attribution du chromosome porteur."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import validate_json_document, validate_tsv_table
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.phasing.shapeit5 import (
    SHAPEIT5_CONTRACT,
    Shapeit5AdapterConfig,
    Shapeit5ContractError,
    build_phase_common_command,
    build_phase_rare_command,
    probe_shapeit5,
)


CommandRunner = Callable[[Sequence[str], float], subprocess.CompletedProcess[str]]
CARRIER_COLUMNS = (
    "SAMPLE_ORDER", "SAMPLE_ID", "TARGET_VARIANT_ID", "EXPLICIT_GENOTYPE",
    "PHASED_GT", "ALT_COPY_COUNT", "CARRIER_HAPLOTYPE", "PHASE_CONFIDENCE",
    "CONFIDENCE_STATUS", "RELIABILITY_STATUS", "UNRELIABLE_REASON",
)
TRANSMISSION_COLUMNS = (
    "CHILD_SAMPLE_ID", "FATHER_SAMPLE_ID", "MOTHER_SAMPLE_ID",
    "CHILD_PHASED_GT", "TRANSMISSION_STATUS", "PATERNAL_CHILD_HAPLOTYPE",
    "MATERNAL_CHILD_HAPLOTYPE",
)
UNRELIABLE_COLUMNS = (
    "CHROMOSOME", "START_BP", "END_BP", "VARIANT_ID", "REASON",
    "AFFECTED_SAMPLE_COUNT",
)


class Shapeit5ExecutionError(ValueError):
    """Signale une entrée ou un résultat SHAPEIT5 invalide."""


class Shapeit5ExecutionExternalError(RuntimeError):
    """Signale l'échec d'un exécutable SHAPEIT5 ou bcftools."""


class Shapeit5ExecutionBlockError(RuntimeError):
    """Signale un résultat scientifiquement incohérent."""


@dataclass(frozen=True)
class PhasedVariant:
    chromosome: str
    position_bp: int
    variant_id: str
    ref: str
    alt: str
    genotypes: tuple[str, ...]
    confidences: tuple[float | None, ...]


@dataclass(frozen=True)
class Shapeit5PhasingResult:
    output_dir: Path
    common_bcf_path: Path
    common_index_path: Path
    final_bcf_path: Path
    final_index_path: Path
    common_log_path: Path
    rare_log_path: Path
    carrier_haplotypes_path: Path
    transmissions_path: Path
    unreliable_regions_path: Path
    manifest_path: Path
    carrier_count: int
    reliable_carrier_count: int


def _default_runner(command: Sequence[str], timeout: float) -> subprocess.CompletedProcess[str]:
    return subprocess.run(list(command), capture_output=True, text=True, timeout=timeout, check=False)


def _run(runner: CommandRunner, command: Sequence[str], timeout: float, operation: str) -> subprocess.CompletedProcess[str]:
    try:
        result = runner(command, timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise Shapeit5ExecutionExternalError(f"{operation}_failed") from error
    if result.returncode != 0:
        raise Shapeit5ExecutionExternalError(f"{operation}_failed:{result.returncode}")
    return result


def _samples(executable: str, path: Path, runner: CommandRunner, timeout: float) -> list[str]:
    result = _run(runner, [executable, "query", "--list-samples", str(path)], timeout, "bcftools_list_phasing_samples")
    samples = [line for line in result.stdout.splitlines() if line]
    if not samples or len(samples) != len(set(samples)):
        raise Shapeit5ExecutionBlockError("phasing_sample_ids_invalid")
    return samples


def _variants(executable: str, path: Path, sample_count: int, with_confidence: bool, runner: CommandRunner, timeout: float) -> list[PhasedVariant]:
    sample_format = "[\\t%GT\\t%PP]" if with_confidence else "[\\t%GT]"
    result = _run(
        runner,
        [executable, "query", "--format", f"%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT{sample_format}\\n", str(path)],
        timeout,
        "bcftools_query_phasing_variants",
    )
    variants: list[PhasedVariant] = []
    width = 5 + sample_count * (2 if with_confidence else 1)
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != width:
            raise Shapeit5ExecutionBlockError("phasing_variant_record_malformed")
        if with_confidence:
            genotypes = tuple(fields[5::2])
            try:
                confidences = tuple(None if value in {"", "."} else float(value) for value in fields[6::2])
            except ValueError as error:
                raise Shapeit5ExecutionBlockError("phasing_confidence_invalid") from error
            if any(value is not None and not 0.5 <= value <= 1 for value in confidences):
                raise Shapeit5ExecutionBlockError("phasing_confidence_invalid")
        else:
            genotypes = tuple(fields[5:])
            confidences = tuple(None for _ in genotypes)
        variants.append(PhasedVariant(fields[0], int(fields[1]), fields[2], fields[3], fields[4], genotypes, confidences))
    if not variants:
        raise Shapeit5ExecutionBlockError("phasing_output_has_no_variants")
    return variants


def _alleles(genotype: str) -> tuple[str, str]:
    separator = "|" if "|" in genotype else "/"
    values = genotype.split(separator)
    if len(values) != 2 or any(value not in {"0", "1"} for value in values):
        raise Shapeit5ExecutionBlockError("unsupported_phasing_genotype")
    return values[0], values[1]


def _mendel_errors(variants: list[PhasedVariant], samples: list[str], pedigree: list[tuple[str, str | None, str | None]]) -> int:
    sample_index = {sample: index for index, sample in enumerate(samples)}
    errors = 0
    for child, father, mother in pedigree:
        for variant in variants:
            child_alleles = _alleles(variant.genotypes[sample_index[child]])
            father_alleles = set(_alleles(variant.genotypes[sample_index[father]])) if father else None
            mother_alleles = set(_alleles(variant.genotypes[sample_index[mother]])) if mother else None
            if father_alleles is not None and mother_alleles is not None:
                compatible = any(first in father_alleles and second in mother_alleles for first, second in (child_alleles, child_alleles[::-1]))
            else:
                known = father_alleles or mother_alleles
                compatible = known is not None and any(allele in known for allele in child_alleles)
            errors += not compatible
    return errors


def _read_pedigree(path: Path, samples: set[str]) -> list[tuple[str, str | None, str | None]]:
    rows: list[tuple[str, str | None, str | None]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 3:
            raise Shapeit5ExecutionError("shapeit5_pedigree_malformed")
        child, father, mother = fields
        normalized = (child, None if father == "NA" else father, None if mother == "NA" else mother)
        if child not in samples or any(parent not in samples for parent in normalized[1:] if parent):
            raise Shapeit5ExecutionError("shapeit5_pedigree_sample_mismatch")
        rows.append(normalized)
    return rows


def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: "" if row[column] is None else row[column] for column in columns})


def _file(path: Path) -> dict[str, Any]:
    return {"filename": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def _carrier_rows(target: PhasedVariant, samples: list[str], genotype_rows: list[dict[str, Any]], threshold: float) -> list[dict[str, Any]]:
    explicit = {row["SAMPLE_ID"]: row["GENOTYPE"] for row in genotype_rows}
    if set(explicit) != set(samples):
        raise Shapeit5ExecutionBlockError("explicit_target_genotype_sample_mismatch")
    rows: list[dict[str, Any]] = []
    for order, (sample, gt, confidence) in enumerate(zip(samples, target.genotypes, target.confidences), start=1):
        first, second = _alleles(gt)
        observed = sorted(target.ref if allele == "0" else target.alt for allele in (first, second))
        if observed != sorted(explicit[sample].split("/")):
            raise Shapeit5ExecutionBlockError("phased_target_genotype_discordant")
        copies = int(first == "1") + int(second == "1")
        haplotype = "NONE" if copies == 0 else "BOTH" if copies == 2 else "H1" if first == "1" else "H2"
        if copies != 1:
            confidence_status, reliability, reason = "NOT_APPLICABLE_HOMOZYGOUS", "PASS", None
        elif confidence is None:
            confidence_status, reliability, reason = "NOT_AVAILABLE", "UNRELIABLE", "target_phase_confidence_not_available"
        elif confidence < threshold:
            confidence_status, reliability, reason = "SCORED_LOW", "UNRELIABLE", "target_phase_confidence_below_threshold"
        else:
            confidence_status, reliability, reason = "SCORED_PASS", "PASS", None
        rows.append({
            "SAMPLE_ORDER": order, "SAMPLE_ID": sample, "TARGET_VARIANT_ID": target.variant_id,
            "EXPLICIT_GENOTYPE": explicit[sample], "PHASED_GT": gt, "ALT_COPY_COUNT": copies,
            "CARRIER_HAPLOTYPE": haplotype, "PHASE_CONFIDENCE": confidence,
            "CONFIDENCE_STATUS": confidence_status, "RELIABILITY_STATUS": reliability,
            "UNRELIABLE_REASON": reason,
        })
    return rows


def _transmissions(target: PhasedVariant, samples: list[str], pedigree: list[tuple[str, str | None, str | None]]) -> list[dict[str, Any]]:
    index = {sample: position for position, sample in enumerate(samples)}
    rows: list[dict[str, Any]] = []
    for child, father, mother in pedigree:
        child_gt = target.genotypes[index[child]]
        child_alleles = _alleles(child_gt)
        father_set = set(_alleles(target.genotypes[index[father]])) if father else None
        mother_set = set(_alleles(target.genotypes[index[mother]])) if mother else None
        direct = father_set is not None and mother_set is not None and child_alleles[0] in father_set and child_alleles[1] in mother_set
        swapped = father_set is not None and mother_set is not None and child_alleles[1] in father_set and child_alleles[0] in mother_set
        if father_set is None or mother_set is None:
            status, paternal, maternal = "DUO_COMPATIBLE", None, None
        elif direct and not swapped:
            status, paternal, maternal = "DIRECT", "H1", "H2"
        elif swapped and not direct:
            status, paternal, maternal = "SWAPPED", "H2", "H1"
        else:
            status, paternal, maternal = "AMBIGUOUS", "UNRESOLVED", "UNRESOLVED"
        rows.append({"CHILD_SAMPLE_ID": child, "FATHER_SAMPLE_ID": father, "MOTHER_SAMPLE_ID": mother, "CHILD_PHASED_GT": child_gt, "TRANSMISSION_STATUS": status, "PATERNAL_CHILD_HAPLOTYPE": paternal, "MATERNAL_CHILD_HAPLOTYPE": maternal})
    return rows


def run_shapeit5_phasing(
    *, input_manifest_path: Path, study_vcf_path: Path, reference_vcf_path: Path,
    genetic_map_path: Path, pedigree_path: Path, target_genotype_audit_path: Path,
    output_dir: Path, adapter_config: Shapeit5AdapterConfig, bcftools_command: str,
    threads: int = 1, seed: int = SHAPEIT5_CONTRACT.default_seed,
    effective_size: int = SHAPEIT5_CONTRACT.default_effective_size,
    minimum_phase_confidence: float = 0.9, timeout_seconds: float = 7200,
    command_runner: CommandRunner = _default_runner,
) -> Shapeit5PhasingResult:
    """Exécute les deux passes et attribue H1/H2 depuis le GT cible phasé."""
    if output_dir.exists() or not 0.5 <= minimum_phase_confidence <= 1:
        raise Shapeit5ExecutionError("invalid_shapeit5_execution_parameters")
    manifest = read_json(input_manifest_path)
    validate_json_document(manifest, "shapeit5_inputs_manifest.schema.json")
    input_files = manifest["files"]
    expected_hashes = (
        (study_vcf_path, input_files["study_vcf"]["sha256"]),
        (reference_vcf_path, input_files["reference_vcf"]["sha256"]),
        (genetic_map_path, input_files["genetic_map"]["sha256"]),
        (pedigree_path, input_files["pedigree"]["sha256"]),
    )
    if any(not path.is_file() or path.is_symlink() or sha256_file(path) != expected for path, expected in expected_hashes):
        raise Shapeit5ExecutionError("shapeit5_execution_input_modified")
    genotype_table = validate_tsv_table(target_genotype_audit_path, "target_genotype_audit.schema.json")
    bcftools = shutil.which(bcftools_command)
    if bcftools is None:
        raise Shapeit5ExecutionExternalError("bcftools_not_found")
    try:
        probe = probe_shapeit5(adapter_config)
    except Shapeit5ContractError as error:
        raise Shapeit5ExecutionExternalError("shapeit5_probe_failed") from error
    samples = _samples(bcftools, study_vcf_path, command_runner, timeout_seconds)
    input_variants = _variants(bcftools, study_vcf_path, len(samples), False, command_runner, timeout_seconds)
    pedigree = _read_pedigree(pedigree_path, set(samples))
    if _mendel_errors(input_variants, samples, pedigree):
        raise Shapeit5ExecutionBlockError("mendel_errors_before_phasing")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    published = False
    try:
        common_bcf, final_bcf = staging / "common.phased.bcf", staging / "target.phased.bcf"
        common_log, rare_log = staging / "common.phase.log", staging / "rare.phase.log"
        pedigree_argument = pedigree_path if pedigree else None
        common_result = _run(command_runner, build_phase_common_command(probe, input_path=study_vcf_path, reference_path=reference_vcf_path, genetic_map_path=genetic_map_path, region=manifest["scaffold_region"], output_path=common_bcf, log_path=common_log, pedigree_path=pedigree_argument, threads=threads, seed=seed), timeout_seconds, "shapeit5_phase_common")
        if not common_log.exists():
            common_log.write_text(common_result.stdout + common_result.stderr, encoding="utf-8")
        rare_result = _run(command_runner, build_phase_rare_command(probe, input_path=study_vcf_path, scaffold_path=common_bcf, genetic_map_path=genetic_map_path, input_region=manifest["input_region"], scaffold_region=manifest["scaffold_region"], output_path=final_bcf, pedigree_path=pedigree_argument, threads=threads, seed=seed, effective_size=effective_size, score_singletons=True), timeout_seconds, "shapeit5_phase_rare")
        rare_log.write_text(rare_result.stdout + rare_result.stderr, encoding="utf-8")
        for bcf_path in (common_bcf, final_bcf):
            if not bcf_path.is_file():
                raise Shapeit5ExecutionExternalError("shapeit5_output_missing")
            index_path = Path(f"{bcf_path}.csi")
            if not index_path.is_file():
                _run(command_runner, [bcftools, "index", "--csi", str(bcf_path)], timeout_seconds, "bcftools_index_shapeit5_output")
        common_samples = _samples(bcftools, common_bcf, command_runner, timeout_seconds)
        final_samples = _samples(bcftools, final_bcf, command_runner, timeout_seconds)
        if common_samples != samples or final_samples != samples:
            raise Shapeit5ExecutionBlockError("shapeit5_output_sample_order_mismatch")
        common_variants = _variants(bcftools, common_bcf, len(samples), False, command_runner, timeout_seconds)
        final_variants = _variants(bcftools, final_bcf, len(samples), True, command_runner, timeout_seconds)
        input_keys = [(v.chromosome, v.position_bp, v.variant_id, v.ref, v.alt) for v in input_variants]
        final_keys = [(v.chromosome, v.position_bp, v.variant_id, v.ref, v.alt) for v in final_variants]
        if input_keys != final_keys or len(common_variants) != manifest["common_variant_count"]:
            raise Shapeit5ExecutionBlockError("shapeit5_output_variant_mismatch")
        common_keys = [(v.chromosome, v.position_bp, v.variant_id, v.ref, v.alt) for v in common_variants]
        if len(common_keys) != len(set(common_keys)) or any(key not in input_keys for key in common_keys):
            raise Shapeit5ExecutionBlockError("shapeit5_common_scaffold_variant_mismatch")
        if any("|" not in genotype for variant in common_variants for genotype in variant.genotypes):
            raise Shapeit5ExecutionBlockError("shapeit5_common_scaffold_unphased")
        for before, after in zip(input_variants, final_variants):
            if any(sorted(_alleles(left)) != sorted(_alleles(right)) for left, right in zip(before.genotypes, after.genotypes)):
                raise Shapeit5ExecutionBlockError("shapeit5_genotype_not_preserved")
            if any("|" not in genotype for genotype in after.genotypes):
                raise Shapeit5ExecutionBlockError("shapeit5_unphased_output_genotype")
        if _mendel_errors(final_variants, samples, pedigree):
            raise Shapeit5ExecutionBlockError("mendel_errors_after_phasing")
        targets = [variant for variant in final_variants if variant.variant_id == manifest["target_variant_id"]]
        if len(targets) != 1:
            raise Shapeit5ExecutionBlockError("phased_target_missing_or_ambiguous")
        target = targets[0]
        carrier_rows = _carrier_rows(target, samples, genotype_table.rows, minimum_phase_confidence)
        transmission_rows = _transmissions(target, samples, pedigree)
        unreliable_counts: dict[str, int] = {}
        for row in carrier_rows:
            if row["UNRELIABLE_REASON"]:
                unreliable_counts[row["UNRELIABLE_REASON"]] = unreliable_counts.get(row["UNRELIABLE_REASON"], 0) + 1
        unreliable_rows = [{"CHROMOSOME": manifest["chromosome"], "START_BP": target.position_bp, "END_BP": target.position_bp, "VARIANT_ID": target.variant_id, "REASON": reason, "AFFECTED_SAMPLE_COUNT": count} for reason, count in sorted(unreliable_counts.items())]
        carrier_path, transmission_path, unreliable_path = staging / "carrier_haplotypes.tsv", staging / "phasing_transmissions.tsv", staging / "phasing_unreliable_regions.tsv"
        _write_tsv(carrier_path, CARRIER_COLUMNS, carrier_rows)
        _write_tsv(transmission_path, TRANSMISSION_COLUMNS, transmission_rows)
        _write_tsv(unreliable_path, UNRELIABLE_COLUMNS, unreliable_rows)
        validate_tsv_table(carrier_path, "carrier_haplotypes.schema.json")
        validate_tsv_table(transmission_path, "phasing_transmissions.schema.json")
        validate_tsv_table(unreliable_path, "phasing_unreliable_regions.schema.json")
        files = {path.stem.replace(".", "_"): _file(path) for path in (common_bcf, Path(f"{common_bcf}.csi"), final_bcf, Path(f"{final_bcf}.csi"), common_log, rare_log, carrier_path, transmission_path, unreliable_path)}
        carrier_count = sum(row["ALT_COPY_COUNT"] > 0 for row in carrier_rows)
        reliable_count = sum(row["ALT_COPY_COUNT"] > 0 and row["RELIABILITY_STATUS"] == "PASS" for row in carrier_rows)
        result_manifest = {"schema_version": "1.0.0", "created_at": utc_now(), "method_id": "shapeit5_common_rare_carrier_assignment_v1", "adapter_id": manifest["adapter_id"], "software_version": probe.phase_common_version, "assembly": manifest["assembly"], "chromosome": manifest["chromosome"], "input_region": manifest["input_region"], "scaffold_region": manifest["scaffold_region"], "seed": seed, "threads": threads, "effective_size": effective_size, "minimum_phase_confidence": minimum_phase_confidence, "sample_count": len(samples), "variant_count": len(final_variants), "common_variant_count": len(common_variants), "carrier_count": carrier_count, "reliable_carrier_count": reliable_count, "pedigree_record_count": len(pedigree), "target_variant_id": target.variant_id, "target_variant_role": manifest["target_variant_role"], "files": files, "checks": {"input_integrity": "PASS", "common_phase": "PASS", "rare_phase": "PASS", "sample_order": "PASS", "genotype_preservation": "PASS", "target_preservation": "PASS", "explicit_target_genotypes": "PASS", "mendel_before": "PASS", "mendel_after": "PASS", "carrier_assignment": "PASS"}}
        manifest_path = staging / "shapeit5_phasing_manifest.json"
        validate_json_document(result_manifest, "shapeit5_phasing_manifest.schema.json")
        atomic_write_json(manifest_path, result_manifest)
        os.replace(staging, output_dir)
        published = True
        return Shapeit5PhasingResult(output_dir, output_dir / common_bcf.name, output_dir / f"{common_bcf.name}.csi", output_dir / final_bcf.name, output_dir / f"{final_bcf.name}.csi", output_dir / common_log.name, output_dir / rare_log.name, output_dir / carrier_path.name, output_dir / transmission_path.name, output_dir / unreliable_path.name, output_dir / manifest_path.name, carrier_count, reliable_count)
    finally:
        if not published:
            shutil.rmtree(staging, ignore_errors=True)
