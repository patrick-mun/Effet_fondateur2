"""Extraction vérifiée d'une fenêtre depuis un panel phasé mis en cache."""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import DocumentValidationError, validate_json_document
from effet_fondateur.orchestrator.state import utc_now
from effet_fondateur.references.cache import CachedReferencePanel
from effet_fondateur.references.catalog import ResolvedReferencePanel


CommandRunner = Callable[[Sequence[str], float], subprocess.CompletedProcess[str]]
GRCH38_AUTOSOME_LENGTHS = {
    1: 248_956_422,
    2: 242_193_529,
    3: 198_295_559,
    4: 190_214_555,
    5: 181_538_259,
    6: 170_805_979,
    7: 159_345_973,
    8: 145_138_636,
    9: 138_394_717,
    10: 133_797_422,
    11: 135_086_622,
    12: 133_275_309,
    13: 114_364_328,
    14: 107_043_718,
    15: 101_991_189,
    16: 90_338_345,
    17: 83_257_441,
    18: 80_373_285,
    19: 58_617_616,
    20: 64_444_167,
    21: 46_709_983,
    22: 50_818_468,
}
CONTIG_PATTERN = re.compile(r"^##contig=<ID=([^,>]+),length=([0-9]+)(?:,.*)?>$")


class ReferenceWindowError(RuntimeError):
    """Signale une configuration ou une exécution d'extraction invalide."""


class ReferenceWindowIntegrityError(ReferenceWindowError):
    """Signale une référence ou une fenêtre extraite incohérente."""


@dataclass(frozen=True)
class ExtractedReferenceWindow:
    """Décrit une fenêtre VCF bgzip, son index et son manifeste vérifiés."""

    output_dir: Path
    manifest_path: Path
    vcf_path: Path
    index_path: Path
    chromosome: int
    start_bp: int
    end_bp: int
    sample_count: int
    variant_count: int
    vcf_sha256: str
    index_sha256: str


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ReferenceWindowError(f"invalid_reference_window_parameter:{name}")
    return value


def _positive_number(value: float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ReferenceWindowError(f"invalid_reference_window_parameter:{name}")
    return float(value)


def _default_runner(
    command: Sequence[str], timeout_seconds: float
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )


def _resolve_executable(command: str) -> str:
    if not isinstance(command, str) or not command.strip():
        raise ReferenceWindowError("invalid_reference_window_parameter:bcftools_command")
    resolved = shutil.which(command)
    if resolved is None:
        raise ReferenceWindowError("bcftools_not_found")
    return resolved


def _run_checked(
    runner: CommandRunner,
    command: Sequence[str],
    timeout_seconds: float,
    operation: str,
) -> subprocess.CompletedProcess[str]:
    try:
        completed = runner(command, timeout_seconds)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ReferenceWindowError(f"bcftools_{operation}_failed") from error
    if completed.returncode != 0:
        raise ReferenceWindowError(f"bcftools_{operation}_failed:{completed.returncode}")
    return completed


def _cache_manifest(
    cached: CachedReferencePanel, resolved: ResolvedReferencePanel
) -> dict[str, Any]:
    if cached.entry_dir.is_symlink() or not cached.entry_dir.is_dir():
        raise ReferenceWindowIntegrityError("reference_cache_entry_invalid")
    expected_paths = (
        cached.manifest_path,
        cached.vcf_path,
        cached.index_path,
        cached.readme_path,
        cached.source_manifest_path,
    )
    if any(path.is_symlink() or not path.is_file() for path in expected_paths):
        raise ReferenceWindowIntegrityError("reference_cache_file_invalid")
    try:
        manifest = read_json(cached.manifest_path)
        validate_json_document(manifest, "reference_cache_manifest.schema.json")
    except (ValueError, DocumentValidationError) as error:
        raise ReferenceWindowIntegrityError("reference_cache_manifest_invalid") from error
    expected_identity = {
        "cache_entry_key": cached.cache_entry_key,
        "catalog_id": resolved.catalog_id,
        "catalog_sha256": resolved.catalog_sha256,
        "panel_id": resolved.panel_id,
        "provider": resolved.provider,
        "release_id": resolved.release_id,
        "assembly": resolved.assembly,
        "chromosome": resolved.chromosome,
        "sample_count": resolved.sample_count,
        "sample_scope": resolved.sample_scope,
        "population_scope": resolved.population_scope,
    }
    for field, expected in expected_identity.items():
        if manifest.get(field) != expected:
            raise ReferenceWindowIntegrityError(f"reference_cache_identity_mismatch:{field}")
    path_by_role = {
        "vcf": cached.vcf_path,
        "index": cached.index_path,
        "readme": cached.readme_path,
        "source_manifest": cached.source_manifest_path,
    }
    for role, path in path_by_role.items():
        record = manifest["files"][role]
        if path.parent != cached.entry_dir or path.name != record["filename"]:
            raise ReferenceWindowIntegrityError(f"reference_cache_path_mismatch:{role}")
        if path.stat().st_size != record["size_bytes"]:
            raise ReferenceWindowIntegrityError(f"reference_cache_size_mismatch:{role}")
        if sha256_file(path) != record["local_sha256"]:
            raise ReferenceWindowIntegrityError(f"reference_cache_sha256_mismatch:{role}")
    return manifest


def _sample_ids(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
    operation: str,
) -> list[str]:
    completed = _run_checked(
        runner,
        [executable, "query", "--list-samples", str(vcf_path)],
        timeout_seconds,
        operation,
    )
    sample_ids = [line for line in completed.stdout.splitlines() if line]
    if len(sample_ids) != len(set(sample_ids)):
        raise ReferenceWindowIntegrityError("reference_duplicate_sample_id")
    return sample_ids


def _header(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
    operation: str,
) -> str:
    return _run_checked(
        runner,
        [executable, "view", "--header-only", str(vcf_path)],
        timeout_seconds,
        operation,
    ).stdout


def _validate_grch38_header(header: str, chromosome: int) -> None:
    expected_length = GRCH38_AUTOSOME_LENGTHS[chromosome]
    expected_contig = f"chr{chromosome}"
    observed_lengths: list[int] = []
    for line in header.splitlines():
        match = CONTIG_PATTERN.fullmatch(line)
        if match is not None and match.group(1) == expected_contig:
            observed_lengths.append(int(match.group(2)))
    if observed_lengths != [expected_length]:
        raise ReferenceWindowIntegrityError("reference_grch38_contig_mismatch")
    if not any(line.startswith("##FORMAT=<ID=GT,") for line in header.splitlines()):
        raise ReferenceWindowIntegrityError("reference_gt_format_missing")


def _variant_count(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
) -> int:
    completed = _run_checked(
        runner,
        [executable, "index", "--nrecords", str(vcf_path)],
        timeout_seconds,
        "count_variants",
    )
    try:
        count = int(completed.stdout.strip())
    except ValueError as error:
        raise ReferenceWindowIntegrityError("reference_variant_count_invalid") from error
    if count <= 0:
        raise ReferenceWindowIntegrityError("reference_window_empty")
    return count


def _validate_phased_genotypes(
    executable: str,
    vcf_path: Path,
    runner: CommandRunner,
    timeout_seconds: float,
) -> None:
    completed = _run_checked(
        runner,
        [
            executable,
            "query",
            "--include",
            'GT!="mis" && GT~"/"',
            "--format",
            "%CHROM\\t%POS\\n",
            str(vcf_path),
        ],
        timeout_seconds,
        "check_phase",
    )
    if completed.stdout.strip():
        raise ReferenceWindowIntegrityError("reference_unphased_genotype_detected")


def _result_from_manifest(output_dir: Path, manifest: dict[str, Any]) -> ExtractedReferenceWindow:
    return ExtractedReferenceWindow(
        output_dir=output_dir,
        manifest_path=output_dir / "reference_window_manifest.json",
        vcf_path=output_dir / manifest["files"]["vcf"]["filename"],
        index_path=output_dir / manifest["files"]["index"]["filename"],
        chromosome=manifest["chromosome"],
        start_bp=manifest["region"]["start_bp"],
        end_bp=manifest["region"]["end_bp"],
        sample_count=manifest["sample_count"],
        variant_count=manifest["variant_count"],
        vcf_sha256=manifest["files"]["vcf"]["sha256"],
        index_sha256=manifest["files"]["index"]["sha256"],
    )


def extract_reference_window(
    resolved: ResolvedReferencePanel,
    cached: CachedReferencePanel,
    output_dir: Path,
    *,
    start_bp: int,
    end_bp: int,
    bcftools_command: str = "bcftools",
    threads: int = 1,
    timeout_seconds: float = 7_200,
    command_runner: CommandRunner = _default_runner,
) -> ExtractedReferenceWindow:
    """Extrait et valide une fenêtre en conservant tous les échantillons du panel."""
    start = _positive_integer(start_bp, "start_bp")
    end = _positive_integer(end_bp, "end_bp")
    resolved_threads = _positive_integer(threads, "threads")
    timeout = _positive_number(timeout_seconds, "timeout_seconds")
    if resolved.assembly != "GRCh38":
        raise ReferenceWindowError("unsupported_reference_assembly")
    if resolved.chromosome not in GRCH38_AUTOSOME_LENGTHS:
        raise ReferenceWindowError("unsupported_reference_chromosome")
    if start > end or end > GRCH38_AUTOSOME_LENGTHS[resolved.chromosome]:
        raise ReferenceWindowError("invalid_reference_window_bounds")
    if output_dir.exists() or output_dir.is_symlink():
        raise ReferenceWindowError("reference_window_output_exists")
    cache_manifest = _cache_manifest(cached, resolved)
    executable = _resolve_executable(bcftools_command)
    version_result = _run_checked(
        command_runner, [executable, "--version"], timeout, "version"
    )
    version_line = version_result.stdout.splitlines()[0] if version_result.stdout else ""
    if not version_line.startswith("bcftools "):
        raise ReferenceWindowIntegrityError("bcftools_version_invalid")
    source_header = _header(
        executable, cached.vcf_path, command_runner, timeout, "read_source_header"
    )
    _validate_grch38_header(source_header, resolved.chromosome)
    source_samples = _sample_ids(
        executable, cached.vcf_path, command_runner, timeout, "read_source_samples"
    )
    if len(source_samples) != resolved.sample_count:
        raise ReferenceWindowIntegrityError("reference_source_sample_count_mismatch")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    published = False
    try:
        vcf_name = f"reference.chr{resolved.chromosome}.{start}-{end}.vcf.gz"
        vcf_path = staging_dir / vcf_name
        index_path = staging_dir / f"{vcf_name}.tbi"
        region = f"chr{resolved.chromosome}:{start}-{end}"
        _run_checked(
            command_runner,
            [
                executable,
                "view",
                "--no-update",
                "--regions",
                region,
                "--output-type",
                "z",
                "--threads",
                str(resolved_threads),
                "--output",
                str(vcf_path),
                str(cached.vcf_path),
            ],
            timeout,
            "extract_window",
        )
        if vcf_path.is_symlink() or not vcf_path.is_file():
            raise ReferenceWindowIntegrityError("reference_window_vcf_missing")
        _run_checked(
            command_runner,
            [executable, "index", "--tbi", "--threads", str(resolved_threads), str(vcf_path)],
            timeout,
            "index_window",
        )
        if index_path.is_symlink() or not index_path.is_file():
            raise ReferenceWindowIntegrityError("reference_window_index_missing")
        output_header = _header(
            executable, vcf_path, command_runner, timeout, "read_window_header"
        )
        _validate_grch38_header(output_header, resolved.chromosome)
        output_samples = _sample_ids(
            executable, vcf_path, command_runner, timeout, "read_window_samples"
        )
        if output_samples != source_samples:
            raise ReferenceWindowIntegrityError("reference_window_sample_set_or_order_mismatch")
        variant_count = _variant_count(executable, vcf_path, command_runner, timeout)
        _validate_phased_genotypes(executable, vcf_path, command_runner, timeout)
        manifest = {
            "schema_version": "1.0.0",
            "created_at": utc_now(),
            "method_id": "bcftools_all_samples_region_extract_v1",
            "assembly": resolved.assembly,
            "chromosome": resolved.chromosome,
            "region": {"start_bp": start, "end_bp": end},
            "sample_count": len(output_samples),
            "variant_count": variant_count,
            "sample_scope": resolved.sample_scope,
            "population_scope": resolved.population_scope,
            "source": {
                "catalog_id": resolved.catalog_id,
                "catalog_sha256": resolved.catalog_sha256,
                "panel_id": resolved.panel_id,
                "provider": resolved.provider,
                "release_id": resolved.release_id,
                "cache_entry_key": cached.cache_entry_key,
                "cache_manifest_sha256": sha256_file(cached.manifest_path),
                "vcf_sha256": cache_manifest["files"]["vcf"]["local_sha256"],
                "index_sha256": cache_manifest["files"]["index"]["local_sha256"],
            },
            "tool": {"name": "bcftools", "version": version_line.removeprefix("bcftools ")},
            "files": {
                "vcf": {
                    "filename": vcf_path.name,
                    "sha256": sha256_file(vcf_path),
                    "size_bytes": vcf_path.stat().st_size,
                },
                "index": {
                    "filename": index_path.name,
                    "sha256": sha256_file(index_path),
                    "size_bytes": index_path.stat().st_size,
                },
            },
            "checks": {
                "assembly": "PASS",
                "chromosome": "PASS",
                "sample_count": "PASS",
                "sample_set_and_order": "PASS",
                "phased_genotypes": "PASS",
                "tabix_index": "PASS",
            },
        }
        validate_json_document(manifest, "reference_window_manifest.schema.json")
        atomic_write_json(staging_dir / "reference_window_manifest.json", manifest)
        os.replace(staging_dir, output_dir)
        published = True
        return _result_from_manifest(output_dir, manifest)
    except (DocumentValidationError, OSError) as error:
        raise ReferenceWindowError("reference_window_publication_failed") from error
    finally:
        if not published:
            shutil.rmtree(staging_dir, ignore_errors=True)
