"""Matérialisation privée D11 depuis un panel 1000G phasé vérifié."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

from effet_fondateur.calibration.partitioning import membership_sha256
from effet_fondateur.calibration.simulation.base_haplotypes import (
    BASE_HAPLOTYPE_METHOD_ID,
    BASE_HAPLOTYPE_METHOD_VERSION,
    BaseHaplotypeAudit,
    BaseHaplotypeError,
    canonical_provenance_bytes,
    pair_donor_haplotypes,
    parse_genetic_map,
    parse_population_assignments,
    transform_phased_vcf,
    validate_base_haplotype_manifest,
)


_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_ALLOWED_PREFIX_BY_ROLE = {
    "CALIBRATION": "SYN_CAL_",
    "VALIDATION": "SYN_VAL_",
}
_PARTITION_ROLE_BY_PANEL = {
    "CALIBRATION": "DONOR_CALIBRATION",
    "VALIDATION": "DONOR_VALIDATION",
}
_TECHNICAL_GENETIC_MAP_ID = "shapeit4_hapmap_grch38_7a0cab7"
_CHECKS = [
    "SOURCE_HASHES_MATCHED",
    "PARTITION_BINDING_MATCHED",
    "HAPLOTYPE_COPY_BIJECTION",
    "SOURCE_PAIRS_DISJOINT",
    "ALLELE_FREQUENCIES_EXACT",
    "HAPLOTYPE_LD_EXACT_BY_PERMUTATION",
    "GENOTYPE_HAPLOTYPE_CONCORDANCE_100_PERCENT",
    "GENETIC_MAP_COORDINATES_VALID",
    "RAW_DONOR_IDENTIFIERS_NOT_IN_MANIFEST",
]
CommandRunner = Callable[[list[str]], None]


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise BaseHaplotypeError("base_haplotype_json_duplicate_key")
        result[key] = value
    return result


def _load_json(path: Path, error_code: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise BaseHaplotypeError(error_code)
    try:
        loaded = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BaseHaplotypeError(error_code) from error
    if not isinstance(loaded, dict):
        raise BaseHaplotypeError(error_code)
    return loaded


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_public_file(path: Path, error_code: str) -> None:
    if not path.is_file() or path.is_symlink():
        raise BaseHaplotypeError(error_code)
    lowered_parts = {part.casefold() for part in path.resolve().parts}
    if "dock6_real" in lowered_parts or {"data", "input"} <= lowered_parts:
        raise BaseHaplotypeError("base_haplotype_real_input_forbidden")
    current = path.parent
    while current != current.parent:
        if current.is_symlink():
            raise BaseHaplotypeError("base_haplotype_source_symlink_forbidden")
        current = current.parent


def _catalog_source(catalog: Mapping[str, Any], source_id: str) -> Mapping[str, Any]:
    sources = catalog.get("sources")
    if not isinstance(sources, list):
        raise BaseHaplotypeError("base_haplotype_public_catalog_invalid")
    matches = [
        source
        for source in sources
        if isinstance(source, Mapping) and source.get("source_id") == source_id
    ]
    if len(matches) != 1:
        raise BaseHaplotypeError("base_haplotype_public_catalog_source_missing")
    return matches[0]


def _technical_panel(catalog: Mapping[str, Any]) -> Mapping[str, Any]:
    panels = catalog.get("panels")
    if not isinstance(panels, list):
        raise BaseHaplotypeError("base_haplotype_technical_catalog_invalid")
    matches = [
        panel
        for panel in panels
        if isinstance(panel, Mapping)
        and panel.get("panel_id") == "1kg_3202_high_coverage_20220422"
    ]
    if len(matches) != 1:
        raise BaseHaplotypeError("base_haplotype_technical_panel_missing")
    return matches[0]


def _verify_source_files(
    source_vcf_path: Path,
    source_index_path: Path,
    technical_catalog: Mapping[str, Any],
    chromosome: int,
) -> tuple[str, str]:
    _validate_public_file(source_vcf_path, "base_haplotype_source_vcf_invalid")
    _validate_public_file(source_index_path, "base_haplotype_source_index_invalid")
    panel = _technical_panel(technical_catalog)
    expected_filename = str(panel.get("vcf_filename_template", "")).format(
        chromosome=chromosome
    )
    expected_index_filename = str(panel.get("index_filename_template", "")).format(
        chromosome=chromosome
    )
    if (
        source_vcf_path.name != expected_filename
        or source_index_path.name != expected_index_filename
    ):
        raise BaseHaplotypeError("base_haplotype_source_filename_mismatch")
    chromosomes = panel.get("chromosomes")
    if not isinstance(chromosomes, list):
        raise BaseHaplotypeError("base_haplotype_technical_catalog_invalid")
    records = [
        record
        for record in chromosomes
        if isinstance(record, Mapping) and record.get("chromosome") == chromosome
    ]
    if len(records) != 1:
        raise BaseHaplotypeError("base_haplotype_chromosome_catalog_missing")
    if _md5_file(source_vcf_path) != records[0].get("vcf_md5"):
        raise BaseHaplotypeError("base_haplotype_source_vcf_md5_mismatch")
    if _md5_file(source_index_path) != records[0].get("index_md5"):
        raise BaseHaplotypeError("base_haplotype_source_index_md5_mismatch")
    return _sha256_file(source_vcf_path), _sha256_file(source_index_path)


def _verify_catalog_bindings(
    public_catalog: Mapping[str, Any],
    technical_catalog_path: Path,
) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    panel_source = _catalog_source(
        public_catalog, "igsr_1kg_3202_phased_grch38_20220422"
    )
    map_source = _catalog_source(
        public_catalog, "shapeit4_hapmap_grch38_git_7a0cab7"
    )
    panel_binding = panel_source.get("catalog_binding")
    if (
        not isinstance(panel_binding, Mapping)
        or _sha256_file(technical_catalog_path) != panel_binding.get("sha256")
    ):
        raise BaseHaplotypeError("base_haplotype_technical_catalog_sha256_mismatch")
    return panel_source, map_source


def _verify_genetic_map(
    genetic_map_path: Path,
    genetic_map_manifest_path: Path,
    map_source: Mapping[str, Any],
    chromosome: int,
) -> str:
    _validate_public_file(genetic_map_path, "base_haplotype_genetic_map_invalid")
    manifest = _load_json(
        genetic_map_manifest_path, "base_haplotype_genetic_map_manifest_invalid"
    )
    fingerprint = map_source.get("fingerprint")
    catalog_binding = map_source.get("catalog_binding")
    archive = manifest.get("archive")
    maps = manifest.get("maps")
    if (
        manifest.get("map_id") != _TECHNICAL_GENETIC_MAP_ID
        or not isinstance(fingerprint, Mapping)
        or not isinstance(catalog_binding, Mapping)
        or not isinstance(archive, Mapping)
        or archive.get("sha256") != fingerprint.get("sha256")
        or manifest.get("catalog_sha256") != catalog_binding.get("sha256")
        or not isinstance(maps, Mapping)
    ):
        raise BaseHaplotypeError("base_haplotype_genetic_map_binding_mismatch")
    map_record = maps.get(str(chromosome))
    map_sha256 = _sha256_file(genetic_map_path)
    if (
        not isinstance(map_record, Mapping)
        or map_record.get("path") != f"normalized/chr{chromosome}.grch38.tsv"
        or map_record.get("sha256") != map_sha256
    ):
        raise BaseHaplotypeError("base_haplotype_genetic_map_sha256_mismatch")
    return map_sha256


def _read_membership(
    partition_root: Path,
    panel_role: str,
) -> tuple[dict[str, Any], tuple[str, ...], str, str]:
    manifest_path = partition_root / "private_truth/reference_partition_manifest.json"
    manifest = _load_json(manifest_path, "base_haplotype_partition_manifest_invalid")
    partitions = manifest.get("partitions")
    if not isinstance(partitions, list):
        raise BaseHaplotypeError("base_haplotype_partition_manifest_invalid")
    partition_role = _PARTITION_ROLE_BY_PANEL[panel_role]
    matches = [
        partition
        for partition in partitions
        if isinstance(partition, dict) and partition.get("role") == partition_role
    ]
    if len(matches) != 1:
        raise BaseHaplotypeError("base_haplotype_partition_role_missing")
    partition = matches[0]
    relative_path = partition.get("membership_artifact_path")
    if not isinstance(relative_path, str):
        raise BaseHaplotypeError("base_haplotype_membership_path_invalid")
    membership_path = partition_root / relative_path
    if (
        not membership_path.is_file()
        or membership_path.is_symlink()
        or partition_root.resolve() not in membership_path.resolve().parents
    ):
        raise BaseHaplotypeError("base_haplotype_membership_path_invalid")
    identifiers = tuple(
        line.strip()
        for line in membership_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    calculated_sha256 = membership_sha256(identifiers)
    if (
        calculated_sha256 != partition.get("membership_sha256")
        or len(identifiers) != partition.get("member_count")
    ):
        raise BaseHaplotypeError("base_haplotype_membership_binding_mismatch")
    return manifest, identifiers, _sha256_file(manifest_path), calculated_sha256


def _run_command(command: list[str]) -> None:
    try:
        completed = subprocess.run(
            command,
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            timeout=3600,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise BaseHaplotypeError("base_haplotype_external_command_failed") from error
    if completed.returncode != 0:
        executable = Path(command[0]).name
        raise BaseHaplotypeError(
            "base_haplotype_external_command_failed",
            f"{executable} exit={completed.returncode}",
        )


def _write_deterministic_gzip(path: Path, payload: bytes) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            mtime=0,
        ) as gzip_handle:
            gzip_handle.write(payload)


def _artifact(
    artifact_id: str,
    role: str,
    relative_path: str,
    absolute_path: Path,
    record_count: int,
) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "role": role,
        "path": relative_path,
        "sha256": _sha256_file(absolute_path),
        "size_bytes": absolute_path.stat().st_size,
        "record_count": record_count,
    }


def materialize_base_haplotypes(
    *,
    public_catalog_path: Path,
    technical_catalog_path: Path,
    source_vcf_path: Path,
    source_index_path: Path,
    population_metadata_path: Path,
    genetic_map_path: Path,
    genetic_map_manifest_path: Path,
    partition_root: Path,
    workspace_root: Path,
    run_id: str,
    panel_role: str,
    chromosome: int,
    start_bp: int,
    end_bp: int,
    random_seed: int,
    created_at: str,
    bcftools_executable: str = "bcftools",
    command_runner: CommandRunner = _run_command,
) -> BaseHaplotypeAudit:
    """Produit atomiquement les quatre artefacts privés et le manifest D11.

    Les VCF publics restent à leur emplacement explicite. Le manifest final ne
    mémorise que leur nom de fichier et leurs empreintes ; les identifiants
    donneurs sont confinés au fichier de provenance privé.
    """
    if _STABLE_IDENTIFIER.fullmatch(run_id) is None:
        raise BaseHaplotypeError("base_haplotype_run_id_invalid")
    if panel_role not in _ALLOWED_PREFIX_BY_ROLE:
        raise BaseHaplotypeError("base_haplotype_panel_role_invalid")
    if not isinstance(chromosome, int) or not 1 <= chromosome <= 22:
        raise BaseHaplotypeError("base_haplotype_chromosome_invalid")
    if (
        not isinstance(start_bp, int)
        or not isinstance(end_bp, int)
        or start_bp < 1
        or end_bp < start_bp
    ):
        raise BaseHaplotypeError("base_haplotype_domain_invalid")
    if workspace_root.name != "calibration_workspace" or workspace_root.is_symlink():
        raise BaseHaplotypeError("base_haplotype_workspace_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise BaseHaplotypeError("base_haplotype_run_exists")

    public_catalog = _load_json(
        public_catalog_path, "base_haplotype_public_catalog_invalid"
    )
    technical_catalog = _load_json(
        technical_catalog_path, "base_haplotype_technical_catalog_invalid"
    )
    panel_source, map_source = _verify_catalog_bindings(
        public_catalog, technical_catalog_path
    )
    map_source = _catalog_source(
        public_catalog, "shapeit4_hapmap_grch38_git_7a0cab7"
    )
    metadata_source = _catalog_source(
        public_catalog, "igsr_1kg_3202_population_metadata_20130606"
    )
    _validate_public_file(
        population_metadata_path, "base_haplotype_population_metadata_invalid"
    )
    metadata_fingerprint = metadata_source.get("fingerprint")
    if (
        not isinstance(metadata_fingerprint, Mapping)
        or _sha256_file(population_metadata_path) != metadata_fingerprint.get("sha256")
    ):
        raise BaseHaplotypeError("base_haplotype_population_metadata_sha256_mismatch")
    genetic_map_sha256 = _verify_genetic_map(
        genetic_map_path,
        genetic_map_manifest_path,
        map_source,
        chromosome,
    )
    source_vcf_sha256, source_index_sha256 = _verify_source_files(
        source_vcf_path,
        source_index_path,
        technical_catalog,
        chromosome,
    )
    partition_manifest, donor_identifiers, partition_manifest_sha256, membership_sha = (
        _read_membership(partition_root, panel_role)
    )
    with population_metadata_path.open(encoding="utf-8") as file_handle:
        populations = parse_population_assignments(file_handle)
    pairings = pair_donor_haplotypes(
        donor_identifiers,
        populations,
        random_seed=random_seed,
        synthetic_identifier_prefix=_ALLOWED_PREFIX_BY_ROLE[panel_role],
    )
    with genetic_map_path.open(encoding="utf-8") as file_handle:
        genetic_map_points = parse_genetic_map(
            file_handle, expected_chromosome=chromosome
        )

    temporary_root = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private_root = temporary_root / "private_truth/d11"
        work_root = temporary_root / ".work"
        private_root.mkdir(parents=True, mode=0o700)
        work_root.mkdir(mode=0o700)
        membership_path = work_root / "donors.ids.txt"
        membership_path.write_text(
            "".join(f"{identifier}\n" for identifier in donor_identifiers),
            encoding="utf-8",
        )
        subset_path = work_root / "source.filtered.vcf"
        region = f"chr{chromosome}:{start_bp}-{end_bp}"
        command_runner(
            [
                bcftools_executable,
                "view",
                "--samples-file",
                str(membership_path),
                "--regions",
                region,
                "--min-alleles",
                "2",
                "--max-alleles",
                "2",
                "--types",
                "snps",
                "--apply-filters",
                "PASS,.",
                "--exclude",
                "F_MISSING>0",
                "--output-type",
                "v",
                "--output-file",
                str(subset_path),
                str(source_vcf_path),
            ]
        )
        if not subset_path.is_file():
            raise BaseHaplotypeError("base_haplotype_subset_missing")

        phased_plain = work_root / "base_haplotypes.phased.vcf"
        unphased_plain = work_root / "base_genotypes.unphased.vcf"
        map_output = private_root / "genetic_map.tsv"
        with (
            subset_path.open(encoding="utf-8") as source_handle,
            phased_plain.open("w", encoding="utf-8", newline="") as phased_handle,
            unphased_plain.open("w", encoding="utf-8", newline="") as unphased_handle,
            map_output.open("w", encoding="utf-8", newline="") as map_handle,
        ):
            transformation_stats = transform_phased_vcf(
                source_handle,
                phased_handle,
                unphased_handle,
                map_handle,
                pairings=pairings,
                genetic_map_points=genetic_map_points,
                chromosome=chromosome,
                start_bp=start_bp,
                end_bp=end_bp,
            )
        phased_output = private_root / "base_haplotypes.phased.vcf.gz"
        unphased_output = private_root / "base_genotypes.unphased.vcf.gz"
        for plain_path, compressed_path in (
            (phased_plain, phased_output),
            (unphased_plain, unphased_output),
        ):
            command_runner(
                [
                    bcftools_executable,
                    "view",
                    "--output-type",
                    "z",
                    "--output-file",
                    str(compressed_path),
                    str(plain_path),
                ]
            )
            if not compressed_path.is_file():
                raise BaseHaplotypeError("base_haplotype_compressed_output_missing")
        provenance_output = private_root / "donor_provenance.tsv.gz"
        _write_deterministic_gzip(
            provenance_output, canonical_provenance_bytes(pairings)
        )
        for path in (
            phased_output,
            unphased_output,
            provenance_output,
            map_output,
        ):
            path.chmod(0o600)

        partition_records = partition_manifest.get("partitions", [])
        partition_role = _PARTITION_ROLE_BY_PANEL[panel_role]
        partition_record = next(
            record
            for record in partition_records
            if isinstance(record, Mapping) and record.get("role") == partition_role
        )
        manifest = {
            "schema_version": "1.0.0",
            "base_haplotype_manifest_id": f"{run_id}.d11",
            "manifest_state": "FINALIZED",
            "created_at": created_at,
            "panel": {
                "panel_id": f"{run_id}.{panel_role.casefold()}",
                "role": panel_role,
                "synthetic_identifier_prefix": _ALLOWED_PREFIX_BY_ROLE[panel_role],
                "individual_count": len(pairings),
                "haplotype_count": 2 * len(pairings),
            },
            "source": {
                "source_id": panel_source["source_id"],
                "release_id": panel_source["version"],
                "assembly": "GRCh38",
                "chromosome": chromosome,
                "vcf_filename": source_vcf_path.name,
                "vcf_sha256": source_vcf_sha256,
                "index_sha256": source_index_sha256,
                "genetic_map_source_id": map_source["source_id"],
                "genetic_map_sha256": genetic_map_sha256,
            },
            "partition": {
                "partition_manifest_id": partition_manifest["partition_manifest_id"],
                "partition_manifest_sha256": partition_manifest_sha256,
                "partition_role": partition_role,
                "membership_sha256": membership_sha,
                "donor_count": len(donor_identifiers),
            },
            "method": {
                "method_id": BASE_HAPLOTYPE_METHOD_ID,
                "method_version": BASE_HAPLOTYPE_METHOD_VERSION,
                "random_seed": random_seed,
                "pairing_unit": "PUBLIC_SOURCE_POPULATION",
                "haplotype_usage": "EACH_DONOR_HAPLOTYPE_EXACTLY_ONCE",
                "same_source_pair_policy": "FORBIDDEN",
                "variant_policy": "BIALLELIC_SNV_COMPLETE_PHASED",
                "genotype_derivation": "UNORDERED_PAIR_OF_TRUE_HAPLOTYPES",
            },
            "domain": {
                "chromosome": chromosome,
                "start_bp": transformation_stats.first_position_bp,
                "end_bp": transformation_stats.last_position_bp,
                "base_position_unit": "BP_1_BASED_CLOSED",
                "genetic_distance_unit": "cM",
                "genetic_map_interpolation": "PIECEWISE_LINEAR_NO_EXTRAPOLATION",
            },
            "confidentiality": {
                "raw_donor_identifiers_in_manifest": False,
                "local_source_paths_in_manifest": False,
                "donor_provenance_private": True,
                "artifacts_version_control_allowed": False,
                "aggregate_counts_publishable": True,
            },
            "artifacts": [
                _artifact(
                    f"{run_id}.phased",
                    "PHASED_HAPLOTYPES",
                    "private_truth/d11/base_haplotypes.phased.vcf.gz",
                    phased_output,
                    transformation_stats.variant_count,
                ),
                _artifact(
                    f"{run_id}.unphased",
                    "UNPHASED_GENOTYPES",
                    "private_truth/d11/base_genotypes.unphased.vcf.gz",
                    unphased_output,
                    transformation_stats.variant_count,
                ),
                _artifact(
                    f"{run_id}.provenance",
                    "DONOR_PROVENANCE",
                    "private_truth/d11/donor_provenance.tsv.gz",
                    provenance_output,
                    2 * len(pairings),
                ),
                _artifact(
                    f"{run_id}.map",
                    "GENETIC_MAP",
                    "private_truth/d11/genetic_map.tsv",
                    map_output,
                    transformation_stats.variant_count,
                ),
            ],
            "integrity": {
                "source_donor_count": len(donor_identifiers),
                "synthetic_individual_count": len(pairings),
                "source_haplotype_count": 2 * len(donor_identifiers),
                "synthetic_haplotype_count": 2 * len(pairings),
                "variant_count": transformation_stats.variant_count,
                "population_stratum_count": len(
                    {pairing.source_population for pairing in pairings}
                ),
                "haplotype_copy_error_count": 0,
                "same_source_pair_count": 0,
                "allele_count_difference": transformation_stats.allele_count_difference,
                "genotype_haplotype_discordance_count": transformation_stats.genotype_haplotype_discordance_count,
                "genetic_map_coordinate_error_count": transformation_stats.genetic_map_coordinate_error_count,
                "genetic_map_domain_covered": True,
                "mendelian_evaluation": "DEFERRED_TO_D13_NO_PEDIGREE",
                "checks": _CHECKS,
            },
        }
        audit = validate_base_haplotype_manifest(
            manifest,
            pairings=pairings,
            transformation_stats=transformation_stats,
        )
        manifest_path = private_root / "base_haplotype_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        manifest_path.chmod(0o600)
        shutil.rmtree(work_root)
        temporary_root.rename(target_root)
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    return audit


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Matérialise les haplotypes de base D11 sans afficher d'identifiant."
    )
    parser.add_argument("--public-catalog", type=Path, required=True)
    parser.add_argument("--technical-catalog", type=Path, required=True)
    parser.add_argument("--source-vcf", type=Path, required=True)
    parser.add_argument("--source-index", type=Path, required=True)
    parser.add_argument("--population-metadata", type=Path, required=True)
    parser.add_argument("--genetic-map", type=Path, required=True)
    parser.add_argument("--genetic-map-manifest", type=Path, required=True)
    parser.add_argument("--partition-root", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--panel-role", choices=sorted(_ALLOWED_PREFIX_BY_ROLE), required=True)
    parser.add_argument("--chromosome", type=int, required=True)
    parser.add_argument("--start-bp", type=int, required=True)
    parser.add_argument("--end-bp", type=int, required=True)
    parser.add_argument("--random-seed", type=int, required=True)
    parser.add_argument("--bcftools", default="bcftools")
    return parser


def main() -> int:
    """Point d'entrée D11 ; imprime seulement les dimensions et contrôles."""
    arguments = _parser().parse_args()
    created_at = datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )
    audit = materialize_base_haplotypes(
        public_catalog_path=arguments.public_catalog,
        technical_catalog_path=arguments.technical_catalog,
        source_vcf_path=arguments.source_vcf,
        source_index_path=arguments.source_index,
        population_metadata_path=arguments.population_metadata,
        genetic_map_path=arguments.genetic_map,
        genetic_map_manifest_path=arguments.genetic_map_manifest,
        partition_root=arguments.partition_root,
        workspace_root=arguments.workspace_root,
        run_id=arguments.run_id,
        panel_role=arguments.panel_role,
        chromosome=arguments.chromosome,
        start_bp=arguments.start_bp,
        end_bp=arguments.end_bp,
        random_seed=arguments.random_seed,
        created_at=created_at,
        bcftools_executable=arguments.bcftools,
    )
    print(
        json.dumps(
            {
                "base_haplotype_manifest_id": audit.base_haplotype_manifest_id,
                "panel_role": audit.panel_role,
                "individual_count": audit.individual_count,
                "haplotype_count": audit.haplotype_count,
                "variant_count": audit.variant_count,
                "population_stratum_count": audit.population_stratum_count,
                "checks": list(audit.checks),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
