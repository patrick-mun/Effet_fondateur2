"""Matérialisation privée D11b de l'univers observable CytoScan 750K Accel."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, TextIO

from effet_fondateur.calibration.simulation.chip_markers import (
    CHIP_ASSEMBLY,
    CHIP_DESIGN_ID,
    CHIP_DESIGN_VERSION,
    ChipMarkerError,
    filter_vcf_to_chip_markers,
    harmonize_chip_markers,
    parse_chip_marker_manifest,
    parse_reference_vcf_variants,
    write_marker_map,
)
from effet_fondateur.calibration.simulation.chip_package import (
    NORMALIZED_MARKER_FILENAME,
    normalize_official_chip_package,
)


_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_PROVIDER_URL = "https://www.thermofisher.com/order/catalog/product/952539"
_CHIP_SOURCE_ID = "thermofisher_cytoscan750k_accel_na36_r1"
_CHECKS = (
    "OFFICIAL_CHIP_SOURCE_HASH_MATCHED",
    "OFFICIAL_PACKAGE_CONTENT_VALIDATED",
    "DENSE_D11A_BINDING_MATCHED",
    "ASSEMBLY_GRCH38_MATCHED",
    "MARKER_COORDINATES_AND_ALLELES_HARMONIZED",
    "AMBIGUOUS_MARKERS_EXCLUDED",
    "OBSERVED_VCF_MARKER_COUNT_MATCHED",
    "NO_REAL_SAMPLE_SOURCE_OR_IDENTIFIER",
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ChipMarkerError("chip_harmonization_json_duplicate_key")
        result[key] = value
    return result


def _load_json(path: Path, code: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise ChipMarkerError(code)
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ChipMarkerError(code) from error
    if not isinstance(value, dict):
        raise ChipMarkerError(code)
    return value


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_chip_catalog_source(
    public_catalog_path: Path,
    package_sha256: str,
) -> None:
    catalog = _load_json(public_catalog_path, "chip_public_catalog_invalid")
    sources = catalog.get("sources")
    scope = catalog.get("scope")
    if (
        catalog.get("catalog_state") != "FROZEN"
        or not isinstance(scope, dict)
        or scope.get("real_study_data_allowed") is not False
        or not isinstance(sources, list)
    ):
        raise ChipMarkerError("chip_public_catalog_invalid")
    matches = [
        source
        for source in sources
        if isinstance(source, dict) and source.get("source_id") == _CHIP_SOURCE_ID
    ]
    if len(matches) != 1:
        raise ChipMarkerError("chip_public_catalog_source_missing")
    source = matches[0]
    origin = source.get("origin")
    assembly = source.get("assembly")
    fingerprint = source.get("fingerprint")
    rights = source.get("rights")
    acquisition = source.get("acquisition")
    allowed_roles = source.get("allowed_roles")
    if not all(
        isinstance(value, dict)
        for value in (origin, assembly, fingerprint, rights, acquisition)
    ) or not isinstance(allowed_roles, list):
        raise ChipMarkerError("chip_public_catalog_source_invalid")
    if (
        source.get("source_type") != "ARRAY_DESIGN_MANIFEST"
        or source.get("version") != CHIP_DESIGN_VERSION
        or "OBSERVATION_MARKER_UNIVERSE" not in allowed_roles
        or origin.get("provider") != "THERMO_FISHER"
        or origin.get("canonical_url") != _PROVIDER_URL
        or assembly.get("name") != CHIP_ASSEMBLY
        or fingerprint.get("algorithm") != "SHA-256"
        or fingerprint.get("sha256") != package_sha256
        or fingerprint.get("verification") != "PROJECT_VERIFIED"
        or rights.get("review_status") != "TERMS_VERIFIED"
        or rights.get("use_decision") != "ALLOWED_LOCAL_CALIBRATION"
        or rights.get("redistribution_allowed_by_project") is not False
        or acquisition.get("local_sha256_required_before_use") is not True
    ):
        raise ChipMarkerError("chip_public_catalog_source_invalid")


def _dense_artifact(
    dense_root: Path,
) -> tuple[dict[str, Any], Path, dict[str, Any], str]:
    if (
        dense_root.is_symlink()
        or not dense_root.is_dir()
        or dense_root.parent.name != "calibration_workspace"
    ):
        raise ChipMarkerError("chip_dense_root_invalid")
    manifest_path = dense_root / "private_truth/d11/base_haplotype_manifest.json"
    manifest = _load_json(manifest_path, "chip_dense_manifest_invalid")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise ChipMarkerError("chip_dense_manifest_invalid")
    matches = [
        artifact
        for artifact in artifacts
        if isinstance(artifact, dict) and artifact.get("role") == "UNPHASED_GENOTYPES"
    ]
    if len(matches) != 1:
        raise ChipMarkerError("chip_dense_genotype_artifact_missing")
    artifact = matches[0]
    relative_path = artifact.get("path")
    if not isinstance(relative_path, str):
        raise ChipMarkerError("chip_dense_genotype_path_invalid")
    artifact_path = dense_root / relative_path
    if (
        not artifact_path.is_file()
        or artifact_path.is_symlink()
        or dense_root.resolve() not in artifact_path.resolve().parents
    ):
        raise ChipMarkerError("chip_dense_genotype_path_invalid")
    if _sha256_file(artifact_path) != artifact.get("sha256"):
        raise ChipMarkerError("chip_dense_genotype_sha256_mismatch")
    return manifest, manifest_path, artifact, _sha256_file(manifest_path)


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _write_deterministic_gzip(path: Path, payload: bytes) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            mtime=0,
        ) as gzip_handle:
            gzip_handle.write(payload)


def _validate_synthetic_vcf_header(
    path: Path,
    *,
    panel_role: str,
    expected_individual_count: int,
) -> None:
    expected_prefix = {
        "CALIBRATION": "SYN_CAL_",
        "VALIDATION": "SYN_VAL_",
    }.get(panel_role)
    if expected_prefix is None or type(expected_individual_count) is not int:
        raise ChipMarkerError("chip_dense_panel_invalid")
    with _open_text(path) as handle:
        header = next((line for line in handle if line.startswith("#CHROM\t")), None)
    if header is None:
        raise ChipMarkerError("chip_dense_vcf_header_missing")
    samples = header.rstrip("\n").split("\t")[9:]
    if (
        len(samples) != expected_individual_count
        or len(samples) != len(set(samples))
        or any(not sample.startswith(expected_prefix) for sample in samples)
    ):
        raise ChipMarkerError("chip_dense_non_synthetic_identifier_forbidden")


def _artifact(
    run_id: str,
    suffix: str,
    role: str,
    relative_path: str,
    absolute_path: Path,
    record_count: int,
) -> dict[str, Any]:
    return {
        "artifact_id": f"{run_id}.{suffix}",
        "role": role,
        "path": relative_path,
        "sha256": _sha256_file(absolute_path),
        "size_bytes": absolute_path.stat().st_size,
        "record_count": record_count,
    }


def materialize_chip_observation(
    *,
    public_catalog_path: Path,
    dense_root: Path,
    official_chip_package_path: Path,
    expected_official_package_sha256: str,
    workspace_root: Path,
    run_id: str,
    created_at: str,
) -> dict[str, Any]:
    """Lie D11a à la puce et produit un VCF observable dans un run neuf.

    Le manifeste de marqueurs est dérivé directement de la base d'annotations
    du paquet officiel NA36.r1. Les entrées d'étude et les dossiers connus
    d'exports individuels sont refusés avant lecture.
    """
    if _STABLE_IDENTIFIER.fullmatch(run_id) is None:
        raise ChipMarkerError("chip_harmonization_run_id_invalid")
    if workspace_root.name != "calibration_workspace" or workspace_root.is_symlink():
        raise ChipMarkerError("chip_harmonization_workspace_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise ChipMarkerError("chip_harmonization_run_exists")
    dense_manifest, dense_manifest_path, dense_artifact, dense_manifest_sha256 = (
        _dense_artifact(dense_root)
    )
    domain = dense_manifest.get("domain")
    panel = dense_manifest.get("panel")
    integrity = dense_manifest.get("integrity")
    source = dense_manifest.get("source")
    method = dense_manifest.get("method")
    if not all(
        isinstance(value, dict)
        for value in (domain, panel, integrity, source, method)
    ):
        raise ChipMarkerError("chip_dense_manifest_invalid")
    if (
        source.get("source_id") != "igsr_1kg_3202_phased_grch38_20220422"
        or source.get("assembly") != CHIP_ASSEMBLY
        or method.get("variant_policy") != "BIALLELIC_SNV_COMPLETE_PHASED"
    ):
        raise ChipMarkerError("chip_dense_source_contract_mismatch")
    chromosome = domain.get("chromosome")
    start_bp = domain.get("start_bp")
    end_bp = domain.get("end_bp")
    if not all(type(value) is int for value in (chromosome, start_bp, end_bp)):
        raise ChipMarkerError("chip_dense_domain_invalid")
    dense_path = dense_root / dense_artifact["path"]
    _validate_synthetic_vcf_header(
        dense_path,
        panel_role=panel.get("role"),
        expected_individual_count=panel.get("individual_count"),
    )
    with _open_text(dense_path) as dense_handle:
        variants = parse_reference_vcf_variants(
            dense_handle,
            chromosome=chromosome,
            start_bp=start_bp,
            end_bp=end_bp,
        )
    if (
        len(variants) != integrity.get("variant_count")
        or len(variants) != dense_artifact.get("record_count")
    ):
        raise ChipMarkerError("chip_dense_variant_count_mismatch")
    temporary_root = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private_root = temporary_root / "private_truth/d11b"
        private_root.mkdir(parents=True, mode=0o700)
        normalized_marker_path = private_root / NORMALIZED_MARKER_FILENAME
        normalization = normalize_official_chip_package(
            package_path=official_chip_package_path,
            expected_package_sha256=expected_official_package_sha256,
            output_path=normalized_marker_path,
        )
        _validate_chip_catalog_source(
            public_catalog_path,
            normalization.package_sha256,
        )
        with normalized_marker_path.open(encoding="utf-8") as marker_handle:
            markers = parse_chip_marker_manifest(marker_handle)
        if len(markers) != normalization.autosomal_marker_count:
            raise ChipMarkerError("chip_normalized_marker_count_mismatch")
        matches, audit = harmonize_chip_markers(
            markers,
            variants,
            chromosome=chromosome,
            start_bp=start_bp,
            end_bp=end_bp,
        )
        if not matches:
            raise ChipMarkerError("chip_harmonization_no_marker_matched")

        marker_map_path = private_root / "chip_reference_marker_map.tsv"
        with marker_map_path.open("w", encoding="utf-8", newline="") as handle:
            map_count = write_marker_map(matches, handle)

        observed_buffer_path = temporary_root / "observed_chip_genotypes.vcf"
        with (
            _open_text(dense_path) as dense_handle,
            observed_buffer_path.open("w", encoding="utf-8", newline="") as output,
        ):
            observed_count = filter_vcf_to_chip_markers(dense_handle, output, matches)
        observed_path = private_root / "observed_chip_genotypes.vcf.gz"
        _write_deterministic_gzip(observed_path, observed_buffer_path.read_bytes())
        observed_buffer_path.unlink()

        audit_payload = {
            "package_annotation_row_count": normalization.annotation_row_count,
            "package_snp_marker_count": normalization.snp_marker_count,
            "normalized_autosomal_marker_count": (
                normalization.autosomal_marker_count
            ),
            "package_non_autosomal_snp_count": (
                normalization.non_autosomal_snp_count
            ),
            "chip_marker_count_in_domain": audit.chip_marker_count,
            "reference_variant_count_in_domain": audit.reference_variant_count,
            "matched_marker_count": audit.matched_marker_count,
            "excluded_marker_count": audit.excluded_marker_count,
            "exclusion_counts": dict(audit.exclusion_counts),
            "relation_counts": dict(audit.relation_counts),
        }
        audit_path = private_root / "harmonization_audit.json"
        audit_path.write_text(
            json.dumps(audit_payload, ensure_ascii=False, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        if map_count != len(matches) or observed_count != len(matches):
            raise ChipMarkerError("chip_observed_marker_count_mismatch")
        for path in (
            normalized_marker_path,
            marker_map_path,
            observed_path,
            audit_path,
        ):
            path.chmod(0o600)

        manifest = {
            "schema_version": "1.0.0",
            "chip_harmonization_manifest_id": f"{run_id}.d11b",
            "manifest_state": "FINALIZED",
            "created_at": created_at,
            "dense_truth": {
                "base_haplotype_manifest_id": dense_manifest[
                    "base_haplotype_manifest_id"
                ],
                "base_haplotype_manifest_sha256": dense_manifest_sha256,
                "panel_role": panel["role"],
                "dense_genotype_artifact_id": dense_artifact["artifact_id"],
                "dense_genotype_sha256": dense_artifact["sha256"],
                "dense_variant_count": integrity["variant_count"],
            },
            "chip_design": {
                "array_design": CHIP_DESIGN_ID,
                "design_version": CHIP_DESIGN_VERSION,
                "assembly": CHIP_ASSEMBLY,
                "source_kind": "OFFICIAL_LIBRARY_DERIVED",
                "public_catalog_source_id": _CHIP_SOURCE_ID,
                "provider": "Thermo Fisher Scientific",
                "provider_url": _PROVIDER_URL,
                "official_package_filename": official_chip_package_path.name,
                "official_package_sha256": normalization.package_sha256,
                "normalized_marker_manifest_filename": NORMALIZED_MARKER_FILENAME,
                "normalized_marker_manifest_sha256": (
                    normalization.normalized_marker_sha256
                ),
                "sample_derived_marker_source": False,
            },
            "domain": {
                "chromosome": chromosome,
                "start_bp": start_bp,
                "end_bp": end_bp,
                "position_unit": "BP_1_BASED_CLOSED",
            },
            "policy": {
                "marker_class": "BIALLELIC_SNP_ONLY",
                "coordinate_match": "EXACT_GRCH38_CHROMOSOME_AND_POSITION",
                "allele_match": (
                    "DIRECT_SWAPPED_COMPLEMENT_OR_COMPLEMENT_SWAPPED"
                ),
                "palindromic_unknown_strand": "EXCLUDE",
                "duplicate_chip_position": "EXCLUDE_ALL",
                "multiple_reference_variants": "EXCLUDE",
                "output_allele_orientation": "PRESERVE_1000G_REFERENCE",
                "pipeline_marker_mode": "INTERSECTION",
            },
            "confidentiality": {
                "real_sample_files_read": False,
                "real_sample_identifiers_present": False,
                "local_source_paths_in_manifest": False,
                "artifacts_version_control_allowed": False,
            },
            "artifacts": [
                _artifact(
                    run_id,
                    "normalized_markers",
                    "NORMALIZED_CHIP_MARKERS",
                    (
                        "private_truth/d11b/"
                        f"{NORMALIZED_MARKER_FILENAME}"
                    ),
                    normalized_marker_path,
                    normalization.autosomal_marker_count,
                ),
                _artifact(
                    run_id,
                    "marker_map",
                    "CHIP_REFERENCE_MARKER_MAP",
                    "private_truth/d11b/chip_reference_marker_map.tsv",
                    marker_map_path,
                    len(matches),
                ),
                _artifact(
                    run_id,
                    "observed_genotypes",
                    "CHIP_OBSERVED_GENOTYPES",
                    "private_truth/d11b/observed_chip_genotypes.vcf.gz",
                    observed_path,
                    observed_count,
                ),
                _artifact(
                    run_id,
                    "harmonization_audit",
                    "HARMONIZATION_AUDIT",
                    "private_truth/d11b/harmonization_audit.json",
                    audit_path,
                    1,
                ),
            ],
            "integrity": {
                **audit_payload,
                "observed_vcf_variant_count": observed_count,
                "checks": list(_CHECKS),
            },
        }
        manifest_path = private_root / "chip_harmonization_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        manifest_path.chmod(0o600)
        temporary_root.rename(target_root)
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Matérialise l'intersection D11b CytoScan 750K Accel–1000G."
    )
    parser.add_argument("--public-catalog", type=Path, required=True)
    parser.add_argument("--dense-root", type=Path, required=True)
    parser.add_argument("--official-chip-package", type=Path, required=True)
    parser.add_argument("--expected-official-package-sha256", required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    return parser


def main() -> int:
    """Point d'entrée D11b ; n'affiche que les comptes agrégés."""
    arguments = _parser().parse_args()
    manifest = materialize_chip_observation(
        public_catalog_path=arguments.public_catalog,
        dense_root=arguments.dense_root,
        official_chip_package_path=arguments.official_chip_package,
        expected_official_package_sha256=(
            arguments.expected_official_package_sha256
        ),
        workspace_root=arguments.workspace_root,
        run_id=arguments.run_id,
        created_at=datetime.now(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z"),
    )
    print(json.dumps(manifest["integrity"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
