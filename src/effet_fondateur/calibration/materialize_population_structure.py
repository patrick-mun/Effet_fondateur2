"""Matérialisation privée D12 de la structure de population synthétique."""

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

from effet_fondateur.calibration.simulation.population_structure import (
    POPULATION_STRUCTURE_METHOD_ID,
    POPULATION_STRUCTURE_METHOD_VERSION,
    PopulationStructureError,
    build_population_panel,
    canonical_ancestry_truth_bytes,
    canonical_population_provenance_bytes,
    parse_d11_population_provenance,
    parse_mixture_specifications,
    parse_observed_variant_keys,
    transform_population_vcf,
)


_CHECKS = (
    "D11A_D11B_BINDING_MATCHED",
    "PANEL_AND_DOMAIN_MATCHED",
    "HAPLOTYPE_COPY_BIJECTION",
    "ANCESTRY_PROPORTIONS_VALID",
    "ANCESTRY_PROPORTIONS_SUM_TO_ONE",
    "CALIBRATION_VALIDATION_ROLE_SEPARATED",
    "OBSERVED_MARKER_UNIVERSE_PRESERVED",
    "BLINDED_PANEL_HAS_NO_ANCESTRY_LABEL",
)
_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PopulationStructureError("population_structure_json_duplicate_key")
        result[key] = value
    return result


def _load_json(path: Path, code: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise PopulationStructureError(code)
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PopulationStructureError(code) from error
    if not isinstance(value, dict):
        raise PopulationStructureError(code)
    return value


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _write_deterministic_gzip(path: Path, payload: bytes) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            gzip_handle.write(payload)


def _artifact_by_role(
    root: Path,
    manifest: dict[str, Any],
    role: str,
    error_prefix: str,
) -> tuple[dict[str, Any], Path]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise PopulationStructureError(f"{error_prefix}_manifest_invalid")
    matches = [
        artifact
        for artifact in artifacts
        if isinstance(artifact, dict) and artifact.get("role") == role
    ]
    if len(matches) != 1:
        raise PopulationStructureError(f"{error_prefix}_artifact_missing")
    artifact = matches[0]
    relative_path = artifact.get("path")
    if not isinstance(relative_path, str):
        raise PopulationStructureError(f"{error_prefix}_artifact_path_invalid")
    path = root / relative_path
    if (
        not path.is_file()
        or path.is_symlink()
        or root.resolve() not in path.resolve().parents
    ):
        raise PopulationStructureError(f"{error_prefix}_artifact_path_invalid")
    if _sha256_file(path) != artifact.get("sha256"):
        raise PopulationStructureError(f"{error_prefix}_artifact_sha256_mismatch")
    return artifact, path


def _load_inputs(
    d11a_root: Path,
    d11b_root: Path,
) -> tuple[
    dict[str, Any],
    Path,
    dict[str, Any],
    Path,
    Path,
    Path,
    Path,
]:
    for root, code in (
        (d11a_root, "population_structure_d11a_root_invalid"),
        (d11b_root, "population_structure_d11b_root_invalid"),
    ):
        if (
            root.is_symlink()
            or not root.is_dir()
            or root.parent.name != "calibration_workspace"
        ):
            raise PopulationStructureError(code)
    d11a_manifest_path = (
        d11a_root / "private_truth/d11/base_haplotype_manifest.json"
    )
    d11b_manifest_path = (
        d11b_root / "private_truth/d11b/chip_harmonization_manifest.json"
    )
    d11a = _load_json(d11a_manifest_path, "population_structure_d11a_manifest_invalid")
    d11b = _load_json(d11b_manifest_path, "population_structure_d11b_manifest_invalid")
    if (
        d11a.get("manifest_state") != "FINALIZED"
        or d11b.get("manifest_state") != "FINALIZED"
    ):
        raise PopulationStructureError("population_structure_input_not_finalized")
    panel = d11a.get("panel")
    d11a_domain = d11a.get("domain")
    d11b_domain = d11b.get("domain")
    dense_binding = d11b.get("dense_truth")
    if not all(
        isinstance(value, dict)
        for value in (panel, d11a_domain, d11b_domain, dense_binding)
    ):
        raise PopulationStructureError("population_structure_input_contract_invalid")
    d11a_manifest_sha256 = _sha256_file(d11a_manifest_path)
    if (
        dense_binding.get("base_haplotype_manifest_id")
        != d11a.get("base_haplotype_manifest_id")
        or dense_binding.get("base_haplotype_manifest_sha256")
        != d11a_manifest_sha256
        or dense_binding.get("panel_role") != panel.get("role")
        or d11b_domain.get("chromosome") != d11a_domain.get("chromosome")
        or d11b_domain.get("start_bp") != d11a_domain.get("start_bp")
        or d11b_domain.get("end_bp") != d11a_domain.get("end_bp")
    ):
        raise PopulationStructureError("population_structure_d11_binding_mismatch")
    phased_artifact, phased_path = _artifact_by_role(
        d11a_root, d11a, "PHASED_HAPLOTYPES", "population_structure_d11a_phased"
    )
    _, provenance_path = _artifact_by_role(
        d11a_root, d11a, "DONOR_PROVENANCE", "population_structure_d11a_provenance"
    )
    observed_artifact, observed_path = _artifact_by_role(
        d11b_root,
        d11b,
        "CHIP_OBSERVED_GENOTYPES",
        "population_structure_d11b_observed",
    )
    integrity = d11a.get("integrity")
    d11b_integrity = d11b.get("integrity")
    if not isinstance(integrity, dict) or not isinstance(d11b_integrity, dict):
        raise PopulationStructureError("population_structure_input_contract_invalid")
    if (
        phased_artifact.get("record_count") != integrity.get("variant_count")
        or observed_artifact.get("record_count")
        != d11b_integrity.get("observed_vcf_variant_count")
    ):
        raise PopulationStructureError("population_structure_input_count_mismatch")
    return (
        d11a,
        d11a_manifest_path,
        d11b,
        d11b_manifest_path,
        phased_path,
        provenance_path,
        observed_path,
    )


def _observed_sample_ids(path: Path) -> tuple[str, ...]:
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                return tuple(line.rstrip("\n\r").split("\t")[9:])
    raise PopulationStructureError("population_structure_observed_header_missing")


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


def materialize_population_structure(
    *,
    d11a_root: Path,
    d11b_root: Path,
    population_config_path: Path,
    workspace_root: Path,
    run_id: str,
    created_at: str,
) -> dict[str, Any]:
    """Matérialise atomiquement D12 après contrôle des deux couches D11."""
    if _STABLE_IDENTIFIER.fullmatch(run_id) is None:
        raise PopulationStructureError("population_structure_run_identifier_invalid")
    if (
        workspace_root.is_symlink()
        or workspace_root.name != "calibration_workspace"
        or not workspace_root.is_dir()
    ):
        raise PopulationStructureError("population_structure_workspace_invalid")
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise PopulationStructureError("population_structure_run_exists")
    configuration = _load_json(
        population_config_path, "population_structure_config_invalid"
    )
    random_seed, mixture_specifications = parse_mixture_specifications(configuration)
    (
        d11a,
        d11a_manifest_path,
        d11b,
        d11b_manifest_path,
        phased_path,
        provenance_path,
        observed_path,
    ) = _load_inputs(d11a_root, d11b_root)
    panel_contract = d11a["panel"]
    panel_role = panel_contract["role"]
    prefix = {
        "CALIBRATION": "SYN_CAL_D12_",
        "VALIDATION": "SYN_VAL_D12_",
    }.get(panel_role)
    if prefix is None:
        raise PopulationStructureError("population_structure_panel_role_invalid")
    with _open_text(provenance_path) as handle:
        source_diploids = parse_d11_population_provenance(handle)
    if len(source_diploids) != panel_contract.get("individual_count"):
        raise PopulationStructureError("population_structure_source_count_mismatch")
    expected_source_ids = {source.sample_id for source in source_diploids}
    observed_sample_ids = _observed_sample_ids(observed_path)
    expected_prefix = "SYN_CAL_" if panel_role == "CALIBRATION" else "SYN_VAL_"
    if (
        len(observed_sample_ids) != len(set(observed_sample_ids))
        or set(observed_sample_ids) != expected_source_ids
        or any(
            not source_id.startswith(expected_prefix)
            for source_id in expected_source_ids
        )
    ):
        raise PopulationStructureError("population_structure_observed_sample_mismatch")
    population_panel = build_population_panel(
        source_diploids,
        mixture_specifications,
        random_seed=random_seed,
        synthetic_identifier_prefix=prefix,
    )
    with _open_text(observed_path) as handle:
        observed_variant_keys = parse_observed_variant_keys(handle)

    temporary_root = Path(
        tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root)
    )
    try:
        private_root = temporary_root / "private_truth/d12"
        private_root.mkdir(parents=True, mode=0o700)
        ancestry_path = private_root / "ancestry_truth.tsv.gz"
        provenance_output_path = private_root / "d11_haplotype_provenance.tsv.gz"
        _write_deterministic_gzip(
            ancestry_path, canonical_ancestry_truth_bytes(population_panel)
        )
        _write_deterministic_gzip(
            provenance_output_path,
            canonical_population_provenance_bytes(population_panel),
        )

        dense_buffer = temporary_root / "population_panel.phased.vcf"
        observed_buffer = temporary_root / "population_panel.observed.vcf"
        domain = d11a["domain"]
        with (
            _open_text(phased_path) as source_handle,
            dense_buffer.open("w", encoding="utf-8", newline="") as dense_handle,
            observed_buffer.open("w", encoding="utf-8", newline="") as observed_handle,
        ):
            stats = transform_population_vcf(
                source_handle,
                dense_handle,
                observed_handle,
                panel=population_panel,
                observed_variant_keys=observed_variant_keys,
                chromosome=domain["chromosome"],
            )
        dense_output_path = private_root / "population_panel.phased.vcf.gz"
        observed_output_path = private_root / "blinded_observed_panel.vcf.gz"
        _write_deterministic_gzip(dense_output_path, dense_buffer.read_bytes())
        _write_deterministic_gzip(observed_output_path, observed_buffer.read_bytes())
        dense_buffer.unlink()
        observed_buffer.unlink()
        for path in (
            ancestry_path,
            provenance_output_path,
            dense_output_path,
            observed_output_path,
        ):
            path.chmod(0o600)

        pure_count = sum(
            len(diploid.mixture_proportions) == 1 for diploid in population_panel
        )
        admixed_count = len(population_panel) - pure_count
        manifest = {
            "schema_version": "1.0.0",
            "population_structure_manifest_id": f"{run_id}.d12",
            "manifest_state": "FINALIZED",
            "created_at": created_at,
            "inputs": {
                "d11a": {
                    "manifest_id": d11a["base_haplotype_manifest_id"],
                    "manifest_sha256": _sha256_file(d11a_manifest_path),
                },
                "d11b": {
                    "manifest_id": d11b["chip_harmonization_manifest_id"],
                    "manifest_sha256": _sha256_file(d11b_manifest_path),
                },
                "panel_role": panel_role,
            },
            "configuration": {
                "configuration_id": configuration[
                    "population_structure_config_id"
                ],
                "configuration_sha256": _sha256_file(population_config_path),
                "random_seed": random_seed,
            },
            "panel": {
                "panel_id": f"{run_id}.{panel_role.lower()}",
                "role": panel_role,
                "synthetic_identifier_prefix": prefix,
                "individual_count": len(population_panel),
                "haplotype_count": 2 * len(population_panel),
                "pure_individual_count": pure_count,
                "admixed_individual_count": admixed_count,
                "source_stratum_count": len(
                    {source.stratum_id for source in source_diploids}
                ),
                "output_stratum_count": len(
                    {diploid.stratum_id for diploid in population_panel}
                ),
                "mixture_profile_count": len(mixture_specifications),
            },
            "domain": {
                "chromosome": domain["chromosome"],
                "start_bp": domain["start_bp"],
                "end_bp": domain["end_bp"],
                "position_unit": "BP_1_BASED_CLOSED",
            },
            "method": {
                "method_id": POPULATION_STRUCTURE_METHOD_ID,
                "method_version": POPULATION_STRUCTURE_METHOD_VERSION,
                "mixture_unit": "RECIPROCAL_WHOLE_HAPLOTYPE_PAIR",
                "supported_proportions": "PURE_OR_TWO_WAY_EQUAL",
                "unequal_mixtures": "DEFERRED_TO_D14_RECOMBINATION_MOSAICS",
            },
            "confidentiality": {
                "real_sample_files_read": False,
                "raw_donor_identifiers_in_d12_artifacts": False,
                "ancestry_labels_in_blinded_panel": False,
                "local_source_paths_in_manifest": False,
                "artifacts_version_control_allowed": False,
            },
            "artifacts": [
                _artifact(
                    run_id,
                    "dense_phased",
                    "DENSE_PHASED_POPULATION_PANEL",
                    "private_truth/d12/population_panel.phased.vcf.gz",
                    dense_output_path,
                    stats.dense_variant_count,
                ),
                _artifact(
                    run_id,
                    "ancestry_truth",
                    "ANCESTRY_TRUTH",
                    "private_truth/d12/ancestry_truth.tsv.gz",
                    ancestry_path,
                    len(population_panel),
                ),
                _artifact(
                    run_id,
                    "d11_provenance",
                    "D11_HAPLOTYPE_PROVENANCE",
                    "private_truth/d12/d11_haplotype_provenance.tsv.gz",
                    provenance_output_path,
                    2 * len(population_panel),
                ),
                _artifact(
                    run_id,
                    "blinded_observed",
                    "BLINDED_OBSERVED_POPULATION_PANEL",
                    "private_truth/d12/blinded_observed_panel.vcf.gz",
                    observed_output_path,
                    stats.observed_variant_count,
                ),
            ],
            "integrity": {
                "dense_variant_count": stats.dense_variant_count,
                "observed_variant_count": stats.observed_variant_count,
                "individual_count": stats.individual_count,
                "source_haplotype_count": stats.source_haplotype_count,
                "output_haplotype_count": stats.output_haplotype_count,
                "haplotype_copy_error_count": stats.haplotype_copy_error_count,
                "allele_count_difference": stats.allele_count_difference,
                "genotype_haplotype_discordance_count": (
                    stats.genotype_haplotype_discordance_count
                ),
                "invalid_proportion_count": 0,
                "proportion_sum_error_count": 0,
                "panel_role_collision_count": 0,
                "checks": list(_CHECKS),
            },
        }
        manifest_path = private_root / "population_structure_manifest.json"
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
        description="Matérialise les strates et mélanges D12 reproductibles."
    )
    parser.add_argument("--d11a-root", type=Path, required=True)
    parser.add_argument("--d11b-root", type=Path, required=True)
    parser.add_argument("--population-config", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    return parser


def main() -> int:
    """Point d'entrée D12 ; n'affiche que des comptes agrégés."""
    arguments = _parser().parse_args()
    manifest = materialize_population_structure(
        d11a_root=arguments.d11a_root,
        d11b_root=arguments.d11b_root,
        population_config_path=arguments.population_config,
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
