"""Matérialisation privée D13 des pedigrees et génotypes descendants."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, TextIO

from effet_fondateur.calibration.simulation.pedigrees import (
    PEDIGREE_METHOD_ID,
    PEDIGREE_METHOD_VERSION,
    PedigreeSimulationError,
    build_pedigree_panel,
    canonical_ancestry_bytes,
    canonical_kinship_bytes,
    canonical_pedigree_bytes,
    canonical_provenance_bytes,
    canonical_transmission_bytes,
    parse_d12_source_truth,
    parse_pedigree_configuration,
    transform_pedigree_vcf,
)
from effet_fondateur.calibration.simulation.population_structure import (
    parse_observed_variant_keys,
)


_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_CHECKS = (
    "D12_BINDING_MATCHED",
    "SOURCE_SAMPLE_UNIVERSE_MATCHED",
    "FOUNDER_HAPLOTYPE_BIJECTION",
    "FAMILY_FOUNDERS_DISJOINT",
    "MENDELIAN_TRANSMISSIONS_EXACT",
    "TRIO_COUNT_MATCHED",
    "KINSHIP_CLASSES_COVERED",
    "ANCESTRY_PROPAGATED",
    "BLINDED_PANEL_HAS_NO_TRUTH_LABEL",
)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PedigreeSimulationError("pedigree_json_duplicate_key")
        result[key] = value
    return result


def _load_json(path: Path, code: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise PedigreeSimulationError(code)
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PedigreeSimulationError(code) from error
    if not isinstance(value, dict):
        raise PedigreeSimulationError(code)
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _write_gzip(path: Path, payload: bytes) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_handle, mtime=0
        ) as gzip_handle:
            gzip_handle.write(payload)


def _artifact_by_role(
    root: Path,
    manifest: dict[str, Any],
    role: str,
) -> tuple[dict[str, Any], Path]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise PedigreeSimulationError("pedigree_d12_manifest_invalid")
    matches = [
        artifact
        for artifact in artifacts
        if isinstance(artifact, dict) and artifact.get("role") == role
    ]
    if len(matches) != 1:
        raise PedigreeSimulationError("pedigree_d12_artifact_missing", role)
    artifact = matches[0]
    relative_path = artifact.get("path")
    if not isinstance(relative_path, str):
        raise PedigreeSimulationError("pedigree_d12_artifact_path_invalid")
    path = root / relative_path
    if (
        not path.is_file()
        or path.is_symlink()
        or root.resolve() not in path.resolve().parents
    ):
        raise PedigreeSimulationError("pedigree_d12_artifact_path_invalid")
    if _sha256(path) != artifact.get("sha256"):
        raise PedigreeSimulationError("pedigree_d12_artifact_sha256_mismatch")
    return artifact, path


def _load_d12(
    d12_root: Path,
) -> tuple[
    dict[str, Any],
    Path,
    dict[str, Any],
    Path,
    dict[str, Any],
    Path,
    dict[str, Any],
    Path,
    dict[str, Any],
    Path,
]:
    if (
        d12_root.is_symlink()
        or not d12_root.is_dir()
        or d12_root.parent.name != "calibration_workspace"
    ):
        raise PedigreeSimulationError("pedigree_d12_root_invalid")
    manifest_path = d12_root / "private_truth/d12/population_structure_manifest.json"
    manifest = _load_json(manifest_path, "pedigree_d12_manifest_invalid")
    if manifest.get("manifest_state") != "FINALIZED":
        raise PedigreeSimulationError("pedigree_d12_manifest_not_finalized")
    panel = manifest.get("panel")
    domain = manifest.get("domain")
    if not isinstance(panel, dict) or not isinstance(domain, dict):
        raise PedigreeSimulationError("pedigree_d12_manifest_invalid")
    dense_artifact, dense_path = _artifact_by_role(
        d12_root, manifest, "DENSE_PHASED_POPULATION_PANEL"
    )
    observed_artifact, observed_path = _artifact_by_role(
        d12_root, manifest, "BLINDED_OBSERVED_POPULATION_PANEL"
    )
    ancestry_artifact, ancestry_path = _artifact_by_role(
        d12_root, manifest, "ANCESTRY_TRUTH"
    )
    provenance_artifact, provenance_path = _artifact_by_role(
        d12_root, manifest, "D11_HAPLOTYPE_PROVENANCE"
    )
    if (
        dense_artifact.get("record_count")
        != manifest.get("integrity", {}).get("dense_variant_count")
        or observed_artifact.get("record_count")
        != manifest.get("integrity", {}).get("observed_variant_count")
        or ancestry_artifact.get("record_count") != panel.get("individual_count")
        or provenance_artifact.get("record_count")
        != 2 * panel.get("individual_count", 0)
    ):
        raise PedigreeSimulationError("pedigree_d12_artifact_count_mismatch")
    return (
        manifest,
        manifest_path,
        dense_artifact,
        dense_path,
        observed_artifact,
        observed_path,
        ancestry_artifact,
        ancestry_path,
        provenance_artifact,
        provenance_path,
    )


def _vcf_samples(path: Path) -> tuple[str, ...]:
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                return tuple(line.rstrip("\n\r").split("\t")[9:])
    raise PedigreeSimulationError("pedigree_vcf_header_missing")


def _artifact(
    run_id: str,
    suffix: str,
    role: str,
    path: Path,
    relative_path: str,
    record_count: int,
) -> dict[str, Any]:
    return {
        "artifact_id": f"{run_id}.{suffix}",
        "role": role,
        "path": relative_path,
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
        "record_count": record_count,
    }


def materialize_pedigrees(
    *,
    d12_root: Path,
    pedigree_config_path: Path,
    workspace_root: Path,
    run_id: str,
    created_at: str,
) -> dict[str, Any]:
    """Matérialise atomiquement D13 à partir d'un unique panel D12 finalisé."""
    if _STABLE_IDENTIFIER.fullmatch(run_id) is None:
        raise PedigreeSimulationError("pedigree_run_identifier_invalid")
    if (
        workspace_root.is_symlink()
        or workspace_root.name != "calibration_workspace"
        or not workspace_root.is_dir()
    ):
        raise PedigreeSimulationError("pedigree_workspace_invalid")
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise PedigreeSimulationError("pedigree_run_exists")
    configuration = _load_json(pedigree_config_path, "pedigree_config_invalid")
    random_seed, family_specs, unrelated_pair_count = (
        parse_pedigree_configuration(configuration)
    )
    (
        d12_manifest,
        d12_manifest_path,
        dense_input_artifact,
        dense_input_path,
        observed_input_artifact,
        observed_input_path,
        ancestry_input_artifact,
        ancestry_input_path,
        provenance_input_artifact,
        provenance_input_path,
    ) = _load_d12(d12_root)
    panel_contract = d12_manifest["panel"]
    panel_role = panel_contract["role"]
    output_prefix = {
        "CALIBRATION": "SYN_CAL_D13_",
        "VALIDATION": "SYN_VAL_D13_",
    }.get(panel_role)
    source_prefix = {
        "CALIBRATION": "SYN_CAL_D12_",
        "VALIDATION": "SYN_VAL_D12_",
    }.get(panel_role)
    if output_prefix is None or source_prefix is None:
        raise PedigreeSimulationError("pedigree_panel_role_invalid")
    with (
        _open_text(ancestry_input_path) as ancestry_handle,
        _open_text(provenance_input_path) as provenance_handle,
    ):
        source_individuals = parse_d12_source_truth(
            ancestry_handle, provenance_handle
        )
    source_ids = {source.sample_id for source in source_individuals}
    observed_samples = _vcf_samples(observed_input_path)
    dense_samples = _vcf_samples(dense_input_path)
    if (
        len(source_individuals) != panel_contract.get("individual_count")
        or len(observed_samples) != len(set(observed_samples))
        or len(dense_samples) != len(set(dense_samples))
        or set(observed_samples) != source_ids
        or set(dense_samples) != source_ids
        or any(not value.startswith(source_prefix) for value in source_ids)
    ):
        raise PedigreeSimulationError("pedigree_source_sample_universe_mismatch")
    panel = build_pedigree_panel(
        source_individuals,
        family_specs,
        unrelated_pair_count=unrelated_pair_count,
        random_seed=random_seed,
        synthetic_identifier_prefix=output_prefix,
    )
    with _open_text(observed_input_path) as handle:
        observed_keys = parse_observed_variant_keys(handle)

    temporary_root = Path(
        tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root)
    )
    try:
        private_root = temporary_root / "private_truth/d13"
        private_root.mkdir(parents=True, mode=0o700)
        payloads = {
            "pedigree_truth.tsv.gz": canonical_pedigree_bytes(panel),
            "kinship_truth.tsv.gz": canonical_kinship_bytes(panel),
            "mendelian_transmission_truth.tsv.gz": canonical_transmission_bytes(
                panel
            ),
            "d12_haplotype_provenance.tsv.gz": canonical_provenance_bytes(panel),
            "ancestry_truth.tsv.gz": canonical_ancestry_bytes(panel),
        }
        for filename, payload in payloads.items():
            _write_gzip(private_root / filename, payload)

        dense_buffer = temporary_root / "pedigree_panel.phased.vcf"
        observed_buffer = temporary_root / "pedigree_panel.observed.vcf"
        domain = d12_manifest["domain"]
        with (
            _open_text(dense_input_path) as source_handle,
            dense_buffer.open("w", encoding="utf-8", newline="") as dense_handle,
            observed_buffer.open(
                "w", encoding="utf-8", newline=""
            ) as observed_handle,
        ):
            stats = transform_pedigree_vcf(
                source_handle,
                dense_handle,
                observed_handle,
                panel=panel,
                observed_variant_keys=observed_keys,
                chromosome=domain["chromosome"],
            )
        dense_output_path = private_root / "pedigree_panel.phased.vcf.gz"
        observed_output_path = private_root / "blinded_observed_panel.vcf.gz"
        _write_gzip(dense_output_path, dense_buffer.read_bytes())
        _write_gzip(observed_output_path, observed_buffer.read_bytes())
        dense_buffer.unlink()
        observed_buffer.unlink()
        for path in private_root.iterdir():
            path.chmod(0o600)

        relationship_counts = Counter(
            record.relationship_class for record in panel.kinship_records
        )
        descendants = len(panel.individuals) - panel.source_individual_count
        manifest = {
            "schema_version": "1.0.0",
            "pedigree_manifest_id": f"{run_id}.d13",
            "manifest_state": "FINALIZED",
            "created_at": created_at,
            "input": {
                "manifest_id": d12_manifest["population_structure_manifest_id"],
                "manifest_sha256": _sha256(d12_manifest_path),
                "panel_role": panel_role,
                "dense_artifact_id": dense_input_artifact["artifact_id"],
                "dense_artifact_sha256": dense_input_artifact["sha256"],
                "observed_artifact_id": observed_input_artifact["artifact_id"],
                "observed_artifact_sha256": observed_input_artifact["sha256"],
                "ancestry_artifact_id": ancestry_input_artifact["artifact_id"],
                "ancestry_artifact_sha256": ancestry_input_artifact["sha256"],
                "provenance_artifact_id": provenance_input_artifact[
                    "artifact_id"
                ],
                "provenance_artifact_sha256": provenance_input_artifact[
                    "sha256"
                ],
            },
            "configuration": {
                "configuration_id": configuration["pedigree_config_id"],
                "configuration_sha256": _sha256(pedigree_config_path),
                "random_seed": random_seed,
            },
            "panel": {
                "panel_id": f"{run_id}.{panel_role.lower()}",
                "role": panel_role,
                "synthetic_identifier_prefix": output_prefix,
                "source_individual_count": panel.source_individual_count,
                "output_individual_count": len(panel.individuals),
                "descendant_count": descendants,
                "family_count": panel.family_count,
                "trio_count": panel.trio_count,
                "unrelated_pair_count": panel.unrelated_pair_count,
            },
            "domain": {
                "chromosome": domain["chromosome"],
                "start_bp": domain["start_bp"],
                "end_bp": domain["end_bp"],
                "position_unit": "BP_1_BASED_CLOSED",
            },
            "method": {
                "method_id": PEDIGREE_METHOD_ID,
                "method_version": PEDIGREE_METHOD_VERSION,
                "family_template": "FIRST_COUSIN_PEDIGREE",
                "source_selection": "SAME_STRATUM_WITHOUT_REPLACEMENT",
                "transmission_unit": "WHOLE_PARENTAL_HAPLOTYPE",
                "kinship_definition": "PEDIGREE_COANCESTRY_COEFFICIENT",
                "recombination_policy": "FORBIDDEN_DEFERRED_TO_D14",
            },
            "confidentiality": {
                "real_sample_files_read": False,
                "d12_identifiers_in_blinded_panel": False,
                "pedigree_labels_in_blinded_panel": False,
                "ancestry_labels_in_blinded_panel": False,
                "local_source_paths_in_manifest": False,
                "artifacts_version_control_allowed": False,
            },
            "artifacts": [
                _artifact(
                    run_id,
                    "dense_phased",
                    "DENSE_PHASED_PEDIGREE_PANEL",
                    dense_output_path,
                    "private_truth/d13/pedigree_panel.phased.vcf.gz",
                    stats.dense_variant_count,
                ),
                _artifact(
                    run_id,
                    "blinded_observed",
                    "BLINDED_OBSERVED_PEDIGREE_PANEL",
                    observed_output_path,
                    "private_truth/d13/blinded_observed_panel.vcf.gz",
                    stats.observed_variant_count,
                ),
                _artifact(
                    run_id,
                    "pedigree_truth",
                    "PEDIGREE_TRUTH",
                    private_root / "pedigree_truth.tsv.gz",
                    "private_truth/d13/pedigree_truth.tsv.gz",
                    len(panel.individuals),
                ),
                _artifact(
                    run_id,
                    "kinship_truth",
                    "KINSHIP_TRUTH",
                    private_root / "kinship_truth.tsv.gz",
                    "private_truth/d13/kinship_truth.tsv.gz",
                    len(panel.kinship_records),
                ),
                _artifact(
                    run_id,
                    "transmission_truth",
                    "MENDELIAN_TRANSMISSION_TRUTH",
                    private_root / "mendelian_transmission_truth.tsv.gz",
                    "private_truth/d13/mendelian_transmission_truth.tsv.gz",
                    len(panel.transmissions),
                ),
                _artifact(
                    run_id,
                    "d12_provenance",
                    "D12_HAPLOTYPE_PROVENANCE",
                    private_root / "d12_haplotype_provenance.tsv.gz",
                    "private_truth/d13/d12_haplotype_provenance.tsv.gz",
                    2 * len(panel.individuals),
                ),
                _artifact(
                    run_id,
                    "ancestry_truth",
                    "ANCESTRY_TRUTH",
                    private_root / "ancestry_truth.tsv.gz",
                    "private_truth/d13/ancestry_truth.tsv.gz",
                    len(panel.individuals),
                ),
            ],
            "integrity": {
                "dense_variant_count": stats.dense_variant_count,
                "observed_variant_count": stats.observed_variant_count,
                "transmission_record_count": len(panel.transmissions),
                "kinship_record_count": len(panel.kinship_records),
                "relationship_class_counts": {
                    key: relationship_counts[key]
                    for key in (
                        "UNRELATED",
                        "FIRST_DEGREE",
                        "SECOND_DEGREE",
                        "THIRD_DEGREE",
                    )
                },
                "founder_haplotype_copy_error_count": (
                    stats.founder_haplotype_copy_error_count
                ),
                "founder_allele_count_difference": (
                    stats.founder_allele_count_difference
                ),
                "mendelian_error_count": stats.mendelian_error_count,
                "family_founder_overlap_count": 0,
                "invalid_parent_reference_count": 0,
                "checks": list(_CHECKS),
            },
        }
        manifest_path = private_root / "pedigree_manifest.json"
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
        description="Matérialise les pedigrees et transmissions D13."
    )
    parser.add_argument("--d12-root", type=Path, required=True)
    parser.add_argument("--pedigree-config", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    return parser


def main() -> int:
    """Point d'entrée D13 ; n'affiche que des compteurs agrégés."""
    arguments = _parser().parse_args()
    manifest = materialize_pedigrees(
        d12_root=arguments.d12_root,
        pedigree_config_path=arguments.pedigree_config,
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
