"""Matérialisation contrôlée des partitions D10 dans le workspace ignoré."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.calibration.partitioning import (
    CalibrationPartitionError,
    PartitionAudit,
    build_partition_manifest,
    canonical_membership_bytes,
    parse_population_pedigree,
    parse_unrelated_sequence_index,
    partition_reference_identifiers,
    validate_partition_manifest,
)


_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise CalibrationPartitionError("partition_catalog_duplicate_key")
        result[key] = value
    return result


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _catalog_source(catalog: dict[str, Any], source_id: str) -> dict[str, Any]:
    sources = catalog.get("sources")
    if not isinstance(sources, list):
        raise CalibrationPartitionError("partition_catalog_invalid")
    matches = [source for source in sources if source.get("source_id") == source_id]
    if len(matches) != 1 or not isinstance(matches[0], dict):
        raise CalibrationPartitionError("partition_catalog_source_missing")
    return matches[0]


def _verify_input_hash(
    path: Path,
    catalog: dict[str, Any],
    source_id: str,
) -> None:
    source = _catalog_source(catalog, source_id)
    fingerprint = source.get("fingerprint")
    if not isinstance(fingerprint, dict):
        raise CalibrationPartitionError("partition_catalog_source_invalid")
    if not path.is_file() or path.is_symlink():
        raise CalibrationPartitionError("partition_input_file_invalid")
    if _sha256_file(path) != fingerprint.get("sha256"):
        raise CalibrationPartitionError("partition_input_sha256_mismatch")


def materialize_reference_partitions(
    *,
    public_catalog_path: Path,
    population_metadata_path: Path,
    unrelated_index_path: Path,
    workspace_root: Path,
    partition_run_id: str,
    created_at: str,
) -> PartitionAudit:
    """Produit atomiquement les listes privées et leur manifest agrégé.

    Le dossier cible doit être un répertoire ``calibration_workspace``. Un run
    existant n'est jamais remplacé. Les fichiers de membres sont créés avec des
    permissions utilisateur uniquement et aucun identifiant n'est envoyé sur la
    sortie standard.
    """
    if _STABLE_IDENTIFIER.fullmatch(partition_run_id) is None:
        raise CalibrationPartitionError("partition_run_id_invalid")
    if workspace_root.name != "calibration_workspace" or workspace_root.is_symlink():
        raise CalibrationPartitionError("partition_workspace_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / partition_run_id
    if target_root.exists() or target_root.is_symlink():
        raise CalibrationPartitionError("partition_run_exists")
    if not public_catalog_path.is_file() or public_catalog_path.is_symlink():
        raise CalibrationPartitionError("partition_catalog_file_invalid")

    catalog = json.loads(
        public_catalog_path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_keys,
    )
    catalog_sha256 = _sha256_file(public_catalog_path)
    _verify_input_hash(
        population_metadata_path,
        catalog,
        "igsr_1kg_3202_population_metadata_20130606",
    )
    _verify_input_hash(
        unrelated_index_path,
        catalog,
        "igsr_1kg_2504_unrelated_sequence_index",
    )
    with population_metadata_path.open(encoding="utf-8") as file_handle:
        source_identifiers, independence_groups = parse_population_pedigree(
            file_handle
        )
    with unrelated_index_path.open(encoding="utf-8") as file_handle:
        candidates = parse_unrelated_sequence_index(file_handle)
    panel = _catalog_source(
        catalog,
        "igsr_1kg_3202_phased_grch38_20220422",
    )
    panel_fingerprint = panel.get("fingerprint")
    if not isinstance(panel_fingerprint, dict):
        raise CalibrationPartitionError("partition_catalog_source_invalid")
    partitions = partition_reference_identifiers(
        source_identifiers=source_identifiers,
        phasing_reference_candidate_identifiers=candidates,
        independence_group_by_identifier=independence_groups,
        source_fingerprint_sha256=str(panel_fingerprint.get("sha256", "")),
    )
    manifest, memberships = build_partition_manifest(
        partition_manifest_id=f"{partition_run_id}.manifest",
        created_at=created_at,
        partitions=partitions,
        source_identifiers=source_identifiers,
        independence_group_by_identifier=independence_groups,
        public_reference_catalog=catalog,
        public_reference_catalog_sha256=catalog_sha256,
    )

    temporary_root = Path(
        tempfile.mkdtemp(prefix=f".{partition_run_id}.", dir=workspace_root)
    )
    try:
        private_root = temporary_root / "private_truth"
        membership_root = private_root / "partitions"
        membership_root.mkdir(parents=True, mode=0o700)
        for partition_id, identifiers in memberships.items():
            membership_path = membership_root / f"{partition_id}.ids.txt"
            membership_path.write_bytes(canonical_membership_bytes(identifiers))
            membership_path.chmod(0o600)
        manifest_path = private_root / "reference_partition_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        manifest_path.chmod(0o600)
        audit = validate_partition_manifest(
            manifest,
            memberships_by_partition_id=memberships,
            source_identifiers=source_identifiers,
            independence_group_by_identifier=independence_groups,
            public_reference_catalog=catalog,
            public_reference_catalog_sha256=catalog_sha256,
        )
        temporary_root.rename(target_root)
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    return audit


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Matérialise les partitions privées D10 sans afficher d'identifiant."
    )
    parser.add_argument("--public-catalog", type=Path, required=True)
    parser.add_argument("--population-metadata", type=Path, required=True)
    parser.add_argument("--unrelated-index", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--partition-run-id", required=True)
    return parser


def main() -> int:
    """Point d'entrée CLI D10 ; imprime uniquement un résumé agrégé."""
    arguments = _parser().parse_args()
    created_at = datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )
    audit = materialize_reference_partitions(
        public_catalog_path=arguments.public_catalog,
        population_metadata_path=arguments.population_metadata,
        unrelated_index_path=arguments.unrelated_index,
        workspace_root=arguments.workspace_root,
        partition_run_id=arguments.partition_run_id,
        created_at=created_at,
    )
    print(
        json.dumps(
            {
                "partition_manifest_id": audit.partition_manifest_id,
                "source_identifier_count": audit.source_identifier_count,
                "donor_calibration_count": audit.donor_calibration_count,
                "donor_validation_count": audit.donor_validation_count,
                "phasing_reference_count": audit.phasing_reference_count,
                "independence_group_count": audit.independence_group_count,
                "checks": list(audit.checks),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
