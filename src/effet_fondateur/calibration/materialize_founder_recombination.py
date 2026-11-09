"""Matérialisation privée D14 du fondateur et des recombinaisons."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, TextIO

from effet_fondateur.calibration.simulation.founder_recombination import (
    FOUNDER_RECOMBINATION_METHOD_ID,
    FOUNDER_RECOMBINATION_METHOD_VERSION,
    FounderRecombinationError,
    build_founder_plan,
    canonical_carrier_bytes,
    canonical_founder_segments_bytes,
    canonical_local_ancestry_bytes,
    canonical_target_truth_bytes,
    interpolate_cm,
    parse_d12_haplotype_strata,
    parse_d13_origins,
    parse_founder_configuration,
    parse_genetic_map,
    parse_pedigree_truth,
    transform_founder_vcf,
)
from effet_fondateur.calibration.simulation.population_structure import (
    parse_observed_variant_keys,
)


_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_CHECKS = (
    "D13_D12_D11_LINEAGE_MATCHED",
    "GENETIC_MAP_HASH_MATCHED",
    "TARGET_VARIANT_MATCHED",
    "FAMILIES_DISJOINT",
    "FOUNDER_TEMPLATE_EXTERNAL_TO_FAMILIES",
    "BREAKPOINTS_WITHIN_DOMAIN",
    "CM_BP_COORDINATES_COHERENT",
    "MENDELIAN_TRANSMISSIONS_EXACT",
    "BLINDED_PANEL_HAS_NO_TRUTH_LABEL",
)


def _duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise FounderRecombinationError("founder_json_duplicate_key")
        result[key] = value
    return result


def _load_json(path: Path, code: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise FounderRecombinationError(code)
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise FounderRecombinationError(code) from error
    if not isinstance(value, dict):
        raise FounderRecombinationError(code)
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


def _validate_root(root: Path, stage: str) -> None:
    if (
        root.is_symlink()
        or not root.is_dir()
        or root.parent.name != "calibration_workspace"
    ):
        raise FounderRecombinationError(f"founder_{stage}_root_invalid")


def _artifact(
    root: Path,
    manifest: dict[str, Any],
    role: str,
) -> tuple[dict[str, Any], Path]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise FounderRecombinationError("founder_input_manifest_invalid")
    matches = [
        item
        for item in artifacts
        if isinstance(item, dict) and item.get("role") == role
    ]
    if len(matches) != 1:
        raise FounderRecombinationError("founder_input_artifact_missing", role)
    item = matches[0]
    relative = item.get("path")
    if not isinstance(relative, str):
        raise FounderRecombinationError("founder_input_artifact_path_invalid")
    path = root / relative
    if (
        not path.is_file()
        or path.is_symlink()
        or root.resolve() not in path.resolve().parents
    ):
        raise FounderRecombinationError("founder_input_artifact_path_invalid")
    if _sha256(path) != item.get("sha256"):
        raise FounderRecombinationError("founder_input_artifact_sha256_mismatch")
    return item, path


def _transmissions(path: Path) -> tuple[tuple[str, str, int, int], ...]:
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "family_unit_id",
            "child_id",
            "parent_id",
            "parental_haplotype_index",
            "child_haplotype_index",
            "recombination_count",
        ):
            raise FounderRecombinationError(
                "founder_transmission_columns_invalid"
            )
        records: list[tuple[str, str, int, int]] = []
        for row in reader:
            try:
                parent_index = int(row["parental_haplotype_index"])
                child_index = int(row["child_haplotype_index"])
                recombination_count = int(row["recombination_count"])
            except ValueError as error:
                raise FounderRecombinationError(
                    "founder_transmission_value_invalid"
                ) from error
            if (
                parent_index not in (1, 2)
                or child_index not in (1, 2)
                or recombination_count != 0
            ):
                raise FounderRecombinationError(
                    "founder_transmission_value_invalid"
                )
            records.append(
                (
                    row["child_id"],
                    row["parent_id"],
                    parent_index,
                    child_index,
                )
            )
    return tuple(records)


def _output_artifact(
    run_id: str,
    suffix: str,
    role: str,
    path: Path,
    record_count: int,
) -> dict[str, Any]:
    return {
        "artifact_id": f"{run_id}.{suffix}",
        "role": role,
        "path": f"private_truth/d14/{path.name}",
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
        "record_count": record_count,
    }


def materialize_founder_recombination(
    *,
    d13_root: Path,
    d12_root: Path,
    d11a_root: Path,
    founder_config_path: Path,
    workspace_root: Path,
    run_id: str,
    created_at: str,
) -> dict[str, Any]:
    """Matérialise atomiquement D14 et lie toute la chaîne D11–D13."""
    if _IDENTIFIER.fullmatch(run_id) is None:
        raise FounderRecombinationError("founder_run_identifier_invalid")
    if (
        workspace_root.is_symlink()
        or not workspace_root.is_dir()
        or workspace_root.name != "calibration_workspace"
    ):
        raise FounderRecombinationError("founder_workspace_invalid")
    for root, stage in ((d13_root, "d13"), (d12_root, "d12"), (d11a_root, "d11a")):
        _validate_root(root, stage)
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise FounderRecombinationError("founder_run_exists")
    configuration_document = _load_json(
        founder_config_path, "founder_config_invalid"
    )
    configuration = parse_founder_configuration(configuration_document)
    d13_manifest_path = d13_root / "private_truth/d13/pedigree_manifest.json"
    d12_manifest_path = (
        d12_root / "private_truth/d12/population_structure_manifest.json"
    )
    d11_manifest_path = d11a_root / "private_truth/d11/base_haplotype_manifest.json"
    d13_manifest = _load_json(d13_manifest_path, "founder_d13_manifest_invalid")
    d12_manifest = _load_json(d12_manifest_path, "founder_d12_manifest_invalid")
    d11_manifest = _load_json(d11_manifest_path, "founder_d11_manifest_invalid")
    if any(
        manifest.get("manifest_state") != "FINALIZED"
        for manifest in (d13_manifest, d12_manifest, d11_manifest)
    ):
        raise FounderRecombinationError("founder_input_manifest_not_finalized")
    if (
        d13_manifest.get("input", {}).get("manifest_id")
        != d12_manifest.get("population_structure_manifest_id")
        or d13_manifest.get("input", {}).get("manifest_sha256")
        != _sha256(d12_manifest_path)
        or d12_manifest.get("inputs", {}).get("d11a", {}).get("manifest_id")
        != d11_manifest.get("base_haplotype_manifest_id")
        or d12_manifest.get("inputs", {}).get("d11a", {}).get("manifest_sha256")
        != _sha256(d11_manifest_path)
    ):
        raise FounderRecombinationError("founder_input_lineage_mismatch")
    domains = (
        d13_manifest.get("domain"),
        d12_manifest.get("domain"),
        d11_manifest.get("domain"),
    )
    if not all(isinstance(domain, dict) for domain in domains):
        raise FounderRecombinationError("founder_input_domain_invalid")
    domain_keys = ("chromosome", "start_bp", "end_bp")
    if any(
        tuple(domain.get(key) for key in domain_keys)
        != tuple(domains[0].get(key) for key in domain_keys)
        for domain in domains[1:]
    ):
        raise FounderRecombinationError("founder_input_domain_mismatch")
    if configuration.chromosome != domains[0]["chromosome"]:
        raise FounderRecombinationError("founder_target_chromosome_mismatch")

    dense_item, dense_path = _artifact(
        d13_root, d13_manifest, "DENSE_PHASED_PEDIGREE_PANEL"
    )
    observed_item, observed_path = _artifact(
        d13_root, d13_manifest, "BLINDED_OBSERVED_PEDIGREE_PANEL"
    )
    pedigree_item, pedigree_path = _artifact(
        d13_root, d13_manifest, "PEDIGREE_TRUTH"
    )
    transmission_item, transmission_path = _artifact(
        d13_root, d13_manifest, "MENDELIAN_TRANSMISSION_TRUTH"
    )
    provenance_item, provenance_path = _artifact(
        d13_root, d13_manifest, "D12_HAPLOTYPE_PROVENANCE"
    )
    d12_provenance_item, d12_provenance_path = _artifact(
        d12_root, d12_manifest, "D11_HAPLOTYPE_PROVENANCE"
    )
    map_item, map_path = _artifact(d11a_root, d11_manifest, "GENETIC_MAP")
    with (
        _open_text(pedigree_path) as pedigree_handle,
        _open_text(provenance_path) as provenance_handle,
        _open_text(d12_provenance_path) as d12_provenance_handle,
        _open_text(map_path) as map_handle,
    ):
        pedigree = parse_pedigree_truth(pedigree_handle)
        origins = parse_d13_origins(provenance_handle)
        strata = parse_d12_haplotype_strata(d12_provenance_handle)
        genetic_map = parse_genetic_map(map_handle, configuration.chromosome)
    plan = build_founder_plan(
        configuration=configuration,
        genetic_map=genetic_map,
        pedigree=pedigree,
        sample_origins=origins,
        background_strata=strata,
    )
    with _open_text(observed_path) as handle:
        observed_keys = parse_observed_variant_keys(handle)
    transmission_records = _transmissions(transmission_path)
    temporary_root = Path(
        tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root)
    )
    try:
        private = temporary_root / "private_truth/d14"
        private.mkdir(parents=True, mode=0o700)
        generated_payloads = {
            "founder_segments.tsv.gz": canonical_founder_segments_bytes(plan),
            "founder_carriers.tsv.gz": canonical_carrier_bytes(plan),
            "target_variant_truth.tsv.gz": canonical_target_truth_bytes(plan),
            "local_ancestry_mosaic_truth.tsv.gz": (
                canonical_local_ancestry_bytes(plan)
            ),
        }
        for filename, payload in generated_payloads.items():
            _write_gzip(private / filename, payload)
        copied = {
            pedigree_path: private / "pedigree_truth.tsv.gz",
            transmission_path: private / "mendelian_transmission_truth.tsv.gz",
            provenance_path: private / "background_d12_haplotype_provenance.tsv.gz",
        }
        for source, destination in copied.items():
            shutil.copyfile(source, destination)
        dense_buffer = temporary_root / "founder_panel.phased.vcf"
        observed_buffer = temporary_root / "blinded_observed_founder_panel.vcf"
        with (
            _open_text(dense_path) as source_handle,
            dense_buffer.open("w", encoding="utf-8", newline="") as dense_handle,
            observed_buffer.open(
                "w", encoding="utf-8", newline=""
            ) as observed_handle,
        ):
            stats = transform_founder_vcf(
                source_handle,
                dense_handle,
                observed_handle,
                plan=plan,
                observed_variant_keys=observed_keys,
                transmissions=transmission_records,
            )
        dense_output = private / "founder_panel.phased.vcf.gz"
        observed_output = private / "blinded_observed_founder_panel.vcf.gz"
        _write_gzip(dense_output, dense_buffer.read_bytes())
        _write_gzip(observed_output, observed_buffer.read_bytes())
        dense_buffer.unlink()
        observed_buffer.unlink()
        for path in private.iterdir():
            path.chmod(0o600)
        target_cm = interpolate_cm(genetic_map, configuration.target_position_bp)
        left_censored = sum(segment.left_censored for segment in plan.segments)
        right_censored = sum(segment.right_censored for segment in plan.segments)
        recombination_events = sum(
            segment.recombination_event_count for segment in plan.segments
        )
        panel = d13_manifest["panel"]
        artifacts = [
            _output_artifact(
                run_id,
                "dense_phased",
                "DENSE_PHASED_FOUNDER_PANEL",
                dense_output,
                stats.dense_variant_count,
            ),
            _output_artifact(
                run_id,
                "blinded_observed",
                "BLINDED_OBSERVED_FOUNDER_PANEL",
                observed_output,
                stats.observed_variant_count,
            ),
            _output_artifact(
                run_id,
                "segments",
                "FOUNDER_SEGMENT_TRUTH",
                private / "founder_segments.tsv.gz",
                len(plan.segments),
            ),
            _output_artifact(
                run_id,
                "carriers",
                "FOUNDER_CARRIER_TRUTH",
                private / "founder_carriers.tsv.gz",
                len(plan.carriers),
            ),
            _output_artifact(
                run_id,
                "target_truth",
                "TARGET_VARIANT_TRUTH",
                private / "target_variant_truth.tsv.gz",
                len(origins),
            ),
            _output_artifact(
                run_id,
                "local_ancestry",
                "LOCAL_ANCESTRY_MOSAIC_TRUTH",
                private / "local_ancestry_mosaic_truth.tsv.gz",
                len(plan.segments),
            ),
            _output_artifact(
                run_id,
                "pedigree",
                "PEDIGREE_TRUTH",
                private / "pedigree_truth.tsv.gz",
                pedigree_item["record_count"],
            ),
            _output_artifact(
                run_id,
                "transmissions",
                "MENDELIAN_TRANSMISSION_TRUTH",
                private / "mendelian_transmission_truth.tsv.gz",
                transmission_item["record_count"],
            ),
            _output_artifact(
                run_id,
                "background_provenance",
                "BACKGROUND_D12_HAPLOTYPE_PROVENANCE",
                private / "background_d12_haplotype_provenance.tsv.gz",
                provenance_item["record_count"],
            ),
        ]
        manifest = {
            "schema_version": "1.0.0",
            "founder_simulation_manifest_id": f"{run_id}.d14",
            "manifest_state": "FINALIZED",
            "created_at": created_at,
            "inputs": {
                "d13_manifest_id": d13_manifest["pedigree_manifest_id"],
                "d13_manifest_sha256": _sha256(d13_manifest_path),
                "d12_manifest_id": d12_manifest[
                    "population_structure_manifest_id"
                ],
                "d12_manifest_sha256": _sha256(d12_manifest_path),
                "d11a_manifest_id": d11_manifest["base_haplotype_manifest_id"],
                "d11a_manifest_sha256": _sha256(d11_manifest_path),
                "dense_sha256": dense_item["sha256"],
                "observed_sha256": observed_item["sha256"],
                "pedigree_sha256": pedigree_item["sha256"],
                "d13_provenance_sha256": provenance_item["sha256"],
                "transmission_sha256": transmission_item["sha256"],
                "d12_provenance_sha256": d12_provenance_item["sha256"],
                "genetic_map_sha256": map_item["sha256"],
            },
            "configuration": {
                "configuration_id": configuration.configuration_id,
                "configuration_sha256": _sha256(founder_config_path),
                "replication_id": configuration.replication_id,
                "random_seed": configuration.random_seed,
                "age_generations": configuration.age_generations,
            },
            "panel": {
                "role": panel["role"],
                "individual_count": len(origins),
                "family_count": len(plan.segments),
                "carrier_sample_count": stats.carrier_sample_count,
                "carrier_haplotype_count": stats.carrier_haplotype_count,
            },
            "domain": {
                "chromosome": configuration.chromosome,
                "start_bp": genetic_map[0].position_bp,
                "end_bp": genetic_map[-1].position_bp,
                "start_cm": genetic_map[0].position_cm,
                "end_cm": genetic_map[-1].position_cm,
                "target_position_bp": configuration.target_position_bp,
                "target_position_cm": target_cm,
                "position_unit": "BP_1_BASED_CLOSED",
                "genetic_distance_unit": "cM",
            },
            "method": {
                "method_id": FOUNDER_RECOMBINATION_METHOD_ID,
                "method_version": FOUNDER_RECOMBINATION_METHOD_VERSION,
                "breakpoint_distribution": "EXPONENTIAL_INDEPENDENT_EACH_SIDE",
                "rate_per_cm": configuration.age_generations / 100.0,
                "edge_policy": "RIGHT_CENSOR_AT_MAP_DOMAIN",
                "recent_pedigree_policy": "PROPAGATE_BY_D12_HAPLOTYPE_ORIGIN",
                "target_variant_policy": "RESET_EXISTING_SNV_TO_SYNTHETIC_MUTATION",
            },
            "confidentiality": {
                "real_sample_files_read": False,
                "truth_labels_in_blinded_panel": False,
                "local_source_paths_in_manifest": False,
                "artifacts_version_control_allowed": False,
            },
            "artifacts": artifacts,
            "integrity": {
                "dense_variant_count": stats.dense_variant_count,
                "observed_variant_count": stats.observed_variant_count,
                "target_variant_count": stats.target_variant_count,
                "family_count": len(plan.segments),
                "left_censored_count": left_censored,
                "right_censored_count": right_censored,
                "recombination_event_count": recombination_events,
                "changed_haplotype_allele_count": (
                    stats.changed_haplotype_allele_count
                ),
                "mendelian_error_count": stats.mendelian_error_count,
                "family_overlap_count": 0,
                "checks": list(_CHECKS),
            },
        }
        manifest_path = private / "founder_simulation_manifest.json"
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
        description="Matérialise le fondateur et les recombinaisons D14."
    )
    parser.add_argument("--d13-root", type=Path, required=True)
    parser.add_argument("--d12-root", type=Path, required=True)
    parser.add_argument("--d11a-root", type=Path, required=True)
    parser.add_argument("--founder-config", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    return parser


def main() -> int:
    """Point d'entrée D14 ; n'affiche que des compteurs agrégés."""
    arguments = _parser().parse_args()
    manifest = materialize_founder_recombination(
        d13_root=arguments.d13_root,
        d12_root=arguments.d12_root,
        d11a_root=arguments.d11a_root,
        founder_config_path=arguments.founder_config,
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
