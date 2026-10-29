"""Matérialisation atomique D14b des fenêtres génétiques de simulation."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

from effet_fondateur.calibration.simulation.chip_markers import (
    CHIP_ASSEMBLY,
    CHIP_DESIGN_ID,
    CHIP_DESIGN_VERSION,
    MarkerMatch,
    ReferenceVariant,
    harmonize_chip_markers,
    parse_chip_marker_manifest,
    parse_reference_vcf_variants,
)
from effet_fondateur.calibration.simulation.windows import (
    SIMULATION_WINDOW_METHOD_ID,
    SIMULATION_WINDOW_METHOD_VERSION,
    DesignedWindow,
    SimulationWindowConfiguration,
    SimulationWindowError,
    design_window,
    interpolate_window_cm,
    parse_public_genetic_map,
    parse_window_configuration,
)


_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_PUBLIC_REFERENCE_SOURCE_ID = "igsr_1kg_3202_phased_grch38_20220422"
_GENETIC_MAP_SOURCE_ID = "shapeit4_hapmap_grch38_git_7a0cab7"
_SOURCE_FILENAME = re.compile(
    r"^1kGP_high_coverage_Illumina\.chr(?P<chromosome>[1-9]|1[0-9]|2[0-2])"
    r"\.filtered\.SNV_INDEL_SV_phased_panel\.vcf\.gz$"
)
_CHECKS = (
    "PUBLIC_REFERENCE_HASHES_RECORDED",
    "GENETIC_MAP_MONOTONIC_AND_COVERING",
    "NO_GENETIC_MAP_EXTRAPOLATION",
    "TARGET_VARIANT_EXACTLY_MATCHED",
    "PRIMARY_ARM_CENSOR_PROBABILITY_WITHIN_LIMIT",
    "CHALLENGE_WINDOWS_SEPARATE",
    "IBD_MARKER_THRESHOLDS_MET",
    "NO_OUTCOME_FILTERING",
    "NO_REAL_STUDY_DATA",
)
CommandRunner = Callable[[list[str]], None]


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise SimulationWindowError("simulation_window_json_duplicate_key")
        value[key] = item
    return value


def _load_json(path: Path, code: str) -> dict[str, Any]:
    if not path.is_file() or path.is_symlink():
        raise SimulationWindowError(code)
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SimulationWindowError(code) from error
    if not isinstance(value, dict):
        raise SimulationWindowError(code)
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_public_or_derived_file(path: Path, code: str) -> None:
    if not path.is_file() or path.is_symlink():
        raise SimulationWindowError(code)
    resolved_parts = tuple(part.casefold() for part in path.resolve().parts)
    if "dock6_real" in resolved_parts:
        raise SimulationWindowError("simulation_window_real_input_forbidden")
    for left, right in zip(resolved_parts, resolved_parts[1:], strict=False):
        if left == "data" and right == "input":
            raise SimulationWindowError("simulation_window_real_input_forbidden")


def _run_command(command: list[str]) -> None:
    try:
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
    except OSError as error:
        raise SimulationWindowError("simulation_window_external_tool_unavailable") from error
    if completed.returncode != 0:
        executable = Path(command[0]).name
        raise SimulationWindowError(
            "simulation_window_external_tool_failed", f"{executable}:{completed.returncode}"
        )


def _extract_reference_window(
    *,
    source_vcf_path: Path,
    chromosome: int,
    start_bp: int,
    end_bp: int,
    destination: Path,
    bcftools_executable: str,
    command_runner: CommandRunner,
) -> None:
    command_runner(
        [
            bcftools_executable,
            "view",
            "--regions",
            f"chr{chromosome}:{start_bp}-{end_bp}",
            "--min-alleles",
            "2",
            "--max-alleles",
            "2",
            "--types",
            "snps",
            "--apply-filters",
            "PASS,.",
            "--output-type",
            "v",
            "--output-file",
            str(destination),
            str(source_vcf_path),
        ]
    )
    if not destination.is_file() or destination.stat().st_size == 0:
        raise SimulationWindowError("simulation_window_reference_extract_missing")


def _variant_in_window(
    variant: ReferenceVariant, window: DesignedWindow
) -> bool:
    return (
        variant.chromosome == window.specification.chromosome
        and window.start_bp <= variant.position_bp <= window.end_bp
    )


def _match_in_window(match: MarkerMatch, window: DesignedWindow) -> bool:
    return (
        match.marker.chromosome == window.specification.chromosome
        and window.start_bp <= match.marker.position_bp <= window.end_bp
    )


def _target_variant(
    variants: Sequence[ReferenceVariant], window: DesignedWindow
) -> ReferenceVariant:
    target = window.specification.target
    matches = [
        variant
        for variant in variants
        if variant.chromosome == window.specification.chromosome
        and variant.position_bp == target.position_bp
        and variant.reference_allele == target.reference_allele
        and variant.alternate_allele == target.alternate_allele
    ]
    if len(matches) != 1:
        raise SimulationWindowError(
            "simulation_window_target_variant_mismatch", window.specification.window_id
        )
    return matches[0]


def _artifact(
    run_id: str,
    suffix: str,
    role: str,
    relative_path: str,
    path: Path,
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


def materialize_simulation_windows(
    *,
    window_config_path: Path,
    source_vcf_path: Path,
    source_index_path: Path,
    genetic_map_path: Path,
    normalized_chip_markers_path: Path,
    workspace_root: Path,
    run_id: str,
    created_at: str,
    bcftools_executable: str = "bcftools",
    command_runner: CommandRunner = _run_command,
) -> dict[str, Any]:
    """Crée le manifest D14b et ses audits agrégés dans un nouveau run."""
    if _IDENTIFIER.fullmatch(run_id) is None:
        raise SimulationWindowError("simulation_window_run_id_invalid")
    if workspace_root.name != "calibration_workspace" or workspace_root.is_symlink():
        raise SimulationWindowError("simulation_window_workspace_invalid")
    workspace_root.mkdir(parents=True, exist_ok=True)
    target_root = workspace_root / run_id
    if target_root.exists() or target_root.is_symlink():
        raise SimulationWindowError("simulation_window_run_exists")
    for path, code in (
        (window_config_path, "simulation_window_config_file_invalid"),
        (source_vcf_path, "simulation_window_source_vcf_invalid"),
        (source_index_path, "simulation_window_source_index_invalid"),
        (genetic_map_path, "simulation_window_genetic_map_invalid"),
        (normalized_chip_markers_path, "simulation_window_chip_markers_invalid"),
    ):
        _validate_public_or_derived_file(path, code)
    source_match = _SOURCE_FILENAME.fullmatch(source_vcf_path.name)
    if source_match is None or source_index_path.name != f"{source_vcf_path.name}.tbi":
        raise SimulationWindowError("simulation_window_source_filename_invalid")
    configuration_document = _load_json(
        window_config_path, "simulation_window_config_file_invalid"
    )
    configuration = parse_window_configuration(configuration_document)
    chromosomes = {spec.chromosome for spec in configuration.window_specs}
    chromosome = next(iter(chromosomes))
    if int(source_match.group("chromosome")) != chromosome:
        raise SimulationWindowError("simulation_window_source_chromosome_mismatch")
    with genetic_map_path.open(encoding="utf-8") as handle:
        map_points = parse_public_genetic_map(
            handle, expected_chromosome=chromosome
        )
    designed_windows = tuple(
        design_window(
            specification,
            map_points,
            maximum_arm_censor_probability=(
                configuration.maximum_arm_censor_probability
            ),
        )
        for specification in configuration.window_specs
    )
    envelope_start = min(window.start_bp for window in designed_windows)
    envelope_end = max(window.end_bp for window in designed_windows)
    temporary_root = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        private_root = temporary_root / "private_truth/d14b"
        work_root = temporary_root / ".work"
        private_root.mkdir(parents=True, mode=0o700)
        work_root.mkdir(mode=0o700)
        reference_extract = work_root / "reference_window.vcf"
        _extract_reference_window(
            source_vcf_path=source_vcf_path,
            chromosome=chromosome,
            start_bp=envelope_start,
            end_bp=envelope_end,
            destination=reference_extract,
            bcftools_executable=bcftools_executable,
            command_runner=command_runner,
        )
        with reference_extract.open(encoding="utf-8") as handle:
            variants = parse_reference_vcf_variants(
                handle,
                chromosome=chromosome,
                start_bp=envelope_start,
                end_bp=envelope_end,
            )
        with normalized_chip_markers_path.open(encoding="utf-8") as handle:
            chip_markers = parse_chip_marker_manifest(handle)
        marker_matches, _ = harmonize_chip_markers(
            chip_markers,
            variants,
            chromosome=chromosome,
            start_bp=envelope_start,
            end_bp=envelope_end,
        )
        window_records: list[dict[str, Any]] = []
        audit_rows: list[tuple[object, ...]] = []
        not_evaluable_count = 0
        for window in designed_windows:
            _target_variant(variants, window)
            window_variants = tuple(
                variant for variant in variants if _variant_in_window(variant, window)
            )
            window_matches = tuple(
                match for match in marker_matches if _match_in_window(match, window)
            )
            if not window_variants or not window_matches:
                raise SimulationWindowError(
                    "simulation_window_marker_universe_empty",
                    window.specification.window_id,
                )
            positions_cm = {
                match.marker.position_bp: interpolate_window_cm(
                    map_points, match.marker.position_bp
                )
                for match in window_matches
            }
            evaluations: list[dict[str, Any]] = []
            if window.specification.purpose == "PRIMARY":
                for requirement in configuration.ibd_requirements:
                    segment_start_cm = (
                        window.target_cm - requirement.segment_length_cm / 2.0
                    )
                    segment_end_cm = (
                        window.target_cm + requirement.segment_length_cm / 2.0
                    )
                    inside_count = sum(
                        segment_start_cm <= position_cm <= segment_end_cm
                        for position_cm in positions_cm.values()
                    )
                    left_flank_count = sum(
                        window.start_cm <= position_cm < segment_start_cm
                        for position_cm in positions_cm.values()
                    )
                    right_flank_count = sum(
                        segment_end_cm < position_cm <= window.end_cm
                        for position_cm in positions_cm.values()
                    )
                    status = (
                        "EVALUABLE"
                        if segment_start_cm > window.start_cm
                        and segment_end_cm < window.end_cm
                        and inside_count >= requirement.minimum_observed_markers
                        and left_flank_count
                        >= configuration.minimum_flank_markers_each_side
                        and right_flank_count
                        >= configuration.minimum_flank_markers_each_side
                        else "NOT_EVALUABLE"
                    )
                    if status == "NOT_EVALUABLE":
                        not_evaluable_count += 1
                    evaluation = {
                        "scenario_id": requirement.scenario_id,
                        "segment_length_cm": requirement.segment_length_cm,
                        "minimum_observed_markers": requirement.minimum_observed_markers,
                        "observed_marker_count": inside_count,
                        "minimum_flank_markers_each_side": (
                            configuration.minimum_flank_markers_each_side
                        ),
                        "left_flank_marker_count": left_flank_count,
                        "right_flank_marker_count": right_flank_count,
                        "status": status,
                    }
                    evaluations.append(evaluation)
                    audit_rows.append(
                        (
                            window.specification.window_id,
                            "IBD_SCENARIO",
                            requirement.scenario_id,
                            requirement.segment_length_cm,
                            requirement.minimum_observed_markers,
                            inside_count,
                            configuration.minimum_flank_markers_each_side,
                            left_flank_count,
                            right_flank_count,
                            status,
                        )
                    )
            audit_rows.append(
                (
                    window.specification.window_id,
                    "WINDOW_SUMMARY",
                    "NA",
                    "NA",
                    "NA",
                    len(window_matches),
                    "NA",
                    sum(position < window.target_cm for position in positions_cm.values()),
                    sum(position > window.target_cm for position in positions_cm.values()),
                    "EVALUABLE",
                )
            )
            window_records.append(
                {
                    "window_id": window.specification.window_id,
                    "purpose": window.specification.purpose,
                    "age_generations": window.specification.age_generations,
                    "chromosome": chromosome,
                    "target": {
                        "variant_id": window.specification.target.variant_id,
                        "position_bp": window.specification.target.position_bp,
                        "position_cm": window.target_cm,
                        "reference_allele": window.specification.target.reference_allele,
                        "alternate_allele": window.specification.target.alternate_allele,
                    },
                    "domain": {
                        "start_bp": window.start_bp,
                        "end_bp": window.end_bp,
                        "start_cm": window.start_cm,
                        "end_cm": window.end_cm,
                    },
                    "requested_half_width_cm": window.requested_half_width_cm,
                    "realized_left_width_cm": window.realized_left_width_cm,
                    "realized_right_width_cm": window.realized_right_width_cm,
                    "theoretical_arm_censor_probability": (
                        window.theoretical_arm_censor_probability
                    ),
                    "dense_variant_count": len(window_variants),
                    "observed_marker_count": len(window_matches),
                    "left_observed_marker_count": sum(
                        position < window.target_cm for position in positions_cm.values()
                    ),
                    "right_observed_marker_count": sum(
                        position > window.target_cm for position in positions_cm.values()
                    ),
                    "ibd_scenarios": evaluations,
                }
            )
        if not_evaluable_count:
            raise SimulationWindowError(
                "simulation_window_ibd_marker_threshold_not_met",
                str(not_evaluable_count),
            )
        definitions_path = private_root / "window_definitions.tsv"
        with definitions_path.open("w", encoding="utf-8", newline="") as handle:
            handle.write(
                "WINDOW_ID\tPURPOSE\tAGE_GENERATIONS\tCHROMOSOME\tSTART_BP\t"
                "TARGET_BP\tEND_BP\tSTART_CM\tTARGET_CM\tEND_CM\t"
                "REQUESTED_HALF_WIDTH_CM\tARM_CENSOR_PROBABILITY\t"
                "DENSE_VARIANT_COUNT\tOBSERVED_MARKER_COUNT\n"
            )
            for record in window_records:
                domain = record["domain"]
                target = record["target"]
                handle.write(
                    "\t".join(
                        str(value)
                        for value in (
                            record["window_id"],
                            record["purpose"],
                            record["age_generations"],
                            record["chromosome"],
                            domain["start_bp"],
                            target["position_bp"],
                            domain["end_bp"],
                            domain["start_cm"],
                            target["position_cm"],
                            domain["end_cm"],
                            record["requested_half_width_cm"],
                            record["theoretical_arm_censor_probability"],
                            record["dense_variant_count"],
                            record["observed_marker_count"],
                        )
                    )
                    + "\n"
                )
        audit_path = private_root / "observable_marker_audit.tsv"
        with audit_path.open("w", encoding="utf-8", newline="") as handle:
            handle.write(
                "WINDOW_ID\tRECORD_TYPE\tSCENARIO_ID\tSEGMENT_LENGTH_CM\t"
                "MINIMUM_OBSERVED_MARKERS\tOBSERVED_MARKER_COUNT\t"
                "MINIMUM_FLANK_MARKERS_EACH_SIDE\tLEFT_FLANK_MARKER_COUNT\t"
                "RIGHT_FLANK_MARKER_COUNT\tSTATUS\n"
            )
            for row in audit_rows:
                handle.write("\t".join(str(value) for value in row) + "\n")
        definitions_path.chmod(0o600)
        audit_path.chmod(0o600)
        artifacts = [
            _artifact(
                run_id,
                "windows",
                "WINDOW_DEFINITION_TABLE",
                "private_truth/d14b/window_definitions.tsv",
                definitions_path,
                len(window_records),
            ),
            _artifact(
                run_id,
                "markers",
                "OBSERVABLE_MARKER_AUDIT",
                "private_truth/d14b/observable_marker_audit.tsv",
                audit_path,
                len(audit_rows),
            ),
        ]
        manifest = {
            "schema_version": "1.0.0",
            "simulation_window_manifest_id": f"{run_id}.d14b",
            "manifest_state": "FINALIZED",
            "created_at": created_at,
            "configuration": {
                "configuration_id": configuration.configuration_id,
                "configuration_sha256": _sha256(window_config_path),
            },
            "panel": {"role": configuration.panel_role},
            "inputs": {
                "assembly": CHIP_ASSEMBLY,
                "public_reference_source_id": _PUBLIC_REFERENCE_SOURCE_ID,
                "source_vcf_filename": source_vcf_path.name,
                "source_vcf_sha256": _sha256(source_vcf_path),
                "source_index_sha256": _sha256(source_index_path),
                "genetic_map_source_id": _GENETIC_MAP_SOURCE_ID,
                "genetic_map_sha256": _sha256(genetic_map_path),
                "chip_design_id": CHIP_DESIGN_ID,
                "chip_design_version": CHIP_DESIGN_VERSION,
                "normalized_chip_markers_sha256": _sha256(
                    normalized_chip_markers_path
                ),
            },
            "method": {
                "method_id": SIMULATION_WINDOW_METHOD_ID,
                "method_version": SIMULATION_WINDOW_METHOD_VERSION,
                "maximum_arm_censor_probability": (
                    configuration.maximum_arm_censor_probability
                ),
                "half_width_formula": (
                    "-100*ln(maximum_arm_censor_probability)/age_generations"
                ),
                "genetic_map_interpolation": (
                    "PIECEWISE_LINEAR_NO_EXTRAPOLATION"
                ),
                "outcome_filtering": "FORBIDDEN",
            },
            "windows": window_records,
            "confidentiality": {
                "real_sample_files_read": False,
                "sample_identifiers_recorded": False,
                "local_source_paths_in_manifest": False,
                "artifacts_version_control_allowed": False,
            },
            "artifacts": artifacts,
            "integrity": {
                "window_count": len(window_records),
                "primary_window_count": sum(
                    record["purpose"] == "PRIMARY" for record in window_records
                ),
                "challenge_window_count": sum(
                    record["purpose"] == "CENSORING_CHALLENGE"
                    for record in window_records
                ),
                "infeasible_window_count": 0,
                "ibd_scenario_evaluation_count": sum(
                    len(record["ibd_scenarios"]) for record in window_records
                ),
                "ibd_not_evaluable_count": 0,
                "checks": list(_CHECKS),
            },
        }
        manifest_path = private_root / "simulation_window_manifest.json"
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
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Matérialise les fenêtres génétiques D14b sans donnée d'étude."
    )
    parser.add_argument("--window-config", type=Path, required=True)
    parser.add_argument("--source-vcf", type=Path, required=True)
    parser.add_argument("--source-index", type=Path, required=True)
    parser.add_argument("--genetic-map", type=Path, required=True)
    parser.add_argument("--normalized-chip-markers", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--bcftools", default="bcftools")
    return parser


def main() -> int:
    """Point d'entrée D14b ; n'affiche que les contrôles agrégés."""
    arguments = _parser().parse_args()
    manifest = materialize_simulation_windows(
        window_config_path=arguments.window_config,
        source_vcf_path=arguments.source_vcf,
        source_index_path=arguments.source_index,
        genetic_map_path=arguments.genetic_map,
        normalized_chip_markers_path=arguments.normalized_chip_markers,
        workspace_root=arguments.workspace_root,
        run_id=arguments.run_id,
        created_at=datetime.now(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z"),
        bcftools_executable=arguments.bcftools,
    )
    print(json.dumps(manifest["integrity"], ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
