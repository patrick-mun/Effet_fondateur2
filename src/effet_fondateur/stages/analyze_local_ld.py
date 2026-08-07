"""Étape 15 : description secondaire du LD dans la région cible."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    TableValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_cohorts_frozen,
    validate_json_document,
    validate_tsv_table,
)
from effet_fondateur.ld import LocalLdError, publish_local_ld
from effet_fondateur.orchestrator.state import utc_now


COHORT_IDS = (
    "controls_unrelated",
    "target_carriers_independent",
    "family_noncarriers",
)


class AnalyzeLocalLdInputError(ValueError):
    """Signale une configuration ou un artefact d'entrée invalide."""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="analyze-local-ld")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise AnalyzeLocalLdInputError(f"invalid_parameter:{name}")
    return value


def _positive_number(parameters: dict[str, Any], name: str, default: float) -> float:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AnalyzeLocalLdInputError(f"invalid_parameter:{name}")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise AnalyzeLocalLdInputError(f"invalid_parameter:{name}")
    return result


def _maf(parameters: dict[str, Any], name: str, default: float) -> float:
    result = _positive_number(parameters, name, default)
    if result > 0.5:
        raise AnalyzeLocalLdInputError(f"invalid_parameter:{name}")
    return result


def _parameters(parameters: dict[str, Any]) -> dict[str, Any]:
    if parameters.get("method", "plink19_local_ld_secondary_v1") != "plink19_local_ld_secondary_v1":
        raise AnalyzeLocalLdInputError("invalid_parameter:method")
    minimum_samples = _positive_integer(parameters, "minimum_samples", 5)
    primary_samples = _positive_integer(parameters, "minimum_primary_samples", 20)
    if primary_samples < minimum_samples:
        raise AnalyzeLocalLdInputError("invalid_parameter:minimum_primary_samples")
    max_distance_cm = _positive_number(parameters, "max_pair_distance_cm", 1.0)
    raw_bins = parameters.get("distance_bins_cm", [0.1, 0.25, 0.5, 1.0])
    if (
        not isinstance(raw_bins, list) or not raw_bins
        or any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in raw_bins)
    ):
        raise AnalyzeLocalLdInputError("invalid_parameter:distance_bins_cm")
    bins = [float(value) for value in raw_bins]
    if (
        any(not math.isfinite(value) or value <= 0 for value in bins)
        or bins != sorted(set(bins))
        or not math.isclose(bins[-1], max_distance_cm, rel_tol=0, abs_tol=1e-12)
    ):
        raise AnalyzeLocalLdInputError("invalid_parameter:distance_bins_cm")
    return {
        "method": "plink19_local_ld_secondary_v1",
        "minimum_samples": minimum_samples,
        "minimum_primary_samples": primary_samples,
        "minimum_called_samples_per_pair": _positive_integer(parameters, "minimum_called_samples_per_pair", 5),
        "minimum_variant_maf": _maf(parameters, "minimum_variant_maf", 0.05),
        "minimum_pairs_per_bin": _positive_integer(parameters, "minimum_pairs_per_bin", 10),
        "max_pair_distance_bp": _positive_integer(parameters, "max_pair_distance_bp", 1_000_000),
        "max_pair_distance_cm": max_distance_cm,
        "distance_bins_cm": bins,
        "max_pair_count": _positive_integer(parameters, "max_pair_count", 1_000_000),
        "plink_timeout_seconds": _positive_integer(parameters, "plink_timeout_seconds", 300),
    }


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [artifact for artifact in stage_inputs["artifacts"] if artifact["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise AnalyzeLocalLdInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _validated_path(artifact: dict[str, Any], run_dir: Path) -> Path:
    path = Path(artifact["path"])
    if not path.is_absolute() and PurePosixPath(artifact["path"]).parts[:1] == ("stages",):
        path = run_dir / path
    elif not path.is_absolute():
        path = Path.cwd() / path
    if not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise AnalyzeLocalLdInputError(f"artifact_missing_or_modified:{artifact['artifact_id']}")
    return path


def _read_fam(path: Path) -> list[tuple[str, str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 6 for row in rows):
        raise AnalyzeLocalLdInputError("malformed_target_region_fam")
    pairs = [(row[0], row[1]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise AnalyzeLocalLdInputError("duplicate_target_region_sample")
    return pairs


def _read_bim(path: Path) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if len(rows) < 2 or any(len(row) != 6 for row in rows):
        raise AnalyzeLocalLdInputError("malformed_target_region_bim")
    if len({row[1] for row in rows}) != len(rows):
        raise AnalyzeLocalLdInputError("duplicate_target_region_variant")
    if any(row[4] == "0" or row[5] == "0" or row[4] == row[5] for row in rows):
        raise AnalyzeLocalLdInputError("non_biallelic_or_unresolved_region_variant")
    return rows


def _read_keep(path: Path) -> list[tuple[str, str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if any(len(row) != 2 for row in rows):
        raise AnalyzeLocalLdInputError(f"malformed_keep:{path.name}")
    pairs = [(row[0], row[1]) for row in rows]
    if len(pairs) != len(set(pairs)):
        raise AnalyzeLocalLdInputError(f"duplicate_keep_sample:{path.name}")
    return pairs


def _validate_region_dataset(
    paths: dict[str, Path], artifacts: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], list[tuple[str, str]], list[list[str]]]:
    descriptor = read_json(paths["target_region_dataset"])
    validate_json_document(descriptor, "plink_dataset.schema.json")
    if paths["target_region_bed"].read_bytes()[:3] != b"\x6c\x1b\x01":
        raise AnalyzeLocalLdInputError("invalid_target_region_bed")
    prefixes = {paths[f"target_region_{suffix}"].with_suffix("") for suffix in ("bed", "bim", "fam")}
    if len(prefixes) != 1:
        raise AnalyzeLocalLdInputError("target_region_prefix_mismatch")
    fam = _read_fam(paths["target_region_fam"])
    bim = _read_bim(paths["target_region_bim"])
    if descriptor["sample_count"] != len(fam) or descriptor["variant_count"] != len(bim):
        raise AnalyzeLocalLdInputError("target_region_descriptor_count_mismatch")
    for suffix in ("bed", "bim", "fam"):
        artifact = artifacts[f"target_region_{suffix}"]
        if descriptor["files"][suffix]["sha256"] != artifact["sha256"]:
            raise AnalyzeLocalLdInputError("target_region_descriptor_checksum_mismatch")
    return descriptor, fam, bim


def _variant_universe(paths: dict[str, Path], bim_rows: list[list[str]]) -> list[dict[str, Any]]:
    genetic_map = validate_tsv_table(paths["target_genetic_map"], "target_genetic_map.schema.json")
    map_by_id = {row["VARIANT_ID"]: row for row in genetic_map.rows}
    qc_rows = validate_tsv_table(paths["variant_qc_final"], "variant_qc_final.schema.json").rows
    retained_qc = {
        row["VARIANT_ID"] for row in qc_rows
        if row["DATASET_ID"] == "target_chromosome_all_qc" and row["PRESENT_AFTER_QC"]
    }
    bim_ids = [row[1] for row in bim_rows]
    if set(bim_ids) != set(map_by_id) or not set(bim_ids) <= retained_qc:
        raise AnalyzeLocalLdInputError("local_variant_universe_mismatch")
    variants = []
    for order, bim in enumerate(bim_rows, start=1):
        map_row = map_by_id[bim[1]]
        if (
            int(map_row["VARIANT_ORDER"]) != order
            or map_row["CHROMOSOME"] != bim[0]
            or int(map_row["POSITION_BP"]) != int(bim[3])
            or not math.isclose(float(map_row["POSITION_CM"]), float(bim[2]), rel_tol=0, abs_tol=1e-8)
        ):
            raise AnalyzeLocalLdInputError("local_map_and_bim_mismatch")
        variants.append({
            "VARIANT_ORDER": order, "VARIANT_ID": bim[1], "CHROMOSOME": bim[0],
            "POSITION_BP": int(bim[3]), "POSITION_CM": float(map_row["POSITION_CM"]),
            "A1": bim[4], "A2": bim[5], "IS_TARGET_VARIANT": map_row["IS_TARGET_VARIANT"],
        })
    return variants


def _cohort_context(
    paths: dict[str, Path], fam_pairs: list[tuple[str, str]],
) -> tuple[dict[str, int], dict[str, bool]]:
    cohort_table = validate_cohorts_frozen(paths["cohorts_frozen"])
    rows_by_cohort: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in cohort_table.rows:
        if row["COHORT_ID"] in COHORT_IDS and row["INCLUDED"]:
            rows_by_cohort[row["COHORT_ID"]].append(row)
    fam_set = set(fam_pairs)
    counts: dict[str, int] = {}
    independent: dict[str, bool] = {}
    for cohort_id in COHORT_IDS:
        keep_pairs = _read_keep(paths[f"cohort_keep_{cohort_id}"])
        if not set(keep_pairs) <= fam_set or len(keep_pairs) != len(rows_by_cohort[cohort_id]):
            raise AnalyzeLocalLdInputError(f"cohort_keep_mismatch:{cohort_id}")
        units = [row["INDEPENDENT_UNIT_ID"] for row in rows_by_cohort[cohort_id]]
        counts[cohort_id] = len(keep_pairs)
        independent[cohort_id] = all(unit is not None for unit in units) and len(units) == len(set(units))
    return counts, independent


def _resolve_plink(command: str | None) -> tuple[str, str]:
    if command is None or shutil.which(command) is None:
        raise LocalLdError("plink_not_available")
    executable = shutil.which(command)
    assert executable is not None
    try:
        completed = subprocess.run([executable, "--version"], capture_output=True, check=False, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise LocalLdError("plink_version_failed") from error
    version = (completed.stdout or completed.stderr).strip().splitlines()
    if completed.returncode != 0 or not version or "1.9" not in version[0]:
        raise LocalLdError("plink19_required")
    return executable, version[0][:200]


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Publie le LD secondaire sans lire les sorties des étapes 13 ou 14."""
    started_at = utc_now()
    started_clock = monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    parameters = _parameters(stage_inputs["parameters"])
    artifact_ids = (
        "target_region_bed", "target_region_bim", "target_region_fam",
        "target_region_dataset", "target_genetic_map", "variant_qc_final",
        "cohorts_frozen", "cohort_keep_controls_unrelated",
        "cohort_keep_target_carriers_independent", "cohort_keep_family_noncarriers",
    )
    artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in artifact_ids}
    paths = {artifact_id: _validated_path(artifact, run_dir) for artifact_id, artifact in artifacts.items()}
    descriptor, fam_rows, bim_rows = _validate_region_dataset(paths, artifacts)
    variants = _variant_universe(paths, bim_rows)
    cohort_counts, cohort_independent = _cohort_context(paths, fam_rows)
    executable, plink_version = _resolve_plink(config["tools"]["plink"])
    publication = publish_local_ld(
        plink_executable=executable, dataset_prefix=paths["target_region_bed"].with_suffix(""),
        variants=variants,
        cohort_keep_paths={cohort_id: paths[f"cohort_keep_{cohort_id}"] for cohort_id in COHORT_IDS},
        cohort_sample_counts=cohort_counts, cohort_independent=cohort_independent,
        output_dir=output_dir / "local_ld", minimum_samples=parameters["minimum_samples"],
        primary_samples=parameters["minimum_primary_samples"],
        minimum_called_samples_per_pair=parameters["minimum_called_samples_per_pair"],
        minimum_variant_maf=parameters["minimum_variant_maf"],
        minimum_pairs_per_bin=parameters["minimum_pairs_per_bin"],
        max_pair_distance_bp=parameters["max_pair_distance_bp"],
        max_pair_distance_cm=parameters["max_pair_distance_cm"],
        distance_bins_cm=parameters["distance_bins_cm"], max_pair_count=parameters["max_pair_count"],
        timeout_seconds=parameters["plink_timeout_seconds"],
    )
    specifications = [
        ("local_ld_variants", publication.variants_path, "text/tab-separated-values", "local_ld_variants.schema.json", "sensitive_genetic"),
        ("local_ld_pairs", publication.pairs_path, "text/tab-separated-values", "local_ld_pairs.schema.json", "sensitive_genetic"),
        ("local_ld_summary", publication.summary_path, "text/tab-separated-values", "local_ld_summary.schema.json", "internal"),
        ("local_ld_analysis_summary", publication.analysis_summary_path, "application/json", "local_ld_analysis_summary.schema.json", "internal"),
    ]
    for native_path in publication.native_paths:
        specifications.append((f"local_ld_native_{native_path.name.replace('.', '_')}", native_path, "text/plain", None, "sensitive_genetic"))
    output_artifacts = [
        build_file_artifact(
            physical_path=path,
            published_path=f"{stage_inputs['published_output_dir']}/{path.relative_to(output_dir).as_posix()}",
            artifact_id=artifact_id, artifact_type=artifact_id, media_type=media_type,
            producer_stage=stage_inputs["stage_name"], producer_signature=stage_inputs["signature"],
            schema_name=schema_name, schema_version="1.0.0" if schema_name else None,
            assembly=descriptor["assembly"], sample_set_id=None,
            variant_set_id=descriptor["variant_set_id"], sensitivity=sensitivity,
        )
        for artifact_id, path, media_type, schema_name, sensitivity in specifications
    ]
    stage_outputs = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "signature": stage_inputs["signature"], "artifacts": output_artifacts,
    }
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    warnings = []
    for cohort_id, status in publication.cohort_statuses.items():
        if status != "DESCRIPTIVE_PRIMARY":
            warnings.append({"code": f"{cohort_id}_{status.lower()}", "count": 1})
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"],
        "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"],
        "method_id": "plink19_local_ld_secondary_v1", "signature": stage_inputs["signature"],
        "started_at": started_at, "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock,
        "inputs": list(artifacts.values()), "outputs": output_artifacts,
        "parameters": parameters,
        "tools": [{"tool": "plink", "configured": config["tools"]["plink"], "version": plink_version}],
        "counts": {
            "region_variants": len(variants), "pair_rows": publication.pair_count,
            "evaluated_pairs": publication.evaluated_pair_count,
            "cohort_samples": cohort_counts,
        },
        "metrics": {
            "analysis_role": "SECONDARY_DESCRIPTIVE",
            "cohort_statuses": publication.cohort_statuses,
            "r2_definition": "GENOTYPE_ALLELE_COUNT_CORRELATION_SQUARED",
            "dprime_definition": "ABSOLUTE_HAPLOTYPE_ML",
            "feeds_founder_haplotype": False, "feeds_variant_age": False,
        },
        "exclusions": [], "warnings": warnings,
        "checks": [
            {"check": "input_artifact_integrity", "status": "PASS"},
            {"check": "common_variant_universe", "status": "PASS"},
            {"check": "frozen_independent_cohorts", "status": "PASS"},
            {"check": "separate_native_r2_and_dprime", "status": "PASS"},
            {"check": "stage_13_14_outputs_not_consumed", "status": "PASS"},
        ],
        "known_limits": [
            "Le LD local est descriptif et ne constitue pas une preuve IBD, d'apparentement ou d'origine fondatrice unique.",
            "Le D-prime est instable aux faibles effectifs et pour les allèles rares ; les cohortes sous les seuils restent exploratoires ou non évaluées.",
            "Les paires impliquant le variant cible sont conditionnées par la définition porteur/non-porteur des cohortes.",
            "La densité ACPA, le QC et la résolution de la carte limitent la structure de LD observable.",
            "Aucun test massif, bloc haplotypique ou ajustement causal n'est produit par cette étape.",
        ],
        "expected_visualizations": ["plot_local_ld_matrices", "plot_local_ld_decay", "plot_local_ld_target_zoom"],
        "manual_validation_required": True,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text(
        "".join(
            f"{artifact['sha256']}  {artifact['path'].removeprefix(stage_inputs['published_output_dir'] + '/')}\n"
            for artifact in output_artifacts
        ), encoding="utf-8",
    )
    return 0


def _classify_error(error: Exception) -> int:
    if isinstance(error, LocalLdError):
        message = str(error)
        if message in {"configured_ld_pair_limit_exceeded", "insufficient_region_variants"}:
            return 4
        if message in {
            "duplicate_region_variant",
            "target_variant_missing_or_ambiguous",
        }:
            return 2
        # Les autres erreurs de ce module apparaissent après un appel PLINK :
        # binaire/version, rapport absent, colonne invalide ou ensemble altéré.
        return 3
    if isinstance(error, (AnalyzeLocalLdInputError, DocumentValidationError, TableValidationError, OSError, ValueError, json.JSONDecodeError)):
        return 2
    return 5


def main(arguments: Sequence[str] | None = None) -> int:
    """Point d'entrée autonome respectant les codes de retour V2."""
    parsed = _build_parser().parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except Exception as error:
        sys.stderr.write(f"{error}\n")
        return _classify_error(error)


if __name__ == "__main__":
    raise SystemExit(main())
