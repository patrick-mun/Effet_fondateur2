"""Construit le panel autosomal homogène utilisé par le phasage et l'IBD."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path, PurePosixPath
from time import monotonic
from typing import Any, Sequence

from effet_fondateur.audit import atomic_write_json, read_json, sha256_file
from effet_fondateur.contracts import (
    DocumentValidationError,
    build_file_artifact,
    load_pipeline_config,
    validate_json_document,
)
from effet_fondateur.orchestrator.state import utc_now


class AutosomalPanelInputError(ValueError):
    """Signale une entrée ou un paramètre invalide."""


class AutosomalPanelExternalError(RuntimeError):
    """Signale un échec de PLINK."""


class AutosomalPanelBlockError(RuntimeError):
    """Signale un panel insuffisant pour une analyse autosomale comparable."""


def _artifact_by_id(stage_inputs: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    matches = [item for item in stage_inputs["artifacts"] if item["artifact_id"] == artifact_id]
    if len(matches) != 1:
        raise AutosomalPanelInputError(f"{artifact_id}_missing_or_ambiguous")
    return matches[0]


def _path(artifact: dict[str, Any], run_dir: Path) -> Path:
    candidate = Path(artifact["path"])
    path = run_dir / candidate if PurePosixPath(artifact["path"]).parts[:1] == ("stages",) else candidate
    if not path.is_absolute():
        path = Path.cwd() / path
    if path.is_symlink() or not path.is_file() or sha256_file(path) != artifact["sha256"]:
        raise AutosomalPanelInputError("declared_input_modified")
    return path


def _positive_integer(parameters: dict[str, Any], name: str, default: int) -> int:
    value = parameters.get(name, default)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise AutosomalPanelInputError(f"invalid_parameter:{name}")
    return value


def _run(command: list[str], timeout: int, code: str) -> subprocess.CompletedProcess[str]:
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise AutosomalPanelExternalError(code) from error
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "")[-500:].replace("\n", " ")
        raise AutosomalPanelExternalError(f"{code}:{completed.returncode}:{detail}")
    return completed


def _rows(path: Path, columns: int, code: str) -> list[list[str]]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != columns for row in rows):
        raise AutosomalPanelInputError(code)
    return rows


def _set_id(prefix: str, values: Sequence[str]) -> str:
    digest = hashlib.sha256("".join(f"{value}\n" for value in values).encode()).hexdigest()
    return f"{prefix}_{digest[:12]}"


def _dataset_artifacts(
    stage_inputs: dict[str, Any], prefix: Path, descriptor_path: Path, descriptor: dict[str, Any]
) -> list[dict[str, Any]]:
    published = PurePosixPath(stage_inputs["published_output_dir"])
    producer = f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}"
    artifacts: list[dict[str, Any]] = []
    for suffix, media_type in (("bed", "application/octet-stream"), ("bim", "text/plain"), ("fam", "text/plain")):
        path = prefix.with_suffix(f".{suffix}")
        artifacts.append(build_file_artifact(
            physical_path=path, published_path=str(published / path.name),
            artifact_id=f"autosomal_phasing_panel_{suffix}", artifact_type=f"plink_{suffix}",
            media_type=media_type, producer_stage=producer,
            producer_signature=stage_inputs["signature"], assembly=descriptor["assembly"],
            sample_set_id=descriptor["sample_set_id"], variant_set_id=descriptor["variant_set_id"],
            sensitivity="sensitive_genetic",
        ))
    artifacts.append(build_file_artifact(
        physical_path=descriptor_path, published_path=str(published / descriptor_path.name),
        artifact_id="autosomal_phasing_panel_dataset", artifact_type="plink_dataset_descriptor",
        media_type="application/json", producer_stage=producer,
        producer_signature=stage_inputs["signature"], schema_name="plink_dataset.schema.json",
        schema_version="1.0.0", assembly=descriptor["assembly"],
        sample_set_id=descriptor["sample_set_id"], variant_set_id=descriptor["variant_set_id"],
        sensitivity="sensitive_genetic",
    ))
    return artifacts


def execute(stage_inputs_path: Path, output_dir: Path) -> int:
    """Filtre une seule fois les 22 autosomes pour les 75 mêmes individus."""
    started_at, started_clock = utc_now(), monotonic()
    stage_inputs = read_json(stage_inputs_path)
    validate_json_document(stage_inputs, "stage_inputs.schema.json")
    parameters = stage_inputs["parameters"]
    timeout = _positive_integer(parameters, "plink_timeout_seconds", 600)
    minimum_per_autosome = _positive_integer(parameters, "minimum_variants_per_autosome", 100)
    minimum_total = _positive_integer(parameters, "minimum_total_variants", 10_000)
    run_dir = output_dir.parent.parent
    config = load_pipeline_config(run_dir / "config.resolved.yaml")
    ids = (
        "genomewide_pre_qc_bed", "genomewide_pre_qc_bim", "genomewide_pre_qc_fam",
        "genomewide_pre_qc_dataset", "cohort_keep_target_chromosome_all_qc",
    )
    artifacts = {artifact_id: _artifact_by_id(stage_inputs, artifact_id) for artifact_id in ids}
    paths = {artifact_id: _path(artifact, run_dir) for artifact_id, artifact in artifacts.items()}
    source_descriptor = read_json(paths["genomewide_pre_qc_dataset"])
    validate_json_document(source_descriptor, "plink_dataset.schema.json")
    keep_rows = _rows(paths["cohort_keep_target_chromosome_all_qc"], 2, "cohort_keep_malformed")
    if len(keep_rows) != len({tuple(row) for row in keep_rows}):
        raise AutosomalPanelInputError("cohort_keep_duplicate")
    plink = shutil.which(config["tools"]["plink"])
    if plink is None:
        raise AutosomalPanelExternalError("plink_not_found")
    version_result = _run([plink, "--version"], timeout, "plink_version_failed")
    plink_version = (version_result.stdout or version_result.stderr).splitlines()[0]
    output_dir.mkdir(parents=True, exist_ok=True)
    source_prefix = paths["genomewide_pre_qc_bed"].with_suffix("")
    panel_prefix = output_dir / "autosomal_phasing_panel"
    _run([
        plink, "--bfile", str(source_prefix), "--keep", str(paths["cohort_keep_target_chromosome_all_qc"]),
        "--autosome", "--geno", "0", "--maf", "0.000001", "--make-bed", "--out", str(panel_prefix),
    ], timeout, "plink_autosomal_panel_failed")
    fam_rows = _rows(panel_prefix.with_suffix(".fam"), 6, "autosomal_panel_fam_malformed")
    bim_rows = _rows(panel_prefix.with_suffix(".bim"), 6, "autosomal_panel_bim_malformed")
    observed_samples = {(row[0], row[1]) for row in fam_rows}
    if observed_samples != {tuple(row) for row in keep_rows}:
        raise AutosomalPanelBlockError("autosomal_panel_sample_universe_mismatch")
    chromosomes = Counter(int(row[0].removeprefix("chr")) for row in bim_rows)
    if set(chromosomes) != set(range(1, 23)):
        raise AutosomalPanelBlockError("autosomal_panel_missing_chromosome")
    if min(chromosomes.values()) < minimum_per_autosome or len(bim_rows) < minimum_total:
        raise AutosomalPanelBlockError("autosomal_panel_marker_density_insufficient")
    if len({row[1] for row in bim_rows}) != len(bim_rows):
        raise AutosomalPanelBlockError("autosomal_panel_duplicate_variant_id")
    descriptor = {
        "schema_version": "1.0.0", "dataset_id": "autosomal_phasing_panel",
        "assembly": source_descriptor["assembly"], "scope": "autosomal_genomewide",
        "target_chromosome": None,
        "sample_set_id": _set_id("autosomal_phasing_samples", [f"{row[0]}\t{row[1]}" for row in fam_rows]),
        "variant_set_id": _set_id("autosomal_phasing_variants", [row[1] for row in bim_rows]),
        "sample_count": len(fam_rows), "variant_count": len(bim_rows),
        "marker_mode": source_descriptor["marker_mode"],
        "source_format": "PLINK_AUTOSOMAL_PHASING_PANEL",
        "allele_orientation": source_descriptor["allele_orientation"],
        "files": {suffix: {"name": panel_prefix.with_suffix(f".{suffix}").name, "sha256": sha256_file(panel_prefix.with_suffix(f".{suffix}"))} for suffix in ("bed", "bim", "fam")},
    }
    validate_json_document(descriptor, "plink_dataset.schema.json")
    descriptor_path = output_dir / "autosomal_phasing_panel.dataset.json"
    atomic_write_json(descriptor_path, descriptor)
    summary = {
        "schema_version": "1.0.0", "method_id": "complete_genotype_autosomal_panel_v1",
        "sample_count": len(fam_rows), "variant_count": len(bim_rows),
        "chromosome_variant_counts": {str(chromosome): chromosomes[chromosome] for chromosome in range(1, 23)},
        "filters": {"sample_universe": "target_chromosome_all_qc", "autosomes_only": True, "maximum_variant_missingness": 0.0, "minimum_maf": 0.000001},
    }
    summary_path = output_dir / "autosomal_phasing_panel_summary.json"
    atomic_write_json(summary_path, summary)
    outputs = _dataset_artifacts(stage_inputs, panel_prefix, descriptor_path, descriptor)
    outputs.append(build_file_artifact(
        physical_path=summary_path,
        published_path=str(PurePosixPath(stage_inputs["published_output_dir"]) / summary_path.name),
        artifact_id="autosomal_phasing_panel_summary", artifact_type="autosomal_phasing_panel_summary",
        media_type="application/json", producer_stage=f"{stage_inputs['stage_id']}_{stage_inputs['stage_name']}",
        producer_signature=stage_inputs["signature"], assembly=descriptor["assembly"],
        sample_set_id=descriptor["sample_set_id"], variant_set_id=descriptor["variant_set_id"], sensitivity="internal",
    ))
    stage_outputs = {"schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"], "stage_name": stage_inputs["stage_name"], "signature": stage_inputs["signature"], "artifacts": outputs}
    validate_json_document(stage_outputs, "stage_outputs.schema.json")
    atomic_write_json(output_dir / "stage_outputs.json", stage_outputs)
    audit = {
        "schema_version": "1.0.0", "run_id": stage_inputs["run_id"], "stage_id": stage_inputs["stage_id"],
        "stage_name": stage_inputs["stage_name"], "method_id": "complete_genotype_autosomal_panel_v1",
        "signature": stage_inputs["signature"], "started_at": started_at, "completed_at": utc_now(),
        "duration_seconds": monotonic() - started_clock, "inputs": list(artifacts.values()), "outputs": outputs,
        "parameters": {"plink_timeout_seconds": timeout, "minimum_variants_per_autosome": minimum_per_autosome, "minimum_total_variants": minimum_total},
        "tools": [{"tool": "plink", "configured": config["tools"]["plink"], "version": plink_version}],
        "counts": {"samples": len(fam_rows), "variants": len(bim_rows), "autosomes": len(chromosomes)},
        "metrics": {"chromosome_variant_counts": summary["chromosome_variant_counts"]}, "exclusions": [], "warnings": [],
        "checks": [{"check": name, "status": "PASS"} for name in ("declared_input_integrity", "single_sample_universe", "complete_genotypes", "all_22_autosomes", "polymorphic_variants", "marker_density")],
        "known_limits": ["Le panel conserve uniquement les génotypes observés complets; aucune imputation de génotype n'est publiée."],
        "expected_visualizations": [], "manual_validation_required": False,
    }
    validate_json_document(audit, "stage_audit.schema.json")
    atomic_write_json(output_dir / "audit.json", audit)
    (output_dir / "checksums.sha256").write_text("".join(f"{item['sha256']}  {Path(str(item['path'])).name}\n" for item in outputs), encoding="utf-8")
    return 0


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="prepare-autosomal-phasing-panel")
    parser.add_argument("--stage-inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parsed = parser.parse_args(arguments)
    try:
        return execute(parsed.stage_inputs, parsed.output_dir)
    except AutosomalPanelBlockError as error:
        sys.stderr.write(f"{error}\n")
        return 4
    except AutosomalPanelExternalError as error:
        sys.stderr.write(f"{error}\n")
        return 3
    except (AutosomalPanelInputError, DocumentValidationError, OSError, ValueError, json.JSONDecodeError) as error:
        sys.stderr.write(f"{error}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
