"""Banc D26 des moteurs réels de LD local et de ROH."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import tempfile
from typing import Any, Iterable

import numpy as np

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.founder_haplotype_oracle import wilson_interval
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document
from effet_fondateur.ld import publish_local_ld
from effet_fondateur.roh import publish_roh


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
PROTOCOL = REPOSITORY / "calibration/protocoles/D26_LD_ROH_PRESPECIFICATION.md"


def _workspace_path(path: Path, *, directory: bool = False) -> Path:
    resolved = path.absolute()
    folded = [part.casefold() for part in resolved.parts]
    valid = resolved.is_dir() if directory else resolved.is_file()
    if (
        not valid or resolved.is_symlink()
        or "calibration_workspace" not in folded
        or "dock6_real" in folded or "cache" in folded
        or any(left == "data" and right == "input" for left, right in zip(folded, folded[1:]))
    ):
        raise ValueError("d26_source_outside_workspace")
    return resolved


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge et contrôle le contrat D26 fermé."""
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(document, "ld_roh_bench_config.schema.json", schemas_dir=SCHEMAS)
    if document["ld"]["distance_bins_cm"][-1] != document["ld"]["max_pair_distance_cm"]:
        raise ValueError("d26_ld_bins_do_not_end_at_maximum")
    return document


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _manifest_sha256(path: Path) -> str:
    return sha256_file(_workspace_path(path))


def _source_paths(d16_root: Path, d17_root: Path) -> dict[str, Path]:
    d16 = _workspace_path(d16_root, directory=True)
    d17 = _workspace_path(d17_root, directory=True)
    return {
        "d16_manifest": _workspace_path(d16 / "roh_simulation_manifest.json"),
        "d16_dense": _workspace_path(d16 / "roh_panel.phased.vcf.gz"),
        "d16_truth": _workspace_path(d16 / "roh_truth.json"),
        "d17_manifest": _workspace_path(d17 / "measurement_error_manifest.json"),
        "d17_root": d17,
    }


def _verify_sources(config: dict[str, Any], paths: dict[str, Path]) -> None:
    if _manifest_sha256(paths["d16_manifest"]) != config["source_d16_manifest_sha256"]:
        raise ValueError("d26_d16_manifest_hash_mismatch")
    if _manifest_sha256(paths["d17_manifest"]) != config["source_d17_manifest_sha256"]:
        raise ValueError("d26_d17_manifest_hash_mismatch")
    manifest = _read_json(paths["d17_manifest"])
    artifact_by_profile = {
        row["profile_id"]: row for row in manifest["artifacts"]
        if row["role"] == "OBSERVED_MEASUREMENT_PANEL"
    }
    for profile in [*config["profiles"], *config["nonphased_aliases"]]:
        artifact = artifact_by_profile.get(profile)
        if artifact is None:
            raise ValueError(f"d26_profile_missing:{profile}")
        source = _workspace_path(paths["d17_root"] / f"profiles/{profile}/observed.vcf.gz")
        if sha256_file(source) != artifact["sha256"]:
            raise ValueError(f"d26_profile_hash_mismatch:{profile}")
    control_hash = artifact_by_profile["control"]["sha256"]
    for alias, control in config["nonphased_aliases"].items():
        if control != "control" or artifact_by_profile[alias]["sha256"] != control_hash:
            raise ValueError(f"d26_nonphased_alias_not_identical:{alias}")


def _plink_version(plink: str) -> str:
    executable = shutil.which(plink)
    if executable is None:
        raise ValueError("d26_plink_not_available")
    completed = subprocess.run([executable, "--version"], capture_output=True, text=True, check=False, timeout=10)
    version = (completed.stdout or completed.stderr).strip().splitlines()
    if completed.returncode or not version or "1.9" not in version[0]:
        raise ValueError("d26_plink19_required")
    return version[0]


def _run(command: list[str], *, timeout: int) -> None:
    completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout)
    if completed.returncode:
        tail = (completed.stderr or completed.stdout).strip().splitlines()[-1:]
        raise RuntimeError(f"d26_external_command_failed:{command[0]}:{' '.join(tail)}")


def _vcf_rows(path: Path) -> tuple[list[str], list[dict[str, Any]]]:
    samples: list[str] | None = None
    rows: list[dict[str, Any]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                rows.append({
                    "chromosome": fields[0].removeprefix("chr"),
                    "position_bp": int(fields[1]), "variant_id": fields[2],
                    "ref": fields[3], "alt": fields[4], "genotypes": fields[9:],
                })
    if samples is None or not rows or any(len(row["genotypes"]) != len(samples) for row in rows):
        raise ValueError(f"d26_malformed_vcf:{path.name}")
    return samples, rows


def _genotype(raw: str) -> str:
    value = raw.split(":", 1)[0].replace("|", "/")
    if value in {"0/0", "0/1", "1/0", "1/1"}:
        return "0/1" if value == "1/0" else value
    return "./."


def _dosage(raw: str) -> float:
    value = _genotype(raw)
    return {"0/0": 0.0, "0/1": 1.0, "1/1": 2.0}.get(value, math.nan)


def _read_map(path: Path) -> tuple[np.ndarray, np.ndarray]:
    positions: list[int] = []
    centimorgans: list[float] = []
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            positions.append(int(row["POSITION_BP"]))
            centimorgans.append(float(row["POSITION_CM"]))
    if len(positions) < 2 or positions != sorted(positions):
        raise ValueError("d26_invalid_genetic_map")
    return np.asarray(positions), np.asarray(centimorgans)


def _cm(position_bp: int, map_positions: np.ndarray, map_cm: np.ndarray) -> float:
    return float(np.interp(position_bp, map_positions, map_cm))


def _make_bed(plink: str, vcf: Path, prefix: Path, timeout: int) -> None:
    _run([
        plink, "--vcf", str(vcf), "--double-id", "--allow-no-sex",
        "--make-bed", "--out", str(prefix),
    ], timeout=timeout)


def _read_fam(path: Path) -> list[tuple[str, str]]:
    return [(row[0], row[1]) for row in (line.split() for line in path.read_text().splitlines()) if row]


def _write_keep(path: Path, pairs: Iterable[tuple[str, str]]) -> None:
    path.write_text("".join(f"{fid}\t{iid}\n" for fid, iid in pairs), encoding="utf-8")


def _variants(rows: list[dict[str, Any]], target_id: str, map_positions: np.ndarray, map_cm: np.ndarray) -> list[dict[str, Any]]:
    return [{
        "VARIANT_ORDER": index, "VARIANT_ID": row["variant_id"],
        "CHROMOSOME": row["chromosome"], "POSITION_BP": row["position_bp"],
        "POSITION_CM": _cm(row["position_bp"], map_positions, map_cm),
        "A1": row["ref"], "A2": row["alt"],
        "IS_TARGET_VARIANT": row["variant_id"] == target_id,
    } for index, row in enumerate(rows, 1)]


def _target_genotypes(
    samples: list[str], rows: list[dict[str, Any]], target_id: str,
    target_position_bp: int,
) -> dict[str, str]:
    selected = [row for row in rows if row["variant_id"] == target_id]
    if not selected:
        selected = [row for row in rows if row["position_bp"] == target_position_bp]
    if len(selected) != 1:
        raise ValueError("d26_target_missing_or_ambiguous")
    variant = selected[0]
    alleles = {"0": variant["ref"], "1": variant["alt"]}
    result: dict[str, str] = {}
    for sample, raw in zip(samples, variant["genotypes"], strict=True):
        genotype = _genotype(raw)
        result[sample] = (
            "./." if genotype == "./."
            else "/".join(alleles[index] for index in genotype.split("/"))
        )
    return result


def _file_inventory(root: Path) -> list[dict[str, Any]]:
    return [{
        "path": path.relative_to(root).as_posix(), "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    } for path in sorted(root.rglob("*")) if path.is_file() and path.name != "run_manifest.json"]


def run_bench(
    *, configuration_path: Path, d16_root: Path, d17_root: Path,
    genetic_map_path: Path, workspace_root: Path, run_id: str, plink: str = "plink",
) -> Path:
    """Exécute et scelle les appels observés sans ouvrir la vérité dense."""
    config = load_configuration(configuration_path)
    paths = _source_paths(d16_root, d17_root)
    _verify_sources(config, paths)
    map_path = _workspace_path(genetic_map_path)
    map_positions, map_cm = _read_map(map_path)
    version = _plink_version(plink)
    executable = shutil.which(plink)
    assert executable is not None
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        for profile in config["profiles"]:
            source = _workspace_path(paths["d17_root"] / f"profiles/{profile}/observed.vcf.gz")
            samples, rows = _vcf_rows(source)
            profile_root = temporary / "profiles" / profile
            with tempfile.TemporaryDirectory(prefix=f"d26_{profile}_", dir=temporary) as work_raw:
                work = Path(work_raw)
                prefix = work / "observed"
                _make_bed(executable, source, prefix, config["ld"]["plink_timeout_seconds"])
                fam_pairs = _read_fam(prefix.with_suffix(".fam"))
                if [iid for _, iid in fam_pairs] != samples:
                    raise ValueError(f"d26_sample_order_mismatch:{profile}")
                keep_all = work / "all.keep"
                keep_empty = work / "empty.keep"
                _write_keep(keep_all, fam_pairs)
                _write_keep(keep_empty, [])
                publish_local_ld(
                    plink_executable=executable, dataset_prefix=prefix,
                    variants=_variants(rows, config["target"]["variant_id"], map_positions, map_cm),
                    cohort_keep_paths={
                        "controls_unrelated": keep_all,
                        "target_carriers_independent": keep_empty,
                        "family_noncarriers": keep_empty,
                    },
                    cohort_sample_counts={"controls_unrelated": len(samples), "target_carriers_independent": 0, "family_noncarriers": 0},
                    cohort_independent={"controls_unrelated": True, "target_carriers_independent": True, "family_noncarriers": True},
                    output_dir=profile_root / "ld",
                    minimum_samples=config["ld"]["minimum_samples"],
                    primary_samples=config["ld"]["minimum_primary_samples"],
                    minimum_called_samples_per_pair=config["ld"]["minimum_called_samples_per_pair"],
                    minimum_variant_maf=config["ld"]["minimum_variant_maf"],
                    minimum_pairs_per_bin=config["ld"]["minimum_pairs_per_bin"],
                    max_pair_distance_bp=config["ld"]["max_pair_distance_bp"],
                    max_pair_distance_cm=config["ld"]["max_pair_distance_cm"],
                    distance_bins_cm=config["ld"]["distance_bins_cm"],
                    max_pair_count=config["ld"]["max_pair_count"],
                    timeout_seconds=config["ld"]["plink_timeout_seconds"],
                )
                sample_by_pair = {pair: pair[1] for pair in fam_pairs}
                target_genotypes = _target_genotypes(
                    samples, rows, config["target"]["variant_id"],
                    config["target"]["position_bp"],
                )
                target_pairs = [pair for pair in fam_pairs if target_genotypes[pair[1]] != "./."]
                if not target_pairs:
                    raise ValueError(f"d26_no_called_target_genotypes:{profile}")
                target_keep = work / "target.keep"
                target_prefix = work / "target_qc"
                _write_keep(target_keep, target_pairs)
                _run([
                    executable, "--bfile", str(prefix), "--keep", str(target_keep),
                    "--make-bed", "--out", str(target_prefix),
                ], timeout=config["roh"]["plink_timeout_seconds"])
                publish_roh(
                    plink_executable=executable,
                    genomewide_datasets={
                        "controls_unrelated": {
                            "dataset_id": f"d26_{profile}_regional_context",
                            "prefix": prefix, "fam_pairs": fam_pairs,
                        }
                    },
                    target_dataset={"dataset_id": f"d26_{profile}", "prefix": target_prefix, "fam_pairs": target_pairs},
                    common_variant_ids=[], common_autosome_count=0,
                    target_variant_count=len(rows), sample_by_pair=sample_by_pair,
                    target_genotypes=target_genotypes,
                    roles_by_sample={sample: {"controls_unrelated"} for sample in samples},
                    target_variant_id=config["target"]["variant_id"],
                    target_chromosome=config["target"]["chromosome"],
                    target_bp=config["target"]["position_bp"],
                    output_dir=profile_root / "roh",
                    parameters={key: value for key, value in config["roh"].items() if key != "minimum_reference_overlap_fraction"},
                )
            write_json(profile_root / "source.json", {
                "profile_id": profile, "observed_vcf_sha256": sha256_file(source),
                "sample_count": len(samples), "marker_count": len(rows),
            })
        inventory = _file_inventory(temporary)
        manifest = {
            "schema_version": "1.0.0", "method_id": "d26_observed_ld_roh_calls_v1",
            "state": "CLOSED_BEFORE_DENSE_TRUTH", "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "source_d16_manifest_sha256": config["source_d16_manifest_sha256"],
            "source_d17_manifest_sha256": config["source_d17_manifest_sha256"],
            "genetic_map_sha256": sha256_file(map_path), "plink_version": version,
            "profiles": config["profiles"], "nonphased_aliases_verified": True,
            "dense_truth_opened": False, "artifacts": inventory,
            "technical_status": "PASS", "calibrated": False, "founder_effect_proven": False,
        }
        write_json(temporary / "run_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _exact_r2(first: np.ndarray, second: np.ndarray) -> float | None:
    called = np.isfinite(first) & np.isfinite(second)
    x, y = first[called], second[called]
    if len(x) < 2 or np.var(x) == 0 or np.var(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def _ld_reference(path: Path) -> dict[str, np.ndarray]:
    _, rows = _vcf_rows(path)
    return {row["variant_id"]: np.asarray([_dosage(value) for value in row["genotypes"]]) for row in rows}


def _percentile(values: list[float], percentile: float) -> float | None:
    return None if not values else float(np.percentile(values, percentile))


def _ld_result(profile: str, pairs_path: Path, reference: dict[str, np.ndarray]) -> dict[str, Any]:
    evaluated = [row for row in _read_tsv(pairs_path) if row["COHORT_ID"] == "controls_unrelated" and row["PAIR_STATUS"] == "EVALUATED"]
    errors: list[float] = []
    missing_reference = 0
    for row in evaluated:
        first, second = reference.get(row["VARIANT_A"]), reference.get(row["VARIANT_B"])
        truth = None if first is None or second is None else _exact_r2(first, second)
        if truth is None:
            missing_reference += 1
        else:
            errors.append(abs(float(row["R2_GENOTYPE"]) - truth))
    return {
        "profile_id": profile, "evaluated_pair_count": len(evaluated),
        "reference_pair_count": len(errors), "missing_reference_pair_count": missing_reference,
        "r_squared_mean_absolute_error": None if not errors else float(np.mean(errors)),
        "r_squared_median_absolute_error": None if not errors else float(np.median(errors)),
        "r_squared_p95_absolute_error": _percentile(errors, 95),
        "r_squared_maximum_absolute_error": None if not errors else max(errors),
    }


def _dense_reference(
    *, dense_vcf: Path, output_dir: Path, config: dict[str, Any],
    map_positions: np.ndarray, map_cm: np.ndarray, plink: str,
) -> Path:
    samples, rows = _vcf_rows(dense_vcf)
    with tempfile.TemporaryDirectory(prefix="d26_dense_", dir=output_dir.parent) as work_raw:
        work = Path(work_raw)
        prefix = work / "dense"
        _make_bed(plink, dense_vcf, prefix, config["roh"]["plink_timeout_seconds"])
        fam_pairs = _read_fam(prefix.with_suffix(".fam"))
        publish_roh(
            plink_executable=plink,
            genomewide_datasets={
                "controls_unrelated": {
                    "dataset_id": "d26_dense_regional_context",
                    "prefix": prefix, "fam_pairs": fam_pairs,
                }
            },
            target_dataset={"dataset_id": "d26_dense_reference", "prefix": prefix, "fam_pairs": fam_pairs},
            common_variant_ids=[], common_autosome_count=0, target_variant_count=len(rows),
            sample_by_pair={pair: pair[1] for pair in fam_pairs},
            target_genotypes=_target_genotypes(
                samples, rows, config["target"]["variant_id"],
                config["target"]["position_bp"],
            ),
            roles_by_sample={sample: {"controls_unrelated"} for sample in samples},
            target_variant_id=config["target"]["variant_id"], target_chromosome=config["target"]["chromosome"],
            target_bp=config["target"]["position_bp"], output_dir=output_dir,
            parameters={key: value for key, value in config["roh"].items() if key != "minimum_reference_overlap_fraction"},
        )
    return output_dir / "roh_segments.tsv"


def _segments(path: Path) -> list[dict[str, Any]]:
    return [{
        "sample_id": row["SAMPLE_ID"], "start_bp": int(row["POS1_BP"]),
        "end_bp": int(row["POS2_BP"]), "length_kb": float(row["LENGTH_KB"]),
    } for row in _read_tsv(path) if row["SCOPE"] == "TARGET_CHROMOSOME"]


def _intersection(first: dict[str, Any], second: dict[str, Any]) -> int:
    return max(0, min(first["end_bp"], second["end_bp"]) - max(first["start_bp"], second["start_bp"]) + 1)


def _match_segments(observed: list[dict[str, Any]], reference: list[dict[str, Any]], fraction: float) -> list[tuple[int, int]]:
    candidates: list[tuple[int, int, int]] = []
    for observed_index, candidate in enumerate(observed):
        for reference_index, truth in enumerate(reference):
            if candidate["sample_id"] != truth["sample_id"]:
                continue
            overlap = _intersection(candidate, truth)
            reference_length = truth["end_bp"] - truth["start_bp"] + 1
            if overlap / reference_length >= fraction:
                candidates.append((overlap, observed_index, reference_index))
    matches: list[tuple[int, int]] = []
    used_observed: set[int] = set()
    used_reference: set[int] = set()
    for _, observed_index, reference_index in sorted(candidates, reverse=True):
        if observed_index not in used_observed and reference_index not in used_reference:
            used_observed.add(observed_index)
            used_reference.add(reference_index)
            matches.append((observed_index, reference_index))
    return matches


def _roh_result(
    profile: str, observed: list[dict[str, Any]], reference: list[dict[str, Any]],
    samples: list[str], truth_rows: list[dict[str, Any]], overlap_fraction: float,
    map_positions: np.ndarray, map_cm: np.ndarray,
) -> dict[str, Any]:
    matches = _match_segments(observed, reference, overlap_fraction)
    reference_positive_samples = {row["sample_id"] for row in reference}
    observed_positive_samples = {row["sample_id"] for row in observed}
    negative_samples = set(samples) - reference_positive_samples
    false_positive_samples = negative_samples & observed_positive_samples
    length_errors: list[float] = []
    for observed_index, reference_index in matches:
        candidate, truth = observed[observed_index], reference[reference_index]
        observed_cm = _cm(candidate["end_bp"], map_positions, map_cm) - _cm(candidate["start_bp"], map_positions, map_cm)
        truth_cm = _cm(truth["end_bp"], map_positions, map_cm) - _cm(truth["start_bp"], map_positions, map_cm)
        length_errors.append(abs(observed_cm - truth_cm))
    evaluated_samples = set(samples)
    injected = [
        row for row in truth_rows
        if row["is_positive"] and row["sample_id"] in evaluated_samples
    ]
    injection_detected = 0
    for truth in injected:
        truth_segment = {"start_bp": truth["start_bp"], "end_bp": truth["end_bp"]}
        truth_length = truth["end_bp"] - truth["start_bp"] + 1
        if any(row["sample_id"] == truth["sample_id"] and _intersection(row, truth_segment) / truth_length >= overlap_fraction for row in observed):
            injection_detected += 1
    sensitivity_low, sensitivity_high = wilson_interval(len(matches), len(reference)) if reference else (None, None)
    fp_low, fp_high = wilson_interval(len(false_positive_samples), len(negative_samples)) if negative_samples else (None, None)
    return {
        "profile_id": profile, "observed_segment_count": len(observed),
        "dense_reference_segment_count": len(reference), "matched_segment_count": len(matches),
        "roh_detection_sensitivity": None if not reference else len(matches) / len(reference),
        "sensitivity_ci95_low": sensitivity_low, "sensitivity_ci95_high": sensitivity_high,
        "reference_negative_sample_count": len(negative_samples),
        "false_positive_sample_count": len(false_positive_samples),
        "roh_false_positive_rate": None if not negative_samples else len(false_positive_samples) / len(negative_samples),
        "false_positive_rate_ci95_low": fp_low, "false_positive_rate_ci95_high": fp_high,
        "unmatched_observed_segment_fraction": None if not observed else (len(observed) - len(matches)) / len(observed),
        "roh_length_mean_absolute_error_cm": None if not length_errors else float(np.mean(length_errors)),
        "roh_length_median_absolute_error_cm": None if not length_errors else float(np.median(length_errors)),
        "injected_positive_count": len(injected), "injected_positive_detected_count": injection_detected,
        "injected_positive_sensitivity": None if not injected else injection_detected / len(injected),
    }


def evaluate_bench(
    *, configuration_path: Path, run_root: Path, d16_root: Path, d17_root: Path,
    genetic_map_path: Path, workspace_root: Path, run_id: str, plink: str = "plink",
) -> Path:
    """Ouvre la vérité après scellement et produit CAL-REP-15/16 en JSON."""
    config = load_configuration(configuration_path)
    paths = _source_paths(d16_root, d17_root)
    _verify_sources(config, paths)
    run_root = _workspace_path(run_root, directory=True)
    run_manifest_path = _workspace_path(run_root / "run_manifest.json")
    run_manifest = _read_json(run_manifest_path)
    if run_manifest["configuration_sha256"] != sha256_file(configuration_path) or run_manifest["protocol_sha256"] != sha256_file(PROTOCOL):
        raise ValueError("d26_run_binding_invalid")
    for artifact in run_manifest["artifacts"]:
        path = _workspace_path(run_root / artifact["path"])
        if sha256_file(path) != artifact["sha256"]:
            raise ValueError(f"d26_run_artifact_modified:{artifact['path']}")
    map_path = _workspace_path(genetic_map_path)
    if sha256_file(map_path) != run_manifest["genetic_map_sha256"]:
        raise ValueError("d26_map_binding_invalid")
    map_positions, map_cm = _read_map(map_path)
    executable = shutil.which(plink)
    if executable is None:
        raise ValueError("d26_plink_not_available")
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        control_vcf = _workspace_path(paths["d17_root"] / "profiles/control/observed.vcf.gz")
        samples, _ = _vcf_rows(control_vcf)
        ld_reference = _ld_reference(control_vcf)
        ld_results = [
            _ld_result(profile, run_root / f"profiles/{profile}/ld/local_ld_pairs.tsv", ld_reference)
            for profile in config["profiles"]
        ]
        dense_output = temporary / "dense_roh_reference"
        dense_segments_path = _dense_reference(
            dense_vcf=paths["d16_dense"], output_dir=dense_output, config=config,
            map_positions=map_positions, map_cm=map_cm, plink=executable,
        )
        reference_segments = _segments(dense_segments_path)
        truth_rows = _read_json(paths["d16_truth"])
        roh_results = []
        for profile in config["profiles"]:
            profile_samples = [
                row["SAMPLE_ID"] for row in _read_tsv(
                    run_root / f"profiles/{profile}/roh/target_in_roh.tsv"
                )
            ]
            profile_sample_set = set(profile_samples)
            profile_reference = [
                row for row in reference_segments
                if row["sample_id"] in profile_sample_set
            ]
            roh_results.append(_roh_result(
                profile,
                _segments(run_root / f"profiles/{profile}/roh/roh_segments.tsv"),
                profile_reference, profile_samples, truth_rows,
                config["roh"]["minimum_reference_overlap_fraction"],
                map_positions, map_cm,
            ))
        common = {
            "schema_version": "1.0.0", "evaluation_state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "run_manifest_sha256": sha256_file(run_manifest_path),
            "source_d16_manifest_sha256": config["source_d16_manifest_sha256"],
            "source_d17_manifest_sha256": config["source_d17_manifest_sha256"],
            "profiles": config["profiles"], "paired_profiles": True,
            "performance_thresholds_sealed": False,
            "technical_status": "PASS", "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
            "production_promotion": False, "calibrated": False, "founder_effect_proven": False,
        }
        ld_document = {
            **common, "report_id": "CAL-REP-15", "method_id": "d26_local_ld_evaluation_v1",
            "r_squared_definition": "GENOTYPE_DOSAGE_PEARSON_CORRELATION_SQUARED",
            "results_by_profile": ld_results,
            "criteria": [
                {"criterion_id": "CAL-LD-001", "metric": "r_squared_mean_absolute_error", "value": ld_results[0]["r_squared_mean_absolute_error"], "threshold": None, "status": "REVIEW_REQUIRED"},
                {"criterion_id": "CAL-LD-002", "metric": "target_chromosome_in_background_count", "value": 0, "threshold": 0, "status": "PASS"},
            ],
            "limits": ["Le contrôle intact est une référence de dosage, pas une vérité populationnelle externe.", "Les profils sont appariés et ne sont pas six réplications indépendantes.", "Le panel régional chr5 ne mesure pas le LD génomique ni DOCK6."],
        }
        roh_document = {
            **common, "report_id": "CAL-REP-16", "method_id": "d26_roh_evaluation_v1",
            "reference_role": "DENSE_D16_PLINK_OPERATIONAL_REFERENCE",
            "dense_reference_segment_count": len(reference_segments),
            "results_by_profile": roh_results,
            "criteria": [
                {"criterion_id": "CAL-ROH-001", "metric": "roh_detection_sensitivity", "value": roh_results[0]["roh_detection_sensitivity"], "threshold": None, "status": "REVIEW_REQUIRED"},
                {"criterion_id": "CAL-ROH-002", "metric": "roh_false_positive_rate", "value": roh_results[0]["roh_false_positive_rate"], "threshold": None, "status": "REVIEW_REQUIRED"},
                {"criterion_id": "CAL-ROH-003", "metric": "roh_length_mean_absolute_error", "value": roh_results[0]["roh_length_mean_absolute_error_cm"], "threshold": None, "status": "REVIEW_REQUIRED"},
                {"criterion_id": "CAL-ROH-004", "metric": "roh_used_as_direct_gamma_length_count", "value": 0, "threshold": 0, "status": "PASS"},
            ],
            "limits": ["La référence dense réutilise l'algorithme PLINK et n'est pas une vérité biologique indépendante.", "Les trois bornes injectées ne sont pas des limites maximales d'appel.", "Le banc régional ne qualifie ni F_ROH génomique ni d'autres classes de longueur."],
        }
        validate_json_document(ld_document, "ld_roh_evaluation.schema.json", schemas_dir=SCHEMAS)
        validate_json_document(roh_document, "ld_roh_evaluation.schema.json", schemas_dir=SCHEMAS)
        write_json(temporary / "CAL-REP-15.json", ld_document)
        write_json(temporary / "CAL-REP-16.json", roh_document)
        inventory = _file_inventory(temporary)
        write_json(temporary / "evaluation_manifest.json", {
            "schema_version": "1.0.0", "state": "CLOSED", "artifacts": inventory,
            "technical_status": "PASS", "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
            "calibrated": False, "founder_effect_proven": False,
        })
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "evaluate"):
        command = subparsers.add_parser(name)
        command.add_argument("--configuration", type=Path, required=True)
        command.add_argument("--d16-root", type=Path, required=True)
        command.add_argument("--d17-root", type=Path, required=True)
        command.add_argument("--genetic-map", type=Path, required=True)
        command.add_argument("--workspace-root", type=Path, required=True)
        command.add_argument("--run-id", required=True)
        command.add_argument("--plink", default="plink")
        if name == "evaluate":
            command.add_argument("--run-root", type=Path, required=True)
    return parser


def main(arguments: list[str] | None = None) -> int:
    parsed = _parser().parse_args(arguments)
    common = dict(
        configuration_path=parsed.configuration, d16_root=parsed.d16_root,
        d17_root=parsed.d17_root, genetic_map_path=parsed.genetic_map,
        workspace_root=parsed.workspace_root, run_id=parsed.run_id, plink=parsed.plink,
    )
    result = run_bench(**common) if parsed.command == "run" else evaluate_bench(run_root=parsed.run_root, **common)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
