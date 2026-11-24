"""Construit et exécute le banc D24 contre le moteur réel de l’étape 13."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
from typing import Any, Iterable

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.evaluate_founder_haplotype import load_configuration
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document, validate_tsv_table
from effet_fondateur.founder.local_ibs import infer_target_centered_ibs


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
TARGET_BP = 84_991_179
TARGET_ID = "SYNTHETIC_TARGET_D24_CHR5"
TARGET_REF = "T"
TARGET_ALT = "A"
PROFILES = ("MAXIMUM", "BASELINE", "DENSITY_50", "MISSING_1PCT")


def _workspace_path(path: Path, *, directory: bool = False) -> Path:
    resolved = path.absolute()
    folded = [part.casefold() for part in resolved.parts]
    valid = resolved.is_dir() if directory else resolved.is_file()
    if (
        not valid or resolved.is_symlink() or "calibration_workspace" not in folded
        or "dock6_real" in folded or "cache" in folded
        or any(left == "data" and right == "input" for left, right in zip(folded, folded[1:]))
    ):
        raise ValueError("d24_bench_source_outside_workspace")
    return resolved


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("d24_bench_json_invalid")
    return value


def _rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.name.endswith(".gz") else Path.open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, columns: tuple[str, ...], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _marker_positions(path: Path, start_bp: int, end_bp: int) -> set[int]:
    positions: set[int] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            if fields[0] == "chr5":
                position = int(fields[1])
                if start_bp <= position <= end_bp:
                    positions.add(position)
    return positions


def _profile_positions(
    maximum: set[int], baseline: set[int], seed: int
) -> dict[str, set[int]]:
    if TARGET_BP not in maximum:
        raise ValueError("d24_target_missing_from_maximum_marker_universe")
    baseline = (baseline & maximum) | {TARGET_BP}
    candidates = sorted(maximum - {TARGET_BP})
    ordered = sorted(
        candidates,
        key=lambda position: hashlib.sha256(f"{seed}:{position}".encode()).digest(),
    )
    density = set(ordered[: round(len(candidates) * 0.5)]) | {TARGET_BP}
    return {
        "MAXIMUM": set(maximum),
        "BASELINE": baseline,
        "DENSITY_50": density,
        "MISSING_1PCT": set(maximum),
    }


def _pedigree(root: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    path = root / "private_truth/d13/pedigree_truth.tsv.gz"
    records = _rows(path)
    by_sample = {row["sample_id"]: row for row in records}
    families = sorted({
        row["family_unit_id"] for row in records
        if not row["family_unit_id"].startswith("UNRELATED_")
    })
    if len(families) != 3:
        raise ValueError("d24_bench_requires_exactly_three_families")
    return by_sample, families


def _positive_source(root: Path, pedigree: dict[str, dict[str, str]]) -> tuple[Path, list[dict[str, Any]], list[dict[str, Any]], str]:
    manifest_path = root / "private_truth/d14/founder_simulation_manifest.json"
    manifest = _read_json(manifest_path)
    if manifest.get("manifest_state") != "FINALIZED" or manifest.get("panel", {}).get("family_count") != 3:
        raise ValueError("d24_positive_manifest_invalid")
    segment_rows = _rows(root / "private_truth/d14/founder_segments.tsv.gz")
    carrier_rows = _rows(root / "private_truth/d14/founder_carriers.tsv.gz")
    selected: list[dict[str, Any]] = []
    for family in sorted({row["family_id"] for row in segment_rows}):
        candidates = [row for row in carrier_rows if row["family_id"] == family]
        maximum_generation = max(int(pedigree[row["sample_id"]]["generation"]) for row in candidates)
        terminal = sorted(
            (row for row in candidates if int(pedigree[row["sample_id"]]["generation"]) == maximum_generation),
            key=lambda row: row["sample_id"],
        )
        row = terminal[0]
        selected.append({
            "family_id": family,
            "sample_id": row["sample_id"],
            "haplotype_index": int(row["carrier_haplotype_index"]),
        })
    segments = [{
        "family_id": row["family_id"],
        "left_boundary_bp": int(row["left_boundary_bp"]),
        "right_boundary_bp": int(row["right_boundary_bp"]),
        "left_boundary_cm": float(row["left_boundary_cm"]),
        "right_boundary_cm": float(row["right_boundary_cm"]),
        "left_censored": row["left_censored"] == "true",
        "right_censored": row["right_censored"] == "true",
    } for row in segment_rows]
    return (
        root / "private_truth/d14/founder_panel.phased.vcf.gz",
        selected,
        segments,
        sha256_file(manifest_path),
    )


def _negative_source(
    d13_root: Path,
    pedigree: dict[str, dict[str, str]],
    families: list[str],
    seed: int,
) -> tuple[Path, list[dict[str, Any]], list[dict[str, Any]], str]:
    generator = random.Random(seed)
    selected = []
    for family in families:
        candidates = [
            row for row in pedigree.values()
            if row["family_unit_id"] == family and int(row["generation"]) == 2
        ]
        candidate = sorted(candidates, key=lambda row: row["sample_id"])[generator.randrange(len(candidates))]
        selected.append({
            "family_id": family,
            "sample_id": candidate["sample_id"],
            "haplotype_index": generator.choice((1, 2)),
        })
    manifest_path = d13_root / "private_truth/d13/pedigree_manifest.json"
    digest = hashlib.sha256(json.dumps({
        "manifest": sha256_file(manifest_path), "seed": seed, "selected": selected,
    }, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return (
        d13_root / "private_truth/d13/pedigree_panel.phased.vcf.gz",
        selected,
        [],
        digest,
    )


def _map_positions(path: Path) -> dict[int, float]:
    return {int(row["POSITION_BP"]): float(row["POSITION_CM"]) for row in _rows(path)}


def _missing(scenario_id: str, sample_id: str, position: int, seed: int) -> bool:
    value = int.from_bytes(
        hashlib.sha256(f"{seed}:{scenario_id}:{sample_id}:{position}".encode()).digest()[:8],
        "big",
    )
    return value / 2**64 < 0.01


def _write_profile_vcf(
    *, source: Path, output: Path, positions: set[int], selected: list[dict[str, Any]],
    negative: bool, profile: str, scenario_id: str, missing_seed: int,
    bcftools: str,
) -> tuple[list[str], dict[int, str], dict[str, str], int]:
    selected_by_sample = {row["sample_id"]: row["haplotype_index"] for row in selected}
    raw = output.with_suffix("")
    variant_ids: dict[int, str] = {}
    target_genotypes: dict[str, str] = {}
    missing_count = 0
    sample_ids: list[str] = []
    observed: set[int] = set()
    with gzip.open(source, "rt", encoding="utf-8") as reader, raw.open("w", encoding="utf-8") as writer:
        for line in reader:
            if line.startswith("##"):
                writer.write(line)
                continue
            if line.startswith("#CHROM"):
                sample_ids = line.rstrip("\n").split("\t")[9:]
                writer.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            position = int(fields[1])
            if position not in positions:
                continue
            observed.add(position)
            if position == TARGET_BP:
                fields[2:5] = [TARGET_ID, TARGET_REF, TARGET_ALT]
            variant_ids[position] = fields[2]
            for index, sample_id in enumerate(sample_ids, 9):
                parts = fields[index].split(":")
                if position == TARGET_BP and negative:
                    haplotype = selected_by_sample.get(sample_id)
                    parts[0] = "1|0" if haplotype == 1 else "0|1" if haplotype == 2 else "0|0"
                if profile == "MISSING_1PCT" and position != TARGET_BP and _missing(
                    scenario_id, sample_id, position, missing_seed
                ):
                    parts[0] = ".|."
                    missing_count += 1
                fields[index] = ":".join(parts)
                if position == TARGET_BP:
                    target_genotypes[sample_id] = parts[0]
            writer.write("\t".join(fields) + "\n")
    if observed != positions:
        raise ValueError("d24_profile_marker_universe_not_observed")
    completed = subprocess.run(
        [bcftools, "view", "-Oz", "-o", str(output), str(raw)],
        check=False, capture_output=True, text=True,
    )
    raw.unlink()
    if completed.returncode != 0:
        raise ValueError(f"d24_bcftools_compression_failed:{completed.stderr.strip()}")
    subprocess.run([bcftools, "index", "--tbi", str(output)], check=True, capture_output=True)
    return sample_ids, variant_ids, target_genotypes, missing_count


def _auxiliary_inputs(
    *, case: Path, sample_ids: list[str], pedigree: dict[str, dict[str, str]],
    selected: list[dict[str, Any]], target_genotypes: dict[str, str],
    variant_ids: dict[int, str], positions_cm: dict[int, float],
) -> tuple[Path, Path, Path, Path]:
    selected_by_sample = {row["sample_id"]: row for row in selected}
    samples_path = case / "private_truth/samples_master.tsv"
    sample_columns = (
        "SAMPLE_ID", "SOURCE_FILE", "FID", "IID", "PID", "MID", "SEX",
        "CLINICAL_STATUS", "GROUP_LABEL", "TARGET_GENOTYPE", "TARGET_GENOTYPE_SOURCE",
        "ARRAY_BATCH", "INCLUDE_GENOMEWIDE", "INCLUDE_TARGET_CHROMOSOME", "NOTES_CODE",
    )
    sample_rows = []
    for sample_id in sample_ids:
        ped = pedigree[sample_id]
        genotype = target_genotypes[sample_id].replace("0", TARGET_REF).replace("1", TARGET_ALT).replace("|", "/")
        sample_rows.append({
            "SAMPLE_ID": sample_id, "SOURCE_FILE": "d24_synthetic.vcf.gz",
            "FID": ped["family_unit_id"], "IID": sample_id,
            "PID": ped["father_id"], "MID": ped["mother_id"], "SEX": ped["sex"],
            "CLINICAL_STATUS": "UNKNOWN", "GROUP_LABEL": "D24_SYNTHETIC",
            "TARGET_GENOTYPE": genotype, "TARGET_GENOTYPE_SOURCE": "D24_EXPLICIT",
            "ARRAY_BATCH": "D24", "INCLUDE_GENOMEWIDE": "true",
            "INCLUDE_TARGET_CHROMOSOME": "true", "NOTES_CODE": "D24_SYNTHETIC",
        })
    _write_tsv(samples_path, sample_columns, sample_rows)
    carriers_path = case / "private_truth/carrier_haplotypes.tsv"
    carrier_columns = (
        "SAMPLE_ORDER", "SAMPLE_ID", "TARGET_VARIANT_ID", "EXPLICIT_GENOTYPE",
        "PHASED_GT", "ALT_COPY_COUNT", "CARRIER_HAPLOTYPE", "PHASE_CONFIDENCE",
        "CONFIDENCE_STATUS", "RELIABILITY_STATUS", "UNRELIABLE_REASON",
    )
    carrier_rows = []
    for order, sample_id in enumerate(sample_ids, 1):
        if sample_id not in selected_by_sample:
            continue
        haplotype = selected_by_sample[sample_id]["haplotype_index"]
        carrier_rows.append({
            "SAMPLE_ORDER": order, "SAMPLE_ID": sample_id, "TARGET_VARIANT_ID": TARGET_ID,
            "EXPLICIT_GENOTYPE": "T/A", "PHASED_GT": "1|0" if haplotype == 1 else "0|1",
            "ALT_COPY_COUNT": 1, "CARRIER_HAPLOTYPE": f"H{haplotype}",
            "PHASE_CONFIDENCE": "1.0", "CONFIDENCE_STATUS": "SCORED_PASS",
            "RELIABILITY_STATUS": "PASS", "UNRELIABLE_REASON": "",
        })
    _write_tsv(carriers_path, carrier_columns, carrier_rows)
    cohorts_path = case / "private_truth/cohorts_frozen.tsv"
    cohort_columns = (
        "COHORT_ID", "SAMPLE_ID", "ROLE", "INCLUDED", "EXCLUSION_CODE",
        "DECISION_SOURCE", "INDEPENDENT_UNIT_ID", "FAMILY_REPRESENTATIVE",
    )
    controls = [
        sample_id for sample_id in sample_ids
        if pedigree[sample_id]["family_unit_id"].startswith("UNRELATED_")
        and sample_id not in selected_by_sample
    ][:20]
    cohort_rows = []
    for row in selected:
        cohort_rows.extend((
            {"COHORT_ID": "target_carriers_independent", "SAMPLE_ID": row["sample_id"], "ROLE": "CARRIER", "INCLUDED": "true", "EXCLUSION_CODE": "", "DECISION_SOURCE": "D24_TRUTH", "INDEPENDENT_UNIT_ID": row["family_id"], "FAMILY_REPRESENTATIVE": "true"},
            {"COHORT_ID": "target_chromosome_all_qc", "SAMPLE_ID": row["sample_id"], "ROLE": "CARRIER", "INCLUDED": "true", "EXCLUSION_CODE": "", "DECISION_SOURCE": "D24_TRUTH", "INDEPENDENT_UNIT_ID": row["family_id"], "FAMILY_REPRESENTATIVE": "true"},
        ))
    for sample_id in controls:
        cohort_rows.extend((
            {"COHORT_ID": "controls_unrelated", "SAMPLE_ID": sample_id, "ROLE": "CONTROL", "INCLUDED": "true", "EXCLUSION_CODE": "", "DECISION_SOURCE": "D24_BACKGROUND", "INDEPENDENT_UNIT_ID": pedigree[sample_id]["family_unit_id"], "FAMILY_REPRESENTATIVE": "false"},
            {"COHORT_ID": "target_chromosome_all_qc", "SAMPLE_ID": sample_id, "ROLE": "CONTROL", "INCLUDED": "true", "EXCLUSION_CODE": "", "DECISION_SOURCE": "D24_BACKGROUND", "INDEPENDENT_UNIT_ID": pedigree[sample_id]["family_unit_id"], "FAMILY_REPRESENTATIVE": "false"},
        ))
    _write_tsv(cohorts_path, cohort_columns, cohort_rows)
    map_path = case / "private_truth/target_genetic_map.tsv"
    map_columns = (
        "DATASET_ID", "VARIANT_ORDER", "VARIANT_ID", "CHROMOSOME", "POSITION_BP",
        "POSITION_CM", "MAP_STATUS", "LEFT_MAP_BP", "LEFT_MAP_CM", "RIGHT_MAP_BP",
        "RIGHT_MAP_CM", "LOCAL_RATE_CM_PER_MB", "IS_TARGET_VARIANT",
    )
    map_rows = []
    for order, position in enumerate(sorted(variant_ids), 1):
        cm = positions_cm[position]
        map_rows.append({
            "DATASET_ID": "D24_CHR5", "VARIANT_ORDER": order,
            "VARIANT_ID": variant_ids[position], "CHROMOSOME": 5,
            "POSITION_BP": position, "POSITION_CM": f"{cm:.12g}", "MAP_STATUS": "EXACT",
            "LEFT_MAP_BP": position, "LEFT_MAP_CM": f"{cm:.12g}",
            "RIGHT_MAP_BP": position, "RIGHT_MAP_CM": f"{cm:.12g}",
            "LOCAL_RATE_CM_PER_MB": "", "IS_TARGET_VARIANT": "true" if position == TARGET_BP else "false",
        })
    _write_tsv(map_path, map_columns, map_rows)
    for path, schema in (
        (samples_path, "samples_master.schema.json"),
        (carriers_path, "carrier_haplotypes.schema.json"),
        (cohorts_path, "cohorts_frozen.schema.json"),
        (map_path, "target_genetic_map.schema.json"),
    ):
        validate_tsv_table(path, schema)
    return carriers_path, cohorts_path, samples_path, map_path


def run_bench(
    *, d13_root: Path, d11a_root: Path, positive_roots: dict[str, Path],
    baseline_marker_vcf: Path, maximum_marker_vcf: Path, configuration_path: Path,
    workspace_root: Path, run_id: str, bcftools: str = "bcftools",
) -> dict[str, Any]:
    """Exécute les 40 cas D24 et publie un reçu agrégé sans identifiant."""
    d13_root = _workspace_path(d13_root, directory=True)
    d11a_root = _workspace_path(d11a_root, directory=True)
    baseline_marker_vcf = _workspace_path(baseline_marker_vcf)
    maximum_marker_vcf = _workspace_path(maximum_marker_vcf)
    config = load_configuration(configuration_path)
    target = new_workspace_run(workspace_root, run_id)
    d13_manifest = _read_json(d13_root / "private_truth/d13/pedigree_manifest.json")
    domain = d13_manifest["domain"]
    maximum = _marker_positions(maximum_marker_vcf, domain["start_bp"], domain["end_bp"])
    baseline = _marker_positions(baseline_marker_vcf, domain["start_bp"], domain["end_bp"])
    profiles = _profile_positions(maximum, baseline, config["profile_seed"])
    map_path = d11a_root / "private_truth/d11/genetic_map.tsv"
    positions_cm = _map_positions(map_path)
    pedigree, families = _pedigree(d13_root)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    cases_manifest = []
    try:
        for scenario in config["scenarios"]:
            scenario_id = scenario["scenario_id"]
            positive = scenario["truth_role"] == "POSITIVE"
            if positive:
                root = _workspace_path(positive_roots[scenario_id], directory=True)
                source, selected, segments, source_sha = _positive_source(root, pedigree)
            else:
                source, selected, segments, source_sha = _negative_source(
                    d13_root, pedigree, families, scenario["random_seed"]
                )
            left = max((row["left_boundary_bp"] for row in segments), default=TARGET_BP)
            right = min((row["right_boundary_bp"] for row in segments), default=TARGET_BP)
            detectable = positive and (
                sum(left <= position < TARGET_BP for position in maximum) >= config["minimum_flank_markers"]
                and sum(TARGET_BP < position <= right for position in maximum) >= config["minimum_flank_markers"]
            )
            for profile in PROFILES:
                case = temporary / "cases" / scenario_id / profile
                (case / "private_truth").mkdir(parents=True)
                stage = case / "stage13"
                vcf = case / "private_truth/candidate.phased.vcf.gz"
                sample_ids, variant_ids, target_genotypes, missing_count = _write_profile_vcf(
                    source=source, output=vcf, positions=profiles[profile], selected=selected,
                    negative=not positive, profile=profile, scenario_id=scenario_id,
                    missing_seed=config["missing_seed"], bcftools=bcftools,
                )
                carriers, cohorts, samples, genetic_map = _auxiliary_inputs(
                    case=case, sample_ids=sample_ids, pedigree=pedigree, selected=selected,
                    target_genotypes=target_genotypes, variant_ids=variant_ids,
                    positions_cm=positions_cm,
                )
                result = infer_target_centered_ibs(
                    phased_bcf_path=vcf, carrier_haplotypes_path=carriers,
                    cohorts_path=cohorts, samples_master_path=samples,
                    genetic_map_path=genetic_map, output_dir=stage,
                    bcftools_command=bcftools, timeout_seconds=300,
                    minimum_independent_carriers=config["minimum_independent_units_per_scenario"],
                    minimum_flank_markers=config["minimum_flank_markers"],
                )
                for path, schema in (
                    (result.segments_path, "founder_segments.schema.json"),
                    (result.consensus_path, "founder_consensus.schema.json"),
                    (result.sharing_path, "founder_sharing_matrix.schema.json"),
                    (result.discordances_path, "founder_discordances.schema.json"),
                ):
                    validate_tsv_table(path, schema)
                validate_json_document(_read_json(result.summary_path), "founder_analysis_summary.schema.json")
                truth = {
                    "schema_version": "1.0.0", "scenario_id": scenario_id,
                    "profile": profile, "truth_role": scenario["truth_role"],
                    "detectable": detectable, "age_generations": scenario["age_generations"],
                    "independent_unit_count": 3,
                    "target": {"variant_id": TARGET_ID, "chromosome": 5, "position_bp": TARGET_BP, "position_cm": positions_cm[TARGET_BP]},
                    "source_manifest_sha256": source_sha, "segments": segments,
                }
                validate_json_document(truth, "founder_haplotype_case.schema.json", schemas_dir=SCHEMAS)
                truth_path = case / "private_truth/founder_haplotype_case.json"
                write_json(truth_path, truth)
                cases_manifest.append({
                    "scenario_id_hash": hashlib.sha256(scenario_id.encode()).hexdigest(),
                    "profile": profile, "truth_role": scenario["truth_role"],
                    "marker_count": len(profiles[profile]), "missing_genotype_count": missing_count,
                    "stage13_status": result.status,
                    "case_sha256": hashlib.sha256("".join(sha256_file(path) for path in (
                        truth_path, result.summary_path, result.consensus_path, result.segments_path,
                    )).encode()).hexdigest(),
                })
        manifest = {
            "schema_version": "1.0.0", "run_id": run_id, "state": "CLOSED",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "configuration_sha256": sha256_file(configuration_path),
            "d13_manifest_sha256": sha256_file(d13_root / "private_truth/d13/pedigree_manifest.json"),
            "d11a_manifest_sha256": sha256_file(d11a_root / "private_truth/d11/base_haplotype_manifest.json"),
            "baseline_marker_source_sha256": sha256_file(baseline_marker_vcf),
            "maximum_marker_source_sha256": sha256_file(maximum_marker_vcf),
            "profile_marker_counts": {profile: len(positions) for profile, positions in profiles.items()},
            "case_count": len(cases_manifest), "cases": cases_manifest,
            "calibrated": False, "founder_effect_proven": False,
        }
        write_json(temporary / "bench_manifest.json", manifest)
        (temporary / "checksums.sha256").write_text(
            f"{sha256_file(temporary / 'bench_manifest.json')}  bench_manifest.json\n",
            encoding="utf-8",
        )
        temporary.rename(target)
        return manifest
    except Exception:
        shutil.rmtree(temporary)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description="Banc réel D24 de l’étape 13")
    parser.add_argument("--d13-root", type=Path, required=True)
    parser.add_argument("--d11a-root", type=Path, required=True)
    parser.add_argument("--positive-root", action="append", required=True)
    parser.add_argument("--baseline-marker-vcf", type=Path, required=True)
    parser.add_argument("--maximum-marker-vcf", type=Path, required=True)
    parser.add_argument("--configuration-path", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--bcftools", default="bcftools")
    arguments = vars(parser.parse_args())
    roots = {}
    for item in arguments.pop("positive_root"):
        scenario_id, separator, path = item.partition("=")
        if not separator or scenario_id in roots:
            raise ValueError("d24_positive_root_argument_invalid")
        roots[scenario_id] = Path(path)
    arguments["positive_roots"] = roots
    manifest = run_bench(**arguments)
    print(manifest["state"], manifest["case_count"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
