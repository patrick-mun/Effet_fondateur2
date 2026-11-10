"""Banc D28 du cœur d'enrichissement haplotypique de l'étape 16B."""

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
import tempfile
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document
from effet_fondateur.founder_enrichment import (
    HaplotypeProfile,
    Marker,
    evaluate_exact_sharing,
    sample_external_null,
    summarize_null,
)
from effet_fondateur.founder_enrichment.observed import (
    validate_independent_family_profiles,
)
from effet_fondateur.stages.evaluate_founder_haplotype_enrichment import (
    _classification,
)


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
PROTOCOL = REPOSITORY / "calibration/protocoles/D28_ENRICHMENT_PRESPECIFICATION.md"
D24_CONFIGURATION = REPOSITORY / "calibration/configurations/d24_founder_haplotype_bench_v1.json"
PROFILES = ("MAXIMUM", "BASELINE", "DENSITY_50", "MISSING_1PCT")


def _workspace_path(path: Path, *, directory: bool = False) -> Path:
    resolved = path.absolute()
    folded = [part.casefold() for part in resolved.parts]
    valid = resolved.is_dir() if directory else resolved.is_file()
    if (
        not valid
        or resolved.is_symlink()
        or "calibration_workspace" not in folded
        or "dock6_real" in folded
        or "cache" in folded
        or any(left == "data" and right == "input" for left, right in zip(folded, folded[1:]))
    ):
        raise ValueError("d28_source_outside_workspace")
    return resolved


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge le contrat D28 et refuse toute classification de production."""
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(document, "enrichment_bench_config.schema.json", schemas_dir=SCHEMAS)
    if (
        document["empirical_classification_threshold"] is not None
        or document["performance_thresholds_sealed"]
        or document["production_promotion"]
    ):
        raise ValueError("d28_unexpected_scientific_promotion")
    if document["null_max_attempts"] < document["null_evaluable_draws"]:
        raise ValueError("d28_null_attempt_count_invalid")
    return document


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _case_hash(case_name: str) -> str:
    return hashlib.sha256(case_name.encode("utf-8")).hexdigest()


def _markers(path: Path) -> tuple[Marker, ...]:
    rows = _rows(path)
    markers = tuple(
        Marker(
            row["VARIANT_ID"], int(row["POSITION_BP"]), float(row["POSITION_CM"]),
            row["IS_TARGET_VARIANT"].casefold() == "true",
        )
        for row in rows
    )
    if len(markers) < 5 or sum(marker.is_target for marker in markers) != 1:
        raise ValueError("d28_marker_map_invalid")
    return markers


def _phased_alleles(
    path: Path, markers: tuple[Marker, ...],
) -> dict[str, tuple[tuple[str | None, ...], tuple[str | None, ...]]]:
    expected = {marker.variant_id for marker in markers}
    samples: tuple[str, ...] | None = None
    by_variant: dict[str, tuple[tuple[str | None, str | None], ...]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = tuple(line.rstrip("\n").split("\t")[9:])
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] not in expected:
                continue
            if samples is None or len(fields[9:]) != len(samples):
                raise ValueError("d28_vcf_shape_invalid")
            pairs: list[tuple[str | None, str | None]] = []
            for sample_field in fields[9:]:
                genotype = sample_field.split(":", 1)[0]
                alleles = genotype.split("|")
                if len(alleles) != 2 or any(value not in {"0", "1", "."} for value in alleles):
                    pairs.append((None, None))
                else:
                    pairs.append(tuple(None if value == "." else value for value in alleles))  # type: ignore[arg-type]
            if fields[2] in by_variant:
                raise ValueError("d28_duplicate_variant")
            by_variant[fields[2]] = tuple(pairs)
    if samples is None or set(by_variant) != expected:
        raise ValueError("d28_vcf_marker_set_mismatch")
    result: dict[str, tuple[tuple[str | None, ...], tuple[str | None, ...]]] = {}
    for sample_index, sample_id in enumerate(samples):
        ordered = [by_variant[marker.variant_id][sample_index] for marker in markers]
        result[sample_id] = (
            tuple(pair[0] for pair in ordered), tuple(pair[1] for pair in ordered)
        )
    return result


def _profile(
    alleles: dict[str, tuple[tuple[str | None, ...], tuple[str | None, ...]]],
    sample_id: str,
    haplotype_id: str,
    family_id: str | None = None,
) -> HaplotypeProfile:
    if sample_id not in alleles or haplotype_id not in {"H1", "H2"}:
        raise ValueError("d28_haplotype_selection_invalid")
    copy_index = 0 if haplotype_id == "H1" else 1
    return HaplotypeProfile(sample_id, haplotype_id, alleles[sample_id][copy_index], family_id)


def _source_digest(paths: list[Path]) -> str:
    hashes = [sha256_file(path) for path in paths]
    return hashlib.sha256(json.dumps(hashes, separators=(",", ":")).encode()).hexdigest()


def _case_call(
    profile_root: Path, config: dict[str, Any], *, seed: int,
) -> tuple[dict[str, Any], str]:
    private = profile_root / "private_truth"
    stage13 = profile_root / "stage13"
    paths = [
        private / "target_genetic_map.tsv", private / "candidate.phased.vcf.gz",
        private / "carrier_haplotypes.tsv", private / "cohorts_frozen.tsv",
        stage13 / "founder_segments.tsv", stage13 / "founder_consensus.tsv",
        stage13 / "founder_analysis_summary.json",
    ]
    if any(not path.is_file() or path.is_symlink() for path in paths):
        raise ValueError("d28_case_input_missing")
    summary = json.loads(paths[-1].read_text(encoding="utf-8"))
    base = {
        "stage13_status": summary["status"], "evaluation_status": "NOT_EVALUATED",
        "non_evaluable_reason": None, "observed_total_shared_cm": None,
        "attempted_null_draws": 0, "evaluable_null_draws": 0,
        "non_evaluable_null_draws": 0, "exceedance_count": 0,
        "empirical_probability": None, "monte_carlo_interval_low": None,
        "monte_carlo_interval_high": None, "benchmark_positive_call": False,
        "production_status": "NOT_EVALUATED",
    }
    if summary["status"] != "SUPPORTED_IBS_CANDIDATE":
        base["non_evaluable_reason"] = "UPSTREAM_STAGE13_NOT_SUPPORTED"
        return base, _source_digest(paths)

    markers = _markers(paths[0])
    alleles = _phased_alleles(paths[1], markers)
    selected = [row for row in _rows(paths[4]) if row["SEGMENT_STATUS"] == "INCLUDED"]
    representatives = tuple(
        _profile(alleles, row["SAMPLE_ID"], row["CARRIER_HAPLOTYPE_ID"], row["FAMILY_ID"])
        for row in selected
    )
    validate_independent_family_profiles(representatives, config["minimum_independent_units"])
    observed = evaluate_exact_sharing(
        markers, representatives, minimum_flank_markers=config["minimum_flank_markers"]
    )
    if observed.evaluation_status != "EVALUATED" or observed.total_shared_cm is None:
        base["non_evaluable_reason"] = observed.non_evaluable_reason or "OBSERVED_NOT_EVALUABLE"
        return base, _source_digest(paths)

    consensus = _rows(paths[5])
    if len(consensus) != 1 or abs(observed.total_shared_cm - (
        float(consensus[0]["LEFT_LENGTH_CM"]) + float(consensus[0]["RIGHT_LENGTH_CM"])
    )) > 1e-8:
        raise ValueError("d28_stage13_statistic_not_reproduced")

    control_ids = {
        row["SAMPLE_ID"] for row in _rows(paths[3])
        if row["COHORT_ID"] == "controls_unrelated" and row["INCLUDED"].casefold() == "true"
    }
    background = tuple(
        _profile(alleles, sample_id, haplotype_id)
        for sample_id in sorted(control_ids) for haplotype_id in ("H1", "H2")
    )
    draws = sample_external_null(
        markers, background, unit_count=len(representatives),
        evaluable_draws=config["null_evaluable_draws"],
        max_attempts=config["null_max_attempts"], random_seed=seed,
        minimum_flank_markers=config["minimum_flank_markers"], stratum="D24_SYNTHETIC_BACKGROUND",
    )
    null = summarize_null(draws, observed.total_shared_cm, exhaustive=False)
    if null.evaluable_draws < config["null_evaluable_draws"] or null.empirical_probability is None:
        base.update({
            "non_evaluable_reason": "INSUFFICIENT_EVALUABLE_NULL_DRAWS",
            "observed_total_shared_cm": observed.total_shared_cm,
            "attempted_null_draws": null.attempted_draws,
            "evaluable_null_draws": null.evaluable_draws,
            "non_evaluable_null_draws": null.non_evaluable_draws,
        })
        return base, _source_digest(paths)

    production_status = _classification(
        null.empirical_probability, config["empirical_classification_threshold"], True, False
    )
    base.update({
        "evaluation_status": "EVALUATED", "observed_total_shared_cm": observed.total_shared_cm,
        "attempted_null_draws": null.attempted_draws,
        "evaluable_null_draws": null.evaluable_draws,
        "non_evaluable_null_draws": null.non_evaluable_draws,
        "exceedance_count": null.exceedance_count,
        "empirical_probability": null.empirical_probability,
        "monte_carlo_interval_low": null.interval_low,
        "monte_carlo_interval_high": null.interval_high,
        "benchmark_positive_call": null.empirical_probability <= config["benchmark_decision_alpha"],
        "production_status": production_status,
    })
    return base, _source_digest(paths)


def run_calls(
    *, configuration_path: Path, d24_bench_root: Path, workspace_root: Path, run_id: str,
) -> Path:
    """Calcule et scelle les probabilités avant toute lecture de la vérité D24."""
    config = load_configuration(configuration_path)
    bench = _workspace_path(d24_bench_root, directory=True)
    manifest_path = _workspace_path(bench / "bench_manifest.json")
    if sha256_file(manifest_path) != config["sources"]["d24_bench_manifest_sha256"]:
        raise ValueError("d28_d24_manifest_hash_mismatch")
    if sha256_file(D24_CONFIGURATION) != config["sources"]["d24_configuration_sha256"]:
        raise ValueError("d28_d24_configuration_hash_mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {(row["scenario_id_hash"], row["profile"]) for row in manifest["cases"]}
    calls: list[dict[str, Any]] = []
    source_digests: list[str] = []
    case_roots = sorted((bench / "cases").iterdir(), key=lambda path: path.name)
    for case_index, case_root in enumerate(case_roots):
        scenario_hash = _case_hash(case_root.name)
        for profile_index, profile in enumerate(config["profiles"]):
            if (scenario_hash, profile) not in expected:
                raise ValueError("d28_case_not_declared_in_d24_manifest")
            result, source_digest = _case_call(
                case_root / profile, config,
                seed=config["random_seed"] + case_index * len(PROFILES) + profile_index,
            )
            calls.append({"scenario_id_hash": scenario_hash, "profile": profile, **result})
            source_digests.append(source_digest)
    if len(calls) != 40 or len(expected) != 40:
        raise ValueError("d28_case_count_invalid")

    destination = new_workspace_run(workspace_root.absolute(), run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        calls_path = temporary / "sealed_calls.json"
        write_json(calls_path, {"schema_version": "1.0.0", "calls": calls})
        call_manifest = {
            "schema_version": "1.0.0", "run_id": run_id,
            "method_id": config["method_id"], "call_state": "SEALED_BEFORE_TRUTH",
            "truth_opened": False, "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "d24_bench_manifest_sha256": sha256_file(manifest_path),
            "case_count": len(calls), "independent_scenario_count": len(case_roots),
            "source_digest_set_sha256": hashlib.sha256("".join(sorted(source_digests)).encode()).hexdigest(),
            "artifact": {"path": calls_path.name, "sha256": sha256_file(calls_path), "size_bytes": calls_path.stat().st_size},
            "checks": ["NO_TRUTH_FILE_READ_DURING_CALLS", "TARGET_EXCLUDED_BY_PRODUCTION_PRIMITIVE", "PRODUCTION_THRESHOLD_NULL", "NO_CACHE_OR_REAL_DATA"],
        }
        write_json(temporary / "calls_manifest.json", call_manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _wilson(successes: int, trials: int) -> tuple[float | None, float | None]:
    if trials == 0:
        return None, None
    z = 1.959963984540054
    proportion = successes / trials
    denominator = 1 + z * z / trials
    centre = (proportion + z * z / (2 * trials)) / denominator
    half = z * math.sqrt(proportion * (1 - proportion) / trials + z * z / (4 * trials * trials)) / denominator
    return max(0.0, centre - half), min(1.0, centre + half)


def _profile_result(rows: list[dict[str, Any]], profile: str) -> dict[str, Any]:
    selected = [row for row in rows if row["profile"] == profile]
    positives = [row for row in selected if row["truth_role"] == "POSITIVE"]
    nulls = [row for row in selected if row["truth_role"] == "NEGATIVE"]
    positive_calls = sum(bool(row["benchmark_positive_call"]) for row in positives)
    null_calls = sum(bool(row["benchmark_positive_call"]) for row in nulls)
    power_ci = _wilson(positive_calls, len(positives))
    fpr_ci = _wilson(null_calls, len(nulls))
    positive_probabilities = [float(row["empirical_probability"]) for row in positives if row["empirical_probability"] is not None]
    null_probabilities = [float(row["empirical_probability"]) for row in nulls if row["empirical_probability"] is not None]
    return {
        "profile": profile, "positive_total": len(positives), "null_total": len(nulls),
        "evaluable_positive": sum(row["evaluation_status"] == "EVALUATED" for row in positives),
        "evaluable_null": sum(row["evaluation_status"] == "EVALUATED" for row in nulls),
        "positive_calls": positive_calls, "null_calls": null_calls,
        "power": positive_calls / len(positives), "power_ci95_low": power_ci[0], "power_ci95_high": power_ci[1],
        "null_false_positive_rate": null_calls / len(nulls), "null_fpr_ci95_low": fpr_ci[0], "null_fpr_ci95_high": fpr_ci[1],
        "median_positive_probability": statistics.median(positive_probabilities) if positive_probabilities else None,
        "median_null_probability": statistics.median(null_probabilities) if null_probabilities else None,
        "non_evaluable_count": sum(row["evaluation_status"] != "EVALUATED" for row in selected),
    }


def evaluate(
    *, configuration_path: Path, d24_bench_root: Path, calls_root: Path, output_root: Path,
) -> Path:
    """Ouvre la vérité D24 après scellement et publie les agrégats D28."""
    config = load_configuration(configuration_path)
    bench = _workspace_path(d24_bench_root, directory=True)
    calls_dir = _workspace_path(calls_root, directory=True)
    manifest_path = _workspace_path(calls_dir / "calls_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    calls_path = _workspace_path(calls_dir / manifest["artifact"]["path"])
    if (
        manifest["call_state"] != "SEALED_BEFORE_TRUTH" or manifest["truth_opened"]
        or manifest["configuration_sha256"] != sha256_file(configuration_path)
        or manifest["protocol_sha256"] != sha256_file(PROTOCOL)
        or manifest["artifact"]["sha256"] != sha256_file(calls_path)
    ):
        raise ValueError("d28_calls_not_validly_sealed")
    calls = json.loads(calls_path.read_text(encoding="utf-8"))["calls"]
    by_key = {(row["scenario_id_hash"], row["profile"]): row for row in calls}
    if len(by_key) != len(calls):
        raise ValueError("d28_duplicate_sealed_call")

    combined: list[dict[str, Any]] = []
    truth_roles: dict[str, str] = {}
    for case_root in sorted((bench / "cases").iterdir(), key=lambda path: path.name):
        scenario_hash = _case_hash(case_root.name)
        for profile in config["profiles"]:
            truth_path = _workspace_path(case_root / profile / "private_truth/founder_haplotype_case.json")
            truth = json.loads(truth_path.read_text(encoding="utf-8"))
            if _case_hash(truth["scenario_id"]) != scenario_hash or truth["profile"] != profile:
                raise ValueError("d28_truth_identity_mismatch")
            truth_roles.setdefault(scenario_hash, truth["truth_role"])
            if truth_roles[scenario_hash] != truth["truth_role"]:
                raise ValueError("d28_truth_role_profile_mismatch")
            combined.append({**by_key[(scenario_hash, profile)], "truth_role": truth["truth_role"], "age_generations": truth["age_generations"]})
    if len(truth_roles) != 10 or list(truth_roles.values()).count("POSITIVE") != 5:
        raise ValueError("d28_truth_balance_invalid")

    label_by_hash: dict[str, str] = {}
    for role, prefix in (("POSITIVE", "P"), ("NEGATIVE", "N")):
        hashes = sorted(scenario_hash for scenario_hash, value in truth_roles.items() if value == role)
        label_by_hash.update({scenario_hash: f"{prefix}{index}" for index, scenario_hash in enumerate(hashes, 1)})
    case_results = [{
        "scenario_label": label_by_hash[row["scenario_id_hash"]],
        "truth_role": row["truth_role"], "age_generations": row["age_generations"],
        "profile": row["profile"], "evaluation_status": row["evaluation_status"],
        "non_evaluable_reason": row["non_evaluable_reason"],
        "observed_total_shared_cm": row["observed_total_shared_cm"],
        "attempted_null_draws": row["attempted_null_draws"],
        "evaluable_null_draws": row["evaluable_null_draws"],
        "non_evaluable_null_draws": row["non_evaluable_null_draws"],
        "exceedance_count": row["exceedance_count"],
        "empirical_probability": row["empirical_probability"],
        "monte_carlo_interval_low": row["monte_carlo_interval_low"],
        "monte_carlo_interval_high": row["monte_carlo_interval_high"],
        "benchmark_positive_call": row["benchmark_positive_call"],
        "production_status": row["production_status"],
    } for row in combined]

    profile_results = [_profile_result(combined, profile) for profile in config["profiles"]]
    primary = next(row for row in profile_results if row["profile"] == config["primary_profile"])
    emitted = sum(
        row["production_status"] in {"UNUSUAL_TARGET_CENTERED_SHARING", "NO_UNUSUAL_SHARING_DETECTED"}
        for row in combined
    )
    reasons: dict[tuple[str, str, str], int] = {}
    for row in combined:
        if row["evaluation_status"] != "EVALUATED":
            key = (row["profile"], row["truth_role"], row["non_evaluable_reason"] or "UNSPECIFIED")
            reasons[key] = reasons.get(key, 0) + 1
    non_evaluations = [
        {"profile": profile, "truth_role": role, "reason": reason, "count": count}
        for (profile, role, reason), count in sorted(reasons.items())
    ]
    criteria = [
        {"criterion_id": "CAL-ENRICH-001", "metric": "classification_emitted_without_prespecified_threshold_count", "value": emitted, "denominator": 40, "threshold": 0, "lifecycle": "fixed_contract", "status": "PASS" if emitted == 0 else "FAIL"},
        {"criterion_id": "CAL-ENRICH-002", "metric": "null_false_positive_rate", "value": primary["null_false_positive_rate"], "denominator": primary["null_total"], "threshold": None, "lifecycle": "pilot_then_seal", "status": "REVIEW_REQUIRED"},
        {"criterion_id": "CAL-ENRICH-003", "metric": "positive_scenario_power", "value": primary["power"], "denominator": primary["positive_total"], "threshold": None, "lifecycle": "pilot_then_seal", "status": "REVIEW_REQUIRED"},
    ]
    technical_status = "FAIL" if emitted else "PASS_CORE_EVALUATOR"
    document = {
        "schema_version": "1.0.0", "evaluation_state": "CLOSED",
        "report_id": "CAL-REP-16B", "method_id": "d28_enrichment_evaluation_v1",
        "configuration_sha256": sha256_file(configuration_path),
        "protocol_sha256": sha256_file(PROTOCOL), "calls_manifest_sha256": sha256_file(manifest_path),
        "analysis_scope": "D24_SYNTHETIC_CHR5_THREE_FAMILY_SCENARIOS",
        "stage_execution_scope": "CORE_16B_PRIMITIVES_NOT_END_TO_END",
        "primary_profile": config["primary_profile"],
        "scenario_counts": {"independent_total": 10, "positive": 5, "null": 5, "profiled_case_total": 40},
        "benchmark_decision_alpha": config["benchmark_decision_alpha"],
        "production_classification_threshold": None,
        "primary_metrics": {
            "null_false_positive_rate": primary["null_false_positive_rate"],
            "positive_scenario_power": primary["power"],
            "classification_emitted_without_prespecified_threshold_count": emitted,
        },
        "profile_results": profile_results, "case_results": case_results,
        "non_evaluations": non_evaluations,
        "criteria": criteria, "technical_status": technical_status,
        "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
        "production_promotion": False, "calibrated": False, "founder_effect_proven": False,
        "checks": [
            "CALLS_SEALED_BEFORE_TRUTH_OPENED", "PLUS_ONE_EMPIRICAL_PROBABILITY",
            "TARGET_EXCLUDED_FROM_SIGNATURE_BY_SHARED_PRIMITIVE", "NON_EVALUATIONS_PRESERVED",
            "PRIMARY_COUNTS_TEN_INDEPENDENT_SCENARIOS", "PRODUCTION_THRESHOLD_REMAINS_NULL",
        ],
        "limits": [
            "Cinq scénarios par rôle donnent des intervalles très larges.",
            "Tous les scénarios comportent exactement trois familles indépendantes.",
            "Le fond synthétique D24 est partagé entre scénarios et ne remplace pas 1000 Genomes.",
            "Le phasage-vérité D24 n'évalue pas la propagation des erreurs SHAPEIT5.",
            "Le banc appelle les primitives de 16B, pas son orchestration ni bcftools.",
            "Le seuil alpha du banc n'est ni un seuil de production ni un seuil de qualification.",
        ],
    }
    validate_json_document(document, "enrichment_evaluation.schema.json", schemas_dir=SCHEMAS)
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        result_path = temporary / "CAL-REP-16B.json"
        write_json(result_path, document)
        write_json(temporary / "evaluation_manifest.json", {
            "schema_version": "1.0.0", "evaluation_state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL), "calls_manifest_sha256": sha256_file(manifest_path),
            "artifacts": [{"path": result_path.name, "sha256": sha256_file(result_path), "size_bytes": result_path.stat().st_size}],
        })
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--configuration-path", type=Path, required=True)
    common.add_argument("--d24-bench-root", type=Path, required=True)
    run_parser = subparsers.add_parser("run", parents=[common])
    run_parser.add_argument("--workspace-root", type=Path, required=True)
    run_parser.add_argument("--run-id", required=True)
    evaluate_parser = subparsers.add_parser("evaluate", parents=[common])
    evaluate_parser.add_argument("--calls-root", type=Path, required=True)
    evaluate_parser.add_argument("--output-root", type=Path, required=True)
    arguments = vars(parser.parse_args())
    command = arguments.pop("command")
    path = run_calls(**arguments) if command == "run" else evaluate(**arguments)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
