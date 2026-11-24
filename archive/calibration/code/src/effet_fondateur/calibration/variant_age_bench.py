"""Banc D25 fermé pour la calibration de la datation Gamma."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
import tempfile
from typing import Any, Iterable

import numpy as np

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.founder_haplotype_oracle import wilson_interval
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document, validate_tsv_table
from effet_fondateur.dating import (
    GammaAgeError,
    chance_sharing_correction_cm,
    estimate_correlated_gamma,
    estimate_independent_gamma,
    publish_variant_age,
)


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
PROTOCOL = REPOSITORY / "calibration/protocoles/D25_VARIANT_AGE_PRESPECIFICATION.md"


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
        raise ValueError("d25_source_outside_workspace")
    return resolved


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge le contrat D25 et vérifie les invariants non exprimés en JSON Schema."""
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(document, "variant_age_bench_config.schema.json", schemas_dir=SCHEMAS)
    profiles = [row["profile"] for row in document["profiles"]]
    scenarios = [row["scenario_id"] for row in document["d24_positive_scenarios"]]
    if len(profiles) != len(set(profiles)) or len(scenarios) != len(set(scenarios)):
        raise ValueError("d25_configuration_duplicates")
    if document["simulation_seed"] == document["bootstrap_seed"]:
        raise ValueError("d25_seeds_not_independent")
    return document


def _write_jsonl_gz(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            for row in rows:
                compressed.write(
                    (json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode()
                )


def _read_jsonl_gz(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _case_seed(master: int, age: int, families: int, replication: int) -> int:
    digest = hashlib.sha256(f"{master}:{age}:{families}:{replication}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def materialize_bench(
    *, configuration_path: Path, workspace_root: Path, run_id: str
) -> Path:
    """Matérialise les longueurs aveugles et conserve l’âge vrai à part."""
    config = load_configuration(configuration_path)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    blinded: list[dict[str, Any]] = []
    truth: list[dict[str, Any]] = []
    try:
        for age in config["ages_generations"]:
            rate_per_cm = age / 100.0
            for family_count in config["family_counts"]:
                for replication in range(1, config["replications_per_cell"] + 1):
                    seed = _case_seed(config["simulation_seed"], age, family_count, replication)
                    generator = random.Random(seed)
                    case_id = "D25_" + hashlib.sha256(str(seed).encode()).hexdigest()[:20]
                    left = [generator.expovariate(rate_per_cm) for _ in range(family_count)]
                    right = [generator.expovariate(rate_per_cm) for _ in range(family_count)]
                    blinded.append({
                        "case_id": case_id,
                        "family_count": family_count,
                        "left_lengths_cm": left,
                        "right_lengths_cm": right,
                    })
                    truth.append({
                        "case_id": case_id,
                        "age_generations": age,
                        "family_count": family_count,
                        "replication": replication,
                        "case_seed_sha256": hashlib.sha256(str(seed).encode()).hexdigest(),
                    })
        blinded_path = temporary / "blinded_inputs/age_lengths.jsonl.gz"
        truth_path = temporary / "private_truth/age_truth.jsonl.gz"
        _write_jsonl_gz(blinded_path, blinded)
        _write_jsonl_gz(truth_path, truth)
        manifest = {
            "schema_version": "1.0.0",
            "method_id": "d25_exponential_age_grid_v1",
            "state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "case_count": len(blinded),
            "age_count": len(config["ages_generations"]),
            "family_count_strata": config["family_counts"],
            "replications_per_cell": config["replications_per_cell"],
            "blinded_inputs_sha256": sha256_file(blinded_path),
            "private_truth_sha256": sha256_file(truth_path),
            "truth_separated": True,
            "calibrated": False,
            "founder_effect_proven": False,
        }
        write_json(temporary / "bench_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _profile_lengths(
    left: list[float], right: list[float], profile: dict[str, Any]
) -> tuple[list[float], list[float]]:
    kind, value = profile["kind"], float(profile["value"])
    if kind == "EXACT":
        return list(left), list(right)
    if kind == "ARM_CAP_CM":
        return [min(value, item) for item in left], [min(value, item) for item in right]
    if kind == "MAP_SCALE":
        return [value * item for item in left], [value * item for item in right]
    raise ValueError("d25_profile_kind_invalid")


def _estimate(
    model: str, left: list[float], right: list[float], confidence_level: float
) -> dict[str, Any]:
    try:
        if model == "CORRELATED":
            estimate = estimate_correlated_gamma(left, right, confidence_level=confidence_level)
        elif model == "INDEPENDENT":
            estimate = estimate_independent_gamma(left, right, confidence_level=confidence_level)
        else:
            raise ValueError("d25_model_invalid")
    except GammaAgeError as error:
        return {
            "status": "NOT_ESTIMATED",
            "exclusion_code": str(error),
            "estimate_generations": None,
            "ci_lower_generations": None,
            "ci_upper_generations": None,
            "effective_n": None,
            "rho": None,
        }
    return {
        "status": "ESTIMATED",
        "exclusion_code": None,
        "estimate_generations": estimate.generations,
        "ci_lower_generations": estimate.confidence_lower,
        "ci_upper_generations": estimate.confidence_upper,
        "effective_n": estimate.effective_sample_size,
        "rho": estimate.correlation,
    }


def _lofo(
    left: list[float], right: list[float], full_estimate: float, confidence_level: float
) -> list[float]:
    if len(left) - 1 < 3:
        return []
    values = []
    for omitted in range(len(left)):
        retained_left = left[:omitted] + left[omitted + 1 :]
        retained_right = right[:omitted] + right[omitted + 1 :]
        estimate = estimate_correlated_gamma(
            retained_left, retained_right, confidence_level=confidence_level
        )
        values.append(estimate.generations - full_estimate)
    return values


def _d24_rows(
    *, d24_root: Path, temporary: Path, config: dict[str, Any]
) -> list[dict[str, Any]]:
    d24_root = _workspace_path(d24_root, directory=True)
    rows: list[dict[str, Any]] = []
    for scenario in config["d24_positive_scenarios"]:
        scenario_id = scenario["scenario_id"]
        for profile in config["d24_profiles"]:
            stage13 = d24_root / "cases" / scenario_id / profile / "stage13"
            segments = _workspace_path(stage13 / "founder_segments.tsv")
            consensus = _workspace_path(stage13 / "founder_consensus.tsv")
            summary = _workspace_path(stage13 / "founder_analysis_summary.json")
            output = temporary / "d24_stage14" / scenario_id / profile
            publication = publish_variant_age(
                founder_segments_path=segments,
                founder_consensus_path=consensus,
                founder_summary_path=summary,
                output_dir=output,
                confidence_level=config["confidence_level"],
                minimum_units_for_estimate=config["minimum_units_for_estimate"],
                minimum_units_for_primary=config["minimum_units_for_primary"],
                generation_years=config["generation_years"],
                leave_one_family_out=True,
                chance_sharing_parameters=None,
            )
            estimates = validate_tsv_table(
                publication.estimates_path, "variant_age_estimates.schema.json"
            ).rows
            correlated = next(row for row in estimates if row["MODEL"] == "CORRELATED")
            rows.append({
                "source": "D24_PROPAGATED",
                "case_id": hashlib.sha256(scenario_id.encode()).hexdigest()[:20],
                "scenario_key": scenario_id,
                "profile": profile,
                "model": "CORRELATED",
                "family_count": publication.included_unit_count,
                "status": "ESTIMATED" if correlated["ESTIMATE_GENERATIONS"] is not None else "NOT_ESTIMATED",
                "exclusion_code": correlated["EXCLUSION_CODE"],
                "estimate_generations": correlated["ESTIMATE_GENERATIONS"],
                "ci_lower_generations": correlated["CI_LOWER_GENERATIONS"],
                "ci_upper_generations": correlated["CI_UPPER_GENERATIONS"],
                "effective_n": correlated["EFFECTIVE_N"],
                "rho": correlated["RHO"],
                "lofo_deltas": [],
                "source_sha256": hashlib.sha256(
                    "".join(sha256_file(path) for path in (segments, consensus, summary)).encode()
                ).hexdigest(),
            })
    return rows


def run_bench(
    *, bench_root: Path, configuration_path: Path, d24_root: Path,
    workspace_root: Path, run_id: str,
) -> Path:
    """Produit les estimations avant toute ouverture des âges vrais."""
    config = load_configuration(configuration_path)
    bench_root = _workspace_path(bench_root, directory=True)
    manifest_path = _workspace_path(bench_root / "bench_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    blinded_path = _workspace_path(bench_root / "blinded_inputs/age_lengths.jsonl.gz")
    if (
        manifest["configuration_sha256"] != sha256_file(configuration_path)
        or manifest["blinded_inputs_sha256"] != sha256_file(blinded_path)
    ):
        raise ValueError("d25_bench_binding_invalid")
    blinded = _read_jsonl_gz(blinded_path)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    rows: list[dict[str, Any]] = []
    try:
        for case in blinded:
            for profile in config["profiles"]:
                left, right = _profile_lengths(
                    case["left_lengths_cm"], case["right_lengths_cm"], profile
                )
                for model in (config["primary_model"], config["reference_model"]):
                    result = _estimate(model, left, right, config["confidence_level"])
                    deltas: list[float] = []
                    if (
                        profile["profile"] == "ORACLE_EXACT"
                        and model == "CORRELATED"
                        and result["estimate_generations"] is not None
                    ):
                        try:
                            deltas = _lofo(
                                left, right, float(result["estimate_generations"]),
                                config["confidence_level"],
                            )
                        except GammaAgeError:
                            deltas = []
                    rows.append({
                        "source": "SIMULATED_GRID",
                        "case_id": case["case_id"],
                        "scenario_key": None,
                        "profile": profile["profile"],
                        "model": model,
                        "family_count": case["family_count"],
                        **result,
                        "lofo_deltas": deltas,
                        "source_sha256": manifest["blinded_inputs_sha256"],
                    })
        rows.extend(_d24_rows(d24_root=d24_root, temporary=temporary, config=config))
        output_path = temporary / "estimates.jsonl.gz"
        _write_jsonl_gz(output_path, rows)
        run_manifest = {
            "schema_version": "1.0.0",
            "method_id": "gamma_gandolfo_2014_v1",
            "state": "CLOSED_BEFORE_TRUTH",
            "configuration_sha256": sha256_file(configuration_path),
            "bench_manifest_sha256": sha256_file(manifest_path),
            "blinded_inputs_sha256": sha256_file(blinded_path),
            "estimate_row_count": len(rows),
            "estimates_sha256": sha256_file(output_path),
            "truth_opened": False,
            "calibrated": False,
            "founder_effect_proven": False,
        }
        write_json(temporary / "run_manifest.json", run_manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _bootstrap_summary(errors: list[float], draws: int, seed: int) -> dict[str, Any]:
    values = np.asarray(errors, dtype=float)
    if values.size == 0:
        return {name: None for name in (
            "bias", "bias_ci95_low", "bias_ci95_high", "mae", "mae_ci95_low",
            "mae_ci95_high", "rmse", "rmse_ci95_low", "rmse_ci95_high",
        )}
    generator = np.random.default_rng(seed)
    bias_draws = np.empty(draws)
    mae_draws = np.empty(draws)
    rmse_draws = np.empty(draws)
    for index in range(draws):
        sampled = values[generator.integers(0, values.size, values.size)]
        bias_draws[index] = sampled.mean()
        mae_draws[index] = np.abs(sampled).mean()
        rmse_draws[index] = np.sqrt(np.square(sampled).mean())
    return {
        "bias": float(values.mean()),
        "bias_ci95_low": float(np.quantile(bias_draws, 0.025)),
        "bias_ci95_high": float(np.quantile(bias_draws, 0.975)),
        "mae": float(np.abs(values).mean()),
        "mae_ci95_low": float(np.quantile(mae_draws, 0.025)),
        "mae_ci95_high": float(np.quantile(mae_draws, 0.975)),
        "rmse": float(np.sqrt(np.square(values).mean())),
        "rmse_ci95_low": float(np.quantile(rmse_draws, 0.025)),
        "rmse_ci95_high": float(np.quantile(rmse_draws, 0.975)),
    }


def _summary(
    rows: list[dict[str, Any]], truth: dict[str, float], *, draws: int, seed: int
) -> dict[str, Any]:
    estimated = [row for row in rows if row["estimate_generations"] is not None]
    errors = [float(row["estimate_generations"]) - truth[row["case_id"]] for row in estimated]
    covered = sum(
        float(row["ci_lower_generations"]) <= truth[row["case_id"]] <= float(row["ci_upper_generations"])
        for row in estimated
    )
    low, high = wilson_interval(covered, len(estimated)) if estimated else (None, None)
    return {
        "total_count": len(rows),
        "estimated_count": len(estimated),
        "not_evaluable_count": len(rows) - len(estimated),
        **_bootstrap_summary(errors, draws, seed),
        "coverage_numerator": covered,
        "coverage_denominator": len(estimated),
        "coverage": None if not estimated else covered / len(estimated),
        "coverage_ci95_low": low,
        "coverage_ci95_high": high,
    }


def _seed(base: int, label: str) -> int:
    return base + int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], "big")


def evaluate_bench(
    *, bench_root: Path, run_root: Path, configuration_path: Path,
    workspace_root: Path, run_id: str,
) -> Path:
    """Ouvre la vérité après les estimations et publie uniquement des agrégats."""
    config = load_configuration(configuration_path)
    bench_root = _workspace_path(bench_root, directory=True)
    run_root = _workspace_path(run_root, directory=True)
    bench_manifest_path = _workspace_path(bench_root / "bench_manifest.json")
    run_manifest_path = _workspace_path(run_root / "run_manifest.json")
    estimates_path = _workspace_path(run_root / "estimates.jsonl.gz")
    bench_manifest = json.loads(bench_manifest_path.read_text(encoding="utf-8"))
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    if (
        run_manifest["configuration_sha256"] != sha256_file(configuration_path)
        or run_manifest["bench_manifest_sha256"] != sha256_file(bench_manifest_path)
        or run_manifest["estimates_sha256"] != sha256_file(estimates_path)
    ):
        raise ValueError("d25_evaluation_binding_invalid")
    estimates = _read_jsonl_gz(estimates_path)  # ouvert avant la vérité
    truth_path = _workspace_path(bench_root / "private_truth/age_truth.jsonl.gz")
    if bench_manifest["private_truth_sha256"] != sha256_file(truth_path):
        raise ValueError("d25_truth_checksum_invalid")
    truth_rows = _read_jsonl_gz(truth_path)
    truth = {row["case_id"]: float(row["age_generations"]) for row in truth_rows}
    family_truth = {row["case_id"]: int(row["family_count"]) for row in truth_rows}
    d24_truth = {
        hashlib.sha256(row["scenario_id"].encode()).hexdigest()[:20]: float(row["age_generations"])
        for row in config["d24_positive_scenarios"]
    }
    grid_rows = [row for row in estimates if row["source"] == "SIMULATED_GRID"]
    d24_rows = [row for row in estimates if row["source"] == "D24_PROPAGATED"]
    primary_rows = [
        row for row in grid_rows
        if row["profile"] == "ORACLE_EXACT" and row["model"] == "CORRELATED"
    ]
    primary = _summary(
        primary_rows, truth, draws=config["bootstrap_draws"],
        seed=_seed(config["bootstrap_seed"], "primary"),
    )
    by_age = []
    for age in config["ages_generations"]:
        selected = [row for row in primary_rows if truth[row["case_id"]] == age]
        by_age.append({
            "age_generations": age,
            **_summary(selected, truth, draws=config["bootstrap_draws"], seed=_seed(config["bootstrap_seed"], f"age:{age}")),
        })
    by_family = []
    for family_count in config["family_counts"]:
        selected = [row for row in primary_rows if family_truth[row["case_id"]] == family_count]
        by_family.append({
            "family_count": family_count,
            **_summary(selected, truth, draws=config["bootstrap_draws"], seed=_seed(config["bootstrap_seed"], f"family:{family_count}")),
        })
    profiles = []
    for profile in [row["profile"] for row in config["profiles"]]:
        for model in (config["primary_model"], config["reference_model"]):
            selected = [row for row in grid_rows if row["profile"] == profile and row["model"] == model]
            profiles.append({
                "profile": profile,
                "model": model,
                **_summary(selected, truth, draws=config["bootstrap_draws"], seed=_seed(config["bootstrap_seed"], f"profile:{profile}:{model}")),
            })
    d24_profiles = []
    for profile in config["d24_profiles"]:
        selected = [row for row in d24_rows if row["profile"] == profile]
        d24_profiles.append({
            "profile": profile,
            **_summary(selected, d24_truth, draws=config["bootstrap_draws"], seed=_seed(config["bootstrap_seed"], f"d24:{profile}")),
        })
    instabilities = [
        max(abs(float(value)) for value in row["lofo_deltas"])
        for row in primary_rows if row["lofo_deltas"]
    ]
    official_left = [1.2, 3.4, 1.6, 4.8, 4.8, 2.3]
    official_right = [2.0, 3.8, 2.7, 2.1, 3.8, 1.8]
    correction = chance_sharing_correction_cm(0.7, 68_378, 290)
    historical_correlated = estimate_correlated_gamma(
        official_left, official_right, chance_sharing_correction_cm=correction
    )
    historical_independent = estimate_independent_gamma(
        official_left, official_right, chance_sharing_correction_cm=correction
    )
    historical_error = max(
        abs(historical_correlated.generations - 22.3),
        abs(historical_independent.generations - 27.5),
    )
    criteria = [
        {"criterion_id": "CAL-AGE-001", "metric": "signed_age_bias", "value": primary["bias"], "threshold": None, "status": "REVIEW_REQUIRED"},
        {"criterion_id": "CAL-AGE-002", "metric": "age_mean_absolute_error", "value": primary["mae"], "threshold": None, "status": "REVIEW_REQUIRED"},
        {"criterion_id": "CAL-AGE-003", "metric": "age_rmse", "value": primary["rmse"], "threshold": None, "status": "REVIEW_REQUIRED"},
        {"criterion_id": "CAL-AGE-004", "metric": "nominal_95_interval_coverage", "value": primary["coverage"], "threshold": None, "status": "REVIEW_REQUIRED"},
        {"criterion_id": "CAL-AGE-005", "metric": "leave_one_family_out_instability", "value": float(np.mean(instabilities)), "threshold": None, "status": "REVIEW_REQUIRED"},
        {"criterion_id": "CAL-AGE-006", "metric": "historical_gamma_reference_example_error", "value": historical_error, "threshold": 0.1, "status": "PASS" if historical_error <= 0.1 else "FAIL"},
    ]
    document = {
        "schema_version": "1.0.0",
        "report_id": "CAL-REP-14",
        "method_id": "gamma_age_truth_bench_v1",
        "evaluation_state": "CLOSED",
        "configuration_sha256": sha256_file(configuration_path),
        "protocol_sha256": sha256_file(PROTOCOL),
        "bench_manifest_sha256": sha256_file(bench_manifest_path),
        "run_manifest_sha256": sha256_file(run_manifest_path),
        "technical_status": "PASS" if historical_error <= 0.1 else "FAIL",
        "performance_qualification": "REVIEW_REQUIRED",
        "design": {
            "ages_generations": config["ages_generations"],
            "family_counts": config["family_counts"],
            "replications_per_cell": config["replications_per_cell"],
            "independent_scenario_count": len(truth),
            "primary_model": config["primary_model"],
            "primary_profile": "ORACLE_EXACT",
            "confidence_level": config["confidence_level"],
        },
        "primary_metrics": primary,
        "results_by_age": by_age,
        "results_by_family_count": by_family,
        "profile_results": profiles,
        "d24_propagation_results": d24_profiles,
        "leave_one_family_out": {
            "evaluable_scenario_count": len(instabilities),
            "mean_max_absolute_delta_generations": float(np.mean(instabilities)),
            "p95_max_absolute_delta_generations": float(np.quantile(instabilities, 0.95)),
            "three_family_status": "NOT_EVALUABLE_AFTER_ONE_OMISSION",
        },
        "historical_reference": {
            "correlated_estimate_generations": historical_correlated.generations,
            "independent_estimate_generations": historical_independent.generations,
            "maximum_absolute_reference_error_generations": historical_error,
        },
        "phase_error_evaluation": {
            "status": "NOT_EVALUABLE",
            "reason": config["phase_error_policy"],
        },
        "criteria": criteria,
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(document, "variant_age_evaluation.schema.json", schemas_dir=SCHEMAS)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        output = temporary / "CAL-REP-14.json"
        write_json(output, document)
        (temporary / "checksums.sha256").write_text(
            f"{sha256_file(output)}  CAL-REP-14.json\n", encoding="utf-8"
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def main() -> int:
    """Expose les trois phases D25 sans mélanger vérité et estimation."""
    parser = argparse.ArgumentParser(description="Banc fermé D25 de datation Gamma")
    subparsers = parser.add_subparsers(dest="command", required=True)
    materialize = subparsers.add_parser("materialize")
    materialize.add_argument("--configuration-path", type=Path, required=True)
    materialize.add_argument("--workspace-root", type=Path, required=True)
    materialize.add_argument("--run-id", required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--bench-root", type=Path, required=True)
    run.add_argument("--configuration-path", type=Path, required=True)
    run.add_argument("--d24-root", type=Path, required=True)
    run.add_argument("--workspace-root", type=Path, required=True)
    run.add_argument("--run-id", required=True)
    evaluate = subparsers.add_parser("evaluate")
    evaluate.add_argument("--bench-root", type=Path, required=True)
    evaluate.add_argument("--run-root", type=Path, required=True)
    evaluate.add_argument("--configuration-path", type=Path, required=True)
    evaluate.add_argument("--workspace-root", type=Path, required=True)
    evaluate.add_argument("--run-id", required=True)
    arguments = vars(parser.parse_args())
    command = arguments.pop("command")
    result = {
        "materialize": materialize_bench,
        "run": run_bench,
        "evaluate": evaluate_bench,
    }[command](**arguments)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
