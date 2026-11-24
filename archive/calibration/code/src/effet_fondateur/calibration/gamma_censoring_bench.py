"""Banc D25b de vraisemblance Gamma avec censure de type II."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
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
from effet_fondateur.calibration.variant_age_bench import _bootstrap_summary, _seed
from effet_fondateur.contracts import validate_json_document
from effet_fondateur.dating import (
    GammaAgeError,
    estimate_correlated_gamma_r_faithful,
    estimate_independent_gamma_r_faithful,
    estimate_independent_gamma_type2_censored,
)


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
PROTOCOL = (
    REPOSITORY
    / "calibration/protocoles/D25B_GAMMA_CENSORING_PRESPECIFICATION.md"
)
METHODS = (
    "TYPE2_CENSORED_INDEPENDENT_V3",
    "R_FAITHFUL_INDEPENDENT_V2",
    "R_FAITHFUL_CORRELATED_V2",
    "FULL_INDEPENDENT_ORACLE_V1",
)


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
        or any(
            left == "data" and right == "input"
            for left, right in zip(folded, folded[1:])
        )
    ):
        raise ValueError("d25b_source_outside_workspace")
    return resolved


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge la configuration D25b et contrôle ses invariants croisés."""
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(
        document, "gamma_censoring_config.schema.json", schemas_dir=SCHEMAS
    )
    if document["simulation_seed"] == document["bootstrap_seed"]:
        raise ValueError("d25b_seeds_not_independent")
    fixture = document["analytic_fixture"]
    if len(fixture["left_lengths_cm"]) != len(fixture["right_lengths_cm"]):
        raise ValueError("d25b_fixture_arm_count_mismatch")
    return document


def _write_jsonl_gz(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0
        ) as compressed:
            for row in rows:
                payload = json.dumps(row, sort_keys=True, separators=(",", ":"))
                compressed.write((payload + "\n").encode())


def _read_jsonl_gz(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _case_seed(
    master_seed: int,
    genealogy_profile: str,
    age_generations: int,
    family_count: int,
    replication: int,
) -> int:
    key = (
        f"{master_seed}:{genealogy_profile}:{age_generations}:"
        f"{family_count}:{replication}"
    )
    return int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")


def _simulate_arm_lengths_cm(
    generator: random.Random,
    *,
    age_generations: int,
    family_count: int,
    shared_path_fraction: float,
) -> list[float]:
    total_rate_per_cm = age_generations / 100.0
    if shared_path_fraction == 0:
        return [
            generator.expovariate(total_rate_per_cm)
            for _ in range(family_count)
        ]
    shared_breakpoint = generator.expovariate(
        total_rate_per_cm * shared_path_fraction
    )
    private_rate = total_rate_per_cm * (1.0 - shared_path_fraction)
    return [
        min(shared_breakpoint, generator.expovariate(private_rate))
        for _ in range(family_count)
    ]


def pairwise_type2_observation(lengths_cm: list[float]) -> list[float]:
    """Remplace l'unique maximum vrai par le seuil du deuxième maximum."""
    if len(lengths_cm) < 2:
        raise ValueError("d25b_pairwise_observation_requires_two_arms")
    observed = list(lengths_cm)
    longest_index = max(range(len(observed)), key=observed.__getitem__)
    observed[longest_index] = sorted(observed)[-2]
    return observed


def materialize_bench(
    *, configuration_path: Path, workspace_root: Path, run_id: str
) -> Path:
    """Matérialise le panel D25b avec âges séparés des longueurs aveugles."""
    config = load_configuration(configuration_path)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    blinded: list[dict[str, Any]] = []
    truth: list[dict[str, Any]] = []
    try:
        for genealogy in config["genealogy_profiles"]:
            genealogy_profile = genealogy["profile"]
            shared_fraction = float(genealogy["shared_path_fraction"])
            for age in config["ages_generations"]:
                for family_count in config["family_counts"]:
                    for replication in range(
                        1, config["replications_per_cell"] + 1
                    ):
                        seed = _case_seed(
                            config["simulation_seed"], genealogy_profile,
                            age, family_count, replication,
                        )
                        generator = random.Random(seed)
                        left_full = _simulate_arm_lengths_cm(
                            generator,
                            age_generations=age,
                            family_count=family_count,
                            shared_path_fraction=shared_fraction,
                        )
                        right_full = _simulate_arm_lengths_cm(
                            generator,
                            age_generations=age,
                            family_count=family_count,
                            shared_path_fraction=shared_fraction,
                        )
                        case_id = (
                            "D25B_"
                            + hashlib.sha256(str(seed).encode()).hexdigest()[:20]
                        )
                        blinded.append({
                            "case_id": case_id,
                            "genealogy_profile": genealogy_profile,
                            "family_count": family_count,
                            "full_left_lengths_cm": left_full,
                            "full_right_lengths_cm": right_full,
                            "observed_left_lengths_cm": (
                                pairwise_type2_observation(left_full)
                            ),
                            "observed_right_lengths_cm": (
                                pairwise_type2_observation(right_full)
                            ),
                        })
                        truth.append({
                            "case_id": case_id,
                            "age_generations": age,
                            "genealogy_profile": genealogy_profile,
                            "shared_path_fraction": shared_fraction,
                            "family_count": family_count,
                            "replication": replication,
                            "case_seed_sha256": hashlib.sha256(
                                str(seed).encode()
                            ).hexdigest(),
                        })
        blinded_path = temporary / "blinded_inputs/gamma_lengths.jsonl.gz"
        truth_path = temporary / "private_truth/gamma_age_truth.jsonl.gz"
        _write_jsonl_gz(blinded_path, blinded)
        _write_jsonl_gz(truth_path, truth)
        manifest = {
            "schema_version": "1.0.0",
            "method_id": "d25b_gamma_type2_grid_v1",
            "state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "case_count": len(blinded),
            "blinded_inputs_sha256": sha256_file(blinded_path),
            "private_truth_sha256": sha256_file(truth_path),
            "truth_separated": True,
            "prior_panel_seed_reused": False,
            "calibrated": False,
            "founder_effect_proven": False,
        }
        write_json(temporary / "bench_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _estimate(
    method: str,
    *,
    full_left: list[float],
    full_right: list[float],
    observed_left: list[float],
    observed_right: list[float],
    confidence_level: float,
) -> dict[str, Any]:
    try:
        if method == "TYPE2_CENSORED_INDEPENDENT_V3":
            estimate = estimate_independent_gamma_type2_censored(
                observed_left, observed_right,
                confidence_level=confidence_level,
            )
        elif method == "R_FAITHFUL_INDEPENDENT_V2":
            estimate = estimate_independent_gamma_r_faithful(
                observed_left, observed_right,
                confidence_level=confidence_level,
                apply_boundary_correction=True,
            )
        elif method == "R_FAITHFUL_CORRELATED_V2":
            estimate = estimate_correlated_gamma_r_faithful(
                observed_left, observed_right,
                confidence_level=confidence_level,
                apply_boundary_correction=True,
            )
        elif method == "FULL_INDEPENDENT_ORACLE_V1":
            estimate = estimate_independent_gamma_r_faithful(
                full_left, full_right,
                confidence_level=confidence_level,
                apply_boundary_correction=False,
            )
        else:
            raise ValueError("d25b_method_invalid")
    except GammaAgeError as error:
        return {
            "estimate_generations": None,
            "ci_lower_generations": None,
            "ci_upper_generations": None,
            "exclusion_code": str(error),
        }
    return {
        "estimate_generations": estimate.generations,
        "ci_lower_generations": estimate.confidence_lower,
        "ci_upper_generations": estimate.confidence_upper,
        "exclusion_code": None,
    }


def _lofo_max_delta(
    *,
    full_left: list[float],
    full_right: list[float],
    full_estimate: float,
    confidence_level: float,
) -> float | None:
    if len(full_left) - 1 < 3:
        return None
    deltas: list[float] = []
    for omitted in range(len(full_left)):
        left = full_left[:omitted] + full_left[omitted + 1:]
        right = full_right[:omitted] + full_right[omitted + 1:]
        result = _estimate(
            "TYPE2_CENSORED_INDEPENDENT_V3",
            full_left=left,
            full_right=right,
            observed_left=pairwise_type2_observation(left),
            observed_right=pairwise_type2_observation(right),
            confidence_level=confidence_level,
        )
        if result["estimate_generations"] is None:
            return None
        deltas.append(
            abs(float(result["estimate_generations"]) - full_estimate)
        )
    return max(deltas)


def run_bench(
    *, bench_root: Path, configuration_path: Path,
    workspace_root: Path, run_id: str,
) -> Path:
    """Ferme les estimations D25b avant toute ouverture de la vérité d'âge."""
    config = load_configuration(configuration_path)
    bench_root = _workspace_path(bench_root, directory=True)
    manifest_path = _workspace_path(bench_root / "bench_manifest.json")
    blinded_path = _workspace_path(
        bench_root / "blinded_inputs/gamma_lengths.jsonl.gz"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest["configuration_sha256"] != sha256_file(configuration_path)
        or manifest["protocol_sha256"] != sha256_file(PROTOCOL)
        or manifest["blinded_inputs_sha256"] != sha256_file(blinded_path)
    ):
        raise ValueError("d25b_bench_binding_invalid")
    blinded = _read_jsonl_gz(blinded_path)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    rows: list[dict[str, Any]] = []
    try:
        for case in blinded:
            candidate_result: dict[str, Any] | None = None
            for method in METHODS:
                result = _estimate(
                    method,
                    full_left=case["full_left_lengths_cm"],
                    full_right=case["full_right_lengths_cm"],
                    observed_left=case["observed_left_lengths_cm"],
                    observed_right=case["observed_right_lengths_cm"],
                    confidence_level=config["confidence_level"],
                )
                if method == "TYPE2_CENSORED_INDEPENDENT_V3":
                    candidate_result = result
                rows.append({
                    "case_id": case["case_id"],
                    "genealogy_profile": case["genealogy_profile"],
                    "family_count": case["family_count"],
                    "method": method,
                    "lofo_max_absolute_delta_generations": None,
                    **result,
                })
            if (
                candidate_result is not None
                and candidate_result["estimate_generations"] is not None
            ):
                rows[-4]["lofo_max_absolute_delta_generations"] = (
                    _lofo_max_delta(
                        full_left=case["full_left_lengths_cm"],
                        full_right=case["full_right_lengths_cm"],
                        full_estimate=float(
                            candidate_result["estimate_generations"]
                        ),
                        confidence_level=config["confidence_level"],
                    )
                )
        estimates_path = temporary / "estimates.jsonl.gz"
        _write_jsonl_gz(estimates_path, rows)
        run_manifest = {
            "schema_version": "1.0.0",
            "method_id": "d25b_gamma_type2_method_matrix_v1",
            "state": "CLOSED_BEFORE_TRUTH",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "bench_manifest_sha256": sha256_file(manifest_path),
            "blinded_inputs_sha256": sha256_file(blinded_path),
            "estimate_row_count": len(rows),
            "estimates_sha256": sha256_file(estimates_path),
            "truth_opened": False,
            "production_promotion": False,
            "calibrated": False,
            "founder_effect_proven": False,
        }
        write_json(temporary / "run_manifest.json", run_manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _summary(
    rows: list[dict[str, Any]], truth: dict[str, float],
    *, draws: int, seed: int,
) -> dict[str, Any]:
    estimated = [row for row in rows if row["estimate_generations"] is not None]
    errors = [
        float(row["estimate_generations"]) - truth[row["case_id"]]
        for row in estimated
    ]
    covered = sum(
        float(row["ci_lower_generations"]) <= truth[row["case_id"]]
        <= float(row["ci_upper_generations"])
        for row in estimated
    )
    coverage_low, coverage_high = (
        wilson_interval(covered, len(estimated))
        if estimated else (None, None)
    )
    interval_lengths = [
        float(row["ci_upper_generations"])
        - float(row["ci_lower_generations"])
        for row in estimated
    ]
    return {
        "total_count": len(rows),
        "estimated_count": len(estimated),
        "not_evaluable_count": len(rows) - len(estimated),
        **_bootstrap_summary(errors, draws, seed),
        "coverage_numerator": covered,
        "coverage_denominator": len(estimated),
        "coverage": None if not estimated else covered / len(estimated),
        "coverage_ci95_low": coverage_low,
        "coverage_ci95_high": coverage_high,
        "mean_interval_length_generations": (
            None if not interval_lengths
            else float(np.mean(interval_lengths))
        ),
    }


def evaluate_bench(
    *, bench_root: Path, run_root: Path, configuration_path: Path,
    workspace_root: Path, run_id: str,
) -> Path:
    """Ouvre la vérité D25b après fermeture et publie les agrégats."""
    config = load_configuration(configuration_path)
    bench_root = _workspace_path(bench_root, directory=True)
    run_root = _workspace_path(run_root, directory=True)
    bench_manifest_path = _workspace_path(bench_root / "bench_manifest.json")
    run_manifest_path = _workspace_path(run_root / "run_manifest.json")
    estimates_path = _workspace_path(run_root / "estimates.jsonl.gz")
    truth_path = _workspace_path(
        bench_root / "private_truth/gamma_age_truth.jsonl.gz"
    )
    bench_manifest = json.loads(bench_manifest_path.read_text(encoding="utf-8"))
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    if (
        bench_manifest["configuration_sha256"]
        != sha256_file(configuration_path)
        or bench_manifest["protocol_sha256"] != sha256_file(PROTOCOL)
        or bench_manifest["private_truth_sha256"] != sha256_file(truth_path)
        or run_manifest["bench_manifest_sha256"]
        != sha256_file(bench_manifest_path)
        or run_manifest["estimates_sha256"] != sha256_file(estimates_path)
    ):
        raise ValueError("d25b_evaluation_binding_invalid")
    rows = _read_jsonl_gz(estimates_path)
    truth_rows = _read_jsonl_gz(truth_path)
    truth = {
        row["case_id"]: float(row["age_generations"])
        for row in truth_rows
    }
    truth_age = {
        row["case_id"]: int(row["age_generations"])
        for row in truth_rows
    }
    truth_family = {
        row["case_id"]: int(row["family_count"])
        for row in truth_rows
    }
    genealogies = [row["profile"] for row in config["genealogy_profiles"]]

    method_results = []
    for genealogy in genealogies:
        for method in METHODS:
            selected = [
                row for row in rows
                if row["genealogy_profile"] == genealogy
                and row["method"] == method
            ]
            method_results.append({
                "genealogy_profile": genealogy,
                "method": method,
                **_summary(
                    selected, truth,
                    draws=config["bootstrap_draws"],
                    seed=_seed(
                        config["bootstrap_seed"], f"method:{genealogy}:{method}"
                    ),
                ),
            })

    candidate_rows = [
        row for row in rows
        if row["method"] == "TYPE2_CENSORED_INDEPENDENT_V3"
    ]
    independent_rows = [
        row for row in candidate_rows
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    lofo = [
        float(row["lofo_max_absolute_delta_generations"])
        for row in independent_rows
        if row["lofo_max_absolute_delta_generations"] is not None
    ]
    primary_result = {
        "genealogy_profile": "INDEPENDENT",
        "method": "TYPE2_CENSORED_INDEPENDENT_V3",
        **_summary(
            independent_rows, truth,
            draws=config["bootstrap_draws"],
            seed=_seed(config["bootstrap_seed"], "primary"),
        ),
        "lofo_evaluable_count": len(lofo),
        "lofo_mean_max_absolute_delta_generations": (
            None if not lofo else float(np.mean(lofo))
        ),
        "lofo_p95_max_absolute_delta_generations": (
            None if not lofo else float(np.quantile(lofo, 0.95))
        ),
    }

    results_by_age = []
    for genealogy in genealogies:
        for age in config["ages_generations"]:
            selected = [
                row for row in candidate_rows
                if row["genealogy_profile"] == genealogy
                and truth_age[row["case_id"]] == age
            ]
            results_by_age.append({
                "genealogy_profile": genealogy,
                "method": "TYPE2_CENSORED_INDEPENDENT_V3",
                "age_generations": age,
                **_summary(
                    selected, truth,
                    draws=config["bootstrap_draws"],
                    seed=_seed(
                        config["bootstrap_seed"], f"age:{genealogy}:{age}"
                    ),
                ),
            })

    results_by_family = []
    for genealogy in genealogies:
        for family_count in config["family_counts"]:
            selected = [
                row for row in candidate_rows
                if row["genealogy_profile"] == genealogy
                and truth_family[row["case_id"]] == family_count
            ]
            results_by_family.append({
                "genealogy_profile": genealogy,
                "method": "TYPE2_CENSORED_INDEPENDENT_V3",
                "family_count": family_count,
                **_summary(
                    selected, truth,
                    draws=config["bootstrap_draws"],
                    seed=_seed(
                        config["bootstrap_seed"],
                        f"family:{genealogy}:{family_count}",
                    ),
                ),
            })

    fixture_config = config["analytic_fixture"]
    fixture_estimate = estimate_independent_gamma_type2_censored(
        fixture_config["left_lengths_cm"],
        fixture_config["right_lengths_cm"],
        confidence_level=config["confidence_level"],
    )
    analytic_error = abs(
        fixture_estimate.generations
        - fixture_config["expected_unbiased_estimate_generations"]
    )
    analytic_fixture = {
        "event_count": fixture_config["expected_event_count"],
        "total_exposure_morgans": (
            fixture_config["expected_total_exposure_morgans"]
        ),
        "expected_generations": (
            fixture_config["expected_unbiased_estimate_generations"]
        ),
        "observed_generations": fixture_estimate.generations,
        "absolute_error_generations": analytic_error,
    }

    independent_family_results = [
        row for row in results_by_family
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    maximum_family_bias = max(
        abs(float(row["bias"])) for row in independent_family_results
    )
    maximum_family_coverage_deviation = max(
        abs(float(row["coverage"]) - 0.95)
        for row in independent_family_results
    )
    thresholds = config["criteria"]
    maximum_coverage_deviation = max(
        0.95 - thresholds["minimum_family_coverage"],
        thresholds["maximum_family_coverage"] - 0.95,
    )
    criteria = [
        {
            "criterion_id": "CAL-AGE-B01",
            "metric": "analytic_fixture_absolute_error",
            "value": analytic_error, "lower_threshold": None,
            "upper_threshold": thresholds[
                "maximum_analytic_fixture_error_generations"
            ],
            "status": "PASS" if analytic_error <= thresholds[
                "maximum_analytic_fixture_error_generations"
            ] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-B02",
            "metric": "absolute_primary_bias",
            "value": abs(float(primary_result["bias"])),
            "lower_threshold": None,
            "upper_threshold": thresholds[
                "maximum_absolute_primary_bias_generations"
            ],
            "status": "PASS" if abs(float(primary_result["bias"]))
            <= thresholds["maximum_absolute_primary_bias_generations"]
            else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-B03",
            "metric": "primary_coverage",
            "value": primary_result["coverage"],
            "lower_threshold": thresholds["minimum_primary_coverage"],
            "upper_threshold": thresholds["maximum_primary_coverage"],
            "status": "PASS" if thresholds["minimum_primary_coverage"]
            <= float(primary_result["coverage"])
            <= thresholds["maximum_primary_coverage"] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-B04",
            "metric": "maximum_absolute_family_bias",
            "value": maximum_family_bias,
            "lower_threshold": None,
            "upper_threshold": thresholds[
                "maximum_absolute_family_bias_generations"
            ],
            "status": "PASS" if maximum_family_bias <= thresholds[
                "maximum_absolute_family_bias_generations"
            ] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-B05",
            "metric": "maximum_family_coverage_deviation_from_0_95",
            "value": maximum_family_coverage_deviation,
            "lower_threshold": None,
            "upper_threshold": maximum_coverage_deviation,
            "status": "PASS" if maximum_family_coverage_deviation
            <= maximum_coverage_deviation else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-B06",
            "metric": "primary_not_estimated_count",
            "value": primary_result["not_evaluable_count"],
            "lower_threshold": None,
            "upper_threshold": thresholds[
                "maximum_primary_not_estimated_count"
            ],
            "status": "PASS" if primary_result["not_evaluable_count"]
            <= thresholds["maximum_primary_not_estimated_count"]
            else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-B07",
            "metric": "correlated_stress_performance",
            "value": None, "lower_threshold": None, "upper_threshold": None,
            "status": "REVIEW_REQUIRED",
        },
    ]
    technical_status = "PASS" if criteria[0]["status"] == "PASS" else "FAIL"
    independent_pass = all(
        row["status"] == "PASS" for row in criteria[1:6]
    )
    document = {
        "schema_version": "1.0.0",
        "report_id": "CAL-REP-14B",
        "method_id": "d25b_gamma_type2_censoring_evaluation_v1",
        "evaluation_state": "CLOSED",
        "configuration_sha256": sha256_file(configuration_path),
        "protocol_sha256": sha256_file(PROTOCOL),
        "bench_manifest_sha256": sha256_file(bench_manifest_path),
        "run_manifest_sha256": sha256_file(run_manifest_path),
        "technical_status": technical_status,
        "performance_qualification": (
            "PASS_INDEPENDENT_ONLY" if independent_pass else "FAIL"
        ),
        "independent_contract_qualification": (
            "PASS" if independent_pass else "FAIL"
        ),
        "correlated_stress_qualification": "REVIEW_REQUIRED",
        "production_promotion": False,
        "design": {
            "ages_generations": config["ages_generations"],
            "family_counts": config["family_counts"],
            "genealogy_profiles": genealogies,
            "replications_per_cell": config["replications_per_cell"],
            "scenario_count": len(truth_rows),
            "estimate_row_count": len(rows),
            "confidence_level": config["confidence_level"],
            "observation_contract": config["observation_contract"]["id"],
        },
        "analytic_fixture": analytic_fixture,
        "primary_result": primary_result,
        "method_results": method_results,
        "results_by_age": results_by_age,
        "results_by_family_count": results_by_family,
        "criteria": criteria,
        "limits": {
            "phase_error": "NOT_EVALUABLE",
            "marker_resolution": "NOT_EVALUABLE",
            "chance_sharing": "NOT_MODELLED",
            "d24_false_positives": "NOT_ADDRESSED",
        },
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(
        document, "gamma_censoring_evaluation.schema.json", schemas_dir=SCHEMAS
    )
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        output = temporary / "CAL-REP-14B.json"
        write_json(output, document)
        (temporary / "checksums.sha256").write_text(
            f"{sha256_file(output)}  CAL-REP-14B.json\n", encoding="utf-8"
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Banc fermé D25b de censure Gamma"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    materialize = subparsers.add_parser("materialize")
    materialize.add_argument("--configuration-path", type=Path, required=True)
    materialize.add_argument("--workspace-root", type=Path, required=True)
    materialize.add_argument("--run-id", required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--bench-root", type=Path, required=True)
    run.add_argument("--configuration-path", type=Path, required=True)
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
    if command == "materialize":
        materialize_bench(**arguments)
    elif command == "run":
        run_bench(**arguments)
    else:
        evaluate_bench(**arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
