"""Banc D25a : fidélité R et domaine d'observation de la méthode Gamma."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import shutil
import subprocess
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
    chance_sharing_correction_cm,
    estimate_correlated_gamma,
    estimate_correlated_gamma_r_faithful,
    estimate_independent_gamma,
    estimate_independent_gamma_r_faithful,
)


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
PROTOCOL = REPOSITORY / "calibration/protocoles/D25A_GAMMA_AUDIT_PRESPECIFICATION.md"


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
        raise ValueError("d25a_source_outside_workspace")
    return resolved


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge D25a et contrôle les invariants croisés de la préspecification."""
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(document, "gamma_audit_config.schema.json", schemas_dir=SCHEMAS)
    genealogy_names = [row["profile"] for row in document["genealogy_profiles"]]
    fixture_names = [row["fixture_id"] for row in document["r_reference_fixtures"]]
    if len(genealogy_names) != len(set(genealogy_names)):
        raise ValueError("d25a_genealogy_profile_duplicates")
    if len(fixture_names) != len(set(fixture_names)):
        raise ValueError("d25a_r_fixture_duplicates")
    if document["simulation_seed"] == document["bootstrap_seed"]:
        raise ValueError("d25a_seeds_not_independent")
    for fixture in document["r_reference_fixtures"]:
        if len(fixture["left_lengths_cm"]) != len(fixture["right_lengths_cm"]):
            raise ValueError("d25a_r_fixture_arm_count_mismatch")
    return document


def _write_jsonl_gz(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            for row in rows:
                payload = json.dumps(row, sort_keys=True, separators=(",", ":"))
                compressed.write((payload + "\n").encode())


def _read_jsonl_gz(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _case_seed(
    master_seed: int, genealogy_profile: str, age: int,
    family_count: int, replication: int,
) -> int:
    key = f"{master_seed}:{genealogy_profile}:{age}:{family_count}:{replication}"
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
        return [generator.expovariate(total_rate_per_cm) for _ in range(family_count)]
    shared_rate = total_rate_per_cm * shared_path_fraction
    private_rate = total_rate_per_cm * (1.0 - shared_path_fraction)
    shared_breakpoint = generator.expovariate(shared_rate)
    return [
        min(shared_breakpoint, generator.expovariate(private_rate))
        for _ in range(family_count)
    ]


def _pairwise_shared_observation(lengths_cm: list[float]) -> list[float]:
    """Tronque l'unique maximum au deuxième maximum observable par partage."""
    if len(lengths_cm) < 2:
        raise ValueError("d25a_pairwise_observation_requires_two_arms")
    observed = list(lengths_cm)
    longest_index = max(range(len(observed)), key=observed.__getitem__)
    second_longest = sorted(observed)[-2]
    observed[longest_index] = second_longest
    return observed


def materialize_bench(
    *, configuration_path: Path, workspace_root: Path, run_id: str,
) -> Path:
    """Matérialise 18 000 scénarios nouveaux avec vérité d'âge séparée."""
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
                    for replication in range(1, config["replications_per_cell"] + 1):
                        seed = _case_seed(
                            config["simulation_seed"], genealogy_profile,
                            age, family_count, replication,
                        )
                        generator = random.Random(seed)
                        left_full = _simulate_arm_lengths_cm(
                            generator, age_generations=age,
                            family_count=family_count,
                            shared_path_fraction=shared_fraction,
                        )
                        right_full = _simulate_arm_lengths_cm(
                            generator, age_generations=age,
                            family_count=family_count,
                            shared_path_fraction=shared_fraction,
                        )
                        case_id = "D25A_" + hashlib.sha256(str(seed).encode()).hexdigest()[:20]
                        blinded.append({
                            "case_id": case_id,
                            "genealogy_profile": genealogy_profile,
                            "family_count": family_count,
                            "full_left_lengths_cm": left_full,
                            "full_right_lengths_cm": right_full,
                            "truncated_left_lengths_cm": _pairwise_shared_observation(left_full),
                            "truncated_right_lengths_cm": _pairwise_shared_observation(right_full),
                        })
                        truth.append({
                            "case_id": case_id,
                            "age_generations": age,
                            "genealogy_profile": genealogy_profile,
                            "shared_path_fraction": shared_fraction,
                            "family_count": family_count,
                            "replication": replication,
                            "case_seed_sha256": hashlib.sha256(str(seed).encode()).hexdigest(),
                        })
        blinded_path = temporary / "blinded_inputs/gamma_lengths.jsonl.gz"
        truth_path = temporary / "private_truth/gamma_age_truth.jsonl.gz"
        _write_jsonl_gz(blinded_path, blinded)
        _write_jsonl_gz(truth_path, truth)
        manifest = {
            "schema_version": "1.0.0",
            "method_id": "d25a_correlated_path_gamma_grid_v1",
            "state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "case_count": len(blinded),
            "blinded_inputs_sha256": sha256_file(blinded_path),
            "private_truth_sha256": sha256_file(truth_path),
            "truth_separated": True,
            "d25_seed_reused": False,
            "calibrated": False,
            "founder_effect_proven": False,
        }
        write_json(temporary / "bench_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _estimate_result(
    method: str,
    left_lengths_cm: list[float],
    right_lengths_cm: list[float],
    confidence_level: float,
    *,
    boundary_correction: bool,
) -> dict[str, Any]:
    try:
        if method == "LEGACY_INDEPENDENT_V1":
            estimate = estimate_independent_gamma(
                left_lengths_cm, right_lengths_cm,
                confidence_level=confidence_level,
            )
        elif method == "LEGACY_CORRELATED_V1":
            estimate = estimate_correlated_gamma(
                left_lengths_cm, right_lengths_cm,
                confidence_level=confidence_level,
            )
        elif method in {"EXACT_INDEPENDENT_ORACLE_V1", "R_FAITHFUL_INDEPENDENT_V2"}:
            estimate = estimate_independent_gamma_r_faithful(
                left_lengths_cm, right_lengths_cm,
                confidence_level=confidence_level,
                apply_boundary_correction=boundary_correction,
            )
        elif method == "R_FAITHFUL_CORRELATED_V2":
            estimate = estimate_correlated_gamma_r_faithful(
                left_lengths_cm, right_lengths_cm,
                confidence_level=confidence_level,
                apply_boundary_correction=boundary_correction,
            )
        else:
            raise ValueError("d25a_method_invalid")
    except GammaAgeError as error:
        return {
            "status": "NOT_ESTIMATED", "exclusion_code": str(error),
            "estimate_generations": None, "ci_lower_generations": None,
            "ci_upper_generations": None, "rho": None,
            "point_effective_n": None, "interval_effective_n": None,
        }
    return {
        "status": "ESTIMATED", "exclusion_code": None,
        "estimate_generations": estimate.generations,
        "ci_lower_generations": estimate.confidence_lower,
        "ci_upper_generations": estimate.confidence_upper,
        "rho": estimate.correlation,
        "point_effective_n": estimate.point_effective_sample_size,
        "interval_effective_n": estimate.effective_sample_size,
    }


def _methods_for_observation(observation_profile: str) -> tuple[str, ...]:
    if observation_profile == "FULL_BREAKPOINTS":
        return (
            "LEGACY_INDEPENDENT_V1", "LEGACY_CORRELATED_V1",
            "EXACT_INDEPENDENT_ORACLE_V1", "R_FAITHFUL_CORRELATED_V2",
        )
    if observation_profile == "PAIRWISE_SHARED_TRUNCATED":
        return (
            "LEGACY_INDEPENDENT_V1", "LEGACY_CORRELATED_V1",
            "R_FAITHFUL_INDEPENDENT_V2", "R_FAITHFUL_CORRELATED_V2",
        )
    raise ValueError("d25a_observation_profile_invalid")


def _primary_method(genealogy_profile: str, observation_profile: str) -> str:
    if genealogy_profile == "INDEPENDENT":
        return (
            "EXACT_INDEPENDENT_ORACLE_V1"
            if observation_profile == "FULL_BREAKPOINTS"
            else "R_FAITHFUL_INDEPENDENT_V2"
        )
    return "R_FAITHFUL_CORRELATED_V2"


def _lofo_max_delta(
    *,
    genealogy_profile: str,
    observation_profile: str,
    full_left: list[float],
    full_right: list[float],
    full_estimate: float,
    confidence_level: float,
) -> float | None:
    if len(full_left) - 1 < 3:
        return None
    method = _primary_method(genealogy_profile, observation_profile)
    deltas: list[float] = []
    for omitted in range(len(full_left)):
        left = full_left[:omitted] + full_left[omitted + 1:]
        right = full_right[:omitted] + full_right[omitted + 1:]
        if observation_profile == "PAIRWISE_SHARED_TRUNCATED":
            left = _pairwise_shared_observation(left)
            right = _pairwise_shared_observation(right)
        result = _estimate_result(
            method, left, right, confidence_level,
            boundary_correction=observation_profile == "PAIRWISE_SHARED_TRUNCATED",
        )
        if result["estimate_generations"] is None:
            return None
        deltas.append(abs(float(result["estimate_generations"]) - full_estimate))
    return max(deltas)


def run_bench(
    *, bench_root: Path, configuration_path: Path,
    workspace_root: Path, run_id: str,
) -> Path:
    """Ferme les estimations candidates et historiques avant ouverture des âges."""
    config = load_configuration(configuration_path)
    bench_root = _workspace_path(bench_root, directory=True)
    manifest_path = _workspace_path(bench_root / "bench_manifest.json")
    blinded_path = _workspace_path(bench_root / "blinded_inputs/gamma_lengths.jsonl.gz")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest["configuration_sha256"] != sha256_file(configuration_path)
        or manifest["blinded_inputs_sha256"] != sha256_file(blinded_path)
    ):
        raise ValueError("d25a_bench_binding_invalid")
    blinded = _read_jsonl_gz(blinded_path)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    rows: list[dict[str, Any]] = []
    try:
        for case in blinded:
            full_left = case["full_left_lengths_cm"]
            full_right = case["full_right_lengths_cm"]
            for observation_profile in config["observation_profiles"]:
                if observation_profile == "FULL_BREAKPOINTS":
                    left, right = full_left, full_right
                    boundary_correction = False
                else:
                    left = case["truncated_left_lengths_cm"]
                    right = case["truncated_right_lengths_cm"]
                    boundary_correction = True
                primary_method = _primary_method(
                    case["genealogy_profile"], observation_profile
                )
                for method in _methods_for_observation(observation_profile):
                    result = _estimate_result(
                        method, left, right, config["confidence_level"],
                        boundary_correction=boundary_correction,
                    )
                    lofo = None
                    if method == primary_method and result["estimate_generations"] is not None:
                        lofo = _lofo_max_delta(
                            genealogy_profile=case["genealogy_profile"],
                            observation_profile=observation_profile,
                            full_left=full_left,
                            full_right=full_right,
                            full_estimate=float(result["estimate_generations"]),
                            confidence_level=config["confidence_level"],
                        )
                    rows.append({
                        "case_id": case["case_id"],
                        "genealogy_profile": case["genealogy_profile"],
                        "observation_profile": observation_profile,
                        "method": method,
                        "family_count": case["family_count"],
                        "boundary_correction_applied": (
                            boundary_correction
                            if method.startswith("R_FAITHFUL")
                            or method == "EXACT_INDEPENDENT_ORACLE_V1"
                            else True
                        ),
                        "is_contract_primary": method == primary_method,
                        "lofo_max_absolute_delta_generations": lofo,
                        **result,
                    })
        estimates_path = temporary / "estimates.jsonl.gz"
        _write_jsonl_gz(estimates_path, rows)
        run_manifest = {
            "schema_version": "1.0.0",
            "method_id": "d25a_gamma_method_matrix_v1",
            "state": "CLOSED_BEFORE_TRUTH",
            "configuration_sha256": sha256_file(configuration_path),
            "bench_manifest_sha256": sha256_file(manifest_path),
            "blinded_inputs_sha256": sha256_file(blinded_path),
            "estimate_row_count": len(rows),
            "estimates_sha256": sha256_file(estimates_path),
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


def _r_vector(values: list[float]) -> str:
    return "c(" + ",".join(format(value, ".17g") for value in values) + ")"


def _r_wrapper(config: dict[str, Any], official_r_script: Path) -> str:
    script_literal = str(official_r_script).replace("\\", "\\\\").replace('"', '\\"')
    lines = ["options(warn=2)"]
    for fixture in config["r_reference_fixtures"]:
        chance = fixture["chance_sharing"]
        lines.extend((
            f'chance.sharing.correction <- {"TRUE" if chance["enabled"] else "FALSE"}',
            f'median.allele.frequency <- {chance["median_allele_frequency"]}',
            f'markers.on.chromosome <- {chance["markers_on_chromosome"]}',
            f'length.of.chromosome <- {chance["chromosome_length_cm"]}',
            f'confidence.coefficient <- {config["confidence_level"]}',
            f'l.lengths <- {_r_vector(fixture["left_lengths_cm"])}',
            f'r.lengths <- {_r_vector(fixture["right_lengths_cm"])}',
            f'invisible(capture.output(source("{script_literal}")))',
            "cat(paste(" + ",".join((
                f'"{fixture["fixture_id"]}"',
                "format(i.tau.hat, digits=17)", "format(i.l, digits=17)",
                "format(i.u, digits=17)", "format(c.tau.hat, digits=17)",
                "format(c.l, digits=17)", "format(c.u, digits=17)",
                "format(rho.hat, digits=17)", "format(n.star, digits=17)",
            )) + ', sep="\\t"), "\\n")',
        ))
    return "\n".join(lines) + "\n"


def _rho_region(rho: float, sample_size: int) -> str:
    if rho < -2.0 / (sample_size - 1):
        return "I1"
    if rho < -1.0 / (sample_size - 1):
        return "I2"
    if rho < 0:
        return "I3"
    return "NONNEGATIVE"


def verify_r_reference(
    *, configuration_path: Path, official_r_script: Path,
    rscript_executable: str, workspace_root: Path, run_id: str,
) -> Path:
    """Exécute le script R officiel et compare les cinq fixtures à Python v2."""
    config = load_configuration(configuration_path)
    official_r_script = _workspace_path(official_r_script)
    if sha256_file(official_r_script) != config["official_r_script_sha256"]:
        raise ValueError("d25a_official_r_checksum_invalid")
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        wrapper = temporary / "verify_reference.R"
        wrapper.write_text(_r_wrapper(config, official_r_script), encoding="utf-8")
        completed = subprocess.run(
            [rscript_executable, str(wrapper)], check=True, capture_output=True,
            text=True, timeout=60,
        )
        version = subprocess.run(
            [rscript_executable, "--version"], check=True, capture_output=True,
            text=True, timeout=10,
        )
        r_rows: dict[str, list[float]] = {}
        for line in completed.stdout.splitlines():
            fields = line.strip().split("\t")
            if len(fields) == 9:
                r_rows[fields[0]] = [float(value) for value in fields[1:]]
        fixture_results = []
        maximum_error = 0.0
        for fixture in config["r_reference_fixtures"]:
            fixture_id = fixture["fixture_id"]
            if fixture_id not in r_rows:
                raise ValueError("d25a_r_fixture_output_missing")
            chance = fixture["chance_sharing"]
            correction = (
                chance_sharing_correction_cm(
                    chance["median_allele_frequency"],
                    chance["markers_on_chromosome"],
                    chance["chromosome_length_cm"],
                ) if chance["enabled"] else 0.0
            )
            independent = estimate_independent_gamma_r_faithful(
                fixture["left_lengths_cm"], fixture["right_lengths_cm"],
                chance_sharing_correction_cm=correction,
            )
            correlated = estimate_correlated_gamma_r_faithful(
                fixture["left_lengths_cm"], fixture["right_lengths_cm"],
                chance_sharing_correction_cm=correction,
            )
            python_values = [
                independent.generations, independent.confidence_lower,
                independent.confidence_upper, correlated.generations,
                correlated.confidence_lower, correlated.confidence_upper,
                float(correlated.correlation), correlated.effective_sample_size,
            ]
            r_values = r_rows[fixture_id]
            errors = [abs(left - right) for left, right in zip(python_values, r_values)]
            maximum_error = max(maximum_error, max(errors))
            observed_region = _rho_region(
                float(correlated.correlation), len(fixture["left_lengths_cm"])
            )
            fixture_results.append({
                "fixture_id": fixture_id,
                "expected_rho_region": fixture["expected_rho_region"],
                "observed_rho_region": observed_region,
                "rho": correlated.correlation,
                "maximum_absolute_error": max(errors),
                "region_matches": observed_region == fixture["expected_rho_region"],
            })
        document = {
            "schema_version": "1.0.0",
            "state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "official_r_script_sha256": sha256_file(official_r_script),
            "r_version": (version.stdout or version.stderr).strip(),
            "fixture_count": len(fixture_results),
            "maximum_absolute_error_generations": maximum_error,
            "all_regions_match": all(row["region_matches"] for row in fixture_results),
            "fixtures": fixture_results,
        }
        output = temporary / "r_parity.json"
        write_json(output, document)
        (temporary / "checksums.sha256").write_text(
            f"{sha256_file(output)}  r_parity.json\n", encoding="utf-8"
        )
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
    errors = [float(row["estimate_generations"]) - truth[row["case_id"]] for row in estimated]
    covered = sum(
        float(row["ci_lower_generations"]) <= truth[row["case_id"]]
        <= float(row["ci_upper_generations"])
        for row in estimated
    )
    low, high = wilson_interval(covered, len(estimated)) if estimated else (None, None)
    interval_lengths = [
        float(row["ci_upper_generations"]) - float(row["ci_lower_generations"])
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
        "coverage_ci95_low": low,
        "coverage_ci95_high": high,
        "mean_interval_length_generations": (
            None if not interval_lengths else float(np.mean(interval_lengths))
        ),
    }


def _contract_primary_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if row["is_contract_primary"]]


def evaluate_bench(
    *, bench_root: Path, run_root: Path, r_parity_root: Path,
    configuration_path: Path, workspace_root: Path, run_id: str,
) -> Path:
    """Ouvre les âges après fermeture et publie les seuls agrégats D25a."""
    config = load_configuration(configuration_path)
    bench_root = _workspace_path(bench_root, directory=True)
    run_root = _workspace_path(run_root, directory=True)
    r_parity_root = _workspace_path(r_parity_root, directory=True)
    bench_manifest_path = _workspace_path(bench_root / "bench_manifest.json")
    run_manifest_path = _workspace_path(run_root / "run_manifest.json")
    estimates_path = _workspace_path(run_root / "estimates.jsonl.gz")
    truth_path = _workspace_path(bench_root / "private_truth/gamma_age_truth.jsonl.gz")
    parity_path = _workspace_path(r_parity_root / "r_parity.json")
    bench_manifest = json.loads(bench_manifest_path.read_text(encoding="utf-8"))
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    parity = json.loads(parity_path.read_text(encoding="utf-8"))
    if (
        run_manifest["configuration_sha256"] != sha256_file(configuration_path)
        or run_manifest["bench_manifest_sha256"] != sha256_file(bench_manifest_path)
        or run_manifest["estimates_sha256"] != sha256_file(estimates_path)
        or bench_manifest["private_truth_sha256"] != sha256_file(truth_path)
        or parity["configuration_sha256"] != sha256_file(configuration_path)
    ):
        raise ValueError("d25a_evaluation_binding_invalid")
    rows = _read_jsonl_gz(estimates_path)
    truth_rows = _read_jsonl_gz(truth_path)
    truth = {row["case_id"]: float(row["age_generations"]) for row in truth_rows}
    truth_age = {row["case_id"]: int(row["age_generations"]) for row in truth_rows}
    truth_family = {row["case_id"]: int(row["family_count"]) for row in truth_rows}
    primary_rows = _contract_primary_rows(rows)

    contract_results = []
    for genealogy in [row["profile"] for row in config["genealogy_profiles"]]:
        for observation in config["observation_profiles"]:
            selected = [
                row for row in primary_rows
                if row["genealogy_profile"] == genealogy
                and row["observation_profile"] == observation
            ]
            label = f"contract:{genealogy}:{observation}"
            lofo = [
                float(row["lofo_max_absolute_delta_generations"])
                for row in selected
                if row["lofo_max_absolute_delta_generations"] is not None
            ]
            contract_results.append({
                "genealogy_profile": genealogy,
                "observation_profile": observation,
                "method": _primary_method(genealogy, observation),
                **_summary(
                    selected, truth, draws=config["bootstrap_draws"],
                    seed=_seed(config["bootstrap_seed"], label),
                ),
                "lofo_evaluable_count": len(lofo),
                "lofo_mean_max_absolute_delta_generations": (
                    None if not lofo else float(np.mean(lofo))
                ),
                "lofo_p95_max_absolute_delta_generations": (
                    None if not lofo else float(np.quantile(lofo, 0.95))
                ),
            })

    method_results = []
    for genealogy in [row["profile"] for row in config["genealogy_profiles"]]:
        for observation in config["observation_profiles"]:
            for method in _methods_for_observation(observation):
                selected = [
                    row for row in rows
                    if row["genealogy_profile"] == genealogy
                    and row["observation_profile"] == observation
                    and row["method"] == method
                ]
                label = f"method:{genealogy}:{observation}:{method}"
                method_results.append({
                    "genealogy_profile": genealogy,
                    "observation_profile": observation,
                    "method": method,
                    **_summary(
                        selected, truth, draws=config["bootstrap_draws"],
                        seed=_seed(config["bootstrap_seed"], label),
                    ),
                })

    results_by_age = []
    for contract in contract_results:
        for age in config["ages_generations"]:
            selected = [
                row for row in primary_rows
                if row["genealogy_profile"] == contract["genealogy_profile"]
                and row["observation_profile"] == contract["observation_profile"]
                and truth_age[row["case_id"]] == age
            ]
            label = (
                f"age:{contract['genealogy_profile']}:"
                f"{contract['observation_profile']}:{age}"
            )
            results_by_age.append({
                "genealogy_profile": contract["genealogy_profile"],
                "observation_profile": contract["observation_profile"],
                "method": contract["method"], "age_generations": age,
                **_summary(
                    selected, truth, draws=config["bootstrap_draws"],
                    seed=_seed(config["bootstrap_seed"], label),
                ),
            })

    results_by_family = []
    for contract in contract_results:
        for family_count in config["family_counts"]:
            selected = [
                row for row in primary_rows
                if row["genealogy_profile"] == contract["genealogy_profile"]
                and row["observation_profile"] == contract["observation_profile"]
                and truth_family[row["case_id"]] == family_count
            ]
            label = (
                f"family:{contract['genealogy_profile']}:"
                f"{contract['observation_profile']}:{family_count}"
            )
            results_by_family.append({
                "genealogy_profile": contract["genealogy_profile"],
                "observation_profile": contract["observation_profile"],
                "method": contract["method"], "family_count": family_count,
                **_summary(
                    selected, truth, draws=config["bootstrap_draws"],
                    seed=_seed(config["bootstrap_seed"], label),
                ),
            })

    rho_counts: list[dict[str, Any]] = []
    for genealogy in [row["profile"] for row in config["genealogy_profiles"]]:
        for observation in config["observation_profiles"]:
            selected = [
                row for row in rows
                if row["genealogy_profile"] == genealogy
                and row["observation_profile"] == observation
                and row["method"] == "R_FAITHFUL_CORRELATED_V2"
                and row["rho"] is not None
            ]
            counts = {region: 0 for region in ("I1", "I2", "I3", "NONNEGATIVE")}
            for row in selected:
                counts[_rho_region(float(row["rho"]), int(row["family_count"]))] += 1
            rho_counts.append({
                "genealogy_profile": genealogy,
                "observation_profile": observation,
                "total_count": len(selected),
                "counts": counts,
            })

    independent_full = next(
        row for row in contract_results
        if row["genealogy_profile"] == "INDEPENDENT"
        and row["observation_profile"] == "FULL_BREAKPOINTS"
    )
    independent_truncated = next(
        row for row in contract_results
        if row["genealogy_profile"] == "INDEPENDENT"
        and row["observation_profile"] == "PAIRWISE_SHARED_TRUNCATED"
    )
    thresholds = config["criteria"]
    criteria = [
        {
            "criterion_id": "CAL-AGE-A01", "metric": "maximum_r_python_error",
            "value": parity["maximum_absolute_error_generations"],
            "lower_threshold": None,
            "upper_threshold": thresholds["maximum_r_parity_error_generations"],
            "status": "PASS" if (
                parity["maximum_absolute_error_generations"]
                <= thresholds["maximum_r_parity_error_generations"]
                and parity["all_regions_match"]
            ) else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-A02", "metric": "absolute_oracle_bias",
            "value": abs(independent_full["bias"]), "lower_threshold": None,
            "upper_threshold": thresholds["maximum_absolute_oracle_bias_generations"],
            "status": "PASS" if abs(independent_full["bias"]) <= thresholds["maximum_absolute_oracle_bias_generations"] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-A03", "metric": "oracle_coverage",
            "value": independent_full["coverage"],
            "lower_threshold": thresholds["minimum_oracle_coverage"],
            "upper_threshold": thresholds["maximum_oracle_coverage"],
            "status": "PASS" if thresholds["minimum_oracle_coverage"] <= independent_full["coverage"] <= thresholds["maximum_oracle_coverage"] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-A04", "metric": "absolute_truncated_independent_bias",
            "value": abs(independent_truncated["bias"]), "lower_threshold": None,
            "upper_threshold": thresholds["maximum_absolute_truncated_independent_bias_generations"],
            "status": "PASS" if abs(independent_truncated["bias"]) <= thresholds["maximum_absolute_truncated_independent_bias_generations"] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-A05", "metric": "truncated_independent_coverage",
            "value": independent_truncated["coverage"],
            "lower_threshold": thresholds["minimum_truncated_independent_coverage"],
            "upper_threshold": thresholds["maximum_truncated_independent_coverage"],
            "status": "PASS" if thresholds["minimum_truncated_independent_coverage"] <= independent_truncated["coverage"] <= thresholds["maximum_truncated_independent_coverage"] else "FAIL",
        },
        {
            "criterion_id": "CAL-AGE-A06", "metric": "correlated_performance",
            "value": None, "lower_threshold": None, "upper_threshold": None,
            "status": "REVIEW_REQUIRED",
        },
    ]
    technical_failure = criteria[0]["status"] == "FAIL"
    performance_failure = any(
        row["status"] == "FAIL" for row in criteria[1:5]
    )
    document = {
        "schema_version": "1.0.0",
        "report_id": "CAL-REP-14A",
        "method_id": "d25a_gamma_domain_audit_v1",
        "evaluation_state": "CLOSED",
        "configuration_sha256": sha256_file(configuration_path),
        "protocol_sha256": sha256_file(PROTOCOL),
        "bench_manifest_sha256": sha256_file(bench_manifest_path),
        "run_manifest_sha256": sha256_file(run_manifest_path),
        "r_parity_sha256": sha256_file(parity_path),
        "technical_status": "FAIL" if technical_failure else "PASS",
        "performance_qualification": (
            "FAIL" if performance_failure else "REVIEW_REQUIRED"
        ),
        "design": {
            "ages_generations": config["ages_generations"],
            "family_counts": config["family_counts"],
            "genealogy_profiles": [row["profile"] for row in config["genealogy_profiles"]],
            "observation_profiles": config["observation_profiles"],
            "replications_per_cell": config["replications_per_cell"],
            "independent_scenario_count": len(truth),
            "estimate_row_count": len(rows),
            "confidence_level": config["confidence_level"],
        },
        "r_parity": parity,
        "contract_results": contract_results,
        "method_results": method_results,
        "results_by_age": results_by_age,
        "results_by_family_count": results_by_family,
        "rho_region_counts": rho_counts,
        "criteria": criteria,
        "phase_error_evaluation": {
            "status": "NOT_EVALUABLE",
            "reason": "D25A_LENGTH_MODEL_ONLY_NO_MATCHED_SHAPEIT5_OUTPUT",
        },
        "calibrated": False,
        "founder_effect_proven": False,
    }
    validate_json_document(document, "gamma_audit_evaluation.schema.json", schemas_dir=SCHEMAS)
    destination = new_workspace_run(workspace_root, run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        output = temporary / "CAL-REP-14A.json"
        write_json(output, document)
        (temporary / "checksums.sha256").write_text(
            f"{sha256_file(output)}  CAL-REP-14A.json\n", encoding="utf-8"
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description="Banc fermé D25a de datation Gamma")
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
    parity = subparsers.add_parser("verify-r")
    parity.add_argument("--configuration-path", type=Path, required=True)
    parity.add_argument("--official-r-script", type=Path, required=True)
    parity.add_argument("--rscript-executable", default="Rscript")
    parity.add_argument("--workspace-root", type=Path, required=True)
    parity.add_argument("--run-id", required=True)
    evaluate = subparsers.add_parser("evaluate")
    evaluate.add_argument("--bench-root", type=Path, required=True)
    evaluate.add_argument("--run-root", type=Path, required=True)
    evaluate.add_argument("--r-parity-root", type=Path, required=True)
    evaluate.add_argument("--configuration-path", type=Path, required=True)
    evaluate.add_argument("--workspace-root", type=Path, required=True)
    evaluate.add_argument("--run-id", required=True)
    arguments = vars(parser.parse_args())
    command = arguments.pop("command")
    function = {
        "materialize": materialize_bench,
        "run": run_bench,
        "verify-r": verify_r_reference,
        "evaluate": evaluate_bench,
    }[command]
    result = function(**arguments)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
