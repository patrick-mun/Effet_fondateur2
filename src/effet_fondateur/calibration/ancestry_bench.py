"""Banc D27 du cœur PCA/projection de l’étape d’ascendance 16A."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import shutil
import tempfile
from typing import Any, Iterable

import numpy as np

from effet_fondateur.ancestry import GenotypePanel, Variant, harmonize_alt_dosages
from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document


REPOSITORY = Path(__file__).resolve().parents[3]
SCHEMAS = REPOSITORY / "calibration/schemas"
PROTOCOL = REPOSITORY / "calibration/protocoles/D27_ANCESTRY_PRESPECIFICATION.md"


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
        raise ValueError("d27_source_outside_workspace")
    return resolved


def load_configuration(path: Path) -> dict[str, Any]:
    """Charge le contrat fermé D27 et interdit toute promotion implicite."""
    document = json.loads(path.read_text(encoding="utf-8"))
    validate_json_document(document, "ancestry_bench_config.schema.json", schemas_dir=SCHEMAS)
    if document["performance_thresholds_sealed"] or document["production_promotion"]:
        raise ValueError("d27_unexpected_scientific_promotion")
    return document


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _read_vcf(
    path: Path, *, ploidy: int, loci: set[tuple[str, int]] | None = None
) -> GenotypePanel:
    samples: tuple[str, ...] | None = None
    variants: list[Variant] = []
    columns: list[list[float]] = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = tuple(line.rstrip("\n").split("\t")[9:])
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chromosome = fields[0].removeprefix("chr")
            position = int(fields[1])
            if loci is not None and (chromosome, position) not in loci:
                continue
            if len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                continue
            if samples is None or len(fields[9:]) != len(samples):
                raise ValueError(f"d27_malformed_vcf:{path.name}")
            genotypes = [value.split(":", 1)[0] for value in fields[9:]]
            if ploidy == 2:
                values = []
                for genotype in genotypes:
                    alleles = genotype.replace("|", "/").split("/")
                    values.append(
                        math.nan if len(alleles) != 2 or "." in alleles
                        else float(sum(allele == "1" for allele in alleles))
                    )
            else:
                values = []
                for genotype in genotypes:
                    alleles = genotype.split("|")
                    if len(alleles) != 2:
                        raise ValueError(f"d27_unphased_local_input:{path.name}")
                    values.extend(math.nan if allele == "." else float(allele) for allele in alleles)
            variants.append(Variant(chromosome, position, fields[2], fields[3], fields[4]))
            columns.append(values)
    if samples is None or not variants:
        raise ValueError(f"d27_empty_vcf:{path.name}")
    if ploidy == 1:
        entity_ids = tuple(f"{sample}:H{copy}" for sample in samples for copy in (1, 2))
    else:
        entity_ids = samples
    return GenotypePanel(
        entity_ids, tuple(variants), np.asarray(columns, dtype=float).T, ploidy
    )


def _loci(panel: GenotypePanel) -> set[tuple[str, int]]:
    return {variant.locus for variant in panel.variants}


def _fit(reference: GenotypePanel, study: GenotypePanel, config: dict[str, Any]):
    return harmonize_alt_dosages(
        reference, study,
        requested_components=config["requested_components"],
        minimum_variants=config["minimum_variants"],
        minimum_reference_call_rate=config["minimum_reference_call_rate"],
    )


def _source_paths(d12_root: Path, d14_root: Path, d17_root: Path, reference_root: Path) -> dict[str, Path]:
    d12 = _workspace_path(d12_root, directory=True)
    d14 = _workspace_path(d14_root, directory=True)
    d17 = _workspace_path(d17_root, directory=True)
    reference = _workspace_path(reference_root, directory=True)
    return {
        "d12_manifest": _workspace_path(d12 / "population_structure_manifest.json"),
        "d12_observed": _workspace_path(d12 / "blinded_observed_panel.vcf.gz"),
        "d12_truth": _workspace_path(d12 / "ancestry_truth.tsv.gz"),
        "d14_manifest": _workspace_path(d14 / "founder_simulation_manifest.json"),
        "d14_carriers": _workspace_path(d14 / "founder_carriers.tsv.gz"),
        "d14_local_truth": _workspace_path(d14 / "local_ancestry_mosaic_truth.tsv.gz"),
        "d17_manifest": _workspace_path(d17 / "measurement_error_manifest.json"),
        "d17_root": d17,
        "reference_manifest": _workspace_path(reference / "preparation_manifest.json"),
        "reference_vcf": _workspace_path(reference / "reference_assets/chr5.reference.vcf.gz"),
        "reference_metadata": _workspace_path(reference / "reference_assets/population_metadata.tsv"),
    }


def _verify_sources(config: dict[str, Any], paths: dict[str, Path]) -> None:
    expected = config["sources"]
    checks = {
        "d12_manifest_sha256": paths["d12_manifest"],
        "d14_manifest_sha256": paths["d14_manifest"],
        "d17_manifest_sha256": paths["d17_manifest"],
        "reference_preparation_manifest_sha256": paths["reference_manifest"],
        "reference_vcf_sha256": paths["reference_vcf"],
        "reference_metadata_sha256": paths["reference_metadata"],
    }
    for key, path in checks.items():
        if sha256_file(path) != expected[key]:
            raise ValueError(f"d27_source_hash_mismatch:{key}")
    d12_manifest = json.loads(paths["d12_manifest"].read_text(encoding="utf-8"))
    d12_artifacts = {row["role"]: row for row in d12_manifest["artifacts"]}
    for role, key in (("BLINDED_OBSERVED_POPULATION_PANEL", "d12_observed"), ("ANCESTRY_TRUTH", "d12_truth")):
        if d12_artifacts[role]["sha256"] != sha256_file(paths[key]):
            raise ValueError(f"d27_d12_artifact_hash_mismatch:{role}")
    d14_manifest = json.loads(paths["d14_manifest"].read_text(encoding="utf-8"))
    d14_artifacts = {row["role"]: row for row in d14_manifest["artifacts"]}
    for role, key in (("FOUNDER_CARRIER_TRUTH", "d14_carriers"), ("LOCAL_ANCESTRY_MOSAIC_TRUTH", "d14_local_truth")):
        if d14_artifacts[role]["sha256"] != sha256_file(paths[key]):
            raise ValueError(f"d27_d14_artifact_hash_mismatch:{role}")
    reference_manifest = json.loads(paths["reference_manifest"].read_text(encoding="utf-8"))
    if (
        reference_manifest["configuration"]["reference_rule"]
        != "EXACT_D10_PHASING_REFERENCE_EXCLUDE_BOTH_DONOR_PARTITIONS"
        or reference_manifest["counts"]["reference_samples"] != config["reference_entity_count"]
        or reference_manifest["counts"]["calibration_donors"] != 900
    ):
        raise ValueError("d27_reference_independence_not_demonstrated")
    d17_manifest = json.loads(paths["d17_manifest"].read_text(encoding="utf-8"))
    artifacts = {(row["role"], row.get("profile_id")): row for row in d17_manifest["artifacts"]}
    for profile in config["profiles"]:
        for role, filename in (
            ("OBSERVED_MEASUREMENT_PANEL", "observed.vcf.gz"),
            ("CANDIDATE_PHASED_MEASUREMENT_PANEL", "candidate.phased.vcf.gz"),
        ):
            path = _workspace_path(paths["d17_root"] / f"profiles/{profile}/{filename}")
            if artifacts[(role, profile)]["sha256"] != sha256_file(path):
                raise ValueError(f"d27_profile_hash_mismatch:{profile}:{role}")


def _save_result(payload: dict[str, np.ndarray], prefix: str, result: Any) -> None:
    payload[f"{prefix}_study_ids"] = np.asarray(result.study_sample_ids)
    payload[f"{prefix}_reference_ids"] = np.asarray(result.reference_sample_ids)
    payload[f"{prefix}_study_scores"] = result.study_scores
    payload[f"{prefix}_reference_scores"] = result.model.reference_scores
    payload[f"{prefix}_eigenvalues"] = result.model.eigenvalues
    payload[f"{prefix}_informative_variant_count"] = np.asarray([len(result.variants)])


def run_calls(
    *, configuration_path: Path, d12_root: Path, d14_root: Path, d17_root: Path,
    reference_root: Path, workspace_root: Path, run_id: str,
) -> Path:
    """Projette les entrées aveugles et scelle les coordonnées avant vérité."""
    config = load_configuration(configuration_path)
    paths = _source_paths(d12_root, d14_root, d17_root, reference_root)
    _verify_sources(config, paths)
    d12_study = _read_vcf(paths["d12_observed"], ploidy=2)
    control_local_path = _workspace_path(paths["d17_root"] / "profiles/control/candidate.phased.vcf.gz")
    control_local = _read_vcf(control_local_path, ploidy=1)
    selected_loci = _loci(d12_study) | _loci(control_local)
    reference_global = _read_vcf(paths["reference_vcf"], ploidy=2, loci=selected_loci)
    reference_local = _read_vcf(paths["reference_vcf"], ploidy=1, loci=selected_loci)
    if len(reference_global.sample_ids) != config["reference_entity_count"]:
        raise ValueError("d27_reference_sample_count_mismatch")

    arrays: dict[str, np.ndarray] = {}
    _save_result(arrays, "mixture", _fit(reference_global, d12_study, config))
    for profile in config["profiles"]:
        profile_root = paths["d17_root"] / "profiles" / profile
        global_study = _read_vcf(_workspace_path(profile_root / "observed.vcf.gz"), ploidy=2)
        local_study = _read_vcf(_workspace_path(profile_root / "candidate.phased.vcf.gz"), ploidy=1)
        _save_result(arrays, f"global_{profile}", _fit(reference_global, global_study, config))
        _save_result(arrays, f"local_{profile}", _fit(reference_local, local_study, config))

    destination = new_workspace_run(workspace_root.absolute(), run_id)
    temporary = Path(tempfile.mkdtemp(prefix=f".{run_id}.", dir=workspace_root))
    try:
        coordinates = temporary / "sealed_coordinates.npz"
        np.savez_compressed(coordinates, **arrays)
        manifest = {
            "schema_version": "1.0.0", "run_id": run_id,
            "method_id": config["method_id"], "call_state": "SEALED_BEFORE_TRUTH",
            "truth_opened": False, "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL),
            "reference_partition_role": config["reference_partition_role"],
            "reference_entity_count": config["reference_entity_count"],
            "profiles": config["profiles"],
            "scope": "CHR5_REGIONAL_PROXY_AND_LOCAL_CONDITIONAL_ON_CANDIDATE_PHASE",
            "source_hashes": config["sources"],
            "artifact": {"path": coordinates.name, "sha256": sha256_file(coordinates), "size_bytes": coordinates.stat().st_size},
            "checks": ["REFERENCE_ONLY_AXES", "D10_REFERENCE_DONOR_PARTITIONS_DISJOINT", "NO_TRUTH_READ_DURING_CALLS", "NO_CACHE_OR_REAL_DATA"],
        }
        write_json(temporary / "calls_manifest.json", manifest)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return destination


def _metadata(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        return {row["SampleID"]: row["Population"] for row in reader}


def _centroid(scores: np.ndarray, ids: Iterable[str], labels: dict[str, str], group: str) -> np.ndarray:
    indexes = [index for index, entity in enumerate(ids) if labels[entity.split(":", 1)[0]] == group]
    if not indexes:
        raise ValueError(f"d27_reference_group_missing:{group}")
    return np.mean(scores[indexes], axis=0)


def _whiten(scores: np.ndarray, eigenvalues: np.ndarray) -> np.ndarray:
    return scores / np.sqrt(eigenvalues)


def _line_proportion(scores: np.ndarray, first: np.ndarray, second: np.ndarray) -> np.ndarray:
    axis = first - second
    if float(axis @ axis) == 0:
        raise ValueError("d27_identical_reference_centroids")
    return np.clip((scores - second) @ axis / float(axis @ axis), 0.0, 1.0)


def _align(
    profile_reference: np.ndarray, profile_study: np.ndarray, profile_eigenvalues: np.ndarray,
    control_reference: np.ndarray,
) -> np.ndarray:
    profile_reference_z = _whiten(profile_reference, profile_eigenvalues)
    profile_study_z = _whiten(profile_study, profile_eigenvalues)
    left, _, right = np.linalg.svd(profile_reference_z.T @ control_reference, full_matrices=False)
    return profile_study_z @ (left @ right)


def _distances(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean((first - second) ** 2, axis=1))


def _local_sample_distances(ids: list[str], aligned: np.ndarray, control: np.ndarray) -> tuple[np.ndarray, int]:
    if len(ids) % 2:
        raise ValueError("d27_local_entity_order_invalid")
    distances: list[float] = []
    swaps = 0
    for index in range(0, len(ids), 2):
        if ids[index].split(":", 1)[0] != ids[index + 1].split(":", 1)[0]:
            raise ValueError("d27_local_haplotype_pair_invalid")
        direct = _distances(aligned[index:index + 2], control[index:index + 2])
        swapped = _distances(aligned[index:index + 2], control[index:index + 2][::-1])
        if float(np.sum(swapped)) < float(np.sum(direct)):
            distances.append(float(np.mean(swapped)))
            swaps += 1
        else:
            distances.append(float(np.mean(direct)))
    return np.asarray(distances), swaps


def _summary(profile: str, scope: str, values: np.ndarray, informative: int, swaps: int | None = None) -> dict[str, Any]:
    return {
        "profile_id": profile, "scope": scope,
        "entity_count": int(values.size), "informative_variant_count": informative,
        "mean_standardized_distance": float(np.mean(values)),
        "median_standardized_distance": float(np.median(values)),
        "p95_standardized_distance": float(np.quantile(values, 0.95)),
        "maximum_standardized_distance": float(np.max(values)),
        "optimal_haplotype_swap_count": swaps,
    }


def _truth_rows(path: Path) -> list[dict[str, str]]:
    with _open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evaluate(
    *, configuration_path: Path, calls_root: Path, d12_root: Path, d14_root: Path,
    d17_root: Path, reference_root: Path, output_root: Path,
) -> Path:
    """Ouvre les vérités après scellement et publie l’évaluation agrégée."""
    config = load_configuration(configuration_path)
    paths = _source_paths(d12_root, d14_root, d17_root, reference_root)
    _verify_sources(config, paths)
    calls = _workspace_path(calls_root, directory=True)
    manifest_path = _workspace_path(calls / "calls_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    coordinates_path = _workspace_path(calls / manifest["artifact"]["path"])
    if (
        manifest["call_state"] != "SEALED_BEFORE_TRUTH" or manifest["truth_opened"]
        or manifest["configuration_sha256"] != sha256_file(configuration_path)
        or manifest["artifact"]["sha256"] != sha256_file(coordinates_path)
    ):
        raise ValueError("d27_calls_not_validly_sealed")
    data = np.load(coordinates_path, allow_pickle=False)
    labels = _metadata(paths["reference_metadata"])

    mixture_reference = _whiten(data["mixture_reference_scores"], data["mixture_eigenvalues"])
    mixture_study = _whiten(data["mixture_study_scores"], data["mixture_eigenvalues"])
    reference_ids = data["mixture_reference_ids"].tolist()
    first_source = config["mixture_strata"]["REFERENCE_STRATUM_A"]
    second_source = config["mixture_strata"]["REFERENCE_STRATUM_B"]
    first_centroid = _centroid(mixture_reference, reference_ids, labels, first_source)
    second_centroid = _centroid(mixture_reference, reference_ids, labels, second_source)
    estimates = _line_proportion(mixture_study, first_centroid, second_centroid)
    estimates_by_id = dict(zip(data["mixture_study_ids"].tolist(), estimates, strict=True))
    mixture_rows = []
    for row in _truth_rows(paths["d12_truth"]):
        proportions = json.loads(row["mixture_proportions"])
        if set(proportions) <= {first_source, second_source}:
            truth = float(proportions.get(first_source, 0.0))
            estimate = float(estimates_by_id[row["sample_id"]])
            mixture_rows.append((row["sample_id"], truth, estimate, abs(estimate - truth)))
    if not mixture_rows or not {truth for _, truth, _, _ in mixture_rows} >= {0.0, 0.5, 1.0}:
        raise ValueError("d27_two_way_mixture_truth_incomplete")
    errors = np.asarray([row[3] for row in mixture_rows])
    generator = np.random.default_rng(config["bootstrap_seed"])
    bootstrap = np.mean(
        errors[generator.integers(0, len(errors), size=(config["bootstrap_replicates"], len(errors)))], axis=1
    )
    by_truth = []
    for truth in sorted({row[1] for row in mixture_rows}):
        selected = [row for row in mixture_rows if row[1] == truth]
        by_truth.append({
            "true_stratum_a_proportion": truth, "sample_count": len(selected),
            "estimated_mean": float(np.mean([row[2] for row in selected])),
            "mean_absolute_error": float(np.mean([row[3] for row in selected])),
        })

    stability: list[dict[str, Any]] = []
    for scope in ("global", "local"):
        prefix = f"{scope}_control"
        control_ids = data[f"{prefix}_study_ids"].tolist()
        control_reference = _whiten(data[f"{prefix}_reference_scores"], data[f"{prefix}_eigenvalues"])
        control_study = _whiten(data[f"{prefix}_study_scores"], data[f"{prefix}_eigenvalues"])
        for profile in config["profiles"]:
            profile_prefix = f"{scope}_{profile}"
            if data[f"{profile_prefix}_study_ids"].tolist() != control_ids:
                raise ValueError(f"d27_profile_entity_mismatch:{scope}:{profile}")
            aligned = _align(
                data[f"{profile_prefix}_reference_scores"], data[f"{profile_prefix}_study_scores"],
                data[f"{profile_prefix}_eigenvalues"], control_reference,
            )
            if scope == "local":
                distances, swaps = _local_sample_distances(control_ids, aligned, control_study)
            else:
                distances, swaps = _distances(aligned, control_study), None
            stability.append(_summary(
                profile, scope.upper(), distances,
                int(data[f"{profile_prefix}_informative_variant_count"][0]), swaps,
            ))

    local_reference = _whiten(data["local_control_reference_scores"], data["local_control_eigenvalues"])
    local_study = _whiten(data["local_control_study_scores"], data["local_control_eigenvalues"])
    local_ids = data["local_control_study_ids"].tolist()
    background = _centroid(local_reference, data["local_control_reference_ids"].tolist(), labels, first_source)
    founder = _centroid(local_reference, data["local_control_reference_ids"].tolist(), labels, config["local_founder_stratum"])
    carrier_rows = _truth_rows(paths["d14_carriers"])
    local_truth_rows = _truth_rows(paths["d14_local_truth"])
    if (
        {row["background_stratum_id"] for row in local_truth_rows} != {first_source}
        or {row["founder_stratum_id"] for row in local_truth_rows} != {config["local_founder_stratum"]}
    ):
        raise ValueError("d27_local_ancestry_truth_not_prespecified")
    carrier_proportions = []
    for row in carrier_rows:
        entity = f'{row["sample_id"]}:H{row["carrier_haplotype_index"]}'
        index = local_ids.index(entity)
        estimate = float(_line_proportion(local_study[index:index + 1], founder, background)[0])
        carrier_proportions.append({"entity_role": "FOUNDER_CARRIER_COPY", "founder_reference_axis_proportion": estimate})

    primary_instability = max(
        row["mean_standardized_distance"] for row in stability if row["profile_id"] != "control"
    )
    document = {
        "schema_version": "1.0.0", "evaluation_state": "CLOSED",
        "report_id": "CAL-REP-16A", "method_id": "d27_ancestry_evaluation_v1",
        "configuration_sha256": sha256_file(configuration_path),
        "protocol_sha256": sha256_file(PROTOCOL), "calls_manifest_sha256": sha256_file(manifest_path),
        "reference_role": config["reference_partition_role"],
        "reference_entity_count": config["reference_entity_count"],
        "analysis_scope": "CHR5_REGIONAL_PROXY_NOT_GENOMEWIDE",
        "stage_execution_scope": "16A_CORE_PCA_HARMONIZATION_PROJECTION_NOT_END_TO_END",
        "profiles": config["profiles"], "component_count": config["requested_components"],
        "mixture": {
            "stratum_a_alias": "REFERENCE_STRATUM_A", "stratum_b_alias": "REFERENCE_STRATUM_B",
            "sample_count": len(mixture_rows), "mean_absolute_error": float(np.mean(errors)),
            "median_absolute_error": float(np.median(errors)), "p95_absolute_error": float(np.quantile(errors, 0.95)),
            "bootstrap_ci95_low": float(np.quantile(bootstrap, 0.025)),
            "bootstrap_ci95_high": float(np.quantile(bootstrap, 0.975)),
            "results_by_true_proportion": by_truth,
        },
        "stability_by_profile": stability,
        "repeated_position_instability": primary_instability,
        "local_founder_check": {
            "carrier_copy_count": len(carrier_proportions), "independent_founder_haplotype_count": 1,
            "mean_founder_reference_axis_proportion": float(np.mean([row["founder_reference_axis_proportion"] for row in carrier_proportions])),
            "closer_to_founder_reference_count": sum(row["founder_reference_axis_proportion"] > 0.5 for row in carrier_proportions),
        },
        "criteria": [
            {"criterion_id": "CAL-ANC-001", "metric": "simulated_mixture_proportion_mean_absolute_error", "value": float(np.mean(errors)), "threshold": None, "status": "REVIEW_REQUIRED"},
            {"criterion_id": "CAL-ANC-002", "metric": "repeated_position_instability", "value": primary_instability, "threshold": None, "status": "REVIEW_REQUIRED"},
        ],
        "technical_status": "PASS_CORE_EVALUATOR",
        "performance_qualification": "REVIEW_REQUIRED_UNSEALED",
        "production_promotion": False, "calibrated": False, "founder_effect_proven": False,
        "checks": ["SOURCE_HASHES_MATCHED", "REFERENCE_DONOR_PARTITIONS_DISJOINT", "CALLS_SEALED_BEFORE_TRUTH", "REFERENCE_ONLY_AXIS_ALIGNMENT", "EIGHT_PAIRED_PROFILES_EVALUATED", "NO_GEOGRAPHIC_OR_ETHNIC_ASSIGNMENT", "NO_CACHE_OR_REAL_DATA"],
        "limits": [
            "Le volet diploïde utilise seulement 487 marqueurs régionaux chr5 et ne qualifie pas la PCA génomique.",
            "Le volet local est conditionnel à la phase candidate D17 et ne qualifie pas le parcours SHAPEIT5 puis 16A bout-en-bout.",
            "Les centroïdes sont des repères de simulation larges, pas des identités ethniques ou géographiques.",
            "Les six copies locales partagent un seul haplotype fondateur et ne sont pas six réplications indépendantes.",
            "Les seuils scientifiques n’étaient pas scellés avant D27 et ne peuvent pas être ajustés après résultat."
        ],
    }
    validate_json_document(document, "ancestry_evaluation.schema.json", schemas_dir=SCHEMAS)
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        result_path = temporary / "CAL-REP-16A.json"
        write_json(result_path, document)
        evaluation_manifest = {
            "schema_version": "1.0.0", "evaluation_state": "CLOSED",
            "configuration_sha256": sha256_file(configuration_path),
            "protocol_sha256": sha256_file(PROTOCOL), "calls_manifest_sha256": sha256_file(manifest_path),
            "artifacts": [{"path": result_path.name, "sha256": sha256_file(result_path), "size_bytes": result_path.stat().st_size}],
        }
        write_json(temporary / "evaluation_manifest.json", evaluation_manifest)
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
    common.add_argument("--d12-root", type=Path, required=True)
    common.add_argument("--d14-root", type=Path, required=True)
    common.add_argument("--d17-root", type=Path, required=True)
    common.add_argument("--reference-root", type=Path, required=True)
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
