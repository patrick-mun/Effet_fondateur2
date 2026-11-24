"""Matérialisation déterministe des shards scientifiques aveugles R13.

Le matérialiseur est le seul composant autorisé à ouvrir les graines privées.
Il écrit séparément une vue aveugle consommable par les outils et une vérité
privée destinée uniquement à l'évaluation post-clôture. Le worker d'appel ne
reçoit jamais le chemin de la vérité.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
from hashlib import sha256
import json
from math import exp
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
from typing import Any, Mapping, Sequence

import numpy as np

from effet_fondateur.calibration.contracts import (
    canonical_sha256,
    file_sha256,
    load_json,
    load_yaml,
)
from effet_fondateur.calibration.r13_scientific_runner import (
    RELEASE_ID,
    build_external_call_schedule,
)


CONFIGURATION_PATH = Path("calibration/configurations/scientific_materialization_v1.yaml")
MATERIALIZER_STATE = "FULL_CALL_MATERIALIZATION_READY"
DENSITY_STEPS = {"D_MAX": 2, "D_MID": 4, "D_LOW": 12}


class R13MaterializationError(RuntimeError):
    """Contrat, entrée privée ou artefact matérialisé invalide."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R13MaterializationError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _safe_relative(path: Path, root: Path, code: str) -> str:
    try:
        return path.resolve(strict=True).relative_to(root.resolve(strict=True)).as_posix()
    except ValueError as error:
        raise R13MaterializationError(code, str(path)) from error


def _artifact(path: Path, root: Path) -> dict[str, Any]:
    return {
        "relative_path": _safe_relative(path, root, "materialized_path_outside_root"),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def _staged_artifact(path: Path, final_path: Path, root: Path) -> dict[str, Any]:
    """Empreinte un fichier staging en publiant son futur chemin final."""
    final_root = root.resolve(strict=True)
    try:
        relative = final_path.resolve(strict=False).relative_to(final_root).as_posix()
    except ValueError as error:
        raise R13MaterializationError(
            "materialized_path_outside_root", str(final_path)
        ) from error
    return {
        "relative_path": relative,
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
    }


def load_materializer_configuration(repository_root: Path) -> dict[str, Any]:
    configuration = load_yaml(repository_root.resolve(strict=True) / CONFIGURATION_PATH)
    if configuration.get("release_id") != RELEASE_ID:
        _fail("materializer_release_mismatch")
    if configuration.get("synthetic_only") is not True:
        _fail("materializer_synthetic_only_required")
    return configuration


def _chromosome_records(configuration: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    records = configuration.get("chromosomes")
    if not isinstance(records, list) or len(records) != 22:
        _fail("materializer_chromosome_contract_invalid")
    result: dict[int, dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, Mapping):
            _fail("materializer_chromosome_contract_invalid")
        chromosome = record.get("chromosome")
        if (
            isinstance(chromosome, bool)
            or not isinstance(chromosome, int)
            or chromosome not in range(1, 23)
            or chromosome in result
            or not isinstance(record.get("length_bp"), int)
            or not isinstance(record.get("length_cm"), (int, float))
        ):
            _fail("materializer_chromosome_contract_invalid")
        result[chromosome] = dict(record)
    if set(result) != set(range(1, 23)):
        _fail("materializer_chromosome_coverage_invalid")
    return result


def dense_variant_allocation(configuration: Mapping[str, Any]) -> dict[int, int]:
    """Alloue 1,2 M variants par unités de 12 pour garantir l'emboîtement."""
    chromosomes = _chromosome_records(configuration)
    dense_total = configuration.get("model", {}).get("dense_autosomal_variants")
    if dense_total != 1_200_000 or dense_total % 12:
        _fail("materializer_dense_total_invalid")
    total_units = dense_total // 12
    minimum_units = 1_000
    remaining = total_units - 22 * minimum_units
    total_bp = sum(record["length_bp"] for record in chromosomes.values())
    raw = {
        chromosome: remaining * record["length_bp"] / total_bp
        for chromosome, record in chromosomes.items()
    }
    units = {chromosome: minimum_units + int(value) for chromosome, value in raw.items()}
    missing = total_units - sum(units.values())
    order = sorted(
        chromosomes,
        key=lambda chromosome: (raw[chromosome] - int(raw[chromosome]), -chromosome),
        reverse=True,
    )
    for chromosome in order[:missing]:
        units[chromosome] += 1
    allocation = {chromosome: count * 12 for chromosome, count in units.items()}
    if sum(allocation.values()) != dense_total:
        _fail("materializer_dense_allocation_invalid")
    return allocation


def resolve_shapeit_chromosome(
    call_contract: Mapping[str, Any], configuration: Mapping[str, Any]
) -> int:
    value = call_contract.get("autosome")
    if isinstance(value, int) and not isinstance(value, bool) and value in range(1, 23):
        return value
    sentinels = configuration.get("sentinels")
    if not isinstance(value, str) or not isinstance(sentinels, Mapping):
        _fail("materializer_shapeit_autosome_invalid", str(value))
    chromosome = sentinels.get(value)
    if not isinstance(chromosome, int) or chromosome not in range(1, 23):
        _fail("materializer_shapeit_sentinel_invalid", value)
    return chromosome


def _validate_call_contract(call_contract: Mapping[str, Any]) -> str:
    call_id = call_contract.get("call_id")
    if not isinstance(call_id, str) or not call_id:
        _fail("materializer_call_id_invalid")
    schedule = build_external_call_schedule()
    matches = [
        call
        for family in ("shapeit", "ibd")
        for call in schedule[family]
        if call["call_id"] == call_id
    ]
    if len(matches) != 1 or dict(call_contract) != matches[0]:
        _fail("materializer_call_contract_mismatch", call_id)
    return "shapeit" if call_id.startswith("SHAPEIT-") else "ibd"


def materializer_contract_audit(repository_root: Path) -> dict[str, Any]:
    configuration = load_materializer_configuration(repository_root)
    allocation = dense_variant_allocation(configuration)
    observed_totals = {
        density: sum(count // step for count in allocation.values())
        for density, step in DENSITY_STEPS.items()
    }
    expected = {"D_MAX": 600_000, "D_MID": 300_000, "D_LOW": 100_000}
    if observed_totals != expected:
        _fail("materializer_density_totals_invalid")
    schedule = build_external_call_schedule()
    for call in schedule["shapeit"]:
        resolve_shapeit_chromosome(call, configuration)
    for call in schedule["ibd"]:
        if call.get("density") not in DENSITY_STEPS or call.get("pair_count") != 100:
            _fail("materializer_ibd_contract_invalid", str(call.get("call_id")))
    return {
        "state": MATERIALIZER_STATE,
        "materializer_id": configuration["materializer_id"],
        "configuration_sha256": canonical_sha256(configuration),
        "dense_variant_count": sum(allocation.values()),
        "observed_variant_totals": observed_totals,
        "shapeit_call_count": len(schedule["shapeit"]),
        "ibd_call_count": len(schedule["ibd"]),
        "supported_call_count": len(schedule["shapeit"]) + len(schedule["ibd"]),
        "sentinels": dict(configuration["sentinels"]),
        "truth_separated_from_worker": True,
    }


def _truth_seed(release_root: Path, truth_id: str) -> str:
    manifest = load_json(release_root / "private_truth/master_truth_manifest.json")
    records = [
        record
        for record in manifest.get("truth_artifacts", [])
        if isinstance(record, Mapping) and record.get("truth_id") == truth_id
    ]
    if len(records) != 1:
        _fail("materializer_truth_binding_missing", truth_id)
    relative = records[0].get("relative_path")
    if not isinstance(relative, str):
        _fail("materializer_truth_path_invalid", truth_id)
    path = release_root / relative
    truth = load_json(path)
    seed = truth.get("generator_seed_sha256")
    if not isinstance(seed, str) or len(seed) != 64:
        _fail("materializer_truth_seed_invalid", truth_id)
    return seed


def _validate_materialization_authorization(
    repository_root: Path,
    release_root: Path,
) -> None:
    receipt = load_json(repository_root / "calibration/readiness/R12C_G5_READY.json")
    frozen_path = release_root / "governance/r12c_frozen_manifest.json"
    checksums_path = release_root / "release_checksums.sha256"
    frozen = load_json(frozen_path)
    evidence = receipt.get("evidence_sha256")
    if (
        receipt.get("gate_g5_passed") is not True
        or receipt.get("materializer_state") != MATERIALIZER_STATE
        or receipt.get("single_call_worker_state")
        != "FULL_SINGLE_CALL_EXECUTION_READY"
        or not isinstance(evidence, Mapping)
        or evidence.get("r12c_frozen_manifest") != file_sha256(frozen_path)
        or evidence.get("release_checksums") != file_sha256(checksums_path)
        or frozen.get("materializer_state") != MATERIALIZER_STATE
    ):
        _fail("materializer_g5_authorization_invalid")
    for line in checksums_path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        digest, relative = line.split("  ", 1)
        path = release_root / relative
        if not path.is_file() or file_sha256(path) != digest:
            _fail("materializer_release_checksum_mismatch", relative)
    for binding in frozen.get("frozen_components", []):
        relative = binding.get("relative_path") if isinstance(binding, Mapping) else None
        if (
            not isinstance(relative, str)
            or not (repository_root / relative).is_file()
            or file_sha256(repository_root / relative) != binding.get("sha256")
        ):
            _fail("materializer_frozen_component_drift", str(relative))


def _rng(seed_text: str, namespace: str) -> np.random.Generator:
    digest = sha256(f"{seed_text}|{namespace}".encode()).digest()
    seed = int.from_bytes(digest[:16], "big")
    return np.random.Generator(np.random.PCG64DXSM(seed))


def _population_weights(count: int) -> np.ndarray:
    states = np.asarray(
        (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.25, 0.75, 0.0),
            (0.50, 0.50, 0.0),
            (0.75, 0.25, 0.0),
            (0.0, 0.25, 0.75),
            (0.0, 0.50, 0.50),
            (0.0, 0.75, 0.25),
        ),
        dtype=np.float64,
    )
    return states[np.arange(count) % len(states)]


def _allele_frequencies(dense_rank: int) -> np.ndarray:
    base = 0.01 + ((dense_rank * 7_919) % 4_700) / 10_000
    shifts = np.asarray((0.0, 0.08, -0.06), dtype=np.float64)
    return np.clip(base + shifts, 0.01, 0.49)


def _update_haplotypes(
    rng: np.random.Generator,
    previous: np.ndarray | None,
    probabilities: np.ndarray,
    refresh_probability: float,
) -> np.ndarray:
    draws = rng.random((len(probabilities), 2)) < probabilities[:, None]
    if previous is None:
        return draws.astype(np.uint8)
    refresh = rng.random((len(probabilities), 2)) < refresh_probability
    return np.where(refresh, draws, previous).astype(np.uint8)


def _vcf_header(
    *,
    source: str,
    chromosome: int,
    length_bp: int,
    sample_ids: Sequence[str],
) -> str:
    return (
        "##fileformat=VCFv4.2\n"
        f"##source={source}\n"
        "##reference=GRCh38_SYNTHETIC_COORDINATE_FRAME\n"
        f"##contig=<ID={chromosome},length={length_bp}>\n"
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_ids)
        + "\n"
    )


def _genotype_text(haplotypes: np.ndarray, separator: str) -> str:
    lookup = np.asarray(
        (f"0{separator}0", f"0{separator}1", f"1{separator}0", f"1{separator}1")
    )
    codes = haplotypes[:, 0] * 2 + haplotypes[:, 1]
    return "\t".join(lookup[codes].tolist())


def _compress_and_index_vcf(source: Path, output: Path, bcftools: Path) -> None:
    completed = subprocess.run(
        [str(bcftools), "view", "-Oz", "-o", str(output), str(source)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        _fail("materializer_bcftools_compression_failed", completed.stderr[-1000:])
    indexed = subprocess.run(
        [str(bcftools), "index", "--tbi", str(output)],
        capture_output=True,
        text=True,
        check=False,
    )
    if indexed.returncode != 0:
        _fail("materializer_bcftools_index_failed", indexed.stderr[-1000:])


def _shapeit_paths(run_root: Path, chromosome: int, density: str) -> Path:
    return run_root / "blinded_materialized/shapeit" / f"chr{chromosome:02d}" / density


def _shapeit_manifest_path(run_root: Path, chromosome: int, density: str) -> Path:
    return _shapeit_paths(run_root, chromosome, density) / "manifest.json"


def _validate_artifacts(root: Path, artifacts: Any) -> None:
    if not isinstance(artifacts, list) or not artifacts:
        _fail("materializer_artifact_manifest_invalid")
    for artifact in artifacts:
        if not isinstance(artifact, Mapping):
            _fail("materializer_artifact_manifest_invalid")
        relative = artifact.get("relative_path")
        if not isinstance(relative, str):
            _fail("materializer_artifact_path_invalid")
        pure = PurePosixPath(relative)
        if pure.is_absolute() or ".." in pure.parts:
            _fail("materializer_artifact_path_invalid")
        path = root / relative
        if (
            not path.is_file()
            or path.stat().st_size != artifact.get("size_bytes")
            or file_sha256(path) != artifact.get("sha256")
        ):
            _fail("materializer_artifact_integrity_mismatch", relative)


@dataclass(frozen=True)
class SmokeDimensions:
    """Dimensions réduites réservées aux tests unitaires explicites."""

    study_count: int
    reference_count: int
    complete_trios: int
    dense_variant_count: int


def materialize_shapeit_shard(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    truth_root: Path,
    call_contract: Mapping[str, Any],
    *,
    test_dimensions: SmokeDimensions | None = None,
) -> dict[str, Any]:
    """Crée une entrée SHAPEIT partagée entre toutes les graines d'une cellule."""
    if _validate_call_contract(call_contract) != "shapeit":
        _fail("materializer_shapeit_call_required")
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    run_root = run_root.resolve(strict=True)
    truth_root = truth_root.resolve(strict=False)
    configuration = load_materializer_configuration(repository_root)
    if test_dimensions is None:
        _validate_materialization_authorization(repository_root, release_root)
    chromosome = resolve_shapeit_chromosome(call_contract, configuration)
    density = str(call_contract["density"])
    manifest_path = _shapeit_manifest_path(run_root, chromosome, density)
    if manifest_path.is_file():
        manifest = load_json(manifest_path)
        _validate_artifacts(run_root, manifest.get("artifacts"))
        return manifest
    chromosome_record = _chromosome_records(configuration)[chromosome]
    allocation = dense_variant_allocation(configuration)
    production_dense_count = allocation[chromosome]
    model = configuration["model"]
    if test_dimensions is None:
        study_count = int(model["study_identities"])
        reference_count = int(model["phase_reference_identities"])
        trio_count = int(model["complete_trios"])
        dense_count = production_dense_count
    else:
        study_count = test_dimensions.study_count
        reference_count = test_dimensions.reference_count
        trio_count = test_dimensions.complete_trios
        dense_count = test_dimensions.dense_variant_count
        if dense_count % 12 or study_count < trio_count * 3 or min(
            study_count, reference_count, dense_count
        ) < 1:
            _fail("materializer_test_dimensions_invalid")
    registry = load_json(
        release_root / "blinded_inputs/VALIDATION_BLIND/study_registry.json"
    )
    registry_ids = [row.get("sample_id") for row in registry.get("subjects", [])]
    if len(registry_ids) < study_count or any(not isinstance(value, str) for value in registry_ids):
        _fail("materializer_study_registry_insufficient")
    study_ids = [str(value) for value in registry_ids[:study_count]]
    reference_ids = [f"SYN-R13-REF-{index:05d}" for index in range(1, reference_count + 1)]
    seed_text = _truth_seed(release_root, "T-PHASE")
    rng = _rng(seed_text, f"shapeit|chr{chromosome}")
    step = DENSITY_STEPS[density]
    expected_variants = dense_count // step
    blind_root = _shapeit_paths(run_root, chromosome, density)
    staging = blind_root.with_name(f".{blind_root.name}.building")
    if staging.exists() or blind_root.exists():
        _fail("materializer_shapeit_incomplete_or_colliding", str(blind_root))
    staging.mkdir(parents=True)
    private_root = truth_root / "shapeit" / f"chr{chromosome:02d}" / density
    if private_root.exists():
        _fail("materializer_private_truth_collision", str(private_root))
    private_root.mkdir(parents=True)
    study_source = staging / "study.source.vcf.gz"
    reference_source = staging / "reference.source.vcf.gz"
    truth_source = private_root / "study_truth.source.vcf.gz"
    pedigree_path = staging / "pedigree.tsv"
    map_path = staging / "genetic_map.tsv.gz"
    father_indices = np.arange(trio_count)
    mother_indices = np.arange(trio_count, 2 * trio_count)
    child_indices = np.arange(2 * trio_count, 3 * trio_count)
    transmissions = rng.integers(0, 2, size=(trio_count, 2), dtype=np.uint8)
    pedigree_path.write_text(
        "".join(
            f"{study_ids[child]}\t{study_ids[father]}\t{study_ids[mother]}\n"
            for father, mother, child in zip(
                father_indices, mother_indices, child_indices, strict=True
            )
        ),
        encoding="utf-8",
    )
    study_weights = _population_weights(study_count)
    reference_weights = _population_weights(reference_count)
    previous_study: np.ndarray | None = None
    previous_reference: np.ndarray | None = None
    selected = 0
    length_bp = int(chromosome_record["length_bp"])
    length_cm = float(chromosome_record["length_cm"])
    delta_cm = length_cm / (dense_count + 1)
    refresh = 1.0 - exp(-delta_cm / float(model["ld_refresh_scale_cm"]))
    recombination = 1.0 - exp(-delta_cm / 100.0)
    with (
        gzip.open(study_source, "wt", encoding="utf-8") as study_handle,
        gzip.open(reference_source, "wt", encoding="utf-8") as reference_handle,
        gzip.open(truth_source, "wt", encoding="utf-8") as truth_handle,
        gzip.open(map_path, "wt", encoding="utf-8") as map_handle,
    ):
        study_handle.write(
            _vcf_header(
                source="effet_fondateur_R13_scientific_blind",
                chromosome=chromosome,
                length_bp=length_bp,
                sample_ids=study_ids,
            )
        )
        truth_handle.write(
            _vcf_header(
                source="effet_fondateur_R13_private_phase_truth",
                chromosome=chromosome,
                length_bp=length_bp,
                sample_ids=study_ids,
            )
        )
        reference_handle.write(
            _vcf_header(
                source="effet_fondateur_R13_synthetic_reference",
                chromosome=chromosome,
                length_bp=length_bp,
                sample_ids=reference_ids,
            )
        )
        map_handle.write("pos\tchr\tcM\n")
        for dense_rank in range(1, dense_count + 1):
            frequencies = _allele_frequencies(dense_rank)
            study_probabilities = study_weights @ frequencies
            reference_probabilities = reference_weights @ frequencies
            previous_study = _update_haplotypes(
                rng, previous_study, study_probabilities, refresh
            )
            previous_reference = _update_haplotypes(
                rng, previous_reference, reference_probabilities, refresh
            )
            if trio_count:
                switches = rng.random((trio_count, 2)) < recombination
                transmissions = np.bitwise_xor(transmissions, switches.astype(np.uint8))
                previous_study[child_indices, 0] = previous_study[
                    father_indices, transmissions[:, 0]
                ]
                previous_study[child_indices, 1] = previous_study[
                    mother_indices, transmissions[:, 1]
                ]
            if dense_rank % step:
                continue
            selected += 1
            position = max(1, int(dense_rank * length_bp / (dense_count + 1)))
            position_cm = dense_rank * length_cm / (dense_count + 1)
            variant_id = f"R13-CHR{chromosome:02d}-{dense_rank:07d}"
            study_ac = int(previous_study.sum())
            reference_ac = int(previous_reference.sum())
            prefix = f"{chromosome}\t{position}\t{variant_id}\tA\tG\t.\tPASS\t"
            study_handle.write(
                prefix
                + f"AC={study_ac};AN={2 * study_count}\tGT\t"
                + _genotype_text(previous_study, "/")
                + "\n"
            )
            truth_handle.write(
                prefix
                + f"AC={study_ac};AN={2 * study_count}\tGT\t"
                + _genotype_text(previous_study, "|")
                + "\n"
            )
            reference_handle.write(
                prefix
                + f"AC={reference_ac};AN={2 * reference_count}\tGT\t"
                + _genotype_text(previous_reference, "|")
                + "\n"
            )
            map_handle.write(f"{position}\t{chromosome}\t{position_cm:.8f}\n")
    if selected != expected_variants:
        _fail("materializer_shapeit_variant_count_mismatch")
    bcftools_name = shutil.which("bcftools")
    if bcftools_name is None:
        _fail("materializer_bcftools_missing")
    bcftools = Path(bcftools_name).resolve(strict=True)
    study_vcf = staging / "study.vcf.gz"
    reference_vcf = staging / "reference.vcf.gz"
    truth_vcf = private_root / "study_truth.vcf.gz"
    _compress_and_index_vcf(study_source, study_vcf, bcftools)
    _compress_and_index_vcf(reference_source, reference_vcf, bcftools)
    _compress_and_index_vcf(truth_source, truth_vcf, bcftools)
    study_source.unlink()
    reference_source.unlink()
    truth_source.unlink()
    blind_artifacts = [
        study_vcf,
        Path(f"{study_vcf}.tbi"),
        reference_vcf,
        Path(f"{reference_vcf}.tbi"),
        pedigree_path,
        map_path,
    ]
    private_manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "truth_id": "T-PHASE",
        "chromosome": chromosome,
        "density": density,
        "variant_count": selected,
        "study_count": study_count,
        "artifacts": [
            _artifact(truth_vcf, truth_root),
            _artifact(Path(f"{truth_vcf}.tbi"), truth_root),
        ],
    }
    _write_json(private_root / "manifest.json", private_manifest)
    manifest = {
        "schema_version": "1.0.0",
        "materializer_id": configuration["materializer_id"],
        "release_id": RELEASE_ID,
        "partition_id": "VALIDATION_BLIND",
        "family": "SHAPEIT",
        "chromosome": chromosome,
        "density": density,
        "variant_count": selected,
        "study_count": study_count,
        "reference_count": reference_count,
        "complete_trios": trio_count,
        "synthetic_only": True,
        "truth_available_to_worker": False,
        "source_truth_binding_sha256": sha256(seed_text.encode()).hexdigest(),
        "artifacts": [
            _staged_artifact(path, blind_root / path.name, run_root)
            for path in blind_artifacts
        ],
    }
    _write_json(staging / "manifest.json", manifest)
    staging.rename(blind_root)
    return load_json(manifest_path)


def _ibd_paths(run_root: Path, world_id: str, density: str) -> Path:
    return run_root / "blinded_materialized/ibd" / world_id / density


def materialize_ibd_world(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    truth_root: Path,
    call_contract: Mapping[str, Any],
    *,
    test_dense_variant_count: int | None = None,
) -> dict[str, Any]:
    """Crée un monde IBD partagé par les deux outils pour une densité."""
    if _validate_call_contract(call_contract) != "ibd":
        _fail("materializer_ibd_call_required")
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    run_root = run_root.resolve(strict=True)
    truth_root = truth_root.resolve(strict=False)
    configuration = load_materializer_configuration(repository_root)
    if test_dense_variant_count is None:
        _validate_materialization_authorization(repository_root, release_root)
    world_id = str(call_contract["world_id"])
    density = str(call_contract["density"])
    world_class = str(call_contract["world_class"])
    chromosome = (
        int(configuration["ibd"]["other_autosome"])
        if call_contract.get("location") == "OTHER_AUTOSOME"
        else int(configuration["ibd"]["default_autosome"])
    )
    manifest_path = _ibd_paths(run_root, world_id, density) / "manifest.json"
    if manifest_path.is_file():
        manifest = load_json(manifest_path)
        _validate_artifacts(run_root, manifest.get("artifacts"))
        return manifest
    allocation = dense_variant_allocation(configuration)
    dense_count = test_dense_variant_count or allocation[chromosome]
    if dense_count % 12:
        _fail("materializer_ibd_dense_count_invalid")
    step = DENSITY_STEPS[density]
    sample_count = int(configuration["ibd"]["samples_per_world"])
    pair_count = int(call_contract["pair_count"])
    if sample_count != pair_count * 2:
        _fail("materializer_ibd_pair_contract_invalid")
    root = _ibd_paths(run_root, world_id, density)
    staging = root.with_name(f".{root.name}.building")
    if staging.exists() or root.exists():
        _fail("materializer_ibd_incomplete_or_colliding", str(root))
    staging.mkdir(parents=True)
    private_root = truth_root / "ibd" / world_id / density
    if private_root.exists():
        _fail("materializer_private_truth_collision", str(private_root))
    private_root.mkdir(parents=True)
    seed_text = _truth_seed(release_root, "T-IBD")
    rng = _rng(seed_text, f"ibd|{world_id}|chr{chromosome}")
    chromosome_record = _chromosome_records(configuration)[chromosome]
    length_bp = int(chromosome_record["length_bp"])
    length_cm = float(chromosome_record["length_cm"])
    delta_cm = length_cm / (dense_count + 1)
    refresh = 1.0 - exp(
        -delta_cm / float(configuration["model"]["ld_refresh_scale_cm"])
    )
    sample_ids = [f"SYN-R13-{world_id}-S{index:03d}" for index in range(1, sample_count + 1)]
    source = staging / "phased.source.vcf.gz"
    map_path = staging / "ibd.map"
    truth_segments: list[dict[str, Any]] = []
    length = float(call_contract.get("length_cm", 0.0))
    positive = world_class != "PRIMARY_NEGATIVE"
    center = 0.75 * length_cm if call_contract.get("location") == "OFF_TARGET_SAME_AUTOSOME" else 0.50 * length_cm
    segment_starts = np.zeros(pair_count, dtype=np.float64)
    segment_ends = np.zeros(pair_count, dtype=np.float64)
    if positive:
        jitter = rng.uniform(-0.10, 0.10, size=pair_count) * max(length, 1.0)
        segment_starts = np.clip(center - length / 2 + jitter, 0.0, length_cm - length)
        segment_ends = segment_starts + length
        truth_segments = [
            {
                "pair_id": f"PAIR-{index + 1:03d}",
                "sample_1": sample_ids[2 * index],
                "sample_2": sample_ids[2 * index + 1],
                "chromosome": chromosome,
                "start_cm": float(segment_starts[index]),
                "end_cm": float(segment_ends[index]),
                "length_cm": length,
                "copies": 2 if world_class == "COMPOUND" else 1,
            }
            for index in range(pair_count)
        ]
    previous: np.ndarray | None = None
    weights = _population_weights(sample_count)
    selected = 0
    with (
        gzip.open(source, "wt", encoding="utf-8") as vcf_handle,
        map_path.open("w", encoding="utf-8") as map_handle,
    ):
        vcf_handle.write(
            _vcf_header(
                source="effet_fondateur_R13_scientific_ibd_world",
                chromosome=chromosome,
                length_bp=length_bp,
                sample_ids=sample_ids,
            )
        )
        for dense_rank in range(1, dense_count + 1):
            frequencies = _allele_frequencies(dense_rank)
            previous = _update_haplotypes(rng, previous, weights @ frequencies, refresh)
            position_cm = dense_rank * length_cm / (dense_count + 1)
            if positive:
                active = (segment_starts <= position_cm) & (position_cm <= segment_ends)
                pair_indices = np.flatnonzero(active)
                previous[2 * pair_indices + 1, 0] = previous[2 * pair_indices, 0]
                if world_class == "COMPOUND":
                    previous[2 * pair_indices + 1, 1] = previous[2 * pair_indices, 1]
            if dense_rank % step:
                continue
            selected += 1
            position = max(1, int(dense_rank * length_bp / (dense_count + 1)))
            variant_id = f"R13-IBD-CHR{chromosome:02d}-{dense_rank:07d}"
            allele_count = int(previous.sum())
            vcf_handle.write(
                f"{chromosome}\t{position}\t{variant_id}\tA\tG\t.\tPASS\t"
                f"AC={allele_count};AN={2 * sample_count}\tGT\t"
                + _genotype_text(previous, "|")
                + "\n"
            )
            map_handle.write(
                f"{chromosome}\t{variant_id}\t{position_cm:.8f}\t{position}\n"
            )
    bcftools_name = shutil.which("bcftools")
    if bcftools_name is None:
        _fail("materializer_bcftools_missing")
    vcf = staging / "phased.vcf.gz"
    _compress_and_index_vcf(source, vcf, Path(bcftools_name).resolve(strict=True))
    source.unlink()
    private_manifest = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "truth_id": "T-IBD",
        "world_id": world_id,
        "world_class": world_class,
        "density": density,
        "segments": truth_segments,
    }
    _write_json(private_root / "manifest.json", private_manifest)
    artifacts = [vcf, Path(f"{vcf}.tbi"), map_path]
    manifest = {
        "schema_version": "1.0.0",
        "materializer_id": configuration["materializer_id"],
        "release_id": RELEASE_ID,
        "partition_id": "VALIDATION_BLIND",
        "family": "IBD",
        "world_id": world_id,
        "world_class": world_class,
        "density": density,
        "chromosome": chromosome,
        "variant_count": selected,
        "sample_count": sample_count,
        "pair_count": pair_count,
        "synthetic_only": True,
        "truth_available_to_worker": False,
        "source_truth_binding_sha256": sha256(seed_text.encode()).hexdigest(),
        "artifacts": [
            _staged_artifact(path, root / path.name, run_root) for path in artifacts
        ],
    }
    _write_json(staging / "manifest.json", manifest)
    staging.rename(root)
    return load_json(manifest_path)


def materialize_call(
    repository_root: Path,
    release_root: Path,
    run_root: Path,
    truth_root: Path,
    call_contract: Mapping[str, Any],
) -> dict[str, Any]:
    family = _validate_call_contract(call_contract)
    if family == "shapeit":
        return materialize_shapeit_shard(
            repository_root, release_root, run_root, truth_root, call_contract
        )
    return materialize_ibd_world(
        repository_root, release_root, run_root, truth_root, call_contract
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("audit", "materialize"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--truth-root", type=Path)
    parser.add_argument("--call-id")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.action == "audit":
        document = materializer_contract_audit(arguments.repository_root)
    else:
        if not all(
            (
                arguments.release_root,
                arguments.run_root,
                arguments.truth_root,
                arguments.call_id,
            )
        ):
            _fail("materializer_cli_arguments_missing")
        schedule = build_external_call_schedule()
        matches = [
            call
            for family in ("shapeit", "ibd")
            for call in schedule[family]
            if call["call_id"] == arguments.call_id
        ]
        if len(matches) != 1:
            _fail("materializer_call_unknown", str(arguments.call_id))
        document = materialize_call(
            arguments.repository_root,
            arguments.release_root,
            arguments.run_root,
            arguments.truth_root,
            matches[0],
        )
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
