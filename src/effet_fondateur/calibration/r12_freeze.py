"""Correction ciblée et gel R12 du candidat de calibration synthétique.

La release R11 reste immuable. R12 dérive uniquement ses entrées régionales,
ajoute les tags ``AC/AN`` avec bcftools, reprend les étapes 12, 12A et 16C,
puis scelle les méthodes, critères, outils et compatibilités avant R13.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import uuid
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.pilot_panel import (
    CONFIGURATION,
    CRITERIA,
    QUALIFICATION_FLAGS,
    CommandMeasure,
    _checksum_index,
    _configuration,
    _load_yaml,
    _resolve_tools,
    _run_measured,
    _write_json,
    validate_pilot_release,
)
from effet_fondateur.phasing.shapeit5 import (
    build_phase_common_command,
    build_phase_rare_command,
)


SOURCE_RELEASE_DIRECTORY = "pilot_panel_generic_v1_r11_calibration_attempt02"
RELEASE_ID = "panel_master_generic_v1_r12_frozen_candidate"
RELEASE_DIRECTORY = RELEASE_ID
COMPATIBILITY = Path("calibration/configurations/compatibility_matrix_v1.yaml")
AFFECTED_STAGES = ("12", "12A", "16C")


@dataclass(frozen=True)
class FrozenCandidateAudit:
    release_id: str
    source_release_id: str
    corrected_vcf_count: int
    shapeit_seed_count: int
    successful_shapeit_seed_count: int
    explicit_ibd_tool_count: int
    failed_command_count: int
    affected_markdown_count: int
    frozen: bool
    validation_blind_opened: bool
    html_pdf_count: int
    logical_bytes: int
    allocated_bytes: int


class R12FreezeError(RuntimeError):
    """Signale une correction incomplète ou un candidat impossible à geler."""


def _validate_compatibility(
    compatibility: Mapping[str, Any], tools: Mapping[str, Any]
) -> None:
    """Vérifie que les outils résolus appartiennent à la matrice gelée."""
    if (
        compatibility.get("freeze_phase") != "R12"
        or compatibility.get("qualification_flags") != QUALIFICATION_FLAGS
    ):
        raise R12FreezeError("r12_compatibility_matrix_invalid")
    expected_fragments = {
        "plink": "PLINK v1.9.0-b.7.11",
        "king": "KING 2.3.2",
        "bcftools": "bcftools 1.21",
        "java": 'openjdk version "17.',
        "hap_ibd": "version 1.0, 15Jun23.92f",
        "refined_ibd": "refined-ibd.17Jan20.102.jar",
    }
    for name, fragment in expected_fragments.items():
        if fragment not in tools["probe_text"][name]:
            raise R12FreezeError(f"r12_tool_version_incompatible:{name}")
    if (
        tools["shapeit_probe"].phase_common_version != "5.1.1"
        or tools["shapeit_probe"].phase_rare_version != "5.1.1"
    ):
        raise R12FreezeError("r12_tool_version_incompatible:shapeit5")
    expected_hashes = {
        "hap_ibd": tools["hap_ibd_jar"],
        "refined_ibd": tools["refined_ibd_jar"],
    }
    for name, path in expected_hashes.items():
        if sha256_file(Path(path)) != compatibility["tools"][name]["sha256"]:
            raise R12FreezeError(f"r12_tool_checksum_incompatible:{name}")


def _copy_checked(source: Path, destination: Path) -> None:
    if source.is_symlink() or not source.is_file():
        raise R12FreezeError(f"r12_source_artifact_invalid:{source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    if sha256_file(source) != sha256_file(destination):
        raise R12FreezeError(f"r12_copy_integrity_failed:{source.name}")


def _source_binding(repository_root: Path, source_root: Path) -> dict[str, Any]:
    validate_pilot_release(repository_root, source_root)
    source_manifest = json.loads(
        (source_root / "governance/panel_release_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    if source_manifest["release_id"] != SOURCE_RELEASE_DIRECTORY:
        raise R12FreezeError("r12_source_release_id_mismatch")
    if (source_root / "blinded_inputs/VALIDATION_BLIND").exists():
        raise R12FreezeError("validation_blind_already_materialized")
    return {
        "source_release_id": source_manifest["release_id"],
        "source_release_manifest_sha256": sha256_file(
            source_root / "governance/panel_release_manifest.json"
        ),
        "source_release_checksums_sha256": sha256_file(
            source_root / "governance/release_checksums.sha256"
        ),
        "source_configuration_sha256": sha256_file(
            source_root / "governance/pilot_panel_v1.yaml"
        ),
        "source_criteria_sha256": sha256_file(
            source_root / "governance/stage_criteria_v1.yaml"
        ),
        "validation_blind_opened": False,
    }


def _prepare_inputs(
    *, source_root: Path, staging: Path, tools: Mapping[str, Any],
    measures: list[CommandMeasure], timeout_seconds: int,
) -> dict[str, Path]:
    source_regional = source_root / "blinded_inputs/CALIBRATION/regional"
    destination = staging / "blinded_inputs/CALIBRATION/regional"
    destination.mkdir(parents=True, exist_ok=False)
    for name in ("shapeit.genetic_map.tsv.gz", "ibd.map"):
        _copy_checked(source_regional / name, destination / name)
    truth_destination = staging / "private_truth/regional_haplotypes.npz"
    _copy_checked(
        source_root / "private_truth/regional_haplotypes.npz", truth_destination
    )
    relationship_destination = staging / "private_truth/relationship_truth.json"
    _copy_checked(
        source_root / "private_truth/relationship_truth.json",
        relationship_destination,
    )

    corrected: dict[str, Path] = {}
    source_names = {
        "study": "study.unphased.vcf.gz",
        "reference": "reference.phased.vcf.gz",
    }
    for role, source_name in source_names.items():
        output_path = destination / f"{role}.ac_an.vcf.gz"
        fill = _run_measured(
            command_id=f"bcftools_fill_tags_{role}",
            tool="bcftools",
            command=[
                tools["bcftools"],
                "+fill-tags",
                str(source_regional / source_name),
                "--output-type",
                "z",
                "--output",
                str(output_path),
                "--",
                "--tags",
                "AC,AN",
            ],
            run_root=staging / "runs/RUN-R12-CORRECTION",
            timeout_seconds=timeout_seconds,
        )
        measures.append(fill)
        if fill.returncode != 0 or not output_path.is_file():
            raise R12FreezeError(f"r12_fill_tags_failed:{role}")
        index = _run_measured(
            command_id=f"bcftools_index_{role}",
            tool="bcftools",
            command=[tools["bcftools"], "index", "--tbi", str(output_path)],
            run_root=staging / "runs/RUN-R12-CORRECTION",
            timeout_seconds=timeout_seconds,
        )
        measures.append(index)
        if index.returncode != 0 or not Path(f"{output_path}.tbi").is_file():
            raise R12FreezeError(f"r12_corrected_index_failed:{role}")
        query = subprocess.run(
            [tools["bcftools"], "query", "-f", "%INFO/AC\t%INFO/AN\n", str(output_path)],
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
        if query.returncode != 0 or not query.stdout:
            raise R12FreezeError(f"r12_corrected_tags_unreadable:{role}")
        for line in query.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 2 or any(value in {"", "."} for value in fields):
                raise R12FreezeError(f"r12_corrected_tags_missing:{role}")
        corrected[role] = output_path
    corrected["shapeit_map"] = destination / "shapeit.genetic_map.tsv.gz"
    corrected["ibd_map"] = destination / "ibd.map"
    corrected["truth"] = truth_destination
    corrected["relationship_truth"] = relationship_destination
    return corrected


def _run_shapeit(
    *, staging: Path, inputs: Mapping[str, Path], tools: Mapping[str, Any],
    measures: list[CommandMeasure], configuration: Mapping[str, Any],
) -> list[Path]:
    output_root = staging / "outputs/shapeit5"
    output_root.mkdir(parents=True, exist_ok=False)
    final_outputs: list[Path] = []
    timeout = int(configuration["resources"]["call_timeout_seconds"]["R2_REGIONAL"])
    for seed_index in range(configuration["counts"]["external_tool_seeds"]):
        seed = 15_052_011 + seed_index
        seed_root = output_root / f"seed_{seed_index + 1:02d}"
        seed_root.mkdir()
        common_bcf = seed_root / "common.bcf"
        final_bcf = seed_root / "final.bcf"
        common = _run_measured(
            command_id=f"shapeit_common_seed_{seed_index + 1:02d}",
            tool="shapeit5_phase_common",
            command=build_phase_common_command(
                tools["shapeit_probe"],
                input_path=inputs["study"],
                reference_path=inputs["reference"],
                genetic_map_path=inputs["shapeit_map"],
                region="7:10000000-20000000",
                output_path=common_bcf,
                log_path=seed_root / "common.shapeit.log",
                threads=1,
                seed=seed,
            ),
            run_root=staging / "runs/RUN-R12-CORRECTION",
            timeout_seconds=timeout,
        )
        measures.append(common)
        if (
            common.returncode != 0
            or not common_bcf.is_file()
            or not Path(f"{common_bcf}.csi").is_file()
        ):
            continue
        rare = _run_measured(
            command_id=f"shapeit_rare_seed_{seed_index + 1:02d}",
            tool="shapeit5_phase_rare",
            command=build_phase_rare_command(
                tools["shapeit_probe"],
                input_path=inputs["study"],
                scaffold_path=common_bcf,
                genetic_map_path=inputs["shapeit_map"],
                input_region="7:10000000-20000000",
                scaffold_region="7:10000000-20000000",
                output_path=final_bcf,
                threads=1,
                seed=seed,
                effective_size=15_000,
                score_singletons=True,
            ),
            run_root=staging / "runs/RUN-R12-CORRECTION",
            timeout_seconds=timeout,
        )
        measures.append(rare)
        if (
            rare.returncode == 0
            and final_bcf.is_file()
            and Path(f"{final_bcf}.csi").is_file()
        ):
            final_outputs.append(final_bcf)
    return final_outputs


def _run_ibd(
    *, staging: Path, final_bcf: Path, inputs: Mapping[str, Path],
    tools: Mapping[str, Any], measures: list[CommandMeasure],
    configuration: Mapping[str, Any],
) -> dict[str, Path]:
    output_root = staging / "outputs/explicit_ibd"
    output_root.mkdir(parents=True, exist_ok=False)
    timeout = int(
        configuration["resources"]["call_timeout_seconds"]["R2_AUTOSOMAL_IBD"]
    )
    phased_vcf = output_root / "phased_seed01.vcf.gz"
    conversion = _run_measured(
        command_id="bcftools_phased_vcf_for_ibd",
        tool="bcftools",
        command=[tools["bcftools"], "view", "-Oz", "-o", str(phased_vcf), str(final_bcf)],
        run_root=staging / "runs/RUN-R12-CORRECTION",
        timeout_seconds=timeout,
    )
    measures.append(conversion)
    if conversion.returncode != 0 or not phased_vcf.is_file():
        raise R12FreezeError("r12_ibd_vcf_conversion_failed")
    outputs: dict[str, Path] = {}
    commands = {
        "hap_ibd": [
            tools["java"],
            "-Xmx4096m",
            "-jar",
            tools["hap_ibd_jar"],
            f"gt={phased_vcf}",
            f"map={inputs['ibd_map']}",
            f"out={output_root / 'hap_ibd'}",
            "min-seed=1.0",
            "min-output=1.0",
            "min-markers=20",
            f"nthreads={min(8, os.cpu_count() or 1)}",
        ],
        "refined_ibd": [
            tools["java"],
            "-Xmx4096m",
            "-jar",
            tools["refined_ibd_jar"],
            f"gt={phased_vcf}",
            f"map={inputs['ibd_map']}",
            f"out={output_root / 'refined_ibd'}",
            "length=1.0",
            f"nthreads={min(8, os.cpu_count() or 1)}",
        ],
    }
    for tool_name, command in commands.items():
        measure = _run_measured(
            command_id=tool_name,
            tool=tool_name,
            command=command,
            run_root=staging / "runs/RUN-R12-CORRECTION",
            timeout_seconds=timeout,
        )
        measures.append(measure)
        output_path = output_root / f"{tool_name}.ibd.gz"
        if measure.returncode != 0 or not output_path.is_file():
            raise R12FreezeError(f"r12_explicit_ibd_failed:{tool_name}")
        outputs[tool_name] = output_path
    return outputs


def _read_phased_calls(
    bcftools: str, bcf_path: Path, sample_ids: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    observed_samples = subprocess.run(
        [bcftools, "query", "--list-samples", str(bcf_path)],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    if observed_samples.returncode != 0 or observed_samples.stdout.splitlines() != list(sample_ids):
        raise R12FreezeError("r12_phased_sample_order_mismatch")
    query = subprocess.run(
        [bcftools, "query", "-f", "%POS[\t%GT]\n", str(bcf_path)],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    if query.returncode != 0:
        raise R12FreezeError("r12_phased_output_query_failed")
    positions: list[int] = []
    calls: list[np.ndarray] = []
    for line in query.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != len(sample_ids) + 1:
            raise R12FreezeError("r12_phased_output_width_mismatch")
        positions.append(int(fields[0]))
        row = np.empty((len(sample_ids), 2), dtype=np.int8)
        for index, genotype in enumerate(fields[1:]):
            if "|" not in genotype:
                raise R12FreezeError("r12_unphased_output_genotype")
            alleles = genotype.split("|")
            if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
                raise R12FreezeError("r12_non_biallelic_output_genotype")
            row[index] = (int(alleles[0]), int(alleles[1]))
        calls.append(row)
    return np.asarray(positions, dtype=np.int64), np.stack(calls, axis=1)


def _switch_error(
    truth: np.ndarray, observed: np.ndarray, variant_slice: slice
) -> tuple[int, int, float | None]:
    switches = 0
    transitions = 0
    subset_truth = truth[:, :, variant_slice]
    subset_observed = observed[:, variant_slice, :]
    for sample_index in range(truth.shape[0]):
        truth_sample = subset_truth[sample_index]
        observed_sample = subset_observed[sample_index]
        informative = np.flatnonzero(
            (truth_sample[0] + truth_sample[1] == 1)
            & (observed_sample[:, 0] + observed_sample[:, 1] == 1)
        )
        if informative.size < 2:
            continue
        orientation = (
            observed_sample[informative, 0] != truth_sample[0, informative]
        ).astype(np.int8)
        switches += int(np.count_nonzero(np.diff(orientation)))
        transitions += int(informative.size - 1)
    return switches, transitions, switches / transitions if transitions else None


def _evaluate_phasing(
    *, inputs: Mapping[str, Path], final_outputs: Sequence[Path],
    bcftools: str, configuration: Mapping[str, Any],
) -> dict[str, Any]:
    with np.load(inputs["truth"], allow_pickle=False) as archive:
        truth = archive["study_haplotypes"].astype(np.int8)
        truth_positions = archive["positions"].astype(np.int64)
        sample_ids = [str(value) for value in archive["study_ids"]]
    common_count = int(configuration["variant_design"]["regional_common_variants"])
    per_seed: list[dict[str, Any]] = []
    for index, final_bcf in enumerate(final_outputs, start=1):
        positions, calls = _read_phased_calls(bcftools, final_bcf, sample_ids)
        if not np.array_equal(positions, truth_positions):
            raise R12FreezeError("r12_phased_variant_order_mismatch")
        truth_genotypes = truth[:, 0, :] + truth[:, 1, :]
        observed_genotypes = calls[:, :, 0] + calls[:, :, 1]
        concordance = float(np.mean(truth_genotypes == observed_genotypes))
        common = _switch_error(truth, calls, slice(0, common_count))
        rare = _switch_error(truth, calls, slice(common_count, truth.shape[2]))
        all_variants = _switch_error(truth, calls, slice(0, truth.shape[2]))
        per_seed.append(
            {
                "seed_index": index,
                "genotype_concordance": concordance,
                "common_switches": common[0],
                "common_transitions": common[1],
                "common_ser": common[2],
                "rare_switches": rare[0],
                "rare_transitions": rare[1],
                "rare_ser": rare[2],
                "all_switches": all_variants[0],
                "all_transitions": all_variants[1],
                "all_ser": all_variants[2],
            }
        )
    return {
        "seed_count": len(per_seed),
        "sample_count": truth.shape[0],
        "variant_count": truth.shape[2],
        "per_seed": per_seed,
        "genotype_concordance_min": min(row["genotype_concordance"] for row in per_seed),
        "common_ser_mean": float(np.mean([row["common_ser"] for row in per_seed])),
        "common_ser_min": min(row["common_ser"] for row in per_seed),
        "common_ser_max": max(row["common_ser"] for row in per_seed),
        "rare_ser_mean": float(np.mean([row["rare_ser"] for row in per_seed])),
        "all_ser_mean": float(np.mean([row["all_ser"] for row in per_seed])),
    }


def _ibd_pairs(path: Path) -> tuple[int, set[tuple[str, str]]]:
    count = 0
    pairs: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise R12FreezeError(f"r12_ibd_output_malformed:{path.name}")
            count += 1
            pairs.add(tuple(sorted((fields[0], fields[2]))))
    return count, pairs


def _evaluate_ibd(
    outputs: Mapping[str, Path], relationship_truth_path: Path
) -> dict[str, Any]:
    truth = json.loads(relationship_truth_path.read_text(encoding="utf-8"))
    duplicate_pairs = {tuple(sorted(pair)) for pair in truth["duplicate_pairs"]}
    results: dict[str, Any] = {}
    observed_sets: dict[str, set[tuple[str, str]]] = {}
    for tool_name, path in outputs.items():
        segment_count, pairs = _ibd_pairs(path)
        observed_sets[tool_name] = pairs
        detected = duplicate_pairs & pairs
        results[tool_name] = {
            "segment_count": segment_count,
            "pair_count": len(pairs),
            "designated_duplicate_detected": len(detected),
            "designated_duplicate_expected": len(duplicate_pairs),
            "designated_duplicate_recall": len(detected) / len(duplicate_pairs),
        }
    intersection = observed_sets["hap_ibd"] & observed_sets["refined_ibd"]
    union = observed_sets["hap_ibd"] | observed_sets["refined_ibd"]
    return {
        "tools": results,
        "pair_jaccard": len(intersection) / len(union) if union else None,
        "ancestry_control_status": "NOT_EVALUABLE",
        "ancestry_control_reason": "OPTIONAL_CELL_ABSENT_BY_DESIGN",
        "length_grid_complete": False,
    }


def _phase_svg(path: Path, phasing: Mapping[str, Any]) -> None:
    values = [float(row["common_ser"]) for row in phasing["per_seed"]]
    width, height = 820, 380
    left, top, chart_width, chart_height = 70, 45, 700, 260
    maximum = max(max(values, default=0.0), 0.01)
    points = []
    for index, value in enumerate(values):
        x = left + index * chart_width / max(1, len(values) - 1)
        y = top + chart_height * (1 - value / maximum)
        points.append((x, y, value))
    polyline = " ".join(f"{x:.2f},{y:.2f}" for x, y, _ in points)
    circles = "".join(
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="#386cb0"/>'
        for x, y, _ in points
    )
    labels = "".join(
        f'<text x="{x:.2f}" y="{top + chart_height + 24}" text-anchor="middle" font-size="12">{index + 1}</text>'
        for index, (x, _, _) in enumerate(points)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
        '<rect width="100%" height="100%" fill="white"/>'
        '<text x="20" y="24" font-size="18" font-weight="bold">SER commun par graine SHAPEIT5 — R12</text>'
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_height}" stroke="black"/>'
        f'<line x1="{left}" y1="{top + chart_height}" x2="{left + chart_width}" y2="{top + chart_height}" stroke="black"/>'
        f'<polyline points="{polyline}" fill="none" stroke="#386cb0" stroke-width="2"/>{circles}{labels}'
        f'<text x="20" y="{top + 12}" font-size="12">{maximum:.4f}</text>'
        f'<text x="20" y="{top + chart_height}" font-size="12">0</text>'
        f'<text x="{left + chart_width / 2}" y="{height - 20}" text-anchor="middle" font-size="13">Graine indexée (N=10)</text>'
        '</svg>\n',
        encoding="utf-8",
    )


def _ibd_svg(path: Path, ibd: Mapping[str, Any]) -> None:
    names = ("hap_ibd", "refined_ibd")
    counts = [ibd["tools"][name]["segment_count"] for name in names]
    maximum = max(counts, default=1)
    bars = []
    for index, (name, count) in enumerate(zip(names, counts, strict=True)):
        x = 130 + index * 250
        height = 220 * count / maximum if maximum else 0
        y = 290 - height
        bars.append(f'<rect x="{x}" y="{y:.2f}" width="120" height="{height:.2f}" fill="#7fc97f"/>')
        bars.append(f'<text x="{x + 60}" y="{y - 8:.2f}" text-anchor="middle" font-size="14">{count}</text>')
        bars.append(f'<text x="{x + 60}" y="315" text-anchor="middle" font-size="13">{name}</text>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="650" height="360">'
        '<rect width="100%" height="100%" fill="white"/>'
        '<text x="20" y="25" font-size="18" font-weight="bold">Segments IBD appelés — reprise R12</text>'
        '<line x1="80" y1="45" x2="80" y2="290" stroke="black"/>'
        '<line x1="80" y1="290" x2="600" y2="290" stroke="black"/>'
        + "".join(bars)
        + '</svg>\n',
        encoding="utf-8",
    )


def _write_affected_reports(
    *, staging: Path, phasing: Mapping[str, Any], ibd: Mapping[str, Any],
    command_failures: Sequence[str],
) -> dict[str, str]:
    result_root = staging / "calibration/resultats"
    evaluation_root = staging / "evaluations"
    result_root.mkdir(parents=True, exist_ok=True)
    evaluation_root.mkdir(parents=True, exist_ok=True)
    if command_failures:
        statuses = {stage: "FAIL" for stage in AFFECTED_STAGES}
    else:
        statuses = {stage: "NOT_EVALUABLE" for stage in AFFECTED_STAGES}
    reasons = {
        "12": "OPTIONAL_CELL_ABSENT_BY_DESIGN",
        "12A": "OPTIONAL_CELL_ABSENT_BY_DESIGN",
        "16C": "OPTIONAL_CELL_ABSENT_BY_DESIGN",
    }
    for stage in AFFECTED_STAGES:
        suffix = stage.lower()
        metrics = {
            "stage_id": stage,
            "status": statuses[stage],
            "reason": "external_command_failed" if command_failures else reasons[stage],
            "command_failures": list(command_failures),
            "phasing": phasing if stage in {"12", "12A"} else None,
            "explicit_ibd": ibd if stage == "16C" else None,
            "qualification_flags": QUALIFICATION_FLAGS,
        }
        metrics_path = evaluation_root / f"ETAPE_{suffix}.r12.metrics.json"
        _write_json(metrics_path, metrics)
        if stage == "12":
            explanation = (
                f"Les 10 graines ont produit 1 200 variants et 500 sujets. La concordance "
                f"génotypique minimale vaut {phasing['genotype_concordance_min']:.6f}; "
                f"le SER commun moyen vaut {phasing['common_ser_mean']:.6f}. La grille "
                "densité × fréquence × pedigree n'est pas complète dans ce pilote."
            )
            figure = "../../figures/shapeit_common_ser.svg"
        elif stage == "12A":
            explanation = (
                "La correction régionale est démontrée, mais la reprise R12 n'ouvre pas "
                "les 22 autosomes : cette cellule reste structurellement absente du pilote."
            )
            figure = "../../figures/shapeit_common_ser.svg"
        else:
            explanation = (
                f"Hap-IBD et Refined IBD ont produit leurs sorties. Le Jaccard des paires "
                f"appelées vaut {ibd['pair_jaccard'] if ibd['pair_jaccard'] is not None else 'NA'}. "
                "La grille 1/1,5/2/3/5 cM et le contrôle 16A ne sont pas complets."
            )
            figure = "../../figures/ibd_segment_counts.svg"
        markdown = (
            f"# Étape {stage} — reprise et gel R12\n\n"
            f"- Statut : `{statuses[stage]}`\n"
            f"- Motif : `{metrics['reason']}`\n"
            "- Partition : `CALIBRATION`\n\n"
            "## Résultat\n\n"
            f"{explanation}\n\n"
            f"![Mesure R12]({figure})\n\n"
            "## Équations et unités\n\n"
            "Le taux d'erreur de commutation est `SER = N_SWITCH / N_TRANSITION`, "
            "calculé sur les transitions hétérozygotes informatives. Le Jaccard IBD "
            "est `|P_HAP ∩ P_REFINED| / |P_HAP ∪ P_REFINED|`, au niveau des paires. "
            "Les longueurs IBD sont en cM et les positions en bp.\n\n"
            "## Limites\n\n"
            "R12 ferme la défaillance technique `AC/AN`; elle ne transforme pas le "
            "pilote en qualification. Les cellules absentes restent `NOT_EVALUABLE` "
            "et aucun seuil R08 n'est modifié.\n\n"
            "## Provenance\n\n"
            f"Métriques : `{metrics_path.relative_to(staging).as_posix()}` "
            f"(SHA-256 `{sha256_file(metrics_path)}`).\n\n"
            "`calibrated=false` ; `founder_effect_proven=false` ; "
            "`production_promotion=false`.\n"
        )
        (result_root / f"ETAPE_{suffix}.md").write_text(markdown, encoding="utf-8")
    return statuses


def _freeze_manifest(
    *, repository_root: Path, staging: Path, source_binding: Mapping[str, Any],
    configuration: Mapping[str, Any], tools: Mapping[str, Any],
    measures: Sequence[CommandMeasure], phasing: Mapping[str, Any],
    ibd: Mapping[str, Any], statuses: Mapping[str, str],
) -> dict[str, Any]:
    criteria_path = repository_root / CRITERIA
    configuration_path = repository_root / CONFIGURATION
    compatibility_path = repository_root / COMPATIBILITY
    production_files = (
        repository_root / "src/effet_fondateur/phasing/shapeit5.py",
        repository_root / "src/effet_fondateur/phasing/execution.py",
        repository_root / "src/effet_fondateur/explicit_ibd/execution.py",
    )
    command_failures = [measure.command_id for measure in measures if measure.returncode != 0]
    return {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "release_design": "panel_master_generic_v1",
        "phase": "R12",
        "state": "FROZEN",
        "frozen_at": "2026-11-12",
        "source_binding": dict(source_binding),
        "correction": {
            "scope": "CALIBRATION_HARNESS_ONLY",
            "cause": "R11_SYNTHETIC_VCF_MISSING_INFO_AC_AN",
            "method": "bcftools_fill_tags_AC_AN",
            "production_software_modified": False,
            "scientific_method_modified": False,
            "criteria_amended": False,
            "affected_stages_rerun": list(AFFECTED_STAGES),
        },
        "frozen_hashes": {
            "pilot_configuration": sha256_file(configuration_path),
            "stage_criteria": sha256_file(criteria_path),
            "compatibility_matrix": sha256_file(compatibility_path),
            "r12_freeze_module": sha256_file(
                repository_root / "src/effet_fondateur/calibration/r12_freeze.py"
            ),
            "production_modules": {
                path.relative_to(repository_root).as_posix(): sha256_file(path)
                for path in production_files
            },
            "shapeit5_phase_common": sha256_file(Path(tools["shapeit_probe"].phase_common_path)),
            "shapeit5_phase_rare": sha256_file(Path(tools["shapeit_probe"].phase_rare_path)),
            "hap_ibd": sha256_file(Path(tools["hap_ibd_jar"])),
            "refined_ibd": sha256_file(Path(tools["refined_ibd_jar"])),
        },
        "dimensions_unchanged": {
            "study_identities": configuration["counts"]["study_identities"],
            "phasing_reference": configuration["counts"]["phasing_reference"],
            "observed_densities": configuration["variant_design"]["observed_densities"],
            "regional_variants": configuration["variant_design"]["regional_variants"],
            "external_tool_seeds": configuration["counts"]["external_tool_seeds"],
        },
        "tool_versions": tools["probe_text"],
        "execution": {
            "command_count": len(measures),
            "failed_commands": command_failures,
            "successful_shapeit_seed_count": phasing["seed_count"],
            "explicit_ibd_tools": sorted(ibd["tools"]),
            "stage_statuses": dict(statuses),
        },
        "validation_blind_opened": False,
        "qualification_flags": QUALIFICATION_FLAGS,
    }


def build_frozen_candidate(
    repository_root: Path, release_root: Path, source_root: Path
) -> FrozenCandidateAudit:
    """Construit un candidat R12 neuf puis le scelle sans ouvrir R13."""
    repository_root = repository_root.resolve(strict=True)
    source_root = source_root.resolve(strict=True)
    if release_root.exists() or release_root.is_symlink():
        raise R12FreezeError("r12_release_already_exists")
    configuration = _configuration(repository_root)
    compatibility = _load_yaml(repository_root / COMPATIBILITY)
    tools = _resolve_tools(repository_root, configuration)
    _validate_compatibility(compatibility, tools)
    source_binding = _source_binding(repository_root, source_root)
    if source_binding["source_configuration_sha256"] != sha256_file(repository_root / CONFIGURATION):
        raise R12FreezeError("r12_configuration_changed_since_r11")
    if source_binding["source_criteria_sha256"] != sha256_file(repository_root / CRITERIA):
        raise R12FreezeError("r12_criteria_changed_since_r11")

    release_root.parent.mkdir(parents=True, exist_ok=True)
    staging = release_root.parent / f".{release_root.name}.building-{uuid.uuid4().hex}"
    staging.mkdir(parents=False, exist_ok=False)
    measures: list[CommandMeasure] = []
    try:
        run_root = staging / "runs/RUN-R12-CORRECTION"
        run_root.mkdir(parents=True)
        governance = staging / "governance"
        governance.mkdir()
        _copy_checked(repository_root / CONFIGURATION, governance / "pilot_panel_v1.yaml")
        _copy_checked(repository_root / CRITERIA, governance / "stage_criteria_v1.yaml")
        _copy_checked(repository_root / COMPATIBILITY, governance / "compatibility_matrix_v1.yaml")
        _write_json(governance / "source_binding.json", source_binding)
        inputs = _prepare_inputs(
            source_root=source_root,
            staging=staging,
            tools=tools,
            measures=measures,
            timeout_seconds=int(configuration["resources"]["call_timeout_seconds"]["R0_R3"]),
        )
        final_outputs = _run_shapeit(
            staging=staging,
            inputs=inputs,
            tools=tools,
            measures=measures,
            configuration=configuration,
        )
        if len(final_outputs) != configuration["counts"]["external_tool_seeds"]:
            raise R12FreezeError(f"r12_shapeit_seed_outputs_incomplete:{len(final_outputs)}")
        ibd_outputs = _run_ibd(
            staging=staging,
            final_bcf=final_outputs[0],
            inputs=inputs,
            tools=tools,
            measures=measures,
            configuration=configuration,
        )
        _write_json(run_root / "command_measures.json", [asdict(item) for item in measures])
        phasing = _evaluate_phasing(
            inputs=inputs,
            final_outputs=final_outputs,
            bcftools=tools["bcftools"],
            configuration=configuration,
        )
        ibd = _evaluate_ibd(ibd_outputs, inputs["relationship_truth"])
        _write_json(staging / "evaluations/phasing_r12.json", phasing)
        _write_json(staging / "evaluations/explicit_ibd_r12.json", ibd)
        _phase_svg(staging / "figures/shapeit_common_ser.svg", phasing)
        _ibd_svg(staging / "figures/ibd_segment_counts.svg", ibd)
        failures = [measure.command_id for measure in measures if measure.returncode != 0]
        statuses = _write_affected_reports(
            staging=staging,
            phasing=phasing,
            ibd=ibd,
            command_failures=failures,
        )
        manifest = _freeze_manifest(
            repository_root=repository_root,
            staging=staging,
            source_binding=source_binding,
            configuration=configuration,
            tools=tools,
            measures=measures,
            phasing=phasing,
            ibd=ibd,
            statuses=statuses,
        )
        _write_json(governance / "frozen_candidate_manifest.json", manifest)
        _checksum_index(staging)
        staging.rename(release_root)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return validate_frozen_candidate(repository_root, release_root)


def validate_frozen_candidate(
    repository_root: Path, release_root: Path
) -> FrozenCandidateAudit:
    """Vérifie le gel R12, les empreintes et l'absence d'ouverture R13."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    required = (
        "governance/frozen_candidate_manifest.json",
        "governance/source_binding.json",
        "governance/release_checksums.sha256",
        "runs/RUN-R12-CORRECTION/command_measures.json",
        "evaluations/phasing_r12.json",
        "evaluations/explicit_ibd_r12.json",
        "figures/shapeit_common_ser.svg",
        "figures/ibd_segment_counts.svg",
    )
    for relative in required:
        if not (release_root / relative).is_file():
            raise R12FreezeError(f"r12_required_artifact_missing:{relative}")
    if (release_root / "blinded_inputs/VALIDATION_BLIND").exists():
        raise R12FreezeError("r12_validation_blind_opened")
    if list(release_root.rglob("*.html")) or list(release_root.rglob("*.pdf")):
        raise R12FreezeError("r12_forbidden_report_created")
    manifest = json.loads(
        (release_root / "governance/frozen_candidate_manifest.json").read_text(encoding="utf-8")
    )
    if (
        manifest.get("state") != "FROZEN"
        or manifest.get("validation_blind_opened") is not False
        or manifest.get("qualification_flags") != QUALIFICATION_FLAGS
        or manifest["correction"]["production_software_modified"] is not False
        or manifest["correction"]["criteria_amended"] is not False
        or manifest["execution"]["failed_commands"]
    ):
        raise R12FreezeError("r12_frozen_manifest_invalid")
    frozen_hashes = manifest["frozen_hashes"]
    if frozen_hashes["pilot_configuration"] != sha256_file(repository_root / CONFIGURATION):
        raise R12FreezeError("r12_frozen_configuration_changed")
    if frozen_hashes["stage_criteria"] != sha256_file(repository_root / CRITERIA):
        raise R12FreezeError("r12_frozen_criteria_changed")
    if frozen_hashes["compatibility_matrix"] != sha256_file(repository_root / COMPATIBILITY):
        raise R12FreezeError("r12_frozen_compatibility_changed")
    if frozen_hashes["r12_freeze_module"] != sha256_file(
        repository_root / "src/effet_fondateur/calibration/r12_freeze.py"
    ):
        raise R12FreezeError("r12_frozen_harness_changed")
    for relative, expected_digest in frozen_hashes["production_modules"].items():
        production_path = repository_root / relative
        if not production_path.is_file() or sha256_file(production_path) != expected_digest:
            raise R12FreezeError(f"r12_frozen_production_module_changed:{relative}")
    index_path = release_root / "governance/release_checksums.sha256"
    declared: dict[str, str] = {}
    for line in index_path.read_text(encoding="utf-8").splitlines():
        digest, separator, relative = line.partition("  ")
        if not separator or relative in declared:
            raise R12FreezeError("r12_checksum_index_invalid")
        declared[relative] = digest
    actual = {
        path.relative_to(release_root).as_posix(): path
        for path in release_root.rglob("*")
        if path.is_file() and path != index_path
    }
    if set(declared) != set(actual):
        raise R12FreezeError("r12_checksum_membership_mismatch")
    for relative, path in actual.items():
        if declared[relative] != sha256_file(path):
            raise R12FreezeError(f"r12_checksum_mismatch:{relative}")
    measures = json.loads(
        (release_root / "runs/RUN-R12-CORRECTION/command_measures.json").read_text(encoding="utf-8")
    )
    phasing = json.loads(
        (release_root / "evaluations/phasing_r12.json").read_text(encoding="utf-8")
    )
    ibd = json.loads(
        (release_root / "evaluations/explicit_ibd_r12.json").read_text(encoding="utf-8")
    )
    markdowns = list((release_root / "calibration/resultats").glob("ETAPE_*.md"))
    if len(markdowns) != 3:
        raise R12FreezeError("r12_affected_markdown_count_mismatch")
    files = [path for path in release_root.rglob("*") if path.is_file()]
    logical = sum(path.stat().st_size for path in files)
    allocated = sum(path.stat().st_blocks * 512 for path in files)
    return FrozenCandidateAudit(
        release_id=RELEASE_ID,
        source_release_id=SOURCE_RELEASE_DIRECTORY,
        corrected_vcf_count=2,
        shapeit_seed_count=10,
        successful_shapeit_seed_count=int(phasing["seed_count"]),
        explicit_ibd_tool_count=len(ibd["tools"]),
        failed_command_count=sum(int(row["returncode"] != 0) for row in measures),
        affected_markdown_count=len(markdowns),
        frozen=True,
        validation_blind_opened=False,
        html_pdf_count=0,
        logical_bytes=logical,
        allocated_bytes=allocated,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "validate"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path)
    parser.add_argument("--source-root", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    repository_root = arguments.repository_root.resolve(strict=True)
    release_root = arguments.release_root or (
        repository_root / "calibration_workspace/releases" / RELEASE_DIRECTORY
    )
    source_root = arguments.source_root or (
        repository_root / "calibration_workspace/releases" / SOURCE_RELEASE_DIRECTORY
    )
    if arguments.action == "build":
        audit = build_frozen_candidate(repository_root, release_root, source_root)
    else:
        audit = validate_frozen_candidate(repository_root, release_root)
    print(json.dumps(asdict(audit), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
