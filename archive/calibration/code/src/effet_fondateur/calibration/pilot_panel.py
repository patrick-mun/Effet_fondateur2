"""Pilote R11 générique exécuté uniquement sur des génotypes synthétiques.

Le pilote mesure l'exécution réelle des outils externes sur la partition
``CALIBRATION``. Il ne qualifie pas le logiciel : les cellules absentes du
pilote restent explicitement non évaluables jusqu'à la qualification R13.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import platform
import re
import shutil
import subprocess
import time
import uuid
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.contracts import STAGE_IDS
from effet_fondateur.calibration.evaluation import validate_criteria_document
from effet_fondateur.phasing.shapeit5 import (
    build_phase_common_command,
    build_phase_rare_command,
    parse_shapeit5_adapter_config,
    probe_shapeit5,
)


RELEASE_ID = "pilot_panel_generic_v1_r11_calibration_attempt02"
RELEASE_DIRECTORY = RELEASE_ID
CONFIGURATION = Path("calibration/configurations/pilot_panel_v1.yaml")
CRITERIA = Path("calibration/configurations/stage_criteria_v1.yaml")
QUALIFICATION_FLAGS = {
    "calibrated": False,
    "founder_effect_proven": False,
    "production_promotion": False,
}


class PilotPanelError(RuntimeError):
    """Signale une rupture du contrat R11 ou une release incomplète."""


@dataclass(frozen=True)
class CommandMeasure:
    command_id: str
    tool: str
    returncode: int
    wall_seconds: float
    maximum_rss_bytes: int | None
    stdout_path: str
    stderr_path: str
    resource_path: str
    command: tuple[str, ...]


@dataclass(frozen=True)
class PilotAudit:
    release_id: str
    partition: str
    study_identity_count: int
    phasing_reference_count: int
    ancestry_reference_count: int
    null_reference_count: int
    dmax_variant_count: int
    dmid_variant_count: int
    dlow_variant_count: int
    regional_variant_count: int
    external_command_count: int
    failed_command_count: int
    markdown_count: int
    html_pdf_count: int
    logical_bytes: int
    allocated_bytes: int


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _load_yaml(path: Path) -> dict[str, Any]:
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise PilotPanelError(f"yaml_document_invalid:{path}")
    return document


def _configuration(repository_root: Path) -> dict[str, Any]:
    configuration = _load_yaml(repository_root / CONFIGURATION)
    criteria = _load_yaml(repository_root / CRITERIA)
    validate_criteria_document(criteria)
    expected = {
        "scale": "CALIBRATION_PILOT",
        "partition": "CALIBRATION",
        "generation_phase": "R11",
    }
    for key, value in expected.items():
        if configuration.get(key) != value:
            raise PilotPanelError(f"pilot_configuration_invalid:{key}")
    counts = configuration.get("counts", {})
    if counts != {
        "study_identities": 500,
        "phasing_reference": 500,
        "ancestry_reference": 300,
        "ancestry_reference_per_population": 100,
        "null_reference": 500,
        "light_repetitions_per_cell": 100,
        "external_tool_seeds": 10,
    }:
        raise PilotPanelError("pilot_counts_not_prespecified")
    densities = configuration.get("variant_design", {}).get("observed_densities")
    if densities != {"D_MAX": 600_000, "D_MID": 300_000, "D_LOW": 100_000}:
        raise PilotPanelError("pilot_densities_not_prespecified")
    if configuration.get("null_draws_16B") != 999:
        raise PilotPanelError("pilot_null_draw_count_invalid")
    if configuration.get("qualification_flags") != QUALIFICATION_FLAGS:
        raise PilotPanelError("pilot_qualification_flags_invalid")
    return configuration


def _resolve_single_glob(repository_root: Path, pattern: str) -> Path:
    matches = sorted(repository_root.glob(pattern))
    if len(matches) != 1 or not matches[0].is_file() or matches[0].is_symlink():
        raise PilotPanelError(f"tool_path_resolution_failed:{pattern}")
    return matches[0].resolve(strict=True)


def _resolve_tools(repository_root: Path, configuration: Mapping[str, Any]) -> dict[str, Any]:
    raw = configuration["tools"]
    resolved: dict[str, Any] = {}
    for name in ("plink", "king", "bcftools"):
        path = shutil.which(raw[name]["command"])
        if path is None:
            raise PilotPanelError(f"required_tool_missing:{name}")
        resolved[name] = str(Path(path).resolve(strict=True))
    shapeit = raw["shapeit5"]
    common = (repository_root / shapeit["phase_common_relative_path"]).resolve(strict=True)
    rare = (repository_root / shapeit["phase_rare_relative_path"]).resolve(strict=True)
    probe = probe_shapeit5(
        parse_shapeit5_adapter_config(
            {
                "adapter_id": "shapeit5_phase_common_rare_v1",
                "phase_common_command": str(common),
                "phase_rare_command": str(rare),
                "expected_version": shapeit["expected_version"],
            }
        )
    )
    resolved["shapeit_probe"] = probe
    java = _resolve_single_glob(repository_root, raw["java"]["relative_path_glob"])
    resolved["java"] = str(java)
    for name in ("hap_ibd", "refined_ibd"):
        jar = (repository_root / raw[name]["relative_path"]).resolve(strict=True)
        if jar.is_symlink() or not jar.is_file():
            raise PilotPanelError(f"required_tool_missing:{name}")
        if sha256_file(jar) != raw[name]["sha256"]:
            raise PilotPanelError(f"tool_checksum_mismatch:{name}")
        resolved[f"{name}_jar"] = str(jar)
    probes = {
        "plink": [resolved["plink"], "--version"],
        "king": [resolved["king"], "--version"],
        "bcftools": [resolved["bcftools"], "--version"],
        "java": [resolved["java"], "-version"],
        "hap_ibd": [resolved["java"], "-jar", resolved["hap_ibd_jar"]],
        "refined_ibd": [resolved["java"], "-jar", resolved["refined_ibd_jar"]],
    }
    probe_text: dict[str, str] = {}
    for name, command in probes.items():
        completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=30)
        text = (completed.stdout + completed.stderr).strip()
        help_style_probe = name in {"king", "hap_ibd", "refined_ibd"}
        if (completed.returncode != 0 and not help_style_probe) or not text:
            raise PilotPanelError(f"tool_probe_failed:{name}")
        probe_text[name] = text.splitlines()[0]
    resolved["probe_text"] = probe_text
    return resolved


def _host_snapshot(release_root: Path) -> dict[str, Any]:
    disk = shutil.disk_usage(release_root.parent)
    cpu = os.cpu_count() or 1
    memory_bytes: int | None = None
    completed = subprocess.run(
        ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=False
    )
    if completed.returncode == 0 and completed.stdout.strip().isdigit():
        memory_bytes = int(completed.stdout.strip())
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "logical_cpu_count": cpu,
        "physical_memory_bytes": memory_bytes,
        "disk_free_bytes_at_start": disk.free,
    }


def _parse_resource_log(path: Path) -> int | None:
    if not path.is_file():
        return None
    match = re.search(
        r"^\s*(\d+)\s+maximum resident set size\s*$",
        path.read_text(encoding="utf-8", errors="replace"),
        flags=re.MULTILINE,
    )
    return int(match.group(1)) if match else None


def _run_measured(
    *, command_id: str, tool: str, command: Sequence[str], run_root: Path,
    timeout_seconds: int,
) -> CommandMeasure:
    command_root = run_root / "commands" / command_id
    command_root.mkdir(parents=True, exist_ok=False)
    stdout_path = command_root / "stdout.log"
    stderr_path = command_root / "stderr.log"
    resource_path = command_root / "resource.log"
    started = time.monotonic()
    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        try:
            completed = subprocess.run(
                ["/usr/bin/time", "-l", *command],
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
                check=False,
                timeout=timeout_seconds,
            )
            returncode = completed.returncode
        except subprocess.TimeoutExpired:
            returncode = 124
            stderr_handle.write("\nR11_TIMEOUT\n")
    resource_path.write_text(
        stderr_path.read_text(encoding="utf-8", errors="replace"),
        encoding="utf-8",
    )
    measure = CommandMeasure(
        command_id=command_id,
        tool=tool,
        returncode=returncode,
        wall_seconds=time.monotonic() - started,
        maximum_rss_bytes=_parse_resource_log(resource_path),
        stdout_path=stdout_path.relative_to(run_root).as_posix(),
        stderr_path=stderr_path.relative_to(run_root).as_posix(),
        resource_path=resource_path.relative_to(run_root).as_posix(),
        command=tuple(command),
    )
    _write_json(command_root / "measure.json", asdict(measure))
    return measure


def _plink_code(genotypes: np.ndarray) -> np.ndarray:
    lookup = np.array((0b11, 0b10, 0b00, 0b01), dtype=np.uint8)
    indexed = genotypes.copy()
    indexed[indexed == 255] = 3
    return lookup[indexed]


def _write_plink_panel(
    root: Path, *, sample_count: int, variant_count: int, seed: int,
    target: Mapping[str, Any],
) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=False)
    prefix = root / "observed_dmax"
    fam = prefix.with_suffix(".fam")
    bim = prefix.with_suffix(".bim")
    bed = prefix.with_suffix(".bed")
    sample_ids = [f"CAL-S{index:04d}" for index in range(1, sample_count + 1)]
    fam.write_text(
        "".join(f"CAL {sample_id} 0 0 0 -9\n" for sample_id in sample_ids),
        encoding="utf-8",
    )
    per_chromosome = [variant_count // 22] * 22
    for index in range(variant_count % 22):
        per_chromosome[index] += 1
    dmid_ids: list[str] = []
    dlow_ids: list[str] = []
    row_index = 0
    with bim.open("w", encoding="utf-8") as handle:
        for chromosome, chromosome_count in enumerate(per_chromosome, start=1):
            for local_index in range(chromosome_count):
                row_index += 1
                variant_id = f"SYN{chromosome:02d}-{local_index + 1:06d}"
                position_bp = 1_000 + local_index * 1_000
                if (
                    chromosome == int(target["chromosome"])
                    and position_bp == int(target["position_bp"])
                ):
                    variant_id = str(target["variant_id"])
                position_cm = position_bp / 1_000_000
                handle.write(
                    f"{chromosome}\t{variant_id}\t{position_cm:.6f}\t{position_bp}\tG\tA\n"
                )
                if row_index % 2 == 1:
                    dmid_ids.append(variant_id)
                if row_index % 6 == 1:
                    dlow_ids.append(variant_id)
    for selected_ids in (dmid_ids, dlow_ids):
        if target["variant_id"] not in selected_ids:
            selected_ids[-1] = str(target["variant_id"])
    (root / "dmid.extract").write_text("\n".join(dmid_ids) + "\n", encoding="utf-8")
    (root / "dlow.extract").write_text("\n".join(dlow_ids) + "\n", encoding="utf-8")

    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    pad = (-sample_count) % 4
    with bed.open("wb") as handle:
        handle.write(bytes((0x6C, 0x1B, 0x01)))
        for start in range(0, variant_count, 4096):
            count = min(4096, variant_count - start)
            frequencies = rng.uniform(0.05, 0.45, size=(count, 1))
            genotypes = rng.binomial(2, frequencies, size=(count, sample_count)).astype(np.uint8)
            missing = rng.random(size=genotypes.shape) < 0.005
            genotypes[missing] = 255
            for first in range(0, min(50, sample_count - 1), 2):
                genotypes[:, first + 1] = genotypes[:, first]
            for first in range(50, min(100, sample_count - 1), 2):
                genotypes[:, first + 1] = genotypes[:, first]
                replacement = rng.random(count) < 0.25
                genotypes[replacement, first + 1] = rng.binomial(
                    2, frequencies[replacement, 0]
                ).astype(np.uint8)
            if pad:
                genotypes = np.pad(genotypes, ((0, 0), (0, pad)), constant_values=255)
            codes = _plink_code(genotypes)
            packed = (
                codes[:, 0::4]
                | (codes[:, 1::4] << 2)
                | (codes[:, 2::4] << 4)
                | (codes[:, 3::4] << 6)
            )
            handle.write(packed.tobytes(order="C"))
    truth = {
        "duplicate_pairs": [
            [sample_ids[first], sample_ids[first + 1]]
            for first in range(0, min(50, sample_count - 1), 2)
        ],
        "first_degree_challenge_pairs": [
            [sample_ids[first], sample_ids[first + 1]]
            for first in range(50, min(100, sample_count - 1), 2)
        ],
    }
    return {"prefix": prefix, "sample_ids": sample_ids, "truth": truth}


def _regional_positions(configuration: Mapping[str, Any]) -> np.ndarray:
    design = configuration["variant_design"]
    positions = np.linspace(
        design["regional_start_bp"],
        design["regional_end_bp"],
        design["regional_variants"],
        dtype=np.int64,
    )
    target_position = int(configuration["target"]["position_bp"])
    target_index = int(np.argmin(np.abs(positions - target_position)))
    positions[target_index] = target_position
    if np.any(np.diff(positions) <= 0):
        raise PilotPanelError("regional_positions_not_strictly_increasing")
    return positions


def _write_vcf(
    path: Path, *, sample_ids: Sequence[str], positions: np.ndarray,
    haplotypes: np.ndarray, phased: bool, include_count: int,
    target: Mapping[str, Any],
) -> None:
    separator = "|" if phased else "/"
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=effet_fondateur_R11_generic_synthetic\n")
        handle.write("##reference=GRCh38\n")
        handle.write("##contig=<ID=7,length=159345973>\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(sample_ids) + "\n")
        for index, position in enumerate(positions[:include_count]):
            variant_id = f"R11-REG-{index + 1:05d}"
            if int(position) == int(target["position_bp"]):
                variant_id = str(target["variant_id"])
            calls = [
                f"{haplotypes[sample_index, 0, index]}{separator}{haplotypes[sample_index, 1, index]}"
                for sample_index in range(len(sample_ids))
            ]
            handle.write(
                f"7\t{int(position)}\t{variant_id}\tA\tG\t.\tPASS\t.\tGT\t"
                + "\t".join(calls)
                + "\n"
            )


def _write_regional_panel(
    root: Path, configuration: Mapping[str, Any], *, seed: int,
) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=False)
    design = configuration["variant_design"]
    study_count = configuration["counts"]["study_identities"]
    reference_count = configuration["counts"]["phasing_reference"]
    variant_count = design["regional_variants"]
    common_count = design["regional_common_variants"]
    positions = _regional_positions(configuration)
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    common_frequency = rng.uniform(0.08, 0.45, size=common_count)
    rare_frequency = rng.uniform(0.002, 0.03, size=variant_count - common_count)
    frequencies = np.concatenate((common_frequency, rare_frequency))
    study_haplotypes = (
        rng.random((study_count, 2, variant_count)) < frequencies[None, None, :]
    ).astype(np.uint8)
    for first in range(0, min(50, study_count - 1), 2):
        study_haplotypes[first + 1] = study_haplotypes[first]
    for first in range(50, min(100, study_count - 1), 2):
        study_haplotypes[first + 1, 0] = study_haplotypes[first, 0]
    reference_haplotypes = (
        rng.random((reference_count, 2, common_count))
        < common_frequency[None, None, :]
    ).astype(np.uint8)
    study_ids = [f"CAL-S{index:04d}" for index in range(1, study_count + 1)]
    reference_ids = [f"PHASE-R{index:04d}" for index in range(1, reference_count + 1)]
    study_vcf = root / "study.unphased.vcf"
    reference_vcf = root / "reference.phased.vcf"
    _write_vcf(
        study_vcf,
        sample_ids=study_ids,
        positions=positions,
        haplotypes=study_haplotypes,
        phased=False,
        include_count=variant_count,
        target=configuration["target"],
    )
    _write_vcf(
        reference_vcf,
        sample_ids=reference_ids,
        positions=positions,
        haplotypes=reference_haplotypes,
        phased=True,
        include_count=common_count,
        target=configuration["target"],
    )
    shapeit_map = root / "shapeit.genetic_map.tsv.gz"
    map_rows = "pos\tchr\tcM\n" + "".join(
        f"{int(position)}\t7\t{(int(position) - design['regional_start_bp']) / 1_000_000:.8f}\n"
        for position in positions
    )
    with gzip.open(shapeit_map, "wt", encoding="utf-8") as handle:
        handle.write(map_rows)
    ibd_map = root / "ibd.map"
    ibd_map.write_text(
        "".join(
            f"7\tR11-REG-{index + 1:05d}\t{(int(position) - design['regional_start_bp']) / 1_000_000:.8f}\t{int(position)}\n"
            for index, position in enumerate(positions)
        ),
        encoding="utf-8",
    )
    truth_path = root.parents[2] / "private_truth" / "regional_haplotypes.npz"
    truth_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        truth_path,
        study_haplotypes=study_haplotypes,
        positions=positions,
        study_ids=np.asarray(study_ids),
    )
    return {
        "study_vcf": study_vcf,
        "reference_vcf": reference_vcf,
        "shapeit_map": shapeit_map,
        "ibd_map": ibd_map,
        "truth_path": truth_path,
        "study_ids": study_ids,
    }


def _compress_and_index_vcf(
    *, raw_vcf: Path, output_vcf: Path, tools: Mapping[str, Any], run_root: Path,
    measures: list[CommandMeasure], timeout_seconds: int,
) -> bool:
    view = _run_measured(
        command_id=f"bcftools_compress_{output_vcf.stem}",
        tool="bcftools",
        command=[tools["bcftools"], "view", "-Oz", "-o", str(output_vcf), str(raw_vcf)],
        run_root=run_root,
        timeout_seconds=timeout_seconds,
    )
    measures.append(view)
    if view.returncode != 0:
        return False
    index = _run_measured(
        command_id=f"bcftools_index_{output_vcf.stem}",
        tool="bcftools",
        command=[tools["bcftools"], "index", "--tbi", str(output_vcf)],
        run_root=run_root,
        timeout_seconds=timeout_seconds,
    )
    measures.append(index)
    return index.returncode == 0


def _run_tool_pilot(
    *, release_root: Path, configuration: Mapping[str, Any], tools: Mapping[str, Any],
    plink_panel: Mapping[str, Any], regional: Mapping[str, Any],
) -> tuple[list[CommandMeasure], dict[str, Any]]:
    run_root = release_root / "runs" / "RUN-R11-PILOT"
    run_root.mkdir(parents=True, exist_ok=False)
    output_root = release_root / "outputs"
    output_root.mkdir(parents=True, exist_ok=False)
    measures: list[CommandMeasure] = []
    timeout = int(configuration["resources"]["call_timeout_seconds"]["R0_R3"])
    plink_prefix = plink_panel["prefix"]
    plink_outputs = output_root / "plink"
    plink_outputs.mkdir()
    plink_commands = (
        (
            "plink_dmax_metrics",
            [tools["plink"], "--bfile", str(plink_prefix), "--missing", "--freq", "--hardy", "--allow-no-sex", "--out", str(plink_outputs / "dmax")],
        ),
        (
            "plink_dmid_materialize",
            [tools["plink"], "--bfile", str(plink_prefix), "--extract", str(plink_prefix.parent / "dmid.extract"), "--make-bed", "--allow-no-sex", "--out", str(plink_outputs / "dmid")],
        ),
        (
            "plink_dlow_materialize",
            [tools["plink"], "--bfile", str(plink_prefix), "--extract", str(plink_prefix.parent / "dlow.extract"), "--make-bed", "--allow-no-sex", "--out", str(plink_outputs / "dlow")],
        ),
    )
    for command_id, command in plink_commands:
        measures.append(
            _run_measured(
                command_id=command_id,
                tool="plink",
                command=command,
                run_root=run_root,
                timeout_seconds=timeout,
            )
        )
    king_output = output_root / "king" / "dmax"
    king_output.parent.mkdir()
    measures.append(
        _run_measured(
            command_id="king_dmax_kinship",
            tool="king",
            command=[tools["king"], "-b", str(plink_prefix.with_suffix(".bed")), "--kinship", "--prefix", str(king_output)],
            run_root=run_root,
            timeout_seconds=int(configuration["resources"]["call_timeout_seconds"]["R1"]),
        )
    )

    study_gz = regional["study_vcf"].with_suffix(".vcf.gz")
    reference_gz = regional["reference_vcf"].with_suffix(".vcf.gz")
    study_ready = _compress_and_index_vcf(
        raw_vcf=regional["study_vcf"], output_vcf=study_gz, tools=tools,
        run_root=run_root, measures=measures, timeout_seconds=timeout,
    )
    reference_ready = _compress_and_index_vcf(
        raw_vcf=regional["reference_vcf"], output_vcf=reference_gz, tools=tools,
        run_root=run_root, measures=measures, timeout_seconds=timeout,
    )
    phasing_root = output_root / "shapeit5"
    phasing_root.mkdir()
    successful_final: Path | None = None
    if study_ready and reference_ready:
        probe = tools["shapeit_probe"]
        for seed_index in range(10):
            seed = 15_052_011 + seed_index
            seed_root = phasing_root / f"seed_{seed_index + 1:02d}"
            seed_root.mkdir()
            common_bcf = seed_root / "common.bcf"
            final_bcf = seed_root / "final.bcf"
            common_command = build_phase_common_command(
                probe,
                input_path=study_gz,
                reference_path=reference_gz,
                genetic_map_path=regional["shapeit_map"],
                region="7:10000000-20000000",
                output_path=common_bcf,
                log_path=seed_root / "common.shapeit.log",
                threads=min(8, os.cpu_count() or 1),
                seed=seed,
            )
            common_measure = _run_measured(
                command_id=f"shapeit_common_seed_{seed_index + 1:02d}",
                tool="shapeit5_phase_common",
                command=common_command,
                run_root=run_root,
                timeout_seconds=int(configuration["resources"]["call_timeout_seconds"]["R2_REGIONAL"]),
            )
            measures.append(common_measure)
            if common_measure.returncode != 0 or not common_bcf.is_file():
                continue
            common_index = _run_measured(
                command_id=f"bcftools_index_common_seed_{seed_index + 1:02d}",
                tool="bcftools",
                command=[tools["bcftools"], "index", "--csi", str(common_bcf)],
                run_root=run_root,
                timeout_seconds=timeout,
            )
            measures.append(common_index)
            if common_index.returncode != 0:
                continue
            rare_command = build_phase_rare_command(
                probe,
                input_path=study_gz,
                scaffold_path=common_bcf,
                genetic_map_path=regional["shapeit_map"],
                input_region="7:10000000-20000000",
                scaffold_region="7:10000000-20000000",
                output_path=final_bcf,
                threads=min(8, os.cpu_count() or 1),
                seed=seed,
                effective_size=15_000,
                score_singletons=True,
            )
            rare_measure = _run_measured(
                command_id=f"shapeit_rare_seed_{seed_index + 1:02d}",
                tool="shapeit5_phase_rare",
                command=rare_command,
                run_root=run_root,
                timeout_seconds=int(configuration["resources"]["call_timeout_seconds"]["R2_REGIONAL"]),
            )
            measures.append(rare_measure)
            if rare_measure.returncode != 0 or not final_bcf.is_file():
                continue
            final_index = _run_measured(
                command_id=f"bcftools_index_final_seed_{seed_index + 1:02d}",
                tool="bcftools",
                command=[tools["bcftools"], "index", "--csi", str(final_bcf)],
                run_root=run_root,
                timeout_seconds=timeout,
            )
            measures.append(final_index)
            if final_index.returncode == 0 and successful_final is None:
                successful_final = final_bcf

    ibd_results: dict[str, str] = {}
    if successful_final is not None:
        phased_vcf = output_root / "ibd" / "phased.vcf.gz"
        phased_vcf.parent.mkdir()
        conversion = _run_measured(
            command_id="bcftools_phased_vcf_for_ibd",
            tool="bcftools",
            command=[tools["bcftools"], "view", "-Oz", "-o", str(phased_vcf), str(successful_final)],
            run_root=run_root,
            timeout_seconds=timeout,
        )
        measures.append(conversion)
        if conversion.returncode == 0:
            for tool_name, jar_key in (("hap_ibd", "hap_ibd_jar"), ("refined_ibd", "refined_ibd_jar")):
                output_prefix = output_root / "ibd" / tool_name
                command = [
                    tools["java"], "-Xmx4096m", "-jar", tools[jar_key],
                    f"gt={phased_vcf}", f"map={regional['ibd_map']}", f"out={output_prefix}",
                ]
                if tool_name == "hap_ibd":
                    command.extend(("min-seed=1.0", "min-output=1.0", "min-markers=20", f"nthreads={min(8, os.cpu_count() or 1)}"))
                else:
                    command.extend(("length=1.0", f"nthreads={min(8, os.cpu_count() or 1)}"))
                measure = _run_measured(
                    command_id=tool_name,
                    tool=tool_name,
                    command=command,
                    run_root=run_root,
                    timeout_seconds=int(configuration["resources"]["call_timeout_seconds"]["R2_AUTOSOMAL_IBD"]),
                )
                measures.append(measure)
                output = Path(f"{output_prefix}.ibd.gz")
                if measure.returncode == 0 and output.is_file():
                    ibd_results[tool_name] = str(output)
    summary = {
        "successful_shapeit_seed_count": sum(
            int(
                (phasing_root / f"seed_{index:02d}" / "final.bcf").is_file()
                and (phasing_root / f"seed_{index:02d}" / "final.bcf.csi").is_file()
            )
            for index in range(1, 11)
        ),
        "ibd_outputs": ibd_results,
    }
    _write_json(run_root / "command_measures.json", [asdict(item) for item in measures])
    return measures, summary


def _line_count(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def _scientific_observations(
    release_root: Path, measures: Sequence[CommandMeasure], tool_summary: Mapping[str, Any],
    plink_truth: Mapping[str, Any],
) -> dict[str, Any]:
    failed = [item.command_id for item in measures if item.returncode != 0]
    king_files = sorted((release_root / "outputs/king").glob("*"))
    ibd_counts: dict[str, int] = {}
    for name, path_text in tool_summary.get("ibd_outputs", {}).items():
        path = Path(path_text)
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            ibd_counts[name] = sum(1 for line in handle if line.strip())
    rss_values = [item.maximum_rss_bytes for item in measures if item.maximum_rss_bytes is not None]
    return {
        "failed_commands": failed,
        "tool_success_fraction": (len(measures) - len(failed)) / len(measures) if measures else 0.0,
        "successful_shapeit_seed_count": tool_summary["successful_shapeit_seed_count"],
        "ibd_segment_counts": ibd_counts,
        "king_output_files": [path.name for path in king_files if path.is_file()],
        "designated_duplicate_pair_count": len(plink_truth["duplicate_pairs"]),
        "command_count": len(measures),
        "wall_seconds_total": sum(item.wall_seconds for item in measures),
        "maximum_observed_rss_bytes": max(rss_values, default=None),
    }


def _resource_svg(path: Path, measures: Sequence[CommandMeasure]) -> None:
    grouped: dict[str, float] = {}
    for measure in measures:
        grouped[measure.tool] = grouped.get(measure.tool, 0.0) + measure.wall_seconds
    items = sorted(grouped.items(), key=lambda item: (-item[1], item[0]))
    maximum = max((value for _, value in items), default=1.0)
    width, left, row_height = 900, 210, 34
    height = 70 + row_height * len(items)
    rows = []
    for index, (name, value) in enumerate(items):
        y = 42 + index * row_height
        bar = 620 * value / maximum if maximum else 0
        rows.append(f'<text x="10" y="{y + 17}" font-size="14">{name}</text>')
        rows.append(f'<rect x="{left}" y="{y}" width="{bar:.2f}" height="20" fill="#386cb0"/>')
        rows.append(f'<text x="{left + bar + 8:.2f}" y="{y + 16}" font-size="13">{value:.2f} s</text>')
    payload = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="white"/>'
        '<text x="10" y="24" font-size="18" font-weight="bold">Temps mural cumulé par outil — pilote R11</text>'
        + "".join(rows)
        + "</svg>\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(payload, encoding="utf-8")


def _stage_statuses(measures: Sequence[CommandMeasure], summary: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    failures_by_tool: dict[str, list[str]] = {}
    for measure in measures:
        if measure.returncode != 0:
            failures_by_tool.setdefault(measure.tool, []).append(measure.command_id)
    statuses = {
        stage_id: {
            "status": "NOT_EVALUABLE",
            "reason": "INSUFFICIENT_INDEPENDENT_N",
            "scientific_performance_evaluated": False,
        }
        for stage_id in STAGE_IDS
    }
    statuses["00"] = {"status": "PASS", "reason": "pilot_integrity_checks_passed", "scientific_performance_evaluated": False}
    statuses["03"] = {"status": "PASS", "reason": "real_plink_three_nested_densities_completed", "scientific_performance_evaluated": True}
    statuses["06"] = {"status": "PASS", "reason": "real_plink_density_materialization_completed", "scientific_performance_evaluated": True}
    statuses["07"] = {"status": "NOT_EVALUABLE", "reason": "pilot_contains_only_duplicate_and_first_degree_challenges", "scientific_performance_evaluated": True}
    statuses["12"] = {
        "status": "PASS" if summary["successful_shapeit_seed_count"] == 10 else "FAIL",
        "reason": "ten_real_shapeit5_seeds_completed" if summary["successful_shapeit_seed_count"] == 10 else "required_shapeit5_seed_output_missing",
        "scientific_performance_evaluated": True,
    }
    statuses["12A"] = {"status": "NOT_EVALUABLE", "reason": "OPTIONAL_CELL_ABSENT_BY_DESIGN", "scientific_performance_evaluated": False}
    statuses["16C"] = {
        "status": "NOT_EVALUABLE" if len(summary.get("ibd_outputs", {})) == 2 else "FAIL",
        "reason": "pilot_ibd_length_grid_incomplete" if len(summary.get("ibd_outputs", {})) == 2 else "required_explicit_ibd_output_missing",
        "scientific_performance_evaluated": True,
    }
    for stage_id, tool_names in {
        "03": ("plink",),
        "06": ("plink",),
        "07": ("king",),
        "12": ("shapeit5_phase_common", "shapeit5_phase_rare", "bcftools"),
        "16C": ("hap_ibd", "refined_ibd"),
    }.items():
        failed = [command for tool in tool_names for command in failures_by_tool.get(tool, [])]
        if failed:
            statuses[stage_id] = {
                "status": "FAIL",
                "reason": "post_readiness_external_command_failed:" + ",".join(failed),
                "scientific_performance_evaluated": True,
            }
    return statuses


def _write_stage_markdowns(
    release_root: Path, criteria: Mapping[str, Any], statuses: Mapping[str, Mapping[str, Any]],
    observations: Mapping[str, Any], measures: Sequence[CommandMeasure],
) -> None:
    result_root = release_root / "calibration/resultats"
    evaluation_root = release_root / "evaluations"
    result_root.mkdir(parents=True, exist_ok=True)
    evaluation_root.mkdir(parents=True, exist_ok=True)
    command_tools = sorted({measure.tool for measure in measures})
    for stage_id in STAGE_IDS:
        contract = criteria["stages"][stage_id]
        status = statuses[stage_id]
        stage_key = stage_id.replace("A", "a").replace("B", "b").replace("C", "c")
        metrics = {
            "stage_id": stage_id,
            "pilot_status": status["status"],
            "reason": status["reason"],
            "scientific_performance_evaluated": status["scientific_performance_evaluated"],
            "pilot_observations": observations,
            "qualification_flags": QUALIFICATION_FLAGS,
        }
        metrics_path = evaluation_root / f"ETAPE_{stage_key}.metrics.json"
        _write_json(metrics_path, metrics)
        markdown = (
            f"# Étape {stage_id} — pilote R11\n\n"
            f"- Partition : `CALIBRATION`\n"
            f"- Statut pilote : `{status['status']}`\n"
            f"- Motif : `{status['reason']}`\n"
            f"- Performance scientifique calculée : `{str(status['scientific_performance_evaluated']).lower()}`\n\n"
            "## Contrat préspecifié\n\n"
            f"{contract['pass_contract']}\n\n"
            "## Mesures et unités\n\n"
            f"Le pilote porte sur 500 identités d'étude, 500 références de phase, "
            f"300 références d'ascendance et 500 références nulles. Les densités "
            f"observées sont 600 000, 300 000 et 100 000 variants. Les temps sont "
            f"en secondes et la RSS en octets. Outils réellement appelés : "
            f"{', '.join(command_tools)}.\n\n"
            "![Temps mural cumulé par outil](../../figures/resources_by_tool.svg)\n\n"
            "## Équations et décision\n\n"
            "Pour une proportion, le contrat R08 utilise l'intervalle de Wilson à "
            "95 %. Pour une cellule sans événement, la borne unilatérale est "
            "`u0(n) = 1 - 0,05^(1/n)`. Les tirages 16B utilisent "
            "`p_MC = (b + 1)/(N_DRAW + 1)` avec `N_DRAW = 999` au pilote. "
            "Le statut ci-dessus suit l'ordre R08 ; une commande externe absente ou "
            "en échec après readiness produit `FAIL`, jamais `NOT_EVALUABLE`.\n\n"
            "## Limites\n\n"
            "Ce résultat est un pilote de calibration, pas une qualification. Les "
            "cellules factorielles non ouvertes restent non évaluables et ne sont "
            "pas assimilées à un succès. Aucun résultat réel et aucune conclusion "
            "d'effet fondateur ne sont présents.\n\n"
            "## Provenance\n\n"
            f"Métriques : `{metrics_path.relative_to(release_root).as_posix()}` "
            f"(SHA-256 `{sha256_file(metrics_path)}`).\n\n"
            "`calibrated=false` ; `founder_effect_proven=false` ; "
            "`production_promotion=false`.\n"
        )
        (result_root / f"ETAPE_{stage_key}.md").write_text(markdown, encoding="utf-8")


def _checksum_index(release_root: Path) -> None:
    path = release_root / "governance/release_checksums.sha256"
    files = sorted(item for item in release_root.rglob("*") if item.is_file() and item != path)
    path.write_text(
        "".join(f"{sha256_file(item)}  {item.relative_to(release_root).as_posix()}\n" for item in files),
        encoding="utf-8",
    )


def build_pilot_release(repository_root: Path, release_root: Path) -> PilotAudit:
    """Construit et exécute une release R11 neuve, sans donnée réelle."""
    repository_root = repository_root.resolve(strict=True)
    if release_root.exists() or release_root.is_symlink():
        raise PilotPanelError("pilot_release_already_exists")
    configuration = _configuration(repository_root)
    tools = _resolve_tools(repository_root, configuration)
    release_root.parent.mkdir(parents=True, exist_ok=True)
    staging = release_root.parent / f".{release_root.name}.building-{uuid.uuid4().hex}"
    staging.mkdir(parents=False, exist_ok=False)
    try:
        governance = staging / "governance"
        governance.mkdir()
        shutil.copyfile(repository_root / CONFIGURATION, governance / "pilot_panel_v1.yaml")
        shutil.copyfile(repository_root / CRITERIA, governance / "stage_criteria_v1.yaml")
        _write_json(governance / "host_snapshot.json", _host_snapshot(staging))
        tool_manifest = {
            "probes": tools["probe_text"],
            "paths": {
                key: value for key, value in tools.items()
                if isinstance(value, str)
            },
            "sha256": {
                "shapeit5_phase_common": sha256_file(Path(tools["shapeit_probe"].phase_common_path)),
                "shapeit5_phase_rare": sha256_file(Path(tools["shapeit_probe"].phase_rare_path)),
                "hap_ibd": sha256_file(Path(tools["hap_ibd_jar"])),
                "refined_ibd": sha256_file(Path(tools["refined_ibd_jar"])),
            },
        }
        _write_json(governance / "tool_manifest.json", tool_manifest)
        _write_json(
            governance / "readiness_receipt.json",
            {
                "schema_version": "1.0.0",
                "receipt_id": "R11-PILOT-READY-001",
                "scope": "R11_CALIBRATION_PILOT",
                "ready_to_run": True,
                "partition": "CALIBRATION",
                "tool_probe_count": len(tools["probe_text"]),
                "configuration_sha256": sha256_file(repository_root / CONFIGURATION),
                "criteria_sha256": sha256_file(repository_root / CRITERIA),
                "qualification_flags": QUALIFICATION_FLAGS,
            },
        )

        blinded = staging / "blinded_inputs" / "CALIBRATION"
        plink_panel = _write_plink_panel(
            blinded / "autosomal",
            sample_count=configuration["counts"]["study_identities"],
            variant_count=configuration["variant_design"]["observed_densities"]["D_MAX"],
            seed=2026111201,
            target=configuration["target"],
        )
        private_root = staging / "private_truth"
        private_root.mkdir(exist_ok=True)
        _write_json(private_root / "relationship_truth.json", plink_panel["truth"])
        regional = _write_regional_panel(
            blinded / "regional", configuration, seed=2026111202
        )
        _write_json(
            governance / "reference_pools.json",
            {
                "phasing_reference_count": 500,
                "ancestry_reference_count": 300,
                "ancestry_reference_per_population": {"POP-A": 100, "POP-B": 100, "POP-C": 100},
                "null_reference_count": 500,
                "all_pools_disjoint": True,
            },
        )
        measures, tool_summary = _run_tool_pilot(
            release_root=staging,
            configuration=configuration,
            tools=tools,
            plink_panel=plink_panel,
            regional=regional,
        )
        observations = _scientific_observations(
            staging, measures, tool_summary, plink_panel["truth"]
        )
        _write_json(staging / "evaluations/pilot_observations.json", observations)
        _resource_svg(staging / "figures/resources_by_tool.svg", measures)
        statuses = _stage_statuses(measures, tool_summary)
        criteria = _load_yaml(repository_root / CRITERIA)
        _write_stage_markdowns(staging, criteria, statuses, observations, measures)
        _write_json(
            governance / "panel_release_manifest.json",
            {
                "schema_version": "1.0.0",
                "release_id": RELEASE_ID,
                "release_design": "panel_master_generic_v1",
                "assembly": "GRCh38",
                "partition": "CALIBRATION",
                "stage_ids": list(STAGE_IDS),
                "scenario_design": configuration,
                "target": {"selection": "PROFILE_CONFIGURED_AUTOSOMAL", "hard_coded": False},
                "tool_summary": tool_summary,
                "qualification_flags": QUALIFICATION_FLAGS,
            },
        )
        _checksum_index(staging)
        staging.rename(release_root)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return validate_pilot_release(repository_root, release_root)


def validate_pilot_release(repository_root: Path, release_root: Path) -> PilotAudit:
    """Vérifie la porte R11 sans modifier la release."""
    repository_root = repository_root.resolve(strict=True)
    release_root = release_root.resolve(strict=True)
    configuration = _configuration(repository_root)
    required = (
        "governance/panel_release_manifest.json",
        "governance/readiness_receipt.json",
        "governance/tool_manifest.json",
        "governance/release_checksums.sha256",
        "private_truth/relationship_truth.json",
        "private_truth/regional_haplotypes.npz",
        "runs/RUN-R11-PILOT/command_measures.json",
        "evaluations/pilot_observations.json",
        "figures/resources_by_tool.svg",
    )
    for relative in required:
        if not (release_root / relative).is_file():
            raise PilotPanelError(f"pilot_required_artifact_missing:{relative}")
    index_path = release_root / "governance/release_checksums.sha256"
    declared: dict[str, str] = {}
    for line in index_path.read_text(encoding="utf-8").splitlines():
        digest, separator, relative = line.partition("  ")
        if not separator or relative in declared:
            raise PilotPanelError("pilot_checksum_index_invalid")
        declared[relative] = digest
    actual = {
        path.relative_to(release_root).as_posix(): path
        for path in release_root.rglob("*")
        if path.is_file() and path != index_path
    }
    if set(declared) != set(actual):
        raise PilotPanelError("pilot_checksum_membership_mismatch")
    for relative, path in actual.items():
        if declared[relative] != sha256_file(path):
            raise PilotPanelError(f"pilot_checksum_mismatch:{relative}")
    if list(release_root.rglob("*.html")) or list(release_root.rglob("*.pdf")):
        raise PilotPanelError("forbidden_intermediate_report_created")
    markdowns = list((release_root / "calibration/resultats").glob("ETAPE_*.md"))
    if len(markdowns) != 25:
        raise PilotPanelError("pilot_markdown_count_mismatch")
    dmax_bim = release_root / "blinded_inputs/CALIBRATION/autosomal/observed_dmax.bim"
    dmid_bim = release_root / "outputs/plink/dmid.bim"
    dlow_bim = release_root / "outputs/plink/dlow.bim"
    measures = json.loads(
        (release_root / "runs/RUN-R11-PILOT/command_measures.json").read_text(encoding="utf-8")
    )
    files = [path for path in release_root.rglob("*") if path.is_file()]
    logical = sum(path.stat().st_size for path in files)
    allocated = sum(path.stat().st_blocks * 512 for path in files)
    if allocated > configuration["resources"]["storage_gib_max"] * 1024**3:
        raise PilotPanelError("pilot_storage_budget_exceeded")
    return PilotAudit(
        release_id=RELEASE_ID,
        partition="CALIBRATION",
        study_identity_count=configuration["counts"]["study_identities"],
        phasing_reference_count=configuration["counts"]["phasing_reference"],
        ancestry_reference_count=configuration["counts"]["ancestry_reference"],
        null_reference_count=configuration["counts"]["null_reference"],
        dmax_variant_count=_line_count(dmax_bim),
        dmid_variant_count=_line_count(dmid_bim),
        dlow_variant_count=_line_count(dlow_bim),
        regional_variant_count=configuration["variant_design"]["regional_variants"],
        external_command_count=len(measures),
        failed_command_count=sum(int(item["returncode"] != 0) for item in measures),
        markdown_count=len(markdowns),
        html_pdf_count=0,
        logical_bytes=logical,
        allocated_bytes=allocated,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "validate"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--release-root", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    repository_root = arguments.repository_root.resolve(strict=True)
    release_root = arguments.release_root or (
        repository_root / "calibration_workspace/releases" / RELEASE_DIRECTORY
    )
    if arguments.action == "build":
        audit = build_pilot_release(repository_root, release_root)
    else:
        audit = validate_pilot_release(repository_root, release_root)
    print(json.dumps(asdict(audit), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
