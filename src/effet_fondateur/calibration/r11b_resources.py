"""Sondes R11b synthétiques prescrites par le dimensionnement R07 v2.

Le harnais prépare des entrées ``RESOURCE_ONLY``, mesure trois appels SHAPEIT5
et dix-huit appels IBD, puis produit une projection conservatrice. Il ne lit
aucune donnée réelle, aucune vérité de validation et ne mesure aucune
performance scientifique.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import json
from math import ceil, exp, log
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.stats import t as student_t

from effet_fondateur.calibration.contracts import file_sha256, load_json, load_yaml
from effet_fondateur.calibration.pilot_panel import CommandMeasure, _run_measured
from effet_fondateur.calibration.r07_revision import (
    CONFIGURATION_PATH as R07_CONFIGURATION_PATH,
    RECEIPT_PATH as R07_RECEIPT_PATH,
    validate_tracked_receipt,
)


WORKSPACE_RELATIVE = Path("calibration_workspace/resource_probes/r11b_scientific_resources_v1")
HARNESS_RELATIVE = Path("src/effet_fondateur/calibration/r11b_resources.py")
SOURCE_PROBE_RELATIVE = Path(
    "calibration_workspace/resource_probes/r11_shapeit_median_chr_attempt02"
)
N250_PROBE_RELATIVE = Path(
    "calibration_workspace/resource_probes/r11_shapeit_sample_scale_n250_attempt01"
)
OLD_PROJECTION_RELATIVE = Path(
    "calibration_workspace/resource_probes/r11_shapeit_sample_scaling_v1/projection.json"
)
SHAPEIT_RELATIVE = Path("calibration_workspace/toolchains/r11/bin/SHAPEIT5_phase_common")
JAVA_RELATIVE = Path(
    "calibration_workspace/toolchains/r11/java/"
    "jdk-17.0.21.0.2+1/Contents/Home/bin/java"
)
HAP_IBD_RELATIVE = Path("calibration_workspace/toolchains/r11/jars/hap-ibd.jar")
REFINED_IBD_RELATIVE = Path(
    "calibration_workspace/toolchains/r11/jars/refined-ibd.17Jan20.102.jar"
)
REGION = "7:10000000-20000000"
VARIANT_COUNT = 26_000
QUALIFICATION_FLAGS = {
    "calibrated": False,
    "founder_effect_proven": False,
    "production_promotion": False,
}
PROBE_DIMENSIONS = {
    "dlow": (250, 250, "D_LOW"),
    "bridge": (2_000, 1_500, "D_MAX"),
    "target": (4_000, 3_000, "D_MAX"),
}
IBD_DENSITY_STEP = {"D_MAX": 1, "D_MID": 2, "D_LOW": 6}
IBD_TOOLS = ("HAP_IBD", "REFINED_IBD")
IBD_REPLICATES = 3
ORIGINAL_HARNESS_SHA256 = "0313d7bdf82e1ddd692dfed3c3438f909d240a5f71119e7fe9a7cfb5d5057502"
AMENDMENT_001_HARNESS_SHA256 = (
    "753ce65bcb6c1d9dd3a608910045b00069f4469d828113c5552caf55f9ce9a16"
)
SHAPEIT_ATTEMPT_NUMBER = {"dlow": 1, "bridge": 2, "target": 1}


class R11bResourceError(RuntimeError):
    """Erreur stable du harnais de ressources R11b."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R11bResourceError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _resolve(repository_root: Path, relative: Path) -> Path:
    return (repository_root / relative).resolve(strict=True)


def _tool_paths(repository_root: Path) -> dict[str, Path]:
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        _fail("r11b_bcftools_missing")
    paths = {
        "bcftools": Path(bcftools).resolve(strict=True),
        "shapeit": _resolve(repository_root, SHAPEIT_RELATIVE),
        "java": _resolve(repository_root, JAVA_RELATIVE),
        "hap_ibd": _resolve(repository_root, HAP_IBD_RELATIVE),
        "refined_ibd": _resolve(repository_root, REFINED_IBD_RELATIVE),
    }
    expected_hashes = {
        "hap_ibd": "846ddf08db5590339c58072fa9944fff8466f9c2825f9c02612cb9d38492752c",
        "refined_ibd": "a2c2d7ee6c1dff5c06831a5987ba3a3ce831cb5798a13d413ab0c10a577db919",
    }
    for name, expected in expected_hashes.items():
        if file_sha256(paths[name]) != expected:
            _fail("r11b_tool_hash_mismatch", name)
    return paths


def _run_checked(command: Sequence[str], *, stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
    else:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                command, stdout=handle, stderr=subprocess.PIPE, text=True, check=False
            )
    if completed.returncode != 0:
        _fail("r11b_command_failed", f"{command[0]}:{completed.stderr[-1000:]}")


def _index_vcf(bcftools: Path, path: Path) -> None:
    index = Path(f"{path}.tbi")
    if index.is_file():
        return
    _run_checked([str(bcftools), "index", "--tbi", str(path)])


def _sample_ids(prefix: str, count: int) -> list[str]:
    return [f"{prefix}{index:05d}" for index in range(1, count + 1)]


def _variant_positions() -> np.ndarray:
    full = np.linspace(10_000_000, 20_000_000, 27_000, dtype=np.int64)
    return full[:VARIANT_COUNT]


def _stream_synthetic_vcf(
    *,
    output_path: Path,
    sample_ids: Sequence[str],
    phased: bool,
    seed: int,
    bcftools: Path,
) -> None:
    """Génère un VCF synthétique bgzip sans conserver la matrice en mémoire."""
    if output_path.is_file() and Path(f"{output_path}.tbi").is_file():
        return
    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial = output_path.with_name(output_path.name + ".partial")
    partial.unlink(missing_ok=True)
    separator = "|" if phased else "/"
    command = [str(bcftools), "view", "-Oz", "-o", str(partial), "-"]
    process = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if process.stdin is None:
        _fail("r11b_vcf_stream_unavailable")
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    try:
        handle = process.stdin
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=effet_fondateur_R11b_resource_only\n")
        handle.write("##reference=GRCh38\n")
        handle.write("##contig=<ID=7,length=159345973>\n")
        handle.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n')
        handle.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(sample_ids) + "\n")
        for index, position in enumerate(_variant_positions(), start=1):
            frequency = 0.08 + ((index * 7_919) % 3_700) / 10_000
            alleles = (
                rng.random((len(sample_ids), 2)) < frequency
            ).astype(np.uint8)
            allele_count = int(alleles.sum())
            calls = (
                str(left) + separator + str(right)
                for left, right in alleles
            )
            handle.write(
                f"7\t{int(position)}\tR11B-REG-{index:05d}\tA\tG\t.\tPASS\t"
                f"AC={allele_count};AN={2 * len(sample_ids)}\tGT\t"
            )
            handle.write("\t".join(calls) + "\n")
        handle.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        returncode = process.wait()
    except BaseException:
        process.kill()
        process.wait()
        partial.unlink(missing_ok=True)
        raise
    if returncode != 0:
        partial.unlink(missing_ok=True)
        _fail("r11b_vcf_generation_failed", stderr[-1000:])
    partial.rename(output_path)
    _index_vcf(bcftools, output_path)


def _write_sample_list(path: Path, sample_ids: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sample_ids) + "\n", encoding="utf-8")


def _subset_vcf(
    *,
    source: Path,
    output: Path,
    bcftools: Path,
    sample_list: Path | None = None,
    regions_file: Path | None = None,
) -> None:
    if output.is_file() and Path(f"{output}.tbi").is_file():
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    partial.unlink(missing_ok=True)
    command = [str(bcftools), "view"]
    if sample_list is not None:
        command.extend(("-S", str(sample_list)))
    if regions_file is not None:
        command.extend(("-R", str(regions_file)))
    command.extend(("-Oz", "-o", str(partial), str(source)))
    _run_checked(command)
    partial.rename(output)
    _index_vcf(bcftools, output)


def _query_lines(bcftools: Path, path: Path, format_string: str) -> list[str]:
    completed = subprocess.run(
        [str(bcftools), "query", "-f", format_string, str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        _fail("r11b_bcftools_query_failed", completed.stderr[-1000:])
    return [line for line in completed.stdout.splitlines() if line]


def _sample_names(bcftools: Path, path: Path) -> list[str]:
    completed = subprocess.run(
        [str(bcftools), "query", "-l", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        _fail("r11b_bcftools_samples_failed", completed.stderr[-1000:])
    return [line for line in completed.stdout.splitlines() if line]


def _write_regions(path: Path, rows: Sequence[str], step: int) -> None:
    selected = rows[::step]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(selected) + "\n", encoding="utf-8")


def _write_shapeit_map(path: Path) -> None:
    if path.is_file():
        return
    rows = "pos\tchr\tcM\n" + "".join(
        f"{int(position)}\t7\t{(int(position) - 10_000_000) / 1_000_000:.8f}\n"
        for position in np.linspace(10_000_000, 20_000_000, 27_000, dtype=np.int64)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(rows)


def _prepare_design(repository_root: Path, workspace: Path) -> dict[str, Any]:
    receipt = validate_tracked_receipt(repository_root)
    r07_configuration = _resolve(repository_root, R07_CONFIGURATION_PATH)
    r07_receipt = _resolve(repository_root, R07_RECEIPT_PATH)
    old_projection = _resolve(repository_root, OLD_PROJECTION_RELATIVE)
    tools = _tool_paths(repository_root)
    design = {
        "schema_version": "1.0.0",
        "design_id": "r11b_scientific_resources_v1",
        "scope": "RESOURCE_CONFIRMATION_ONLY",
        "partition": "RESOURCE_ONLY",
        "synthetic_only": True,
        "scientific_performance_measured": False,
        "shapeit_probes": {
            name: {
                "study_samples": study,
                "reference_samples": reference,
                "density": density,
                "variants": VARIANT_COUNT if density == "D_MAX" else ceil(VARIANT_COUNT / 6),
                "threads": 1,
                "seed": 15_052_101,
            }
            for name, (study, reference, density) in PROBE_DIMENSIONS.items()
        },
        "ibd_probes": {
            "samples": 200,
            "pairs_represented": 100,
            "densities": list(IBD_DENSITY_STEP),
            "tools": list(IBD_TOOLS),
            "replicates_per_cell": IBD_REPLICATES,
            "total_calls": len(IBD_DENSITY_STEP) * len(IBD_TOOLS) * IBD_REPLICATES,
        },
        "projection_policy": {
            "one_sided_prediction_level": 0.95,
            "safety_factor": 2.0,
            "shapeit_schedule_calls": 640,
            "ibd_worlds": 90,
            "ibd_tool_density_calls": 540,
        },
        "source_receipt_status": receipt["status"],
        "evidence_sha256": {
            R07_CONFIGURATION_PATH.as_posix(): file_sha256(r07_configuration),
            R07_RECEIPT_PATH.as_posix(): file_sha256(r07_receipt),
            OLD_PROJECTION_RELATIVE.as_posix(): file_sha256(old_projection),
            HARNESS_RELATIVE.as_posix(): file_sha256(
                _resolve(repository_root, HARNESS_RELATIVE)
            ),
            SHAPEIT_RELATIVE.as_posix(): file_sha256(tools["shapeit"]),
            HAP_IBD_RELATIVE.as_posix(): file_sha256(tools["hap_ibd"]),
            REFINED_IBD_RELATIVE.as_posix(): file_sha256(tools["refined_ibd"]),
        },
        "qualification_flags": QUALIFICATION_FLAGS,
    }
    design_path = workspace / "governance/design.json"
    if design_path.is_file():
        frozen_design = load_json(design_path)
        frozen_evidence = dict(frozen_design.get("evidence_sha256", {}))
        current_evidence = dict(design["evidence_sha256"])
        frozen_harness = frozen_evidence.pop(HARNESS_RELATIVE.as_posix(), None)
        current_harness = current_evidence.pop(HARNESS_RELATIVE.as_posix(), None)
        frozen_without_harness = {**frozen_design, "evidence_sha256": frozen_evidence}
        current_without_harness = {**design, "evidence_sha256": current_evidence}
        if frozen_without_harness != current_without_harness:
            _fail("r11b_design_drift")
        if frozen_harness != current_harness:
            if frozen_harness != ORIGINAL_HARNESS_SHA256:
                _fail("r11b_unrecognized_harness_drift")
            amendment_001 = {
                "schema_version": "1.0.0",
                "amendment_id": "R11B-HARNESS-AMENDMENT-001",
                "reason": "POSTPROCESS_INDEX_IDEMPOTENCE_AND_RECEIPT_RESUME_ONLY",
                "scientific_parameters_changed": False,
                "probe_inputs_changed": False,
                "completed_measurements_reused": True,
                "previous_harness_sha256": frozen_harness,
                "current_harness_sha256": AMENDMENT_001_HARNESS_SHA256,
                "qualification_flags": QUALIFICATION_FLAGS,
            }
            amendment_001_path = workspace / "governance/HARNESS_AMENDMENT_001.json"
            if (
                amendment_001_path.is_file()
                and load_json(amendment_001_path) != amendment_001
            ):
                _fail("r11b_harness_amendment_drift")
            _write_json(amendment_001_path, amendment_001)
            if current_harness != AMENDMENT_001_HARNESS_SHA256:
                amendment_002 = {
                    "schema_version": "1.0.0",
                    "amendment_id": "R11B-HARNESS-AMENDMENT-002",
                    "reason": "INTERRUPTED_BRIDGE_NEW_ATTEMPT_ROUTING_ONLY",
                    "scientific_parameters_changed": False,
                    "probe_inputs_changed": False,
                    "completed_measurements_reused": True,
                    "interrupted_measurements_reused": False,
                    "previous_harness_sha256": AMENDMENT_001_HARNESS_SHA256,
                    "current_harness_sha256": current_harness,
                    "qualification_flags": QUALIFICATION_FLAGS,
                }
                amendment_002_path = workspace / "governance/HARNESS_AMENDMENT_002.json"
                if (
                    amendment_002_path.is_file()
                    and load_json(amendment_002_path) != amendment_002
                ):
                    _fail("r11b_harness_amendment_drift")
                _write_json(amendment_002_path, amendment_002)
        return frozen_design
    else:
        _write_json(design_path, design)
    return design


def prepare_workspace(repository_root: Path, workspace: Path) -> dict[str, Any]:
    """Prépare les entrées synthétiques emboîtées et les manifestes de R11b."""
    repository_root = repository_root.resolve(strict=True)
    workspace.mkdir(parents=True, exist_ok=True)
    design = _prepare_design(repository_root, workspace)
    tools = _tool_paths(repository_root)
    bcftools = tools["bcftools"]
    inputs = workspace / "blinded_inputs/RESOURCE_ONLY"
    target = inputs / "shapeit_target"
    _stream_synthetic_vcf(
        output_path=target / "study.vcf.gz",
        sample_ids=_sample_ids("R11B-S", 4_000),
        phased=False,
        seed=91_001,
        bcftools=bcftools,
    )
    _stream_synthetic_vcf(
        output_path=target / "reference.vcf.gz",
        sample_ids=_sample_ids("R11B-R", 3_000),
        phased=True,
        seed=91_002,
        bcftools=bcftools,
    )
    _write_shapeit_map(target / "shapeit.genetic_map.tsv.gz")

    bridge = inputs / "shapeit_bridge"
    study_list = bridge / "study.samples.txt"
    reference_list = bridge / "reference.samples.txt"
    _write_sample_list(study_list, _sample_ids("R11B-S", 2_000))
    _write_sample_list(reference_list, _sample_ids("R11B-R", 1_500))
    _subset_vcf(
        source=target / "study.vcf.gz",
        output=bridge / "study.vcf.gz",
        bcftools=bcftools,
        sample_list=study_list,
    )
    _subset_vcf(
        source=target / "reference.vcf.gz",
        output=bridge / "reference.vcf.gz",
        bcftools=bcftools,
        sample_list=reference_list,
    )
    _write_shapeit_map(bridge / "shapeit.genetic_map.tsv.gz")

    n250 = _resolve(repository_root, N250_PROBE_RELATIVE) / "blinded_inputs/CALIBRATION/regional"
    dlow = inputs / "shapeit_dlow_n250"
    reference_rows = _query_lines(
        bcftools, n250 / "reference.ac_an.vcf.gz", "%CHROM\t%POS\n"
    )
    dlow_regions = dlow / "dlow.regions.tsv"
    _write_regions(dlow_regions, reference_rows, 6)
    _subset_vcf(
        source=n250 / "study.ac_an.vcf.gz",
        output=dlow / "study.vcf.gz",
        bcftools=bcftools,
        regions_file=dlow_regions,
    )
    _subset_vcf(
        source=n250 / "reference.ac_an.vcf.gz",
        output=dlow / "reference.vcf.gz",
        bcftools=bcftools,
        regions_file=dlow_regions,
    )
    _write_shapeit_map(dlow / "shapeit.genetic_map.tsv.gz")

    ibd_source = _resolve(repository_root, SOURCE_PROBE_RELATIVE) / "outputs/shapeit5/common.bcf"
    ibd_root = inputs / "ibd_batch_100"
    source_samples = _sample_names(bcftools, ibd_source)
    if len(source_samples) < 200:
        _fail("r11b_ibd_source_sample_count")
    ibd_samples = ibd_root / "samples.txt"
    _write_sample_list(ibd_samples, source_samples[:200])
    source_positions = _query_lines(bcftools, ibd_source, "%CHROM\t%POS\n")
    for density, step in IBD_DENSITY_STEP.items():
        density_root = ibd_root / density
        regions = density_root / "regions.tsv"
        _write_regions(regions, source_positions, step)
        vcf = density_root / "phased.vcf.gz"
        _subset_vcf(
            source=ibd_source,
            output=vcf,
            bcftools=bcftools,
            sample_list=ibd_samples,
            regions_file=regions,
        )
        variants = _query_lines(bcftools, vcf, "%CHROM\t%ID\t%POS\n")
        map_path = density_root / "ibd.map"
        map_path.write_text(
            "".join(
                f"{chromosome}\t{variant_id}\t{(int(position) - 10_000_000) / 1_000_000:.8f}\t{position}\n"
                for chromosome, variant_id, position in (row.split("\t") for row in variants)
            ),
            encoding="utf-8",
        )

    inventory_paths = sorted(
        path for path in inputs.rglob("*") if path.is_file()
    )
    inventory = {
        "schema_version": "1.0.0",
        "design_sha256": file_sha256(workspace / "governance/design.json"),
        "artifacts": [
            {
                "relative_path": path.relative_to(workspace).as_posix(),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in inventory_paths
        ],
        "qualification_flags": QUALIFICATION_FLAGS,
    }
    _write_json(workspace / "governance/input_inventory.json", inventory)
    return {"design": design, "artifact_count": len(inventory_paths)}


def _parse_time_cpu(resource_log: Path) -> tuple[float, float]:
    text = resource_log.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"([0-9.]+) real\s+([0-9.]+) user\s+([0-9.]+) sys",
        text,
    )
    if match is None:
        _fail("r11b_resource_cpu_missing", str(resource_log))
    return float(match.group(2)), float(match.group(3))


def _shapeit_probe_root(workspace: Path, probe_name: str) -> Path:
    attempt = SHAPEIT_ATTEMPT_NUMBER[probe_name]
    return workspace / "probes" / f"shapeit_{probe_name}_attempt{attempt:02d}"


def _shapeit_inputs(workspace: Path, probe_name: str) -> Path:
    mapping = {
        "target": "shapeit_target",
        "bridge": "shapeit_bridge",
        "dlow": "shapeit_dlow_n250",
    }
    return workspace / "blinded_inputs/RESOURCE_ONLY" / mapping[probe_name]


def _count_vcf(bcftools: Path, path: Path) -> tuple[int, int]:
    samples = _sample_names(bcftools, path)
    positions = _query_lines(bcftools, path, "%POS\n")
    return len(samples), len(positions)


def run_shapeit_probe(
    repository_root: Path,
    workspace: Path,
    probe_name: str,
) -> dict[str, Any]:
    """Exécute une sonde SHAPEIT5 unique et reprend un résultat déjà clos."""
    if probe_name not in PROBE_DIMENSIONS:
        _fail("r11b_shapeit_probe_unknown", probe_name)
    repository_root = repository_root.resolve(strict=True)
    _prepare_design(repository_root, workspace)
    tools = _tool_paths(repository_root)
    probe_root = _shapeit_probe_root(workspace, probe_name)
    receipt_path = probe_root / "governance/measurement.json"
    if receipt_path.is_file():
        return load_json(receipt_path)
    inputs = _shapeit_inputs(workspace, probe_name)
    for name in ("study.vcf.gz", "reference.vcf.gz", "shapeit.genetic_map.tsv.gz"):
        if not (inputs / name).is_file():
            _fail("r11b_inputs_not_prepared", name)
    run_root = probe_root / "run"
    output_root = probe_root / "outputs"
    output_bcf = output_root / "common.bcf"
    output_log = output_root / "common.shapeit.log"
    command = [
        str(tools["shapeit"]),
        "--input", str(inputs / "study.vcf.gz"),
        "--reference", str(inputs / "reference.vcf.gz"),
        "--map", str(inputs / "shapeit.genetic_map.tsv.gz"),
        "--region", REGION,
        "--output", str(output_bcf),
        "--output-format", "bcf",
        "--log", str(output_log),
        "--thread", "1",
        "--seed", "15052101",
    ]
    measure_path = run_root / "commands" / f"shapeit_{probe_name}" / "measure.json"
    if probe_root.exists():
        if not measure_path.is_file() or not output_bcf.is_file():
            _fail("r11b_incomplete_probe_requires_new_attempt", probe_name)
        measure = CommandMeasure(**load_json(measure_path))
    else:
        probe_root.mkdir(parents=True, exist_ok=False)
        output_root.mkdir(parents=True)
        measure = _run_measured(
            command_id=f"shapeit_{probe_name}",
            tool="shapeit5_phase_common",
            command=command,
            run_root=run_root,
            timeout_seconds=168 * 3_600,
        )
    if measure.returncode != 0 or not output_bcf.is_file():
        _fail("r11b_shapeit_probe_failed", probe_name)
    if not Path(f"{output_bcf}.csi").is_file():
        _run_checked([str(tools["bcftools"]), "index", "--csi", str(output_bcf)])
    user_cpu, system_cpu = _parse_time_cpu(run_root / measure.resource_path)
    study_count, input_variants = _count_vcf(tools["bcftools"], inputs / "study.vcf.gz")
    reference_count, _ = _count_vcf(tools["bcftools"], inputs / "reference.vcf.gz")
    output_samples, output_variants = _count_vcf(tools["bcftools"], output_bcf)
    receipt = {
        "schema_version": "1.0.0",
        "probe_id": (
            f"r11b_shapeit_{probe_name}_attempt"
            f"{SHAPEIT_ATTEMPT_NUMBER[probe_name]:02d}"
        ),
        "state": "COMPLETED",
        "completed_at": _utc_now(),
        "scope": "RESOURCE_CONFIRMATION_ONLY",
        "partition": "RESOURCE_ONLY",
        "synthetic_only": True,
        "scientific_performance_measured": False,
        "counts": {
            "study_samples": study_count,
            "reference_samples": reference_count,
            "input_variants": input_variants,
            "output_samples": output_samples,
            "output_variants": output_variants,
            "threads": 1,
        },
        "measurement": {
            "wall_seconds": measure.wall_seconds,
            "user_cpu_seconds": user_cpu,
            "system_cpu_seconds": system_cpu,
            "maximum_rss_bytes": measure.maximum_rss_bytes,
        },
        "outputs": [
            {
                "relative_path": path.relative_to(probe_root).as_posix(),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in (output_bcf, Path(f"{output_bcf}.csi"), output_log)
        ],
        "qualification_flags": QUALIFICATION_FLAGS,
    }
    _write_json(receipt_path, receipt)
    return receipt


def _ibd_command(
    *,
    tool: str,
    tools: Mapping[str, Path],
    vcf: Path,
    map_path: Path,
    output_prefix: Path,
) -> list[str]:
    jar = tools["hap_ibd"] if tool == "HAP_IBD" else tools["refined_ibd"]
    command = [
        str(tools["java"]), "-Xmx4096m", "-jar", str(jar),
        f"gt={vcf}", f"map={map_path}", f"out={output_prefix}",
    ]
    if tool == "HAP_IBD":
        command.extend(("min-seed=1.0", "min-output=1.0", "min-markers=20", "nthreads=4"))
    else:
        command.extend(("length=1.0", "nthreads=4"))
    return command


def run_ibd_probes(repository_root: Path, workspace: Path) -> dict[str, Any]:
    """Mesure trois répétitions par outil et densité sur un lot de 100 paires."""
    repository_root = repository_root.resolve(strict=True)
    _prepare_design(repository_root, workspace)
    tools = _tool_paths(repository_root)
    results: list[dict[str, Any]] = []
    for density in IBD_DENSITY_STEP:
        inputs = workspace / "blinded_inputs/RESOURCE_ONLY/ibd_batch_100" / density
        vcf = inputs / "phased.vcf.gz"
        map_path = inputs / "ibd.map"
        if not vcf.is_file() or not map_path.is_file():
            _fail("r11b_inputs_not_prepared", density)
        sample_count, variant_count = _count_vcf(tools["bcftools"], vcf)
        for tool in IBD_TOOLS:
            for replicate in range(1, IBD_REPLICATES + 1):
                probe_id = f"r11b_ibd_{density.lower()}_{tool.lower()}_rep{replicate:02d}"
                probe_root = workspace / "probes" / probe_id
                receipt_path = probe_root / "governance/measurement.json"
                if receipt_path.is_file():
                    results.append(load_json(receipt_path))
                    continue
                probe_root.mkdir(parents=True, exist_ok=False)
                output_root = probe_root / "outputs"
                output_root.mkdir(parents=True)
                output_prefix = output_root / "segments"
                command = _ibd_command(
                    tool=tool,
                    tools=tools,
                    vcf=vcf,
                    map_path=map_path,
                    output_prefix=output_prefix,
                )
                measure = _run_measured(
                    command_id=probe_id,
                    tool=tool,
                    command=command,
                    run_root=probe_root / "run",
                    timeout_seconds=168 * 3_600,
                )
                output_path = Path(f"{output_prefix}.ibd.gz")
                if measure.returncode != 0 or not output_path.is_file():
                    _fail("r11b_ibd_probe_failed", probe_id)
                user_cpu, system_cpu = _parse_time_cpu(
                    probe_root / "run" / measure.resource_path
                )
                receipt = {
                    "schema_version": "1.0.0",
                    "probe_id": probe_id,
                    "state": "COMPLETED",
                    "completed_at": _utc_now(),
                    "scope": "RESOURCE_CONFIRMATION_ONLY",
                    "partition": "RESOURCE_ONLY",
                    "synthetic_only": True,
                    "scientific_performance_measured": False,
                    "tool": tool,
                    "density": density,
                    "replicate": replicate,
                    "counts": {
                        "samples": sample_count,
                        "pairs_represented": 100,
                        "variants": variant_count,
                        "threads": 4,
                    },
                    "measurement": {
                        "wall_seconds": measure.wall_seconds,
                        "user_cpu_seconds": user_cpu,
                        "system_cpu_seconds": system_cpu,
                        "maximum_rss_bytes": measure.maximum_rss_bytes,
                    },
                    "outputs": [
                        {
                            "relative_path": output_path.relative_to(probe_root).as_posix(),
                            "sha256": file_sha256(output_path),
                            "size_bytes": output_path.stat().st_size,
                        }
                    ],
                    "qualification_flags": QUALIFICATION_FLAGS,
                }
                _write_json(receipt_path, receipt)
                results.append(receipt)
    summary = {
        "schema_version": "1.0.0",
        "state": "COMPLETED",
        "probe_count": len(results),
        "expected_probe_count": 18,
        "measurements": [
            {
                "probe_id": item["probe_id"],
                "tool": item["tool"],
                "density": item["density"],
                "replicate": item["replicate"],
                **item["measurement"],
            }
            for item in results
        ],
        "qualification_flags": QUALIFICATION_FLAGS,
    }
    if len(results) != 18:
        _fail("r11b_ibd_probe_count_mismatch")
    _write_json(workspace / "evaluations/ibd_batch_measurements.json", summary)
    return summary


def _one_sided_prediction_bound(values: Sequence[float]) -> float:
    if len(values) < 3 or any(value <= 0 for value in values):
        _fail("r11b_prediction_values_invalid")
    logged = np.log(np.asarray(values, dtype=float))
    standard_deviation = float(logged.std(ddof=1))
    quantile = float(student_t.ppf(0.95, len(values) - 1))
    return exp(float(logged.mean()) + quantile * standard_deviation * (1 + 1 / len(values)) ** 0.5)


def _shapeit_fit_bound(
    totals: Sequence[int],
    values: Sequence[float],
    target_total: int,
) -> float:
    if len(totals) != len(values) or len(values) < 5:
        _fail("r11b_shapeit_fit_invalid")
    x = np.log(np.asarray(totals, dtype=float))
    y = np.log(np.asarray(values, dtype=float))
    matrix = np.column_stack((np.ones_like(x), x))
    coefficients, *_ = np.linalg.lstsq(matrix, y, rcond=None)
    residuals = y - matrix @ coefficients
    degrees = len(y) - 2
    residual_standard_error = float(np.sqrt(np.sum(residuals**2) / degrees))
    target_vector = np.asarray([1.0, log(target_total)])
    leverage = float(target_vector @ np.linalg.inv(matrix.T @ matrix) @ target_vector)
    quantile = float(student_t.ppf(0.95, degrees))
    prediction = float(target_vector @ coefficients)
    return exp(prediction + quantile * residual_standard_error * (1 + leverage) ** 0.5)


def _load_shapeit_receipts(workspace: Path) -> dict[str, dict[str, Any]]:
    receipts: dict[str, dict[str, Any]] = {}
    for name in PROBE_DIMENSIONS:
        path = _shapeit_probe_root(workspace, name) / "governance/measurement.json"
        if not path.is_file():
            _fail("r11b_shapeit_measurement_missing", name)
        receipts[name] = load_json(path)
    return receipts


def _load_ibd_receipts(workspace: Path) -> list[dict[str, Any]]:
    paths = sorted((workspace / "probes").glob("r11b_ibd_*/governance/measurement.json"))
    if len(paths) != 18:
        _fail("r11b_ibd_measurement_missing", str(len(paths)))
    return [load_json(path) for path in paths]


def build_projection(repository_root: Path, workspace: Path) -> dict[str, Any]:
    """Projette le calendrier R07 v2 et décide seulement la porte de ressources."""
    repository_root = repository_root.resolve(strict=True)
    r07 = load_yaml(_resolve(repository_root, R07_CONFIGURATION_PATH))
    previous = load_json(_resolve(repository_root, OLD_PROJECTION_RELATIVE))
    shapeit = _load_shapeit_receipts(workspace)
    ibd = _load_ibd_receipts(workspace)
    old_totals = [item["total_samples"] for item in previous["equal_panel_observations"]]
    bridge_total = sum(PROBE_DIMENSIONS["bridge"][:2])
    target_total = sum(PROBE_DIMENSIONS["target"][:2])
    totals = [*old_totals, bridge_total, target_total]

    def resource_values(field: str) -> list[float]:
        old_field = "cpu_seconds" if field == "cpu_seconds" else field
        old_values = [float(item[old_field]) for item in previous["equal_panel_observations"]]
        for name in ("bridge", "target"):
            measurement = shapeit[name]["measurement"]
            value = (
                measurement["user_cpu_seconds"] + measurement["system_cpu_seconds"]
                if field == "cpu_seconds"
                else measurement[field]
            )
            old_values.append(float(value))
        return old_values

    shapeit_wall_bound = _shapeit_fit_bound(totals, resource_values("wall_seconds"), target_total)
    shapeit_cpu_bound = _shapeit_fit_bound(totals, resource_values("cpu_seconds"), target_total)
    shapeit_rss_bound = _shapeit_fit_bound(
        totals, resource_values("maximum_rss_bytes"), target_total
    )
    safety_factor = 2.0
    dmid_cpu_ratio = float(previous["density_control"]["cpu_seconds_ratio"])
    n250_cpu = float(previous["equal_panel_observations"][2]["cpu_seconds"])
    dlow_measurement = shapeit["dlow"]["measurement"]
    dlow_cpu = dlow_measurement["user_cpu_seconds"] + dlow_measurement["system_cpu_seconds"]
    dlow_cpu_ratio = dlow_cpu / n250_cpu
    expected_calls = r07["phasing_schedule"]["expected_calls"]
    safe_target_cpu_hours = shapeit_cpu_bound * safety_factor / 3_600
    phasing_core_hours = safe_target_cpu_hours * (
        expected_calls["D_MAX"]
        + expected_calls["D_MID"] * dmid_cpu_ratio
        + expected_calls["D_LOW"] * dlow_cpu_ratio
    )
    safe_target_wall_hours = shapeit_wall_bound * safety_factor / 3_600
    phasing_serial_wall_hours = safe_target_wall_hours * (
        expected_calls["D_MAX"]
        + expected_calls["D_MID"] * dmid_cpu_ratio
        + expected_calls["D_LOW"] * dlow_cpu_ratio
    )

    ibd_cells: list[dict[str, Any]] = []
    ibd_core_hours = 0.0
    ibd_serial_wall_hours = 0.0
    ibd_max_rss = 0.0
    for density in IBD_DENSITY_STEP:
        for tool in IBD_TOOLS:
            cell = [item for item in ibd if item["density"] == density and item["tool"] == tool]
            wall_values = [float(item["measurement"]["wall_seconds"]) for item in cell]
            cpu_values = [
                float(item["measurement"]["user_cpu_seconds"])
                + float(item["measurement"]["system_cpu_seconds"])
                for item in cell
            ]
            rss_values = [float(item["measurement"]["maximum_rss_bytes"]) for item in cell]
            wall_bound = _one_sided_prediction_bound(wall_values) * safety_factor
            cpu_bound = _one_sided_prediction_bound(cpu_values) * safety_factor
            rss_bound = _one_sided_prediction_bound(rss_values) * safety_factor
            calls = int(r07["ibd_worlds"]["expected_independent_worlds"]["total"])
            ibd_core_hours += cpu_bound * calls / 3_600
            ibd_serial_wall_hours += wall_bound * calls / 3_600
            ibd_max_rss = max(ibd_max_rss, rss_bound)
            ibd_cells.append(
                {
                    "density": density,
                    "tool": tool,
                    "replicates": len(cell),
                    "safe_wall_seconds_per_call": wall_bound,
                    "safe_cpu_seconds_per_call": cpu_bound,
                    "safe_rss_bytes_per_call": rss_bound,
                    "projected_calls": calls,
                }
            )

    workers = int(r07["resources"]["r2_workers_simultaneous_max"])
    schedule_days = (phasing_serial_wall_hours + ibd_serial_wall_hours) / workers / 24
    safe_shapeit_rss = shapeit_rss_bound * safety_factor
    aggregate_rss_gib = max(safe_shapeit_rss, ibd_max_rss) * workers / 2**30

    inputs_size = sum(
        path.stat().st_size
        for path in (workspace / "blinded_inputs").rglob("*")
        if path.is_file()
    )
    target_output_size = sum(item["size_bytes"] for item in shapeit["target"]["outputs"])
    dlow_output_size = sum(item["size_bytes"] for item in shapeit["dlow"]["outputs"])
    dmax_mid_output_projection = target_output_size * (
        expected_calls["D_MAX"] + expected_calls["D_MID"] * 0.5
    )
    dlow_output_projection = dlow_output_size * expected_calls["D_LOW"] * target_total / 500
    ibd_output_projection = sum(item["outputs"][0]["size_bytes"] for item in ibd) / 3 * 90
    storage_gib = (
        2 * (inputs_size * 22 + dmax_mid_output_projection + dlow_output_projection + ibd_output_projection)
        / 2**30
    )

    resources = r07["resources"]
    checks = [
        {
            "check": "SHAPEIT_CORE_HOURS",
            "projected": phasing_core_hours,
            "ceiling": resources["core_hour_allocations"]["phasing"],
            "status": "PASS" if phasing_core_hours <= resources["core_hour_allocations"]["phasing"] else "RESOURCE_BLOCKED",
        },
        {
            "check": "IBD_CORE_HOURS",
            "projected": ibd_core_hours,
            "ceiling": resources["core_hour_allocations"]["explicit_ibd"],
            "status": "PASS" if ibd_core_hours <= resources["core_hour_allocations"]["explicit_ibd"] else "RESOURCE_BLOCKED",
        },
        {
            "check": "R2_AGGREGATE_RSS_GIB_AT_24_WORKERS",
            "projected": aggregate_rss_gib,
            "ceiling": resources["aggregate_rss_gib_max"],
            "status": "PASS" if aggregate_rss_gib <= resources["aggregate_rss_gib_max"] else "RESOURCE_BLOCKED",
        },
        {
            "check": "QUALIFICATION_SCHEDULE_DAYS_AT_24_WORKERS",
            "projected": schedule_days,
            "ceiling": resources["qualification_calendar_days_max"],
            "status": "PASS" if schedule_days <= resources["qualification_calendar_days_max"] else "RESOURCE_BLOCKED",
        },
        {
            "check": "FULL_RELEASE_STORAGE_GIB",
            "projected": storage_gib,
            "ceiling": resources["storage_gib_max"],
            "status": "PASS" if storage_gib <= resources["storage_gib_max"] else "RESOURCE_BLOCKED",
        },
    ]
    resource_ready = all(item["status"] == "PASS" for item in checks)
    decision = {
        "schema_version": "1.0.0",
        "decision_id": "R11B-SCIENTIFIC-RESOURCE-CONFIRMATION-V1",
        "created_at": _utc_now(),
        "scope": "RESOURCE_CONFIRMATION_ONLY",
        "status": "RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE" if resource_ready else "RESOURCE_BLOCKED",
        "scientific_performance_measured": False,
        "shapeit": {
            "safe_wall_hours_per_target_call": safe_target_wall_hours,
            "safe_cpu_hours_per_target_call": safe_target_cpu_hours,
            "safe_rss_gib_per_target_call": safe_shapeit_rss / 2**30,
            "dlow_cpu_ratio": dlow_cpu_ratio,
            "projected_core_hours": phasing_core_hours,
            "projected_serial_wall_hours": phasing_serial_wall_hours,
        },
        "ibd": {
            "cells": ibd_cells,
            "projected_core_hours": ibd_core_hours,
            "projected_serial_wall_hours": ibd_serial_wall_hours,
        },
        "schedule": {
            "workers": workers,
            "projected_days": schedule_days,
            "execution_host_logical_cpus": os.cpu_count() or 1,
            "execution_host_can_supply_24_workers": (os.cpu_count() or 1) >= workers,
        },
        "storage": {"projected_gib_with_factor_2": storage_gib},
        "checks": checks,
        "scientific_g5_authorized": False,
        "next_phase": "BUILD_NEW_SCIENTIFIC_RELEASE" if resource_ready else "REVISE_R07_OR_RESOURCES",
        "restrictions": [
            "RESOURCE_EVIDENCE_ONLY",
            "NEW_RELEASE_FREEZE_AND_BLIND_SCIENTIFIC_VALIDATION_REQUIRED",
            "NO_REAL_DATA_AUTHORIZED",
        ],
        "qualification_flags": QUALIFICATION_FLAGS,
    }
    _write_json(workspace / "evaluations/R11B_RESOURCE_DECISION.json", decision)
    return decision


def status(repository_root: Path, workspace: Path) -> dict[str, Any]:
    """Retourne l'avancement sans modifier le workspace."""
    repository_root = repository_root.resolve(strict=True)
    prepared = (workspace / "governance/input_inventory.json").is_file()
    shapeit = {
        name: (_shapeit_probe_root(workspace, name) / "governance/measurement.json").is_file()
        for name in PROBE_DIMENSIONS
    }
    ibd_count = len(list((workspace / "probes").glob("r11b_ibd_*/governance/measurement.json")))
    decision_path = workspace / "evaluations/R11B_RESOURCE_DECISION.json"
    try:
        display_workspace = workspace.relative_to(repository_root).as_posix()
    except ValueError:
        display_workspace = workspace.as_posix()
    return {
        "workspace": display_workspace,
        "prepared": prepared,
        "shapeit": shapeit,
        "ibd_completed": ibd_count,
        "ibd_expected": 18,
        "decision": load_json(decision_path) if decision_path.is_file() else None,
        "qualification_flags": QUALIFICATION_FLAGS,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--workspace", type=Path)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("prepare")
    shapeit = subparsers.add_parser("run-shapeit")
    shapeit.add_argument("--probe", choices=tuple(PROBE_DIMENSIONS), required=True)
    subparsers.add_parser("run-ibd")
    subparsers.add_parser("project")
    subparsers.add_parser("status")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    repository_root = arguments.repository_root.resolve(strict=True)
    workspace = (
        arguments.workspace.resolve()
        if arguments.workspace is not None
        else repository_root / WORKSPACE_RELATIVE
    )
    if arguments.command == "prepare":
        document = prepare_workspace(repository_root, workspace)
    elif arguments.command == "run-shapeit":
        document = run_shapeit_probe(repository_root, workspace, arguments.probe)
    elif arguments.command == "run-ibd":
        document = run_ibd_probes(repository_root, workspace)
    elif arguments.command == "project":
        document = build_projection(repository_root, workspace)
    else:
        document = status(repository_root, workspace)
    print(json.dumps(document, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
