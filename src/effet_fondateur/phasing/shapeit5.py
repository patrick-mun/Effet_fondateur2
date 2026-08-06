"""Contrat reproductible de l'adaptateur SHAPEIT5 retenu à l'étape 12.0."""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


SHAPEIT5_VERSION_PATTERN = re.compile(r"\*\s+Version\s+:\s+(\d+\.\d+\.\d+)")


class Shapeit5ContractError(ValueError):
    """Signale une configuration ou une installation SHAPEIT5 non conforme."""


@dataclass(frozen=True)
class Shapeit5Contract:
    """Décrit les choix scientifiques et techniques figés pour SHAPEIT5."""

    adapter_id: str
    software: str
    version: str
    license: str
    supported_assemblies: tuple[str, ...]
    supported_chromosomes: tuple[int, ...]
    input_formats: tuple[str, ...]
    reference_formats: tuple[str, ...]
    output_format: str
    pedigree_format: str
    genetic_map_required: bool
    reference_panel_required: bool
    preserves_input_only_variants_in_common_phase: bool
    rare_variant_strategy: str
    default_seed: int
    default_effective_size: int


SHAPEIT5_CONTRACT = Shapeit5Contract(
    adapter_id="shapeit5_phase_common_rare_v1",
    software="SHAPEIT5",
    version="5.1.1",
    license="MIT",
    supported_assemblies=("GRCh38",),
    supported_chromosomes=tuple(range(1, 23)),
    input_formats=("VCF", "VCF_GZ", "BCF"),
    reference_formats=("VCF", "VCF_GZ", "BCF"),
    output_format="BCF",
    pedigree_format="CHILD_FATHER_MOTHER_TSV",
    genetic_map_required=True,
    reference_panel_required=True,
    preserves_input_only_variants_in_common_phase=False,
    rare_variant_strategy="PHASE_RARE_ON_PHASE_COMMON_SCAFFOLD",
    default_seed=15_052_011,
    default_effective_size=15_000,
)


@dataclass(frozen=True)
class Shapeit5AdapterConfig:
    """Commandes et version attendue déclarées dans la configuration du run."""

    adapter_id: str
    phase_common_command: str
    phase_rare_command: str
    expected_version: str


@dataclass(frozen=True)
class Shapeit5Probe:
    """Résultat non sensible du contrôle des deux exécutables SHAPEIT5."""

    phase_common_path: str
    phase_rare_path: str
    phase_common_version: str
    phase_rare_version: str


def parse_shapeit5_adapter_config(
    raw_config: Mapping[str, Any] | None,
) -> Shapeit5AdapterConfig:
    """Valide la déclaration de l'adaptateur contre le contrat figé."""
    if raw_config is None:
        raise Shapeit5ContractError("phasing_adapter_not_configured")
    required_fields = {
        "adapter_id",
        "phase_common_command",
        "phase_rare_command",
        "expected_version",
    }
    if set(raw_config) != required_fields:
        raise Shapeit5ContractError("invalid_phasing_adapter_fields")
    if not all(
        isinstance(raw_config[field], str) and raw_config[field].strip()
        for field in required_fields
    ):
        raise Shapeit5ContractError("invalid_phasing_adapter_value")
    if raw_config["adapter_id"] != SHAPEIT5_CONTRACT.adapter_id:
        raise Shapeit5ContractError("unsupported_phasing_adapter")
    if raw_config["expected_version"] != SHAPEIT5_CONTRACT.version:
        raise Shapeit5ContractError("unsupported_shapeit5_version")
    return Shapeit5AdapterConfig(
        adapter_id=raw_config["adapter_id"],
        phase_common_command=raw_config["phase_common_command"],
        phase_rare_command=raw_config["phase_rare_command"],
        expected_version=raw_config["expected_version"],
    )


def _resolve_executable(command: str, component: str) -> str:
    executable = shutil.which(command)
    if executable is None:
        raise Shapeit5ContractError(f"shapeit5_executable_not_available:{component}")
    return executable


def _read_version(executable: str, component: str, timeout_seconds: int) -> str:
    try:
        completed = subprocess.run(
            [executable, "--help"],
            capture_output=True,
            check=False,
            text=True,
            timeout=timeout_seconds,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise Shapeit5ContractError(
            f"shapeit5_probe_failed:{component}"
        ) from error
    output = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
    match = SHAPEIT5_VERSION_PATTERN.search(output)
    if completed.returncode != 0 or match is None:
        raise Shapeit5ContractError(f"shapeit5_version_unreadable:{component}")
    return match.group(1)


def probe_shapeit5(
    config: Shapeit5AdapterConfig,
    *,
    timeout_seconds: int = 10,
) -> Shapeit5Probe:
    """Résout les exécutables et exige la version exacte du contrat."""
    if isinstance(timeout_seconds, bool) or timeout_seconds <= 0:
        raise Shapeit5ContractError("invalid_shapeit5_probe_timeout")
    phase_common_path = _resolve_executable(
        config.phase_common_command, "phase_common"
    )
    phase_rare_path = _resolve_executable(config.phase_rare_command, "phase_rare")
    common_version = _read_version(
        phase_common_path, "phase_common", timeout_seconds
    )
    rare_version = _read_version(phase_rare_path, "phase_rare", timeout_seconds)
    if common_version != config.expected_version:
        raise Shapeit5ContractError("shapeit5_version_mismatch:phase_common")
    if rare_version != config.expected_version:
        raise Shapeit5ContractError("shapeit5_version_mismatch:phase_rare")
    return Shapeit5Probe(
        phase_common_path=phase_common_path,
        phase_rare_path=phase_rare_path,
        phase_common_version=common_version,
        phase_rare_version=rare_version,
    )


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise Shapeit5ContractError(f"invalid_shapeit5_parameter:{name}")
    return value


def _validate_region(region: str) -> str:
    if not isinstance(region, str) or not re.fullmatch(
        r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2])(?::[1-9][0-9]*-[1-9][0-9]*)?",
        region,
    ):
        raise Shapeit5ContractError("invalid_shapeit5_region")
    return region


def _validate_input_path(path: Path, name: str) -> Path:
    path_text = str(path).lower()
    if not path_text.endswith((".bcf", ".vcf", ".vcf.gz")):
        raise Shapeit5ContractError(f"invalid_shapeit5_format:{name}")
    return path


def _validate_output_path(path: Path) -> Path:
    if path.suffix.lower() != ".bcf":
        raise Shapeit5ContractError("invalid_shapeit5_output_format")
    return path


def build_phase_common_command(
    probe: Shapeit5Probe,
    *,
    input_path: Path,
    reference_path: Path,
    genetic_map_path: Path,
    region: str,
    output_path: Path,
    log_path: Path,
    pedigree_path: Path | None = None,
    threads: int = 1,
    seed: int = SHAPEIT5_CONTRACT.default_seed,
) -> list[str]:
    """Construit la commande du scaffold commun sans l'exécuter."""
    command = [
        probe.phase_common_path,
        "--input",
        str(_validate_input_path(input_path, "input")),
        "--reference",
        str(_validate_input_path(reference_path, "reference")),
        "--map",
        str(genetic_map_path),
        "--region",
        _validate_region(region),
        "--output",
        str(_validate_output_path(output_path)),
        "--output-format",
        "bcf",
        "--log",
        str(log_path),
        "--thread",
        str(_positive_integer(threads, "threads")),
        "--seed",
        str(_positive_integer(seed, "seed")),
    ]
    if pedigree_path is not None:
        command.extend(("--pedigree", str(pedigree_path)))
    return command


def build_phase_rare_command(
    probe: Shapeit5Probe,
    *,
    input_path: Path,
    scaffold_path: Path,
    genetic_map_path: Path,
    input_region: str,
    scaffold_region: str,
    output_path: Path,
    log_path: Path,
    pedigree_path: Path | None = None,
    threads: int = 1,
    seed: int = SHAPEIT5_CONTRACT.default_seed,
    effective_size: int = SHAPEIT5_CONTRACT.default_effective_size,
) -> list[str]:
    """Construit la commande qui replace les variants rares sur le scaffold."""
    command = [
        probe.phase_rare_path,
        "--input",
        str(_validate_input_path(input_path, "input")),
        "--scaffold",
        str(_validate_input_path(scaffold_path, "scaffold")),
        "--map",
        str(genetic_map_path),
        "--input-region",
        _validate_region(input_region),
        "--scaffold-region",
        _validate_region(scaffold_region),
        "--output",
        str(_validate_output_path(output_path)),
        "--log",
        str(log_path),
        "--thread",
        str(_positive_integer(threads, "threads")),
        "--seed",
        str(_positive_integer(seed, "seed")),
        "--effective-size",
        str(_positive_integer(effective_size, "effective_size")),
    ]
    if pedigree_path is not None:
        command.extend(("--pedigree", str(pedigree_path)))
    return command
