import sys
from pathlib import Path

import pytest

from effet_fondateur.contracts import ConfigurationError, load_pipeline_config
from effet_fondateur.orchestrator.environment import build_environment
from effet_fondateur.phasing import (
    SHAPEIT5_CONTRACT,
    Shapeit5ContractError,
    build_phase_common_command,
    build_phase_rare_command,
    parse_shapeit5_adapter_config,
    probe_shapeit5,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def write_shapeit_executable(path: Path, component: str, version: str = "5.1.1") -> None:
    path.write_text(
        f'''#!{sys.executable}
import sys

if sys.argv[1:] == ["--help"]:
    print("[SHAPEIT5] {component}")
    print("  * Version       : {version} / commit = test")
    raise SystemExit(0)
raise SystemExit(2)
''',
        encoding="utf-8",
    )
    path.chmod(0o755)


def adapter_config(common: Path, rare: Path) -> dict[str, str]:
    return {
        "adapter_id": SHAPEIT5_CONTRACT.adapter_id,
        "phase_common_command": str(common),
        "phase_rare_command": str(rare),
        "expected_version": SHAPEIT5_CONTRACT.version,
    }


def probe_fixture(tmp_path: Path):
    common = tmp_path / "SHAPEIT5_phase_common"
    rare = tmp_path / "SHAPEIT5_phase_rare"
    write_shapeit_executable(common, "phase_common")
    write_shapeit_executable(rare, "phase_rare")
    return probe_shapeit5(parse_shapeit5_adapter_config(adapter_config(common, rare)))


def test_example_declares_frozen_shapeit5_contract() -> None:
    config = load_pipeline_config(REPOSITORY_ROOT / "config" / "pipeline.example.yaml")

    adapter = parse_shapeit5_adapter_config(config["tools"]["phasing_adapter"])

    assert adapter.adapter_id == "shapeit5_phase_common_rare_v1"
    assert adapter.expected_version == "5.1.1"
    assert SHAPEIT5_CONTRACT.license == "MIT"
    assert SHAPEIT5_CONTRACT.supported_assemblies == ("GRCh38",)
    assert SHAPEIT5_CONTRACT.genetic_map_required is True
    assert SHAPEIT5_CONTRACT.preserves_input_only_variants_in_common_phase is False


def test_pipeline_schema_rejects_unpinned_shapeit_version(tmp_path: Path) -> None:
    config_text = (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
        encoding="utf-8"
    )
    config_path = tmp_path / "invalid.yaml"
    config_path.write_text(
        config_text.replace("expected_version: 5.1.1", "expected_version: latest"),
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="expected_version"):
        load_pipeline_config(config_path)


def test_probe_requires_exact_versions_for_both_components(tmp_path: Path) -> None:
    common = tmp_path / "SHAPEIT5_phase_common"
    rare = tmp_path / "SHAPEIT5_phase_rare"
    write_shapeit_executable(common, "phase_common")
    write_shapeit_executable(rare, "phase_rare", version="5.1.0")

    with pytest.raises(
        Shapeit5ContractError, match="shapeit5_version_mismatch:phase_rare"
    ):
        probe_shapeit5(parse_shapeit5_adapter_config(adapter_config(common, rare)))


def test_environment_records_both_shapeit_versions(tmp_path: Path) -> None:
    common = tmp_path / "SHAPEIT5_phase_common"
    rare = tmp_path / "SHAPEIT5_phase_rare"
    write_shapeit_executable(common, "phase_common")
    write_shapeit_executable(rare, "phase_rare")
    config = load_pipeline_config(REPOSITORY_ROOT / "config" / "pipeline.example.yaml")
    config["tools"]["phasing_adapter"] = adapter_config(common, rare)

    environment = build_environment(config)

    adapter_environment = environment["tools"]["phasing_adapter"]
    assert adapter_environment["phase_common"]["available"] is True
    assert "5.1.1" in adapter_environment["phase_common"]["version"]
    assert adapter_environment["phase_rare"]["available"] is True
    assert "5.1.1" in adapter_environment["phase_rare"]["version"]


def test_common_command_requires_reference_map_seed_and_pedigree(tmp_path: Path) -> None:
    probe = probe_fixture(tmp_path)

    command = build_phase_common_command(
        probe,
        input_path=Path("study.bcf"),
        reference_path=Path("reference.vcf.gz"),
        genetic_map_path=Path("chr19.b38.gmap.gz"),
        region="chr19:1000-2000",
        output_path=Path("common.phased.bcf"),
        log_path=Path("common.log"),
        pedigree_path=Path("pedigree.tsv"),
        threads=4,
    )

    assert command[0] == probe.phase_common_path
    assert command[command.index("--reference") + 1] == "reference.vcf.gz"
    assert command[command.index("--map") + 1] == "chr19.b38.gmap.gz"
    assert command[command.index("--seed") + 1] == "15052011"
    assert command[command.index("--thread") + 1] == "4"
    assert command[command.index("--pedigree") + 1] == "pedigree.tsv"
    assert command[command.index("--output-format") + 1] == "bcf"


def test_rare_command_preserves_input_only_variant_path(tmp_path: Path) -> None:
    probe = probe_fixture(tmp_path)

    command = build_phase_rare_command(
        probe,
        input_path=Path("study_with_target.bcf"),
        scaffold_path=Path("common.phased.bcf"),
        genetic_map_path=Path("chr19.b38.gmap.gz"),
        input_region="19:1000-2000",
        scaffold_region="19:900-2100",
        output_path=Path("target.phased.bcf"),
        pedigree_path=Path("pedigree.tsv"),
        threads=2,
        score_singletons=True,
    )

    assert command[0] == probe.phase_rare_path
    assert command[command.index("--input") + 1] == "study_with_target.bcf"
    assert command[command.index("--scaffold") + 1] == "common.phased.bcf"
    assert command[command.index("--effective-size") + 1] == "15000"
    assert command[command.index("--pedigree") + 1] == "pedigree.tsv"
    assert "--log" not in command
    assert "--score-singletons" in command


def test_command_rejects_default_map_substitution_and_non_bcf_output(
    tmp_path: Path,
) -> None:
    probe = probe_fixture(tmp_path)

    with pytest.raises(Shapeit5ContractError, match="invalid_shapeit5_output_format"):
        build_phase_common_command(
            probe,
            input_path=Path("study.bcf"),
            reference_path=Path("reference.bcf"),
            genetic_map_path=Path("chr19.b38.gmap.gz"),
            region="19",
            output_path=Path("phased.vcf.gz"),
            log_path=Path("phase.log"),
        )
