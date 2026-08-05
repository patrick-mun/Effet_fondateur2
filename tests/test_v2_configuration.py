from pathlib import Path

import pytest

from effet_fondateur.cli import main
from effet_fondateur.contracts import ConfigurationError, load_pipeline_config


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_generic_example_is_valid() -> None:
    config = load_pipeline_config(REPOSITORY_ROOT / "config" / "pipeline.example.yaml")

    assert config["schema_version"] == "1.0.0"
    assert config["target"]["project_variant_id"] == "target_GRCh38_1_100000_A_G"
    assert config["target"]["rsid"] is None


def test_dock6_profile_accepts_unconfirmed_target_fields() -> None:
    config = load_pipeline_config(
        REPOSITORY_ROOT / "config" / "studies" / "dock6.example.yaml"
    )

    assert config["target"]["gene"] == "DOCK6"
    assert config["target"]["chromosome"] == 19
    assert config["target"]["position_bp"] is None


def test_synthetic_orchestrator_profile_is_valid() -> None:
    config = load_pipeline_config(
        REPOSITORY_ROOT / "config" / "testing" / "synthetic.example.yaml"
    )

    assert config["stages"]["synthetic_stage"]["enabled"] is True
    assert all(command is None for command in config["tools"].values())


def test_unknown_configuration_key_is_rejected(tmp_path: Path) -> None:
    config_path = tmp_path / "invalid.yaml"
    config_path.write_text(
        (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(encoding="utf-8")
        + "\nunexpected: true\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="Additional properties"):
        load_pipeline_config(config_path)


def test_forbidden_clinical_inference_is_rejected(tmp_path: Path) -> None:
    example = (REPOSITORY_ROOT / "config" / "pipeline.example.yaml").read_text(
        encoding="utf-8"
    )
    config_path = tmp_path / "unsafe.yaml"
    config_path.write_text(
        example.replace("allow_clinical_genotype_inference: false", "allow_clinical_genotype_inference: true"),
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="allow_clinical_genotype_inference"):
        load_pipeline_config(config_path)


def test_cli_validates_configuration(capsys: pytest.CaptureFixture[str]) -> None:
    exit_code = main(
        [
            "validate-config",
            str(REPOSITORY_ROOT / "config" / "pipeline.example.yaml"),
        ]
    )

    assert exit_code == 0
    assert "Configuration valide" in capsys.readouterr().out
