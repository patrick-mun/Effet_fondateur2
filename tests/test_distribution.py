from __future__ import annotations

import email
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def built_wheel(tmp_path_factory: pytest.TempPathFactory) -> Path:
    build_root = tmp_path_factory.mktemp("wheel-source")
    shutil.copy2(REPOSITORY_ROOT / "pyproject.toml", build_root / "pyproject.toml")
    shutil.copytree(REPOSITORY_ROOT / "src", build_root / "src")
    wheel_dir = tmp_path_factory.mktemp("wheel-artifact")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            ".",
            "--no-deps",
            "--no-build-isolation",
            "--wheel-dir",
            str(wheel_dir),
        ],
        cwd=build_root,
        check=True,
        capture_output=True,
        text=True,
    )
    wheels = tuple(wheel_dir.glob("effet_fondateur-*.whl"))
    assert len(wheels) == 1
    return wheels[0]


def test_wheel_declares_runtime_dependencies_and_resources(
    built_wheel: Path,
) -> None:
    with zipfile.ZipFile(built_wheel) as archive:
        names = set(archive.namelist())
        metadata_name = next(
            name for name in names if name.endswith(".dist-info/METADATA")
        )
        metadata = email.message_from_bytes(archive.read(metadata_name))

    requirements = metadata.get_all("Requires-Dist", [])
    assert any(requirement.startswith("numpy==1.26.4") for requirement in requirements)
    assert any(requirement.startswith("scipy<1.14,>=1.11") for requirement in requirements)
    assert any(requirement.startswith("fpdf2>=2.7") for requirement in requirements)
    expected_resources = {
        path.relative_to(REPOSITORY_ROOT / "src").as_posix()
        for resource_dir in (
            REPOSITORY_ROOT / "src" / "effet_fondateur" / "resources" / "schemas",
            REPOSITORY_ROOT / "src" / "effet_fondateur" / "resources" / "references",
            REPOSITORY_ROOT / "src" / "effet_fondateur" / "reporting" / "assets",
        )
        for path in resource_dir.rglob("*")
        if path.is_file()
    }
    assert expected_resources <= names


def test_installed_wheel_uses_only_packaged_resources(
    built_wheel: Path,
    tmp_path: Path,
) -> None:
    install_dir = tmp_path / "installed"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-deps",
            "--target",
            str(install_dir),
            str(built_wheel),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    script = """
from pathlib import Path

import effet_fondateur
from effet_fondateur.contracts import ConfigurationError, load_pipeline_config
from effet_fondateur.references import resolve_reference_panel
from effet_fondateur.resources import REPORTING_ASSETS_DIR, SCHEMAS_DIR

install_dir = Path.cwd() / "installed"
assert Path(effet_fondateur.__file__).resolve().is_relative_to(install_dir.resolve())
assert (SCHEMAS_DIR / "pipeline_config.schema.json").is_file()

invalid_config = Path.cwd() / "invalid.yaml"
invalid_config.write_text("{}\\n", encoding="utf-8")
try:
    load_pipeline_config(invalid_config)
except ConfigurationError as error:
    assert "Configuration invalide" in str(error)
else:
    raise AssertionError("invalid configuration unexpectedly accepted")

panel = resolve_reference_panel(
    panel_id="1kg_3202_high_coverage_20220422",
    assembly="GRCh38",
    chromosome=19,
)
assert panel.sample_count == 3202
assert "@page" in (REPORTING_ASSETS_DIR / "report-linear.css").read_text(
    encoding="utf-8"
)
assert (REPORTING_ASSETS_DIR / "vendor" / "paged.polyfill.min.js").is_file()
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(install_dir)
    subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
