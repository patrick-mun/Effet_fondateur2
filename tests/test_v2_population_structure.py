import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from effet_fondateur.orchestrator import StageExecutionError
from effet_fondateur.orchestrator.pipeline import run_pipeline
from effet_fondateur.stages.analyze_population_structure import _fit_and_project
from test_v2_infer_kinship import prepare_kinship_inputs


def write_population_plink(
    path: Path,
    delegated_plink: Path,
    *,
    fail_export: bool = False,
    malformed_header: bool = False,
) -> None:
    script = f'''#!{sys.executable}
import os
import pathlib
import sys

if "--recode" not in sys.argv:
    os.execv({str(delegated_plink)!r}, [{str(delegated_plink)!r}, *sys.argv[1:]])
if {fail_export!r}:
    raise SystemExit(8)
input_prefix = pathlib.Path(sys.argv[sys.argv.index("--bfile") + 1])
output_prefix = pathlib.Path(sys.argv[sys.argv.index("--out") + 1])
fam_rows = [line.split() for line in input_prefix.with_suffix(".fam").read_text().splitlines() if line]
bim_rows = [line.split() for line in input_prefix.with_suffix(".bim").read_text().splitlines() if line]
variant_columns = [f"{{row[1]}}_{{row[4]}}" for row in bim_rows]
if {malformed_header!r}:
    variant_columns[0] = "wrong_variant_A"
lines = ["FID IID PAT MAT SEX PHENOTYPE " + " ".join(variant_columns) + "\\n"]
for sample_index, fam_row in enumerate(fam_rows):
    dosages = []
    for variant_index, _ in enumerate(bim_rows):
        if sample_index == 0:
            dosage = 0
        elif sample_index == 1:
            dosage = 1 if variant_index % 2 else 0
        else:
            dosage = 2 if variant_index % 2 else 0
        dosages.append(str(dosage))
    lines.append(" ".join([*fam_row[:6], *dosages]) + "\\n")
output_prefix.with_suffix(".raw").write_text("".join(lines))
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def prepare_population_inputs(
    tmp_path: Path,
    *,
    related_pair: bool = False,
    fail_export: bool = False,
    malformed_header: bool = False,
    parameter_overrides: dict[str, object] | None = None,
) -> tuple[Path, Path]:
    between_rows = ["F1 I1 F2 I2 22 0.25 0.10 0.25"] if related_pair else None
    config_path, runs_dir = prepare_kinship_inputs(
        tmp_path,
        between_rows=between_rows,
    )
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    delegated_plink = Path(config["tools"]["plink"])
    population_plink = tmp_path / "population_plink"
    write_population_plink(
        population_plink,
        delegated_plink,
        fail_export=fail_export,
        malformed_header=malformed_header,
    )
    config["tools"]["plink"] = str(population_plink)
    parameters: dict[str, object] = {
        "requested_components": 3,
        "outlier_components": 2,
        "min_reference_samples": 2,
        "min_informative_variants": 10,
        "min_reference_call_rate": 0.95,
        "outlier_alpha": 1e-12,
        "plink_timeout_seconds": 30,
    }
    parameters.update(parameter_overrides or {})
    config["stages"]["analyze_population_structure"] = {
        "enabled": True,
        "parameters": parameters,
    }
    config_path.write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    return config_path, runs_dir


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as input_file:
        return list(csv.DictReader(input_file, delimiter="\t"))


def test_population_pca_uses_independent_reference(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(tmp_path)

    run_dir = run_pipeline(config_path, runs_dir)

    stage_dir = run_dir / "stages" / "08_analyze_population_structure"
    scores = read_tsv(stage_dir / "population_scores.tsv")
    assert len(scores) == 3
    assert all(row["REFERENCE_INCLUDED"] == "true" for row in scores)
    assert all(row["PROJECTED"] == "false" for row in scores)
    assert all(row["PC1"] for row in scores)
    assert all(row["PC3"] == "" for row in scores)
    eigenvalues = read_tsv(stage_dir / "population_eigenvalues.tsv")
    assert [row["COMPONENT"] for row in eigenvalues] == ["PC1", "PC2"]
    report = json.loads(
        (stage_dir / "population_structure_report.json").read_text(encoding="utf-8")
    )
    assert report["reference_sample_count"] == 3
    assert report["projected_sample_count"] == 0
    assert report["informative_variant_count"] == 11
    assert report["reference_sample_set_id"].startswith("population_reference_")
    assert scores[0]["SAMPLE_SET_ID"] == report["population_sample_set_id"]
    assert scores[0]["SAMPLE_SET_ID"] != report["reference_sample_set_id"]
    loadings = read_tsv(stage_dir / "population_variant_loadings.tsv")
    assert sum(row["MODEL_STATUS"] == "INCLUDED" for row in loadings) == 11
    assert {
        row["EXCLUSION_CODE"]
        for row in loadings
        if row["MODEL_STATUS"] == "EXCLUDED"
    } == {"MONOMORPHIC_IN_REFERENCE"}
    stage_outputs = json.loads(
        (stage_dir / "stage_outputs.json").read_text(encoding="utf-8")
    )
    artifacts = {
        artifact["artifact_id"]: artifact for artifact in stage_outputs["artifacts"]
    }
    assert artifacts["population_variant_loadings"]["sensitivity"] == "sensitive_genetic"
    assert not list(stage_dir.glob("*.raw"))


def test_related_sample_is_projected_without_influencing_axes(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(tmp_path, related_pair=True)

    run_dir = run_pipeline(config_path, runs_dir)

    scores = read_tsv(
        run_dir
        / "stages"
        / "08_analyze_population_structure"
        / "population_scores.tsv"
    )
    sample_1 = next(row for row in scores if row["SAMPLE_ID"] == "sample_1")
    assert sample_1["REFERENCE_INCLUDED"] == "false"
    assert sample_1["PROJECTED"] == "true"
    assert sum(row["REFERENCE_INCLUDED"] == "true" for row in scores) == 2


def test_population_outliers_require_manual_review(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(
        tmp_path,
        parameter_overrides={"outlier_alpha": 1.0},
    )

    run_dir = run_pipeline(config_path, runs_dir)

    outliers = read_tsv(
        run_dir
        / "stages"
        / "08_analyze_population_structure"
        / "population_outliers.tsv"
    )
    assert any(row["PROPOSED_EXCLUSION"] == "true" for row in outliers)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert "population_structure_exclusion_approval" in manifest["manual_decisions_required"]


def test_insufficient_reference_samples_is_scientific_block(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(
        tmp_path,
        related_pair=True,
        parameter_overrides={"min_reference_samples": 3},
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 4


def test_population_plink_failure_uses_external_tool_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(tmp_path, fail_export=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3


def test_malformed_dosage_header_uses_external_tool_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(tmp_path, malformed_header=True)

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 3


def test_invalid_component_configuration_uses_input_code(tmp_path: Path) -> None:
    config_path, runs_dir = prepare_population_inputs(
        tmp_path,
        parameter_overrides={"requested_components": 11},
    )

    with pytest.raises(StageExecutionError) as error:
        run_pipeline(config_path, runs_dir)

    assert error.value.return_code == 2


def test_projected_dosages_do_not_change_reference_axes() -> None:
    parameters = {
        "requested_components": 2,
        "outlier_components": 2,
        "min_informative_variants": 2,
        "min_reference_call_rate": 1.0,
    }
    first_dosages = np.array(
        [[0, 0, 0], [0, 1, 2], [2, 1, 0], [1, 1, 1]], dtype=float
    )
    second_dosages = first_dosages.copy()
    second_dosages[3] = [2, 2, 2]
    reference_indices = np.array([0, 1, 2])

    first_model = _fit_and_project(first_dosages, reference_indices, parameters)
    second_model = _fit_and_project(second_dosages, reference_indices, parameters)

    np.testing.assert_allclose(first_model["right_vectors"], second_model["right_vectors"])
    np.testing.assert_allclose(first_model["scores"][:3], second_model["scores"][:3])
