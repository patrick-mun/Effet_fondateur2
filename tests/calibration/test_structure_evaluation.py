"""D22 : agrégateur 08 et vrai composant 09 sans aucun outil externe."""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from effet_fondateur.calibration.evaluate_structure import _invariance, _population_expectation, evaluate_freeze, evaluate_population
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_bench import write_json, write_table
from effet_fondateur.calibration.run_structure_bench import component_inputs, fixture_reviews
from effet_fondateur.calibration.structure_bench import CASES, COHORTS, PRODUCER, fixture_tables, generate_cases, registry_rows, upstream_configuration
from effet_fondateur.calibration.technical_evaluation import observation_matrix
from effet_fondateur.contracts import build_file_artifact
from effet_fondateur.stages.freeze_cohorts import CohortFreezeBlockError, CohortFreezeInputError, execute


REPOSITORY = Path(__file__).resolve().parents[2]


@pytest.fixture
def context(tmp_path):
    config = _json(REPOSITORY / "calibration/configurations/d22_structure_bench_v1.json")
    config["markers_per_autosome"] = 32
    target = _json(REPOSITORY / "calibration/configurations/d19_technical_target_v1.json")
    panels, truth = generate_cases(config)
    panel = panels[0]
    registry = registry_rows(panel, CASES[0])
    tables = fixture_tables(panel, registry, config, target["target_metadata"])
    matrix = observation_matrix(panel, registry)
    root = tmp_path / "run"
    fixture = root / "stages/D22_fixture"
    fixture.mkdir(parents=True)
    write_json(root / "manifest.json", dict(run_id="d22_test"))
    (root / "config.resolved.yaml").write_text(yaml.safe_dump(upstream_configuration(config, target)))
    bim = fixture / "panel.bim"
    bim.write_text("".join(f"{chromosome}\t{variant}\t0\t{position}\tG\tA\n" for chromosome, position, variant in matrix.variants))
    model, expected = _population_expectation(matrix, bim, registry, tables["independent_set_proposal"], "fixture_set", config)
    paths, artifacts = {}, []

    def publish(name, rows, path):
        converted = [{key: "" if isinstance(value, float) and not np.isfinite(value) else value for key, value in row.items()} for row in rows]
        columns = tuple(_json(REPOSITORY / "schemas" / (name + ".schema.json"))["x-tsv"]["columns"])
        write_table(path, columns, converted)
        paths[name] = path
        artifacts.append(build_file_artifact(path, path.relative_to(root).as_posix(), name, name, "text/tab-separated-values", PRODUCER, "a" * 64, assembly="GRCh38"))

    for name, rows in {"samples_master": registry, **tables, **expected}.items():
        publish(name, rows, fixture / (name + ".tsv"))
    metadata_path = fixture / "target_metadata.yaml"
    metadata_path.write_text(yaml.safe_dump(target["target_metadata"]))
    artifacts.append(build_file_artifact(metadata_path, metadata_path.relative_to(root).as_posix(), "config_input_target_variant_metadata", "target_metadata", "text/yaml", PRODUCER, "a" * 64, assembly="GRCh38"))
    paths["population_structure_report"] = fixture / "population_structure_report.json"
    write_json(paths["population_structure_report"], dict(panel_sample_count=128, reference_sample_count=108, projected_sample_count=20, panel_variant_count=704,
               informative_variant_count=model["informative_count"], component_count=3, outlier_component_count=2, outlier_count=int(model["outliers"].sum()), manual_validation_required=bool(model["outliers"].any())))
    population = {item["artifact_id"]: item for item in artifacts}
    reviews = fixture_reviews(root, artifacts, population, config["fixture_review_policy"])
    parameters = dict(control_group_labels=["CONTROL"], required_nonempty_cohorts=list(COHORTS), decision_reviews=reviews)
    freeze_artifacts = [item for item in artifacts if item["artifact_id"] in {"samples_master", "target_genotype_audit", "qc_individual_metrics", "independent_set_proposal", "population_scores", "population_outliers", "config_input_target_variant_metadata"}]
    return dict(config=config, root=root, matrix=matrix, paths=paths, registry=registry, tables=tables, descriptor=dict(sample_set_id="fixture_set"),
                bim=bim, truth=truth, parameters=parameters, artifacts=freeze_artifacts)


def test_population_aggregate_has_no_duplicate_metric_error(context):
    result, published = evaluate_population(context["matrix"], context["paths"], context["registry"], context["tables"]["independent_set_proposal"],
                                           context["descriptor"], context["bim"], context["truth"], context["config"])
    assert result["status"] == "PASS" and result["reasons"] == []
    assert result["metrics"]["reference_sample_count"] == 108
    assert published["scores"].shape == (128, 3)


@pytest.mark.parametrize("artifact,field,value", [("population_scores", "PC1", "123"), ("population_outliers", "PROPOSED_EXCLUSION", "true"), ("population_variant_loadings", "MODEL_STATUS", "INCLUDED")])
def test_population_oracle_detects_mutated_outputs(context, artifact, field, value):
    from effet_fondateur.calibration.evaluate_technical import raw_table
    path = context["paths"][artifact]
    rows = raw_table(path)
    rows[0][field] = value
    write_table(path, tuple(rows[0]), rows)
    result, _ = evaluate_population(context["matrix"], context["paths"], context["registry"], context["tables"]["independent_set_proposal"],
                                    context["descriptor"], context["bim"], context["truth"], context["config"])
    assert result["status"] == "FAIL" and result["reasons"]


def test_real_freeze_component_matches_independent_oracle_and_repeats(context):
    directories = []
    for label in ("a", "b"):
        path = component_inputs(context["root"], "09", "freeze_cohorts", "09_test_" + label, context["parameters"], context["artifacts"])
        assert execute(path / "stage_inputs.json", path) == 0
        directories.append(path.relative_to(context["root"]).as_posix())
    result = evaluate_freeze(context["root"], directories, context["paths"], context["registry"], context["tables"], "A", context["config"])
    assert result["status"] == "PASS", result["reasons"]
    assert result["metrics"]["cohort_repeat_mismatch_count"] == 0
    assert result["metrics"]["membership_count"] == 896
    # Une sortie altérée ne devient pas un faux PASS de reproductibilité.
    path = context["root"] / directories[1] / "keep/controls_unrelated.keep"
    path.write_text("")
    with pytest.raises(ValueError, match="artifact_invalid"):
        evaluate_freeze(context["root"], directories, context["paths"], context["registry"], context["tables"], "A", context["config"])


@pytest.mark.parametrize("mutation", ["missing", "sha", "incomplete", "timezone"])
def test_real_freeze_refuses_invalid_fixture_reviews(context, mutation):
    parameters = deepcopy(context["parameters"])
    review = parameters["decision_reviews"]["kinship_exclusion_approval"]
    expected = CohortFreezeBlockError
    if mutation == "missing":
        parameters["decision_reviews"] = {}
    elif mutation == "sha":
        review["reviewed_artifact_sha256"] = "0" * 64
    elif mutation == "incomplete":
        review["approved_exclusion_sample_ids"] = []
    else:
        review["reviewed_at"] = "2026-10-30T07:00:00"
        expected = CohortFreezeInputError
    path = component_inputs(context["root"], "09", "freeze_cohorts", "09_negative", parameters, context["artifacts"])
    with pytest.raises(expected):
        execute(path / "stage_inputs.json", path)
    assert not (path / "stage_outputs.json").exists() and not (path / "cohorts.frozen.tsv").exists()


def test_axis_sign_is_not_a_structural_change_but_score_distortion_is():
    config = _json(REPOSITORY / "calibration/configurations/d22_structure_bench_v1.json")
    first = dict(reference=np.ones(3, bool), scores=np.array([[1., 2], [0., 0], [-1., -2]]), eigenvalues=np.array([2., 1.]))
    second = deepcopy(first)
    second["scores"][:, 0] *= -1
    assert _invariance(first, second, config)["status"] == "PASS"
    second["scores"][0, 0] += .1
    assert _invariance(first, second, config)["status"] == "FAIL"
