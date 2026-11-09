"""D22 : rapports indépendants, verdicts conservés et numérotation continue."""

from pathlib import Path
import re

import pytest

from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.report_structure import build_reports
from effet_fondateur.calibration.structure_bench import CASES, COHORTS, METHOD_ID


def test_three_reports_preserve_limits_and_fail_without_individual_identifiers(tmp_path, monkeypatch):
    repository = Path(__file__).resolve().parents[2]
    config = _json(repository / "calibration/configurations/d22_structure_bench_v1.json")
    source = tmp_path / "calibration_workspace/evaluation"
    source.mkdir(parents=True)
    population, cohorts = [], []
    for index, case in enumerate(CASES):
        detection = dict(tp=4, fn=0, fp=0, tn=104, sensitivity_wilson=[.51, 1.], non_injected_rate_wilson=[0., .036])
        population.append(dict(case=case, status="FAIL" if index == 0 else "PASS", reasons=["injected_test_discordance"] if index == 0 else [],
                               metrics=dict(informative_variant_count=1406, centroid_separation_over_within_sd=40., mixed_coordinate_mean=.5, independent_reference_detection=detection)))
        cohorts.append(dict(case=case, status="PASS", reasons=[], metrics=dict(included_counts={cohort: 10 for cohort in COHORTS}), review_rejections=[]))
    for stage, rows in (("08", population), ("09", cohorts)):
        write_json(source / ("CAL-REP-" + stage + ".json"), dict(evaluations=rows))
    manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, evaluation_state="CLOSED", created_at="2026-10-30T07:00:00Z",
                    configuration=config, tool_versions=["PLINK report fixture"], panel_manifest_sha256="a" * 64, run_index_sha256="b" * 64,
                    scope="ISOLATED_COMPONENTS_08_09_NOT_04_07", kinship_07_issue="OPEN_DEFERRED_UNTIL_GLOBAL_REVIEW_BY_USER",
                    stages=[dict(stage_id="08", status="FAIL", case_count=4, fail_count=1), dict(stage_id="09", status="PASS", case_count=4, fail_count=0)],
                    artifacts=file_inventory(source), calibrated=False, founder_effect_proven=False)
    write_json(source / "evaluation_manifest.json", manifest)
    panel_root, pipeline_root = source.parent / "panel", source.parent / "pipeline"
    for root in (panel_root, pipeline_root):
        root.mkdir()
        write_json(root / "placeholder.json", {})
    shared = dict(schema_version="1.0.0", method_id=METHOD_ID, calibrated=False, founder_effect_proven=False)
    write_json(panel_root / "panel_manifest.json", dict(**shared, manifest_state="FINALIZED", created_at=manifest["created_at"],
               configuration_sha256="c" * 64, bundles=[dict(bundle_id=f"B{index:02d}", sample_count=128, marker_count=1408) for index in range(4)], artifacts=file_inventory(panel_root)))
    panel_sha = _sha(panel_root / "panel_manifest.json")
    write_json(pipeline_root / "run_index.json", dict(**shared, closure_state="CLOSED", panel_manifest_sha256=panel_sha,
               runs=[dict(bundle_id=f"B{index:02d}", run_path=f"fixture/B{index:02d}", manifest_sha256="d" * 64, upstream_returncode=0, population_returncode=0,
                          freeze_directories=["stages/a", "stages/b"], review_rejections=[], scope=manifest["scope"], fixture_producer="D22_UPSTREAM_CONTRACT_FIXTURE") for index in range(4)], artifacts=file_inventory(pipeline_root)))
    centroids = [dict(group=group, role="REFERENCE", count=count, mean=mean, standard_deviation=[.02, .03]) for group, count, mean in
                 (("A", 40, [-.7, 0]), ("B", 40, [.7, 0]), ("MIXED", 24, [0, 0]), ("RARE", 4, [0, 1.7]))]
    visual_case = dict(centroids=centroids, sample_count=128, eigenvalues=[dict(component=f"PC{index+1}", explained_variance_ratio=value) for index, value in enumerate((.42, .025, .012))],
                       distance_histogram=dict(edges=[0., 10., 20., 30.], reference_counts=[104, 0, 4], projected_counts=[16, 0, 4], threshold=9.21034),
                       cohorts={cohort: dict(included=10, not_included=118, exclusion_counts={"not_explicit_control_group": 118}) for cohort in COHORTS})
    visuals = dict(cases=[dict(case=case, **visual_case) for case in CASES], aggregation_only=True, individual_identifiers_exported=False,
                   source_fingerprints=dict(panel_manifest_sha256=panel_sha, run_index_sha256=_sha(pipeline_root / "run_index.json"), evaluation_manifest_sha256=_sha(source / "evaluation_manifest.json")))
    monkeypatch.setattr("effet_fondateur.calibration.report_structure.collect_visual_evidence", lambda **kwargs: visuals)
    template = tmp_path / "template.html"
    template.write_text('<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>')
    output = source.parent / "reports"
    original_manifest = (source / "evaluation_manifest.json").read_bytes()
    paths = build_reports(run_root=source, template_path=template, output_root=output, panel_root=panel_root, pipeline_root=pipeline_root)
    assert {path.parent.name for path in paths} == {"CAL-DEV-D22", "CAL-REP-08", "CAL-REP-09"}
    for path in paths:
        html = path.read_text()
        assert html.count('<section id="s') == 14 and html.count('class="toc-entry"') == 14
        expected_figures = 2 if path.parent.name == "CAL-REP-09" else 4
        assert html.count('class="scientific-figure"') == expected_figures
        assert [int(number) for number in re.findall(r'data-figure="(\d+)"', html)] == list(range(1, expected_figures + 1))
        equations = [int(number) for number in re.findall(r'data-equation="(\d+)"', html)]
        assert len(equations) >= 2 and equations == list(range(1, len(equations) + 1))
        assert 'role="math"' in html and "Équation 1" in html
        captions = [int(number) for number in re.findall(r'<caption>Tableau (\d+)', html)]
        assert captions == list(range(1, len(captions) + 1))
        assert "calibrated=false" in html and "FAIL 07" in html and "calibration_contract_fixture" in html
        assert "[To complete:" not in html and "SYNTHETIC_INDIVIDUAL_KEY" not in html
        assert _json(path.parent / "visual_data.json")["individual_identifiers_exported"] is False
    assert (source / "evaluation_manifest.json").read_bytes() == original_manifest
    assert "Verdict : FAIL" in (output / "CAL-REP-08/index.html").read_text()
    with pytest.raises(ValueError, match="run_exists"):
        build_reports(run_root=source, template_path=template, output_root=output, panel_root=panel_root, pipeline_root=pipeline_root)
