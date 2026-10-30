"""D21 : rapports séparés, verdicts non masqués et agrégats sans identifiants."""

from pathlib import Path

import pytest

from effet_fondateur.calibration.kinship_bench import METHOD_ID
from effet_fondateur.calibration.kinship_oracle import classification_metrics
from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.report_kinship import build_reports
from effet_fondateur.calibration.simulation.kinship_bench import DEGREES


def test_independent_reports_keep_scientific_fail_visible(tmp_path):
    repository = Path(__file__).resolve().parents[2]
    config = _json(repository / "calibration/configurations/d21_kinship_bench_v1.json")
    template = tmp_path / "template.html"
    template.write_text('<!DOCTYPE html><html><head><title>Test</title></head><body></body></html>')
    source = tmp_path / "calibration_workspace/evaluation"
    source.mkdir(parents=True)
    records = [dict(true_degree=degree, observed_degree=degree, stratum_id="ARTIFICIAL_A") for degree in DEGREES]
    records[4]["observed_degree"] = "UNRELATED"
    performance = classification_metrics(records, config)
    panel_rows, kinship_rows = [], []
    for index, profile in enumerate(config["profiles"]):
        common = dict(bundle_id=f"B{index:02d}", profile=profile, context=dict(post_qc_markers=2000))
        panel_rows.append(dict(**common, status="PASS", reasons=[], metrics=dict(retained_markers=1000, autosome_coverage=22, residual_ld_pair_count=3000, residual_ld_max_r2=.1)))
        kinship_rows.append(dict(**common, status="FAIL", functional_concordance="PASS", reasons=["true_related_pair_in_independent_proposal"],
                                 metrics=dict(forbidden_true_pair_count=1), performance=performance, pre_pruning_oracle_performance=performance))
    for stage, rows in (("06", panel_rows), ("07", kinship_rows)):
        write_json(source / f"CAL-REP-{stage}.json", dict(evaluations=rows))
    manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, evaluation_state="CLOSED", created_at="2026-10-30T00:00:00Z",
                    panel_manifest_sha256="a" * 64, run_index_sha256="b" * 64, configuration=config, tool_versions=["KING simulated report test"],
                    stages=[dict(stage_id="06", status="PASS", case_count=3, fail_count=0), dict(stage_id="07", status="FAIL", case_count=3, fail_count=3)],
                    artifacts=file_inventory(source), calibrated=False, founder_effect_proven=False)
    write_json(source / "evaluation_manifest.json", manifest)
    output = source.parent / "reports"
    reports = build_reports(run_root=source, template_path=template, output_root=output)
    assert {path.parent.name for path in reports} == {"CAL-DEV-D21", "CAL-REP-06", "CAL-REP-07"}
    for path in reports:
        html = path.read_text()
        assert html.count('<section id="s') == 14 and html.count('class="toc-entry"') == 14
        assert html.count('<figure>') == 2 and '<table>' in html
        assert "calibrated=false" in html and "kinship_exclusion_approval" in html
        assert "[To complete:" not in html and "SYNTHETIC_INDIVIDUAL_KEY" not in html
    html = (output / "CAL-REP-07/index.html").read_text()
    assert "Verdict mesuré : FAIL" in html and "Confusion complète" in html
    assert "Wilson 95 %" in html and "conditionnel au panel" in html
    with pytest.raises(ValueError, match="run_exists"):
        build_reports(run_root=source, template_path=template, output_root=output)
