"""Contrôle du rendu D20 à partir d'agrégats fictifs, sans outil externe."""

from pathlib import Path

import pytest

from effet_fondateur.calibration.materialize_ibd import _json
from effet_fondateur.calibration.qc_bench import METHOD_ID, file_inventory, write_json
from effet_fondateur.calibration.report_qc import build_reports


def test_three_independent_aggregate_reports_without_individual_identifiers(tmp_path):
    repository = Path(__file__).resolve().parents[2]
    # Le modèle minimal n'essaie pas de simuler la pagination du navigateur :
    # Chrome et le PDF réel sont contrôlés séparément en intégration annoncée.
    template = tmp_path / "template.html"
    template.write_text('<!DOCTYPE html><html><head><title>Template</title><style>body{color:black}</style></head><body></body></html>')
    workspace = tmp_path / "calibration_workspace"
    source = workspace / "evaluation"
    source.mkdir(parents=True)
    config = _json(repository / "calibration/configurations/d20_qc_bench_v1.json")
    for stage in ("05", "10"):
        key = "sample_missingness" if stage == "05" else "controls_unrelated/sample_missingness"
        write_json(source / f"CAL-REP-{stage}.json", dict(evaluations=[dict(reasons=[]) for _ in range(14)], missing_classes=[],
                   rules={key: dict(unit="sample", tp=1, fn=0, fp=0, tn=10)}))
    manifest = dict(schema_version="1.0.0", method_id=METHOD_ID, evaluation_state="CLOSED", created_at="2026-10-30T00:00:00Z",
                    panel_manifest_sha256="a" * 64, run_index_sha256="b" * 64, configuration=config, plink_versions=["PLINK simulated report test"],
                    stages=[dict(stage_id=stage, status="PASS", case_count=14, fail_count=0) for stage in ("05", "10")],
                    artifacts=file_inventory(source), calibrated=False, founder_effect_proven=False)
    write_json(source / "evaluation_manifest.json", manifest)
    output = workspace / "reports"
    reports = build_reports(run_root=source, template_path=template, output_root=output)
    assert {path.parent.name for path in reports} == {"CAL-DEV-D20", "CAL-REP-05", "CAL-REP-10"}
    for path in reports:
        html = path.read_text()
        assert html.count('<section id="s') == 14
        assert html.count('class="toc-entry"') == 14
        assert html.count('<footer class="page-footer"') == 1
        assert html.count("<figure>") == 2 and "<table>" in html
        assert "calibrated=false" in html and "REVIEW_REQUIRED" in html
        assert "[To complete:" not in html and "SYNTHETIC_INDIVIDUAL_KEY" not in html
    with pytest.raises(ValueError, match="run_exists"):
        build_reports(run_root=source, template_path=template, output_root=output)
