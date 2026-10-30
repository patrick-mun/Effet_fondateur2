"""Rendu D23a : agrégats clos uniquement, sans identifiant ni relance d'outil."""

import json
from pathlib import Path
import re

import pytest

from effet_fondateur.calibration import report_phasing_preparation as report
from effet_fondateur.calibration.phasing_preparation import load_configuration


def test_report_is_independent_illustrated_mathematical_and_aggregate_only(tmp_path: Path, monkeypatch) -> None:
    workspace = tmp_path / "calibration_workspace"
    workspace.mkdir()
    panel = workspace / "panel"
    panel.mkdir()
    manifest_path = panel / "preparation_manifest.json"
    manifest_path.write_text("closed synthetic fixture")
    source = manifest_path.read_bytes()
    config_path = Path(__file__).resolve().parents[2] / "calibration/configurations/d23a_phasing_panel_v2.json"
    manifest = dict(schema_version="1.0.0", method_id="chip_autosomal_phasing_preparation_v1", created_at="2026-10-30T12:00:00Z",
                    configuration=load_configuration(config_path), counts=dict(samples=40, markers=22529, founders=32,
                    trios=8, independent_simulation_units=24, source_donor_groups=64, reference_samples=1403,
                    calibration_donors=900, validation_donors=899), audit=dict(compared_genotype_cells=901160, genotype_mismatch_count=0),
                    tool_versions=["mock-bcftools"], partition_manifest_sha256="a" * 64,
                    reference_selection_sha256="b" * 64, source_fingerprints=[], chip_package_sha256="c" * 64,
                    chip_marker_table_sha256="d" * 64, map_archive_sha256="e" * 64, phasing_executed=False,
                    performance_status="NOT_MEASURED_PREPARATION_ONLY", calibrated=False, founder_effect_proven=False,
                    chromosomes=[dict(chromosome=chromosome, retained_chip_marker_count=1024, rare_heterozygote_count=10)
                                 for chromosome in range(1,23)], artifacts=[dict(path="SAMPLE_IDENTIFIER_MUST_NOT_APPEAR")])
    monkeypatch.setattr(report, "verify_preparation", lambda root: manifest)
    template = tmp_path / "template.html"
    template.write_text('<html><head><title>Template</title></head><body></body></html>')
    output = workspace / "report"
    html_path = report.build_report(panel_root=panel, template_path=template, output_root=output)
    html = html_path.read_text()
    assert len(re.findall(r'<section id="s\d+">', html)) == 14
    assert html.count('<figure class="scientific-figure">') == 2
    assert html.count('role="math"') == 2
    assert all(f"Figure {number} —" in html and f"Équation {number} —" in html for number in (1,2))
    assert "aucune précision de phasage" in html and "FAIL 07 toujours ouvert" in html
    data = json.loads((html_path.parent / "report_data.json").read_text())
    assert "artifacts" not in data and "SAMPLE_IDENTIFIER_MUST_NOT_APPEAR" not in html + json.dumps(data)
    assert not data["phasing_executed"] and not data["calibrated"]
    assert manifest_path.read_bytes() == source
    with pytest.raises(ValueError, match="run_exists"):
        report.build_report(panel_root=panel, template_path=template, output_root=output)
