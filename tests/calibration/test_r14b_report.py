from __future__ import annotations

import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.r14b_report import (
    REPORT_ID,
    REPORT_TITLE,
    R14bReportError,
    _quality_control,
    _render_html,
    summarize_ibd_measurements,
    validate_r14b_report,
    wilson_two_sided_95,
)


def test_wilson_intervals_show_small_sample_uncertainty() -> None:
    hap = wilson_two_sided_95(23, 25)
    refined = wilson_two_sided_95(22, 25)

    assert hap == pytest.approx((0.7503389105, 0.9777795987))
    assert refined == pytest.approx((0.7004420608, 0.9583318285))
    with pytest.raises(R14bReportError, match="binomial_counts_invalid"):
        wilson_two_sided_95(26, 25)


def test_ibd_summary_requires_three_replicates_in_six_cells() -> None:
    measurements = []
    for density in ("D_MAX", "D_MID", "D_LOW"):
        for tool in ("HAP_IBD", "REFINED_IBD"):
            for replicate in range(1, 4):
                measurements.append(
                    {
                        "density": density,
                        "tool": tool,
                        "wall_seconds": float(replicate),
                        "maximum_rss_bytes": replicate * 2**20,
                    }
                )

    summaries = summarize_ibd_measurements(measurements)

    assert len(summaries) == 6
    assert summaries[0]["n"] == 3
    assert summaries[0]["wall_mean_seconds"] == 2.0
    assert summaries[0]["rss_mean_mib"] == 2.0


def test_render_replaces_template_body_and_title() -> None:
    rendered = _render_html(
        "<html><head><title>modèle</title></head><body>ancien</body></html>",
        '<body lang="fr"><main>R14b</main></body>',
    )

    assert f"<title>{REPORT_TITLE}</title>" in rendered
    assert "ancien" not in rendered
    assert "R14b" in rendered


def test_quality_control_counts_numbering_and_links() -> None:
    document = """
    <section id="s1"></section><a href="#s1">lien</a>
    <figure><figcaption>Figure 1 : A</figcaption></figure>
    <figure><figcaption>Figure 2 : B</figcaption></figure>
    <figure><figcaption>Figure 3 : C</figcaption></figure>
    <figure><figcaption>Figure 4 : D</figcaption></figure>
    <table><caption>Table 1 : A</caption></table>
    <table><caption>Table 2 : B</caption></table>
    <table><caption>Table 3 : C</caption></table>
    <table><caption>Table 4 : D</caption></table>
    <table><caption>Table 5 : E</caption></table>
    <table><caption>Table 6 : F</caption></table>
    <table><caption>Table 7 : G</caption></table>
    <table><caption>Table 8 : H</caption></table>
    <table><caption>Table 9 : I</caption></table>
    <table><caption>Table 10 : J</caption></table>
    <table><caption>Table 11 : K</caption></table>
    <footer class="page-footer"></footer>
    Paged.js v0.4.3
    """

    qc = _quality_control(document, 15)

    assert qc["figure_numbering_continuous"] is True
    assert qc["table_numbering_continuous"] is True
    assert qc["missing_internal_links"] == []


def test_quality_control_accepts_an_explicit_real_data_exclusion() -> None:
    qc = _quality_control("Aucun fichier sous data/input/ n'a été lu.", 12)

    assert qc["placeholder_count"] == 0


def test_validation_refuses_granted_scientific_qualification(tmp_path: Path) -> None:
    output_root = tmp_path / "r14b"
    output_root.mkdir()
    (output_root / "R14B_REPORT_MANIFEST.json").write_text(
        json.dumps({"report_id": REPORT_ID}), encoding="utf-8"
    )
    (output_root / "R14B_DECISION.json").write_text(
        json.dumps(
            {
                "decision_id": REPORT_ID,
                "shapeit5_integration": "VERIFIED",
                "scientific_qualification": "GRANTED",
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(R14bReportError, match="scientific_qualification_invalid"):
        validate_r14b_report(output_root)
