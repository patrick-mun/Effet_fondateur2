from __future__ import annotations

import json
from pathlib import Path

import pytest

from effet_fondateur.calibration.r14_report import (
    REPORT_ID,
    REPORT_TITLE,
    R14ReportError,
    _render_html,
    _source_table_rows,
    validate_r14_report,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_r14_render_replaces_title_and_body() -> None:
    template = "<html><head><title>Modèle</title></head><body>ancien</body></html>"

    rendered = _render_html(template, '<body lang="fr"><main>nouveau</main></body>')

    assert f"<title>{REPORT_TITLE}</title>" in rendered
    assert '<body lang="fr"><main>nouveau</main></body>' in rendered
    assert "ancien" not in rendered


def test_r14_source_table_uses_short_paths_and_truncated_hashes() -> None:
    rows = _source_table_rows(
        [
            {
                "relative_path": "calibration_workspace/zone/intermediaire/source.json",
                "role": "preuve technique",
                "sha256": "a" * 64,
            }
        ]
    )

    assert "intermediaire/source.json" in rows
    assert "calibration_workspace/zone" not in rows
    assert "aaaaaaaaaaaa…" in rows


def test_r14_validation_refuses_scientific_qualification(tmp_path: Path) -> None:
    output_root = tmp_path / "r14"
    output_root.mkdir()
    (output_root / "R14_REPORT_MANIFEST.json").write_text(
        json.dumps({"report_id": REPORT_ID}),
        encoding="utf-8",
    )
    (output_root / "R14_DECISION.json").write_text(
        json.dumps(
            {
                "decision_id": REPORT_ID,
                "scientific_qualification": "GRANTED",
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(R14ReportError, match="r14_scientific_decision_mismatch"):
        validate_r14_report(REPOSITORY_ROOT, output_root)
