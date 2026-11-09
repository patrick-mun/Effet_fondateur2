"""Visuels D23b : uniquement les agrégats fournis."""
from effet_fondateur.calibration.report_phasing_sensitivity import _contrast_plot, _series_plot


def test_sensitivity_graphics_use_supplied_aggregates_and_are_accessible() -> None:
    rows = [
        {"profile": f"DENSITY_{density}", "estimate": value, "interval_95": [value-.01, value+.01]}
        for density, value in ((25, .2), (50, .15), (75, .1), (100, .05))
    ]
    svg = _series_plot([("Avec pedigree", "circle", rows)], "SER", "SER par densité")
    assert 'role="img"' in svg and "20.00%" in svg and "5.00%" in svg
    contrasts = [dict(contrast_type="WITHOUT_MINUS_WITH_PEDIGREE", metric="switch_error_rate",
                      profile=row["profile"], estimate=.01, interval_95=[-.01, .03]) for row in rows]
    contrast_svg = _contrast_plot(contrasts)
    assert 'role="img"' in contrast_svg and "+1.00 pt" in contrast_svg
