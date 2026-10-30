"""Visuels D23 : mise à l'échelle des seuls agrégats clos, aucun estimateur."""
from __future__ import annotations

from html import escape
from typing import Any


def _text(x: float, y: float, text: str, anchor: str = "start") -> str:
    return f'<text x="{x:g}" y="{y:g}" text-anchor="{anchor}">{escape(text)}</text>'


def _svg(title: str, body: str, height: int) -> str:
    return f'<svg viewBox="0 0 680 {height}" role="img" aria-label="{escape(title, quote=True)}"><title>{escape(title)}</title><style>text{{font:12px Arial,sans-serif;fill:#243641}}</style>{body}</svg>'


def count_bars(labels: list[str], values: list[int | float], title: str, unit: str) -> str:
    """Barres de comptes ou résidus, avec valeurs numériques explicites."""
    maximum = max(max(values, default=0), 1e-15)
    body = _text(8, 18, title) + _text(190, 40, unit)
    for index, (label, value) in enumerate(zip(labels, values, strict=True)):
        y = 60 + 40 * index
        width = 360 * value / maximum
        body += _text(8, y + 14, label) + f'<rect x="190" y="{y}" width="360" height="18" fill="#e7edf0"/><rect x="190" y="{y}" width="{width:g}" height="18" fill="#286d91"/>' + _text(565, y + 14, f"{value:g}")
    return _svg(title, body, 80 + 40 * len(labels))


def rate_intervals(evaluations: list[dict[str, Any]], metric: str, title: str) -> str:
    """Points et IC familiaux : le rendu conserve les estimations et bornes reçues."""
    measures = [row.get("metrics", {}).get(metric, {}).get("ALL") for row in evaluations]
    maxima = [value for measure in measures if measure is not None for value in ([measure["estimate"]] + (measure["interval_95"] or [])) if value is not None]
    maximum = min(max(max(maxima, default=.01) * 1.12, .01), 1)
    body = _text(8, 18, title)
    for fraction in (0, .25, .5, .75, 1):
        x = 230 + 390 * fraction
        body += f'<line x1="{x:g}" y1="37" x2="{x:g}" y2="{55+65*len(evaluations)}" stroke="#dce3e8"/>'
        body += _text(x, 50 + 65 * len(evaluations), f"{100*maximum*fraction:.2f}", "middle")
    for index, (row, measure) in enumerate(zip(evaluations, measures, strict=True)):
        y = 62 + 65 * index
        body += _text(8, y, row["profile"])
        if measure is None or measure["estimate"] is None:
            body += _text(8, y + 18, "Non évaluable : dénominateur nul ou sortie absente")
            continue
        value = measure["estimate"]
        body += _text(8, y + 18, f'{measure["numerator"]}/{measure["denominator"]} ; {100*value:.3f} %')
        if measure["interval_95"] is not None:
            lower, upper = measure["interval_95"]
            start, end = 230 + 390 * lower / maximum, 230 + 390 * upper / maximum
            body += f'<line x1="{start:g}" y1="{y}" x2="{end:g}" y2="{y}" stroke="#286d91" stroke-width="4"/>'
            body += _text(230, y + 20, f'IC 95 % : {100*lower:.3f}–{100*upper:.3f} % ; unités={measure["contributing_units"]}')
        body += f'<circle cx="{230+390*value/maximum:g}" cy="{y}" r="4" fill="#9b3e61"/>'
    body += _text(420, 72 + 65 * len(evaluations), "Taux (%) — IC descriptifs, pas dispersion", "middle")
    return _svg(title, body, 90 + 65 * len(evaluations))


def genetic_map_curve(data: dict[str, Any]) -> str:
    """Carte du cas de base : lignes attendues et observations publiées."""
    positions, expected, observed = data["position_bp"], data["expected_cm"], data["observed_cm"]
    first, last = min(positions), max(positions)
    lower, upper = min(expected + observed), max(expected + observed)
    def x(value: float) -> float:
        return 75 + 550 * (value - first) / max(last - first, 1)
    def y(value: float) -> float:
        return 220 - 160 * (value - lower) / max(upper - lower, 1e-12)
    body = _text(8, 18, "Carte régionale : oracle et sortie 11, cas de base")
    for fraction in (0, .5, 1):
        body += _text(65, 224 - 160 * fraction, f"{lower+(upper-lower)*fraction:.3f}", "end")
        body += _text(75 + 550 * fraction, 240, f"{(first+(last-first)*fraction)/1e6:.2f}", "middle")
    points = " ".join(f"{x(bp):g},{y(cm):g}" for bp, cm in zip(positions, expected, strict=True))
    body += f'<polyline points="{points}" fill="none" stroke="#286d91" stroke-width="2"/>'
    stride = max(len(positions) // 40, 1)
    for index in range(0, len(positions), stride):
        body += f'<circle cx="{x(positions[index]):g}" cy="{y(observed[index]):g}" r="2" fill="#b06a20"/>'
    body += _text(8, 42, "Position génétique (cM)") + _text(350, 262, "Position physique (Mb, GRCh38)", "middle")
    body += _text(70, 286, "Bleu : oracle · ocre : sortie (points espacés pour lisibilité)")
    return _svg("Concordance des cartes en bp et cM", body, 300)
