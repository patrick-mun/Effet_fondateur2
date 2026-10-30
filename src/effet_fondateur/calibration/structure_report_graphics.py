"""Figures SVG D22 : uniquement mise à l'échelle des agrégats empreintés."""

from __future__ import annotations

from html import escape
from typing import Any


COLORS = dict(A="#286d91", B="#b06a20", MIXED="#287c58", RARE="#9b3e61")
GROUP_LABELS = dict(A="A", B="B", MIXED="Mélanges", RARE="Atypiques")
COHORT_LABELS = ("Témoins indépendants", "Tous porteurs", "Porteurs indépendants", "Non-porteurs familiaux", "Atteints non-porteurs", "Chromosome cible", "Référence structure")
REASON_LABELS = dict(not_explicit_control_group="Rôle non témoin", target_variant_carrier="Porteur explicite",
                     approved_kinship_exclusion="Exclusion 07 de test", approved_population_exclusion="Atypique 08 de test",
                     same_family_control_represented="Représentant FID déjà retenu", preliminary_qc_excluded="Absent du panel",
                     genomewide_not_included="Génome entier non autorisé")


def _text(x: float, y: float, label: str, **attributes: str) -> str:
    """Texte toujours échappé ; aucun fragment fourni par les données n'est HTML."""
    extra = " ".join(f'{key.replace("_", "-")}="{escape(value, quote=True)}"' for key, value in attributes.items())
    return f'<text x="{x:.2f}" y="{y:.2f}" {extra}>{escape(label)}</text>'


def _svg(title: str, height: int, body: str) -> str:
    return f'<svg viewBox="0 0 680 {height}" role="img" aria-label="{escape(title, quote=True)}"><title>{escape(title)}</title><style>text{{font:12px Arial,sans-serif;fill:#243641}}line{{stroke:#cad5db}}</style>' + body + '</svg>'


def figure(number: int, svg: str, caption: str) -> str:
    """Enveloppe numérotée et légende source, indissociables à l'impression."""
    return f'<figure class="scientific-figure" data-figure="{number}">{svg}<figcaption>Figure {number} — {escape(caption)}</figcaption></figure>'


def pca_centroids(case: dict[str, Any]) -> str:
    """Graphique agrégé : moyenne ± un écart-type, jamais intervalle de confiance."""
    centers = case["centroids"]
    extrema = [[row["mean"][axis] + direction * row["standard_deviation"][axis] for row in centers for direction in (-1, 1)] for axis in (0, 1)]
    limits = []
    for values in extrema:
        span = max(max(values) - min(values), .1)
        limits.append((min(values) - .15 * span, max(values) + .15 * span))

    def x(value: float) -> float:
        return 65 + 550 * (value - limits[0][0]) / (limits[0][1] - limits[0][0])

    def y(value: float) -> float:
        return 240 - 195 * (value - limits[1][0]) / (limits[1][1] - limits[1][0])

    body = _text(8, 18, "PCA agrégée — cas de base, axes sans unité")
    for index in range(5):
        first = limits[0][0] + index * (limits[0][1] - limits[0][0]) / 4
        second = limits[1][0] + index * (limits[1][1] - limits[1][0]) / 4
        body += f'<line x1="{x(first)}" y1="45" x2="{x(first)}" y2="240"/><line x1="65" y1="{y(second)}" x2="615" y2="{y(second)}"/>'
        body += _text(x(first), 258, f"{first:.2f}", text_anchor="middle") + _text(57, y(second) + 4, f"{second:.2f}", text_anchor="end")
    for row in sorted(centers, key=lambda item: item["role"] == "REFERENCE"):
        first, second = row["mean"]
        sd_x, sd_y = row["standard_deviation"]
        color = COLORS[row["group"]]
        body += f'<g style="color:{color}"><line x1="{x(first-sd_x)}" y1="{y(second)}" x2="{x(first+sd_x)}" y2="{y(second)}" style="stroke:currentColor;stroke-width:2"/><line x1="{x(first)}" y1="{y(second-sd_y)}" x2="{x(first)}" y2="{y(second+sd_y)}" style="stroke:currentColor;stroke-width:2"/>'
        body += f'<circle cx="{x(first)}" cy="{y(second)}" r="4" fill="{color}"/>' if row["role"] == "REFERENCE" else f'<rect x="{x(first)-6}" y="{y(second)-6}" width="12" height="12" fill="none" stroke="{color}" stroke-width="2"/>'
        body += '</g>'
    body += _text(340, 280, "PC1", text_anchor="middle") + _text(8, 37, "PC2")
    for index, (group, label) in enumerate(GROUP_LABELS.items()):
        body += f'<rect x="{30+index*160}" y="298" width="10" height="10" fill="{COLORS[group]}"/>' + _text(46 + index * 160, 308, label)
    body += _text(30, 330, "● Référence : A=40, B=40, mélanges=24, atypiques=4")
    body += _text(30, 348, "□ Projetés : A=8, B=8, atypiques=4 ; barres = dispersion ± 1 ET")
    return _svg("Position des groupes simulés et des sujets projetés", 360, body)


def scree(case: dict[str, Any]) -> str:
    body = _text(8, 18, "Part de variance expliquée — trois axes demandés")
    for value in (0, 25, 50, 75, 100):
        vertical = 145 - value
        body += f'<line x1="60" y1="{vertical}" x2="640" y2="{vertical}"/>' + _text(53, vertical + 4, f"{value} %", text_anchor="end")
    for index, row in enumerate(case["eigenvalues"]):
        percent = 100 * row["explained_variance_ratio"]
        horizontal = 150 + 175 * index
        body += f'<rect x="{horizontal-42}" y="{145-percent}" width="84" height="{percent}" fill="#286d91"/>'
        body += _text(horizontal, 139 - percent, f"{percent:.2f} %", text_anchor="middle") + _text(horizontal, 165, row["component"], text_anchor="middle")
    return _svg("Variance expliquée par PC1 à PC3", 180, body)


def distance_histogram(case: dict[str, Any]) -> str:
    histogram = case["distance_histogram"]
    edges, reference, projected = histogram["edges"], histogram["reference_counts"], histogram["projected_counts"]
    maximum = max(left + right for left, right in zip(reference, projected, strict=True))
    bin_width = 576 / len(reference)
    body = _text(8, 18, "Distances D² — référence et projetés, cas de base")
    for index, (left, right) in enumerate(zip(reference, projected, strict=True)):
        horizontal = 60 + index * bin_width
        height_reference, height_projected = 120 * left / maximum, 120 * right / maximum
        body += f'<rect x="{horizontal}" y="{160-height_reference}" width="{bin_width-5}" height="{height_reference}" fill="#286d91"/><rect x="{horizontal}" y="{160-height_reference-height_projected}" width="{bin_width-5}" height="{height_projected}" fill="#b06a20"/>'
        if left + right:
            body += _text(horizontal + (bin_width - 5) / 2, 154 - height_reference - height_projected, str(left + right), text_anchor="middle")
    position = 60 + 576 * histogram["threshold"] / edges[-1]
    body += f'<line x1="{position}" y1="35" x2="{position}" y2="160" style="stroke:#9b3e61;stroke-width:2;stroke-dasharray:5 4"/>'
    body += _text(position + 5, 47, f'Seuil D² = {histogram["threshold"]:.2f}')
    for fraction in (0, .25, .5, .75, 1):
        body += _text(60 + 576 * fraction, 178, f"{edges[-1]*fraction:.1f}", text_anchor="middle")
    body += _text(12, 36, "Effectif") + _text(600, 197, "D²")
    body += _text(65, 215, "Bleu : référence (108) · ocre : projetés (20) · nombres : sujets par classe")
    return _svg("Histogramme des distances PCA et seuil d’alerte", 230, body)


def detection_intervals(evaluation: dict[str, Any]) -> str:
    data = evaluation["metrics"]["independent_reference_detection"]
    body = _text(8, 18, "Incertitude — références indépendantes, cas de base")
    for horizontal in (0, 25, 50, 75, 100):
        x = 210 + horizontal * 4.1
        body += f'<line x1="{x}" y1="35" x2="{x}" y2="120"/>' + _text(x, 142, f"{horizontal} %", text_anchor="middle")
    for vertical, label, successes, count, interval in ((55, "Sensibilité", data["tp"], data["tp"] + data["fn"], data["sensitivity_wilson"]),
                                                       (105, "Fausses alertes", data["fp"], data["fp"] + data["tn"], data["non_injected_rate_wilson"])):
        start, end = (210 + 410 * value for value in interval)
        estimate = 210 + 410 * successes / count
        body += _text(8, vertical + 4, f"{label} : {successes}/{count}")
        body += f'<line x1="{start}" y1="{vertical}" x2="{end}" y2="{vertical}" style="stroke:#286d91;stroke-width:4"/><circle cx="{estimate}" cy="{vertical}" r="5" fill="#286d91"/>'
        body += _text(210, vertical + 21, f"Wilson 95 % : {100*interval[0]:.2f}–{100*interval[1]:.2f} %")
    return _svg("Rappels et fausses alertes avec intervalles descriptifs", 155, body)


def cohort_counts(case: dict[str, Any]) -> str:
    body = _text(8, 18, "Inclusions — parmi les 128 sujets, cas de base")
    maximum = case["sample_count"]
    for index, (cohort, label) in enumerate(zip(case["cohorts"], COHORT_LABELS, strict=True)):
        row = case["cohorts"][cohort]
        vertical, width = 40 + 30 * index, 360 * row["included"] / maximum
        body += _text(8, vertical + 14, label) + f'<rect x="210" y="{vertical}" width="360" height="18" fill="#e2e7eb"/><rect x="210" y="{vertical}" width="{width}" height="18" fill="#286d91"/>'
        body += _text(585, vertical + 14, f'{row["included"]} / {maximum}')
    body += _text(210, 268, "Bleu : inclus · gris : non inclus dans cette cohorte")
    return _svg("Effectifs des sept cohortes analytiques", 285, body)


def control_exclusions(case: dict[str, Any]) -> str:
    rows = case["cohorts"]["controls_unrelated"]["exclusion_counts"]
    body = _text(8, 18, "Pourquoi un sujet n’entre-t-il pas dans les témoins ?")
    for index, (reason, count) in enumerate(rows.items()):
        vertical = 40 + index * 30
        body += _text(8, vertical + 14, REASON_LABELS.get(reason, reason))
        body += f'<rect x="255" y="{vertical}" width="{330*count/case["sample_count"]}" height="18" fill="#b06a20"/>' + _text(600, vertical + 14, str(count))
    return _svg("Motifs exclusifs de non-inclusion des témoins", 65 + 30 * len(rows), body)
