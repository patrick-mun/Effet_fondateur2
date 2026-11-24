"""Rapport autonome CAL-REP-14B de la censure Gamma D25b."""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document


SCHEMAS = Path(__file__).resolve().parents[3] / "calibration/schemas"
HEADINGS = (
    "Identification et domaine de validité",
    "Résumé décisionnel",
    "Pourquoi une vraisemblance censurée",
    "Autorité, préspecification et aveuglement",
    "Plan expérimental D25b",
    "Modèle biologique et observation",
    "Vraisemblance et intervalle exact",
    "Fixture analytique et contrôles techniques",
    "Résultat primaire indépendant",
    "Résultats selon l’âge",
    "Résultats selon l’effectif familial",
    "Comparateurs et stress corrélé",
    "Critères, décision et usage autorisé",
    "Limites, provenance et empreintes",
)


def _num(value: float | None, digits: int = 2) -> str:
    return "Non évaluable" if value is None else f"{value:.{digits}f}"


def _pct(value: float | None) -> str:
    return "Non évaluable" if value is None else f"{100 * value:.2f} %"


def _figure(svg: str, caption: str, number: int) -> str:
    return (
        '<figure class="scientific-figure">' + svg
        + f"<figcaption>Figure {number} : {escape(caption)}</figcaption></figure>"
    )


def _table(
    headers: tuple[str, ...], rows: list[list[str]], caption: str, number: int,
) -> str:
    head = "".join(f"<th>{escape(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row)
        + "</tr>"
        for row in rows
    )
    return (
        f"<table><caption>Tableau {number} : {escape(caption)}</caption>"
        f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    )


def _method_label(method: str) -> str:
    return {
        "TYPE2_CENSORED_INDEPENDENT_V3": "V3 censurée exacte",
        "R_FAITHFUL_INDEPENDENT_V2": "V2 indépendante",
        "R_FAITHFUL_CORRELATED_V2": "V2 corrélée",
        "FULL_INDEPENDENT_ORACLE_V1": "Oracle complet",
    }[method]


def _genealogy_label(genealogy: str) -> str:
    return {
        "INDEPENDENT": "Indépendante",
        "CORRELATED_25": "Corrélée 25 %",
        "CORRELATED_50": "Corrélée 50 %",
    }[genealogy]


def _method_bias_figure(document: dict[str, Any]) -> str:
    rows = [
        row for row in document["method_results"]
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    maximum = max(6.0, max(abs(float(row["bias"])) for row in rows) * 1.2)
    zero_x, scale = 310.0, 175.0 / maximum
    parts = [
        '<svg viewBox="0 0 540 230" role="img" aria-label="Biais des méthodes indépendantes D25b">',
        '<text x="8" y="18" font-size="13" font-weight="600">Biais sur les généalogies indépendantes</text>',
        '<line x1="310" y1="32" x2="310" y2="200" stroke="#555"/>',
    ]
    for index, row in enumerate(rows):
        y = 48 + index * 38
        bias = float(row["bias"])
        width = abs(bias) * scale
        x = zero_x - width if bias < 0 else zero_x
        color = "#2f7d69" if row["method"] in {
            "TYPE2_CENSORED_INDEPENDENT_V3", "FULL_INDEPENDENT_ORACLE_V1"
        } else "#b24745"
        parts.append(
            f'<text x="8" y="{y + 11}" font-size="9">'
            f'{escape(_method_label(row["method"]))}</text>'
        )
        parts.append(
            f'<rect x="{x:.2f}" y="{y}" width="{width:.2f}" '
            f'height="15" fill="{color}" rx="2"/>'
        )
        label_x = x - 38 if bias < 0 else x + width + 6
        parts.append(
            f'<text x="{label_x:.2f}" y="{y + 11}" font-size="9">'
            f'{bias:.2f}</text>'
        )
    parts.append(
        '<text x="205" y="222" font-size="9">Biais signé (générations)</text></svg>'
    )
    return "".join(parts)


def _age_bias_figure(document: dict[str, Any]) -> str:
    rows = sorted(
        [
            row for row in document["results_by_age"]
            if row["genealogy_profile"] == "INDEPENDENT"
        ],
        key=lambda row: row["age_generations"],
    )
    values = [float(row["bias"]) for row in rows]
    limit = max(2.0, max(abs(value) for value in values) * 1.25)
    x_coord = lambda age: 65 + 430 * float(age) / 150
    y_coord = lambda bias: 135 - 95 * float(bias) / limit
    points = [
        (x_coord(row["age_generations"]), y_coord(row["bias"])) for row in rows
    ]
    parts = [
        '<svg viewBox="0 0 540 280" role="img" aria-label="Biais D25b selon l’âge">',
        '<text x="8" y="18" font-size="13" font-weight="600">Biais V3 selon l’âge vrai</text>',
        '<line x1="65" y1="230" x2="500" y2="230" stroke="#555"/>',
        '<line x1="65" y1="40" x2="65" y2="230" stroke="#555"/>',
        f'<line x1="65" y1="{y_coord(0):.2f}" x2="500" '
        f'y2="{y_coord(0):.2f}" stroke="#777" stroke-dasharray="4 3"/>',
        '<polyline points="' + " ".join(
            f"{x:.2f},{y:.2f}" for x, y in points
        ) + '" fill="none" stroke="#2f7d69" stroke-width="2"/>',
    ]
    for (x, y), row in zip(points, rows, strict=True):
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="#2f7d69"/>')
        parts.append(
            f'<text x="{x - 10:.2f}" y="{y - 8:.2f}" font-size="8">'
            f'{float(row["bias"]):.1f}</text>'
        )
    for age in (0, 50, 100, 150):
        parts.append(
            f'<text x="{x_coord(age) - 8:.1f}" y="248" font-size="9">{age}</text>'
        )
    parts.extend((
        '<text x="360" y="270" font-size="9">Âge vrai (générations)</text>',
        '<text x="14" y="175" font-size="9" transform="rotate(-90 14 175)">Biais (générations)</text>',
        '</svg>',
    ))
    return "".join(parts)


def _family_coverage_figure(document: dict[str, Any]) -> str:
    rows = [
        row for row in document["results_by_family_count"]
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    x0, width = 170.0, 330.0
    minimum, maximum = 0.92, 0.98
    coordinate = lambda value: x0 + width * (float(value) - minimum) / (
        maximum - minimum
    )
    parts = [
        '<svg viewBox="0 0 540 220" role="img" aria-label="Couverture D25b par effectif">',
        '<text x="8" y="18" font-size="13" font-weight="600">Couverture exacte par effectif</text>',
        f'<line x1="{coordinate(.95):.2f}" y1="30" '
        f'x2="{coordinate(.95):.2f}" y2="175" stroke="#b24745" '
        'stroke-dasharray="4 3"/>',
    ]
    for tick in (.92, .94, .95, .96, .98):
        x = coordinate(tick)
        parts.append(f'<text x="{x - 10:.2f}" y="205" font-size="9">{100*tick:.0f}%</text>')
    for index, row in enumerate(rows):
        y = 58 + index * 45
        low = coordinate(row["coverage_ci95_low"])
        high = coordinate(row["coverage_ci95_high"])
        value = coordinate(row["coverage"])
        parts.append(
            f'<text x="20" y="{y + 4}" font-size="10">'
            f'{row["family_count"]} familles</text>'
        )
        parts.append(
            f'<line x1="{low:.2f}" y1="{y}" x2="{high:.2f}" y2="{y}" '
            'stroke="#286d91" stroke-width="3"/>'
        )
        parts.append(f'<circle cx="{value:.2f}" cy="{y}" r="5" fill="#2f7d69"/>')
        parts.append(
            f'<text x="{value + 7:.2f}" y="{y - 7}" font-size="8">'
            f'{100*float(row["coverage"]):.2f}%</text>'
        )
    parts.append('</svg>')
    return "".join(parts)


def _correlated_coverage_figure(document: dict[str, Any]) -> str:
    rows = [
        row for row in document["results_by_family_count"]
        if row["genealogy_profile"] != "INDEPENDENT"
    ]
    series = {
        genealogy: sorted(
            [row for row in rows if row["genealogy_profile"] == genealogy],
            key=lambda row: row["family_count"],
        )
        for genealogy in ("CORRELATED_25", "CORRELATED_50")
    }
    x_coord = {3: 130.0, 5: 280.0, 10: 455.0}
    y_coord = lambda value: 220 - 170 * (float(value) - .60) / .40
    colors = {"CORRELATED_25": "#286d91", "CORRELATED_50": "#b24745"}
    parts = [
        '<svg viewBox="0 0 540 285" role="img" aria-label="Couverture D25b sous corrélation">',
        '<text x="8" y="18" font-size="13" font-weight="600">Rupture de couverture sous dépendance</text>',
        '<line x1="70" y1="220" x2="500" y2="220" stroke="#555"/>',
        '<line x1="70" y1="40" x2="70" y2="220" stroke="#555"/>',
        f'<line x1="70" y1="{y_coord(.95):.2f}" x2="500" '
        f'y2="{y_coord(.95):.2f}" stroke="#777" stroke-dasharray="4 3"/>',
    ]
    for tick in (.60, .70, .80, .90, .95, 1.00):
        y = y_coord(tick)
        parts.append(f'<text x="38" y="{y + 3:.1f}" font-size="9">{100*tick:.0f}%</text>')
    for genealogy, values in series.items():
        points = [
            (x_coord[row["family_count"]], y_coord(row["coverage"]))
            for row in values
        ]
        parts.append(
            '<polyline points="' + " ".join(
                f"{x:.2f},{y:.2f}" for x, y in points
            ) + f'" fill="none" stroke="{colors[genealogy]}" stroke-width="2"/>'
        )
        for x, y in points:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{colors[genealogy]}"/>')
    for family, x in x_coord.items():
        parts.append(f'<text x="{x - 5}" y="240" font-size="9">{family}</text>')
    parts.extend((
        '<rect x="125" y="258" width="12" height="7" fill="#286d91"/><text x="142" y="265" font-size="9">chemin 25 %</text>',
        '<rect x="285" y="258" width="12" height="7" fill="#b24745"/><text x="302" y="265" font-size="9">chemin 50 %</text>',
        '<text x="420" y="265" font-size="9">Familles</text>',
        '</svg>',
    ))
    return "".join(parts)


def _equations() -> str:
    equations = (
        ("Bras vrai", "X_f ~ Exp(g)", "X est en Morgans et g en générations."),
        ("Temps sous risque d’un côté", "T_s = Σ_{i=1}^{n−1} X_(i) + X_(n−1)", "Le dernier terme représente le bras le plus long censuré."),
        ("Vraisemblance", "ℓ(g) = k log(g) − gT + C ; k = 2(n−1)", "Les deux côtés apportent chacun n−1 ruptures."),
        ("Estimateurs", "ĝ_MLE = k/T ; ĝ_sans_biais = (k−1)/T", "D25b publie l’estimateur sans biais."),
        ("Intervalle exact", "IC95 = [G_k⁻¹(0,025)/T ; G_k⁻¹(0,975)/T]", "G_k est une Gamma de forme k et de taux 1."),
        ("Évaluation", "biais = mean(ĝ−g) ; RMSE = sqrt(mean((ĝ−g)²)) ; C = mean(1[g∈IC95])", "Le scénario simulé est l’unité indépendante."),
    )
    return "".join(
        f'<div class="equation"><p><strong>Équation {index} — '
        f'{escape(title)}</strong></p><div class="math" role="math">'
        f'{escape(formula)}</div><p>{escape(detail)}</p></div>'
        for index, (title, formula, detail) in enumerate(equations, 1)
    )


def _sections(document: dict[str, Any], source_sha: str) -> list[str]:
    primary = document["primary_result"]
    independent_methods = [
        row for row in document["method_results"]
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    independent_family = [
        row for row in document["results_by_family_count"]
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    correlated_family = [
        row for row in document["results_by_family_count"]
        if row["genealogy_profile"] != "INDEPENDENT"
    ]
    design_rows = [
        ["Âges", ", ".join(map(str, document["design"]["ages_generations"])), "générations"],
        ["Effectifs", ", ".join(map(str, document["design"]["family_counts"])), "familles"],
        ["Réplications", str(document["design"]["replications_per_cell"]), "par cellule"],
        ["Scénarios", str(document["design"]["scenario_count"]), "unités indépendantes"],
        ["Estimations", str(document["design"]["estimate_row_count"]), "comparaisons appariées"],
    ]
    method_rows = [[
        _method_label(row["method"]), _num(row["bias"]), _num(row["rmse"]),
        _pct(row["coverage"]), _num(row["mean_interval_length_generations"]),
    ] for row in independent_methods]
    family_rows = [[
        str(row["family_count"]), _num(row["bias"]), _num(row["rmse"]),
        _pct(row["coverage"]),
        f'[{_pct(row["coverage_ci95_low"])}, {_pct(row["coverage_ci95_high"])}]',
    ] for row in independent_family]
    correlated_rows = [[
        _genealogy_label(row["genealogy_profile"]), str(row["family_count"]),
        _num(row["bias"]), _num(row["rmse"]), _pct(row["coverage"]),
    ] for row in correlated_family]
    labels = {
        "analytic_fixture_absolute_error": "Erreur fixture analytique",
        "absolute_primary_bias": "|biais| primaire",
        "primary_coverage": "Couverture primaire",
        "maximum_absolute_family_bias": "|biais| maximal par effectif",
        "maximum_family_coverage_deviation_from_0_95": "Écart max à 95 %",
        "primary_not_estimated_count": "Non estimés primaires",
        "correlated_stress_performance": "Stress corrélé",
    }
    criterion_rows = []
    for row in document["criteria"]:
        formatter = (
            (lambda value: "Non évaluable" if value is None else f"{value:.3e}")
            if row["criterion_id"] == "CAL-AGE-B01"
            else (lambda value: "Non évaluable" if value is None else f"{value:.4f}")
        )
        criterion_rows.append([
            row["criterion_id"], labels[row["metric"]], formatter(row["value"]),
            formatter(row["lower_threshold"]), formatter(row["upper_threshold"]),
            "REVUE REQUISE" if row["status"] == "REVIEW_REQUIRED" else row["status"],
        ])
    correlated_25 = next(
        row for row in document["method_results"]
        if row["genealogy_profile"] == "CORRELATED_25"
        and row["method"] == "TYPE2_CENSORED_INDEPENDENT_V3"
    )
    correlated_50 = next(
        row for row in document["method_results"]
        if row["genealogy_profile"] == "CORRELATED_50"
        and row["method"] == "TYPE2_CENSORED_INDEPENDENT_V3"
    )
    return [
        '<p>CAL-REP-14B évalue une nouvelle vraisemblance pour les longueurs ancestrales observables seulement par partage entre au moins deux porteurs. Le résultat qualifiable concerne uniquement les unités familiales à histoires recombinatoires indépendantes.</p><p>La candidate V3 ne remplace ni V1 ni V2. Aucune donnée DOCK6 réelle, clinique ou individuelle n’est utilisée.</p>',
        f'<p>Statut technique : <strong>{document["technical_status"]}</strong>. Contrat indépendant : <strong>{document["independent_contract_qualification"]}</strong>. Qualification globale : <strong>{document["performance_qualification"]}</strong>.</p><p>Le biais primaire vaut {_num(primary["bias"])} génération, IC95 [{_num(primary["bias_ci95_low"])}, {_num(primary["bias_ci95_high"])}], et la couverture atteint {_pct(primary["coverage"])}.</p><aside class="callout">La censure est corrigée pour l’hypothèse indépendante, mais la couverture s’effondre sous forte dépendance. Production : NOT_PROMOTED ; calibrated=false ; founder_effect_proven=false.</aside>',
        '<p>La méthode publiée ajoutait une longueur moyenne aux deux bras les plus longs. D25a a montré que cette approximation surestimait encore l’âge de 5,91 générations sur les longueurs tronquées. D25b utilise plutôt l’information réellement disponible : n−1 ruptures exactes et un bras encore à risque au dernier seuil observé.</p><p>Cette représentation est une censure de type II classique. Elle évite d’inventer la valeur du prolongement invisible.</p>',
        '<p>Le protocole, la configuration, les graines et les sept décisions ont été fermés avant la simulation. Les 36 000 âges vrais sont séparés des entrées aveugles ; les 144 000 estimations sont fermées avant ouverture de la vérité.</p><p>Le rapport lit uniquement le document agrégé CAL-REP-14B.json et ne recalcule aucune métrique scientifique.</p>',
        '<p>La grille reprend les dix âges et les trois effectifs de D25a, double le nombre de réplications et utilise deux nouvelles graines. Les profils corrélés ne participent pas aux critères de succès indépendants.</p>' + _table(("Dimension", "Valeur", "Unité/rôle"), design_rows, "Plan fermé D25b.", 1),
        '<p>Sous le modèle de Haldane, chaque bras est exponentiel de taux égal à l’âge. Le partage entre au moins deux porteurs révèle toutes les ruptures sauf la plus lointaine : le maximum observé est donc répété, une fois comme rupture exacte et une fois comme exposition censurée.</p><p>Le contrat refuse une correction de partage fortuit et toute censure multiple. Une égalité due à la résolution des marqueurs réels ne doit pas être confondue avec ce mécanisme biologique.</p>',
        '<p>Le temps total sous risque suffit pour estimer le taux exponentiel. La forme Gamma k=2(n−1) donne un point sans biais et un intervalle exact, sans approximation asymptotique.</p>' + _equations(),
        f'<p>La fixture contient trois familles, quatre événements et un temps total sous risque de {document["analytic_fixture"]["total_exposure_morgans"]:.2f} Morgan. La valeur fermée est {document["analytic_fixture"]["expected_generations"]:.2f} générations ; la candidate retourne {document["analytic_fixture"]["observed_generations"]:.15f}.</p><p>L’erreur absolue {document["analytic_fixture"]["absolute_error_generations"]:.3e} génération est inférieure au seuil 1e−12. Les données sans maximum répété et les corrections non modélisées sont refusées.</p>',
        f'<p>Sur {primary["estimated_count"]} scénarios indépendants, le biais est {_num(primary["bias"])} génération et la RMSE {_num(primary["rmse"])}. La couverture {_pct(primary["coverage"])} possède un IC95 de Wilson [{_pct(primary["coverage_ci95_low"])}, {_pct(primary["coverage_ci95_high"])}].</p><p>L’intervalle moyen reste large ({_num(primary["mean_interval_length_generations"])} générations) : absence de biais ne signifie pas précision individuelle.</p>' + _figure(_method_bias_figure(document), "Biais des quatre méthodes appariées sur les 12 000 scénarios indépendants. Source : CAL-REP-14B.json.", 1) + _table(("Méthode", "Biais", "RMSE", "Couverture", "Longueur IC"), method_rows, "Comparaison des méthodes en généalogie indépendante ; unités en générations.", 2),
        '<p>Le biais reste fluctuant selon l’âge, mais aucune dérive croissante systématique comparable à D25 n’apparaît. Les valeurs aux âges élevés ont naturellement une variance absolue plus grande, car l’échelle exponentielle est proportionnelle à l’âge.</p>' + _figure(_age_bias_figure(document), "Biais signé de la candidate V3 par âge vrai en généalogie indépendante. Source : CAL-REP-14B.json.", 2),
        f'<p>Les trois effectifs respectent les limites préspecifiées. Le biais maximal vaut {max(abs(float(row["bias"])) for row in independent_family):.2f} génération et l’écart maximal de couverture à 95 % vaut {100*max(abs(float(row["coverage"])-.95) for row in independent_family):.2f} point.</p><p>La précision reste fortement dépendante de l’effectif : la RMSE passe de {_num(independent_family[0]["rmse"])} générations à trois familles à {_num(independent_family[-1]["rmse"])} à dix. Le P95 leave-one-family-out atteint {_num(primary["lofo_p95_max_absolute_delta_generations"])} générations.</p>' + _figure(_family_coverage_figure(document), "Couverture et IC de Wilson de la candidate V3 par nombre de familles. Trait rouge : 95 %. Source : CAL-REP-14B.json.", 3) + _table(("Familles", "Biais", "RMSE", "Couverture", "IC95 Wilson"), family_rows, "Résultats indépendants par effectif ; unités en générations.", 3),
        f'<p>La V2 indépendante conserve un biais de {_num(independent_methods[1]["bias"])} générations, tandis que la V3 rejoint l’oracle complet. Ce gain disparaît dès que l’hypothèse d’indépendance est violée.</p><p>Avec 25 % de chemin partagé, la V3 couvre {_pct(correlated_25["coverage"])} ; avec 50 %, seulement {_pct(correlated_50["coverage"])}. À dix familles et 50 %, la couverture tombe à {_pct(correlated_family[-1]["coverage"])} et le biais atteint {_num(correlated_family[-1]["bias"])} générations.</p>' + _figure(_correlated_coverage_figure(document), "Couverture de la candidate indépendante sous deux niveaux de chemin partagé. Source : CAL-REP-14B.json.", 4) + _table(("Profil", "Familles", "Biais", "RMSE", "Couverture"), correlated_rows, "Stress de dépendance pour la candidate V3 ; unités en générations.", 4),
        f'<p>Les critères B01 à B06 passent ; B07 reste en revue humaine par construction. La qualification est donc {document["performance_qualification"]}, sans promotion de production.</p><p>Usage autorisé : référence ou sensibilité indépendante lorsque les unités familiales peuvent être défendues comme recombinatoirement indépendantes. Usage interdit : estimation principale corrélée ou preuve d’effet fondateur.</p>' + _table(("Critère", "Métrique", "Valeur", "Seuil bas", "Seuil haut", "Décision"), criterion_rows, "Décisions D25b préspecifiées.", 5),
        f'<p>D25b ne mesure ni erreur SHAPEIT5, ni résolution de marqueurs, ni partage fortuit, ni faux positifs D24. Le succès indépendant ne valide donc pas le pipeline complet. Une méthode corrélée censurée nécessiterait un modèle de dépendance identifiable et un nouveau panel.</p><p>Évaluation close : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Banc : <span class="hash">{document["bench_manifest_sha256"]}</span>. Run : <span class="hash">{document["run_manifest_sha256"]}</span>.</p>',
    ]


def build_report(
    *, evaluation_root: Path, template_path: Path, output_root: Path,
) -> Path:
    """Produit CAL-REP-14B depuis la seule évaluation agrégée close."""
    source = evaluation_root / "CAL-REP-14B.json"
    checksum = evaluation_root / "checksums.sha256"
    if (
        not source.is_file() or source.is_symlink()
        or not checksum.is_file() or checksum.is_symlink()
    ):
        raise ValueError("d25b_report_source_missing")
    expected_checksum = f"{sha256_file(source)}  CAL-REP-14B.json"
    if checksum.read_text(encoding="utf-8").strip() != expected_checksum:
        raise ValueError("d25b_report_source_checksum_invalid")
    document = json.loads(source.read_text(encoding="utf-8"))
    validate_json_document(
        document, "gamma_censoring_evaluation.schema.json", schemas_dir=SCHEMAS
    )
    source_sha = sha256_file(source)
    template = template_path.read_text(encoding="utf-8")
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d25b_report_template_invalid")
    style = (
        '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table'
        '{break-inside:avoid}.scientific-figure{margin:1em 0}.equation'
        '{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;'
        'margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}'
        'table{font-size:.68em}main p{orphans:3;widows:3}</style>'
    )
    head = re.sub(
        r"<title>.*?</title>",
        "<title>CAL-REP-14B — Censure Gamma D25b</title>",
        head, flags=re.S,
    ).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}">'
        f'<span class="toc-title">{index}. {escape(title)}</span>'
        '<span class="toc-dots"></span><span class="toc-page-num"></span>'
        '</a></li>'
        for index, title in enumerate(HEADINGS, 1)
    ) + "</ol></nav>"
    body = (
        '<body lang="fr"><section class="cover">'
        '<h1>CAL-REP-14B — Censure Gamma D25b</h1>'
        '<p>Effet_fondateur2 · version 1.0.0 · 10 novembre 2026</p>'
        '<p>Vraisemblance exacte indépendante et stress de dépendance</p>'
        '</section>' + toc + "<main>"
    )
    body += "".join(
        f'<section id="s{index}"><h2>{index}. {escape(title)}</h2>'
        f'{content}</section>'
        for index, (title, content) in enumerate(
            zip(HEADINGS, _sections(document, source_sha), strict=True), 1
        )
    )
    html = (
        head + body
        + '</main><footer class="page-footer"><span>CAL-REP-14B · '
        'Effet_fondateur2</span></footer></body></html>'
    )
    destination = new_workspace_run(
        output_root.absolute().parent, output_root.name
    )
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{output_root.name}.", dir=destination.parent
    ))
    try:
        report_dir = temporary / "CAL-REP-14B"
        report_dir.mkdir()
        html_path = report_dir / "index.html"
        html_path.write_text(html, encoding="utf-8")
        data_path = report_dir / "report_data.json"
        write_json(data_path, {"results": document})
        (report_dir / "checksums.sha256").write_text(
            "".join(
                f"{sha256_file(path)}  {path.name}\n"
                for path in (html_path, data_path)
            ),
            encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    if sha256_file(source) != source_sha:
        raise ValueError("d25b_report_source_changed")
    return destination / "CAL-REP-14B/index.html"


def main() -> int:
    parser = argparse.ArgumentParser(description="Rapport CAL-REP-14B")
    for name in ("evaluation-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    print(build_report(**vars(parser.parse_args())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
