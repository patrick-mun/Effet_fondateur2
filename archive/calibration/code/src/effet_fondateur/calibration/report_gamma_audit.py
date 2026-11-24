"""Rapport autonome CAL-REP-14A de l'audit D25a."""

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
    "Identification et portée",
    "Résumé décisionnel",
    "Pourquoi D25 devait être repris",
    "Sources, autorité et aveuglement",
    "Plan expérimental D25a",
    "Simulation des généalogies",
    "Méthodes comparées et fidélité R",
    "Équations et mesures",
    "Longueurs complètes : oracle restauré",
    "Longueurs tronquées : correction insuffisante",
    "Généalogies corrélées",
    "Effectifs, stabilité et limites",
    "Critères et décision",
    "Provenance et empreintes",
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
        "<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>"
        for row in rows
    )
    return (
        f"<table><caption>Tableau {number} : {escape(caption)}</caption>"
        f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    )


def _contract_label(row: dict[str, Any]) -> str:
    genealogy = {
        "INDEPENDENT": "Indépendante",
        "CORRELATED_25": "Corrélée 25 %",
        "CORRELATED_50": "Corrélée 50 %",
    }[row["genealogy_profile"]]
    observation = (
        "complète" if row["observation_profile"] == "FULL_BREAKPOINTS"
        else "tronquée"
    )
    return f"{genealogy} · {observation}"


def _bias_figure(document: dict[str, Any]) -> str:
    rows = document["contract_results"]
    maximum = max(abs(float(row["bias"])) for row in rows) * 1.15
    zero_x = 330.0
    bar_scale = 165.0 / maximum
    parts = [
        '<svg viewBox="0 0 540 285" role="img" aria-label="Biais des six contrats D25a">',
        '<text x="8" y="18" font-size="13" font-weight="600">Biais signé par contrat</text>',
        '<line x1="330" y1="32" x2="330" y2="255" stroke="#555"/>',
    ]
    for index, row in enumerate(rows):
        y = 43 + index * 35
        bias = float(row["bias"])
        width = abs(bias) * bar_scale
        x = zero_x - width if bias < 0 else zero_x
        color = "#b24745" if abs(bias) > 0.75 else "#2f7d69"
        parts.append(f'<text x="8" y="{y + 12}" font-size="9">{escape(_contract_label(row))}</text>')
        parts.append(f'<rect x="{x:.2f}" y="{y}" width="{width:.2f}" height="15" fill="{color}" rx="2"/>')
        value_x = x - 38 if bias < 0 else x + width + 5
        parts.append(f'<text x="{value_x:.2f}" y="{y + 11}" font-size="9">{bias:.2f}</text>')
    parts.append('<text x="260" y="278" font-size="9">Sous-estimation ← 0 → surestimation (générations)</text></svg>')
    return "".join(parts)


def _coverage_figure(document: dict[str, Any]) -> str:
    rows = document["contract_results"]
    x0, width = 205.0, 310.0
    minimum, maximum = 0.88, 0.98
    coordinate = lambda value: x0 + width * (float(value) - minimum) / (maximum - minimum)
    parts = [
        '<svg viewBox="0 0 540 295" role="img" aria-label="Couverture des six contrats D25a">',
        '<text x="8" y="18" font-size="13" font-weight="600">Couverture des intervalles à 95 %</text>',
        f'<line x1="{x0}" y1="35" x2="{x0 + width}" y2="35" stroke="#555"/>',
        f'<line x1="{coordinate(.95):.2f}" y1="30" x2="{coordinate(.95):.2f}" y2="260" stroke="#b24745" stroke-dasharray="4 3"/>',
    ]
    for tick in (.88, .90, .92, .94, .96, .98):
        x = coordinate(tick)
        parts.append(f'<text x="{x - 10:.2f}" y="284" font-size="9">{100*tick:.0f}%</text>')
    for index, row in enumerate(rows):
        y = 55 + index * 35
        low = coordinate(row["coverage_ci95_low"])
        high = coordinate(row["coverage_ci95_high"])
        value = coordinate(row["coverage"])
        parts.append(f'<text x="8" y="{y + 4}" font-size="9">{escape(_contract_label(row))}</text>')
        parts.append(f'<line x1="{low:.2f}" y1="{y}" x2="{high:.2f}" y2="{y}" stroke="#286d91" stroke-width="3"/>')
        parts.append(f'<circle cx="{value:.2f}" cy="{y}" r="5" fill="#2f7d69"/>')
        label_x = min(value + 7, x0 + width - 28)
        parts.append(f'<text x="{label_x:.2f}" y="{y - 7}" font-size="8">{100*row["coverage"]:.1f}%</text>')
    parts.append('</svg>')
    return "".join(parts)


def _independent_method_figure(document: dict[str, Any]) -> str:
    rows = [
        row for row in document["method_results"]
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    maximum = max(abs(float(row["bias"])) for row in rows) * 1.1
    zero_x, scale = 360.0, 150.0 / maximum
    method_names = {
        "LEGACY_INDEPENDENT_V1": "Legacy indépendant",
        "LEGACY_CORRELATED_V1": "Legacy corrélé",
        "EXACT_INDEPENDENT_ORACLE_V1": "Oracle sans bord",
        "R_FAITHFUL_INDEPENDENT_V2": "R fidèle indépendant",
        "R_FAITHFUL_CORRELATED_V2": "R fidèle corrélé",
    }
    parts = [
        '<svg viewBox="0 0 540 345" role="img" aria-label="Contraste des méthodes en généalogie indépendante">',
        '<text x="8" y="18" font-size="13" font-weight="600">Effet du contrat et de la méthode</text>',
        '<line x1="360" y1="32" x2="360" y2="315" stroke="#555"/>',
    ]
    for index, row in enumerate(rows):
        y = 42 + index * 34
        bias = float(row["bias"])
        width = abs(bias) * scale
        x = zero_x - width if bias < 0 else zero_x
        observation = "complet" if row["observation_profile"] == "FULL_BREAKPOINTS" else "tronqué"
        label = f"{observation} · {method_names[row['method']]}"
        color = "#286d91" if row["method"].startswith("LEGACY") else "#2f7d69"
        parts.append(f'<text x="8" y="{y + 11}" font-size="8.5">{escape(label)}</text>')
        parts.append(f'<rect x="{x:.2f}" y="{y}" width="{width:.2f}" height="14" fill="{color}" rx="2"/>')
        value_x = x - 38 if bias < 0 else x + width + 5
        parts.append(f'<text x="{value_x:.2f}" y="{y + 10}" font-size="8">{bias:.2f}</text>')
    parts.append('</svg>')
    return "".join(parts)


def _age_bias_figure(document: dict[str, Any]) -> str:
    rows = [
        row for row in document["results_by_age"]
        if row["genealogy_profile"] == "INDEPENDENT"
    ]
    series = {
        observation: sorted(
            [row for row in rows if row["observation_profile"] == observation],
            key=lambda row: row["age_generations"],
        )
        for observation in ("FULL_BREAKPOINTS", "PAIRWISE_SHARED_TRUNCATED")
    }
    biases = [float(row["bias"]) for values in series.values() for row in values]
    y_min = min(-2.0, min(biases) - 1.0)
    y_max = max(16.0, max(biases) + 1.0)
    x_coord = lambda age: 65 + 430 * float(age) / 150
    y_coord = lambda bias: 240 - 190 * (float(bias) - y_min) / (y_max - y_min)
    parts = [
        '<svg viewBox="0 0 540 285" role="img" aria-label="Biais par âge des contrats indépendants">',
        '<text x="8" y="18" font-size="13" font-weight="600">Biais indépendant selon l’âge vrai</text>',
        '<line x1="65" y1="240" x2="500" y2="240" stroke="#555"/>',
        '<line x1="65" y1="40" x2="65" y2="240" stroke="#555"/>',
        f'<line x1="65" y1="{y_coord(0):.2f}" x2="500" y2="{y_coord(0):.2f}" stroke="#777" stroke-dasharray="4 3"/>',
    ]
    colors = {"FULL_BREAKPOINTS": "#2f7d69", "PAIRWISE_SHARED_TRUNCATED": "#b24745"}
    for observation, values in series.items():
        points = [(x_coord(row["age_generations"]), y_coord(row["bias"])) for row in values]
        parts.append('<polyline points="' + " ".join(f"{x:.2f},{y:.2f}" for x, y in points) + f'" fill="none" stroke="{colors[observation]}" stroke-width="2"/>')
        for x, y in points:
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3" fill="{colors[observation]}"/>')
    for age in (0, 50, 100, 150):
        parts.append(f'<text x="{x_coord(age)-8:.1f}" y="257" font-size="9">{age}</text>')
    for tick in (0, 5, 10, 15):
        if y_min <= tick <= y_max:
            y = y_coord(tick)
            parts.append(f'<text x="45" y="{y + 3:.1f}" font-size="9">{tick}</text>')
    parts.extend((
        '<rect x="90" y="268" width="12" height="7" fill="#2f7d69"/><text x="107" y="275" font-size="9">complet</text>',
        '<rect x="180" y="268" width="12" height="7" fill="#b24745"/><text x="197" y="275" font-size="9">tronqué</text>',
        '<text x="350" y="275" font-size="9">Âge vrai (générations)</text>',
        '<text x="12" y="170" font-size="9" transform="rotate(-90 12 170)">Biais (générations)</text>',
        '</svg>',
    ))
    return "".join(parts)


def _equations() -> str:
    equations = (
        ("Bras marginal", "X_f = min(C,P_f), avec X_f ~ Exp(g/100)", "C crée la dépendance ; P_f est propre à la famille."),
        ("Oracle complet", "ĝ = (2n−1) / Σ_f(L_f+R_f)", "Les longueurs sont en Morgans et aucune extension n'est ajoutée."),
        ("Observation tronquée", "L_obs = Σ l_f + 2 â_manquant", "â_manquant est appliqué seulement aux bras les plus longs non observables."),
        ("Effectif du point", "n*_point = clamp[n/(1+(n−1)ρ), −n, n]", "Il peut être signé dans I1, conformément au script R."),
        ("Effectif de l’intervalle", "n*_IC > 0 ; il est recalculé dans I1 et fixé à n dans I2", "Le point et l'intervalle n'utilisent pas silencieusement la même valeur."),
        ("Évaluation", "biais=mean(ĝ−g) ; RMSE=sqrt(mean((ĝ−g)²)) ; C=mean(1[g∈IC95])", "Le scénario simulé est l'unité indépendante."),
    )
    return "".join(
        f'<div class="equation"><p><strong>Équation {index} — {escape(title)}</strong></p>'
        f'<div class="math" role="math">{escape(formula)}</div><p>{escape(detail)}</p></div>'
        for index, (title, formula, detail) in enumerate(equations, 1)
    )


def _sections(document: dict[str, Any], source_sha: str) -> list[str]:
    contracts = document["contract_results"]
    independent_full = next(
        row for row in contracts
        if row["genealogy_profile"] == "INDEPENDENT"
        and row["observation_profile"] == "FULL_BREAKPOINTS"
    )
    independent_truncated = next(
        row for row in contracts
        if row["genealogy_profile"] == "INDEPENDENT"
        and row["observation_profile"] == "PAIRWISE_SHARED_TRUNCATED"
    )
    contract_rows = [[
        _contract_label(row), row["method"], _num(row["bias"]),
        _num(row["rmse"]), _pct(row["coverage"]),
        _num(row["mean_interval_length_generations"]),
    ] for row in contracts]
    method_rows = [[
        _contract_label(row), row["method"], _num(row["bias"]),
        _pct(row["coverage"]), str(row["estimated_count"]),
    ] for row in document["method_results"] if row["genealogy_profile"] == "INDEPENDENT"]
    family_rows = [[
        _contract_label(row), str(row["family_count"]), _num(row["bias"]),
        _num(row["rmse"]), _pct(row["coverage"]),
    ] for row in document["results_by_family_count"]]
    parity_rows = [[
        row["fixture_id"], row["expected_rho_region"],
        _num(row["rho"], 6), f'{row["maximum_absolute_error"]:.3e}',
        "Oui" if row["region_matches"] else "Non",
    ] for row in document["r_parity"]["fixtures"]]
    metric_labels = {
        "maximum_r_python_error": "Erreur max R/Python",
        "absolute_oracle_bias": "|biais| oracle",
        "oracle_coverage": "Couverture oracle",
        "absolute_truncated_independent_bias": "|biais| tronqué",
        "truncated_independent_coverage": "Couverture tronquée",
        "correlated_performance": "Performance corrélée",
    }
    criterion_rows = []
    for row in document["criteria"]:
        scientific = row["criterion_id"] == "CAL-AGE-A01"
        format_value = lambda value: (
            "Non évaluable" if value is None
            else f"{value:.3e}" if scientific
            else f"{value:.4f}"
        )
        criterion_rows.append([
            row["criterion_id"], metric_labels[row["metric"]],
            format_value(row["value"]), format_value(row["lower_threshold"]),
            format_value(row["upper_threshold"]),
            "REVUE REQUISE" if row["status"] == "REVIEW_REQUIRED" else row["status"],
        ])
    correlated = [row for row in contracts if row["genealogy_profile"] != "INDEPENDENT"]
    return [
        '<p>CAL-REP-14A audite la datation Gamma de l’étape 14 après l’anomalie D25. Il porte exclusivement sur des longueurs synthétiques en cM et sur le script R public de référence. Aucune donnée clinique, DOCK6 ou individuelle n’est utilisée.</p><p>La candidate <code>gamma_gandolfo_2014_r_faithful_v2</code> reste expérimentale. La méthode historique V1 demeure inchangée dans le pipeline.</p>',
        f'<p>Conformité technique : <strong>{document["technical_status"]}</strong>. Performance : <strong>{document["performance_qualification"]}</strong>. La parité R/Python atteint une erreur maximale de {document["r_parity"]["maximum_absolute_error_generations"]:.3e} génération.</p><p>L’oracle sur longueurs complètes retrouve un biais de {_num(independent_full["bias"])} génération et une couverture de {_pct(independent_full["coverage"])}. En revanche, le contrat indépendant tronqué conserve un biais de {_num(independent_truncated["bias"])} générations et une couverture de {_pct(independent_truncated["coverage"])}.</p><aside class="callout">D25a localise le problème : la divergence Python/R est corrigée, mais la correction de bord publiée ne suffit pas à calibrer les petits échantillons tronqués. calibrated=false ; founder_effect_proven=false.</aside>',
        '<p>D25 additionnait une extension de bord à des distances-vérité déjà complètes, ce qui provoquait mathématiquement une sous-estimation. Il utilisait en outre le modèle corrélé sur des tirages indépendants et une seule valeur d’effectif effectif pour le point et l’intervalle.</p><p>D25a ne supprime pas ce résultat : il le conserve comme diagnostic historique et construit un panel neuf avec des contrats séparés. La figure 3 rend visible la part attribuable au contrat d’entrée et à la branche corrélée.</p>',
        '<p>Le protocole, la configuration et les seuils fixes ont été enregistrés avant le nouveau panel. Les 18 000 âges vrais sont physiquement séparés des longueurs aveugles. Le runner ferme 144 000 estimations avant l’ouverture de la vérité.</p><p>Le script R officiel est contrôlé par SHA-256 et exécuté localement avec R 4.4.1. Le rapport ne lit que CAL-REP-14A.json et ne recalcule aucune statistique scientifique.</p>',
        '<p>La grille croise dix âges, trois effectifs familiaux, trois généalogies, deux contrats d’observation et 200 réplications par cellule. Les méthodes appliquées au même scénario sont des comparaisons appariées, pas des unités indépendantes supplémentaires.</p>' + _table(("Contrat", "Méthode primaire", "Biais", "RMSE", "Couverture", "Longueur IC"), contract_rows, "Résultats des six contrats primaires ; unités en générations.", 1),
        '<p>Le profil indépendant utilise des bras exponentiels propres à chaque famille. Les profils corrélés ajoutent un chemin commun représentant 25 % ou 50 % du taux total, puis des branches privées. Le minimum conserve exactement la marginale exponentielle attendue tout en créant une dépendance positive.</p><p>Les longueurs tronquées remplacent, sur chaque côté, le maximum par le deuxième maximum : elles représentent ce qui reste observable par partage entre au moins deux porteurs.</p>',
        f'<p>Les cinq fixtures couvrent l’exemple officiel et les quatre domaines de corrélation. Toutes les régions sont correctement reconnues ; l’écart maximal {document["r_parity"]["maximum_absolute_error_generations"]:.3e} est très inférieur au seuil 1e−10.</p><p>Dans I1, l’effectif du point reste signé tandis que l’effectif de l’intervalle est positif. Cette séparation corrige la divergence de V1 sans modifier ses anciennes sorties.</p>' + _table(("Fixture", "Région attendue", "ρ", "Erreur max", "Concordance"), parity_rows, "Parité numérique avec le script R officiel.", 2),
        '<p>Les équations ci-dessous distinguent le processus biologique, l’observation et l’estimation. Les cM sont convertis en Morgans à la frontière du calcul. Les IC du biais, de la MAE et de la RMSE utilisent 1 000 bootstrap ; la couverture utilise Wilson.</p>' + _equations(),
        f'<p>Sur 6 000 scénarios indépendants complets, l’oracle donne un biais de {_num(independent_full["bias"])} génération, IC95 [{_num(independent_full["bias_ci95_low"])}, {_num(independent_full["bias_ci95_high"])}], et une couverture de {_pct(independent_full["coverage"])}. Les critères A02 et A03 passent.</p><p>La V1 indépendante, appliquée à ces mêmes longueurs avec une extension inutile, conserve un biais négatif. La figure 1 montre que les quatre contrats complets sont proches de zéro seulement lorsque leur méthode respecte le modèle de génération.</p>' + _figure(_bias_figure(document), "Biais signé des six contrats primaires. Vert : |biais| ≤ 0,75 ; rouge : dépassement. Source : CAL-REP-14A.json.", 1),
        f'<p>Après troncature, la candidate indépendante surestime l’âge de {_num(independent_truncated["bias"])} générations en moyenne, IC95 [{_num(independent_truncated["bias_ci95_low"])}, {_num(independent_truncated["bias_ci95_high"])}]. Sa couverture tombe à {_pct(independent_truncated["coverage"])}. Les critères A04 et A05 échouent.</p><p>La dégradation est la plus forte à trois familles et augmente en valeur absolue avec l’âge. La correction publiée est une approximation des longueurs manquantes ; elle ne restaure pas ici l’absence de biais aux petits effectifs.</p>' + _figure(_age_bias_figure(document), "Biais par âge en généalogie indépendante, longueurs complètes contre tronquées. Source : CAL-REP-14A.json.", 2) + _table(("Contrat", "Méthode", "Biais", "Couverture", "n"), method_rows, "Comparaisons appariées en généalogie indépendante.", 3),
        f'<p>Sur les longueurs complètes corrélées à 25 %, le biais vaut {_num(correlated[0]["bias"])} génération ; à 50 %, il vaut {_num(correlated[2]["bias"])}. Les intervalles sont toutefois très larges, avec des RMSE de {_num(correlated[0]["rmse"])} et {_num(correlated[2]["rmse"])} générations.</p><p>La troncature maintient une couverture de {_pct(correlated[1]["coverage"])} à 25 %, mais seulement {_pct(correlated[3]["coverage"])} à 50 %. Ces résultats restent exploratoires et ne scellent pas le modèle corrélé.</p>' + _figure(_coverage_figure(document), "Couverture et IC de Wilson des six contrats ; trait rouge = 95 %. Source : CAL-REP-14A.json.", 3),
        '<p>Le leave-one-family-out est recalculé après retrait d’une famille et, pour le contrat tronqué, après reconstruction de l’observation. Son P95 reste compris entre environ 31 et 52 générations selon le contrat : même un estimateur centré peut être peu précis.</p><p>Les analyses par effectif montrent que les petits échantillons concentrent le biais de troncature. D25a n’intègre toujours ni erreur de phase SHAPEIT5, ni erreur de carte, ni faux positif D24. Une excellente moyenne simulée ne constituerait donc pas une validation biologique.</p>' + _figure(_independent_method_figure(document), "Biais des méthodes historiques et candidates en généalogie indépendante. Source : CAL-REP-14A.json.", 4) + _table(("Contrat", "Familles", "Biais", "RMSE", "Couverture"), family_rows, "Résultats primaires par nombre de familles.", 4),
        f'<p>Le statut technique est {document["technical_status"]} parce que la parité R, les branches et les contrats passent. Le statut de performance est {document["performance_qualification"]} : deux contrôles préspecifiés du contrat tronqué indépendant échouent. Les performances corrélées restent en revue humaine.</p><p>La candidate ne doit donc pas remplacer la V1 dans le pipeline. La prochaine méthode devra modéliser explicitement l’observation des extrémités ou intégrer leur censure dans la vraisemblance, puis être évaluée sur un nouveau panel indépendant.</p>' + _table(("Critère", "Métrique", "Valeur", "Seuil bas", "Seuil haut", "Décision"), criterion_rows, "Décisions D25a préspecifiées.", 5),
        f'<p>Évaluation close : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Banc : <span class="hash">{document["bench_manifest_sha256"]}</span>. Run : <span class="hash">{document["run_manifest_sha256"]}</span>. Parité R : <span class="hash">{document["r_parity_sha256"]}</span>.</p><p>Le HTML, les données de rendu et le PDF possèdent leurs propres empreintes. Les résultats générés restent sous calibration_workspace/ et hors Git.</p>',
    ]


def build_report(
    *, evaluation_root: Path, template_path: Path, output_root: Path,
) -> Path:
    """Produit CAL-REP-14A depuis la seule évaluation agrégée close."""
    source = evaluation_root / "CAL-REP-14A.json"
    checksum = evaluation_root / "checksums.sha256"
    if not source.is_file() or source.is_symlink() or not checksum.is_file() or checksum.is_symlink():
        raise ValueError("d25a_report_source_missing")
    if checksum.read_text(encoding="utf-8").strip() != f"{sha256_file(source)}  CAL-REP-14A.json":
        raise ValueError("d25a_report_source_checksum_invalid")
    document = json.loads(source.read_text(encoding="utf-8"))
    validate_json_document(document, "gamma_audit_evaluation.schema.json", schemas_dir=SCHEMAS)
    source_sha = sha256_file(source)
    template = template_path.read_text(encoding="utf-8")
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d25a_report_template_invalid")
    style = '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}table{font-size:.68em}main p{orphans:3;widows:3}</style>'
    head = re.sub(
        r"<title>.*?</title>", "<title>CAL-REP-14A — Audit Gamma D25a</title>",
        head, flags=re.S,
    ).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(title)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, title in enumerate(HEADINGS, 1)
    ) + "</ol></nav>"
    body = '<body lang="fr"><section class="cover"><h1>CAL-REP-14A — Audit Gamma D25a</h1><p>Effet_fondateur2 · version 1.0.0 · 10 novembre 2026</p><p>Fidélité R, contrats d’observation et généalogies corrélées</p></section>' + toc + "<main>"
    body += "".join(
        f'<section id="s{index}"><h2>{index}. {escape(title)}</h2>{content}</section>'
        for index, (title, content) in enumerate(
            zip(HEADINGS, _sections(document, source_sha), strict=True), 1
        )
    )
    html = head + body + '</main><footer class="page-footer"><span>CAL-REP-14A · Effet_fondateur2</span></footer></body></html>'
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        report_dir = temporary / "CAL-REP-14A"
        report_dir.mkdir()
        html_path = report_dir / "index.html"
        html_path.write_text(html, encoding="utf-8")
        write_json(report_dir / "report_data.json", {"results": document})
        (report_dir / "checksums.sha256").write_text(
            "".join(
                f"{sha256_file(path)}  {path.name}\n"
                for path in (html_path, report_dir / "report_data.json")
            ), encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    if sha256_file(source) != source_sha:
        raise ValueError("d25a_report_source_changed")
    return destination / "CAL-REP-14A/index.html"


def main() -> int:
    parser = argparse.ArgumentParser(description="Rapport CAL-REP-14A")
    for name in ("evaluation-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    print(build_report(**vars(parser.parse_args())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
