"""Rapport autonome CAL-REP-14 de l’évaluation D25."""

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
    "Question scientifique",
    "Sources et séparation de la vérité",
    "Plan expérimental",
    "Modèle Gamma évalué",
    "Équations et incertitude",
    "Résultats primaires",
    "Effet du nombre de familles",
    "Sensibilités de mesure",
    "Propagation des bornes D24",
    "Critères, anomalies et limites",
    "Conclusion et suite",
    "Provenance et empreintes",
)


def _pct(value: float | None) -> str:
    return "Non évaluable" if value is None else f"{100 * value:.2f} %"


def _num(value: float | None, digits: int = 2) -> str:
    return "Non évaluable" if value is None else f"{value:.{digits}f}"


def _figure(svg: str, caption: str, number: int) -> str:
    return (
        '<figure class="scientific-figure">' + svg
        + f"<figcaption>Figure {number} : {escape(caption)}</figcaption></figure>"
    )


def _table(headers: tuple[str, ...], rows: list[list[str]], caption: str, number: int) -> str:
    head = "".join(f"<th>{escape(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>"
        for row in rows
    )
    return (
        f"<table><caption>Tableau {number} : {escape(caption)}</caption>"
        f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    )


def _age_figure(document: dict[str, Any]) -> str:
    rows = document["results_by_age"]
    parts = [
        '<svg viewBox="0 0 540 280" role="img" aria-label="Âge estimé moyen contre âge vrai">',
        '<text x="8" y="18" font-size="13" font-weight="600">Calibration de l’âge moyen</text>',
        '<line x1="62" y1="235" x2="505" y2="235" stroke="#555"/>',
        '<line x1="62" y1="35" x2="62" y2="235" stroke="#555"/>',
        '<line x1="62" y1="235" x2="505" y2="35" stroke="#999" stroke-dasharray="5 4"/>',
    ]
    for tick in (0, 50, 100, 150):
        x = 62 + 443 * tick / 150
        y = 235 - 200 * tick / 150
        parts.append(f'<text x="{x - 9:.1f}" y="252" font-size="9">{tick}</text>')
        parts.append(f'<text x="34" y="{y + 3:.1f}" font-size="9">{tick}</text>')
    points = []
    for row in rows:
        true_age = row["age_generations"]
        estimate = true_age + row["bias"]
        x = 62 + 443 * true_age / 150
        y = 235 - 200 * max(0, min(150, estimate)) / 150
        points.append((x, y))
        parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="#b24745"/>')
    parts.append('<polyline points="' + " ".join(f"{x:.2f},{y:.2f}" for x, y in points) + '" fill="none" stroke="#b24745" stroke-width="2"/>')
    parts.extend((
        '<text x="230" y="273" font-size="10">Âge vrai (générations)</text>',
        '<text x="8" y="130" font-size="10" transform="rotate(-90 8 130)">Âge estimé moyen</text>',
        '</svg>',
    ))
    return "".join(parts)


def _family_figure(document: dict[str, Any]) -> str:
    rows = document["results_by_family_count"]
    maximum = max(row["rmse"] for row in rows)
    parts = [
        '<svg viewBox="0 0 540 190" role="img" aria-label="Biais et RMSE par nombre de familles">',
        '<text x="8" y="18" font-size="13" font-weight="600">Erreur selon les familles indépendantes</text>',
    ]
    for index, row in enumerate(rows):
        y = 45 + index * 42
        parts.append(f'<text x="8" y="{y + 14}" font-size="11">n={row["family_count"]}</text>')
        bias_width = 370 * abs(row["bias"]) / maximum
        rmse_width = 370 * row["rmse"] / maximum
        parts.append(f'<rect x="65" y="{y}" width="{bias_width:.2f}" height="10" fill="#b24745"/>')
        parts.append(f'<rect x="65" y="{y + 15}" width="{rmse_width:.2f}" height="10" fill="#286d91"/>')
        parts.append(f'<text x="{72 + rmse_width:.2f}" y="{y + 24}" font-size="9">RMSE {row["rmse"]:.1f}</text>')
    parts.extend((
        '<rect x="65" y="174" width="12" height="7" fill="#b24745"/><text x="82" y="181" font-size="9">|biais|</text>',
        '<rect x="140" y="174" width="12" height="7" fill="#286d91"/><text x="157" y="181" font-size="9">RMSE</text>',
        '</svg>',
    ))
    return "".join(parts)


def _coverage_figure(document: dict[str, Any]) -> str:
    rows = document["results_by_family_count"]
    parts = [
        '<svg viewBox="0 0 540 175" role="img" aria-label="Couverture Gamma par nombre de familles">',
        '<text x="8" y="18" font-size="13" font-weight="600">Couverture nominale à 95 %</text>',
        '<line x1="120" y1="40" x2="500" y2="40" stroke="#777"/>',
        '<line x1="481" y1="35" x2="481" y2="155" stroke="#b24745" stroke-dasharray="4 3"/>',
    ]
    for index, row in enumerate(rows):
        y = 70 + index * 36
        low = 120 + 380 * row["coverage_ci95_low"]
        high = 120 + 380 * row["coverage_ci95_high"]
        value = 120 + 380 * row["coverage"]
        parts.append(f'<text x="8" y="{y + 4}" font-size="11">{row["family_count"]} familles</text>')
        parts.append(f'<line x1="{low:.2f}" y1="{y}" x2="{high:.2f}" y2="{y}" stroke="#286d91" stroke-width="3"/>')
        parts.append(f'<circle cx="{value:.2f}" cy="{y}" r="5" fill="#2f7d69"/>')
        parts.append(f'<text x="{min(value + 8, 495):.2f}" y="{y + 4}" font-size="9">{100 * row["coverage"]:.1f}%</text>')
    parts.append('</svg>')
    return "".join(parts)


def _profile_figure(document: dict[str, Any]) -> str:
    rows = [row for row in document["profile_results"] if row["model"] == "CORRELATED"]
    maximum = max(abs(row["bias"]) for row in rows)
    zero_x = 420
    maximum_bar_width = 190
    parts = [
        '<svg viewBox="0 0 540 225" role="img" aria-label="Biais des profils de mesure">',
        '<text x="8" y="18" font-size="13" font-weight="600">Biais corrélé par profil apparié</text>',
        f'<line x1="{zero_x}" y1="35" x2="{zero_x}" y2="195" stroke="#555"/>',
    ]
    for index, row in enumerate(rows):
        y = 50 + index * 38
        value = row["bias"]
        width = maximum_bar_width * abs(value) / maximum
        x = zero_x - width if value < 0 else zero_x
        value_x = x - 40 if value < 0 else zero_x + width + 6
        parts.append(f'<text x="8" y="{y + 12}" font-size="9">{escape(row["profile"])}</text>')
        parts.append(f'<rect x="{x:.2f}" y="{y}" width="{width:.2f}" height="16" fill="#b24745" rx="2"/>')
        parts.append(f'<text x="{value_x:.2f}" y="{y + 12}" font-size="9">{value:.1f}</text>')
    parts.append('<text x="300" y="218" font-size="9">Sous-estimation ← 0 → surestimation</text></svg>')
    return "".join(parts)


def _equations() -> str:
    rows = (
        ("Simulation", "X_f,s ~ Exp(g/100), E[X_f,s] = 100/g cM", "f indexe la famille et s le côté."),
        ("Estimateur indépendant", "ĝ = (2n−1) / Σ_f(L_f+R_f corrigés)", "Les cM sont convertis en Morgans ; la correction de bord est conservée."),
        ("Effectif corrélé", "n_eff = n / [1+(n−1)ρ]", "L’implémentation applique ses branches bornées lorsque ρ est négatif."),
        ("Erreurs", "biais = mean(ĝ−g) ; MAE = mean(|ĝ−g|) ; RMSE = sqrt(mean((ĝ−g)²))", "L’unité est la génération et le scénario est l’unité statistique."),
        ("Couverture", "C = n⁻¹ Σ 1[IC95_inf ≤ g ≤ IC95_sup]", "L’incertitude publiée sur C est l’intervalle de Wilson."),
    )
    return "".join(
        f'<div class="equation"><p><strong>Équation {index} — {escape(title)}</strong></p>'
        f'<div class="math" role="math">{escape(formula)}</div><p>{escape(detail)}</p></div>'
        for index, (title, formula, detail) in enumerate(rows, 1)
    )


def _sections(document: dict[str, Any], source_sha: str) -> list[str]:
    primary = document["primary_metrics"]
    family_rows = [[
        str(row["family_count"]), _num(row["bias"]), _num(row["mae"]),
        _num(row["rmse"]), _pct(row["coverage"]), str(row["estimated_count"]),
    ] for row in document["results_by_family_count"]]
    age_rows = [[
        str(row["age_generations"]), _num(row["bias"]), _num(row["rmse"]),
        _pct(row["coverage"]), str(row["estimated_count"]),
    ] for row in document["results_by_age"]]
    profile_rows = [[
        row["profile"], row["model"], _num(row["bias"]), _num(row["rmse"]),
        _pct(row["coverage"]), str(row["not_evaluable_count"]),
    ] for row in document["profile_results"]]
    d24_rows = [[
        row["profile"], str(row["estimated_count"]), _num(row["bias"]),
        _num(row["rmse"]), _pct(row["coverage"]),
    ] for row in document["d24_propagation_results"]]
    criteria_rows = [[
        row["criterion_id"], row["metric"], _num(row["value"], 4),
        "Non scellé" if row["threshold"] is None else _num(row["threshold"]), row["status"],
    ] for row in document["criteria"]]
    return [
        '<p>CAL-REP-14 évalue la datation de l’étape <code>14_estimate_variant_age</code> dans Effet_fondateur2. Il porte sur une cible synthétique et des longueurs en cM ; aucune donnée clinique, aucune donnée DOCK6 et aucun résultat individuel ne sont publiés.</p><p>La portée primaire est la calibration interne du moteur Gamma sous son modèle de simulation. La propagation des bornes D24 est secondaire. Une datation conditionnelle ne démontre jamais l’existence d’un ancêtre fondateur commun.</p>',
        f'<p>Conformité technique : {document["technical_status"]}. Performance : {document["performance_qualification"]}. Sur {primary["estimated_count"]} scénarios exacts, le biais corrélé est {_num(primary["bias"])} générations, la MAE {_num(primary["mae"])} et la RMSE {_num(primary["rmse"])}. La couverture à 95 % vaut {_pct(primary["coverage"])}.</p><aside class="callout">Le moteur reproduit la référence historique, mais sous-estime systématiquement les âges simulés. Les seuils ne sont pas scellés : calibrated=false ; founder_effect_proven=false.</aside>',
        '<p>D25 demande si l’intervalle Gamma est correctement centré et couvert sur la grille minimale d’âges et d’effectifs familiaux, puis comment la censure, l’échelle de carte et les bornes récupérées modifient cette calibration.</p><p>Le biais, la MAE et la RMSE répondent à des questions différentes. Une couverture proche de 95 % ne suffit pas si l’intervalle est très large ou si des strates se compensent. Les résultats sont donc toujours publiés par âge et par nombre de familles.</p>',
        '<p>La configuration et le protocole sont figés avant résultats. Les 6 000 jeux de longueurs aveugles sont produits séparément de la table des âges vrais. Le runner ferme ensuite les estimations avant que l’évaluateur ouvre cette vérité et ne publie que des agrégats.</p><p>Le rapport ne rouvre ni les longueurs familiales ni les identifiants de cas. Chaque étape est liée par SHA-256. Les sorties D24 sont datées à partir des fichiers réels de l’étape 13, mais leurs cinq âges restent un volet secondaire de faible effectif.</p>',
        '<p>La grille croise dix âges (5 à 150 générations), trois effectifs (3, 5 et 10 familles) et 200 réplications par cellule. Elle fournit 600 réplications indépendantes par âge. Le primaire <code>ORACLE_EXACT</code> utilise les longueurs exactes ; trois profils appariés modifient ensuite uniquement la fenêtre ou l’échelle de carte.</p>' + _table(("Familles", "Biais", "MAE", "RMSE", "Couverture", "n"), family_rows, "Résultats primaires par effectif familial ; source : évaluation close D25.", 1),
        '<p>La méthode évaluée est <code>gamma_gandolfo_2014_v1</code>. Le modèle corrélé est primaire parce que le pipeline l’a fixé ainsi avant D25 ; le modèle indépendant demeure une référence visible et n’est pas remplacé après lecture des résultats.</p><p>La correction de partage fortuit n’est pas appliquée au primaire synthétique. La référence historique l’utilise avec les constantes publiées par le test de six individus. Les intervalles sont conditionnels aux longueurs fournies et n’intègrent pas spontanément l’erreur de phase, de carte ou de sélection du segment.</p>',
        '<p>Le taux exponentiel est exprimé par cM. Le passage en Morgans intervient dans l’estimateur. Les intervalles de biais, MAE et RMSE utilisent 1 000 bootstrap iid ; la couverture utilise Wilson à 95 %. Le leave-one-family-out mesure un changement d’estimation, pas une nouvelle vérité.</p>' + _equations(),
        f'<p>Le biais primaire vaut {_num(primary["bias"])} générations, avec IC bootstrap [{_num(primary["bias_ci95_low"])}, {_num(primary["bias_ci95_high"])}]. L’intervalle exclut zéro. La RMSE atteint {_num(primary["rmse"])} générations. La couverture globale {_pct(primary["coverage"])} possède un IC Wilson [{_pct(primary["coverage_ci95_low"])}, {_pct(primary["coverage_ci95_high"])}].</p><p>La sous-estimation augmente avec l’âge : elle passe de {_num(document["results_by_age"][0]["bias"])} génération à 5 générations vraies à {_num(document["results_by_age"][-1]["bias"])} générations à 150. La figure 1 montre que l’âge moyen s’éloigne progressivement de la diagonale d’identité.</p>' + _figure(_age_figure(document), "Âge estimé moyen par le modèle corrélé contre âge vrai ; diagonale = calibration parfaite. Source : CAL-REP-14.json.", 1) + _table(("Âge vrai", "Biais", "RMSE", "Couverture", "n"), age_rows, "Calibration primaire par âge vrai.", 2),
        '<p>Trois familles donnent la plus forte sous-estimation et la plus faible couverture. Dix familles réduisent l’erreur, mais ne suppriment pas le biais : la couverture de 97,70 % devient alors supérieure au nominal. L’agrégat à 94,12 % masque donc une sous-couverture à trois familles et une sur-couverture à dix.</p><p>Le leave-one-family-out est évaluable sur 4 000 scénarios à 5 ou 10 familles. Le maximum absolu change en moyenne de ' + _num(document["leave_one_family_out"]["mean_max_absolute_delta_generations"]) + ' générations ; son 95e percentile atteint ' + _num(document["leave_one_family_out"]["p95_max_absolute_delta_generations"]) + '. Il est non évaluable à trois familles après retrait d’une unité.</p>' + _figure(_family_figure(document), "Biais absolu et RMSE selon le nombre de familles indépendantes. Source : CAL-REP-14.json.", 2) + _figure(_coverage_figure(document), "Couverture et IC de Wilson par effectif familial ; trait rouge = 95 %. Source : CAL-REP-14.json.", 3),
        '<p>Le plafonnement de chaque bras à 3 cM réduit l’amplitude du biais négatif, mais abaisse la couverture corrélée à 87,98 %. Une erreur d’échelle de carte de −5 % ou +5 % déplace le biais dans les directions attendues, sans résoudre la calibration.</p><p>Le modèle indépendant présente un biais primaire moins négatif que le modèle corrélé, mais une couverture de 89,13 % seulement. Ce contraste interdit de choisir rétrospectivement le modèle indépendant comme nouveau primaire : toute correction devra recevoir un nouvel identifiant et une validation distincte.</p>' + _figure(_profile_figure(document), "Biais du modèle corrélé sous profils appariés de censure et d’échelle de carte. Source : CAL-REP-14.json.", 4) + _table(("Profil", "Modèle", "Biais", "RMSE", "Couverture", "Non évaluables"), profile_rows, "Sensibilités appariées du moteur Gamma.", 3),
        '<p>La reprise des cinq positifs D24 montre l’effet concret des bornes récupérées. MAXIMUM sous-estime en moyenne de 27,61 générations ; BASELINE n’évalue que quatre cas et sous-estime de 79,29 générations ; MISSING_1PCT change même le signe du biais moyen. Avec seulement quatre ou cinq cas, ces valeurs sont descriptives et leurs couvertures ne qualifient rien.</p><p>Ce volet utilise le phasage-vérité D14. Il propage la densité et les erreurs de bornes de l’étape 13, mais pas les erreurs de commutation SHAPEIT5 observées en D23. L’effet de phase reste donc <code>NOT_EVALUABLE</code>.</p>' + _table(("Profil D24", "Estimés", "Biais", "RMSE", "Couverture"), d24_rows, "Propagation secondaire des bornes D24 vers Gamma.", 4),
        '<p>CAL-AGE-006 passe : les estimations historiques corrélée et indépendante valent respectivement ' + _num(document["historical_reference"]["correlated_estimate_generations"], 4) + ' et ' + _num(document["historical_reference"]["independent_estimate_generations"], 4) + ' générations ; l’erreur maximale de référence est ' + _num(document["historical_reference"]["maximum_absolute_reference_error_generations"], 4) + '.</p><p>CAL-AGE-001 à 005 restent REVIEW_REQUIRED faute de seuils scellés. Les limites majeures sont l’hypothèse exponentielle conforme au modèle, l’absence de phase estimée, l’absence de diversité de régions génomiques et le fait que la couverture conditionnelle ne représente pas toute l’incertitude biologique.</p>' + _table(("Critère", "Métrique", "Valeur", "Seuil", "Décision"), criteria_rows, "Raccord aux critères D04.", 5),
        '<p>D25 est techniquement réussi mais sa performance scientifique n’est pas qualifiée. Le biais corrélé, son aggravation avec l’âge, la dépendance au nombre de familles et l’instabilité leave-one-out exigent une revue du modèle avant scellement.</p><p>La suite correcte est de conserver <code>gamma_gandolfo_2014_v1</code> comme référence, de préspecifier une méthode candidate sous un nouvel identifiant, puis de la tester sur des simulations séparées et sur un parcours phasage→bornes→âge apparié. Aucun résultat présent ne justifie <code>calibrated=true</code> ou <code>founder_effect_proven=true</code>.</p>',
        f'<p>Évaluation close : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Banc : <span class="hash">{document["bench_manifest_sha256"]}</span>. Run : <span class="hash">{document["run_manifest_sha256"]}</span>.</p><p>Le HTML, les données de rendu et le PDF reçoivent leurs propres sommes de contrôle. Le rapport consomme seulement CAL-REP-14.json et ne recalcule aucune métrique scientifique.</p>',
    ]


def build_report(*, evaluation_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit CAL-REP-14 depuis la seule évaluation agrégée close."""
    source = evaluation_root / "CAL-REP-14.json"
    checksum = evaluation_root / "checksums.sha256"
    if not source.is_file() or source.is_symlink() or not checksum.is_file() or checksum.is_symlink():
        raise ValueError("d25_report_source_missing")
    if checksum.read_text(encoding="utf-8").strip() != f"{sha256_file(source)}  CAL-REP-14.json":
        raise ValueError("d25_report_source_checksum_invalid")
    document = json.loads(source.read_text(encoding="utf-8"))
    validate_json_document(document, "variant_age_evaluation.schema.json", schemas_dir=SCHEMAS)
    source_sha = sha256_file(source)
    template = template_path.read_text(encoding="utf-8")
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d25_report_template_invalid")
    style = '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}table{font-size:.70em}main p{orphans:3;widows:3}</style>'
    head = re.sub(r"<title>.*?</title>", "<title>CAL-REP-14 — Datation Gamma</title>", head, flags=re.S).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(title)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, title in enumerate(HEADINGS, 1)
    ) + "</ol></nav>"
    body = '<body lang="fr"><section class="cover"><h1>CAL-REP-14 — Datation Gamma</h1><p>Effet_fondateur2 · D25 · version 1.0.0</p><p>Biais, erreur, couverture et sensibilités familiales</p></section>' + toc + "<main>"
    body += "".join(
        f'<section id="s{index}"><h2>{index}. {escape(title)}</h2>{content}</section>'
        for index, (title, content) in enumerate(zip(HEADINGS, _sections(document, source_sha), strict=True), 1)
    )
    html = head + body + '</main><footer class="page-footer"><span>CAL-REP-14 · Effet_fondateur2</span></footer></body></html>'
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        report_dir = temporary / "CAL-REP-14"
        report_dir.mkdir()
        html_path = report_dir / "index.html"
        html_path.write_text(html, encoding="utf-8")
        write_json(report_dir / "report_data.json", {"results": document})
        (report_dir / "checksums.sha256").write_text(
            "".join(f"{sha256_file(path)}  {path.name}\n" for path in (html_path, report_dir / "report_data.json")),
            encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    if sha256_file(source) != source_sha:
        raise ValueError("d25_report_source_changed")
    return destination / "CAL-REP-14/index.html"


def main() -> int:
    parser = argparse.ArgumentParser(description="Rapport CAL-REP-14")
    for name in ("evaluation-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    print(build_report(**vars(parser.parse_args())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
