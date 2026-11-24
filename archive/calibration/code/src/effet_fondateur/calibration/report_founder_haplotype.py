"""Rapport autonome CAL-REP-13 conforme au contrat report-formater."""

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


_SCHEMAS = Path(__file__).resolve().parents[3] / "calibration/schemas"
HEADINGS = (
    "Identification et portée",
    "Résumé décisionnel",
    "Question scientifique",
    "Sources et autorité",
    "Plan expérimental",
    "Vérité indépendante",
    "Méthode évaluée",
    "Équations et incertitude",
    "Résultats quantitatifs",
    "Critères d’acceptation",
    "Anomalies et exclusions",
    "Limites",
    "Conclusion et suite",
    "Provenance et empreintes",
)


def _pct(value: float | None) -> str:
    return "Non évaluable" if value is None else f"{100 * value:.2f} %"


def _num(value: float | None, unit: str = "") -> str:
    return "Non évaluable" if value is None else f"{value:.4g}{unit}"


def _figure(svg: str, caption: str, number: int) -> str:
    return (
        '<figure class="scientific-figure">' + svg
        + f"<figcaption>Figure {number} : {escape(caption)}</figcaption></figure>"
    )


def _confusion_figure(document: dict[str, Any]) -> str:
    confusion = document["confusion"]
    labels = (("Vrai +", "true_positive", "#2f7d69"), ("Faux −", "false_negative", "#c17c2c"),
              ("Vrai −", "true_negative", "#286d91"), ("Faux +", "false_positive", "#b24745"))
    maximum = max(1, *(confusion[key] for _, key, _ in labels))
    bars = []
    for index, (label, key, color) in enumerate(labels):
        y = 38 + index * 35
        width = 360 * confusion[key] / maximum
        bars.append(
            f'<text x="8" y="{y + 14}" font-size="12">{label}</text>'
            f'<rect x="75" y="{y}" width="{width:.2f}" height="19" fill="{color}" rx="2"/>'
            f'<text x="{82 + width:.2f}" y="{y + 14}" font-size="12">{confusion[key]}</text>'
        )
    return '<svg viewBox="0 0 500 185" role="img" aria-label="Matrice de détection agrégée">' \
        '<text x="8" y="18" font-size="13" font-weight="600">Scénarios détectés et non détectés</text>' \
        + ''.join(bars) + '</svg>'


def _rates_figure(document: dict[str, Any]) -> str:
    metrics = document["metrics"]
    rows = (
        ("Sensibilité", metrics["founder_segment_detection_sensitivity"]),
        ("Spécificité", metrics["founder_segment_detection_specificity"]),
        ("Cible contenue", metrics["target_containment_rate_when_detectable"]),
    )
    parts = ['<svg viewBox="0 0 520 165" role="img" aria-label="Taux et intervalles de Wilson">',
             '<text x="8" y="18" font-size="13" font-weight="600">Taux par scénario indépendant</text>',
             '<line x1="160" y1="35" x2="480" y2="35" stroke="#777"/>']
    for tick in range(5):
        x = 160 + tick * 80
        parts.append(f'<text x="{x - 8}" y="31" font-size="9">{tick * 25}%</text>')
    for index, (label, measure) in enumerate(rows):
        y = 62 + index * 38
        parts.append(f'<text x="8" y="{y + 4}" font-size="11">{label}</text>')
        if measure["estimate"] is None:
            parts.append(f'<text x="160" y="{y + 4}" font-size="11">Non évaluable</text>')
            continue
        low, high, estimate = measure["ci95_low"], measure["ci95_high"], measure["estimate"]
        x_low, x_high, x_est = 160 + 320 * low, 160 + 320 * high, 160 + 320 * estimate
        parts.extend((
            f'<line x1="{x_low:.2f}" y1="{y}" x2="{x_high:.2f}" y2="{y}" stroke="#286d91" stroke-width="3"/>',
            f'<circle cx="{x_est:.2f}" cy="{y}" r="5" fill="#b24745"/>',
            f'<text x="{min(475, x_est + 8):.2f}" y="{y + 4}" font-size="10">{100 * estimate:.1f}%</text>',
        ))
    parts.append('</svg>')
    return ''.join(parts)


def _boundary_figure(document: dict[str, Any]) -> str:
    boundary = document["metrics"]["boundary_errors"]
    rows = (
        ("Gauche", boundary["left_mean_absolute_error_cm"], "#286d91"),
        ("Droite", boundary["right_mean_absolute_error_cm"], "#2f7d69"),
    )
    finite = [value for _, value, _ in rows if value is not None]
    maximum = max(finite, default=0.0)
    parts = ['<svg viewBox="0 0 500 130" role="img" aria-label="Erreurs absolues moyennes des bornes">',
             '<text x="8" y="18" font-size="13" font-weight="600">Erreur absolue moyenne des bornes (cM)</text>']
    for index, (label, value, color) in enumerate(rows):
        y = 43 + index * 38
        parts.append(f'<text x="8" y="{y + 15}" font-size="12">{label}</text>')
        if value is None:
            parts.append(f'<text x="90" y="{y + 15}" font-size="11">Non évaluable</text>')
        else:
            width = 330 * value / maximum if maximum else 0
            parts.append(f'<rect x="90" y="{y}" width="{width:.2f}" height="20" fill="{color}" rx="2"/>')
            parts.append(f'<text x="{98 + width:.2f}" y="{y + 15}" font-size="11">{value:.4g} cM</text>')
    parts.append('</svg>')
    return ''.join(parts)


def _profile_figure(document: dict[str, Any]) -> str:
    """Compare les taux par profil sans assimiler les profils à des réplications."""
    rows = document["profile_results"]
    parts = [
        '<svg viewBox="0 0 540 220" role="img" aria-label="Sensibilité et spécificité par profil">',
        '<text x="8" y="18" font-size="13" font-weight="600">Profils appariés : sensibilité et spécificité</text>',
        '<line x1="145" y1="36" x2="505" y2="36" stroke="#777"/>',
    ]
    for tick in range(5):
        x = 145 + tick * 90
        parts.append(f'<text x="{x - 8}" y="32" font-size="9">{tick * 25}%</text>')
    for index, row in enumerate(rows):
        y = 58 + index * 39
        sensitivity = row["sensitivity"]
        specificity = row["specificity"]
        parts.append(f'<text x="8" y="{y + 13}" font-size="10">{escape(row["profile"])}</text>')
        if sensitivity is not None:
            parts.append(
                f'<rect x="145" y="{y}" width="{360 * sensitivity:.2f}" height="10" '
                'fill="#2f7d69" rx="2"/>'
            )
        if specificity is not None:
            parts.append(
                f'<rect x="145" y="{y + 14}" width="{360 * specificity:.2f}" height="10" '
                'fill="#b24745" rx="2"/>'
            )
    parts.extend((
        '<rect x="145" y="205" width="12" height="7" fill="#2f7d69"/>',
        '<text x="162" y="212" font-size="9">Sensibilité</text>',
        '<rect x="235" y="205" width="12" height="7" fill="#b24745"/>',
        '<text x="252" y="212" font-size="9">Spécificité</text>',
        '</svg>',
    ))
    return ''.join(parts)


def _table(headers: tuple[str, ...], rows: list[list[str]], caption: str, number: int) -> str:
    head = ''.join(f'<th>{escape(value)}</th>' for value in headers)
    body = ''.join('<tr>' + ''.join(f'<td>{escape(value)}</td>' for value in row) + '</tr>' for row in rows)
    return f'<table><caption>Tableau {number} : {escape(caption)}</caption><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def _equations() -> str:
    equations = (
        ("Intersection vraie", "L* = max_f(L_f) ; R* = min_f(R_f)", "f indexe les familles indépendantes ; bp et cM sont traités séparément."),
        ("Sensibilité et spécificité", "Se = TP/(TP+FN) ; Sp = TN/(TN+FP)", "L’unité statistique est le scénario, jamais le marqueur."),
        ("Erreur de borne", "e_L = L_prédit−L* ; e_R = R_prédit−R*", "Les erreurs signées et absolues sont publiées pour chaque côté."),
        ("Erreur absolue moyenne", "MAE = Σ_j |e_j| / m", "m est le nombre de bornes vraies non censurées et évaluables ; une borne censurée n’est jamais traitée comme exacte."),
    )
    return ''.join(
        f'<div class="equation"><p><strong>Équation {index} — {escape(title)}</strong></p>'
        f'<div class="math" role="math">{escape(formula)}</div><p>{escape(detail)}</p></div>'
        for index, (title, formula, detail) in enumerate(equations, 1)
    )


def _sections(document: dict[str, Any], source_sha: str) -> list[str]:
    metrics, counts = document["metrics"], document["scenario_counts"]
    sensitivity = metrics["founder_segment_detection_sensitivity"]
    specificity = metrics["founder_segment_detection_specificity"]
    containment = metrics["target_containment_rate_when_detectable"]
    boundary = metrics["boundary_errors"]
    result_rows = [
        ["Sensibilité", str(sensitivity["numerator"]), str(sensitivity["denominator"]), _pct(sensitivity["estimate"])],
        ["Spécificité", str(specificity["numerator"]), str(specificity["denominator"]), _pct(specificity["estimate"])],
        ["Cible contenue", str(containment["numerator"]), str(containment["denominator"]), _pct(containment["estimate"])],
        ["MAE des bornes", "—", str(boundary["evaluated_boundary_count"]), _num(boundary["boundary_mean_absolute_error"], " cM")],
    ]
    criteria_rows = [[row["criterion_id"], row["metric"], "Non mesuré" if row["value"] is None else str(row["value"]), row["status"]] for row in document["criteria"]]
    profile_rows = [
        [
            row["profile"],
            str(row["positive"]),
            str(row["negative"]),
            _pct(row["sensitivity"]),
            _pct(row["specificity"]),
            _num(row["boundary_mean_absolute_error_cm"], " cM"),
        ]
        for row in document["profile_results"]
    ]
    anomalies = '<p>Aucune anomalie agrégée.</p>' if not document["anomalies"] else '<ul>' + ''.join(f'<li>{escape(value)}</li>' for value in document["anomalies"]) + '</ul>'
    return [
        '<p>CAL-REP-13 présente l’évaluation indépendante de l’étape 13, <code>13_infer_founder_haplotype</code>, dans Effet_fondateur2 version 1.0.0. Son objet est limité au candidat haplotypique fondateur autour d’une cible synthétique ; il ne porte ni sur une cohorte clinique ni sur DOCK6.</p><p>Le document est autonome : il décrit le plan, la vérité, les règles de décision, les résultats, l’incertitude et les limites nécessaires à leur lecture. Les données individuelles et la table de vérité restent privées ; seules des statistiques agrégées et leurs empreintes sont publiées.</p>',
        f'<p>Conformité technique : {document["technical_status"]}. Qualification de performance : {document["performance_qualification"]}. Le primaire MAXIMUM comprend {counts["positive"]} scénarios positifs et {counts["negative"]} négatifs ; {counts["profiled_case_total"]} cas appariés ont été contrôlés sur quatre profils.</p><p>Les taux mesurés ne disposent encore d’aucun seuil scientifique scellé. Un résultat techniquement calculable reste donc soumis à revue, même si les contrats, empreintes et dénominateurs sont conformes.</p><aside class="callout">Ce rapport contrôle un banc de développement ; il ne qualifie pas encore la méthode scientifique. calibrated=false ; founder_effect_proven=false.</aside>',
        '<p>La question principale est triple : l’étape 13 détecte-t-elle le segment IBS commun lorsqu’un segment fondateur détectable a été injecté ; évite-t-elle un consensus sur les négatifs ; et, lorsqu’elle fournit des bornes, quelle est leur distance à la vérité ?</p><p>L’inclusion de la cible est évaluée séparément de la simple détection. Cette séparation empêche qu’un segment appelé loin de la cible soit compté comme succès. Enfin, le terme IBS décrit une identité d’état locale : il ne constitue pas une preuve IBD ni une preuve de transmission par un ancêtre commun.</p>',
        '<p>Les autorités normatives sont D24_FOUNDER_HAPLOTYPE_PRESPECIFICATION.md, les critères CAL-FOUND-001 à CAL-FOUND-005 et les contrats de production de l’étape 13. La configuration fixe avant résultats les scénarios, graines, profils, effectifs minimaux et règles de censure.</p><p>L’évaluateur ouvre les sorties closes avant la vérité, contrôle leur schéma et leur cible, puis compare les deux. Chaque cas et chaque document de référence est lié par SHA-256. Aucune donnée <code>dock6_real</code>, aucun cache interne et aucun identifiant individuel ne sont admis dans le rendu.</p>',
        '<p>Le parcours primaire MAXIMUM provient de l’univers de marqueurs à densité maximale exploitable établi en D23b. BASELINE, DENSITY_50 et MISSING_1PCT rejouent exactement les mêmes dix vérités pour décrire la sensibilité à la mesure ; ils ne constituent pas trente réplications biologiques supplémentaires.</p><p>Les cinq positifs couvrent 50, 100, 150 et 200 générations, avec une seconde réplication à 100. Les cinq négatifs placent la cible sur des lignées indépendantes sans copier de segment fondateur commun. Chaque scénario positif contient au moins trois unités familiales indépendantes.</p><p>Les haplotypes soumis à l’étape 13 sont ceux de la vérité phasée D14, ensuite projetés sur chaque profil de marqueurs. Ce choix isole la récupération du segment par l’étape 13 : aucun phasage SHAPEIT5 n’est exécuté dans D24 et les erreurs de commutation mesurées en D23 ne sont pas propagées ici.</p>' + _table(("Profil", "Positifs", "Négatifs", "Sensibilité", "Spécificité", "MAE bornes"), profile_rows, "Effectifs et métriques descriptives par profil apparié ; source : évaluation close D24.", 1),
        '<p>La vérité positive contient, pour chaque famille, les bornes en bp et cM, l’état de censure, la cible et la décision de détectabilité fixée avant les sorties. L’intervalle commun vrai est l’intersection des segments familiaux, non une borne reconstruite depuis l’appel évalué.</p><p>Une borne censurée signifie que l’événement réel se situe au-delà du domaine observable. Elle reste comptée dans l’audit, mais elle n’est pas utilisée comme une position exacte pour calculer une erreur. Les négatifs conservent de vrais porteurs synthétiques de la cible tout en interdisant un segment commun injecté.</p>',
        '<p>La méthode évaluée est le moteur de production <code>target_centered_exact_ibs_v1</code>. Une détection exige simultanément les statuts <code>SUPPORTED_IBS_CANDIDATE</code> et <code>IBS_SHARED_CANDIDATE</code> dans le résumé et le consensus ; une simple présence de fichier ne suffit pas.</p><p>Les paramètres sont contrôlés à trois porteurs indépendants et deux marqueurs flanquants par côté. Les familles absentes ou exclues restent comptées dans les anomalies, sans imputation de borne. La cible, sa position et son identifiant doivent concorder dans le consensus et tous les segments familiaux.</p><p>Le banc compare donc la sortie scientifique réelle de l’étape 13 à une vérité de segments indépendante de cette sortie. Il ne teste toutefois que cette étape conditionnellement à un phasage exact ; il ne constitue pas encore un parcours de bout en bout phasage→segment fondateur.</p>',
        '<p>Les intervalles à 95 % sont ceux de Wilson sur les dix scénarios indépendants du profil primaire. Le containment utilise tous les positifs définis détectables a priori : une absence de segment ne devient pas un succès. Les mesures de sensibilité restent descriptives et ne remplacent pas les seuils qui seront scellés après D34.</p><p>Les erreurs signées indiquent le sens du déplacement ; les erreurs absolues quantifient son amplitude. Les bp et cM sont conservés séparément, car une même distance physique ne représente pas une même distance génétique selon la région. Les bornes censurées sont exclues uniquement du calcul qui exige leur valeur exacte.</p>' + _equations(),
        f'<p>Sur le primaire MAXIMUM, la sensibilité vaut {_pct(sensitivity["estimate"])} sur {sensitivity["denominator"]} positifs détectables et la spécificité {_pct(specificity["estimate"])} sur {specificity["denominator"]} négatifs. Le taux d’inclusion de cible vaut {_pct(containment["estimate"])} ; la MAE des bornes est {_num(boundary["boundary_mean_absolute_error"], " cM")} sur {boundary["evaluated_boundary_count"]} bornes exactes évaluables.</p><p>Le résultat dominant est la présence de quatre faux positifs parmi cinq négatifs au profil MAXIMUM. La sensibilité parfaite sur cinq positifs ne compense pas cette faible capacité à rejeter des porteurs indépendants partageant fortuitement un état allélique autour de la cible. Les intervalles de Wilson restent larges, mais la valeur observée interdit déjà toute qualification favorable.</p><p>La comparaison appariée n’est pas monotone : BASELINE perd un positif mais rejette les cinq négatifs, DENSITY_50 produit deux faux positifs et MISSING_1PCT en produit trois. Ces dix vérités communes décrivent une interaction entre densité et règle IBS ; elles ne permettent pas d’attribuer causalement les différences à la seule densité.</p><p>La figure 1 rend visibles les quatre cellules de détection. La figure 2 associe les taux à leurs intervalles de Wilson. La figure 3 sépare les erreurs gauche et droite. La figure 4 compare les profils sans gonfler les dénominateurs primaires.</p>'
        + _figure(_confusion_figure(document), "Détection agrégée des scénarios positifs détectables et négatifs. Source : CAL-REP-13.json.", 1)
        + _figure(_rates_figure(document), "Sensibilité, spécificité et containment avec intervalles de Wilson à 95 %. Source : CAL-REP-13.json.", 2)
        + _figure(_boundary_figure(document), "Erreur absolue moyenne gauche et droite en cM. Source : CAL-REP-13.json.", 3)
        + _figure(_profile_figure(document), "Sensibilité et spécificité descriptives des quatre profils appliqués aux mêmes vérités. Source : CAL-REP-13.json.", 4)
        + _table(("Mesure", "Numérateur", "Dénominateur", "Résultat"), result_rows, "Résultats principaux sans seuil scientifique scellé.", 2),
        '<p>Les critères <code>pilot_then_seal</code> restent REVIEW_REQUIRED dès qu’ils sont mesurables : D24 fournit les observations nécessaires, mais ne choisit pas rétrospectivement un seuil. Seul CAL-FOUND-005 possède un seuil fixe, égal à zéro.</p><p>Un FAIL de CAL-FOUND-005 signale une revendication IBD non autorisée et rend la conformité technique négative. À l’inverse, un PASS de ce contrôle sémantique ne suffit pas à démontrer une bonne sensibilité, une bonne spécificité ou des bornes exactes.</p>' + _table(("Critère", "Métrique", "Valeur", "Décision"), criteria_rows, "Raccord aux critères approuvés avant pilote.", 3),
        anomalies + '<p>Les trois familles manquantes ou exclues proviennent d’un unique positif BASELINE non détecté : les trois segments familiaux y sont classés <code>INSUFFICIENT</code>. Ce résultat explique la sensibilité BASELINE de 4/5 et ne touche pas le primaire MAXIMUM.</p><p>Les anomalies et les reclassifications IBS→IBD sont conservées comme résultats du banc ; elles ne sont pas supprimées pour améliorer les taux. Seuls des comptes agrégés sont publiés. Les identifiants de famille, de sujet et de scénario restent dans les cas privés.</p>',
        '<p>Les résultats dépendent de la densité, du nombre de familles, de la carte génétique et de la censure. Ils sont conditionnels à un phasage-vérité : D24 n’observe donc pas les conséquences des erreurs de commutation D23b et ne requalifie pas D23. Un banc de bout en bout devra projeter les sorties réelles de phasage sur les mêmes vérités pour mesurer cette propagation.</p><p>Le profil MAXIMUM contient 487 marqueurs dans la région technique, tandis que BASELINE en contient 28 avec la cible protégée. Cette construction renseigne la robustesse locale de l’étape 13, pas toute la densité autosomale de la puce. Le problème scientifique 07 identifié en D21 reste inchangé.</p><p>Le banc utilise une cible technique sur le chromosome 5 et des génomes synthétiques. Il ne valide ni la mutation DOCK6, ni une population humaine, ni une interprétation clinique. Avec cinq scénarios par rôle, les intervalles restent nécessairement larges ; le pilote et la validation aveugle devront fournir l’évidence indépendante.</p>',
        '<p>D24 est techniquement interprétable : les contrats, la matrice appariée, les empreintes et les calculs passent ensemble. Sa performance scientifique n’est pas qualifiée. En particulier, une spécificité primaire observée de 20 % rend la règle actuelle impropre à une validation favorable, même avant le scellement de seuils.</p><p>La prochaine action scientifique n’est donc pas de déclarer D24 réussi, mais d’analyser la cause des appels IBS sur les négatifs et de préspecifier une correction ou un contrôle discriminant. Toute variante devra être évaluée sur des scénarios séparés pour éviter d’optimiser et de valider sur les mêmes dix vérités. D33–D35 restent soumis à ces portes.</p><p>Aucun résultat de ce rapport ne permet de passer <code>calibrated</code> ou <code>founder_effect_proven</code> à vrai.</p>',
        f'<p>Évaluation close : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>.</p><p>Les {len(document["source_case_sha256s"])} empreintes de cas sont conservées dans report_data.json, sans chemin ni identifiant individuel. Le fichier de résultats, le HTML et les données de rendu possèdent leurs propres sommes de contrôle.</p><p>Le moteur de rapport ne rouvre jamais la vérité privée et ne recalcule aucune métrique scientifique. Cette séparation permet de reproduire la présentation à partir de l’évaluation close sans modifier les résultats.</p>',
    ]


def build_report(*, evaluation_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit un rapport français autonome sans rouvrir la vérité privée."""
    source = evaluation_root / "CAL-REP-13.json"
    checksum = evaluation_root / "checksums.sha256"
    if not source.is_file() or source.is_symlink() or not checksum.is_file() or checksum.is_symlink():
        raise ValueError("d24_report_source_missing")
    expected = checksum.read_text(encoding="utf-8").strip()
    if expected != f"{sha256_file(source)}  CAL-REP-13.json":
        raise ValueError("d24_report_source_checksum_invalid")
    document = json.loads(source.read_text(encoding="utf-8"))
    validate_json_document(document, "founder_haplotype_evaluation.schema.json", schemas_dir=_SCHEMAS)
    source_sha = sha256_file(source)
    template = template_path.read_text(encoding="utf-8")
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d24_report_template_invalid")
    style = '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}table{font-size:.72em}main p{orphans:3;widows:3}</style>'
    head = re.sub(r"<title>.*?</title>", "<title>CAL-REP-13 — Haplotype fondateur</title>", head, flags=re.S).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(title)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, title in enumerate(HEADINGS, 1)
    ) + '</ol></nav>'
    body = '<body lang="fr"><section class="cover"><h1>CAL-REP-13 — Haplotype fondateur</h1><p>Effet_fondateur2 · D24 · version 1.0.0</p><p>Comparaison du segment IBS récupéré à la vérité synthétique indépendante</p></section>' + toc
    body += '<main>'
    body += ''.join(
        f'<section id="s{index}"><h2>{index}. {escape(title)}</h2>{content}</section>'
        for index, (title, content) in enumerate(zip(HEADINGS, _sections(document, source_sha), strict=True), 1)
    )
    html = head + body + '</main><footer class="page-footer"><span>CAL-REP-13 · Effet_fondateur2</span></footer></body></html>'
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        report_dir = temporary / "CAL-REP-13"
        report_dir.mkdir()
        html_path = report_dir / "index.html"
        html_path.write_text(html, encoding="utf-8")
        write_json(report_dir / "report_data.json", {"results": document})
        (report_dir / "checksums.sha256").write_text(
            ''.join(f'{sha256_file(path)}  {path.name}\n' for path in (html_path, report_dir / "report_data.json")),
            encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    if sha256_file(source) != source_sha:
        raise ValueError("d24_report_source_changed")
    return destination / "CAL-REP-13/index.html"


def main() -> int:
    """CLI du rendu D24 seul ; aucune statistique n’est recalculée."""
    parser = argparse.ArgumentParser(description="Rapport CAL-REP-13")
    for name in ("evaluation-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    print(build_report(**vars(parser.parse_args())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
