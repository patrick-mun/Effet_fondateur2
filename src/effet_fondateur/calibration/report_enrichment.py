"""Rapport autonome CAL-REP-16B de D28."""

from __future__ import annotations

import argparse
from html import escape
import json
import math
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
    "Identification et domaine de validité", "Résumé décisionnel",
    "Question évaluée", "Séparation appels–vérité",
    "Panel et unités indépendantes", "Moteur réel évalué",
    "Statistique IBS et probabilité empirique", "Probabilités du profil primaire",
    "Calibration sous scénarios nuls", "Puissance sous alternatives",
    "Non-évaluations et densité", "Sensibilités appariées",
    "Critères et décision", "Limites, provenance et empreintes",
)
PROFILE_LABELS = {
    "MAXIMUM": "Maximum", "BASELINE": "Base", "DENSITY_50": "Densité 50 %",
    "MISSING_1PCT": "Manque 1 %",
}


def _pct(value: float) -> str:
    return f"{100 * value:.1f} %"


def _num(value: float | None, digits: int = 4) -> str:
    return "Non évaluable" if value is None else f"{value:.{digits}f}"


def _table(headers: tuple[str, ...], rows: list[list[str]], caption: str, number: int) -> str:
    return (
        f'<table><caption>Tableau {number} : {escape(caption)}</caption><thead><tr>'
        + "".join(f"<th>{escape(value)}</th>" for value in headers)
        + "</tr></thead><tbody>"
        + "".join("<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>" for row in rows)
        + "</tbody></table>"
    )


def _figure(svg: str, caption: str, number: int) -> str:
    return f'<figure class="scientific-figure">{svg}<figcaption>Figure {number} : {escape(caption)}</figcaption></figure>'


def _probability_figure(rows: list[dict[str, Any]], alpha: float) -> str:
    selected = sorted((row for row in rows if row["profile"] == "MAXIMUM"), key=lambda row: (row["truth_role"], row["scenario_label"]))
    parts = ['<svg viewBox="0 0 650 330" role="img" aria-label="Probabilités empiriques du profil primaire">', '<text x="14" y="20" font-size="14" font-weight="600">Probabilité empirique par scénario — profil maximal</text>']
    top, bottom = 45, 260
    for exponent in (0, -1, -2, -3, -4):
        y = top + (0 - exponent) * (bottom - top) / 4
        parts.append(f'<line x1="65" y1="{y:.1f}" x2="630" y2="{y:.1f}" stroke="#d7dde1"/><text x="55" y="{y+4:.1f}" text-anchor="end" font-size="10">10^{exponent}</text>')
    alpha_y = top + (0 - math.log10(alpha)) * (bottom - top) / 4
    parts.append(f'<line x1="65" y1="{alpha_y:.1f}" x2="630" y2="{alpha_y:.1f}" stroke="#c4573a" stroke-width="2" stroke-dasharray="6 4"/><text x="625" y="{alpha_y-5:.1f}" text-anchor="end" font-size="10">α banc = 0,05</text>')
    for index, row in enumerate(selected):
        x = 90 + index * 55
        probability = row["empirical_probability"]
        color = "#286d91" if row["truth_role"] == "POSITIVE" else "#c48b32"
        if probability is None:
            parts.append(f'<path d="M{x-6} 275 L{x+6} 287 M{x+6} 275 L{x-6} 287" stroke="#777" stroke-width="2"/><text x="{x}" y="302" text-anchor="middle" font-size="9">{escape(row["scenario_label"])}</text>')
        else:
            y = top + (0 - max(-4.0, math.log10(float(probability)))) * (bottom - top) / 4
            parts.append(f'<circle cx="{x}" cy="{y:.1f}" r="7" fill="{color}"/><text x="{x}" y="302" text-anchor="middle" font-size="9">{escape(row["scenario_label"])}</text>')
    parts.append('<text x="350" y="324" text-anchor="middle" font-size="11">Scénarios P = positifs ; N = nuls ; croix = non évaluable</text></svg>')
    return "".join(parts)


def _performance_figure(rows: list[dict[str, Any]]) -> str:
    parts = ['<svg viewBox="0 0 650 300" role="img" aria-label="Puissance et faux positifs par profil">', '<text x="14" y="20" font-size="14" font-weight="600">Appels du banc sur cinq scénarios par rôle</text>']
    for index, row in enumerate(rows):
        x = 70 + index * 145
        power_height = 180 * float(row["power"])
        fpr_height = 180 * float(row["null_false_positive_rate"])
        parts.append(f'<rect x="{x}" y="{235-power_height:.1f}" width="42" height="{power_height:.1f}" fill="#286d91"/><rect x="{x+50}" y="{235-fpr_height:.1f}" width="42" height="{fpr_height:.1f}" fill="#c4573a"/><text x="{x+21}" y="{225-power_height:.1f}" text-anchor="middle" font-size="10">{row["positive_calls"]}/5</text><text x="{x+71}" y="{225-fpr_height:.1f}" text-anchor="middle" font-size="10">{row["null_calls"]}/5</text><text x="{x+46}" y="258" text-anchor="middle" font-size="10">{escape(PROFILE_LABELS[row["profile"]])}</text>')
    parts.append('<line x1="55" y1="235" x2="635" y2="235" stroke="#333"/><rect x="220" y="278" width="12" height="12" fill="#286d91"/><text x="238" y="288" font-size="10">positifs appelés</text><rect x="355" y="278" width="12" height="12" fill="#c4573a"/><text x="373" y="288" font-size="10">nuls appelés</text></svg>')
    return "".join(parts)


def _evaluable_figure(rows: list[dict[str, Any]]) -> str:
    parts = ['<svg viewBox="0 0 650 285" role="img" aria-label="Cas évaluables et non évaluables">', '<text x="14" y="20" font-size="14" font-weight="600">Évaluabilité des dix cas par profil</text>']
    for index, row in enumerate(rows):
        evaluable = int(row["evaluable_positive"]) + int(row["evaluable_null"])
        non_evaluable = int(row["non_evaluable_count"])
        x = 65 + index * 145
        eval_height = 18 * evaluable
        non_height = 18 * non_evaluable
        parts.append(f'<rect x="{x}" y="{230-eval_height:.1f}" width="70" height="{eval_height:.1f}" fill="#6b8e23"/><rect x="{x}" y="{230-eval_height-non_height:.1f}" width="70" height="{non_height:.1f}" fill="#b5b8ba"/><text x="{x+35}" y="{220-eval_height:.1f}" text-anchor="middle" font-size="10">{evaluable} évaluables</text><text x="{x+35}" y="250" text-anchor="middle" font-size="10">{escape(PROFILE_LABELS[row["profile"]])}</text>')
    parts.append('<line x1="50" y1="230" x2="630" y2="230" stroke="#333"/><text x="325" y="276" text-anchor="middle" font-size="10">Vert : évaluables · gris : non évaluables</text></svg>')
    return "".join(parts)


def _sharing_figure(rows: list[dict[str, Any]]) -> str:
    selected = [row for row in rows if row["profile"] == "MAXIMUM" and row["empirical_probability"] is not None]
    maximum = max(float(row["observed_total_shared_cm"]) for row in selected)
    parts = ['<svg viewBox="0 0 650 300" role="img" aria-label="Partage IBS et probabilité empirique">', '<text x="14" y="20" font-size="14" font-weight="600">Partage IBS observé et rareté empirique</text>', '<line x1="65" y1="245" x2="625" y2="245" stroke="#333"/><line x1="65" y1="245" x2="65" y2="45" stroke="#333"/>']
    for row in selected:
        x = 65 + 540 * float(row["observed_total_shared_cm"]) / maximum
        y = 45 + 200 * (0 - max(-4.0, math.log10(float(row["empirical_probability"])))) / 4
        color = "#286d91" if row["truth_role"] == "POSITIVE" else "#c48b32"
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{color}"/><text x="{x:.1f}" y="{y-10:.1f}" text-anchor="middle" font-size="9">{escape(row["scenario_label"])}</text>')
    parts.append(f'<text x="345" y="280" text-anchor="middle" font-size="11">Partage total observé (cM), maximum {maximum:.3f}</text><text transform="translate(18 150) rotate(-90)" text-anchor="middle" font-size="11">p empirique (échelle log)</text></svg>')
    return "".join(parts)


def _sections(document: dict[str, Any], source_sha: str) -> list[str]:
    primary = next(row for row in document["profile_results"] if row["profile"] == "MAXIMUM")
    profile_table = [[
        PROFILE_LABELS[row["profile"]], f'{row["evaluable_positive"]}/5', f'{row["evaluable_null"]}/5',
        f'{row["positive_calls"]}/5 ({_pct(row["power"])})', f'{row["null_calls"]}/5 ({_pct(row["null_false_positive_rate"])})',
        _num(row["median_positive_probability"]), _num(row["median_null_probability"]),
    ] for row in document["profile_results"]]
    primary_cases = sorted((row for row in document["case_results"] if row["profile"] == "MAXIMUM"), key=lambda row: (row["truth_role"], row["scenario_label"]))
    case_table = [[
        row["scenario_label"], "Positif" if row["truth_role"] == "POSITIVE" else "Nul",
        "—" if row["age_generations"] is None else str(row["age_generations"]),
        _num(row["observed_total_shared_cm"], 3), _num(row["empirical_probability"], 5),
        "Oui" if row["benchmark_positive_call"] else "Non", row["production_status"],
    ] for row in primary_cases]
    criterion_table = [[row["criterion_id"], row["metric"], _num(float(row["value"]), 3), "0" if row["threshold"] == 0 else "Non scellé", row["status"]] for row in document["criteria"]]
    probability_equation = '<div class="equation"><div class="math">S<sub>obs</sub> = [cM(T)−cM(L)] + [cM(R)−cM(T)]<br>p<sub>emp</sub> = (1+k)/(1+N), avec k = Σ 1(S<sub>b</sub> ≥ S<sub>obs</sub>)</div><p>T est la cible, exclue de la signature ; N est le nombre de tirages nuls évaluables.</p></div>'
    performance_equation = '<div class="equation"><div class="math">FPR = FP / 5 &nbsp;&nbsp; ; &nbsp;&nbsp; puissance = TP / 5</div><p>Une non-évaluation reste au dénominateur et ne devient jamais un appel positif.</p></div>'
    return [
        '<p><strong>Rapport :</strong> CAL-REP-16B · D28. <strong>Domaine :</strong> cœur de l’enrichissement IBS sur chromosome 5 synthétique. <strong>Objet :</strong> calibration sous le nul, puissance descriptive et non-évaluations.</p><aside class="callout">Ce banc évalue des primitives de 16B, pas le parcours externe 1000 Genomes ni l’orchestration bout-en-bout.</aside>',
        f'<p>Le cœur passe ses contrôles techniques : <strong>{document["technical_status"]}</strong>. Sur cinq positifs au profil maximal, {primary["positive_calls"]} sont appelés par la règle du banc ; sur cinq nuls, {primary["null_calls"]} est appelé à tort et {primary["non_evaluable_count"]} n’est pas évaluable.</p><aside class="callout">Puissance {_pct(primary["power"])} (IC95 {_pct(primary["power_ci95_low"])}–{_pct(primary["power_ci95_high"])}), FPR {_pct(primary["null_false_positive_rate"])} (IC95 {_pct(primary["null_fpr_ci95_low"])}–{_pct(primary["null_fpr_ci95_high"])}). Performance : <strong>{document["performance_qualification"]}</strong>.</aside>',
        '<p>D28 demande si une longueur de partage IBS réellement calculée par 16B paraît rare face au fond synthétique, si les scénarios nuls respectent la fréquence attendue d’appels et si les alternatives sont détectées.</p><p>Il ne demande ni de prouver une transmission IBD, ni de conclure à un effet fondateur, ni de choisir après résultat un seuil de production.</p>',
        '<p>Les 40 probabilités ont été calculées sans ouvrir le rôle positif/nul ni le fichier de vérité D24. Le manifeste d’appels a ensuite été scellé par SHA-256. L’évaluation n’a ouvert la vérité qu’après ce scellement.</p><p>Les scénarios sont désignés P1–P5 et N1–N5 dans le rapport ; aucun identifiant individuel n’est publié.</p>',
        '<p>Le primaire comprend dix scénarios biologiques indépendants : cinq positifs à 50, 100, 100, 150 et 200 générations, et cinq nuls sans segment commun injecté. Chacun comporte trois familles indépendantes.</p><p>Les quatre profils sont appariés sur ces mêmes scénarios. Ils ne valent donc que dix unités indépendantes, pas quarante. Les 10 000 tirages nuls par cas mesurent une incertitude Monte-Carlo conditionnelle et ne créent pas de réplications biologiques.</p>',
        '<p>Le banc appelle directement les fonctions de production <code>evaluate_exact_sharing</code>, <code>sample_external_null</code>, <code>summarize_null</code> et la règle de statut. Les représentants familiaux viennent des sorties closes de l’étape 13.</p><p>Le fond est composé des contrôles synthétiques non apparentés D24. bcftools, 1000 Genomes et l’orchestrateur ne sont pas exécutés ; le verdict technique reste donc limité au cœur scientifique.</p>',
        '<p>Le partage exact s’étend de part et d’autre de la cible jusqu’au premier désaccord ou manque. Deux marqueurs partagés sont requis sur chaque flanc. La correction +1 évite une probabilité nulle avec un Monte-Carlo fini.</p>' + probability_equation,
        f'<p>Les cinq scénarios positifs évaluables atteignent le plancher Monte-Carlo 1/10 001 ≈ 0,00010. Parmi les nuls, un scénario passe sous α=0,05, trois restent au-dessus et un est arrêté en amont par l’étape 13.</p><p>La Figure 1 utilise une échelle logarithmique ; une croix représente une non-évaluation, pas une probabilité élevée.</p>' + _figure(_probability_figure(document["case_results"], document["benchmark_decision_alpha"]), "Probabilités empiriques des dix scénarios indépendants du profil maximal ; seuil rouge réservé au banc.", 1) + _table(("Scénario", "Vérité", "Âge", "S obs. cM", "p empirique", "Appel banc", "Statut production"), case_table, "Résultats clos du profil primaire.", 1),
        f'<p>Le taux de faux positifs primaire vaut {primary["null_calls"]}/{primary["null_total"]} = {_pct(primary["null_false_positive_rate"])}. L’IC95 Wilson { _pct(primary["null_fpr_ci95_low"])}–{_pct(primary["null_fpr_ci95_high"])} montre que cinq nuls ne suffisent pas à fixer une calibration précise.</p><p>Le nul non évaluable reste dans le dénominateur selon le protocole. Le taux parmi les seuls évaluables serait descriptif, mais il ne remplace pas le critère primaire.</p>' + performance_equation,
        f'<p>La puissance primaire vaut {primary["positive_calls"]}/{primary["positive_total"]} = {_pct(primary["power"])}. Son IC95 Wilson descend toutefois à {_pct(primary["power_ci95_low"])}. Le résultat est compatible avec une bonne détection sur ces cinq alternatives, mais ne démontre pas une puissance générale.</p><p>Les âges modifient les longueurs partagées, mais tous les positifs proviennent du même schéma à trois familles et du phasage-vérité D24.</p>' + _figure(_sharing_figure(document["case_results"]), "Relation descriptive entre longueur IBS observée et probabilité empirique au profil maximal.", 2),
        '<p>La densité maximale produit neuf cas évaluables sur dix ; le profil avec 1 % de manque en produit huit. À la densité de base et à 50 %, aucun cas n’atteint 10 000 tirages nuls évaluables sous la règle de deux marqueurs par flanc.</p><p>Ces profils ne prouvent donc pas une FPR nulle ni une puissance nulle : ils montrent une perte d’évaluabilité. Les sorties 13 insuffisantes sont distinguées des distributions nulles trop peu évaluables.</p>' + _figure(_evaluable_figure(document["profile_results"]), "Cas évaluables et non évaluables par profil apparié.", 3),
        '<p>La Figure 4 présente les appels bruts du banc, mais doit être lue avec la Figure 3. Les barres nulles des profils de base et à demi-densité résultent de dix non-évaluations, pas de dix décisions négatives fiables.</p><p>Au profil manque 1 %, cinq positifs sur cinq et un nul sur cinq sont appelés, comme au primaire, mais deux nuls sont non évaluables.</p>' + _figure(_performance_figure(document["profile_results"]), "Appels descriptifs positifs et nuls par profil ; chaque rôle compte cinq scénarios.", 4) + _table(("Profil", "Positifs évaluables", "Nuls évaluables", "Appels positifs", "Faux positifs", "p médiane +", "p médiane nulle"), profile_table, "Sensibilités appariées et évaluabilité.", 2),
        '<p>CAL-ENRICH-001 passe : aucune classification positive ou négative de production n’a été émise avec un seuil nul. Les 40 sorties sont <code>NOT_CLASSIFIED</code> lorsqu’elles sont évaluables, ou <code>NOT_EVALUATED</code>.</p><p>CAL-ENRICH-002 et 003 restent <strong>REVIEW_REQUIRED</strong>. Leurs seuils d’acceptation n’étaient pas scellés avant D28 et ne peuvent pas être choisis à partir des valeurs 20 % et 100 %.</p>' + _table(("Critère", "Métrique", "Valeur", "Seuil", "Statut"), criterion_table, "Décision selon les critères préspecifiés.", 3),
        f'<p>Les limites majeures sont n=5 par rôle, trois familles fixes, un fond synthétique partagé, le phasage-vérité et l’absence de la référence externe de production. La forte non-évaluabilité à densité réduite interdit une comparaison naïve des taux d’appels.</p><p>Source agrégée : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Appels scellés : <span class="hash">{document["calls_manifest_sha256"]}</span>.</p><p><strong>Verdicts maintenus :</strong> calibrated=false ; founder_effect_proven=false ; production_promotion=false.</p>',
    ]


def _render(document: dict[str, Any], source_sha: str, template: str) -> str:
    sections = _sections(document, source_sha)
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d28_report_template_invalid")
    style = '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}table{font-size:.65em}main p{orphans:3;widows:3}code{font-size:.9em}</style>'
    head = re.sub(r"<title>.*?</title>", "<title>CAL-REP-16B — Enrichissement haplotypique D28</title>", head, flags=re.S).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + "</ol></nav>"
    body = '<body lang="fr"><section class="cover"><h1>CAL-REP-16B — Calibration de l’enrichissement IBS</h1><p>Effet_fondateur2 · version 1.0.0 · 10 novembre 2026</p><p>D28 · évaluation fermée du cœur de l’étape 16B</p></section>' + toc + '<footer class="page-footer"><span>CAL-REP-16B · Effet_fondateur2</span></footer><main>'
    body += "".join(f'<section><h2 id="s{index}">{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, sections, strict=True), 1))
    toc_sync = '''<script>(function(){
      function sync(){
        document.querySelectorAll('.pagedjs_pages a.toc-entry[href^="#"]').forEach(function(link){
          var target=document.querySelector('.pagedjs_pages '+link.getAttribute('href'));
          var page=target&&target.closest('.pagedjs_page');
          var number=page&&page.getAttribute('data-page-number');
          var cell=link.querySelector('.toc-page-num');
          if(cell&&number&&cell.textContent!==number) cell.textContent=number;
        });
      }
      new MutationObserver(sync).observe(document.documentElement,{childList:true,subtree:true,attributes:true});
      window.addEventListener('load',sync);
    })();</script>'''
    return head + body + '</main>' + toc_sync + '</body></html>'


def build_report(*, evaluation_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit le rapport uniquement depuis le JSON agrégé et empreinté."""
    manifest_path = evaluation_root / "evaluation_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source = evaluation_root / "CAL-REP-16B.json"
    source_sha = sha256_file(source)
    if manifest["artifacts"][0]["sha256"] != source_sha:
        raise ValueError("d28_report_source_checksum_invalid")
    document = json.loads(source.read_text(encoding="utf-8"))
    validate_json_document(document, "enrichment_evaluation.schema.json", schemas_dir=SCHEMAS)
    template = template_path.read_text(encoding="utf-8")
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        report_dir = temporary / "CAL-REP-16B"
        report_dir.mkdir()
        html_path = report_dir / "index.html"
        html_path.write_text(_render(document, source_sha, template), encoding="utf-8")
        data_path = report_dir / "report_data.json"
        write_json(data_path, {"results": document})
        (report_dir / "checksums.sha256").write_text("".join(f"{sha256_file(path)}  {path.name}\n" for path in (html_path, data_path)), encoding="utf-8")
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    if sha256_file(source) != source_sha:
        raise ValueError("d28_report_source_changed")
    return destination / "CAL-REP-16B/index.html"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-root", type=Path, required=True)
    parser.add_argument("--template-path", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    print(build_report(**vars(parser.parse_args())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
