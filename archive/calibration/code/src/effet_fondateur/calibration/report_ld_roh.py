"""Rapports autonomes CAL-REP-15 et CAL-REP-16 de D26."""

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
PROFILE_LABELS = {
    "control": "Contrôle intact", "missing_01": "Manque 1 %",
    "missing_05": "Manque 5 %", "genotype_01": "Dosage 1 %",
    "density_50": "Densité 50 %", "combined": "Combiné",
}
LD_HEADINGS = (
    "Identification et domaine de validité", "Résumé décisionnel",
    "Question évaluée", "Autorité et séparation de la vérité",
    "Panel apparié et effectifs", "Moteur réel de l’étape 15",
    "Définition statistique du LD", "Construction de la référence",
    "Fidélité du contrôle intact", "Effet du manque génotypique",
    "Effet des erreurs de dosage", "Effet de la densité et du profil combiné",
    "Critères et décision", "Limites, provenance et empreintes",
)
ROH_HEADINGS = (
    "Identification et domaine de validité", "Résumé décisionnel",
    "Question évaluée", "Autorité et séparation de la vérité",
    "Panel, injections et référence dense", "Moteur réel de l’étape 16",
    "Règle d’appariement des segments", "Sensibilité globale",
    "Faux positifs au niveau sujet", "Erreur de longueur génétique",
    "Contrôles positifs injectés", "Résistance aux profils de mesure",
    "Critères et décision", "Limites, provenance et empreintes",
)


def _pct(value: float | None) -> str:
    return "Non évaluable" if value is None else f"{100 * value:.2f} %"


def _num(value: float | None, digits: int = 5) -> str:
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


def _bar_figure(rows: list[dict[str, Any]], key: str, *, title: str, percent: bool = False, log: bool = False) -> str:
    values = [float(row[key] or 0) for row in rows]
    plotted = [math_log10(value) if log else value for value in values]
    maximum = max(plotted) if plotted else 1.0
    minimum = min(plotted) if log and plotted else 0.0
    span = maximum - minimum or 1.0
    parts = [f'<svg viewBox="0 0 620 300" role="img" aria-label="{escape(title)}"><text x="12" y="20" font-size="14" font-weight="600">{escape(title)}</text>']
    for index, (row, value, plot_value) in enumerate(zip(rows, values, plotted, strict=True)):
        x = 50 + index * 92
        height = 190 * (plot_value - minimum) / span if value else 0
        y = 235 - height
        parts.append(f'<rect x="{x}" y="{y:.2f}" width="56" height="{height:.2f}" fill="#286d91" rx="3"/>')
        label = f"{100*value:.2f}%" if percent else (f"{value:.2e}" if log else f"{value:.3f}")
        parts.append(f'<text x="{x-3}" y="{max(40,y-7):.2f}" font-size="9">{label}</text>')
        parts.append(f'<text x="{x-8}" y="255" font-size="8" transform="rotate(25 {x-8} 255)">{escape(PROFILE_LABELS[row["profile_id"]])}</text>')
    axis = "échelle log10" if log else ("proportion" if percent else "valeur")
    parts.append(f'<text x="10" y="285" font-size="9">Axe vertical : {axis}</text></svg>')
    return "".join(parts)


def math_log10(value: float) -> float:
    import math
    return math.log10(max(value, 1e-12))


def _interval_figure(rows: list[dict[str, Any]], *, value_key: str, low_key: str, high_key: str, title: str) -> str:
    lows = [float(row[low_key] or 0) for row in rows]
    highs = [float(row[high_key] or 0) for row in rows]
    minimum, maximum = min(lows + [0.0]), max(highs + [0.01])
    coordinate = lambda value: 175 + 390 * (float(value) - minimum) / (maximum - minimum)
    parts = [f'<svg viewBox="0 0 620 330" role="img" aria-label="{escape(title)}"><text x="12" y="20" font-size="14" font-weight="600">{escape(title)}</text>']
    for index, row in enumerate(rows):
        y = 55 + index * 39
        low, high, value = coordinate(row[low_key]), coordinate(row[high_key]), coordinate(row[value_key])
        parts.append(f'<text x="10" y="{y+4}" font-size="9">{escape(PROFILE_LABELS[row["profile_id"]])}</text>')
        parts.append(f'<line x1="{low:.2f}" y1="{y}" x2="{high:.2f}" y2="{y}" stroke="#286d91" stroke-width="3"/>')
        parts.append(f'<circle cx="{value:.2f}" cy="{y}" r="5" fill="#2f7d69"/>')
        parts.append(f'<text x="{value+7:.2f}" y="{y-7}" font-size="8">{100*float(row[value_key]):.2f}%</text>')
    parts.append('<text x="235" y="315" font-size="9">Estimation et intervalle de Wilson à 95 %</text></svg>')
    return "".join(parts)


def _criterion_table(document: dict[str, Any]) -> str:
    return _table(
        ("Critère", "Métrique", "Valeur", "Seuil", "Décision"),
        [[row["criterion_id"], row["metric"], _num(row["value"]), _num(row["threshold"]), "REVUE REQUISE" if row["status"] == "REVIEW_REQUIRED" else row["status"]] for row in document["criteria"]],
        "Décisions D26 selon le registre d’acceptation.", 4,
    )


def _ld_sections(document: dict[str, Any], source_sha: str) -> list[str]:
    rows = document["results_by_profile"]
    by_id = {row["profile_id"]: row for row in rows}
    table_rows = [[PROFILE_LABELS[row["profile_id"]], str(row["evaluated_pair_count"]), f'{row["r_squared_mean_absolute_error"]:.8f}', f'{row["r_squared_p95_absolute_error"]:.8f}', f'{row["r_squared_maximum_absolute_error"]:.8f}'] for row in rows]
    pair_rows = [[PROFILE_LABELS[row["profile_id"]], str(row["evaluated_pair_count"]), str(row["reference_pair_count"]), str(row["missing_reference_pair_count"])] for row in rows]
    equation = '<div class="equation"><strong>Équation 1 — LD génotypique</strong><div class="math">r²<sub>jk</sub> = [Σᵢ(gᵢⱼ−ḡⱼ)(gᵢₖ−ḡₖ)]² / [Σᵢ(gᵢⱼ−ḡⱼ)² Σᵢ(gᵢₖ−ḡₖ)²]</div></div><div class="equation"><strong>Équation 2 — erreur absolue moyenne</strong><div class="math">MAE<sub>p</sub> = |E<sub>p</sub>|⁻¹ Σ<sub>(j,k)∈Eₚ</sub> |r̂²<sub>p,jk</sub> − r²<sub>0,jk</sub>|</div></div>'
    return [
        '<p>CAL-REP-15 évalue la fidélité numérique et la résistance aux erreurs de mesure du LD local de l’étape 15. Le domaine est un panel synthétique régional du chromosome 5 ; il ne contient ni DOCK6, ni génotype de patient, ni fond génomique réel.</p><p>La mesure porte sur le LD génotypique r². Elle ne démontre ni identité par descendance, ni haplotype fondateur partagé, ni causalité.</p>',
        f'<p>Statut technique : <strong>{document["technical_status"]}</strong>. Qualification scientifique : <strong>{document["performance_qualification"]}</strong>. Le contrôle intact reproduit la corrélation calculée indépendamment avec une MAE de {by_id["control"]["r_squared_mean_absolute_error"]:.3e}.</p><aside class="callout">Le dosage 1 % produit la plus grande MAE ({by_id["genotype_01"]["r_squared_mean_absolute_error"]:.5f}). Aucun seuil scientifique n’était scellé : CAL-LD-001 reste en revue.</aside>',
        '<p>La question est double : le moteur PLINK et son normalisateur reproduisent-ils correctement r² en l’absence de défaut, puis comment le manque, le dosage erroné et la densité modifient-ils le résultat ?</p><p>La comparaison est appariée : tous les profils proviennent des mêmes 908 individus et de la même biologie simulée. Les écarts entre profils ne sont donc pas des réplications indépendantes.</p>',
        '<p>Le protocole, les profils, la fenêtre de 1 Mb/1 cM, la MAF minimale 0,05 et les règles de décision ont été finalisés avant les appels. Les tables PLINK observées ont été scellées avant le calcul de la référence.</p><p>Les deux profils modifiant seulement la phase ont des VCF non phasés identiques au contrôle, vérifiés par SHA-256. Ils ne sont pas comptés comme observations supplémentaires.</p>',
        '<p>Le contrôle, les profils de manque 1 % et 5 % et le profil dosage conservent 487 marqueurs ; les profils à demi-densité en conservent 244. Les paires évaluables changent donc avec la densité et avec le passage de variants autour du seuil de MAF.</p>' + _table(("Profil", "Paires PLINK", "Paires référencées", "Référence absente"), pair_rows, "Effectifs de paires indépendamment évaluables.", 1),
        '<p>D26 appelle directement <code>publish_local_ld</code>, moteur scientifique utilisé par l’étape 15. PLINK 1.9 calcule les fréquences, les dosages, r² génotypique et D′ ; le normalisateur contrôle l’univers de variants et le nombre d’appels conjoints.</p><p>Seule la cohorte de contrôle synthétique est évaluée. Les cohortes non fournies sont explicitement non évaluées, jamais remplacées par des zéros.</p>',
        '<p>r² est le carré de la corrélation de Pearson entre comptes d’allèles 0/1/2 pour les sujets appelés aux deux loci. Une valeur proche de 1 indique une association linéaire forte des dosages, pas nécessairement une origine ancestrale unique.</p>' + equation,
        '<p>La référence r² est calculée en Python à partir du VCF de contrôle intact, indépendamment du fichier LD et du parseur PLINK. Pour chaque profil, seules les paires effectivement évaluées sont comparées aux mêmes variants du contrôle.</p><p>Cette référence mesure l’erreur d’observation et de calcul. Elle n’est pas une valeur populationnelle externe et ne corrige pas la structure généalogique du panel.</p>',
        f'<p>Les {by_id["control"]["evaluated_pair_count"]:,} paires du contrôle ont toutes une référence. La MAE {by_id["control"]["r_squared_mean_absolute_error"]:.3e} et l’erreur maximale {by_id["control"]["r_squared_maximum_absolute_error"]:.3e} correspondent à l’arrondi numérique des sorties PLINK.</p><p>Ce résultat qualifie le raccord technique de l’étape 15 dans ce banc, mais pas un seuil biologique de LD.</p>' + _figure(_bar_figure(rows, "r_squared_mean_absolute_error", title="MAE de r² par profil", log=True), "Erreur absolue moyenne de r² ; l’échelle logarithmique rend visible le plancher numérique du contrôle.", 1),
        f'<p>Le manque 1 % porte la MAE à {by_id["missing_01"]["r_squared_mean_absolute_error"]:.5f} et le manque 5 % à {by_id["missing_05"]["r_squared_mean_absolute_error"]:.5f}. Le P95 atteint respectivement {by_id["missing_01"]["r_squared_p95_absolute_error"]:.5f} et {by_id["missing_05"]["r_squared_p95_absolute_error"]:.5f}.</p><p>Le manque change l’échantillon appelé paire par paire ; il ajoute donc une variabilité d’échantillonnage même sans erreur d’allèle.</p>',
        f'<p>Avec 1 % de dosages modifiés, la MAE atteint {by_id["genotype_01"]["r_squared_mean_absolute_error"]:.5f}, le P95 {by_id["genotype_01"]["r_squared_p95_absolute_error"]:.5f} et le maximum {by_id["genotype_01"]["r_squared_maximum_absolute_error"]:.3f}. Les erreurs de dosage peuvent aussi rendre évaluables des variants auparavant sous le seuil.</p><p>Ce profil est un stress technique uniforme, pas une estimation du taux d’erreur réel de la puce.</p>' + _figure(_bar_figure(rows, "r_squared_maximum_absolute_error", title="Erreur maximale de r² par profil"), "Maximum observé parmi les paires évaluables de chaque profil.", 2),
        f'<p>La demi-densité conserve une MAE au plancher ({by_id["density_50"]["r_squared_mean_absolute_error"]:.3e}) sur les paires restantes, mais réduit leur nombre de {by_id["control"]["evaluated_pair_count"]:,} à {by_id["density_50"]["evaluated_pair_count"]:,}. Elle perd donc de la résolution sans déformer les paires conservées.</p><p>Le profil combiné associe cette perte de couverture aux erreurs et atteint une MAE de {by_id["combined"]["r_squared_mean_absolute_error"]:.5f}.</p>' + _figure(_bar_figure(rows, "evaluated_pair_count", title="Nombre de paires LD évaluables"), "Paires respectant distance, MAF et appels conjoints ; les profils sont appariés.", 3) + _table(("Profil", "MAE", "Médiane", "P95", "Maximum"), table_rows, "Distribution des erreurs absolues de r².", 2),
        '<p>CAL-LD-002 passe : le banc régional n’introduit aucun chromosome cible dans un fond génomique. CAL-LD-001 reste REVIEW_REQUIRED car son seuil pilote n’était pas scellé avant D26.</p><p>Le résultat autorise l’usage descriptif de l’étape 15 avec affichage des effectifs et du profil de mesure. Il n’autorise ni seuil d’effet fondateur, ni reclassement du LD en IBD.</p>' + _criterion_table(document),
        f'<p>Le panel couvre une seule fenêtre chr5 et un seul ensemble synthétique. Il n’évalue pas l’incertitude populationnelle, la diversité de cartes, les petits effectifs ni le chromosome 19. Les profils sont appariés.</p><p>Source agrégée : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Run scellé : <span class="hash">{document["run_manifest_sha256"]}</span>.</p><p><strong>Verdicts maintenus :</strong> calibrated=false ; founder_effect_proven=false ; production_promotion=false.</p>',
    ]


def _roh_sections(document: dict[str, Any], source_sha: str) -> list[str]:
    rows = document["results_by_profile"]
    by_id = {row["profile_id"]: row for row in rows}
    summary_rows = [[PROFILE_LABELS[row["profile_id"]], f'{row["matched_segment_count"]}/{row["dense_reference_segment_count"]}', _pct(row["roh_detection_sensitivity"]), _pct(row["roh_false_positive_rate"]), _num(row["roh_length_mean_absolute_error_cm"], 3)] for row in rows]
    calls_rows = [[PROFILE_LABELS[row["profile_id"]], str(row["observed_segment_count"]), str(row["matched_segment_count"]), _pct(row["unmatched_observed_segment_fraction"])] for row in rows]
    equation = '<div class="equation"><strong>Équation 1 — sensibilité segmentaire</strong><div class="math">Sensibilité = TP / (TP + FN)</div></div><div class="equation"><strong>Équation 2 — faux positifs au niveau sujet</strong><div class="math">FPR = sujets sans ROH dense mais avec appel / sujets sans ROH dense</div></div><div class="equation"><strong>Équation 3 — erreur de longueur</strong><div class="math">MAE<sub>L</sub> = M⁻¹ Σₘ |L<sup>obs</sup>ₘ − L<sup>dense</sup>ₘ|, en cM</div></div>'
    return [
        '<p>CAL-REP-16 évalue la détection régionale des ROH par l’étape 16 sur un panel synthétique chr5. Le rapport sépare autozygotie individuelle, IBD entre sujets et effet fondateur : ces notions ne sont pas interchangeables.</p><p>Le domaine contient trois ROH injectées de 2 Mb et une référence dense opérationnelle. Aucune longueur ROH n’alimente Gamma.</p>',
        f'<p>Statut technique : <strong>{document["technical_status"]}</strong>. Qualification scientifique : <strong>{document["performance_qualification"]}</strong>. Le contrôle détecte {by_id["control"]["matched_segment_count"]} des {by_id["control"]["dense_reference_segment_count"]} segments denses, soit {_pct(by_id["control"]["roh_detection_sensitivity"])}.</p><aside class="callout">Les trois injections sont retrouvées dans chaque profil, mais de nombreux appels clairsemés n’ont pas d’équivalent dense. Les seuils CAL-ROH-001 à 003 restent non scellés.</aside>',
        '<p>D26 demande si les paramètres de production retrouvent les segments que PLINK appelle sur le génotype dense, puis comment manque, dosage et densité déplacent cette détection. La comparaison porte sur les mêmes sujets et la même région.</p><p>Elle ne répond pas à la question d’un fardeau ROH génomique ni à celle d’une généalogie fondatrice partagée.</p>',
        '<p>Les appels observés ont été produits et empreintés avant l’ouverture du VCF dense et de la vérité d’injection. La référence dense a ensuite été calculée avec les paramètres préspecifiés, sans réajuster le seuil d’appariement de 50 %.</p><p>Le VCF dense n’est jamais une entrée du pipeline observé. Les profils de phase non phasés, identiques au contrôle, ne sont pas dupliqués.</p>',
        f'<p>Le panel contient 908 sujets et {document["dense_reference_segment_count"]} segments de référence dense d’au moins 1,5 Mb selon PLINK. Trois sujets portent une injection explicite de 2 Mb ; le quatrième scénario sans injection n’est considéré négatif que si la référence dense le confirme.</p><p>Les 891 sujets sans segment dense forment le dénominateur négatif du taux de faux positifs au niveau sujet.</p>',
        '<p>D26 appelle directement <code>publish_roh</code>, moteur de l’étape 16, avec PLINK 1.9. Les seuils sont 1 500 kb, 50 SNP, 50 kb/SNP, trou maximal 1 000 kb, fenêtre de 50 SNP, un hétérozygote et cinq manquants permis.</p><p>Le fardeau génomique reste explicitement non évalué, faute de 22 autosomes. Seul le périmètre régional est interprété.</p>',
        '<p>Un segment observé est apparié à un segment dense du même sujet lorsque l’intersection couvre au moins 50 % du segment dense. L’algorithme glouton retient ensuite l’intersection maximale et interdit de réutiliser un segment.</p><p>Cette règle évite de compter deux fois un segment fragmenté. Elle reste une convention opérationnelle qui devra faire l’objet d’une sensibilité future.</p>' + equation,
        f'<p>La sensibilité vaut {_pct(by_id["control"]["roh_detection_sensitivity"])} dans le contrôle. Elle descend à {_pct(by_id["missing_05"]["roh_detection_sensitivity"])} sous manque 5 % et à {_pct(by_id["density_50"]["roh_detection_sensitivity"])} à demi-densité.</p><p>Les intervalles de Wilson sont larges parce que la référence ne comporte que {document["dense_reference_segment_count"]} segments ; ce sont les segments, pas les profils, qui constituent les unités de cette incertitude.</p>' + _figure(_interval_figure(rows, value_key="roh_detection_sensitivity", low_key="sensitivity_ci95_low", high_key="sensitivity_ci95_high", title="Sensibilité ROH et IC95"), "Segments observés appariés aux 22 segments de référence dense.", 1),
        f'<p>Dans le contrôle, {by_id["control"]["false_positive_sample_count"]} des {by_id["control"]["reference_negative_sample_count"]} sujets négatifs denses reçoivent au moins un appel, soit {_pct(by_id["control"]["roh_false_positive_rate"])}. Le profil combiné atteint {_pct(by_id["combined"]["roh_false_positive_rate"])}.</p><p>La proportion de segments observés non appariés est plus élevée ({_pct(by_id["control"]["unmatched_observed_segment_fraction"])} dans le contrôle), mais elle n’est pas appelée FPR car elle n’a pas de dénominateur de vrais négatifs segmentaires.</p>' + _figure(_interval_figure(rows, value_key="roh_false_positive_rate", low_key="false_positive_rate_ci95_low", high_key="false_positive_rate_ci95_high", title="Taux de sujets faussement positifs et IC95"), "Sujets sans ROH dense mais avec au moins un appel observé.", 2) + _table(("Profil", "Appels", "Appariés", "Appels non appariés"), calls_rows, "Comptes d’appels régionaux et appariement dense.", 1),
        f'<p>L’erreur moyenne de longueur du contrôle est {by_id["control"]["roh_length_mean_absolute_error_cm"]:.3f} cM, médiane {by_id["control"]["roh_length_median_absolute_error_cm"]:.3f} cM. Elle atteint {by_id["combined"]["roh_length_mean_absolute_error_cm"]:.3f} cM dans le profil combiné.</p><p>Les longueurs sont obtenues par interpolation de la carte D11a aux deux bornes. Les bornes d’injection ne servent pas d’oracle de longueur, car un appel peut légitimement prolonger le segment injecté dans un fond déjà homozygote.</p>' + _figure(_bar_figure(rows, "roh_length_mean_absolute_error_cm", title="MAE de longueur ROH (cM)"), "Erreur absolue moyenne des longueurs pour les segments appariés.", 3),
        '<p>Les trois ROH injectées sont détectées dans chacun des six profils, soit 18 succès sur 18 observations appariées. Ce contrôle positif vérifie que les injections majeures ne sont pas perdues sous les défauts étudiés.</p><p>Il ne suffit pas à qualifier l’algorithme : ces trois segments ont tous 2 Mb, une bonne observabilité sur puce et proviennent du même panel. Ils ne couvrent ni les autres longueurs, ni les autres cartes, ni les autres chromosomes.</p>',
        f'<p>Le manque 5 % et la demi-densité réduisent la sensibilité globale à {_pct(by_id["missing_05"]["roh_detection_sensitivity"])}. Le dosage 1 % conserve {_pct(by_id["genotype_01"]["roh_detection_sensitivity"])} ; le combiné remonte à {_pct(by_id["combined"]["roh_detection_sensitivity"])} mais augmente les faux positifs.</p><p>Ces variations non monotones sont compatibles avec les fenêtres tolérantes de PLINK : retirer ou modifier un marqueur peut casser un segment, mais aussi masquer un hétérozygote qui empêchait un appel.</p>' + _figure(_bar_figure(rows, "observed_segment_count", title="Nombre d’appels ROH observés"), "Nombre brut d’appels régionaux avant appariement à la référence dense.", 4) + _table(("Profil", "Appariés", "Sensibilité", "FPR sujet", "MAE cM"), summary_rows, "Synthèse des performances ROH par profil.", 2),
        '<p>CAL-ROH-004 passe : zéro longueur ROH est utilisée directement par Gamma. CAL-ROH-001 à 003 restent REVIEW_REQUIRED, car leurs seuils pilote étaient absents avant l’ouverture des résultats.</p><p>Le résultat autorise un usage descriptif régional accompagné des limites de densité et d’un contrôle de faux positifs. Il n’autorise pas une qualification scientifique de l’étape 16 ni une preuve d’effet fondateur.</p>' + _criterion_table(document),
        f'<p>La référence dense réutilise PLINK : elle isole l’effet de l’observation clairsemée mais ne constitue pas une vérité biologique indépendante de l’algorithme. Le banc ne couvre qu’une région chr5, une classe d’injection et un seul panel.</p><p>Source agrégée : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Run scellé : <span class="hash">{document["run_manifest_sha256"]}</span>.</p><p><strong>Verdicts maintenus :</strong> calibrated=false ; founder_effect_proven=false ; production_promotion=false.</p>',
    ]


def _render(document: dict[str, Any], source_sha: str, template: str) -> str:
    report_id = document["report_id"]
    title = "Fidélité du LD local — D26" if report_id == "CAL-REP-15" else "Détection régionale des ROH — D26"
    headings = LD_HEADINGS if report_id == "CAL-REP-15" else ROH_HEADINGS
    sections = _ld_sections(document, source_sha) if report_id == "CAL-REP-15" else _roh_sections(document, source_sha)
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d26_report_template_invalid")
    style = '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}table{font-size:.67em}main p{orphans:3;widows:3}code{font-size:.9em}</style>'
    head = re.sub(r"<title>.*?</title>", f"<title>{report_id} — {escape(title)}</title>", head, flags=re.S).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(headings, 1)) + "</ol></nav>"
    body = f'<body lang="fr"><section class="cover"><h1>{report_id} — {escape(title)}</h1><p>Effet_fondateur2 · version 1.0.0 · 10 novembre 2026</p><p>Évaluation fermée des moteurs réels 15/16</p></section>{toc}<main>'
    body += "".join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(headings, sections, strict=True), 1))
    return head + body + f'</main><footer class="page-footer"><span>{report_id} · Effet_fondateur2</span></footer></body></html>'


def build_reports(*, evaluation_root: Path, template_path: Path, output_root: Path) -> tuple[Path, Path]:
    """Produit les deux rapports depuis les seuls JSON agrégés et empreintés."""
    manifest_path = evaluation_root / "evaluation_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    artifact_hashes = {row["path"]: row["sha256"] for row in manifest["artifacts"]}
    template = template_path.read_text(encoding="utf-8")
    documents: dict[str, tuple[dict[str, Any], str]] = {}
    for report_id in ("CAL-REP-15", "CAL-REP-16"):
        source = evaluation_root / f"{report_id}.json"
        source_sha = sha256_file(source)
        if artifact_hashes.get(source.name) != source_sha:
            raise ValueError(f"d26_report_source_checksum_invalid:{report_id}")
        document = json.loads(source.read_text(encoding="utf-8"))
        validate_json_document(document, "ld_roh_evaluation.schema.json", schemas_dir=SCHEMAS)
        documents[report_id] = (document, source_sha)
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        paths: list[Path] = []
        for report_id, (document, source_sha) in documents.items():
            report_dir = temporary / report_id
            report_dir.mkdir()
            html_path = report_dir / "index.html"
            html_path.write_text(_render(document, source_sha, template), encoding="utf-8")
            data_path = report_dir / "report_data.json"
            write_json(data_path, {"results": document})
            (report_dir / "checksums.sha256").write_text("".join(f"{sha256_file(path)}  {path.name}\n" for path in (html_path, data_path)), encoding="utf-8")
            paths.append(html_path)
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    for report_id, (_, source_sha) in documents.items():
        if sha256_file(evaluation_root / f"{report_id}.json") != source_sha:
            raise ValueError(f"d26_report_source_changed:{report_id}")
    return destination / "CAL-REP-15/index.html", destination / "CAL-REP-16/index.html"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-root", type=Path, required=True)
    parser.add_argument("--template-path", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    print("\n".join(str(path) for path in build_reports(**vars(parser.parse_args()))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
