"""Quatre rapports D23 : sources closes, graphiques et méthodes explicites."""
from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.evaluate_phasing import STAGES, checked_closure, validate_evaluation
from effet_fondateur.calibration.phasing_report_graphics import count_bars, genetic_map_curve, rate_intervals
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.report_phasing_preparation import column_chart, figure
from effet_fondateur.calibration.report_technical_evaluation import HEADINGS, _table
from effet_fondateur.calibration.technical_scenarios import new_workspace_run


LABELS = {"10A": "Préparation autosomale", "11": "Région cible et carte génétique", "12": "Phasage régional et cible", "12A": "Phasage des 22 autosomes"}


def equations(stage: str) -> str:
    """Représente exclusivement les règles effectivement utilisées par les oracles."""
    if stage == "10A":
        expressions = [("Sélection des marqueurs", "K = {j : tous les G<sub>ij</sub> sont observés et MAF<sub>j</sub> ≥ 10<sup>−6</sup>}", "G est le génotype diploïde observé ; complétude sur les 40 sujets, MAF minoritaire sur les 32 fondateurs, conformément à PLINK sans --nonfounders. Le filtre exclut tout locus incomplet ou monomorphe. Il ne s’agit pas d’une estimation de fréquence dans une population humaine."),
                       ("Concordance exhaustive", "D = ∑<sub>i,j∈K</sub> 𝟙[G<sub>ij,sortie</sub> ≠ G<sub>ij,attendu</sub>]", "Comparer les bases nucléotidiques, sans supposer A1=ALT. Exiger D=0 et les mêmes ensembles et ordres de sujets/loci. Aucun intervalle d’échantillonnage : toutes les cellules de ce banc sont comparées.")]
    elif stage == "11":
        expressions = [("Fenêtre et interpolation", "W = {j : chr<sub>j</sub> = chr<sub>cible</sub>, |bp<sub>j</sub> − bp<sub>cible</sub>| ≤ 5×10<sup>6</sup>}<br>cM(bp) = cM<sub>L</sub> + (cM<sub>R</sub>−cM<sub>L</sub>)(bp−bp<sub>L</sub>)/(bp<sub>R</sub>−bp<sub>L</sub>)", "bp : bases ; cM : centimorgans ; L/R : points encadrants de la carte officielle. Pas d’extrapolation ; un taux bp/cM constant n’est pas supposé."),
                       ("Contrôle numérique", "E = ∑<sub>j∈W</sub> 𝟙[|cM<sub>sortie,j</sub> − cM<sub>oracle,j</sub>| &gt; 10<sup>−8</sup>]", "La tolérance est numérique, pas un seuil de précision biologique. Exiger E=0, cible unique, ordre et génotypes concordants. Les comptes exhaustifs ne requièrent pas d’IC d’échantillonnage.")]
    else:
        expressions = [("Commutations et phase ponctuelle", "o<sub>j</sub> = 𝟙[H<sub>prédit,j,1</sub> ≠ H<sub>vrai,j,1</sub>]<br>SER = ∑ 𝟙[o<sub>j</sub> ≠ o<sub>j−1</sub>] / N<sub>transitions</sub><br>Hamming = ∑ 𝟙[o<sub>j</sub> ≠ b] / N<sub>hets ancrés</sub><br>A<sub>cible</sub> = ∑ 𝟙[o<sub>cible</sub> = b] / N<sub>cibles hétérozygotes ancrables</sub>", "o vaut 0/1 aux hétérozygotes concordants et phasés. b est l’orientation majoritaire hors cible, par sujet/chromosome ou bloc PS. Une inversion globale des deux copies ne crée pas d’erreur. Homozygotes ignorés ; manque, discordance et rupture PS interrompent les transitions. Égalité d’ancrage : phase ponctuelle non évaluable ; A est la précision d’orientation cible, applicable seulement à 12."),
                       ("Phase rare et incertitude familiale", "MAF<sub>j</sub> = min(p<sub>j</sub>, 1−p<sub>j</sub>), p<sub>j</sub> = ∑<sub>32 fondateurs,2 copies</sub> H / 64<br>R<sup>*</sup> = ∑<sub>u tirés avec remise</sub> erreurs<sub>u</sub> / ∑<sub>u tirés avec remise</sub> opportunités<sub>u</sub><br>IC<sub>95 %</sub> = [quantile<sub>0,025</sub>(R<sup>*</sup>), quantile<sub>0,975</sub>(R<sup>*</sup>)]", "La phase rare porte sur MAF ≤5 %. Bootstrap de 24 unités : 8 trios entiers et 16 sujets isolés ; 2 000 tirages, graine 2305. Les chromosomes, marqueurs et trois profils ne sont pas des sujets indépendants. IC descriptif conditionnel à ce petit banc, sans garantie sur l’étude réelle.")]
    return ''.join(f'<div class="equation"><p><strong>Équation {index} — {escape(title)}</strong></p><div class="math" role="math">{formula}</div><p>{description}</p></div>' for index, (title, formula, description) in enumerate(expressions, 1))


def graphics(stage: str, rows: list[dict[str, Any]]) -> str:
    """Graphiques quantitatifs, sans recalcul de métriques ni de décision."""
    available = [row for row in rows if "observed_markers" in row or "metrics" in row]
    if not available:
        return '<p>Aucune sortie quantitative disponible : les processus ont échoué. Aucun graphique de performance inventé.</p>'
    if stage == "10A":
        plots = [(count_bars([row["profile"] for row in available], [row["observed_markers"] for row in available], "Combien de marqueurs sont retenus par 10A ?", "Nombre de SNP"), "Effectifs réellement retenus, profils appariés. Source : 10A.evaluation.json."),
                 (column_chart([dict(chromosome=row["chromosome"], markers=row["observed"]) for row in available[0]["chromosome_counts"]], "markers", "Couverture des 22 autosomes — cas de base", "Nombre de marqueurs retenus"), "Couverture autosomale, pas indépendance des marqueurs. Source : oracle 10A, cas de base.")]
    elif stage == "11":
        plots = [(genetic_map_curve(available[0]["map_plot"]), "Carte et interpolation indépendante : les courbes peuvent se superposer exactement. Source : 11.evaluation.json, map_plot ; 40 sujets par profil."),
                 (count_bars([row["profile"] for row in available], [row["observed_markers"] for row in available], "Quelle densité reste dans la fenêtre cible ?", "Nombre de marqueurs, fenêtre ±5 Mb"), "Effectifs régionaux après le profil de densité ; pas un indicateur de qualité du phasage. Source : oracle 11.")]
    else:
        plots = [(rate_intervals(rows, "switch_error_rate", "Commutations entre hétérozygotes successifs"), "Numérateur : commutations ; dénominateur : transitions évaluables. IC bootstrap familial descriptif à 95 %, 24 unités ; profils appariés. Source : métriques closes " + stage + "."),
                 (rate_intervals(rows, "rare_variant_phase_error_rate", "Erreurs de phase aux variants rares (MAF ≤5 %)"), "Ancrage hors cible, fréquence calculée sur les fondateurs vrais. IC familial descriptif, pas un seuil de qualification. Source : métriques closes " + stage + ".")]
        if stage == "12":
            plots.append((rate_intervals(rows, "target_mutation_orientation_accuracy", "Orientation de la cible avec ancrage indépendant"), "Numérateur : orientations correctes ; dénominateur : hétérozygotes cibles ancrables. Les homozygotes et égalités ne sont pas des succès. Source : oracle 12."))
        else:
            runs = available[0]["error_run_lengths"]
            plots.append((count_bars(["1 hétérozygote", "2–5 hétérozygotes", "6 et plus"], list(runs.values()), "Longueur des séquences de phase discordante", "Nombre de séquences — cas de base"), "Séquences après alignement global ; description de motifs, pas identification causale d’une erreur ponctuelle ou de bloc. Source : oracle 12A."))
    return ''.join(figure(index, svg, caption) for index, (svg, caption) in enumerate(plots, 1))


def report_sections(manifest: dict[str, Any], document: dict[str, Any], source_sha: str) -> list[str]:
    stage, rows = document["stage_id"], document["evaluations"]
    baseline_summary = ""
    if stage in {"12", "12A"}:
        baseline = rows[0].get("metrics", {})
        for name, label in (("switch_error_rate", "Commutations"), ("rare_variant_phase_error_rate", "Erreurs de phase rare")):
            measure = baseline.get(name, {}).get("ALL")
            if measure is not None and measure["estimate"] is not None:
                baseline_summary += f'<p>Profil de base — {label} : {100*measure["estimate"]:.2f} % ({measure["numerator"]}/{measure["denominator"]}).</p>'
        baseline_summary += '<p>Un succès technique n’est pas une preuve de phase fiable. La densité très réduite de ce banc limite son extrapolation à l’IBD et à la puce entière ; une évaluation à plus forte densité reste nécessaire, sans modifier les seuils à partir de ces résultats.</p>'
    table_rows = []
    for row in rows:
        if stage in {"10A", "11"}:
            table_rows.append([row["profile"], row["status"], str(row.get("expected_markers", "Non mesuré")), str(row.get("observed_markers", "Non mesuré"))])
        else:
            measure = row.get("metrics", {}).get("switch_error_rate", {}).get("ALL")
            table_rows.append([row["profile"], row["status"], "Non mesuré" if measure is None else f'{measure["numerator"]}/{measure["denominator"]}',
                               "Non évaluable" if measure is None or measure["estimate"] is None else f'{100*measure["estimate"]:.3f} %'])
    table = _table(("Profil", "Conformité", "Attendu" if stage in {"10A", "11"} else "Switches/transitions", "Observé" if stage in {"10A", "11"} else "SER"), table_rows,
                   "Trois profils appariés, non trois réplications indépendantes. Source : évaluation close " + stage + ".")
    if stage == "12":
        table += '<p>Qualité signalée par la production (porteurs, homozygotes compris) : ' + escape('; '.join(
            f'{row["profile"]} : {row.get("reliable_carrier_count", "non mesuré")} fiables, {row.get("unreliable_carrier_count", "non mesuré")} non fiables ; revue manuelle={row.get("published_manual_review_required", "non mesuré")}' for row in rows)) + '.</p><p>La précision d’orientation est mesurée sur tous les hétérozygotes ancrables, sans retirer a posteriori les appels signalés non fiables. Elle ne vaut pas approbation de leur usage scientifique.</p>'
    details = ""
    if stage in {"12", "12A"}:
        detail_rows = []
        for row in rows:
            for name, title in (("switch_error_rate", "SER"), ("rare_variant_phase_error_rate", "Phase rare")):
                for label in ("MAF_LE_005", "MAF_GT_005", "GAP_LE_100KB", "GAP_GT_100KB"):
                    if name == "rare_variant_phase_error_rate" and label == "MAF_GT_005":
                        continue
                    measure = row.get("metrics", {}).get(name, {}).get(label)
                    if measure is not None:
                        detail_rows.append([row["profile"], title, label, f'{measure["numerator"]}/{measure["denominator"]}', "NE" if measure["estimate"] is None else f'{100*measure["estimate"]:.3f} %'])
        details = _table(("Profil", "Mesure", "Strate", "Erreurs/n", "Taux"), detail_rows,
                         "Strates marginales : fréquence fondatrice et écart au marqueur précédent. Tous les IC et raisons de non-évaluation figurent dans report_data.json.").replace("Tableau 1", "Tableau 2")
    question = {"10A": "Le filtrage complet/polymorphe conserve-t-il exactement les sujets, marqueurs et dosages attendus sur les 22 autosomes ?",
                "11": "La fenêtre cible conserve-t-elle les bons variants et génotypes, avec une carte bp/cM correctement interpolée et la cible unique ?",
                "12": "Quelle phase régionale est reconstruite, et peut-on orienter la cible sans utiliser la cible elle-même pour aligner les copies ?",
                "12A": "Quelle phase est reconstruite sur chaque autosome, avec le même ensemble d’individus, sans donneur de simulation dans la référence ?"}[stage]
    criterion_rows = [[row["profile"], criterion["criterion_id"], "Non mesuré" if criterion["value"] is None else str(criterion["value"]), criterion["status"]]
                      for row in rows for criterion in row["criteria"]]
    criteria_table = _table(("Profil", "Critère D04", "Valeur", "Décision"), criterion_rows,
                            "Raccord explicite aux critères approuvés avant pilote, version 0.2.0 ; seuils de performance non scellés.").replace("Tableau 1", "Tableau 3" if stage in {"12", "12A"} else "Tableau 2")
    execution_scope = ""
    if stage == "12":
        execution_scope = '<p>Parcours réel 12 : phase_common puis phase_rare pour la cible absente de référence, avec pedigree synthétique transmis à l’adaptateur. La fenêtre conserve les marqueurs régionaux, y compris ceux qui sont monomorphes dans la seule étude.</p>'
    elif stage == "12A":
        execution_scope = '<p>Parcours réel 12A : phase_common séparément sur 22 autosomes, après filtre complet/polymorphe et appariement de référence. Aucun argument de pedigree n’est passé à SHAPEIT5 par ce parcours ; la présence du registre en entrée ne signifie pas que les trios y contraignent la phase. La cible ajoutée est traitée séparément par 12.</p>'
    return [
        f'<p>CAL-REP-{stage} — {escape(LABELS[stage])}. Effet_fondateur2, version 1.0.0, {escape(manifest["created_at"][:10])}. Rapport indépendant de développement D23.</p>',
        f'<p>Conformité du banc : {document["status"]}. Trois observations appariées des mêmes 40 sujets, 22 autosomes disponibles, 24 unités de simulation. Qualification de performance : {document["performance_qualification"]}.</p>' + baseline_summary + '<aside class="callout">Ce contrôle ne qualifie ni le pipeline entier, ni les seuils scientifiques, ni un effet fondateur. calibrated=false ; founder_effect_proven=false.</aside>',
        '<p>' + question + '</p>',
        '<p>Autorité : D23_PHASING_PRESPECIFICATION.md, critères D04 approuvés avant pilote, D23a v2 et contrats des composants. Sources publiques épinglées : 1000G 20220422 GRCh38, puce CytoScan 750K Accel NA36.r1, cartes HapMap/SHAPEIT D09. Aucun échantillon d’étude, aucun génotype du panel synthétique de validation.</p>',
        '<p>Phasage : SHAPEIT5 5.1.1, graine 2303, un thread ; masquage 2304 ; bootstrap 2305, 2 000 tirages. Profils BASELINE, un marqueur sur deux DENSITY_50 et masquage Bernoulli diploïde 1 % MISSING_1PCT, cible protégée. Les seuils de production 10A/12A ne sont pas abaissés. Paramètres et empreintes sont conservés dans report_data.json ; commandes et logs restent privés dans les reçus d’intégration.</p>',
        '<p>32 fondateurs issus de 64 donneurs de calibration D10 sans groupe documenté partagé ; 8 enfants issus de méioses Haldane, 8 trios et 16 sujets isolés. Vérité phasée conservée à part. La cible technique APC chr5 est ajoutée explicitement et transmise, sans statut clinique. Sa fréquence est artificielle ; ce n’est pas une calibration d’une mutation biologique rare DOCK6.</p>',
        '<p>Les outils ne consomment que GT non phasés, pedigree synthétique déclaré, génotypes moléculaires explicites, cartes et références publiques D10 sélectionnées. Conversion VCF–PLINK réelle ; QC/cohorte amont portent D23_UPSTREAM_CONTRACT_FIXTURE. 03–09 ne sont pas prétendues exécutées ou qualifiées. La vérité est ouverte après clôture par un évaluateur distinct.</p>',
        '<p>Oracles séparés des résumés PASS : décodeur BED indépendant, comptages de complétude/polymorphisme, fenêtre et interpolation NumPy ; comparaison des BCF aux copies vraies par sujet et locus, recodage REF/ALT explicite. Les mutations d’axes et dosages sont des défauts, pas des erreurs de phase.</p>' + execution_scope + equations(stage),
        '<p>La figure 1 répond à la question principale ; les figures suivantes précisent couverture ou rareté. Les tableaux gardent les numérateurs et dénominateurs exacts.</p>' + graphics(stage, rows) + table + details,
        '<p>Critères fixes : aucune discordance d’univers/dosage ou de carte, cible retenue pour 11/12, référence D10 exacte. Les taux mesurés ne deviennent pas PASS faute de seuil : switch_error_rate, rare_variant_phase_error_rate et target_mutation_orientation_accuracy sont pilot_then_seal, donc qualification REVIEW_REQUIRED.</p>' + criteria_table,
        '<p>Premier essai 11 bloqué : l’adaptateur du banc omettait sample_set_id/variant_set_id des artefacts malgré leur présence dans le descripteur. Premier essai 12 bloqué : colonnes de l’audit cible hors ordre canonique. Corrections limitées aux enveloppes et à l’ordre des colonnes ; valeurs génétiques inchangées. Nouveau reçu régional 11/12, anciens essais conservés. Aucun filtre ou méthode scientifique de production changé.</p>',
        '<p>Sous-échantillon de puce, 40 sujets et seulement 24 unités : limites de puissance et d’extrapolation. Le bootstrap est descriptif, les marqueurs sont corrélés. Les strates sont marginales, pas un plan factoriel indépendant. Le premier locus de chaque chromosome n’a pas d’écart précédent et reste hors strates de densité. La cellule « rare / MAF &gt;5 % » est hors domaine par définition, pas un prérequis manquant. La phase d’un singleton isolé ou d’un ancrage à égalité reste non évaluable ; zéro erreur sur un petit dénominateur n’est pas une garantie.</p>'
        '<p>12A exclut la cible ajoutée absente de sa référence et la traite via 12 ; cela ne constitue pas une perte inattendue. Les deux parcours diffèrent par leurs univers, densités et usage du pedigree : la comparaison de leurs SER n’isole pas un effet causal de la seule densité. Les motifs de points et blocs ne permettent pas d’identifier leur cause. Une phase exacte ne démontre ni IBD, ni datation, ni effet fondateur. Le FAIL scientifique 07 de D21 demeure ouvert.</p>',
        f'<p>{document["status"]} de conformité sur ce banc seulement. Les métriques et intervalles n’autorisent pas de changement automatique des seuils ni de calibrated. Poursuite : D24 (haplotype fondateur, segment récupéré contre segment vrai), après lecture et validation humaine de ces rapports. Le pilote D33 et la validation indépendante restent requis.</p>',
        '<p>Évaluation close : <span class="hash">' + source_sha + '</span>. D23a : <span class="hash">' + manifest["preparation_manifest_sha256"] + '</span>.</p><p>report_data.json contient uniquement agrégats, configuration et empreintes. Les graphiques reprennent les mesures closes sans nouveau test, qualification ou estimation. Les sources et sorties individuelles restent dans le workspace ignoré.</p>',
    ]


def build_reports(*, evaluation_root: Path, template_path: Path, output_root: Path) -> list[Path]:
    """Publie quatre HTML autonomes avec le gabarit report-formater fourni."""
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    manifest = checked_closure(evaluation_root)
    validate_evaluation(manifest, manifest=True)
    source_sha = sha256_file(evaluation_root / "run_index.json")
    head = re.split(r"<body\b", template_path.read_text(), maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d23_report_template_invalid")
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}main section{break-inside:avoid}main #s8,main #s9{break-inside:auto}main #s9{break-before:page}.scientific-figure,.equation{break-inside:avoid}.scientific-figure{margin:1em 0}table{font-size:.7em}th,td{overflow-wrap:anywhere}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.equation p{font-size:.85em}.math{font:1em Georgia,serif;line-height:1.8}main p{orphans:3;widows:3}</style></head>')
    destination.mkdir(mode=0o700)
    outputs = []
    for stage in STAGES:
        document = json.loads((evaluation_root / (stage + ".evaluation.json")).read_text())
        validate_evaluation(document)
        title = "CAL-REP-" + stage + " — " + LABELS[stage]
        toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
        body = f'<body lang="fr"><section class="cover"><h1>{title}</h1><p>Effet_fondateur2 · {manifest["created_at"][:10]} · version 1.0.0</p><p>Évaluation D23 : résultats, graphiques et méthodes mathématiques</p></section>' + toc
        body += f'<footer class="page-footer"><span>CAL-REP-{stage} · Effet_fondateur2</span></footer><main>'
        body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, report_sections(manifest, document, source_sha), strict=True), 1))
        directory = destination / ("CAL-REP-" + stage)
        directory.mkdir(mode=0o700)
        path = directory / "index.html"
        path.write_text(re.sub(r"<title>.*?</title>", f"<title>{escape(title)}</title>", head, flags=re.S) + body + '</main></body></html>')
        write_json(directory / "report_data.json", dict(manifest=manifest, results=document))
        (directory / "checksums.sha256").write_text(''.join(f'{sha256_file(item)}  {item.name}\n' for item in (path, directory / "report_data.json")))
        outputs.append(path)
    if sha256_file(evaluation_root / "run_index.json") != source_sha:
        raise ValueError("d23_report_source_changed")
    return outputs


def main() -> int:
    """CLI du rendu seul ; aucune vérité ni sortie individuelle en argument."""
    parser = argparse.ArgumentParser(description="Quatre rapports D23")
    for name in ("evaluation-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    print("Rapports créés :", len(build_reports(**vars(parser.parse_args()))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
