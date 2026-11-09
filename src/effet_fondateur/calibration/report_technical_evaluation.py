"""Rapports D19 indépendants : CAL-DEV-D19 et un CAL-REP par étape 00–04.

Le rendu ne recalcule aucune métrique. Il publie seulement les agrégats
contrôlés, sans clés individuelles, génotypes ou correspondance de profils.
"""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.calibration.evaluate_technical import contained_file, load_checked_evaluation
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.technical_evaluation import STAGES, TechnicalEvaluationError


HEADINGS = (
    "Identification", "Résumé exécutif", "Question évaluée", "Sources et autorité",
    "Configuration, graines et versions", "Vérité de référence", "Entrées aveugles",
    "Méthode d’évaluation", "Résultats et incertitude", "Critères préspecifiés",
    "Anomalies et écarts", "Limites", "Décision et suite", "Provenance et empreintes",
)
LABELS = dict(zip(STAGES, ("Initialisation", "Sources ACPA", "Registre des échantillons", "Conversion PLINK", "Préparation de la cible"), strict=True))
METHODS = {
    "00": "Comparaison de la configuration résolue à l’entrée validée ; contrôle des empreintes de configuration et d’environnement, de Python 3.12 et du résumé d’initialisation reconstruit indépendamment.",
    "01": "Comparaison de l’inventaire publié à l’inventaire D18 scellé : fichiers, effectifs, ordre, tailles, empreintes, nombre de marqueurs et statut. Une substitution reste détectable même si les empreintes internes du run sont recalculées.",
    "02": "Comparaison de chaque colonne du registre publié au registre neutre D18, avec contrôle des pertes, ajouts, doublons, changements de valeur et d’ordre.",
    "03": "Décodage indépendant des deux triplets BED/BIM/FAM ; contrôle des descripteurs, des axes, de l’alignement et de tous les génotypes observés, manques inclus. Les bases nucléotidiques sont comparées, pas un dosage supposant A1=REF. Un réalignement permet de compter les discordances mais ne masque pas un changement d’ordre.",
    "04": "Comparaison du jeu enrichi à la base observée et à une attente cible synthétique explicite, séparée et empreintée. La cible doit être absente de la base, et les appels de fond doivent être préservés. Aucune inférence depuis le statut clinique ou le groupe n’est autorisée.",
}


def _table(headers: tuple[str, ...], rows: list[list[str]], caption: str) -> str:
    """Table agrégée avec échappement systématique du contenu."""
    return '<table><caption>Tableau 1 — ' + escape(caption) + '</caption><thead><tr>' + ''.join(
        '<th>' + escape(value) + '</th>' for value in headers
    ) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join('<td>' + escape(value) + '</td>' for value in row) + '</tr>' for row in rows) + '</tbody></table>'


def _sections(manifest: dict[str, Any], documents: dict[str, dict[str, Any]], stage: str | None, manifest_sha: str) -> list[str]:
    """Texte scientifique propre au stade ; absence de mesure jamais remplacée par zéro."""
    selected = manifest["stages"] if stage is None else [item for item in manifest["stages"] if item["stage_id"] == stage]
    extended = "scenario" in manifest
    rejection = manifest.get("scenario") == "TARGET_PRESENT_REJECTION"
    report_id = "CAL-DEV-D19" if stage is None else f"CAL-REP-{stage}"
    counts = [sum(item[key] for item in selected) for key in ("pass_count", "fail_count", "not_evaluable_count")]
    cards = '<figure><div class="stat-grid">' + ''.join(
        f'<div class="stat-card"><div class="num">{count}</div><div class="label">{label}</div></div>'
        for count, label in zip(counts, ("PASS techniques", "FAIL techniques", "Non évaluables"), strict=True)
    ) + '</div><figcaption>Figure 1 — Verdicts lot × étape. Source : résultats D19 scellés ; les profils appariés ne sont pas des réplications indépendantes.</figcaption></figure>'
    flow = '<figure><svg viewBox="0 0 680 115" role="img" aria-label="Exécution aveugle puis évaluation séparée">'
    for index, (first, second) in enumerate((("Entrées techniques" if extended else "Entrées D18", "Observations et cible explicite" if extended else "Observations seules"), ("Pipeline terminé", "Reçu CLOSED empreinté"), ("Évaluation séparée", "Référence et sorties"))):
        horizontal = 8 + index * 224
        flow += f'<rect x="{horizontal}" y="30" width="215" height="65" rx="5" fill="#eaf1f5" stroke="#236c75"/><text x="{horizontal+107}" y="54" text-anchor="middle" font-size="13">{first}</text><text x="{horizontal+107}" y="75" text-anchor="middle" font-size="11">{second}</text>'
    flow += '</svg><figcaption>Figure 2 — Séparation procédurale. Source : contrats D18/D19 ; la référence n’est pas fournie au processus pipeline.</figcaption></figure>'
    if stage is None:
        rows = [[item["stage_id"], LABELS[item["stage_id"]], item["status"], str(item["pass_count"]), str(item["fail_count"]), str(item["not_evaluable_count"])] for item in selected]
        table = _table(("Étape", "Contrôle", "Verdict", "PASS", "FAIL", "NE"), rows, "Synthèse des cinq évaluations. Source : manifest D19.")
        method = "Cinq évaluateurs distincts vérifient les fichiers, effectifs, ordres, appels et empreintes. Les tests injectent pertes, doublons, substitutions et permutations, y compris après recalcul des empreintes internes. Le décodeur BED dispose d’un oracle binaire issu de la documentation primaire PLINK, indépendant de l’encodeur de fixture."
    else:
        rows = []
        for evaluation in documents[stage]["evaluations"]:
            metrics = evaluation["metrics"]
            compared = metrics.get("compared_count", metrics.get("regional_base_genotype_compared_count", metrics.get("genotype_compared_count")))
            rows.append([evaluation["bundle_id"], evaluation["status"], "Non mesuré" if compared is None else str(compared), ", ".join(evaluation["reasons"]) or "Aucun défaut détecté"])
        table = _table(("Lot opaque", "Verdict", "Lignes/cellules comparées", "Motif"), rows, "Évaluation par lot. Source : résultat indépendant D19 ; motifs codés, sans information individuelle.")
        method = METHODS[stage]
    decision = "FAIL : défaut technique à corriger avant qualification." if counts[1] else (
        "NOT_EVALUABLE : preuves manquantes ; aucune validation de cette étape n’est revendiquée." if counts[2] else
        "PASS technique pour les contrôles définis et les lots exécutés ; ce verdict ne valide pas l’interprétation scientifique."
    )
    if rejection:
        decision += " Scénario négatif : un PASS en 04 valide le refus de doublon, pas une insertion réussie ; le pipeline reste FAILED/BLOCKED."
        if stage == "04":
            method = "Contrôle négatif préspecifié : la cible est déjà présente dans une copie des observations. Vérification du code V2=4, de l’état FAILED, du motif exact target_variant_already_present empreinté, de la tentative liée à la signature, de l’absence de publication 04 et de la base PLINK conservée. Un échec pour une autre raison ne vaut jamais PASS. Le taux d’erreur d’insertion reste null, car aucune insertion n’est réalisée."
    hashes = '<ul>' + ''.join(f'<li>{escape(label)} : <span class="hash">{escape(value)}</span></li>' for label, value in (
        ("Manifest D19", manifest_sha), ("Manifest D18 source", manifest["d18_manifest_sha256"]),
        ("Index des runs clos", manifest["run_index_sha256"]),
    )) + '</ul>'
    if extended:
        hashes += f'<p>Manifest des entrées techniques : <span class="hash">{manifest["technical_panel_manifest_sha256"]}</span>.</p>'
    if stage is not None:
        artifact = next(item for item in manifest["artifacts"] if item["stage_id"] == stage)
        hashes += f'<p>Résultat {escape(report_id)} : <span class="hash">{artifact["sha256"]}</span>.</p>'
    sections = [
        f'<p>{report_id} · version 1.0 · {escape(manifest["created_at"][:10])}. Responsable : contrôle de développement Effet_fondateur2. Périmètre : technique régional.</p>',
        f'<p>{decision} Les compteurs de la figure 1 portent sur des évaluations, pas sur des individus indépendants.</p>' + cards,
        '<p>Les étapes conservent-elles les entrées attendues sans perte, duplication, substitution ou changement d’ordre non signalé ?</p>',
        '<p>AGENTS.md ; conventions de calibration ; plan D19 ; manifest et inventaire D18 ; schémas du pipeline et contrats D19. '
        'Le format binaire est vérifié contre la <a href="https://www.cog-genomics.org/plink/1.9/formats#bed">spécification primaire BED de PLINK 1.9</a>.</p>',
        f'<p>Méthode {escape(manifest["method_id"])}. Évaluateurs déterministes, aucune nouvelle graine et aucun réglage de seuil. Les versions d’exécution sont scellées dans environment.json et liées au reçu de clôture. Python 3.12 est contrôlé en 00. Le rendu n’appelle aucun outil génétique.</p>',
        '<p>Pour 00–03, référence = configuration, inventaire, registre neutre et VCF observé D18. Les erreurs déjà injectées dans l’observation doivent être conservées, pas corrigées à partir de la vérité biologique. Pour 04, une attente cible indépendante et explicite est obligatoire ; elle n’est pas disponible dans le présent panel.</p>',
        '<p>Le processus pipeline reçoit uniquement les entrées D18. L’évaluation commence après clôture des processus et scellement du reçu ; les correspondances et profils restent privés. Voir figure 2.</p>' + flow,
        f'<p>{method}</p><p>Manifest → audit → sorties → artefacts sont reliés par empreintes et identités. L’évaluateur ne lance aucune étape et vérifie l’immuabilité avant publication atomique.</p>',
        table + '<p>Les métriques sans comparaison sont null en JSON. Aucun intervalle de performance scientifique n’est estimé : il s’agit de concordance technique exacte sur un panel régional de développement, non d’une validation indépendante biologique.</p>',
        '<p>Seuil technique : zéro perte, ajout, doublon, changement de valeur, désordre et discordance alignable ; empreintes et contrats cohérents. FAIL prime sur les autres résultats. Une étape non exécutée est NOT_EVALUABLE ; une sortie annoncée réussie mais absente ou corrompue est FAIL. Absence de vérité cible : NOT_EVALUABLE, jamais PASS.</p>',
        '<p>Le tableau 1 expose les verdicts observés. Les défauts injectés dans les tests restent confinés à des dossiers temporaires. Aucun résultat de production n’est réparé par l’évaluateur et aucune approbation du registre n’est fabriquée.</p>',
        '<p>Le rejeu D19 est limité à 00–02 : les configurations D18 n’autorisent pas 03–04. Le registre attend encore une approbation humaine. La cible est déjà présente parmi les marqueurs de base, alors que le contrat 04 attend une insertion nouvelle : une politique explicite sera nécessaire avant ce test intégré. Les évaluateurs 03–04 sont testés sur de petits jeux contrôlés, sans validation intégrée de ces étapes. La région chr5 ne qualifie pas un pipeline génome entier. Les profils appariés ne sont pas indépendants.</p>'
        '<p>Les verdicts portent sur la concordance D19, pas sur l’ensemble des critères transversaux D04/D08. Le journal système des accès n’est pas instrumenté ici ; aucun compteur d’accès réel égal à zéro n’est revendiqué. Le scellement global D06/D07 et les portes de passage restent requis avant qualification.</p>',
        f'<p>{decision} calibrated=false ; founder_effect_proven=false. Prochaine étape de développement : D20, évaluateurs QC 05 et 10. Une qualification complète exigera d’autres panels, les prérequis de chaque étape et les portes de passage prévues au plan.</p>',
        hashes,
    ]
    if extended:
        target = manifest["technical_target"]
        annotation = f'{escape(target["gene"])} · GRCh38 chr{target["chromosome"]}:{target["position_bp"]} · {escape(target["ref"])}&gt;{escape(target["alt"])} · {escape(target["transcript"])}:{escape(target["hgvs_c"])}'
        sections[3] += f'<p>Annotation de la cible technique : {annotation}. Source primaire : <a href="{escape(target["source_url"], quote=True)}">{escape(target["source_id"])}</a>. Seules référence et annotation sont utilisées, sans classification clinique ni sélection de patients.</p>'
        sections[4] += f'<p>PLINK est exécuté réellement en 03–04, avec version consignée dans les audits. L’approbation du registre neutre est fournie par la décision {escape(target["decision_id"])} ; le runner ne la crée pas.</p>'
        sections[5] = '<p>Les observations et le registre source D18 sont conservés dans le scénario principal. Une cible technique séparée, absente de ces sondes, reçoit une attente synthétique scellée avant exécution. Les appels RR, RA et AA suivent le rang canonique du registre modulo 3. Ils ne viennent ni des patients, ni de leur clinique, ni des donneurs 1000G ; cette distribution de stress n’est pas une fréquence populationnelle. La confirmation du contrat concerne l’annotation publique, pas une mesure expérimentale sur les sujets.</p>'
        if rejection:
            sections[5] += '<p>Le négatif reprend un lot représentatif et ajoute exactement cette cible à une copie du VCF et des ACPA. Toutes les observations préexistantes restent identiques ; aucun fichier D18 n’est modifié.</p>'
        sections[6] = '<p>Le pipeline reçoit uniquement les copies aveugles, l’annotation publique et les génotypes cibles explicitement simulés nécessaires à 04. Aucun chemin de vérité privée ou profil D17 n’entre dans sa configuration. Le reçu CLOSED précède l’évaluation séparée. Voir figure 2.</p>' + flow
        sections[10] = '<p>Les verdicts observés sont au tableau 1. Les anciens runs 00–02 restent immuables. Les nouvelles configurations autorisent explicitement 00–04 ; les sources D18 sont vérifiées avant/après. Le scénario de doublon est évalué séparément et son échec attendu n’est jamais confondu avec une réussite d’insertion.</p>'
        sections[11] = '<p>Ce diagnostic qualifie seulement la concordance technique régionale sur les lots exécutés. La cible technique est distincte de la cible biologique D14–D18 et peut être hors de leur fenêtre : aucun haplotype, âge, LD, ROH ou résultat IBD n’est évalué ici. Le registre conserve sexe/clinique inconnus et filiation non renseignée ; le contrôle Mendel sans trio complet reste NOT_APPLICABLE, pas un test réussi. Les profils sont appariés. Aucun patient positif ni validation biologique externe n’est ajouté.</p><p>Les critères transversaux D04/D08 et le journal système d’accès non instrumenté ici restent hors de ces PASS. Le scellement global D06/D07 et les portes du pilote restent requis ; calibrated=false et founder_effect_proven=false.</p>'
    return sections


def build_reports(*, run_root: Path, template_path: Path, output_root: Path) -> list[Path]:
    """Publie six HTML autonomes, données et empreintes dans un dossier neuf."""
    parts = tuple(part.casefold() for part in output_root.absolute().parts)
    if output_root.absolute().parent.name != "calibration_workspace" or output_root.exists() or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(path.is_symlink() for path in (output_root, *output_root.parents)):
        raise TechnicalEvaluationError("d19_report_output_invalid")
    manifest = load_checked_evaluation(run_root)
    source_sha = _sha(run_root / "private_truth/d19/technical_evaluation_manifest.json")
    documents = {item["stage_id"]: _json(contained_file(run_root, item["path"])) for item in manifest["artifacts"]}
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    if "</head>" not in head:
        raise TechnicalEvaluationError("d19_report_template_invalid")
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}table{font-size:0.78em}th,td{overflow-wrap:anywhere}main p{orphans:3;widows:3}main section{break-inside:avoid}</style></head>')
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=output_root.parent))
    paths = []
    try:
        for stage in (None, *STAGES):
            report_id = "CAL-DEV-D19" if stage is None else f"CAL-REP-{stage}"
            title = "D19 — Évaluateurs techniques" if stage is None else f"{report_id} — {LABELS[stage]}"
            sections = _sections(manifest, documents, stage, source_sha)
            toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(
                f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
                for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
            body = f'<body lang="fr"><section class="cover"><h1>{escape(title)}</h1><p>Effet_fondateur2 · {report_id} · version 1.0 · {escape(manifest["created_at"][:10])}</p></section>' + toc
            body += f'<footer class="page-footer"><span>{report_id} · Effet_fondateur2</span></footer><main>'
            body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, sections, strict=True), 1)) + '</main></body></html>'
            directory = temporary / report_id
            directory.mkdir(mode=0o700)
            output = directory / "index.html"
            output.write_text(re.sub(r"<title>.*?</title>", f"<title>{escape(title)}</title>", head, flags=re.S) + body, encoding="utf-8")
            (directory / "report_data.json").write_text(json.dumps(manifest if stage is None else documents[stage], indent=2) + "\n", encoding="utf-8")
            (directory / "checksums.sha256").write_text(''.join(f"{_sha(path)}  {path.name}\n" for path in (output, directory / "report_data.json")), encoding="utf-8")
            paths.append(output_root / report_id / "index.html")
        if source_sha != _sha(run_root / "private_truth/d19/technical_evaluation_manifest.json"):
            raise TechnicalEvaluationError("d19_report_source_changed")
        load_checked_evaluation(run_root)
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(output_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return paths


def main() -> int:
    """CLI des six rapports indépendants ; PDF après contrôle de pagination."""
    parser = argparse.ArgumentParser(description="Rapports indépendants D19")
    for name in ("run-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    build_reports(**vars(parser.parse_args()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
