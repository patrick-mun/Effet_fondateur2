"""D22 : trois rapports autonomes depuis les agrégats clos, sans recalcul PCA."""

from __future__ import annotations

import argparse
from html import escape
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.report_technical_evaluation import HEADINGS, _table
from effet_fondateur.calibration.structure_bench import COHORTS, checked_document
from effet_fondateur.calibration.structure_report_evidence import collect_visual_evidence
from effet_fondateur.calibration.structure_report_graphics import cohort_counts, control_exclusions, detection_intervals, distance_histogram, figure, pca_centroids, scree
from effet_fondateur.calibration.structure_report_math import cohort_equations, pca_equations, wilson_equation
from effet_fondateur.calibration.technical_scenarios import new_workspace_run


LABELS = dict(zip(COHORTS, ("Témoins indépendants", "Tous porteurs", "Porteurs indépendants", "Non-porteurs familiaux", "Atteints non-porteurs", "Chromosome cible", "Référence de structure"), strict=True))


def report_sections(manifest: dict[str, Any], documents: dict[str, Any], stage: str | None, source_sha: str, visuals: dict[str, Any]) -> list[str]:
    """Préserve les verdicts, limites de dépendance et portée des fixtures."""
    selected = [row for row in manifest["stages"] if stage is None or row["stage_id"] == stage]
    status = "FAIL" if any(row["status"] == "FAIL" for row in selected) else "PASS"
    population, cohorts = documents["08"]["evaluations"], documents["09"]["evaluations"]
    cards = '<aside class="callout">Quatre cas appariés · 108 références et 20 projetés par cas. Les cas ne sont pas des réplications indépendantes.</aside>'
    base = visuals["cases"][0]
    if stage == "08":
        plots = [(pca_centroids(base), "Où se placent groupes et projetés ? Barycentres et ± un écart-type, pas IC ; positions agrégées, jamais sujets individuels. Source : population_scores.tsv clos, visual_data.json, cas de base."),
                 (scree(base), "Quel poids ont les axes ? Seuls les trois axes demandés sont représentés, pas toute la variance. Source : population_eigenvalues.tsv clos, cas de base."),
                 (distance_histogram(base), "Pourquoi des sujets sont-ils signalés ? Effectifs par classe de distance, seuil α=0,01 et D²≈9,21 ; voir équation 2. Source : population_outliers.tsv clos, cas de base."),
                 (detection_intervals(population[0]), "Quelle incertitude reste-t-il malgré 4/4 ? Points = taux observés, segments = IC Wilson 95 % descriptifs. Source : CAL-REP-08.json, références du cas de base ; copies projetées exclues.")]
    elif stage == "09":
        plots = [(cohort_counts(base), "Qui est retenu ? Les sept cohortes se chevauchent : ne pas additionner les barres. Non inclus = complément parmi les 128 sujets, pas seulement parmi les candidats. Source : cohorts.frozen.tsv clos, cas de base."),
                 (control_exclusions(base), "Pourquoi ne pas être témoin ? Motifs prioritaires mutuellement exclusifs ; ce n’est pas le nombre de tous les défauts possibles. Source : cohorts.frozen.tsv clos, controls_unrelated, cas de base.")]
    else:
        plots = [(pca_centroids(base), "Structure visible sans publier de sujet : moyennes ± dispersion, pas intervalles de confiance. Source : population_scores.tsv clos, cas de base."),
                 (detection_intervals(population[0]), "Un taux parfait observé n’est pas une certitude. Source : CAL-REP-08.json, Wilson 95 % descriptif sur références indépendantes."),
                 (cohort_counts(base), "Les cohortes se chevauchent et ne sont pas additionnables. Source : cohorts.frozen.tsv clos, cas de base ; compléments parmi 128 sujets."),
                 (control_exclusions(base), "Chaque non-inclusion témoin a un motif prioritaire. Source : cohorts.frozen.tsv clos, cas de base.")]
    graphics = ''.join(figure(number, svg, caption) for number, (svg, caption) in enumerate(plots, 1))
    if stage is None:
        table = _table(("Étape", "Cas mesurés", "FAIL", "Verdict"), [[row["stage_id"], str(row["case_count"]), str(row["fail_count"]), row["status"]] for row in selected], "Synthèse des composants réellement évalués. Source : manifest D22.")
    elif stage == "08":
        table = _table(("Cas", "Loci informatifs", "Séparation¹", "Mélanges²", "VP/FN/FP/VN³", "Verdict"), [
            [row["case"], str(row["metrics"]["informative_variant_count"]), f'{row["metrics"]["centroid_separation_over_within_sd"]:.3f}',
             f'{row["metrics"]["mixed_coordinate_mean"]:.6f}', "/".join(str(row["metrics"]["independent_reference_detection"][key]) for key in ("tp", "fn", "fp", "tn")), row["status"]] for row in population],
            "Source : CAL-REP-08.json. ¹Distance des centroïdes / dispersion intra-groupe ; ²coordonnée moyenne entre centroïdes (attente descriptive 0,5), pas proportion d’ascendance estimée ; ³référence indépendante uniquement.")
    else:
        table = _table(("Cohorte", "Base", "Projetés modifiés", "Rôles permutés", "Manque 5 %"), [
            [LABELS[cohort], *[str(row["metrics"]["included_counts"][cohort]) for row in cohorts]] for cohort in COHORTS],
            "Effectifs inclus, non données d’étude. Source : CAL-REP-09.json. Les changements de rôles descriptifs peuvent légitimement modifier les cohortes, pas les axes PCA.")
        table += _table(("Contrôle négatif", "Code réel", "Résultat attendu"), [[row["case"], str(row["returncode"]), "Refus sans cohorte publiée"] for row in cohorts[0]["review_rejections"]],
                        "Gardes de revue réellement exécutés dans le cas de base. Source : reçu clos et logs empreintés.").replace("Tableau 1", "Tableau 2")
    intervals = []
    for row in population:
        detection = row["metrics"]["independent_reference_detection"]
        interval = detection["sensitivity_wilson"]
        false_interval = detection["non_injected_rate_wilson"]
        intervals.append(f'{row["case"]} : sensibilité {detection["tp"]}/{detection["tp"]+detection["fn"]}, IC Wilson 95 % {100*interval[0]:.2f}–{100*interval[1]:.2f} % ; alertes négatives {detection["fp"]}/{detection["fp"]+detection["tn"]}, IC {100*false_interval[0]:.2f}–{100*false_interval[1]:.2f} %.')
    reasons = sorted({reason for key in ("08", "09") if stage is None or key == stage for row in documents[key]["evaluations"] for reason in row["reasons"]})
    config = manifest["configuration"]
    sections = [
        f'<p>{"CAL-DEV-D22" if stage is None else "CAL-REP-" + stage} · présentation 2.0 · analyse close du {escape(manifest["created_at"][:10])}. Effet_fondateur2 : évaluation technique de structure et gel. Seul le rendu change, jamais les mesures ou verdicts.</p>',
        f'<p>Verdict : {status} sur le périmètre mesuré. Les calculs et décisions sont confrontés à des oracles distincts. Ce résultat ne qualifie ni l’indépendance de 07, ni une ascendance humaine, ni le pipeline entier.</p>' + cards,
        '<p>08 calcule-t-elle et projette-t-elle une PCA sans influence des sujets projetés ou des étiquettes descriptives ? 09 produit-elle les sept cohortes attendues, reproductibles, et refuse-t-elle les revues invalides ?</p>',
        '<p>Autorité : contrats 08/09, AGENTS.md, conventions et préspecification D22 ; CAL-DEC-006 autorise uniquement ce nouveau banc. '
        'Référence primaire : <a href="https://www.cog-genomics.org/plink/1.9/strat">PLINK 1.9, PCA et projection sur sous-ensemble</a>. La PCA évaluée est celle du module Python, pas la commande --pca de PLINK. L’annotation cible publique épinglée D19 ne constitue pas une mesure sur un individu.</p>',
        f'<p>Graines fond/manque : {config["seed"]}/{config["missingness_seed"]}. 22 autosomes × {config["markers_per_autosome"]} loci IID, trois axes PCA, deux axes d’alerte, alpha {config["population_parameters"]["outlier_alpha"]}. '
        f'Versions : {escape(" ; ".join(manifest["tool_versions"]))}. Paramètres complets, règles de revue et empreintes dans report_data.json.</p>',
        '<p>Référence : 40 sujets A, 40 B, 24 mélanges 50/50 par copie et locus, quatre atypiques rares de contraste orthogonal. Vingt projetés dérivent de copies de référence avant les perturbations d’observation. A/B sont des fréquences artificielles 0,2/0,8, pas des populations géographiques. Les FID sont des groupes techniques, pas des pedigrees biologiques. La cible T/T, T/A ou A/A suit le rang modulo 3, indépendamment de la clinique et du groupe.</p>',
        '<p>Cas base, projetés modifiés, rôles/clinique/lots permutés et manque 5 %. Les vérités de groupe et d’injection restent privées. Les sources 04–07 sont explicitement D22_UPSTREAM_CONTRACT_FIXTURE. Les revues portent le rôle calibration_contract_fixture et sont liées aux sorties après clôture, selon une règle fixée avant le run. Elles ne résolvent aucune décision D21 ou d’étude.</p>',
        '<p>Oracle 08 : lecture indépendante des allèles BED, dosages A1 réels, fréquence et imputation sur la seule référence ; Gram/eigh plutôt que la SVD de production. Les loadings, scores, valeurs propres, distances, comptes et statuts sont comparés. La queue chi-deux à deux degrés de liberté est calculée par exp(-distance²/2). Les invariances comparent les scores publiés après le seul alignement de signe, sans rotation arbitraire.</p>'
        '<p>Oracle 09 : règles de rôle explicites, représentants orientés qualité, motifs, sources, unités, sept keep et décisions vérifiés. Deux gels partagent exactement données et paramètres ; seules les enveloppes de publication et audits horodatés diffèrent.</p>',
        table + ('<p>' + escape(" ".join(intervals)) + '</p><p>Les huit atypiques incluant les copies projetées sont aussi comptés dans le JSON, sans intervalle les assimilant à huit unités indépendantes.</p>' if stage != "09" else '<p>896 lignes de membres et 84 lignes de décisions par gel ; les deux gels de chaque cas sont comparés sur toutes les tables, keep et rapports, pas seulement leurs effectifs. La cohorte de référence reste la sélection d’entrée moins les exclusions techniques acceptées, sans refit PCA après 09.</p>'),
        '<p>Zéro discordance fonctionnelle attendue. CAL-POP-001/002 : sensibilité et alertes des négatifs mesurées, seuils de qualification non scellés. CAL-POP-003 : aucune origine géographique attribuée dans les agrégats publiés ; audit limité à ces textes, pas au système entier. CAL-POP-004 : zéro différence entre les gels scientifiques répétés. Une sortie manquante ou altérée ne devient jamais PASS.</p>',
        '<p>Discordances finales : ' + escape(", ".join(reasons) or "aucune") + '.</p><p>Le premier run a révélé un vrai défaut du schéma population_eigenvalues : REFERENCE_SAMPLE_COUNT refusait les nombres commençant par 1, dont 108. Six cas de régression échouaient avant correction. Le schéma accepte désormais tous les entiers canoniques ≥ 2 ; zéro, un, signes, décimales et zéros initiaux restent refusés. Aucun calcul ni seuil PCA modifié. Le run initial en échec est conservé, les verdicts présents concernent le nouveau run clos.</p>',
        '<p>Banc favorable fortement contrasté, quatre atypiques indépendants et profils appariés : les intervalles restent larges et descriptifs. Aucune qualification populationnelle, calibration des p-values, robustesse sur mélanges humains complexes ou densité réelle de puce n’est démontrée. Pas de données 1000G, LD biologique, phasage, IBD ou datation évalués ici. Les effectifs de familles sont techniques ; aucune vérité de parenté genome-wide n’est qualifiée par ce gel.</p>'
        '<p>La portée est celle de composants 08/09 avec amont simulé, pas de l’intégration réelle 00–09. Aucun journal exhaustif d’accès système n’est instrumenté. Le rôle de revue est une fixture, pas une validation humaine des exclusions. Un PASS de reproductibilité peut reproduire une décision scientifiquement imparfaite.</p>',
        f'<p>{status} de conformité seulement ; performance REVIEW_REQUIRED, calibrated=false et founder_effect_proven=false. Le FAIL 07 de D21 reste ouvert, son étude est reportée après l’ensemble de l’évaluation sur décision de l’utilisateur. Prochaine étape de développement : D23, évaluateurs de phasage 10A/11/12/12A, avec vérification préalable de leurs propres ressources ; aucun lancement automatique.</p>',
        '<p>Méthode : ' + escape(manifest["method_id"]) + '.</p><ul>' + ''.join(f'<li>{label} : <span class="hash">{escape(value)}</span></li>' for label, value in (
            ("Évaluation close", source_sha), ("Entrées", manifest["panel_manifest_sha256"]), ("Reçu des runs", manifest["run_index_sha256"]))) + '</ul>',
    ]
    sections[7] += cohort_equations() if stage == "09" else pca_equations()
    sections[8] = '<p>Les graphiques suivants expliquent les mesures, puis les tableaux conservent les comptes exacts. Les vues détaillées utilisent le cas de base ; les quatre cas restent comparés dans le tableau et le JSON.</p>' + graphics + sections[8]
    if stage != "09":
        sections[9] += wilson_equation()
    sections[13] += '<p>Données graphiques : visual_data.json, avec empreintes des quatre artefacts sources par cas. Aucun score individuel n’y est publié ; les agrégations ne recalculent aucune PCA ni aucun test.</p>'
    return sections


def build_reports(*, run_root: Path, template_path: Path, output_root: Path, panel_root: Path, pipeline_root: Path) -> list[Path]:
    """Publie trois HTML et JSON agrégés, sans écraser un rendu existant."""
    output_root = new_workspace_run(output_root.absolute().parent, output_root.name)
    manifest = checked_document(run_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
    source_sha = _sha(run_root / "evaluation_manifest.json")
    visuals = collect_visual_evidence(panel_root=panel_root, pipeline_root=pipeline_root, evaluation_root=run_root)
    documents = {stage: _json(run_root / ("CAL-REP-" + stage + ".json")) for stage in ("08", "09")}
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d22_report_template_invalid")
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}table{font-size:0.7em}th,td{overflow-wrap:anywhere}main p{orphans:3;widows:3}main section{break-inside:avoid}main #s8,main #s9{break-inside:auto}.scientific-figure,.equation{break-inside:avoid}.scientific-figure{margin:1.2em 0}.equation{padding:0.8em 1em;border-left:3px solid #286d91;background:#f2f6f8;margin:1em 0}.equation p{font-size:0.85em}.equation-title{font-weight:600}.math{font:1.05em Georgia,serif;line-height:1.9}.math>div{margin:0.4em 0}.fraction{display:inline-flex;flex-direction:column;vertical-align:middle;text-align:center;line-height:1.35;margin:0 0.25em}.fraction>span:first-child{border-bottom:1px solid;padding:0 0.2em}.fraction>span:last-child{padding-top:0.1em}</style></head>')
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=output_root.parent))
    paths = []
    try:
        for stage in (None, "08", "09"):
            report_id = "CAL-DEV-D22" if stage is None else "CAL-REP-" + stage
            title = report_id + " — " + ("Structure et gel" if stage is None else "Structure populationnelle" if stage == "08" else "Gel des cohortes")
            toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
            body = f'<body lang="fr"><section class="cover"><h1>{escape(title)}</h1><p>Effet_fondateur2 · présentation 2.0 · analyse du {manifest["created_at"][:10]}</p><p>Graphiques scientifiques et rappels mathématiques</p></section>' + toc
            body += f'<footer class="page-footer"><span>{report_id} · Effet_fondateur2</span></footer><main>'
            body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, report_sections(manifest, documents, stage, source_sha, visuals), strict=True), 1))
            directory = temporary / report_id
            directory.mkdir(mode=0o700)
            path = directory / "index.html"
            path.write_text(re.sub(r"<title>.*?</title>", f"<title>{escape(title)}</title>", head, flags=re.S) + body + '</main></body></html>', encoding="utf-8")
            write_json(directory / "report_data.json", dict(manifest=manifest, results=documents))
            write_json(directory / "visual_data.json", visuals)
            (directory / "checksums.sha256").write_text(''.join(f"{_sha(item)}  {item.name}\n" for item in (path, directory / "report_data.json", directory / "visual_data.json")), encoding="utf-8")
            paths.append(output_root / report_id / "index.html")
        checked_document(run_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
        if _sha(run_root / "evaluation_manifest.json") != source_sha:
            raise ValueError("d22_report_source_changed")
        checked_document(panel_root, "panel_manifest.json", "manifest_state", "FINALIZED")
        checked_document(pipeline_root, "run_index.json", "closure_state", "CLOSED")
        if visuals["source_fingerprints"]["panel_manifest_sha256"] != _sha(panel_root / "panel_manifest.json") or visuals["source_fingerprints"]["run_index_sha256"] != _sha(pipeline_root / "run_index.json"):
            raise ValueError("d22_visual_source_changed")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(output_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return paths


def main() -> int:
    """CLI des rapports français 08/09 et développement D22."""
    parser = argparse.ArgumentParser(description="Rapports indépendants D22")
    for name in ("run-root", "template-path", "output-root", "panel-root", "pipeline-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    for path in build_reports(**vars(parser.parse_args())):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
