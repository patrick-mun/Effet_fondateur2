"""Rapports D21 autonomes : agrégats clos, sans calcul ni clé individuelle."""

from __future__ import annotations

import argparse
from html import escape
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.calibration.kinship_bench import checked_document
from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.report_technical_evaluation import HEADINGS, _table
from effet_fondateur.calibration.simulation.kinship_bench import DEGREES
from effet_fondateur.calibration.technical_scenarios import new_workspace_run


DEGREE_LABELS = ("Clone", "D1", "D2", "D3", "D4", "Non apparenté")


def _percentage(value: float | None) -> str:
    """Formate une proportion mesurée, sans transformer une absence en zéro."""
    return "Non mesuré" if value is None else f"{100 * value:.2f} %"


def _interval(values: list[float] | None) -> str:
    """Formate les bornes déjà calculées par l'évaluateur."""
    return "Non mesuré" if values is None else " – ".join(_percentage(value) for value in values)


def report_sections(manifest: dict[str, Any], documents: dict[str, Any], stage: str | None, manifest_sha: str) -> list[str]:
    """Compose les quatorze sections en distinguant conformité et performance."""
    selected = [row for row in manifest["stages"] if stage is None or row["stage_id"] == stage]
    status = "FAIL" if any(row["status"] == "FAIL" for row in selected) else "PASS"
    config = manifest["configuration"]
    units = len(config["strata"]) * 7 * config["families_per_relation_per_stratum"]
    cards = '<figure><div class="stat-grid">' + ''.join(
        f'<div class="stat-card"><div class="num">{value}</div><div class="label">{label}</div></div>'
        for value, label in ((units, "Unités généalogiques"), (2 * units, "Sujets simulés"), (3, "Profils appariés")))
    cards += '</div><figcaption>Figure 1 — Domaine du banc D21. Source : configuration scellée. Les profils ne sont pas des réplications indépendantes.</figcaption></figure>'
    flow = '<figure><svg viewBox="0 0 680 115" role="img" aria-label="Deux références distinctes pour évaluer 07">'
    for index, (title, detail) in enumerate((("Génotypes observés", "Conformité des calculs 06/07"), ("Pedigree et origines", "Exactitude des parentés"), ("Deux verdicts distincts", "Code conforme ≠ panel qualifié"))):
        horizontal = 8 + index * 224
        flow += f'<rect x="{horizontal}" y="30" width="215" height="65" rx="5" fill="#eaf1f5" stroke="#236c75"/><text x="{horizontal+107}" y="54" text-anchor="middle" font-size="12">{title}</text><text x="{horizontal+107}" y="75" text-anchor="middle" font-size="11">{detail}</text>'
    flow += '</svg><figcaption>Figure 2 — Références complémentaires. Source : préspecification D21. L’évaluateur ouvre la vérité après clôture, sans la transmettre au pipeline.</figcaption></figure>'
    kinship = documents["07"]["evaluations"]
    if stage is None:
        table = _table(("Étape", "Cas mesurés", "Échecs", "Verdict"), [[row["stage_id"], str(row["case_count"]), str(row["fail_count"]), row["status"]] for row in selected], "Verdicts distincts de 06 et 07. Source : manifest D21.")
    elif stage == "06":
        table = _table(("Profil", "Après QC", "Après LD", "Autosomes", "Paires LD", "R² max", "Verdict"), [
            [row["profile"]["profile_id"], str(row["context"]["post_qc_markers"]), str(row["metrics"]["retained_markers"]), str(row["metrics"]["autosome_coverage"]), str(row["metrics"]["residual_ld_pair_count"]), f'{row["metrics"]["residual_ld_max_r2"]:.6f}', row["status"]]
            for row in documents[stage]["evaluations"]], "Panel conservé et LD résiduel dans les fenêtres configurées, pas une estimation génomique globale. Source : CAL-REP-06.json.")
    else:
        table = _table(("Profil", "Conformité", "BA", "IC bootstrap 95 %", "Paires interdites", "Verdict"), [
            [row["profile"]["profile_id"], row["functional_concordance"], _percentage(row["performance"]["balanced_accuracy"]), _interval(row["performance"]["family_bootstrap_interval"]), str(row["metrics"]["forbidden_true_pair_count"]), row["status"]]
            for row in kinship], "BA = moyenne des six rappels ; IC conditionnel au panel sélectionné. Source : CAL-REP-07.json.")
        confusion = []
        for index, row in enumerate(kinship, 1):
            confusion.extend([[str(index), label, *[str(row["performance"]["confusion"][degree][observed]) for observed in (*DEGREES, "MISSING")]] for degree, label in zip(DEGREES, DEGREE_LABELS, strict=True)])
        table += _table(("Profil¹", "Vrai", "Clone", "D1", "D2", "D3", "D4", "NA", "Absent"), confusion,
                        "Confusion complète. ¹1=complet, 2=manque 5 %, 3=densité 25 %. Colonnes=estimé ; NA=non apparenté, pas donnée manquante.").replace("Tableau 1", "Tableau 2")
    intervals = _table(("Vrai degré", "Complet : rappel [IC]", "Manque : rappel [IC]", "Densité : rappel [IC]"), [
        [label, *[f'{row["performance"]["per_degree"][degree]["correct"]}/{row["performance"]["per_degree"][degree]["unit_count"]} [{_interval(row["performance"]["per_degree"][degree]["wilson_interval"])}]' for row in kinship]]
        for degree, label in zip(DEGREES, DEGREE_LABELS, strict=True)], "Rappels par unité indépendante et Wilson 95 %, descriptifs après sélection commune du panel.")
    reasons = sorted({reason for key in ("06", "07") if stage is None or stage == key for row in documents[key]["evaluations"] for reason in row["reasons"]})
    diagnostics = "; ".join(f'{row["profile"]["profile_id"]} : BA oracle avant pruning {_percentage(row["pre_pruning_oracle_performance"]["balanced_accuracy"])}, après KING {_percentage(row["performance"]["balanced_accuracy"])}' for row in kinship)
    return [
        f'<p>{"CAL-DEV-D21" if stage is None else "CAL-REP-" + stage} · version 1.0 · {escape(manifest["created_at"][:10])}. Banc de parenté de recherche, Effet_fondateur2.</p>',
        f'<p>Verdict mesuré : {status}. Les contrôles de calcul et la capacité à exclure toutes les vraies parentés sont deux résultats distincts. Un FAIL reste un résultat évaluable ; il ne doit pas être remplacé par NOT_EVALUABLE.</p>' + cards,
        '<p>06 conserve-t-elle un panel conforme, et 07 retrouve-t-elle les degrés connus et une proposition sans paire réellement apparentée conservée ensemble ?</p>',
        '<p>Autorité : contrats de production, conventions, CAL-DEC-005 et préspecification D21 v1.0.0. Sources primaires : '
        '<a href="https://www.cog-genomics.org/plink/1.9/ld">PLINK 1.9, LD calculée sur les fondateurs</a> ; '
        '<a href="https://www.cog-genomics.org/static/pdf/Manichaikuletal2010.pdf">Manichaikul et al., 2010, KING-robust, équations 9 et 11</a>. Aucune donnée de patient ou référence 1000G n’est utilisée dans ce banc.</p>',
        f'<p>Graines simulation/manque/bootstrap : {config["simulation_seed"]}/{config["missingness_seed"]}/{config["bootstrap_seed"]}. '
        f'{config["markers_per_autosome"]} loci par autosome, 22 autosomes ; {config["bootstrap_replicates"]} tirages. Versions : {escape(" ; ".join(manifest["tool_versions"]))}. Seuils complets dans report_data.json et la configuration empreintée.</p>',
        '<p>Sept relations : clone, parent-enfant, germains, demi-germains, cousins germains, cousins avec une génération de décalage, non-apparentés. Deux strates artificielles de fréquences 0,25/0,65, huit unités par relation et strate. Les origines des copies et le pedigree exact sont privés. La transmission est mendélienne locus par locus, sans recombinaison par segments ni LD biologique.</p>',
        '<p>Trois profils appariés : complet, manque aléatoire 5 %, densité 25 %. Registre neutre, seuls les liens père-enfant observables sont déclarés ; strates et degrés restent privés. Deux clones inter-FID contrôlent la voie cryptique. Les négatifs cibles intra-FID sont effectivement mesurés : une absence dans .kin0 ne vaut jamais UNRELATED. Exécution limitée à 00–03/05–07.</p>' + flow,
        '<p>Oracle 06 : décodage exact des génotypes, MAF des fondateurs, motifs, partition des variants, couverture et LD résiduelle. Le tie-break PLINK n’est pas imposé. Oracle 07 : métriques KING depuis les génotypes, degrés depuis les valeurs publiées, filiation, comptes et proposition gloutonne maximale orientée qualité. Tolérance native KING : 5,1e-5 absolue. Le graphe théorique est contrôlé séparément ; les profils ne sont jamais cumulés comme unités indépendantes.</p>',
        table,
        '<p>CAL-KIN-001 : 22 autosomes. CAL-KIN-003 : zéro vraie paire apparentée conservée ensemble dans la proposition 07, pas dans une cohorte 09 absente. CAL-KIN-002 : performance du pilote REVIEW_REQUIRED, aucun seuil scellé. Les contrôles fonctionnels exigent zéro discordance. Les défauts de couverture ou paires obligatoires absentes ne sont pas masqués.</p>' + (intervals if stage != "06" else '<p>Les contrôles sont exacts sur les sorties exécutées ; aucun intervalle populationnel n’est inventé pour les 22 autosomes synthétiques.</p>'),
        '<p>Motifs : ' + escape(", ".join(reasons) or "aucune discordance") + '.</p><p>' + escape(diagnostics) + '. Le calcul avant pruning est un diagnostic de l’oracle, pas une seconde exécution KING ni une comparaison causale randomisée. Les moyennes φ théorique/réalisé/estimé et les performances par strate sont publiées dans report_data.json. Aucun seuil n’a été ajusté après observation.</p>',
        '<p>Les effectifs sont 16 unités par degré, sauf D1=32, et huit par strate et relation : minimum de qualification de 200 par classe non atteint. Wilson est descriptif ; le bootstrap familial est stratifié par degré et strate, conditionnel au panel commun sélectionné et sans réexécuter son pruning. Il ne couvre ni l’incertitude de sélection ni le domaine biologique réel. Un rappel 100 % avec petit effectif ne prouve pas une absence d’erreur.</p><p>Les fréquences, coordonnées et densités sont artificielles : pas de puce réelle, de LD populationnelle, de phasage, d’IBD explicite ni de datation validés. La diminution du nombre de paires interdites peut aussi refléter davantage de faux positifs et d’exclusions. La robustesse inter-FID n’est couverte ici que par les clones cryptiques, pas tous les degrés. Aucun journal exhaustif d’accès système n’est instrumenté.</p>',
        f'<p>{status} sur le périmètre mesuré seulement. Les exclusions de 07 restent soumises à kinship_exclusion_approval, non approuvées. Aucune cohorte 09 n’est gelée. D21 fournit les mesures ; un échec d’indépendance exige une décision et une étude de sensibilité avant toute qualification. D22 n’est pas exécutée. calibrated=false ; founder_effect_proven=false.</p>',
        '<p>Méthode : ' + escape(manifest["method_id"]) + '.</p><ul>' + ''.join(f'<li>{label} : <span class="hash">{escape(value)}</span></li>' for label, value in (
            ("Évaluation close", manifest_sha), ("Entrées", manifest["panel_manifest_sha256"]), ("Runs clos", manifest["run_index_sha256"]))) + '</ul>',
    ]


def build_reports(*, run_root: Path, template_path: Path, output_root: Path) -> list[Path]:
    """Publie trois rapports dans un dossier neuf du workspace ; PDF après QA."""
    output_root = new_workspace_run(output_root.absolute().parent, output_root.name)
    manifest = checked_document(run_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
    source_sha = _sha(run_root / "evaluation_manifest.json")
    documents = {stage: _json(run_root / f"CAL-REP-{stage}.json") for stage in ("06", "07")}
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d21_report_template_invalid")
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}table{font-size:0.65em}th,td{overflow-wrap:anywhere}main p{orphans:3;widows:3}main section{break-inside:avoid}</style></head>')
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=output_root.parent))
    paths = []
    try:
        for stage in (None, "06", "07"):
            report_id = "CAL-DEV-D21" if stage is None else "CAL-REP-" + stage
            title = report_id + " — " + ("Banc de parenté" if stage is None else "Panel de parenté" if stage == "06" else "Parenté et proposition indépendante")
            toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(
                f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
            body = f'<body lang="fr"><section class="cover"><h1>{escape(title)}</h1><p>Effet_fondateur2 · version 1.0 · {manifest["created_at"][:10]}</p></section>' + toc
            body += f'<footer class="page-footer"><span>{report_id} · Effet_fondateur2</span></footer><main>'
            body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, report_sections(manifest, documents, stage, source_sha), strict=True), 1))
            directory = temporary / report_id
            directory.mkdir(mode=0o700)
            output = directory / "index.html"
            output.write_text(re.sub(r"<title>.*?</title>", f"<title>{escape(title)}</title>", head, flags=re.S) + body + '</main></body></html>', encoding="utf-8")
            write_json(directory / "report_data.json", dict(manifest=manifest, results=documents if stage is None else {stage: documents[stage]}, diagnostic_07=documents["07"] if stage == "06" else None))
            (directory / "checksums.sha256").write_text(''.join(f"{_sha(path)}  {path.name}\n" for path in (output, directory / "report_data.json")), encoding="utf-8")
            paths.append(output_root / report_id / "index.html")
        checked_document(run_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
        if _sha(run_root / "evaluation_manifest.json") != source_sha:
            raise ValueError("d21_report_source_changed")
        for path in temporary.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
        temporary.rename(output_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return paths


def main() -> int:
    """CLI des rapports indépendants 06/07 et développement D21."""
    parser = argparse.ArgumentParser(description="Rapports D21 après évaluation close")
    for name in ("run-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    for path in build_reports(**vars(parser.parse_args())):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
