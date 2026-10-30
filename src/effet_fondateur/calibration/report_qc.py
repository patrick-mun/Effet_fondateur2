"""Rapports D20 autonomes, un document par étape ; aucun calcul statistique.

Les matrices de confusion proviennent exclusivement des JSON évalués clos.
Le modèle report-formater fournit CSS, Paged.js et pagination A4 embarqués.
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

from effet_fondateur.calibration.materialize_ibd import _json, _sha
from effet_fondateur.calibration.qc_bench import checked_inventory, write_json
from effet_fondateur.calibration.report_technical_evaluation import HEADINGS, _table
from effet_fondateur.calibration.technical_scenarios import new_workspace_run


LABELS = {"sample_missingness": "Manque sujets", "variant_missingness": "Manque variants", "heterozygosity_outlier": "Hétérozygotie",
          "low_maf_for_kinship": "MAF initiale", "batch_missingness_delta": "Écart de lots", "potential_technical_duplicate": "Clone injecté",
          "low_maf_controls": "MAF témoins", "hwe_controls": "HWE témoins", "low_maf_target_region": "MAF régionale", "target_exception": "Exception cible"}
COHORT_LABELS = {"controls_unrelated": "Témoins", "target_carriers_independent": "Porteurs", "target_chromosome_all_qc": "Région"}
UNIT_LABELS = {"sample": "Sujets", "variant": "Variants", "pair": "Paires", "batch": "Lots"}


def _contents(manifest: dict[str, Any], documents: dict[str, Any], stage: str | None, manifest_sha: str) -> list[str]:
    """Compose les quatorze sections depuis les seuls agrégats clos et leurs sources."""
    selected = manifest["stages"] if stage is None else [item for item in manifest["stages"] if item["stage_id"] == stage]
    report_id = "CAL-DEV-D20" if stage is None else "CAL-REP-" + stage
    total = sum(item["case_count"] for item in selected)
    failed = sum(item["fail_count"] for item in selected)
    status = "PASS" if all(item["status"] == "PASS" for item in selected) else "FAIL"
    config = manifest["configuration"]
    cards = '<figure><div class="stat-grid">' + ''.join(
        f'<div class="stat-card"><div class="num">{value}</div><div class="label">{label}</div></div>'
        for value, label in ((total, "Cas × étape mesurés"), (failed, "Cas en défaut"), (22, "Autosomes synthétiques")))
    cards += '</div><figcaption>Figure 1 — Domaine exécuté. Source : manifest d’évaluation D20. Les cas sont appariés, pas des réplications indépendantes.</figcaption></figure>'
    flow = '<figure><svg viewBox="0 0 680 115" role="img" aria-label="Exécution limitée puis évaluation indépendante">'
    for index, (title, detail) in enumerate((("Entrées synthétiques", "IID, défauts préspecifiés"), ("Exécution bornée", "00–05 puis 10 isolée"), ("Évaluation après clôture", "Oracle distinct du pipeline"))):
        horizontal = 8 + index * 224
        flow += f'<rect x="{horizontal}" y="30" width="215" height="65" rx="5" fill="#eaf1f5" stroke="#236c75"/><text x="{horizontal+107}" y="54" text-anchor="middle" font-size="12">{title}</text><text x="{horizontal+107}" y="75" text-anchor="middle" font-size="11">{detail}</text>'
    flow += '</svg><figcaption>Figure 2 — Le producteur des cohortes de 10 est D20_COHORT_CONTRACT_FIXTURE ; aucune exécution 09 n’est prétendue. Source : reçus D20.</figcaption></figure>'
    if stage is None:
        table = _table(("Étape", "Cas mesurés", "Défauts", "Verdict"), [[item["stage_id"], str(item["case_count"]), str(item["fail_count"]), item["status"]] for item in selected], "Contrôles de composants. Source : manifest d’évaluation D20.")
    else:
        rows = []
        for key, values in documents[stage]["rules"].items():
            parts = key.split("/")
            cohort = COHORT_LABELS[parts[0]] if len(parts) == 2 else "Initial"
            rows.append([cohort, LABELS[parts[-1]], UNIT_LABELS[values["unit"]], *[str(values[metric]) for metric in ("tp", "fn", "fp", "tn")]])
        table = _table(("Jeu", "Règle", "Unité", "VP", "FN", "FP", "VN"), rows,
                       "Matrices agrégées par règle et unité. Source : JSON CAL-REP-" + stage + ". VP/FN : positifs détectés/manqués ; FP/VN : négatifs alertés/non alertés.")
    method = "Les fractions sont recalculées depuis les génotypes observés. Les codes, statuts et comptes sont comparés exactement, puis les axes et chaque cellule des triplets retenus sont décodés indépendamment."
    if stage == "05":
        method += " En 05, les exclusions de sujets et variants sont simultanées sur le jeu initial. MAF, hétérozygotie, lots et doublons restent des alertes. Les tables d’alertes et de lots sont contrôlées séparément."
    elif stage == "10":
        method += " En 10, le retrait des sujets précède le recalcul des métriques de variants. HWE et MAF témoins ne sont pas appliqués aux porteurs. La cible régionale est retenue par exception documentée si un filtre la vise. Trois jeux sont comparés séparément."
    reasons = sorted({reason for document in documents.values() for row in document["evaluations"] for reason in row["reasons"]})
    coverage = [key for document in documents.values() for key in document["missing_classes"]]
    return [
        f'<p>{report_id} · version 1.0 · {escape(manifest["created_at"][:10])}. Recherche Effet_fondateur2. Banc autosomique de conformité QC, non qualification biologique.</p>',
        f'<p>Verdict technique : {status}. Les règles visées sont mesurées sur des défauts connus et des contrôles négatifs. Ce verdict ne scelle aucun seuil de performance et ne qualifie pas le pipeline entier.</p>' + cards,
        '<p>Les exclusions, alertes, motifs et génotypes retenus par 05 et 10 correspondent-ils à un oracle indépendant, dans le domaine synthétique préspecifié ?</p>',
        '<p>Autorité : contrats de production, AGENTS.md, conventions de calibration et préspecification D20. Statistiques : '
        '<a href="https://www.cog-genomics.org/plink/1.9/basic_stats">documentation primaire PLINK 1.9</a> ; '
        '<a href="https://www.cog-genomics.org/plink/1.9/ibd">PI_HAT</a> ; '
        '<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC1199378/">Wigginton et al., 2005, HWE exact</a>. Les annotations publiques de la cible sont celles de D19 ; aucun appel expérimental n’est revendiqué.</p>',
        f'<p>Graine {config["seed"]} ; {config["sample_count"]} sujets ; {config["markers_per_autosome"]} marqueurs par autosome ; fréquence allélique initiale {config["allele_frequency"]}. '
        f'PLINK : {escape(" ; ".join(manifest["plink_versions"]))}. Les seuils complets et empreintes sont dans le JSON associé et les reçus. Approbation du seul registre technique : CAL-DEC-004.</p>',
        '<p>Quatorze cas partagent le même fond IID mathématique. Perturbations : manques global, régional et par cohorte, hétérozygotie, lots, clone, allèle rare, HWE témoins et cible rare. L’oracle HWE énumère la distribution conditionnelle exacte sans mid-p. Le clone connu teste la détection, pas la justesse d’un estimateur de parenté.</p>',
        '<p>ACPA, VCF, métadonnées neutres, appels cibles synthétiques explicites et fixtures publiques de cohortes sont distincts de la vérité privée. Le processus pipeline ne reçoit aucun chemin de vérité. Les cohortes sont préspecifiées avant exécution, pas déduites des exclusions observées.</p>' + flow,
        '<p>' + method + '</p><p>Tolérance des seules mesures arrondies : 5e-4 relative, 5e-6 absolue ; tolérance absolue nulle pour HWE_P. Aucun assouplissement des décisions ou des génotypes. Les seuils ne sont pas ajustés après observation.</p>',
        table + '<p>Les matrices mesurent la conformité aux décisions de l’oracle ; seule la ligne clone utilise directement une vérité d’injection. Les tests ciblés couvrent les défauts injectés. Des alertes de variation de fond ne sont pas automatiquement des défauts du code. Les taux descriptifs sont conservés dans le JSON ; aucun intervalle populationnel n’est estimé sur des stress appariés.</p>',
        '<p>Réussite technique : aucune discordance de décision, motif, compte ou génotype, et au moins un positif et un négatif par règle agrégée. Les matrices restent séparées par unité et par cohorte. Une sortie manquante ou altérée bloque l’évaluation ; elle n’est jamais comptée comme réussite.</p>',
        '<p>Motifs observés : ' + escape(", ".join(reasons) or "aucune discordance fonctionnelle détectée") + '. Classes manquantes : ' + escape(", ".join(coverage) or "aucune") + '.</p>'
        '<p>La relecture D20 a identifié une perte de précision dans la sérialisation de HWE_P : le formateur de 10 écrasait les valeurs de module inférieur à 5e-16. La correction normalise seulement le vrai zéro signé, sans modifier le test ni les exclusions. Source : préspecification D20 v1.0.1 et test de régression. Les verdicts ci-dessus concernent uniquement le run empreinté ; les anciens diagnostics ne sont jamais réécrits.</p>',
        '<p>Ce banc ne simule ni LD, ni population 1000G, ni densité réelle de puce, ni IBD ou datation. Les deux lots sont artificiels. Dans le cas cible rare, la cohorte de porteurs singleton teste seulement le comportement mono-lot ; la détection entre lots est couverte par les autres profils. Les contrôles sexuels restent NOT_EVALUATED_BY_AUTOSOMAL_CONTRACT. HWE en 05 et chez les porteurs en 10 est exclu par politique. L’intégration de 10 est isolée : 06–09 ne sont pas exécutées ni validées. Aucun journal système exhaustif d’accès n’est revendiqué.</p>',
        f'<p>{status} de composant seulement. Décision de performance scientifique : REVIEW_REQUIRED, seuils du pilote non scellés. calibrated=false ; founder_effect_proven=false. Prochaine étape : D21, parenté 06/07, avant l’intégration réelle des cohortes de D22.</p>',
        '<p>Méthode : ' + escape(manifest["method_id"]) + '.</p><ul>' + ''.join(
            f'<li>{label} : <span class="hash">{escape(value)}</span></li>' for label, value in (
                ("Évaluation", manifest_sha), ("Panel", manifest["panel_manifest_sha256"]), ("Reçu des exécutions", manifest["run_index_sha256"]))) + '</ul>',
    ]


def build_reports(*, run_root: Path, template_path: Path, output_root: Path) -> list[Path]:
    """Publie trois HTML/JSON sans écraser un rendu ; PDF après contrôle Chrome."""
    output_root = new_workspace_run(output_root.absolute().parent, output_root.name)
    manifest = checked_inventory(run_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
    source_sha = _sha(run_root / "evaluation_manifest.json")
    documents = {stage: _json(run_root / f"CAL-REP-{stage}.json") for stage in ("05", "10")}
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d20_report_template_invalid")
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}table{font-size:0.73em}th,td{overflow-wrap:anywhere}main p{orphans:3;widows:3}main section{break-inside:avoid}</style></head>')
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=output_root.parent))
    paths = []
    try:
        for stage in (None, "05", "10"):
            report_id = "CAL-DEV-D20" if stage is None else "CAL-REP-" + stage
            title = report_id + " — " + ("Banc de conformité QC" if stage is None else "QC préliminaire" if stage == "05" else "QC final par cohorte")
            toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(
                f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
                for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
            contents = _contents(manifest, documents, stage, source_sha)
            body = f'<body lang="fr"><section class="cover"><h1>{escape(title)}</h1><p>Effet_fondateur2 · {report_id} · version 1.0 · {manifest["created_at"][:10]}</p></section>' + toc
            body += f'<footer class="page-footer"><span>{report_id} · Effet_fondateur2</span></footer><main>'
            body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, contents, strict=True), 1))
            body += '</main></body></html>'
            directory = temporary / report_id
            directory.mkdir(mode=0o700)
            output = directory / "index.html"
            output.write_text(re.sub(r"<title>.*?</title>", f"<title>{escape(title)}</title>", head, flags=re.S) + body, encoding="utf-8")
            write_json(directory / "report_data.json", dict(manifest=manifest, result=documents[stage] if stage else None))
            (directory / "checksums.sha256").write_text("".join(f"{_sha(path)}  {path.name}\n" for path in (output, directory / "report_data.json")), encoding="utf-8")
            paths.append(output_root / report_id / "index.html")
        checked_inventory(run_root, "evaluation_manifest.json", "evaluation_state", "CLOSED")
        if _sha(run_root / "evaluation_manifest.json") != source_sha:
            raise ValueError("d20_report_source_changed")
        temporary.rename(output_root)
    except Exception:
        shutil.rmtree(temporary)
        raise
    return paths


def main() -> int:
    """CLI des trois rapports français autonomes, sans appel génétique ni PDF implicite."""
    parser = argparse.ArgumentParser(description="Rapports D20 indépendants 05/10 et développement")
    for name in ("run-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    print(json.dumps(dict(reports=[str(path) for path in build_reports(**vars(parser.parse_args()))])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
