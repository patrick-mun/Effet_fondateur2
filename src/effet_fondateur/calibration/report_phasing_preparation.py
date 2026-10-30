"""Rapport D23a agrégé depuis une préparation close, sans recalcul scientifique."""

from __future__ import annotations

import argparse
from html import escape
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.phasing_preparation import verify_preparation, write_json
from effet_fondateur.calibration.report_technical_evaluation import HEADINGS, _table
from effet_fondateur.calibration.technical_scenarios import new_workspace_run


def column_chart(rows: list[dict[str, Any]], field: str, title: str, axis: str) -> str:
    """Barres exactes par chromosome ; aucun profil ni identifiant individuel."""
    maximum = max(max(row[field] for row in rows), 1)
    scale = 170 / maximum
    content = f'<text x="12" y="20">{escape(title)}</text><text x="12" y="43" font-size="11">{escape(axis)}</text>'
    for fraction in (0, 0.5, 1):
        y = 240 - 170 * fraction
        content += f'<line x1="60" x2="665" y1="{y}" y2="{y}" stroke="#dce3e8"/><text x="54" y="{y+4}" text-anchor="end" font-size="10">{maximum*fraction:g}</text>'
    for index, row in enumerate(rows):
        x, height = 67 + index * 27, row[field] * scale
        content += f'<rect x="{x}" y="{240-height}" width="17" height="{height}" fill="#286d91"/><text x="{x+8.5}" y="{230-height}" text-anchor="middle" font-size="9">{row[field]}</text><text x="{x+8.5}" y="258" text-anchor="middle" font-size="10">{row["chromosome"]}</text>'
    content += '<text x="350" y="281" text-anchor="middle" font-size="11">Autosome</text>'
    return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 295" role="img" aria-label="' + escape(title, quote=True) + '">' + content + '</svg>'


def figure(number: int, svg: str, caption: str) -> str:
    return f'<figure class="scientific-figure">{svg}<figcaption>Figure {number} — {escape(caption)}</figcaption></figure>'


def sections(data: dict[str, Any], source_sha: str) -> list[str]:
    """Explique la préparation et ses limites sans annoncer de performance D23."""
    counts, audit = data["counts"], data["audit"]
    configuration = data["configuration"]
    plots = figure(1, column_chart(data["chromosomes"], "retained_chip_marker_count", "Quels marqueurs de puce sont disponibles ?", "Nombre de SNP retenus"),
                   "SNP de la puce harmonisés et retenus par autosome. La cible synthétique supplémentaire n’est pas comptée comme sonde de puce. Source : preparation_manifest.json, chromosomes.")
    plots += figure(2, column_chart(data["chromosomes"], "rare_heterozygote_count", "Où pourra-t-on mesurer la phase des variants rares ?", "Génotypes hétérozygotes, 0 < MAF ≤ 5 %"),
                    "Opportunités de mesure, pas erreurs de phasage. MAF = fréquence de l’allèle minoritaire dans les 40 sujets simulés, enfants compris ; les observations apparentées ne sont pas indépendantes. Source : preparation_manifest.json, chromosomes.")
    total_rare = sum(row["rare_heterozygote_count"] for row in data["chromosomes"])
    source_table = _table(("Partition D10", "Effectif", "Utilisation D23a"), [
        ["Donneurs calibration", str(counts["calibration_donors"]), str(counts["source_donor_groups"]) + " groupes contributeurs"],
        ["Donneurs validation", str(counts["validation_donors"]), "Aucun génotype utilisé"],
        ["Référence de phasage", str(counts["reference_samples"]), "Liste exacte imposée, 22 extraits vérifiés"]],
        "Partitions sans intersection. Source : manifest D10 vérifié et préparation D23a.")
    labels = dict(genotype_mismatch_count="Allèles VCF observés", acpa_mismatch_count="Appels ACPA",
                  founder_copy_mismatch_count="Copies fondatrices", transmission_mismatch_count="Copies transmises",
                  target_mismatch_count="Génotypes cibles", pedigree_mismatch_count="Filiations déclarées")
    audit_table = _table(("Contrôle exhaustif", "Discordances"), [[labels[key], str(value)]
                        for key, value in audit.items() if key.endswith("_mismatch_count")],
                        "Relecture indépendante des fichiers exportés. Source : preparation_manifest.json, audit.").replace("Tableau 1", "Tableau 2")
    selection_math = '<div class="equation"><p><strong>Équation 1 — Séparation et génotype observé</strong></p><div class="math" role="math">R ∩ D<sub>cal</sub> = R ∩ D<sub>val</sub> = D<sub>cal</sub> ∩ D<sub>val</sub> = ∅<br>G<sub>ij</sub> = H<sub>ij1</sub> + H<sub>ij2</sub></div><p>R : références D10 ; D : donneurs des deux partitions. H ∈ {0,1} code REF/ALT ; G ∈ {0,1,2} est le dosage ALT, sans unité. L’ordre des deux copies n’est jamais exporté comme observation. Le contrôle porte sur les donneurs d’origine, pas sur le seul nom synthétique.</p></div>'
    meiosis_math = '<div class="equation"><p><strong>Équation 2 — Recombinaison et transmission</strong></p><div class="math" role="math">N ∼ Poisson(L), L = (cM<sub>fin</sub> − cM<sub>début</sub>)/100<br>b<sub>k</sub> | N ∼ Uniforme(cM<sub>début</sub>, cM<sub>fin</sub>)<br>T<sub>j</sub> = (B + ∑<sub>k</sub> 𝟙[b<sub>k</sub> ≤ cM<sub>j</sub>]) mod 2<br>H<sub>enfant,j,père</sub> = H<sub>père,j,Tj+1</sub></div><p>L est en Morgans ; B est la copie initiale 0/1 tirée avec probabilité 1/2 ; T vaut 0/1 et indexe les copies H1/H2 par T+1 ; b sont les positions de recombinaison en cM. Même règle, sous-flux distinct, pour la mère. Modèle de Haldane sans interférence, pas de carte sexuée ni de recombinaison en bp uniforme.</p></div>'
    return [
        '<p>CAL-DEV-D23A — Préparation autosomique et référence de phasage. Effet_fondateur2, ' + escape(data["created_at"][:10]) + ', version 1.0.0. Rapport de développement, pas rapport de performance 12/12A.</p>',
        f'<p>Préparation close : {counts["samples"]} sujets, 22 autosomes, {counts["markers"]} marqueurs observés ; référence limitée aux {counts["reference_samples"]} sujets D10 autorisés. Les contrôles exhaustifs d’export ne trouvent aucune discordance. SHAPEIT5 et le pipeline complet n’ont pas été exécutés.</p><aside class="callout">Prêt pour développer les évaluateurs D23 ; aucune précision de phasage ni qualification scientifique n’est encore démontrée.</aside>',
        '<p>Les entrées couvrent-elles tous les autosomes aux marqueurs de la puce ? La référence exclut-elle les donneurs des deux panels malgré le renommage synthétique ? Pedigrees et génotypes moléculaires sont-ils explicites, sans fuite de phase vraie ?</p>',
        '<p>Autorité : CAL-DEC-007, préspecification D23a, conventions et AGENTS.md. Données : panel public IGSR 1000G 20220422 GRCh38, annotation officielle CytoScan 750K Accel NA36.r1, archive HapMap/SHAPEIT épinglée D09. Le catalogue D10 historique est conservé à son SHA exact ; la version ultérieure ajoutant la puce ne réécrit pas cet engagement.</p>',
        f'<p>Graines fondateurs/transmissions : {configuration["founder_seed"]}/{configuration["meiosis_seed"]}, privées ; graine distincte du futur phasage : {configuration["phasing_seed"]}. Huit trios et seize sujets supplémentaires. Sélection physique bornée à {configuration["maximum_markers_per_autosome"]} SNP par autosome, plus les sondes dans la fenêtre cible ±{configuration["target_flank_bp"] / 1e6:g} Mb. Minimum préspecifié : {configuration["minimum_markers_per_autosome"]} par autosome et {configuration["minimum_total_markers"]} au total. Outil réel : {escape(" ; ".join(data["tool_versions"]))} ; SHAPEIT5 non lancé.</p>',
        f'<p>{counts["founders"]} fondateurs issus de {counts["source_donor_groups"]} groupes donneurs D10 distincts, une copie par donneur. Les enfants héritent des copies parentales réellement recombinées. Il y a {counts["independent_simulation_units"]} unités de simulation (familles ou sujets isolés), pas 40 observations indépendantes.</p>' + source_table,
        '<p>Observations non phasées ACPA/VCF, registre de pedigree synthétique, génotypes moléculaires explicites, annotation cible et sélection de référence. Clinique UNKNOWN, groupe UNSPECIFIED. Phase, origines donneurs et recombinaisons restent dans private_truth. La validation aveugle n’est pas matérialisée ni inspectée. Le profil fourni s’arrête à 00–02, sans approbation scientifique automatique ; les paramètres de sélection 12/12A sont fournis séparément.</p>',
        '<p>Copies fondatrices contrôlées contre les appels donneurs conservés ; transmissions rejouées par sous-graines ; VCF observé, ACPA, cartes, pedigree et génotypes cibles relus. Les extraits de référence sont contrôlés sur sujets, ordre, identité des variants et phase complète. L’intégrité est ancrée aux empreintes des 22 VCF/index publics et de tous les artefacts générés.</p>' + selection_math + meiosis_math,
        '<p>La figure 1 montre la couverture en marqueurs de puce ; la figure 2 distingue la disponibilité de génotypes rares d’une précision encore inconnue.</p>' + plots + f'<p>Total : {total_rare} génotypes hétérozygotes rares disponibles pour un futur comparateur. Ces dénombrements exhaustifs n’ont pas d’intervalle d’échantillonnage de performance.</p>' + audit_table,
        '<p>Réussite de préparation : 22 autosomes, seuils de couverture atteints, zéro intersection interdite D10, référence exactement sélectionnée, cible unique et zéro discordance des exports. Les seuils de sensibilité/commutation de D23 ne sont pas scellés et ne sont pas évalués ici.</p>',
        '<p>Le contrôle a refusé le catalogue actuel incompatible avec D10 et l’archive locale de cartes non conforme. Le catalogue D10 exact et la carte publique épinglée ont été utilisés sans modifier les sources historiques.</p><p>Le brouillon v1 proposait la graine des méioses pour le futur phasage. V2 impose une graine distincte, conserve les SHA-256 des artefacts génétiques et relit les exports ; seuls trois fichiers de paramètres/chemins changent. V1 est conservé, mais interdit d’usage.</p>',
        '<p>Puce sous-échantillonnée, petit effectif : aucune qualification de densité, KING, IBD, datation ou effet fondateur. D10 sépare les filiations documentées, sans exclure tout apparentement cryptique. Haldane ignore l’interférence et les différences de carte entre sexes.</p><p>Cible APC chr5 : annotation D19, génotypes injectés, pas de mesure expérimentale. La cible supplémentaire est déclarée, pas présentée comme une sonde de puce. Aucun run DOCK6.</p>',
        '<p>D23 reste à faire : quatre évaluateurs 10A/11/12/12A, contrôles positifs/négatifs, intégration bornée et rapports indépendants. Comparer modulo permutation H1/H2, ancrer la cible indépendamment, stratifier fréquence/densité/manque et quantifier l’incertitude par famille. FAIL 07 toujours ouvert ; calibrated=false et founder_effect_proven=false.</p>',
        '<p>Source close : preparation_manifest.json, SHA-256 <span class="hash">' + escape(source_sha) + '</span>. Report_data.json contient seulement les agrégats et empreintes, jamais les listes de sujets, correspondances donneurs ou génotypes. Les figures reprennent directement ces agrégats, sans nouveau calcul de performance.</p>',
    ]


def build_report(*, panel_root: Path, template_path: Path, output_root: Path) -> Path:
    """Publie le rapport indépendant, sans modifier les ressources closes."""
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    manifest = verify_preparation(panel_root)
    source_sha = sha256_file(panel_root / "preparation_manifest.json")
    fields = ("schema_version", "method_id", "created_at", "configuration", "counts", "chromosomes", "audit", "tool_versions",
              "partition_manifest_sha256", "reference_selection_sha256", "source_fingerprints", "chip_package_sha256",
              "chip_marker_table_sha256", "map_archive_sha256", "phasing_executed", "performance_status", "calibrated", "founder_effect_proven")
    data = {key: manifest[key] for key in fields}
    data["preparation_manifest_sha256"] = source_sha
    head = re.split(r"<body\b", template_path.read_text(), maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d23a_template_invalid")
    title = "CAL-DEV-D23A — Préparation du phasage autosomique"
    head = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", head, flags=re.S)
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}main section{break-inside:avoid}main #s8,main #s9{break-inside:auto}.scientific-figure,.equation{break-inside:avoid}.scientific-figure{margin:1em 0}table{font-size:.78em}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.equation p{font-size:.85em}.math{font:1.05em Georgia,serif;line-height:1.8}main p{orphans:3;widows:3}</style></head>')
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
    body = f'<body lang="fr"><section class="cover"><h1>{title}</h1><p>Effet_fondateur2 · {escape(manifest["created_at"][:10])} · version 1.0.0</p><p>Préparation contrôlée, graphiques de disponibilité et modèle mathématique</p></section>' + toc
    body += '<footer class="page-footer"><span>CAL-DEV-D23A · Effet_fondateur2</span></footer><main>'
    body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, sections(data, source_sha), strict=True), 1))
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        directory = temporary / "CAL-DEV-D23A"
        directory.mkdir(mode=0o700)
        html = directory / "index.html"
        html.write_text(head + body + '</main></body></html>')
        write_json(directory / "report_data.json", data)
        (directory / "checksums.sha256").write_text(''.join(f'{sha256_file(path)}  {path.name}\n' for path in (html, directory / "report_data.json")))
        verify_preparation(panel_root)
        if sha256_file(panel_root / "preparation_manifest.json") != source_sha:
            raise ValueError("d23a_report_source_changed")
        temporary.rename(destination)
        return destination / "CAL-DEV-D23A/index.html"
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("panel-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    build_report(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()
