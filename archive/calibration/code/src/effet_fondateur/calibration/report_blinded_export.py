"""Rapport indépendant CAL-DEV-D18 : contrôle des formats et de l'aveuglement."""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
from typing import Any

from effet_fondateur.calibration.export_blinded_inputs import BlindedExportError, load_checked_export
from effet_fondateur.calibration.materialize_ibd import _sha


def _sections(manifest: dict[str, Any], manifest_sha: str) -> list[tuple[str, str]]:
    """Présente uniquement les agrégats scellés ; jamais la correspondance privée."""
    cards = '<figure><div class="stat-grid">' + "".join(
        f'<div class="stat-card"><div class="num">{value:,}</div><div class="label">{label}</div></div>'.replace(",", " ")
        for value, label in ((manifest["bundle_count"], "Lots aveugles"), (manifest["file_count"], "Fichiers inventoriés"), (0, "Discordance d’export"))
    ) + '</div><figcaption>Figure 1 — Export D18. Source : manifest D18 ; inventaire JSON exclu du nombre de fichiers.</figcaption></figure>'
    table = '<p>Le tableau 1 décrit les observations exportées, pas des réplications biologiques indépendantes.</p>'
    table += '<table><caption>Tableau 1 — Concordance des lots. Source : audits scellés du manifest D18.</caption>'
    table += '<thead><tr><th>Lot opaque</th><th>Sujets</th><th>Sondes</th><th>Appels absents</th><th>Discordances</th></tr></thead><tbody>'
    for item in manifest["audits"]:
        table += f'<tr><th>{escape(item["bundle_id"])}</th>' + "".join(f'<td>{item[field]:,}</td>'.replace(",", " ") for field in ("sample_count", "marker_count", "missing_call_count", "genotype_mismatch_count")) + '</tr>'
    table += '</tbody></table>'
    flow = '<figure><svg viewBox="0 0 680 155" role="img" aria-label="Séparation des entrées observées et de la correspondance privée">'
    for x, y, first, second in ((8, 30, "Observation D17", "GT non phasés seulement"), (345, 30, "Entrées aveugles", "ACPA · VCF · registre · YAML"), (345, 100, "Contrôle privé", "Correspondance · profils · graine")):
        flow += f'<rect x="{x}" y="{y}" width="325" height="48" rx="5" fill="#eaf1f5" stroke="#236c75"/><text x="{x+162}" y="{y+19}" text-anchor="middle" font-size="13">{first}</text><text x="{x+162}" y="{y+37}" text-anchor="middle" font-size="11">{second}</text>'
    flow += '</svg><figcaption>Figure 2 — Deux espaces séparés. Source : contrat D18 ; le pipeline reçoit uniquement les entrées aveugles.</figcaption></figure>'
    hashes = '<ul>' + "".join(f'<li>{label} : <span class="hash">{value}</span></li>' for label, value in (
        ("Manifest D18", manifest_sha), ("Manifest D17 source", manifest["source_manifest_sha256"]),
        ("Configuration D18", manifest["configuration_sha256"]), ("Inventaire aveugle", manifest["inventory_sha256"]),
        ("Correspondance privée", manifest["blinding_map_sha256"]))) + '</ul>'
    return [
        ("Identification", f'<p>CAL-DEV-D18 · version 1.0 · {escape(manifest["created_at"][:10])}. '
         'Responsable : contrôle de développement Effet_fondateur2. Périmètre régional de validation des entrées.</p>'),
        ("Résumé exécutif", '<ul><li>Les huit observations D17 sont exportées sous des noms opaques.</li>'
         '<li>ACPA et VCF non phasés, registre neutre et configuration des étapes 00–02.</li>'
         '<li>Formats et appels vérifiés ; aucune vérité biologique ajoutée aux entrées.</li>'
         '<li>Le pedigree et les informations moléculaires nécessaires à une analyse complète ne sont pas inventés.</li></ul>' + cards),
        ("Question évaluée", '<p>Les génotypes réellement observés peuvent-ils entrer dans les lecteurs du pipeline '
         'sans changement de dosage, perte de manque ni exposition des profils et identifiants de simulation ?</p>'),
        ("Sources et autorité", '<p>AGENTS.md ; conventions de calibration ; protocole D03 §5 et §9 ; plan D18 ; '
         'manifest D17 finalisé ; schémas D18 ; contrats samples_master et pipeline_config ; lecteurs ACPA des étapes 01 et 03.</p>'),
        ("Configuration, graines et versions", f'<p>Méthode {escape(manifest["method_id"])}. La graine de blinding '
         'reste privée. Les pseudonymes et permutations sont déterministes par rôle '
         'et lot ; ils ne dépendent ni du génotype, ni de l’ascendance. Aucun outil génétique externe n’est appelé.</p>'),
        ("Vérité simulée", '<p>Les artefacts D17 sont contrôlés par empreinte avant/après. Seuls les VCF '
         'OBSERVED_MEASUREMENT_PANEL sont décodés. Le registre d’erreurs et les candidats phasés ne sont pas interprétés. '
         'Les tables de correspondance restent privées et permettent la comparaison ultérieure par l’évaluateur.</p>'),
        ("Entrées aveugles", '<p>Les champs autorisés sont reconstruits : aucun en-tête privé n’est recopié. '
         'Les noms de lots, fichiers et individus ne portent pas le profil. L’inventaire opérateur ne contient ni '
         'graine, profil, correspondance, empreinte de vérité ni chemin privé. Voir figure 2.</p>' + flow),
        ("Méthode d’évaluation", '<p>Chaque VCF est relu et comparé cellule par cellule après pseudonymisation. '
         'Chaque ACPA est inspecté par le lecteur 01, puis décodé et comparé aux appels observés. Le registre passe '
         'le validateur 02 ; la configuration passe le schéma du pipeline. Un inventaire fermé refuse tout fichier, '
         'répertoire ou lien supplémentaire. Les tests utilisent aussi le lecteur de conversion 03 et un run synthétique '
         'limité à 00–02, sans approbation automatique.</p>'),
        ("Résultats et incertitude", table),
        ("Critères préspecifiés", '<p>Zéro différence de génotype, de manque, de locus ou d’ordre d’échantillon après '
         'réalignement. Fichiers et inventaires complets, aucune référence privée dans les entrées, aucune mutation '
         'des sources, publication atomique et refus d’écrasement. Les contrôles du tableau 1 sont techniques.</p>'),
        ("Anomalies et écarts", '<p>Un appel absent devient explicitement NoCall en ACPA et reste absent après '
         'décodage. Le lecteur 03 existant compte aussi NoCall dans ses appels invalides : ce compteur n’est pas une '
         'nouvelle erreur de simulation. Aucun rsID absent n’est inventé. Les allèles sont sérialisés comme appels '
         'forward synthétiques depuis le VCF GRCh38 ; ni intensité ni confiance ChAS ne sont simulées.</p>'),
        ("Limites", '<p>Les étapes 03 et suivantes ne sont pas activées dans la configuration livrée. Le registre '
         'reste UNKNOWN pour sexe et clinique ; PID/MID=0 signifie filiation non renseignée, pas absence de parenté. '
         'Chaque FID individuel est une clé technique, pas une unité familiale indépendante. Les génotypes cibles '
         'ne sont pas complétés depuis la vérité, aucun HGVS ni statut moléculaire CONFIRMED n’est fabriqué. '
         'Des métadonnées d’observation autorisées, une validation moléculaire et une décision humaine restent '
         'nécessaires avant l’analyse familiale complète. La région ne remplace pas 22 autosomes. L’aveuglement '
         'est procédural, non une protection contre un opérateur possédant le code de simulation et ses graines.</p>'),
        ("Décision et suite", '<p>PASS technique des formats régionaux et du cloisonnement. Performance scientifique : '
         'NOT_EVALUABLE ; calibrated=false ; founder_effect_proven=false. Aucun pipeline complet n’a été exécuté. '
         'Le manifest global D06/D07 et les portes du pilote restent à sceller avant les runs correspondants. '
         'Prochaine étape : D19, évaluateurs techniques 00–04, en conservant explicitement les entrées manquantes.</p>'),
        ("Provenance et empreintes", hashes),
    ]


def build_report(*, run_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit un rapport agrégé, sans liste de sujets ni correspondance privée."""
    parts = tuple(part.casefold() for part in output_root.absolute().parts)
    if output_root.absolute().parent.name != "calibration_workspace" or output_root.exists() or ".." in parts or "dock6_real" in parts or any(
        left == "data" and right == "input" for left, right in zip(parts, parts[1:])
    ) or any(path.is_symlink() for path in (output_root, *output_root.parents)):
        raise BlindedExportError("d18_report_output_invalid")
    manifest = load_checked_export(run_root)
    source_sha = _sha(run_root / "private_truth/d18/blinded_export_manifest.json")
    sections = _sections(manifest, source_sha)
    title = "D18 — Validation de l’export aveugle"
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    if "</head>" not in head:
        raise BlindedExportError("d18_report_template_invalid")
    head = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", head, flags=re.S)
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}table{font-size:0.83em}th,td{overflow-wrap:anywhere}main p{orphans:3;widows:3}#s9,#s12{break-inside:avoid}</style></head>')
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, (heading, _) in enumerate(sections, 1)) + '</ol></nav>'
    body = f'<body lang="fr"><section class="cover"><h1>{title}</h1><p>Effet_fondateur2 · {escape(manifest["created_at"][:10])} · CAL-DEV-D18 · version 1.0</p></section>' + toc
    body += '<footer class="page-footer"><span>D18 · Effet_fondateur2</span></footer><main>'
    body += "".join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(sections, 1)) + '</main></body></html>'
    output_root.mkdir(parents=True, mode=0o700)
    output = output_root / "index.html"
    output.write_text(head + body, encoding="utf-8")
    (output_root / "report_data.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output_root / "checksums.sha256").write_text("".join(f"{_sha(path)}  {path.name}\n" for path in (output, output_root / "report_data.json")), encoding="utf-8")
    return output


def main() -> int:
    """CLI du rendu D18 ; aucune métrique scientifique recalculée."""
    parser = argparse.ArgumentParser(description="Rapport indépendant D18")
    for name in ("run-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    build_report(**vars(parser.parse_args()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
