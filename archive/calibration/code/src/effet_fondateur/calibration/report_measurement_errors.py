"""Rendu CAL-DEV-D17 autonome ; agrégats calculés par le générateur seulement."""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
from typing import Any

from effet_fondateur.calibration.materialize_ibd import _sha
from effet_fondateur.calibration.materialize_measurement_errors import load_checked_run
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementError


def _sections(manifest: dict[str, Any], source_sha: str) -> list[tuple[str, str]]:
    """Mappe les résultats scellés sur les quatorze rubriques du plan D17."""
    profiles = manifest["profiles"]
    cards = '<figure><div class="stat-grid">' + "".join(
        f'<div class="stat-card"><div class="num">{value}</div><div class="label">{label}</div></div>'
        for value, label in ((len(profiles), "Profils appariés"), (manifest["sample_count"], "Sujets synthétiques"),
                             (manifest["source_marker_count"], "Sondes initiales"))
    ) + '</div><figcaption>Figure 1 — Périmètre D17. Source : manifest D17.</figcaption></figure>'
    table = '<p>Le tableau 1 distingue les quotas d’injection ; ce ne sont pas des performances de détection.</p>'
    table += '<table><caption>Tableau 1 — Événements enregistrés. Source : profils du manifest D17.</caption>'
    table += '<thead><tr><th>Profil</th><th>Sondes</th><th>Manque</th><th>Dosage</th><th>Phase ponctuelle</th><th>Frontières</th></tr></thead><tbody>'
    for profile in profiles:
        counts = profile["event_counts"]
        values = [profile["retained_marker_count"], counts["MISSING_GT"], counts["GENOTYPE"],
                  counts["PHASE_POINT"], counts["PHASE_SWITCH_BOUNDARY"]]
        table += f'<tr><th>{escape(profile["profile_id"])}</th>' + "".join(f'<td>{value:,}</td>'.replace(",", " ") for value in values) + '</tr>'
    table += '</tbody></table>'
    switch_table = '<table><caption>Tableau 2 — Frontières et appels inversés. Source : manifest D17.</caption>'
    switch_table += '<thead><tr><th>Profil</th><th>Frontières éligibles</th><th>Frontières injectées</th><th>GT inversés par blocs</th></tr></thead><tbody>'
    for profile in profiles:
        if profile["requested"]["phase_switch_rate"]:
            switch_table += f'<tr><th>{escape(profile["profile_id"])}</th>' + "".join(
                f'<td>{value:,}</td>'.replace(",", " ") for value in (
                    profile["eligible_switch_boundary_count"], profile["event_counts"]["PHASE_SWITCH_BOUNDARY"],
                    profile["event_counts"]["PHASE_SWITCH_GT"])) + '</tr>'
    switch_table += '</tbody></table>'
    flow = '<figure><svg viewBox="0 0 680 100" role="img" aria-label="Séparation vérité observation et phase candidate">'
    for index, (first, second) in enumerate((("Vérité D16", "inchangée"), ("Profils D17", "registre privé"), ("GT non phasés", "futur export D18"), ("Phase candidate", "évaluateur seulement"))):
        x = 8 + index * 169
        flow += f'<rect x="{x}" y="30" width="155" height="52" rx="5" fill="#eaf1f5" stroke="#236c75"/><text x="{x+77}" y="51" text-anchor="middle" font-size="13">{first}</text><text x="{x+77}" y="70" text-anchor="middle" font-size="11">{second}</text>'
    flow += '</svg><figcaption>Figure 2 — Deux observations distinctes de la vérité. Source : contrat D17 ; les deux derniers blocs sont des sorties parallèles.</figcaption></figure>'
    provenance = '<p>SHA-256 du manifest D17 :</p><p class="hash">' + source_sha + '</p>'
    provenance += '<p>Manifest D16 : <span class="hash">' + manifest["source"]["manifest_sha256"] + '</span></p>'
    provenance += '<p>Configuration : <span class="hash">' + manifest["configuration"]["sha256"] + '</span></p>'
    seeds = '<ul>' + "".join(f'<li>{escape(item["profile_id"])} : <span class="hash">' + escape(json.dumps(item["seeds"], sort_keys=True)) + '</span></li>' for item in profiles) + '</ul>'
    return [
        ("Identification", f'<p>CAL-DEV-D17 · version 1.0 · {escape(manifest["created_at"][:10])}. '
         'Responsable : contrôle de développement Effet_fondateur2. Analyse régionale du panel de calibration.</p>'),
        ("Résumé exécutif", '<ul><li>Un contrôle intact et des variantes de mesure appariées au même run D16.</li>'
         '<li>Chaque défaut et chaque appel inversé possèdent un registre avant/après.</li>'
         '<li>Rejeu indépendant et relecture exhaustive des deux exports : PASS technique.</li>'
         '<li>Aucune qualification scientifique ni modèle empirique des erreurs de puce.</li></ul>' + cards),
        ("Question évaluée", '<p>Peut-on dégrader l’observation de façon déterministe et intégralement traçable, '
         'tout en conservant les origines biologiques IBD (identité par descendance), ROH (régions d’homozygotie), '
         'les pedigrees et les transmissions vrais ?</p>'),
        ("Sources et autorité", '<p>AGENTS.md ; conventions de calibration ; protocole D03 §8.5 ; '
         'plan D17 ; trois schémas measurement_error ; configuration D17 et manifest D16 liés par SHA-256. '
         'Les chiffres des tableaux 1 et 2 et de la figure 1 proviennent exclusivement du manifest D17.</p>'),
        ("Configuration, graines et versions", f'<p>Méthode {escape(manifest["method_id"])} ; graine measurement_error '
         f'{manifest["configuration"]["measurement_error_seed"]}. Sous-flux SHA-256 séparés par rôle, identifiant '
         'de profil et opération, indépendants de l’ordre des profils. Python '
         f'{escape(manifest["software"]["python"])} ; NumPy {escape(manifest["software"]["numpy"])}.</p>' + seeds),
        ("Vérité simulée", '<p>La couche dense D16 et ses huit artefacts privés sont vérifiés avant et après. '
         'La vérité phasée est projetée sans imputation sur les sondes ; son dosage doit correspondre au VCF D16. '
         'Les profils sont des répétitions de mesure d’un même panel, jamais des réplications biologiques indépendantes.</p>'),
        ("Entrées aveugles", '<p>Le VCF observable contient seulement GT non phasé, sans graine, profil ou vérité dans '
         'les champs. Le candidat phasé reste privé et réservé aux évaluateurs : il ne remplace ni la phase vraie '
         'ni une sortie SHAPEIT. L’export final et ses métadonnées restent à réaliser en D18.</p>' + flow),
        ("Méthode d’évaluation", '<p>Ordre fixé : densité → erreur de dosage → commutation → point de phase → manque. '
         'Quota = taux × effectif, arrondi au plus proche avec demi supérieur, tirage uniforme sans remise. '
         'Dosage et manque : cellules retenues ; phase ponctuelle : hétérozygotes après erreurs de dosage ; '
         'commutations : frontières entre hétérozygotes successifs d’un même sujet/chromosome. '
         'Une commutation inverse tous les hétérozygotes du bloc jusqu’à la suivante, sans changement du dosage. '
         'Le registre seul reconstruit ensuite tous les GT et retraits ; les fichiers sont relus et comparés cellule par cellule.</p>'),
        ("Résultats et incertitude", table + '<p>Les frontières ne sont pas des nombres de GT inversés : voir tableau 2. '
         'La figure 2 rappelle la séparation des exports.</p>' + switch_table),
        ("Critères préspecifiés", '<p>Zéro changement non enregistré, zéro divergence de rejeu/export, contrôle intact, '
         'quotas exacts, aucune mutation de la vérité biologique et inventaire empreinté complet. '
         'Les compteurs de comparaison des exports sont nuls pour chaque profil ; les seuils du pipeline ne sont pas ajustés.</p>'),
        ("Anomalies et écarts", '<p>La cible est conservée si elle appartient aux sondes initiales ; aucune sonde '
         'absente n’est ajoutée. Son génotype reste exposé aux défauts comme les autres. Les pertes de densité '
         'ne déplacent pas les segments vrais et ne déclenchent pas d’abaissement des seuils IBD/ROH. '
         'Le manque ou une erreur ponctuelle peut masquer ou compenser une erreur antérieure : le registre conserve toutes les opérations.</p>'),
        ("Limites", '<p>Stress technique régional, défauts uniformes indépendants de l’ascendance et de la qualité : '
         'ni intensités de sonde, ni erreurs corrélées de lot, ni calibration empirique de puce. '
         'Les points et frontières de phase représentent des mécanismes injectés, pas le taux d’erreur mesuré '
         'd’un logiciel de phasage. Les erreurs de carte/bornes ne sont pas étudiées ici. '
         'Les critères après dégradation, le genome-wide, les effectifs et les évaluateurs restent à qualifier.</p>'),
        ("Décision et suite", '<p>PASS technique du générateur. Performances du pipeline et intervalles : NOT_EVALUABLE '
         'avant appels réels et évaluateurs ; aucune inférence de sensibilité/spécificité. '
         'calibrated=false ; founder_effect_proven=false ; gamma_inputs_modified=false. Prochaine étape : D18.</p>'),
        ("Provenance et empreintes", provenance),
    ]


def build_report(*, run_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit un HTML autoportant et ses agrégats dans un dossier neuf du workspace."""
    if output_root.absolute().parent.name != "calibration_workspace" or output_root.exists() or any(
        path.is_symlink() for path in (output_root, *output_root.parents)
    ) or ".." in output_root.parts:
        raise MeasurementError("d17_report_output_invalid")
    parts = tuple(part.casefold() for part in output_root.absolute().parts)
    if "dock6_real" in parts or any(left == "data" and right == "input" for left, right in zip(parts, parts[1:])):
        raise MeasurementError("d17_report_output_invalid")
    manifest = load_checked_run(run_root)
    source_sha = _sha(run_root / "private_truth/d17/measurement_error_manifest.json")
    sections = _sections(manifest, source_sha)
    title = "D17 — Validation des imperfections de mesure"
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    if "</head>" not in head:
        raise MeasurementError("d17_report_template_invalid")
    head = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", head, flags=re.S)
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere} table{font-size:0.83em} th,td{overflow-wrap:anywhere}</style></head>')
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span>'
        '<span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, (heading, _) in enumerate(sections, 1)) + '</ol></nav>'
    body = f'<body lang="fr"><section class="cover"><h1>{title}</h1><p>Effet_fondateur2 · '
    body += escape(manifest["created_at"][:10]) + ' · CAL-DEV-D17 · version 1.0</p><p>Contrôle technique régional ; rapport indépendant de développement.</p></section>' + toc
    body += '<footer class="page-footer"><span>D17 · Effet_fondateur2</span></footer><main>'
    body += "".join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(sections, 1))
    body += '</main></body></html>'
    output_root.mkdir(parents=True, mode=0o700)
    output = output_root / "index.html"
    output.write_text(head + body, encoding="utf-8")
    (output_root / "report_data.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output_root / "checksums.sha256").write_text("".join(f"{_sha(path)}  {path.name}\n" for path in (output, output_root / "report_data.json")), encoding="utf-8")
    return output


def main() -> int:
    """CLI de rendu sans réestimation scientifique ni lecture individuelle."""
    parser = argparse.ArgumentParser(description="Rapport de validation D17")
    for name in ("run-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    build_report(**vars(parser.parse_args()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
