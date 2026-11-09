"""Rapport indépendant CAL-DEV-D16 depuis les agrégats scellés du run ROH."""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
from typing import Any

from effet_fondateur.calibration.materialize_ibd import _artifact, _json, _sha
from effet_fondateur.calibration.materialize_roh import _validate_schema


def _sections(
    manifest: dict[str, Any], source_sha: str, truth: list[dict[str, Any]],
) -> list[tuple[str, str]]:
    """Présente les contrôles calculés par le générateur, sans métrique d'appel."""
    integrity = manifest["integrity"]
    lengths = " ; ".join(f'{item["length_bp"]:,} bp'.replace(",", " ")
                         for item in truth if item["is_positive"])
    labels = {
        "dense_variant_count": "Variants denses", "observed_variant_count": "Marqueurs observables",
        "changed_h2_allele_count": "Substitutions H2", "h1_modified_count": "Modifications H1",
        "roh_copy_error_count": "Erreurs de copie ROH", "ibd_copy_error_count": "Erreurs de copie IBD",
        "mendelian_error_count": "Erreurs mendéliennes", "negative_modified_count": "Négatifs modifiés",
        "induced_ibd_copy_error_count": "Erreurs IBD induit sur H2",
    }
    table = '<table><caption>Tableau 1 — Contrôles techniques. Source : manifest D16.</caption><tbody>' + "".join(
        f'<tr><th>{labels[key]}</th><td>{value}</td></tr>' for key, value in integrity.items()
    ) + '</tbody></table>'
    cards = '<figure><div class="stat-grid">' + "".join(
        f'<div class="stat-card"><div class="num">{value}</div><div class="label">{label}</div></div>'
        for value, label in ((manifest["positive_roh_count"], "ROH injectées"),
                             (manifest["negative_roh_count"], "Témoin sans injection"),
                             (manifest["ibd_coupled_roh_count"], "ROH couplée à l’IBD"))
    ) + '</div><figcaption>Figure 1 — Scénarios D16. Source : manifest D16.</figcaption></figure>'
    bars = '<figure><svg viewBox="0 0 650 160" role="img" aria-label="Marqueurs dans les ROH injectées">'
    for index, item in enumerate(manifest["observability"]):
        y = 30 + index * 40
        label = escape(item["scenario_id"].replace("roh_", ""))
        count = item["observed_marker_count"]
        bars += f'<text x="8" y="{y+14}" font-size="11">{label}</text><rect x="240" y="{y}" width="{min(count*2,320)}" height="20" fill="#236c75"/><text x="{250+min(count*2,320)}" y="{y+14}">{count}</text>'
    bars += '</svg><figcaption>Figure 2 — Marqueurs observables. Source : audit de densité D16.</figcaption></figure>'
    provenance = '<p>Manifest D16 :</p><p style="overflow-wrap:anywhere">' + source_sha + '</p><ul>' + "".join(
        f'<li>{escape(key)} : <code style="overflow-wrap:anywhere">{value}</code></li>'
        for key, value in manifest["inputs"].items()
    ) + '</ul>'
    return [
        ("Identification", f'<p>CAL-DEV-D16 · version 1.0 · {escape(manifest["created_at"])}. '
         f'Run {escape(manifest["roh_simulation_manifest_id"])}.</p>'),
        ("Résumé exécutif", '<ul><li>Trois intervalles autozygotes synthétiques et un témoin sans injection.</li>'
         '<li>Copies H1, vérité IBD et transmissions conservées.</li><li>Validation technique régionale ; performances de détection à mesurer en D26.</li></ul>'),
        ("Question évaluée", '<p>Peut-on injecter une origine haplotypique identique sur les deux copies d’un sujet, '
         'tout en distinguant cette autozygotie du partage IBD entre deux sujets ?</p>'),
        ("Sources et autorité", '<p>AGENTS.md ; protocole D03 §8.4 ; plan D16 ; contrat docs/modules/roh.md ; '
         'configuration D16 et manifests D15, D14 et D11a liés par SHA-256.</p>'),
        ("Configuration, graines et versions", f'<p>Méthode {escape(manifest["method_id"])}. '
         f'Configuration {escape(manifest["configuration"]["configuration_id"])}. '
         f'Graine dérivée {manifest["seed"]["derived_seed"]}. SHA-256 de la graine de recombinaison, '
         'du rôle et de la méthode ; sélection avant observation des génotypes.</p>'),
        ("Vérité simulée", '<p>Un témoin sans injection, une ROH éloignée de la cible, une ROH centrée sur la cible '
         'sans IBD injecté et une ROH centrée chez un sujet positif D15. Longueurs physiques enregistrées : '
         f'{lengths}. Leurs longueurs génétiques sont interpolées séparément. '
         'Le témoin possède des bornes nulles.</p>'),
        ("Entrées aveugles", '<p>Le VCF observable contient seulement les sondes de puce et les génotypes non phasés. '
         'Les copies, matrices, scénarios et origines restent dans les vérités privées. L’export final relève de D18.</p>'),
        ("Méthode d’évaluation", '<p>Copie H1 vers H2 dans l’intervalle configuré, contrôle de l’homozygotie et de '
         'l’IBD H1 préexistant, conservation des transmissions. Le sujet couplé appartient à une paire positive D15 ; '
         'les autres sujets sont hors pedigree et exclus des porteurs D14, des paires et des matrices D15. La densité, les trous entre '
         'sondes et les marqueurs flanquants sont vérifiés avant injection.</p>'),
        ("Résultats et incertitude", cards + table + bars),
        ("Critères préspecifiés", '<p>Zéro erreur de copie ROH, d’IBD ou de Mendel. H1 conservée, négatifs inchangés. '
         'Au moins 50 sondes, longueur minimale de 1 500 kb, densité au plus 50 kb par sonde, trou interne au plus '
         '1 000 kb et flancs suffisants. Les performances et leurs intervalles sont NOT_EVALUABLE : PLINK n’a pas '
         'encore analysé ce panel.</p>'),
        ("Anomalies et écarts", '<p>Aucun déplacement d’intervalle ni assouplissement de seuil. La vérité IBD D15 '
         'est copiée à l’identique. Le cas couplé crée en plus une copie IBD H2–H1 dans la même paire, '
         'sur l’intersection ROH/IBD : cette vérité dérivée est enregistrée séparément et ne constitue '
         'pas une répétition indépendante. Les changements H2 sont aussi une surcharge locale d’origine '
         'à intégrer dans les futures évaluations d’ascendance.</p>'),
        ("Limites", '<p>Une origine commune est synthétiquement imposée ; aucune généalogie consanguine complète '
         'n’est reconstruite. L’homozygotie du fond peut prolonger un ROH au-delà des bornes injectées. '
         'Le témoin sans injection peut posséder une ROH historique. TARGET_DISTINCT signifie absence de la cible '
         'dans la ROH, sans garantir une disjonction avec toutes les coordonnées fondatrices d’autres familles. '
         'Un seul chromosome et une grille réduite de longueurs préparent les évaluateurs ; ils ne qualifient pas le fardeau '
         'genome-wide ni les performances par classe de longueur.</p>'),
        ("Décision et suite", '<p>PASS technique si les contrôles sources sont conformes. Statut scientifique : '
         'NOT_EVALUABLE. calibrated=false ; gamma_inputs_modified=false ; founder_effect_proven=false. '
         'Prochaine étape : D17, imperfections de mesure. Le label SIMULATION_CALIBRATED restera réservé à la '
         'qualification après validation aveugle et décision humaine.</p>'),
        ("Provenance et empreintes", provenance),
    ]


def build_report(*, run_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit HTML autonome/JSON/empreintes sous un nouveau dossier du workspace."""
    if output_root.absolute().parent.name != "calibration_workspace" or output_root.exists():
        raise ValueError("d16_report_output_invalid")
    if any(path.is_symlink() for path in (output_root, *output_root.parents)):
        raise ValueError("d16_report_symlink_forbidden")
    manifest_path = run_root / "private_truth/d16/roh_simulation_manifest.json"
    manifest = _json(manifest_path)
    _validate_schema(manifest, "roh_simulation_manifest")
    for item in manifest["artifacts"]:
        _artifact(run_root, manifest, item["role"])
    truth_path = _artifact(run_root, manifest, "ROH_INJECTION_TRUTH")
    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    _validate_schema(truth, "roh_injection_truth")
    sections = _sections(manifest, _sha(manifest_path), truth)
    title = "D16 — Validation des injections ROH"
    head = re.split(r"<body\b", template_path.read_text(encoding="utf-8"), maxsplit=1)[0]
    head = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", head, flags=re.S)
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span>'
        '<span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, (heading, _) in enumerate(sections, 1)
    ) + '</ol></nav>'
    # Le pied courant précède main pour être visible dès la première page de contenu.
    body = f'<body lang="fr"><section class="cover"><h1>{title}</h1><p>Effet_fondateur2 · '
    body += escape(manifest["created_at"][:10]) + ' · Contrôle technique régional</p></section>' + toc
    body += '<footer class="page-footer"><span>D16 · Effet_fondateur2</span></footer><main>'
    body += "".join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>'
                    for index, (heading, content) in enumerate(sections, 1))
    body += '</main></body></html>'
    output_root.mkdir(parents=True)
    output = output_root / "index.html"
    output.write_text(head + body, encoding="utf-8")
    (output_root / "report_data.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output_root / "checksums.sha256").write_text("".join(
        f"{_sha(path)}  {path.name}\n" for path in (output, output_root / "report_data.json")
    ), encoding="utf-8")
    return output


def main() -> int:
    """CLI de rendu du rapport technique D16, sans nouvelle analyse scientifique."""
    parser = argparse.ArgumentParser(description="Rapport de validation technique ROH D16")
    for argument in ("run-root", "template-path", "output-root"):
        parser.add_argument("--" + argument, type=Path, required=True)
    build_report(**vars(parser.parse_args()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
