"""Rapport technique autonome D15, alimenté uniquement par le run finalisé."""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re

import jsonschema

from effet_fondateur.calibration.materialize_ibd import _artifact, _json, _sha


def build_report(*, run_root: Path, template_path: Path, output_root: Path) -> Path:
    """Publie un rapport D15 français depuis les agrégats et empreintes sources.

    Le modèle HTML autonome est fourni explicitement (report-formater).
    Le rendu ne calcule aucune sensibilité ou performance de détection.
    Le répertoire de rapport doit être neuf et rester dans calibration_workspace.
    """
    if output_root.absolute().parent.name != "calibration_workspace" or output_root.exists():
        raise ValueError("d15_report_output_invalid")
    manifest_path = run_root / "private_truth/d15/ibd_simulation_manifest.json"
    manifest = _json(manifest_path)
    schema_path = Path(__file__).resolve().parents[3] / "calibration/schemas/ibd_simulation_manifest.schema.json"
    jsonschema.Draft202012Validator(_json(schema_path)).validate(manifest)
    for item in manifest["artifacts"]:
        _artifact(run_root, manifest, item["role"])
    truth_path = _artifact(run_root, manifest, "IBD_INJECTION_TRUTH")
    truth = json.loads(truth_path.read_text())
    integrity = manifest["integrity"]
    title = "D15 — Validation des injections IBD"
    sections = [
        ("Identification", f"<p>Rapport CAL-DEV-D15 · version 1.0 · {escape(manifest['created_at'])}. "
         f"Run : {escape(manifest['ibd_simulation_manifest_id'])}.</p>"),
        ("Résumé exécutif", "<ul><li>Injection technique des cinq longueurs prévues réalisée.</li>"
         "<li>Paires et matrices disjointes ; transmissions mendéliennes contrôlées.</li>"
         "<li>Le panel régional prépare la validation des appels IBD.</li></ul>"),
        ("Question évaluée", "<p>La copie d’un haplotype commun sur un intervalle préspecifié "
         "produit-elle une vérité traçable sans modifier les paires négatives ni rompre Mendel ?</p>"),
        ("Sources et autorité", "<p>AGENTS.md ; protocole D03 §8.4 ; critères D04 IBD ; "
         "plan D14b–D15 ; manifests D14b, reçu de rejeu, D11a et D14 liés dans le run.</p>"),
        ("Configuration, graines et versions", f"<p>Méthode : {escape(manifest['method_id'])}. "
         f"Rôle : {escape(manifest['panel_role'])}. Graine dérivée : {manifest['seed']['derived_seed']}. "
         "Dérivation SHA-256 explicitement enregistrée ; aucun tirage fondé sur un appel IBD.</p>"),
        ("Vérité simulée", "<p>Une paire positive et une paire négative par longueur : "
         "1 ; 1,5 ; 2 ; 3 ; 5 cM. H1 des deux sujets positifs reçoit H1 d’une même matrice "
         "sur l’intervalle. Les bornes négatives sont nulles. La vérité détaillée reste privée.</p>"),
        ("Entrées aveugles", "<p>Le VCF observable conserve les identifiants synthétiques, "
         "les sondes de puce et les génotypes non phasés. Aucun libellé de scénario ni matrice "
         "n’y est ajouté. L’export vers blinded_inputs/ relève de D18.</p>"),
        ("Méthode d’évaluation", "<p>Contrôles d’empreintes, de disjonction, d’identité injectée, "
         "d’invariance des paires négatives et des transmissions. Les bornes physiques sont "
         "arrondies au bp puis reconverties en cM ; la longueur réalisée peut différer légèrement "
         "de la longueur demandée.</p>"),
        ("Résultats et incertitude", ""),
        ("Critères préspecifiés", "<p>Critères techniques : zéro erreur de copie, de modification "
         "négative ou de transmission ; univers observable intégralement conservé ; seuils "
         "de marqueurs D14b respectés. Les intervalles de performance sont NOT_EVALUABLE : "
         "aucun outil de détection n’a été exécuté.</p>"),
        ("Anomalies et écarts", "<p>Aucun élargissement de fenêtre ni déplacement de cible. "
         "Le run conserve les fondateurs et les recombinaisons D14. Les nouveaux sujets positifs "
         "sont des témoins techniques IBD ; ils ne constituent pas de nouvelles familles mutantes.</p>"),
        ("Limites", "<p>Le fond 1000G peut contenir du partage ancestral non connu et de "
         "l’identité par état. Une paire négative signifie uniquement absence d’injection. "
         "Les limites injectées ne garantissent pas les limites maximales des appels. Les "
         "matrices et sujets distincts ne prouvent pas une indépendance populationnelle complète. "
         "Les 22 autosomes, 3 000 négatifs et répétitions de qualification restent à réaliser.</p>"),
        ("Décision et suite", "<p>PASS technique D15 si tous les contrôles sources sont nuls. "
         "Statut scientifique : NOT_EVALUABLE. calibrated=false et founder_effect_proven=false. "
         "Étape suivante : D16, injection des ROH. La qualification SIMULATION_CALIBRATED "
         "nécessitera les rapports de performance et la décision humaine prévus.</p>"),
        ("Provenance et empreintes", ""),
    ]
    cards = '<div class="stat-grid">' + "".join(
        f'<div class="stat-card"><div class="num">{value}</div><div class="label">{label}</div></div>'
        for value, label in ((manifest["positive_pair_count"], "Paires positives"),
                             (manifest["negative_pair_count"], "Paires négatives"),
                             (integrity["observed_variant_count"], "Marqueurs observables"))
    ) + '</div><figcaption>Figure 1 — Effectifs techniques. Source : manifest D15.</figcaption>'
    table = '<table><caption>Tableau 1 — Contrôles du run (manifest D15).</caption><tbody>' + "".join(
        f"<tr><th>{escape(key)}</th><td>{value}</td></tr>" for key, value in integrity.items()
    ) + '</tbody></table>'
    bars = '<svg viewBox="0 0 600 220" role="img" aria-label="Longueurs IBD demandées">'
    for index, item in enumerate(item for item in truth if item["is_positive"]):
        length = item["requested_length_cm"]
        y = 30 + index * 35
        bars += f'<text x="8" y="{y+14}">{length:g} cM</text><rect x="90" y="{y}" width="{length*80}" height="20" fill="#236c75"/>'
    bars += '</svg><figcaption>Figure 2 — Longueurs demandées. Source : vérité D15.</figcaption>'
    sections[8] = (sections[8][0], cards + table + bars)
    sections[13] = (sections[13][0], "<p>Empreinte du manifest source :</p><p style='overflow-wrap:anywhere'>" + _sha(manifest_path) + "</p><ul>" + "".join(
        f"<li>{escape(key)} : <code style='overflow-wrap:anywhere'>{escape(value)}</code></li>"
        for key, value in manifest["inputs"].items()
    ) + "</ul>")
    template = template_path.read_text(encoding="utf-8")
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    head = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", head, flags=re.S)
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(
        f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>'
        for index, (heading, _) in enumerate(sections, 1)
    ) + '</ol></nav>'
    # Paged.js doit rencontrer l'élément courant avant les pages de contenu.
    # Après main, il ne répète le pied qu'à partir de la dernière page.
    footer = '<footer class="page-footer"><span>D15 · Effet_fondateur2</span></footer>'
    body = '<body lang="fr"><section class="cover"><h1>' + title + '</h1><p>Effet_fondateur2 · 29 octobre 2026 · Contrôle technique régional</p></section>' + toc + footer + '<main>' + "".join(
        f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>'
        for index, (heading, content) in enumerate(sections, 1)
    ) + '</main></body></html>'
    if any(parent.is_symlink() for parent in (output_root, *output_root.parents)):
        raise ValueError("d15_report_symlink_forbidden")
    output_root.mkdir(parents=True)
    output = output_root / "index.html"
    output.write_text(head + body, encoding="utf-8")
    (output_root / "report_data.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output_root / "checksums.sha256").write_text("".join(
        f"{_sha(path)}  {path.name}\n" for path in (output, output_root / "report_data.json")
    ), encoding="utf-8")
    return output


def main() -> int:
    """Rend le rapport depuis un modèle fourni, sans accès aux données d'étude."""
    parser = argparse.ArgumentParser(description="Rapport de validation technique D15")
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--template-path", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    build_report(**vars(parser.parse_args()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
