"""Rapport autonome CAL-REP-16A de D27."""

from __future__ import annotations

import argparse
from html import escape
import json
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any

from effet_fondateur.audit import sha256_file
from effet_fondateur.calibration.qc_bench import write_json
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document


SCHEMAS = Path(__file__).resolve().parents[3] / "calibration/schemas"
HEADINGS = (
    "Identification et domaine de validité", "Résumé décisionnel",
    "Question évaluée", "Autorité et séparation de la vérité",
    "Panel, référence et indépendance", "Moteur réel évalué",
    "Estimateur du mélange simulé", "Résultat du mélange régional",
    "Définition de la stabilité", "Stabilité du positionnement diploïde",
    "Stabilité du positionnement local", "Contrôle des copies fondatrices locales",
    "Critères et décision", "Limites, provenance et empreintes",
)
PROFILE_LABELS = {
    "control": "Contrôle", "missing_01": "Manque 1 %", "missing_05": "Manque 5 %",
    "genotype_01": "Dosage 1 %", "phase_point_02": "Phase ponctuelle 2 %",
    "phase_switch_02": "Commutation 2 %", "density_50": "Densité 50 %",
    "combined": "Combiné",
}


def _pct(value: float) -> str:
    return f"{100 * value:.2f} %"


def _num(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def _table(headers: tuple[str, ...], rows: list[list[str]], caption: str, number: int) -> str:
    return (
        f'<table><caption>Tableau {number} : {escape(caption)}</caption><thead><tr>'
        + "".join(f"<th>{escape(value)}</th>" for value in headers)
        + "</tr></thead><tbody>"
        + "".join("<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>" for row in rows)
        + "</tbody></table>"
    )


def _figure(svg: str, caption: str, number: int) -> str:
    return f'<figure class="scientific-figure">{svg}<figcaption>Figure {number} : {escape(caption)}</figcaption></figure>'


def _mixture_figure(rows: list[dict[str, Any]]) -> str:
    points = []
    for row in rows:
        truth = float(row["true_stratum_a_proportion"])
        estimate = float(row["estimated_mean"])
        x = 75 + 450 * truth
        y = 245 - 190 * estimate
        points.append(f'<line x1="{x}" y1="245" x2="{x}" y2="{y:.1f}" stroke="#9eb7c4"/><circle cx="{x}" cy="{y:.1f}" r="7" fill="#c4573a"/><text x="{x}" y="{y-13:.1f}" text-anchor="middle" font-size="12">{estimate:.3f}</text>')
    return (
        '<svg viewBox="0 0 600 290" role="img" aria-label="Proportion estimée selon la proportion simulée">'
        '<text x="16" y="20" font-size="14" font-weight="600">Moyenne estimée par proportion vraie</text>'
        '<line x1="75" y1="245" x2="525" y2="55" stroke="#286d91" stroke-dasharray="6 4"/>'
        '<line x1="75" y1="245" x2="525" y2="245" stroke="#333"/><line x1="75" y1="245" x2="75" y2="45" stroke="#333"/>'
        + "".join(points)
        + '<text x="300" y="280" text-anchor="middle" font-size="12">Proportion vraie de la strate A</text>'
        '<text transform="translate(18 155) rotate(-90)" text-anchor="middle" font-size="12">Proportion moyenne estimée</text>'
        '<text x="75" y="262" text-anchor="middle" font-size="11">0</text><text x="300" y="262" text-anchor="middle" font-size="11">0,5</text><text x="525" y="262" text-anchor="middle" font-size="11">1</text>'
        '</svg>'
    )


def _stability_figure(rows: list[dict[str, Any]], scope: str) -> str:
    selected = [row for row in rows if row["scope"] == scope and row["profile_id"] != "control"]
    maximum = max(float(row["mean_standardized_distance"]) for row in selected) or 1.0
    parts = [f'<svg viewBox="0 0 640 330" role="img" aria-label="Stabilité {scope.lower()}"><text x="14" y="20" font-size="14" font-weight="600">Distance moyenne standardisée — {escape(scope)}</text>']
    for index, row in enumerate(selected):
        value = float(row["mean_standardized_distance"])
        x = 45 + index * 84
        height = 220 * value / maximum
        parts.append(f'<rect x="{x}" y="{260-height:.1f}" width="54" height="{height:.1f}" fill="#286d91"/><text x="{x+27}" y="{250-height:.1f}" text-anchor="middle" font-size="11">{value:.3f}</text><text transform="translate({x+27} 282) rotate(35)" text-anchor="start" font-size="10">{escape(PROFILE_LABELS[row["profile_id"]])}</text>')
    parts.append('<line x1="35" y1="260" x2="625" y2="260" stroke="#333"/></svg>')
    return "".join(parts)


def _founder_figure(document: dict[str, Any]) -> str:
    local = document["local_founder_check"]
    value = float(local["mean_founder_reference_axis_proportion"])
    x = 70 + 470 * value
    return (
        '<svg viewBox="0 0 610 170" role="img" aria-label="Position des copies fondatrices sur l’axe de référence">'
        '<text x="15" y="22" font-size="14" font-weight="600">Position locale moyenne sur l’axe des deux repères</text>'
        '<line x1="70" y1="90" x2="540" y2="90" stroke="#333" stroke-width="3"/>'
        '<circle cx="70" cy="90" r="8" fill="#7597a8"/><circle cx="540" cy="90" r="8" fill="#c4573a"/>'
        f'<circle cx="{x:.1f}" cy="90" r="11" fill="#6b8e23"/><text x="{x:.1f}" y="65" text-anchor="middle" font-size="13">moyenne {value:.3f}</text>'
        '<text x="70" y="120" text-anchor="middle" font-size="11">repère de fond</text><text x="540" y="120" text-anchor="middle" font-size="11">repère fondateur</text>'
        f'<text x="305" y="150" text-anchor="middle" font-size="11">{local["closer_to_founder_reference_count"]}/{local["carrier_copy_count"]} copies plus proches du repère fondateur ; 1 haplotype indépendant</text></svg>'
    )


def _criterion_table(document: dict[str, Any]) -> str:
    return _table(
        ("Critère", "Métrique", "Valeur", "Seuil", "Statut"),
        [[row["criterion_id"], row["metric"], _num(float(row["value"])), "Non scellé", row["status"]] for row in document["criteria"]],
        "Décision selon les critères préspecifiés.", 3,
    )


def _sections(document: dict[str, Any], source_sha: str) -> list[str]:
    mixture = document["mixture"]
    global_rows = [row for row in document["stability_by_profile"] if row["scope"] == "GLOBAL"]
    local_rows = [row for row in document["stability_by_profile"] if row["scope"] == "LOCAL"]
    global_by_id = {row["profile_id"]: row for row in global_rows}
    local_by_id = {row["profile_id"]: row for row in local_rows}
    mixture_table = [[
        _num(float(row["true_stratum_a_proportion"]), 1), str(row["sample_count"]),
        _num(float(row["estimated_mean"])), _num(float(row["mean_absolute_error"])),
    ] for row in mixture["results_by_true_proportion"]]
    stability_table = [[
        PROFILE_LABELS[profile], _num(float(global_by_id[profile]["mean_standardized_distance"])),
        _num(float(local_by_id[profile]["mean_standardized_distance"])),
        str(global_by_id[profile]["informative_variant_count"]),
    ] for profile in document["profiles"]]
    mixture_equation = '<div class="equation"><div class="math">α̂ᵢ = clip<sub>[0,1]</sub>[(zᵢ−c<sub>B</sub>)ᵀ(c<sub>A</sub>−c<sub>B</sub>) / ‖c<sub>A</sub>−c<sub>B</sub>‖²]</div><p>Les coordonnées z sont les scores PC divisés par √λ. MAE = n⁻¹Σ|α̂ᵢ−αᵢ|.</p></div>'
    stability_equation = '<div class="equation"><div class="math">dᵢ = √{K⁻¹ Σ<sub>k=1…K</sub>(z<sup>aligné</sup><sub>ik,p</sub>−z<sub>ik,0</sub>)²}</div><p>K = 10. L’alignement orthogonal utilise seulement les 1 403 références communes.</p></div>'
    return [
        '<p><strong>Rapport :</strong> CAL-REP-16A · D27. <strong>Domaine :</strong> 487 marqueurs régionaux du chromosome 5, GRCh38. <strong>Objet :</strong> fidélité du cœur PCA/projection de 16A.</p><aside class="callout">Le mot « global » décrit ici le calcul diploïde de 16A appliqué à un panel régional. Il ne signifie pas genome-wide.</aside>',
        f'<p>Le cœur de l’évaluateur passe ses contrôles techniques : <strong>{document["technical_status"]}</strong>. La MAE de proportion simulée vaut <strong>{_num(mixture["mean_absolute_error"])}</strong>, IC95 bootstrap {_num(mixture["bootstrap_ci95_low"])}–{_num(mixture["bootstrap_ci95_high"])}. L’instabilité primaire maximale vaut <strong>{_num(document["repeated_position_instability"])}</strong> distance standardisée.</p><aside class="callout">Qualification scientifique : <strong>{document["performance_qualification"]}</strong>. Aucun seuil n’était scellé ; ces nombres ne produisent donc pas de PASS scientifique.</aside>',
        '<p>D27 demande si des mélanges synthétiques connus occupent la position attendue entre deux repères et si les mêmes individus gardent une position stable lorsque densité, manque, dosage ou phase sont perturbés.</p><p>Il ne demande pas d’assigner une population, un territoire, une ethnicité ou un ancêtre à un sujet. Les strates A/B sont uniquement des codes de simulation.</p>',
        '<p>Les coordonnées sont calculées et scellées avant lecture des vérités D12/D14. L’évaluation ouvre ensuite les proportions simulées et les six copies fondatrices. Les axes et leur alignement ne voient jamais ces vérités.</p><p>Le manifeste d’appels atteste <code>SEALED_BEFORE_TRUTH</code>. Les sorties agrégées conservent les empreintes de la configuration, du protocole et des coordonnées.</p>',
        f'<p>La référence contient {document["reference_entity_count"]} individus de la partition D10 <code>{escape(document["reference_role"])}</code>. La préparation D23b atteste l’exclusion des deux partitions de donneurs, dont les 900 donneurs ayant produit D12.</p><p>Le contraste de mélange comprend {mixture["sample_count"]} sujets : deux groupes purs et 20 mélanges réciproques 50/50. Les huit profils D17 portent sur les mêmes 908 individus.</p>',
        '<p>Le banc appelle directement <code>harmonize_alt_dosages</code>, qui ajuste <code>fit_reference_pca</code> sur la référence seulement, puis projette l’étude avec <code>project_pca</code>. Dix composantes sont demandées ; l’imputation des valeurs manquantes utilise la fréquence de référence.</p><p>Le parcours d’orchestration, le cache, le téléchargement et SHAPEIT5 ne sont pas exécutés. Le verdict technique est donc limité au cœur numérique et à l’évaluateur.</p>',
        '<p>La proportion n’est pas une sortie native de 16A. D27 la dérive de façon préspecifiée par projection orthogonale sur le segment joignant deux centroïdes de référence, après standardisation par les valeurs propres. Les estimations hors segment sont tronquées à [0,1].</p>' + mixture_equation,
        f'<p>La MAE totale est {_num(mixture["mean_absolute_error"])}. Le mélange vrai 0,5 est centré à {_num(mixture["results_by_true_proportion"][1]["estimated_mean"])} en moyenne, mais son erreur absolue reste {_num(mixture["results_by_true_proportion"][1]["mean_absolute_error"])}. Les deux groupes purs sont ramenés vers le centre, ce qui explique une MAE de {_num(mixture["results_by_true_proportion"][0]["mean_absolute_error"])} et {_num(mixture["results_by_true_proportion"][2]["mean_absolute_error"])}.</p><p>L’IC95 bootstrap reflète l’échantillonnage des sujets de ce panel, pas le changement de région, de puce ou de référence.</p>' + _figure(_mixture_figure(mixture["results_by_true_proportion"]), "Moyenne estimée face à la proportion simulée ; la diagonale indique l’accord parfait.", 1) + _table(("Proportion vraie A", "n", "Moyenne estimée", "MAE"), mixture_table, "Erreur du mélange par proportion simulée.", 1),
        '<p>Les PCA réajustées sur des univers de marqueurs différents peuvent subir rotations et inversions. D27 aligne chaque profil au contrôle par Procrustes en n’utilisant que les références communes. La distance est ensuite exprimée en unités de variance de composante.</p><p>Pour les haplotypes locaux, les deux copies sont appariées selon la permutation de distance minimale. Cela retire l’orientation arbitraire H1/H2 sans effacer une commutation interne.</p>' + stability_equation,
        f'<p>Le manque 1 % déplace en moyenne les positions de {_num(global_by_id["missing_01"]["mean_standardized_distance"])} ; le manque 5 % de {_num(global_by_id["missing_05"]["mean_standardized_distance"])} et l’erreur de dosage 1 % de {_num(global_by_id["genotype_01"]["mean_standardized_distance"])}. Les erreurs de phase seules laissent le dosage inchangé et ont une distance numérique nulle.</p><p>La demi-densité atteint {_num(global_by_id["density_50"]["mean_standardized_distance"])} et le combiné {_num(global_by_id["combined"]["mean_standardized_distance"])}. Ils ne conservent que 244 variants informatifs contre 487.</p>' + _figure(_stability_figure(document["stability_by_profile"], "GLOBAL"), "Instabilité moyenne du positionnement diploïde par profil.", 2),
        f'<p>La phase ponctuelle produit une distance locale moyenne de {_num(local_by_id["phase_point_02"]["mean_standardized_distance"])}. Les commutations atteignent {_num(local_by_id["phase_switch_02"]["mean_standardized_distance"])} malgré l’appariement optimal des copies ; {local_by_id["phase_switch_02"]["optimal_haplotype_swap_count"]} sujets requièrent une permutation globale H1/H2.</p><p>La demi-densité vaut {_num(local_by_id["density_50"]["mean_standardized_distance"])} et le combiné {_num(local_by_id["combined"]["mean_standardized_distance"])}. Le combiné fixe ainsi la métrique primaire conservatrice.</p>' + _figure(_stability_figure(document["stability_by_profile"], "LOCAL"), "Instabilité moyenne du positionnement haplotypique local par profil.", 3) + _table(("Profil", "Diploïde", "Local", "Variants diploïdes"), stability_table, "Distances moyennes appariées et densité informative.", 2),
        f'<p>Les six copies D14 portant le segment fondateur sont toutes plus proches du repère fondateur que du repère de fond sur l’axe préspecifié. Leur proportion axiale moyenne vaut {_num(document["local_founder_check"]["mean_founder_reference_axis_proportion"])}.</p><p>Ce contrôle est descriptif : les six copies descendent d’un seul haplotype fondateur. L’effectif indépendant est donc 1, sans intervalle de confiance ni prétention de sensibilité.</p>' + _figure(_founder_figure(document), "Position moyenne des copies porteuses entre les deux repères locaux.", 4),
        '<p>CAL-ANC-001 et CAL-ANC-002 restent <strong>REVIEW_REQUIRED</strong> : leurs seuils <code>pilot_then_seal</code> étaient absents avant l’ouverture des résultats. La stabilité sous manque faible est encourageante, mais la demi-densité, les commutations et le profil combiné déplacent substantiellement les positions.</p><p>Aucune attribution de strate n’est émise, et aucune promotion de production n’est autorisée.</p>' + _criterion_table(document),
        f'<p>Ce banc est régional, conditionnel à la phase candidate et limité à une puce, une région et une référence. La proportion est un estimateur dérivé à deux repères, pas une fonction publiée par 16A. Les résultats ne généralisent ni à une PCA autosomale complète ni à des mélanges multivoies.</p><p>Source agrégée : <span class="hash">{source_sha}</span>. Configuration : <span class="hash">{document["configuration_sha256"]}</span>. Protocole : <span class="hash">{document["protocol_sha256"]}</span>. Appels scellés : <span class="hash">{document["calls_manifest_sha256"]}</span>.</p><p><strong>Verdicts maintenus :</strong> calibrated=false ; founder_effect_proven=false ; production_promotion=false.</p>',
    ]


def _render(document: dict[str, Any], source_sha: str, template: str) -> str:
    sections = _sections(document, source_sha)
    head = re.split(r"<body\b", template, maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d27_report_template_invalid")
    style = '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.7}table{font-size:.68em}main p{orphans:3;widows:3}code{font-size:.9em}</style>'
    head = re.sub(r"<title>.*?</title>", "<title>CAL-REP-16A — Évaluation de l’ascendance D27</title>", head, flags=re.S).replace("</head>", style + "</head>")
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + "".join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + "</ol></nav>"
    body = '<body lang="fr"><section class="cover"><h1>CAL-REP-16A — Positionnement de référence et stabilité</h1><p>Effet_fondateur2 · version 1.0.0 · 10 novembre 2026</p><p>D27 · évaluation fermée du cœur de l’étape 16A</p></section>' + toc + "<main>"
    body += "".join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, sections, strict=True), 1))
    return head + body + '</main><footer class="page-footer"><span>CAL-REP-16A · Effet_fondateur2</span></footer></body></html>'


def build_report(*, evaluation_root: Path, template_path: Path, output_root: Path) -> Path:
    """Produit CAL-REP-16A depuis le seul JSON agrégé et empreinté."""
    manifest_path = evaluation_root / "evaluation_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source = evaluation_root / "CAL-REP-16A.json"
    source_sha = sha256_file(source)
    if manifest["artifacts"][0]["sha256"] != source_sha:
        raise ValueError("d27_report_source_checksum_invalid")
    document = json.loads(source.read_text(encoding="utf-8"))
    validate_json_document(document, "ancestry_evaluation.schema.json", schemas_dir=SCHEMAS)
    template = template_path.read_text(encoding="utf-8")
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_root.name}.", dir=destination.parent))
    try:
        report_dir = temporary / "CAL-REP-16A"
        report_dir.mkdir()
        html_path = report_dir / "index.html"
        html_path.write_text(_render(document, source_sha, template), encoding="utf-8")
        data_path = report_dir / "report_data.json"
        write_json(data_path, {"results": document})
        (report_dir / "checksums.sha256").write_text("".join(f"{sha256_file(path)}  {path.name}\n" for path in (html_path, data_path)), encoding="utf-8")
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary)
        raise
    if sha256_file(source) != source_sha:
        raise ValueError("d27_report_source_changed")
    return destination / "CAL-REP-16A/index.html"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-root", type=Path, required=True)
    parser.add_argument("--template-path", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    print(build_report(**vars(parser.parse_args())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
