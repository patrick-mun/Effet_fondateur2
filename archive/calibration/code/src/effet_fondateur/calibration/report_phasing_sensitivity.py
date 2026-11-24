"""Rapport D23b autonome à partir des seuls agrégats clos."""
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
from effet_fondateur.calibration.qc_bench import file_inventory, write_json
from effet_fondateur.calibration.report_technical_evaluation import HEADINGS
from effet_fondateur.calibration.technical_scenarios import new_workspace_run
from effet_fondateur.contracts import validate_json_document


DENSITY_ORDER = ("DENSITY_25", "DENSITY_50", "DENSITY_75", "DENSITY_100")


def _checked_evaluation(root: Path) -> dict[str, Any]:
    manifest = json.loads((root / "run_index.json").read_text())
    validate_json_document(
        manifest, "phasing_sensitivity_evaluation_manifest.schema.json",
        schemas_dir=Path(__file__).resolve().parents[3] / "calibration/schemas",
    )
    if manifest.get("method_id") != "paired_density_pedigree_phasing_evaluation_v1" or manifest.get("closure_state") != "CLOSED":
        raise ValueError("d23b_report_evaluation_not_closed")
    if manifest.get("calibrated") is not False or manifest.get("founder_effect_proven") is not False:
        raise ValueError("d23b_report_verdict_invalid")
    if file_inventory(root, ("run_index.json",)) != manifest.get("artifacts"):
        raise ValueError("d23b_report_inventory_invalid")
    return manifest


def _table(number: int, headers: tuple[str, ...], rows: list[list[str]], caption: str) -> str:
    head = "".join(f"<th>{escape(value)}</th>" for value in headers)
    body = "".join("<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>" for row in rows)
    return f'<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table><p class="caption">Tableau {number} — {escape(caption)}</p>'


def _figure(number: int, svg: str, caption: str) -> str:
    return f'<figure class="scientific-figure">{svg}<figcaption>Figure {number} — {escape(caption)}</figcaption></figure>'


def _series_plot(series: list[tuple[str, str, list[dict[str, Any]]]], title: str, aria: str) -> str:
    """Courbes exactes avec IC déjà calculés ; aucune estimation dans le rapport."""
    available = [row for _, _, rows in series for row in rows if row["estimate"] is not None]
    maximum = max((max([row["estimate"], *(row["interval_95"] or [])]) for row in available), default=1.0)
    maximum = max(maximum * 1.15, 0.01)
    content = f'<text x="12" y="20" font-weight="600">{escape(title)}</text>'
    for fraction in (0, .25, .5, .75, 1):
        y = 245 - 180 * fraction
        content += f'<line x1="70" x2="650" y1="{y:.1f}" y2="{y:.1f}" stroke="#dce3e8"/><text x="62" y="{y+4:.1f}" text-anchor="end" font-size="10">{100*maximum*fraction:.1f}%</text>'
    colors = ("#286d91", "#b0523a", "#4d7c56")
    for series_index, (label, marker, rows) in enumerate(series):
        points = []
        color = colors[series_index]
        by_profile = {row["profile"]: row for row in rows}
        for index, profile in enumerate(DENSITY_ORDER):
            row = by_profile.get(profile)
            if row is None or row["estimate"] is None:
                continue
            x, y = 110 + index * 165, 245 - 180 * row["estimate"] / maximum
            points.append(f"{x},{y:.2f}")
            if row["interval_95"]:
                low = 245 - 180 * row["interval_95"][0] / maximum
                high = 245 - 180 * row["interval_95"][1] / maximum
                content += f'<line x1="{x}" x2="{x}" y1="{high:.2f}" y2="{low:.2f}" stroke="{color}"/><line x1="{x-5}" x2="{x+5}" y1="{high:.2f}" y2="{high:.2f}" stroke="{color}"/><line x1="{x-5}" x2="{x+5}" y1="{low:.2f}" y2="{low:.2f}" stroke="{color}"/>'
            content += f'<{marker} cx="{x}" cy="{y:.2f}" r="5" fill="{color}"/><text x="{x}" y="{y-9:.2f}" text-anchor="middle" font-size="9">{100*row["estimate"]:.2f}%</text>'
        if len(points) > 1:
            content += f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="2"/>'
        content += f'<rect x="{90+series_index*190}" y="285" width="12" height="12" fill="{color}"/><text x="{108+series_index*190}" y="295" font-size="10">{escape(label)}</text>'
    for index, profile in enumerate(DENSITY_ORDER):
        content += f'<text x="{110+index*165}" y="268" text-anchor="middle" font-size="10">{profile.removeprefix("DENSITY_")} %</text>'
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 315" role="img" aria-label="{escape(aria, quote=True)}">{content}</svg>'


def _count_plot(performance: list[dict[str, Any]]) -> str:
    rows = [row for row in performance if row["scope"] == "AUTOSOMAL" and "observed_loci" in row]
    if not rows:
        return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 120" role="img" aria-label="Effectif non évaluable"><text x="12" y="20">Effectifs non évaluables : sorties de phase absentes.</text></svg>'
    maximum = max(row.get("observed_loci", 0) for row in rows)
    by_profile = {row["profile"]: row for row in rows}
    content = '<text x="12" y="20" font-weight="600">Loci évalués dans le parcours autosomal</text>'
    for index, profile in enumerate(DENSITY_ORDER):
        row = by_profile.get(profile)
        if row is None:
            content += f'<text x="12" y="{63+index*50}" font-size="11">{profile.removeprefix("DENSITY_")} % — non évaluable</text>'
            continue
        width = 500 * row["observed_loci"] / maximum
        y = 48 + index * 50
        content += f'<text x="12" y="{y+15}" font-size="11">{profile.removeprefix("DENSITY_")} %</text><rect x="78" y="{y}" width="{width:.2f}" height="22" fill="#286d91"/><text x="{84+width:.2f}" y="{y+15}" font-size="10">{row["observed_loci"]}</text>'
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 260" role="img" aria-label="Effectif de loci selon la densité">{content}</svg>'


def _contrast_plot(contrasts: list[dict[str, Any]]) -> str:
    rows = [row for row in contrasts if row["contrast_type"] == "WITHOUT_MINUS_WITH_PEDIGREE" and row["metric"] == "switch_error_rate"]
    maximum = max((max(abs(value) for value in [row["estimate"], *(row["interval_95"] or [])]) for row in rows if row["estimate"] is not None), default=.01)
    maximum *= 1.2
    center, scale = 340, 250 / maximum
    content = '<text x="12" y="20" font-weight="600">Retrait du pedigree : différence de SER</text><line x1="340" x2="340" y1="35" y2="245" stroke="#333"/>'
    for index, profile in enumerate(DENSITY_ORDER):
        row = next(row for row in rows if row["profile"] == profile)
        y = 55 + index * 48
        content += f'<text x="12" y="{y+4}" font-size="11">{profile.removeprefix("DENSITY_")} %</text>'
        if row["estimate"] is not None:
            x = center + row["estimate"] * scale
            low, high = row["interval_95"] or (row["estimate"], row["estimate"])
            content += f'<line x1="{center+low*scale:.2f}" x2="{center+high*scale:.2f}" y1="{y}" y2="{y}" stroke="#b0523a" stroke-width="3"/><circle cx="{x:.2f}" cy="{y}" r="5" fill="#b0523a"/><text x="{x+8:.2f}" y="{y-8}" font-size="9">{100*row["estimate"]:+.2f} pt</text>'
    content += '<text x="340" y="270" text-anchor="middle" font-size="10">0 = aucune différence ; positif = davantage d’erreurs sans pedigree</text>'
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 290" role="img" aria-label="Contrastes appariés du retrait du pedigree">{content}</svg>'


def _metric_rows(performance: list[dict[str, Any]], scope: str, mode: str, metric: str) -> list[dict[str, Any]]:
    return [dict(profile=row["profile"], **row["metrics"][metric]["ALL"])
            for row in performance if row["scope"] == scope and row["pedigree_mode"] == mode and metric in row.get("metrics", {})]


def _integer(value: int) -> str:
    return f"{value:,}".replace(",", " ")


def _scope_label(value: str) -> str:
    return {"AUTOSOMAL": "Autosomal", "REGIONAL": "Régional"}.get(value, value)


def _pedigree_label(value: str) -> str:
    return {
        "WITH_PEDIGREE": "Avec",
        "WITHOUT_PEDIGREE": "Sans",
        "WITHOUT_PEDIGREE_BY_CONTRACT": "Sans — contrat 12A",
        "PAIRED": "Apparié",
    }.get(value, value)


def _sections(manifest: dict[str, Any], technical: list[dict[str, Any]], performance: list[dict[str, Any]],
              contrasts: list[dict[str, Any]], preparation: dict[str, Any], source_sha: str) -> list[str]:
    regional_with = _metric_rows(performance, "REGIONAL", "WITH_PEDIGREE", "switch_error_rate")
    regional_without = _metric_rows(performance, "REGIONAL", "WITHOUT_PEDIGREE", "switch_error_rate")
    autosomal = _metric_rows(performance, "AUTOSOMAL", "WITHOUT_PEDIGREE_BY_CONTRACT", "switch_error_rate")
    primary_regional = next((row for row in regional_with if row["profile"] == "DENSITY_100"), None)
    primary_autosomal = next((row for row in autosomal if row["profile"] == "DENSITY_100"), None)
    primary_summary = (
        f'Au maximum, la SER régionale avec pedigree vaut {100*primary_regional["estimate"]:.2f} % ; '
        f'la SER autosomale sans pedigree vaut {100*primary_autosomal["estimate"]:.2f} %.'
        if primary_regional is not None and primary_regional["estimate"] is not None
        and primary_autosomal is not None and primary_autosomal["estimate"] is not None
        else "Au moins un parcours primaire est non évaluable ; aucun taux de remplacement n’est inventé."
    )
    primary_regional_cell = next((row for row in performance if row["scope"] == "REGIONAL" and row["profile"] == "DENSITY_100" and row["pedigree_mode"] == "WITH_PEDIGREE"), None)
    primary_autosomal_cell = next((row for row in performance if row["scope"] == "AUTOSOMAL" and row["profile"] == "DENSITY_100" and row["pedigree_mode"] == "WITHOUT_PEDIGREE_BY_CONTRACT"), None)
    hamming_warning = ""
    if primary_regional_cell is not None and primary_autosomal_cell is not None:
        regional_hamming = primary_regional_cell.get("metrics", {}).get("haplotype_hamming_rate", {}).get("ALL", {}).get("estimate")
        autosomal_hamming = primary_autosomal_cell.get("metrics", {}).get("haplotype_hamming_rate", {}).get("ALL", {}).get("estimate")
        if regional_hamming is not None and autosomal_hamming is not None:
            hamming_warning = f' Le Hamming de phase reste toutefois à {100*regional_hamming:.2f} % en régional et {100*autosomal_hamming:.2f} % en autosomal.'
    figures = (
        _figure(1, _count_plot(performance), "Effectifs de loci phasés par densité. Source : performance_evaluations.json ; la cible technique est exclue de 12A."),
        _figure(2, _series_plot([("Avec pedigree", "circle", regional_with), ("Sans pedigree", "circle", regional_without)], "Erreur de commutation régionale", "SER régionale selon densité et pedigree"), "Taux d’erreur de commutation régional et IC descriptifs à 95 %. Source : oracle D23b fermé."),
        _figure(3, _series_plot([("12A sans pedigree", "circle", autosomal)], "Erreur de commutation autosomale", "SER autosomale selon densité"), "Taux d’erreur de commutation autosomal et IC descriptifs à 95 %. Source : oracle D23b fermé."),
        _figure(4, _contrast_plot(contrasts), "Différence sans pedigree moins avec pedigree, avec IC bootstrap appariés. Source : paired_contrasts.json."),
    )
    performance_rows = []
    for row in sorted(performance, key=lambda value: (value["scope"], value["pedigree_mode"], DENSITY_ORDER.index(value["profile"]))):
        measure = row.get("metrics", {}).get("switch_error_rate", {}).get("ALL")
        if measure is None:
            performance_rows.append([_scope_label(row["scope"]), row["profile"].removeprefix("DENSITY_") + " %", _pedigree_label(row["pedigree_mode"]), row["status"], str(row.get("observed_loci", "NE")), "NE", "NE", "NE"])
            continue
        interval = "NE" if measure["interval_95"] is None else f'{100*measure["interval_95"][0]:.2f}–{100*measure["interval_95"][1]:.2f} %'
        performance_rows.append([_scope_label(row["scope"]), row["profile"].removeprefix("DENSITY_") + " %", _pedigree_label(row["pedigree_mode"]), row["status"], _integer(row["observed_loci"]), f'{_integer(measure["numerator"])}/{_integer(measure["denominator"])}', f'{100*measure["estimate"]:.3f} %', interval])
    metric_labels = {
        "switch_error_rate": "Erreur de commutation",
        "haplotype_hamming_rate": "Hamming de phase",
        "isolated_point_motif_rate": "Motifs isolés",
        "rare_variant_phase_error_rate": "Erreur de phase rare",
        "target_mutation_orientation_accuracy": "Orientation cible correcte",
    }
    primary_rows = []
    primary_cells = [row for row in (primary_regional_cell, primary_autosomal_cell) if row is not None]
    for row in primary_cells:
        for metric, label in metric_labels.items():
            measure = row.get("metrics", {}).get(metric, {}).get("ALL")
            if measure is None:
                continue
            interval = "NE" if measure["interval_95"] is None else f'{100*measure["interval_95"][0]:.2f}–{100*measure["interval_95"][1]:.2f} %'
            estimate = "NE" if measure["estimate"] is None else f'{100*measure["estimate"]:.3f} %'
            primary_rows.append([_scope_label(row["scope"]), label, f'{_integer(measure["numerator"])}/{_integer(measure["denominator"])}', estimate, interval])
    contrast_rows = []
    for row in contrasts:
        if row["metric"] != "switch_error_rate":
            continue
        interval = "NE" if row["interval_95"] is None else f'{100*row["interval_95"][0]:+.2f} à {100*row["interval_95"][1]:+.2f} pt'
        contrast_label = "Densité − 100 %" if row["contrast_type"] == "DENSITY_MINUS_DENSITY_100" else "Sans − avec"
        contrast_rows.append([contrast_label, _scope_label(row["scope"]), row["profile"].removeprefix("DENSITY_") + " %", _pedigree_label(row["pedigree_mode"]), "NE" if row["estimate"] is None else f'{100*row["estimate"]:+.3f} pt', interval])
    technical_rows = [[_scope_label(row["scope"]), row["profile"].removeprefix("DENSITY_") + " %", _pedigree_label(row["pedigree_mode"]), row["stage_id"], row["status"]] for row in technical]
    density_rows = []
    for profile in DENSITY_ORDER:
        autosomal_technical = next(row for row in technical if row["scope"] == "AUTOSOMAL" and row["profile"] == profile)
        autosomal_performance = next(row for row in performance if row["scope"] == "AUTOSOMAL" and row["profile"] == profile)
        regional_technical = next(row for row in technical if row["scope"] == "REGIONAL" and row["profile"] == profile)
        density_rows.append([
            profile.removeprefix("DENSITY_") + " %",
            _integer(autosomal_technical["source_markers"]),
            _integer(autosomal_technical["excluded_monomorphic"]),
            _integer(autosomal_performance["observed_loci"]),
            _integer(regional_technical["observed_markers"]),
        ])
    pedigree_contrasts = {
        row["profile"]: row for row in contrasts
        if row["contrast_type"] == "WITHOUT_MINUS_WITH_PEDIGREE"
        and row["metric"] == "switch_error_rate"
    }
    maximum_pedigree = pedigree_contrasts["DENSITY_100"]
    quarter_pedigree = pedigree_contrasts["DENSITY_25"]
    primary_regional_hamming = primary_regional_cell["metrics"]["haplotype_hamming_rate"]["ALL"]
    primary_autosomal_hamming = primary_autosomal_cell["metrics"]["haplotype_hamming_rate"]["ALL"]
    primary_regional_rare = primary_regional_cell["metrics"]["rare_variant_phase_error_rate"]["ALL"]
    primary_autosomal_rare = primary_autosomal_cell["metrics"]["rare_variant_phase_error_rate"]["ALL"]
    primary_regional_isolated = primary_regional_cell["metrics"]["isolated_point_motif_rate"]["ALL"]
    primary_autosomal_isolated = primary_autosomal_cell["metrics"]["isolated_point_motif_rate"]["ALL"]
    primary_target = primary_regional_cell["metrics"]["target_mutation_orientation_accuracy"]["ALL"]
    provenance_rows = [
        ["Préparation maximale", manifest["preparation_manifest_sha256"], "Panel et vérité simulée"],
        ["Exécution D23b", manifest["source_run_sha256"], "24 composants 10A/11/12/12A"],
        ["Évaluation D23b", source_sha, "Agrégats, IC et contrastes"],
        ["Configuration D23b", manifest["configuration_sha256"], "Densités, graines et paramètres"],
        ["Préspecification", manifest["prespecification_sha256"], "Questions et règles fixées avant résultats"],
        ["Paquet CytoScan", preparation["chip_package_sha256"], "Annotation officielle NA36.r1"],
        ["Archive des cartes", preparation["map_archive_sha256"], "Cartes GRCh38 épinglées"],
    ]
    equations = '<div class="equation"><p><strong>Équation 1 — Taux agrégé</strong></p><div class="math" role="math">R = (Σ<sub>u</sub> e<sub>u</sub>) / (Σ<sub>u</sub> n<sub>u</sub>)</div><p>u indexe 24 unités indépendantes de simulation ; e compte les erreurs, n les opportunités.</p></div>'
    equations += '<div class="equation"><p><strong>Équation 2 — Contrastes préspecifiés</strong></p><div class="math" role="math">Δ<sub>densité</sub>(d,p) = R(d,p) − R(100 %,p)<br>Δ<sub>pedigree</sub>(d) = R(d,sans) − R(d,avec)</div><p>Un contraste positif indique davantage d’erreurs dans le premier terme. Il s’agit de différences en points de proportion.</p></div>'
    equations += '<div class="equation"><p><strong>Équation 3 — Intervalle apparié</strong></p><div class="math" role="math">Δ* = (Σ e<sub>u*</sub>/Σ n<sub>u*</sub>)<sub>1</sub> − (Σ e<sub>u*</sub>/Σ n<sub>u*</sub>)<sub>2</sub><br>IC<sub>95 %</sub> = [q<sub>0,025</sub>(Δ*), q<sub>0,975</sub>(Δ*)]</div><p>Mêmes 24 unités tirées avec remise dans les deux termes ; 2 000 tirages, graine 2311. Un tirage à dénominateur nul est exclu et compté.</p></div>'
    return [
        '<p><strong>Document :</strong> CAL-DEV-D23B — Sensibilité contrôlée de la qualité du phasage à la densité de marqueurs et au pedigree.</p>'
        '<p><strong>Projet :</strong> Effet_fondateur2. <strong>Version du rapport :</strong> 1.1.0. '
        f'<strong>Date des agrégats :</strong> {escape(manifest["created_at"][:10])}. <strong>Statut :</strong> évaluation close, performance REVIEW_REQUIRED.</p>'
        '<h3>Périmètre exact</h3><p>Ce document rend compte des parties 1 à 3 de l’investigation complémentaire D23b : préspecification, détermination de la densité maximale exploitable, puis comparaison appariée des densités et du pedigree. Il ne constitue ni une qualification du pipeline complet, ni une validation d’IBD, ni une preuve d’effet fondateur.</p>'
        '<p>Les calculs portent exclusivement sur un panel synthétique de calibration. La cible du chromosome 5 est un repère technique APC et non la mutation DOCK6 étudiée par le projet. Les données réelles <code>dock6_real</code>, les résultats historiques de production et le cache interne n’ont pas été utilisés.</p>'
        '<aside class="callout"><strong>Lecture du statut.</strong> « CLOSED » signifie que les artefacts attendus sont présents et liés par empreinte. « REVIEW_REQUIRED » signifie que les performances scientifiques ne satisfont à aucun seuil qualifiant scellé. Ces deux notions ne sont pas interchangeables.</aside>',

        '<div class="stat-grid"><div class="stat-card"><div class="num">243 844</div><div class="label">SNP autosomiques harmonisés</div></div>'
        '<div class="stat-card"><div class="num">4,31 %</div><div class="label">SER régionale maximale, avec pedigree</div></div>'
        '<div class="stat-card"><div class="num">4,47 %</div><div class="label">SER autosomale maximale</div></div>'
        '<div class="stat-card"><div class="num">24</div><div class="label">Unités indépendantes rééchantillonnées</div></div></div>'
        f'<p>Le résultat central est une diminution régulière de l’erreur de commutation lorsque la densité augmente. {primary_summary}{hamming_warning} La densité améliore donc nettement la cohérence locale de phase dans ce banc, sans suffire à établir une phase globalement fiable.</p>'
        f'<p>Le retrait du pedigree n’a pas le même effet à toutes les densités. À 25 %, la différence sans moins avec pedigree est de {100*quarter_pedigree["estimate"]:+.2f} points, avec un IC à 95 % de {100*quarter_pedigree["interval_95"][0]:+.2f} à {100*quarter_pedigree["interval_95"][1]:+.2f}. À 100 %, elle n’est plus que de {100*maximum_pedigree["estimate"]:+.2f} point, IC {100*maximum_pedigree["interval_95"][0]:+.2f} à {100*maximum_pedigree["interval_95"][1]:+.2f} : l’effet est alors non résolu par ce panel.</p>'
        '<ul><li><strong>Conclusion autorisée :</strong> dans ces données simulées et avec ces paramètres, davantage de marqueurs est associé à moins de commutations.</li><li><strong>Conclusion non autorisée :</strong> la densité n’est pas démontrée comme cause unique de toutes les erreurs de phase.</li><li><strong>Décision :</strong> conserver <code>calibrated=false</code> et <code>founder_effect_proven=false</code>.</li><li><strong>Suite :</strong> ne pas ouvrir D24 avant un audit séparé de la vérité de segments fondateurs.</li></ul>'
        '<aside class="callout"><strong>Message principal.</strong> Les contrôles techniques passent, mais le Hamming et l’erreur sur variants rares restent trop élevés pour transformer ce succès d’exécution en validation scientifique.</aside>',

        '<h3>Question 1 — densité</h3><p>À sujets, référence, cartes, paramètres SHAPEIT5 et vérité identiques, comment évoluent les erreurs de phase lorsque l’on conserve 25 %, 50 %, 75 % puis 100 % des marqueurs harmonisés ? Le contraste préspecifié compare chaque densité au maximum, sur les mêmes 24 unités.</p>'
        '<h3>Question 2 — pedigree</h3><p>Dans le parcours régional 11/12, quel changement observe-t-on lorsque les liens parentaux sont retirés, tout en conservant exactement les mêmes individus, génotypes, familles nominales, marqueurs et paramètres ? Cette comparaison est réalisée séparément à chaque densité.</p>'
        '<h3>Ce que l’expérience ne cherche pas à mesurer</h3><p>Le banc ne compare pas deux logiciels, ne teste pas un panel clinique, ne mesure pas l’effet d’un manque de génotypes et ne tente pas d’identifier une cause biologique. Il ne compare pas causalement le parcours régional au parcours autosomal : leurs univers de loci et leurs adaptateurs diffèrent.</p>'
        '<p>Aucune hypothèse de succès n’a été transformée en seuil après observation. Les taux et leurs intervalles sont descriptifs. L’absence de zéro dans un intervalle de contraste indique un signal résolu dans ce banc, pas une généralisation à d’autres puces, populations ou logiciels.</p>',

        '<h3>Autorité méthodologique</h3><p>La préspecification <code>D23B_DENSITY_PEDIGREE_PRESPECIFICATION.md</code> et la décision CAL-DEC-009 ont été fixées avant la matérialisation maximale et avant le phasage. Elles définissent les deux questions, les profils emboîtés, les unités, les métriques, les contrastes et la règle de décision.</p>'
        '<h3>Sources techniques</h3><ul><li>préparation maximale D23b close, contenant les génotypes synthétiques et la vérité privée ;</li><li>partition D10, dont 1 403 références publiques sont utilisées pour le phasage et restent distinctes des donneurs de simulation ;</li><li>annotation officielle CytoScan 750K Accel NA36.r1 pour l’univers de puce ;</li><li>VCF publics 1000 Genomes et cartes GRCh38 épinglées pour l’harmonisation et le phasage ;</li><li>sorties réelles des composants 10A, 11, 12 et 12A ;</li><li>oracle D23 et évaluateur D23b pour les métriques et contrastes.</li></ul>'
        '<h3>Règle de priorité</h3><p>La préspecification gouverne l’interprétation ; les manifests clos gouvernent les nombres et empreintes ; les sorties d’outils gouvernent les contrôles techniques. Le rapport ne remplace aucune de ces sources et ne recalcule pas les statistiques.</p>'
        '<p>Les sources explicitement exclues sont les données réelles DOCK6, les génotypes de validation aveugle, le cache interne et tout statut clinique. La mutation technique a été injectée par une règle moléculaire documentée, jamais inférée depuis un groupe ou une clinique.</p>',

        '<h3>Construction du maximum exploitable</h3><p>Le maximum ne correspond pas au nombre commercial inscrit sur le nom de la puce. Il est obtenu après harmonisation : position GRCh38 unique, carte disponible sans extrapolation, variant SNP biallélique présent dans la référence, orientation REF/ALT compatible et génotypes complets chez les 64 groupes donneurs mobilisés.</p>'
        f'<p>Ce filtrage retient {_integer(preparation["counts"]["markers"] - 1)} SNP autosomiques de puce. Une cible technique supplémentaire est ajoutée sur le chromosome 5, soit {_integer(preparation["counts"]["markers"])} marqueurs dans la préparation. Le chromosome le moins couvert en conserve {_integer(min(row["retained_chip_marker_count"] for row in preparation["chromosomes"]))} et le plus couvert {_integer(max(row["retained_chip_marker_count"] for row in preparation["chromosomes"]))}.</p>'
        '<h3>Profils emboîtés</h3><p>Les loci sont ordonnés par position dans chaque chromosome. Les rangs modulo quatre définissent des sous-ensembles strictement emboîtés : 25 % retient le rang 0 ; 50 % les rangs 0 et 2 ; 75 % les rangs 0, 1 et 2 ; 100 % tous les rangs. Aucun locus n’est sélectionné selon la qualité de sa phase observée.</p>'
        + _table(1, ("Densité", "Marqueurs source", "Monomorphes exclus", "Loci autosomaux évalués", "Loci régionaux"), density_rows, "Passage du profil source aux loci effectivement évalués. Source : technical_evaluations.json et performance_evaluations.json.")
        + '<p>Le tableau 1 explique la différence entre densité nominale, marqueurs transmis à 10A et loci réellement évaluables. Les variants monomorphes sont retirés du parcours autosomal ; la cible technique n’entre pas dans 12A. Ces exclusions sont déterministes et ne dépendent pas des erreurs de phase.</p>'
        '<h3>Paramètres constants</h3><p>Les 40 sujets, les 1 403 références, les cartes, la graine SHAPEIT5 2303, le thread unique, les paramètres de fenêtre et la vérité sont identiques dans toutes les cellules. Les intervalles utilisent 2 000 tirages appariés avec la graine 2311.</p>',

        '<h3>Origine de la vérité</h3><p>La vérité de phase est simulée avant et indépendamment des sorties SHAPEIT5. Elle comprend 32 fondateurs, huit trios et seize sujets isolés. Les fondateurs proviennent de 64 groupes donneurs D10 distincts ; les transmissions parent–enfant sont enregistrées explicitement.</p>'
        '<p>Les huit familles sont rééchantillonnées comme huit blocs entiers et les seize isolés comme seize blocs individuels. L’effectif statistique est donc de 24 unités indépendantes, et non de 40 sujets supposés indépendants, ni du nombre de loci.</p>'
        '<h3>Séparation entre vérité et données de phasage</h3><p>Les outils ne reçoivent que les génotypes non phasés, les cartes, la référence et, selon la cellule régionale, les liens parentaux. Les haplotypes vrais, copies fondatrices et transmissions restent privés jusqu’à la clôture des outils. L’évaluateur ouvre ensuite la vérité pour comparer les sorties.</p>'
        f'<h3>Audit exhaustif avant phasage</h3><p>{_integer(preparation["audit"]["compared_genotype_cells"])} cellules génotypiques ont été relues. Les discordances de génotype, ACPA, copie fondatrice, transmission, cible et pedigree sont toutes nulles. Cette concordance garantit que le banc mesure le phasage des données prévues ; elle ne préjuge pas de sa précision.</p>'
        '<h3>Alignement des haplotypes</h3><p>Une permutation globale H1/H2 n’est pas une erreur biologique. L’oracle choisit donc l’orientation globale minimisant le Hamming dans chaque chromosome ou bloc de phase. Les pertes de loci, discordances de dosage, ruptures de bloc et erreurs locales restent pénalisées ou signalées.</p>',

        '<h3>Contenu commun aux cellules</h3><p>Chaque cellule reçoit les mêmes 40 sujets, le même ordre, les mêmes FID, la même sélection de référence et le même sous-ensemble de marqueurs pour une densité donnée. Aucun génotype n’est ajouté, supprimé ou imputé pour fabriquer la comparaison pedigree.</p>'
        '<h3>Manipulation du pedigree</h3><p>En mode <code>WITH_PEDIGREE</code>, les PID/MID synthétiques décrivent huit trios. En mode <code>WITHOUT_PEDIGREE</code>, ces champs sont remplacés par zéro dans le FAM et dans <code>samples_master</code>. Les génotypes et identifiants de sujets restent inchangés. L’audit retrouve huit enregistrements parentaux dans le premier mode et zéro dans le second.</p>'
        '<p>Le parcours autosomal 10A/12A n’accepte pas de pedigree dans son adaptateur de production. Il est donc étiqueté <code>WITHOUT_PEDIGREE_BY_CONTRACT</code>. D23b n’a pas modifié cet adaptateur pour créer une cellule artificielle ; l’effet du pedigree est estimé uniquement dans le parcours régional.</p>'
        '<h3>Aveuglement et confidentialité</h3><p>La vérité privée n’apparaît dans aucune commande ni entrée d’outil. Les résultats publiés sont agrégés ; le rapport n’exporte ni identifiants individuels, ni génotypes, ni liste de références. Les artefacts restent dans <code>calibration_workspace/</code>, ignoré par Git.</p>'
        '<aside class="callout"><strong>Conséquence pour l’interprétation.</strong> La comparaison avec/sans pedigree est appariée et propre dans le parcours régional. La comparaison régional/autosomal ne l’est pas et ne doit pas servir à attribuer un écart au pedigree.</aside>',

        '<h3>Chaîne d’évaluation</h3><ol><li>10A et 11 vérifient échantillons, ordre, variants, cartes, dosages et sorties attendues.</li><li>12 et 12A produisent les haplotypes estimés.</li><li>L’oracle aligne sortie et vérité par chromosome, position, REF/ALT et sujet.</li><li>Les numérateurs et dénominateurs sont calculés par unité indépendante.</li><li>Les contrastes utilisent les mêmes tirages d’unités dans les deux termes.</li></ol>'
        '<h3>Ce que mesurent les cinq indicateurs</h3><ul><li><strong>SER :</strong> fréquence des changements d’orientation entre deux hétérozygotes évaluables successifs. Elle décrit la cohérence locale des blocs.</li><li><strong>Hamming de phase :</strong> part d’hétérozygotes mal orientés après la meilleure inversion globale H1/H2. Il décrit l’exactitude locus par locus.</li><li><strong>Motifs isolés :</strong> orientation différente de deux voisins égaux ; ce motif est descriptif et n’identifie pas sa cause.</li><li><strong>Erreur rare :</strong> Hamming restreint aux hétérozygotes dont la fréquence allélique fondatrice est ≤ 5 %.</li><li><strong>Orientation cible :</strong> proportion d’hétérozygotes cibles correctement orientés lorsque l’ancrage hors cible possède une majorité stricte.</li></ul>'
        '<p>Une SER faible peut coexister avec un Hamming élevé : un long segment peut être localement cohérent tout en étant orienté incorrectement sur de nombreux loci. C’est précisément pourquoi le rapport ne résume pas la qualité du phasage par la seule SER.</p>'
        + equations
        + '<p>Les intervalles à 95 % sont les quantiles 2,5 % et 97,5 % de 2 000 rééchantillonnages appariés. Tous les tirages de contraste SER sont valides et aucun n’a été exclu pour dénominateur nul.</p>',

        '<h3>9.1 Couverture effectivement évaluée</h3><p>La figure 1 montre la progression des loci autosomaux réellement comparés à la vérité, après exclusion des monomorphes et de la cible technique. Le passage de 59 169 à 236 483 loci augmente presque linéairement l’information disponible, sans changer les sujets.</p>'
        + figures[0]
        + '<h3>9.2 Densité et erreur régionale</h3><p>La figure 2 montre une baisse prononcée de la SER régionale entre 25 % et 50 %, puis des gains plus modestes jusqu’au maximum. Avec pedigree, les taux sont 14,95 %, 7,83 %, 5,47 % et 4,31 %. Sans pedigree, ils sont 19,50 %, 8,75 %, 5,09 % et 4,46 %.</p>'
        '<p>Les courbes se rapprochent lorsque la densité augmente. Elles se croisent légèrement à 75 %, ce qui interdit une lecture simpliste selon laquelle le pedigree améliorerait uniformément chaque cellule. Les intervalles appariés de la figure 4 doivent gouverner cette comparaison.</p>'
        + figures[1]
        + '<h3>9.3 Densité et erreur autosomale</h3><p>La figure 3 présente la relation la plus nette du banc : la SER autosomale passe de 17,87 % à 9,55 %, 6,14 % puis 4,47 %. Par rapport au maximum, les excès sont respectivement +13,41, +5,08 et +1,67 points ; leurs intervalles appariés excluent zéro.</p>'
        '<p>Cette amélioration ne doit pas être confondue avec une validation globale. Au maximum, le Hamming autosomal reste de 35,78 % et l’erreur de phase rare de 35,67 %. L’algorithme produit donc des blocs localement plus cohérents sans restituer correctement l’orientation d’une large part des hétérozygotes.</p>'
        + figures[2]
        + '<h3>9.4 Effet descriptif du retrait du pedigree</h3><p>La figure 4 représente « sans pedigree moins avec pedigree ». Une valeur positive signifie davantage de commutations sans pedigree. À 25 %, le contraste est résolu et défavorable au retrait du pedigree. À 50 %, 75 % et 100 %, les intervalles recouvrent zéro : le banc ne résout pas la direction de l’effet à ces densités.</p>'
        + figures[3]
        + '<h3>9.5 Résultats numériques complets</h3><p>Le tableau 2 fournit les numérateurs, dénominateurs et intervalles de toutes les cellules SER. Le tableau 3 complète les deux cellules primaires avec le Hamming, les motifs isolés, l’erreur rare et l’orientation cible. Le tableau 4 donne les contrastes appariés préspecifiés.</p>'
        + _table(2, ("Portée", "Densité", "Pedigree", "Statut", "Loci", "Switches/n", "SER", "IC 95 %"), performance_rows, "Taux d’erreur de commutation, effectifs et intervalles descriptifs. Source : performance_evaluations.json.")
        + _table(3, ("Primaire", "Métrique", "Numérateur/n", "Taux", "IC 95 %"), primary_rows, "Toutes les métriques primaires à densité maximale ; l’orientation cible est hors domaine de 12A. Source : performance_evaluations.json.")
        + f'<p>À densité maximale, les motifs isolés sont rares ({100*primary_regional_isolated["estimate"]:.2f} % en régional et {100*primary_autosomal_isolated["estimate"]:.2f} % en autosomal), mais le Hamming demeure élevé ({100*primary_regional_hamming["estimate"]:.2f} % et {100*primary_autosomal_hamming["estimate"]:.2f} %). L’erreur rare vaut {100*primary_regional_rare["estimate"]:.2f} % et {100*primary_autosomal_rare["estimate"]:.2f} %. La cible régionale est correctement orientée dans {primary_target["numerator"]}/{primary_target["denominator"]} cas évaluables, soit {100*primary_target["estimate"]:.2f} % ; cet effectif est trop petit et artificiel pour une interprétation biologique.</p>'
        + _table(4, ("Contraste", "Portée", "Densité", "Mode", "Différence", "IC 95 %"), contrast_rows, "Contrastes de SER appariés préspecifiés. Source : paired_contrasts.json."),

        '<h3>Contrôles fixes</h3><p>Les ensembles de densité doivent être strictement emboîtés ; la cible doit être conservée en régional ; sujets, ordre, génotypes et cartes doivent correspondre aux sources ; les sorties doivent être liées aux empreintes attendues. En mode pedigree, huit relations familiales sont attendues ; sans pedigree, aucune.</p>'
        '<p>Les douze contrôles techniques 10A/11 passent. Aucun sujet ou variant inattendu, aucune discordance de dosage et aucune divergence de carte n’ont été détectés. Les 24 composants génétiques se sont terminés avec un code de retour nul.</p>'
        + _table(5, ("Portée", "Densité", "Pedigree", "Étape", "Statut"), technical_rows, "Contrôles techniques exhaustifs 10A et 11. Source : technical_evaluations.json.")
        + '<h3 class="subsection-page">Pourquoi le verdict scientifique ne devient pas PASS</h3><p>Les critères D04 de performance sont encore <code>pilot_then_seal</code> : aucun seuil de SER, Hamming ou erreur rare n’est scellé. Le système interdit donc de convertir automatiquement une métrique mesurable en réussite scientifique. Les taux restent <code>REVIEW_REQUIRED</code>.</p>'
        '<ul><li><strong>PASS technique :</strong> fichiers, échantillons, cartes, génotypes et sorties conformes.</li><li><strong>Performance mesurable :</strong> numérateurs, dénominateurs et intervalles disponibles.</li><li><strong>Qualification scientifique :</strong> non acquise, car seuils non scellés et erreurs globales encore importantes.</li></ul>'
        '<aside class="callout"><strong>Règle de décision appliquée.</strong> Le succès de 10A/11 prouve que le banc a exécuté la comparaison prévue. Il ne prouve pas que 12/12A sont assez précis pour soutenir une analyse IBD ou un effet fondateur.</aside>',

        '<h3>Écart 1 — contexte chromosome incomplet</h3><p>La première tentative de préparation maximale a échoué avant publication parce qu’une référence attendue n’était pas transmise au contexte de chromosome. Le dossier temporaire a été nettoyé. Le défaut a été reproduit, corrigé et couvert par un test ciblé avant la préparation v2.</p>'
        '<h3>Écart 2 — chemin de schémas absent à la clôture</h3><p>Après le succès des 24 composants, la publication du reçu final a été refusée parce que le répertoire de schémas n’était pas transmis. Le contrôle a joué son rôle : aucun reçu CLOSED n’a été produit tant que le document ne pouvait pas être validé.</p>'
        '<h3>Écart 3 — racine de cellule erronée dans le finaliseur</h3><p>La première reprise bornée de publication pointait vers une racine de cellule incorrecte. Le garde a refusé un fichier manquant. La racine a été corrigée, puis les 24 reçus de composants ont été relus avant la clôture.</p>'
        '<h3>Impact sur les résultats</h3><p>Ces trois écarts concernent la préparation du contexte ou la publication des reçus, pas les génotypes, les paramètres ou les calculs de phase. Après le succès des outils, aucun composant génétique n’a été relancé : seul le finaliseur a été rejoué. Aucune valeur scientifique ni aucun seuil n’a été modifié en réponse aux résultats.</p>'
        '<p>Les essais et leurs raisons d’échec restent conservés dans le workspace pour distinguer un refus de sécurité d’un succès scientifique. Le rapport final repose exclusivement sur les manifests clos v2/v1 dont les empreintes figurent en section 14.</p>',

        '<h3>Un seul panel et un seul logiciel</h3><p>Les résultats conditionnent l’interprétation à ces 40 sujets synthétiques, à SHAPEIT5 5.1.1 et aux paramètres fixés. Ils ne quantifient ni la variabilité entre panels, ni la robustesse à d’autres références ou logiciels.</p>'
        '<h3>Profils emboîtés et intervalles descriptifs</h3><p>Les densités partagent les mêmes sujets et une grande partie de leurs marqueurs. Elles ne sont donc pas quatre réplications indépendantes. Le bootstrap apparié respecte cette dépendance, mais ses intervalles restent descriptifs et ne constituent pas une garantie de performance externe.</p>'
        '<h3>Portées régionale et autosomale non comparables causalement</h3><p>Le régional travaille dans une fenêtre cible et peut utiliser le pedigree ; l’autosomal couvre 22 chromosomes et son adaptateur n’utilise pas le pedigree. Une différence brute entre leurs taux mélange densité locale, architecture des loci, portée et adaptateur. Elle ne doit pas être étiquetée « effet du pedigree ».</p>'
        '<h3>Indicateurs incomplets pris isolément</h3><p>La SER décrit les transitions locales ; le Hamming décrit l’orientation locus par locus. Leur divergence au maximum montre qu’aucun des deux indicateurs ne suffit seul. Les erreurs rares restent également importantes, ce qui est directement pertinent pour les analyses de segments hérités.</p>'
        '<h3>Cible technique et portée biologique</h3><p>La cible APC du chromosome 5 est injectée à une fréquence artificielle pour tester l’orientation. Elle ne reproduit ni DOCK6, ni sa fréquence, ni son contexte haplotypique. Les 13 orientations correctes sur 17 ne démontrent donc aucun effet fondateur.</p>'
        '<h3>Limites de projet maintenues</h3><p>Le problème scientifique 07 identifié en D21 reste volontairement différé. D23b n’évalue pas l’IBD, les bornes de segments fondateurs ou la datation Gamma. Ces questions appartiennent aux étapes ultérieures et nécessitent leurs propres vérités et contrôles négatifs.</p>',

        '<h3>Décision sur les parties 1 à 3</h3><p>La préspecification, la construction du maximum, l’exécution des sensibilités, l’évaluation et la publication sont closes. Le banc répond à la question de densité et fournit un contraste pedigree correctement apparié dans le régional.</p>'
        '<p>La conclusion opérationnelle est de retenir la densité maximale comme configuration primaire pour les investigations suivantes, tout en conservant les profils plus faibles comme sensibilités. Ce choix repose sur la diminution observée de la SER ; il ne transforme pas le phasage en méthode qualifiée.</p>'
        '<h3>Verdicts maintenus</h3><ul><li><code>performance_qualification=REVIEW_REQUIRED</code> ;</li><li><code>calibrated=false</code> ;</li><li><code>founder_effect_proven=false</code> ;</li><li>cible chr5 classée contrôle technique, jamais preuve DOCK6.</li></ul>'
        '<h3>Condition d’ouverture de D24</h3><p>D24 n’est pas lancé. Avant sa partie 4, il faut démontrer l’existence d’une vérité indépendante de segments fondateurs contenant des bornes vraies, des contrôles négatifs et au moins trois unités indépendantes. Si cette vérité n’existe pas ou n’est pas suffisante, D24 doit rester bloqué plutôt que réutiliser D23 comme substitut.</p>'
        '<h3>Travail volontairement reporté</h3><p>Les parties 4 et 5 devront traiter séparément l’évaluateur d’haplotype fondateur puis la suite de qualification prévue par le plan. Elles ne doivent ni modifier rétroactivement les seuils D23b, ni résoudre silencieusement le problème 07.</p>',

        '<p>Chaque source ci-dessous est identifiée par son SHA-256. Une modification d’un manifest, d’une configuration, du paquet de puce ou des cartes invaliderait la chaîne de reproduction et exigerait un nouveau run, pas une réédition silencieuse du rapport.</p>'
        + _table(6, ("Objet", "SHA-256", "Rôle"), provenance_rows, "Empreintes de la chaîne D23b. Sources : manifests clos et préparation maximale.")
        + '<h3 class="subsection-page">Artefacts du rapport</h3><p><code>report_data.json</code> contient les agrégats effectivement rendus : préparation, contrôles techniques, performances et contrastes. <code>checksums.sha256</code> lie ce fichier à l’HTML autonome. Le PDF possède sa propre empreinte.</p>'
        '<p>Le rapport ne contient aucune liste de sujets, aucun génotype et aucune vérité locus par locus. Les données sources restent dans <code>calibration_workspace/</code>, hors Git. Les modules, schémas, configurations et tests reproductibles sont les seuls éléments destinés au dépôt.</p>'
        '<aside class="callout"><strong>Reproductibilité.</strong> Évaluation close : <span class="hash">' + escape(source_sha) + '</span>. Préparation maximale : <span class="hash">' + escape(manifest["preparation_manifest_sha256"]) + '</span>. Ces empreintes sont les ancres de toute comparaison future.</aside>',
    ]


def build_report(*, evaluation_root: Path, preparation_root: Path, template_path: Path, output_root: Path) -> Path:
    """Construit l'HTML autonome sans modifier les deux sources closes."""
    destination = new_workspace_run(output_root.absolute().parent, output_root.name)
    manifest = _checked_evaluation(evaluation_root)
    source_sha = sha256_file(evaluation_root / "run_index.json")
    preparation = json.loads((preparation_root / "preparation_manifest.json").read_text())
    if sha256_file(preparation_root / "preparation_manifest.json") != manifest["preparation_manifest_sha256"]:
        raise ValueError("d23b_report_preparation_anchor_invalid")
    technical = json.loads((evaluation_root / "technical_evaluations.json").read_text())
    performance = json.loads((evaluation_root / "performance_evaluations.json").read_text())
    contrasts = json.loads((evaluation_root / "paired_contrasts.json").read_text())
    report_data = {
        "manifest": manifest,
        "preparation": {
            "counts": preparation["counts"],
            "chromosomes": preparation["chromosomes"],
            "audit": preparation["audit"],
            "tool_versions": preparation["tool_versions"],
            "performance_status": preparation["performance_status"],
            "configuration_sha256": preparation["configuration_sha256"],
            "chip_package_sha256": preparation["chip_package_sha256"],
            "map_archive_sha256": preparation["map_archive_sha256"],
        },
        "technical": technical,
        "performance": performance,
        "contrasts": contrasts,
    }
    head = re.split(r"<body\b", template_path.read_text(), maxsplit=1)[0]
    if "</head>" not in head:
        raise ValueError("d23b_report_template_invalid")
    title = "CAL-DEV-D23B — Sensibilité densité × pedigree"
    head = re.sub(r"<title>.*?</title>", f"<title>{title}</title>", head, flags=re.S)
    head = head.replace("</head>", '<style>.hash{overflow-wrap:anywhere}.scientific-figure,.equation,table{break-inside:avoid}.scientific-figure{margin:1em 0}.equation{border-left:3px solid #286d91;padding:.7em;background:#f2f6f8;margin:1em 0}.math{font:1em Georgia,serif;line-height:1.8}table{font-size:.68em}.caption{font-size:.78em;color:#46545d}.subsection-page{page-break-before:always!important;break-before:page!important}main,main>section{page-break-inside:auto!important;break-inside:auto!important}main>section>h2{page-break-before:always!important;break-before:page!important}main p{orphans:3;widows:3}</style></head>')
    toc = '<nav class="toc-page"><h2>Sommaire</h2><ol>' + ''.join(f'<li><a class="toc-entry" href="#s{index}"><span class="toc-title">{index}. {escape(heading)}</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>' for index, heading in enumerate(HEADINGS, 1)) + '</ol></nav>'
    body = f'<body lang="fr"><section class="cover"><h1>{title}</h1><p>Effet_fondateur2 · {escape(manifest["created_at"][:10])} · version 1.1.0</p><p>Investigation complémentaire D23 — parties 1 à 3</p><p>Version narrative enrichie : méthode, résultats, interprétation et limites</p></section>' + toc
    body += '<main>'
    body += ''.join(f'<section id="s{index}"><h2>{index}. {escape(heading)}</h2>{content}</section>' for index, (heading, content) in enumerate(zip(HEADINGS, _sections(manifest, technical, performance, contrasts, preparation, source_sha), strict=True), 1))
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        directory = temporary / "CAL-DEV-D23B"
        directory.mkdir(mode=0o700)
        html = directory / "index.html"
        html.write_text(head + body + '</main><footer class="page-footer"><span>CAL-DEV-D23B · Effet_fondateur2</span></footer></body></html>')
        write_json(directory / "report_data.json", report_data)
        (directory / "checksums.sha256").write_text("".join(f"{sha256_file(path)}  {path.name}\n" for path in (html, directory / "report_data.json")))
        if sha256_file(evaluation_root / "run_index.json") != source_sha:
            raise ValueError("d23b_report_source_changed")
        temporary.rename(destination)
        return destination / "CAL-DEV-D23B/index.html"
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("evaluation-root", "preparation-root", "template-path", "output-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    build_report(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()
