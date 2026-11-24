"""Rapport final R14b de l'amendement R13 limité à l'intégration.

Le rapport consolide des preuves synthétiques closes. Il produit un HTML
autonome paginé et son PDF, sans lire de donnée réelle, rouvrir la vérité ou
relancer un outil scientifique. Les statistiques dérivées sont descriptives
et leurs hypothèses sont affichées dans le rapport.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import html
import json
from math import sqrt
from pathlib import Path
import re
import shutil
import statistics
import sys
from typing import Any, Mapping, Sequence

from effet_fondateur.calibration.contracts import file_sha256, load_json
from effet_fondateur.calibration.r14_report import _run_chrome


REPORT_ID = "R14B-INTEGRATION-CLOSURE-001"
REPORT_TITLE = "Rapport final R14b — intégration scientifique des outils et clôture amendée"
REPORT_SCOPE = "R14B_TOOL_INTEGRATION_CLOSURE"
FIXED_DATE = "15 novembre 2026"
READINESS_PATH = Path("calibration/readiness/R14B_INTEGRATION_CLOSURE.json")
R13_AMENDMENT_PATH = Path("calibration/readiness/R13_INTEGRATION_AMENDMENT.json")
R13_CALL_RELATIVE = Path(
    "outputs/external_calls/SHAPEIT-DMAX-CHR01-SEED01/attempt-03/receipt.json"
)
R13_LOG_RELATIVE = Path(
    "outputs/external_calls/SHAPEIT-DMAX-CHR01-SEED01/attempt-03/shapeit.log"
)
QUALIFICATION_FLAGS = {
    "calibrated": False,
    "founder_effect_proven": False,
    "production_promotion": False,
}


class R14bReportError(RuntimeError):
    """Source incohérente, rendu incomplet ou décision hors périmètre."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R14bReportError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _source_record(repository_root: Path, path: Path, role: str) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    try:
        relative = resolved.relative_to(repository_root).as_posix()
    except ValueError:
        relative = resolved.as_posix()
    return {
        "relative_path": relative,
        "role": role,
        "sha256": file_sha256(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def wilson_two_sided_95(successes: int, trials: int) -> tuple[float, float]:
    """Intervalle descriptif de Wilson bilatéral à 95 %.

    L'interprétation binomiale exige des unités indépendantes. R14b l'affiche
    comme hypothèse non démontrée pour les duplicats R12 et n'en tire donc
    aucune décision de qualification.
    """
    if (
        isinstance(successes, bool)
        or isinstance(trials, bool)
        or not isinstance(successes, int)
        or not isinstance(trials, int)
        or trials <= 0
        or successes < 0
        or successes > trials
    ):
        _fail("r14b_binomial_counts_invalid")
    z = 1.959963984540054
    estimate = successes / trials
    denominator = 1.0 + z * z / trials
    center = (estimate + z * z / (2.0 * trials)) / denominator
    half_width = (
        z
        * sqrt(
            estimate * (1.0 - estimate) / trials
            + z * z / (4.0 * trials * trials)
        )
        / denominator
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def summarize_ibd_measurements(
    measurements: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Résume les trois répétitions de chaque cellule outil–densité."""
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for item in measurements:
        grouped[(str(item.get("tool")), str(item.get("density")))].append(item)
    expected = {
        (tool, density)
        for tool in ("HAP_IBD", "REFINED_IBD")
        for density in ("D_MAX", "D_MID", "D_LOW")
    }
    if set(grouped) != expected or any(len(items) != 3 for items in grouped.values()):
        _fail("r14b_ibd_measurement_cells_incomplete")
    summaries = []
    for density in ("D_MAX", "D_MID", "D_LOW"):
        for tool in ("HAP_IBD", "REFINED_IBD"):
            items = grouped[(tool, density)]
            wall = [float(item["wall_seconds"]) for item in items]
            rss = [float(item["maximum_rss_bytes"]) / 2**20 for item in items]
            summaries.append(
                {
                    "tool": tool,
                    "density": density,
                    "n": len(items),
                    "wall_mean_seconds": statistics.fmean(wall),
                    "wall_min_seconds": min(wall),
                    "wall_max_seconds": max(wall),
                    "rss_mean_mib": statistics.fmean(rss),
                    "rss_min_mib": min(rss),
                    "rss_max_mib": max(rss),
                }
            )
    return summaries


def _validate_sources(
    repository_root: Path,
    r13_run_root: Path,
    r12_ibd_root: Path,
    r11b_root: Path,
    template_path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    amendment_path = repository_root / R13_AMENDMENT_PATH
    amendment = load_json(amendment_path)
    if amendment.get("status") != "CLOSED_LIGHTWEIGHT_INTEGRATION_VERIFIED":
        _fail("r14b_r13_amendment_not_closed")
    if amendment.get("shapeit5_integration", {}).get("verdict") != "INTEGRATION_VERIFIED":
        _fail("r14b_shapeit_integration_not_verified")
    if amendment.get("truth_opened") is not False:
        _fail("r14b_truth_opened")
    if amendment.get("qualification_flags") != QUALIFICATION_FLAGS:
        _fail("r14b_amendment_flags_invalid")
    campaign = amendment.get("original_campaign", {})
    expected_campaign = {
        "historical_call_count": 1_180,
        "retained_completed_calls": 1,
        "withdrawn_unexecuted_calls": 1_179,
        "newly_authorized_calls": 0,
        "disposition": "SUPERSEDED_NOT_EXECUTED",
    }
    if any(campaign.get(key) != value for key, value in expected_campaign.items()):
        _fail("r14b_campaign_disposition_invalid")

    call_receipt_path = r13_run_root / R13_CALL_RELATIVE
    call_log_path = r13_run_root / R13_LOG_RELATIVE
    call_receipt = load_json(call_receipt_path)
    if (
        call_receipt.get("state") != "COMPLETED"
        or call_receipt.get("measurement", {}).get("returncode") != 0
        or call_receipt.get("truth_available_to_worker") is not False
    ):
        _fail("r14b_shapeit_call_invalid")
    if amendment["evidence_sha256"].get(
        f"r13_run/{R13_CALL_RELATIVE.as_posix()}"
    ) != file_sha256(call_receipt_path):
        _fail("r14b_shapeit_receipt_hash_mismatch")
    if amendment["evidence_sha256"].get(
        f"r13_run/{R13_LOG_RELATIVE.as_posix()}"
    ) != file_sha256(call_log_path):
        _fail("r14b_shapeit_log_hash_mismatch")

    ibd_evaluation_path = r12_ibd_root / "evaluations/explicit_ibd_r12.json"
    ibd_evaluation = load_json(ibd_evaluation_path)
    ibd_measurements_path = r11b_root / "evaluations/ibd_batch_measurements.json"
    ibd_measurements = load_json(ibd_measurements_path)
    r11b_decision_path = r11b_root / "evaluations/R11B_RESOURCE_DECISION.json"
    r11b_decision = load_json(r11b_decision_path)
    if ibd_measurements.get("probe_count") != 18:
        _fail("r14b_ibd_probe_count_mismatch")
    if r11b_decision.get("status") != "RESOURCE_READY_FOR_NEW_SCIENTIFIC_RELEASE":
        _fail("r14b_r11b_resource_decision_invalid")

    tracked_sources = [
        ("calibration/PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md", "plan normatif R00–R14b"),
        ("calibration/CONVENTIONS_DEVELOPPEMENT.md", "contrat permanent des rapports"),
        (R13_AMENDMENT_PATH.as_posix(), "reçu principal de l'amendement R13"),
        ("calibration/configurations/r13_integration_amendment_v1.yaml", "contrat R13 amendé"),
        ("calibration/readiness/R07_SCIENTIFIC_DIMENSIONING_V2.json", "dimensionnement R07 v2"),
        ("calibration/readiness/R11B_RESOURCE_CONFIRMATION.json", "reçu de ressources R11b"),
        ("calibration/readiness/R12C_G5_READY.json", "reçu G5 R12c historique"),
        ("calibration/readiness/R14_TECHNICAL_CLOSURE.json", "clôture technique R14 antérieure"),
        ("src/effet_fondateur/calibration/r11b_resources.py", "formules de projection de ressources"),
        ("src/effet_fondateur/calibration/r14b_report.py", "constructeur et statistiques R14b"),
    ]
    sources = [
        _source_record(repository_root, repository_root / relative, role)
        for relative, role in tracked_sources
    ]
    sources.extend(
        (
            _source_record(repository_root, call_receipt_path, "reçu SHAPEIT5 grandeur réelle"),
            _source_record(repository_root, call_log_path, "journal SHAPEIT5 grandeur réelle"),
            _source_record(repository_root, ibd_evaluation_path, "contrôles positifs IBD R12"),
            _source_record(repository_root, ibd_measurements_path, "18 mesures de ressources IBD R11b"),
            _source_record(repository_root, r11b_decision_path, "projection et décision de ressources R11b"),
            _source_record(repository_root, template_path, "modèle autonome report-formater"),
        )
    )
    return (
        {
            "amendment": amendment,
            "call_receipt": call_receipt,
            "ibd_evaluation": ibd_evaluation,
            "ibd_measurements": ibd_measurements,
            "r11b_decision": r11b_decision,
        },
        sources,
    )


def _fr(value: float, digits: int = 2) -> str:
    return f"{value:,.{digits}f}".replace(",", "X").replace(".", ",").replace("X", " ")


def _source_rows(sources: Sequence[Mapping[str, Any]]) -> str:
    rows = []
    for index, item in enumerate(sources, start=1):
        path = Path(str(item["relative_path"]))
        display = "/".join(path.parts[-3:]) if len(path.parts) >= 3 else path.as_posix()
        rows.append(
            "<tr>"
            f'<td style="width:2.5rem;white-space:nowrap">{index}</td>'
            f"<td><code>{html.escape(display)}</code></td>"
            f"<td>{html.escape(str(item['role']))}</td>"
            f"<td><code>{item['sha256'][:12]}…</code></td>"
            "</tr>"
        )
    return "".join(rows)


def _ibd_resource_svg(summaries: Sequence[Mapping[str, Any]]) -> str:
    maximum = max(float(item["wall_max_seconds"]) for item in summaries)
    rows = []
    colors = {"HAP_IBD": "#1f4e79", "REFINED_IBD": "#e85d38"}
    names = {"HAP_IBD": "Hap-IBD", "REFINED_IBD": "Refined IBD"}
    y = 48
    for item in summaries:
        start = 180 + 360 * float(item["wall_min_seconds"]) / maximum
        end = 180 + 360 * float(item["wall_max_seconds"]) / maximum
        mean = 180 + 360 * float(item["wall_mean_seconds"]) / maximum
        label = f"{item['density']} · {names[str(item['tool'])]}"
        rows.append(
            f'<text x="10" y="{y + 4}" font-size="10">{label}</text>'
            f'<line x1="{start:.1f}" y1="{y}" x2="{end:.1f}" y2="{y}" '
            f'stroke="{colors[str(item["tool"])]}" stroke-width="4"/>'
            f'<circle cx="{mean:.1f}" cy="{y}" r="5" fill="{colors[str(item["tool"])]}"/>'
            f'<text x="555" y="{y + 4}" font-size="10">'
            f'{_fr(float(item["wall_mean_seconds"]), 2)} s</text>'
        )
        y += 30
    return "".join(rows)


def _report_body(data: Mapping[str, Any], sources: Sequence[Mapping[str, Any]]) -> str:
    amendment = data["amendment"]
    shapeit = amendment["shapeit5_integration"]
    measurement = data["call_receipt"]["measurement"]
    ibd_tools = data["ibd_evaluation"]["tools"]
    summaries = summarize_ibd_measurements(data["ibd_measurements"]["measurements"])
    r11b = data["r11b_decision"]
    hap_low, hap_high = wilson_two_sided_95(
        ibd_tools["hap_ibd"]["designated_duplicate_detected"],
        ibd_tools["hap_ibd"]["designated_duplicate_expected"],
    )
    refined_low, refined_high = wilson_two_sided_95(
        ibd_tools["refined_ibd"]["designated_duplicate_detected"],
        ibd_tools["refined_ibd"]["designated_duplicate_expected"],
    )
    wall_hours = float(measurement["wall_seconds"]) / 3_600
    user_hours = float(measurement["user_cpu_seconds"]) / 3_600
    system_hours = float(measurement["system_cpu_seconds"]) / 3_600
    rss_gib = float(measurement["maximum_rss_native"]) / 2**30
    time_max = max(wall_hours, user_hours, system_hours)
    hap_x = 110 + 420 * 0.92
    hap_x_low = 110 + 420 * hap_low
    hap_x_high = 110 + 420 * hap_high
    refined_x = 110 + 420 * 0.88
    refined_x_low = 110 + 420 * refined_low
    refined_x_high = 110 + 420 * refined_high
    source_rows = _source_rows(sources)
    ibd_rows = "".join(
        "<tr>"
        f"<td>{item['density']}</td>"
        f"<td>{'Hap-IBD' if item['tool'] == 'HAP_IBD' else 'Refined IBD'}</td>"
        f"<td>{item['n']}</td>"
        f"<td>{_fr(item['wall_mean_seconds'], 3)} [{_fr(item['wall_min_seconds'], 3)}–{_fr(item['wall_max_seconds'], 3)}]</td>"
        f"<td>{_fr(item['rss_mean_mib'], 1)} [{_fr(item['rss_min_mib'], 1)}–{_fr(item['rss_max_mib'], 1)}]</td>"
        "</tr>"
        for item in summaries
    )
    resource_rows = "".join(
        "<tr>"
        f"<td>{html.escape(str(item['check']))}</td>"
        f"<td>{_fr(float(item['projected']), 3)}</td>"
        f"<td>{_fr(float(item['ceiling']), 3)}</td>"
        f"<td>{html.escape(str(item['status']))}</td>"
        "</tr>"
        for item in r11b["checks"]
    )
    return f"""
<body lang="fr">
  <section class="cover" id="cover">
    <div class="cover-logo"><span style="font-family:Georgia,serif;font-size:1.45em;font-weight:700;color:var(--accent);letter-spacing:.08em;">EFFET FONDATEUR 2</span></div>
    <div class="cover-accent-band"></div>
    <div class="cover-center">
      <h1>{REPORT_TITLE}</h1>
      <div class="cover-subtitle">Rapport scientifique détaillé · données synthétiques uniquement</div>
    </div>
    <div class="cover-meta"><table>
      <tr><td>Décision</td><td>Intégration des outils vérifiée ; performance complète non évaluée</td></tr>
      <tr><td>Release</td><td><code>{amendment['release_id']}</code></td></tr>
      <tr><td>Date</td><td>{FIXED_DATE}</td></tr>
      <tr><td>Version</td><td>1.0.0 · {REPORT_ID}</td></tr>
      <tr><td>Format</td><td>HTML autonome et PDF A4 · report-formater</td></tr>
    </table></div>
  </section>

  <nav class="toc-page" id="toc"><h2>Sommaire</h2><ol>
    <li><a class="toc-entry" href="#s1"><span class="toc-title">1. Résumé exécutif</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s2"><span class="toc-title">2. Question, portée et décisions</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s3"><span class="toc-title">3. Sources et hiérarchie des preuves</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s4"><span class="toc-title">4. Design R12c–R14b et aveuglement</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s5"><span class="toc-title">5. Méthodes statistiques et équations</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s6"><span class="toc-title">6. Résultat principal : SHAPEIT5</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s7"><span class="toc-title">7. Résultats d'intégration IBD</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s8"><span class="toc-title">8. Ressources et faisabilité</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s9"><span class="toc-title">9. Disposition des checkpoints</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s10"><span class="toc-title">10. Interprétation, risques et limites</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s11"><span class="toc-title">11. Décision finale R14b</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    <li><a class="toc-entry" href="#s12"><span class="toc-title">12. Annexes, références et traçabilité</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
  </ol></nav>

  <main>
    <section id="s1"><div class="summary"><h2>1. Résumé exécutif</h2><ul>
      <li><strong>SHAPEIT5 est correctement intégré au pipeline.</strong> La preuve principale est un appel grandeur réelle terminé sur 4 000 sujets, 3 000 références, 46 524 variants et 200 trios, sans erreur mendélienne.</li>
      <li><strong>Hap-IBD et Refined IBD disposent de preuves locales d'intégration.</strong> Les deux adaptateurs ont produit des segments et détecté respectivement 23/25 et 22/25 duplicats synthétiques désignés.</li>
      <li><strong>La validité intrinsèque des méthodes est référencée, pas redémontrée.</strong> Les publications fondatrices sont liées ; le dépôt démontre localement le raccordement, les paramètres, les sorties et la provenance.</li>
      <li><strong>La campagne exhaustive a été amendée.</strong> Un appel est retenu, 1 179 appels non exécutés sont retirés, et aucun nouvel appel n'est autorisé.</li>
      <li><strong>La portée reste strictement limitée.</strong> La performance R07/R08 complète, l'effet fondateur et la promotion en production ne sont pas établis.</li>
    </ul></div>
    <div class="stat-grid">
      <div class="stat-card"><div class="num">4 000 + 3 000</div><div class="label">sujets + références SHAPEIT5</div></div>
      <div class="stat-card"><div class="num">46 524</div><div class="label">variants produits</div></div>
      <div class="stat-card"><div class="num">0</div><div class="label">erreur mendélienne</div></div>
      <div class="stat-card"><div class="num">1 / 1 180</div><div class="label">appel retenu / plan historique</div></div>
    </div>
    <aside class="callout"><span class="label">Conclusion en une phrase</span>R14b clôt une vérification d'intégration crédible et traçable ; il ne transforme pas une preuve de bon fonctionnement local en validation complète des performances scientifiques.</aside></section>

    <section id="s2"><h2>2. Question, portée et décisions</h2>
      <p>La question R14b n'est pas « SHAPEIT5 ou les méthodes IBD sont-ils valides en général ? », question traitée par la littérature, mais « le dépôt appelle-t-il correctement ces outils dans un cadre reproductible et conforme à ses contrats ? ». La réponse est positive pour l'intégration, avec des niveaux de preuve distincts selon l'outil.</p>
      <table><caption>Table 1 : Questions de décision R14b et réponses. Sources : reçu R13 amendé, contrôles R12 et plan R00–R14b.</caption><thead><tr><th>Question</th><th>Réponse</th><th>Portée</th></tr></thead><tbody>
        <tr><td>Intégration SHAPEIT5</td><td><strong>VÉRIFIÉE</strong></td><td>Appel grandeur réelle, sortie BCF et index intègres</td></tr>
        <tr><td>Intégration Hap-IBD</td><td><strong>SOUTENUE</strong></td><td>Contrôle positif synthétique et sorties non vides</td></tr>
        <tr><td>Intégration Refined IBD</td><td><strong>SOUTENUE</strong></td><td>Contrôle positif synthétique et sorties non vides</td></tr>
        <tr><td>Performance complète R07/R08</td><td><strong>NON ÉVALUÉE</strong></td><td>Pas de grille exhaustive de sensibilité/spécificité</td></tr>
        <tr><td>Effet fondateur</td><td><strong>NON DÉMONTRÉ</strong></td><td>Aucune donnée d'étude réelle analysée</td></tr>
        <tr><td>Promotion en production</td><td><strong>NON AUTORISÉE</strong></td><td>Qualification complète absente</td></tr>
      </tbody></table>
      <aside class="callout"><span class="label">Données utilisées</span>Toutes les preuves sont synthétiques. Aucun fichier sous <code>data/input/</code>, aucun génotype réel et aucune vérité privée n'apparaissent dans ce rapport.</aside>
    </section>

    <section id="s3"><h2>3. Sources et hiérarchie des preuves</h2>
      <p>Les résultats bruts et reçus clos priment sur les résumés de session. Les contrats et le plan déterminent la portée ; les publications déterminent la validité externe des méthodes. La Figure 1 montre cette hiérarchie sans confondre preuve externe et preuve locale.</p>
      <figure role="img" aria-labelledby="fig1cap"><svg viewBox="0 0 680 190" width="100%" aria-label="Hiérarchie des preuves R14b">
        <defs><marker id="a1" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#1f4e79"/></marker></defs>
        <text x="8" y="18" font-size="13" font-weight="700">Deux familles de preuves, une décision bornée</text>
        <rect x="18" y="42" width="190" height="54" rx="6" fill="#eaf0f6" stroke="#1f4e79"/><text x="113" y="63" text-anchor="middle" font-size="11"><tspan x="113">Publications méthodologiques</tspan><tspan x="113" dy="16">validité externe des outils</tspan></text>
        <rect x="246" y="42" width="190" height="54" rx="6" fill="#eaf0f6" stroke="#1f4e79"/><text x="341" y="63" text-anchor="middle" font-size="11"><tspan x="341">Preuves locales closes</tspan><tspan x="341" dy="16">appels, sorties, hashes</tspan></text>
        <rect x="474" y="42" width="190" height="54" rx="6" fill="#dcefe2" stroke="#29633a"/><text x="569" y="63" text-anchor="middle" font-size="11"><tspan x="569">R14b</tspan><tspan x="569" dy="16">intégration vérifiée</tspan></text>
        <line x1="208" y1="69" x2="242" y2="69" stroke="#1f4e79" stroke-width="2" marker-end="url(#a1)"/><line x1="436" y1="69" x2="470" y2="69" stroke="#1f4e79" stroke-width="2" marker-end="url(#a1)"/>
        <rect x="246" y="122" width="418" height="42" rx="6" fill="#fff1c9" stroke="#b07800"/><text x="455" y="139" text-anchor="middle" font-size="10"><tspan x="455">Limite explicite : performance exhaustive, effet fondateur</tspan><tspan x="455" dy="14">et promotion restent hors de la conclusion</tspan></text>
      </svg><figcaption id="fig1cap">Figure 1 : Hiérarchie des preuves de R14b. Sources : contrat d'amendement R13, reçu R13 et références bibliographiques.</figcaption></figure>
      <table><caption>Table 2 : Inventaire des sources et empreintes. Source : inventaire R14b généré au moment du rendu.</caption><thead><tr><th style="width:2.5rem;white-space:nowrap">#</th><th>Source</th><th>Rôle</th><th>SHA-256</th></tr></thead><tbody>{source_rows}</tbody></table>
    </section>

    <section id="s4"><h2>4. Design R12c–R14b et aveuglement</h2>
      <p>R12c a gelé une partition <code>VALIDATION_BLIND</code>, un matérialiseur et un worker unitaire. Le worker reçoit seulement les entrées aveugles, la carte, le pedigree et les paramètres de l'appel. La vérité privée n'est pas transmise. R13 a exécuté un checkpoint, puis l'amendement a fermé l'autorisation de calcul sans réécrire le registre historique.</p>
      <table><caption>Table 3 : Séparation des responsabilités. Sources : manifeste R12c, reçu du worker et overlay R13.</caption><thead><tr><th>Composant</th><th>Accès</th><th>Sortie contrôlée</th></tr></thead><tbody>
        <tr><td>Matérialiseur</td><td>Contrat scellé et graines déterministes</td><td>VCF d'étude, référence, carte et pedigree</td></tr>
        <tr><td>Worker SHAPEIT5</td><td>Vue aveugle uniquement</td><td>BCF, CSI, journal et reçu de ressources</td></tr>
        <tr><td>Registre</td><td>Contrats des 1 180 appels</td><td>Transitions et empreintes des sorties</td></tr>
        <tr><td>Amendement</td><td>Preuves closes, sans vérité</td><td>1 retenu, 1 179 retirés, 0 autorisé</td></tr>
        <tr><td>R14b</td><td>Résultats agrégés seulement</td><td>Décision, HTML et PDF</td></tr>
      </tbody></table>
      <pre><code>preuve_intégration =
    identité_et_version_outil
  ∧ entrées_et_paramètres_conformes
  ∧ exécution_sans_shell_et_code_retour_0
  ∧ contrôles_synthétiques_cohérents
  ∧ sorties_indexées_et_empreintées
  ∧ vérité_indisponible_au_worker</code></pre>
      <p>Cette expression est une règle logique de gouvernance, pas une équation statistique. Chaque terme doit être vrai ; aucun score compensatoire n'est calculé.</p>
    </section>

    <section id="s5"><h2>5. Méthodes statistiques et équations</h2>
      <h3>5.1 Proportion observée et intervalle descriptif</h3>
      <pre><code>p̂ = k / n
centre = (p̂ + z² / 2n) / (1 + z² / n)
demi-largeur =
  z × √[p̂(1−p̂)/n + z²/(4n²)] / (1 + z²/n)
IC95_Wilson = [centre − demi-largeur ; centre + demi-largeur]

k : duplicats désignés détectés
n : duplicats désignés attendus
z : 1,959964 pour un intervalle bilatéral à 95 %</code></pre>
      <p>L'intervalle de Wilson est affiché pour rendre visible l'incertitude liée à <code>n=25</code>. Il suppose des essais binomiaux indépendants. L'indépendance des duplicats R12 n'étant pas établie comme unité scientifique, cet intervalle reste descriptif et ne sert à aucun verdict de qualification.</p>
      <h3>5.2 Résumé des mesures de ressources</h3>
      <pre><code>moyenne = (1/n) × Σ xᵢ
étendue_observée = [min(xᵢ) ; max(xᵢ)]
n = 3 répétitions par cellule outil–densité</code></pre>
      <p>La moyenne et l'étendue décrivent les 18 sondes R11b. Une étendue n'est pas un intervalle de confiance. Les six cellules sont comparées sans supposer de distribution commune entre outils ou densités.</p>
      <h3>5.3 Borne conservatrice R11b</h3>
      <pre><code>B95_log = exp[moyenne(log x)
              + t₀,₉₅,n−1 × s_log × √(1 + 1/n)]
valeur_sûre = 2 × B95_log

x : temps, CPU ou RSS strictement positif
n : 3 répétitions par cellule IBD
2 : facteur de sécurité préspecifié</code></pre>
      <p>Pour SHAPEIT5, R11b emploie une régression log–log sur l'effectif total et une borne prédictive unilatérale, également multipliée par deux. Ces projections décrivent la faisabilité de calcul ; elles ne mesurent ni exactitude de phase ni performance biologique.</p>
    </section>

    <section id="s6"><h2>6. Résultat principal : intégration SHAPEIT5</h2>
      <p>Le checkpoint <code>{shapeit['call_id']}</code>, tentative {shapeit['attempt']}, est la preuve principale de R13 amendé. SHAPEIT5 {shapeit['tool_version']} a lu la cohorte, la référence, la carte génétique et le pedigree, puis a produit un BCF indexé avec les dimensions attendues.</p>
      <div class="stat-grid">
        <div class="stat-card"><div class="num">{_fr(wall_hours, 2)} h</div><div class="label">temps mural observé</div></div>
        <div class="stat-card"><div class="num">{_fr(rss_gib, 2)} Gio</div><div class="label">RSS maximal observé</div></div>
        <div class="stat-card"><div class="num">200</div><div class="label">trios contrôlés</div></div>
        <div class="stat-card"><div class="num">0</div><div class="label">erreur mendélienne</div></div>
      </div>
      <figure role="img" aria-labelledby="fig2cap"><svg viewBox="0 0 680 210" width="100%" aria-label="Temps observés du checkpoint SHAPEIT5">
        <text x="8" y="18" font-size="13" font-weight="700">Temps observés du checkpoint SHAPEIT5</text>
        <g font-size="11"><text x="20" y="63">Mur</text><rect x="105" y="46" width="{430 * wall_hours / time_max:.1f}" height="22" fill="#1f4e79"/><text x="550" y="63">{_fr(wall_hours, 3)} h</text>
        <text x="20" y="108">CPU utilisateur</text><rect x="105" y="91" width="{430 * user_hours / time_max:.1f}" height="22" fill="#3a7ca5"/><text x="550" y="108">{_fr(user_hours, 3)} h</text>
        <text x="20" y="153">CPU système</text><rect x="105" y="136" width="{430 * system_hours / time_max:.1f}" height="22" fill="#e85d38"/><text x="550" y="153">{_fr(system_hours, 3)} h</text></g>
        <text x="20" y="192" font-size="10" fill="#555">Un seul thread ; les valeurs proviennent du reçu de la tentative 3.</text>
      </svg><figcaption id="fig2cap">Figure 2 : Temps mur, CPU utilisateur et CPU système du checkpoint SHAPEIT5. Source : reçu <code>attempt-03</code> ; n = 1 appel.</figcaption></figure>
      <table><caption>Table 4 : Contrôles observables du checkpoint SHAPEIT5. Sources : reçu, journal, BCF et CSI de la tentative 3.</caption><thead><tr><th>Contrôle</th><th>Valeur</th><th>Interprétation</th></tr></thead><tbody>
        <tr><td>Code retour</td><td>0</td><td>Exécution terminée sans erreur externe</td></tr>
        <tr><td>Sujets / références</td><td>{shapeit['study_samples']} / {shapeit['reference_samples']}</td><td>Dimensions conformes au contrat</td></tr>
        <tr><td>Variants entrée / sortie</td><td>{shapeit['variants']} / {shapeit['output_variants']}</td><td>Aucune perte de ligne observée</td></tr>
        <tr><td>Trios / erreurs mendéliennes</td><td>{shapeit['trios']} / {shapeit['mendel_errors']}</td><td>Pedigree accepté et cohérent</td></tr>
        <tr><td>Shell</td><td>Non</td><td>Commande structurée, sans interpolation shell</td></tr>
        <tr><td>Vérité accessible</td><td>Non</td><td>Aveuglement maintenu</td></tr>
        <tr><td>BCF / CSI / journal</td><td>Empreintes valides</td><td>Sorties auditables et non altérées</td></tr>
      </tbody></table>
      <aside class="callout"><span class="label">Ce que cette preuve ne dit pas</span>Un appel réussi ne mesure pas le taux d'erreur de commutation, l'exactitude de phase par fréquence allélique ni la robustesse sur 22 autosomes. Ces performances restent non évaluées.</aside>
    </section>

    <section id="s7"><h2>7. Résultats d'intégration Hap-IBD et Refined IBD</h2>
      <p>Les preuves locales IBD proviennent du candidat R12 clos. Elles montrent que les adaptateurs exécutent les deux outils, normalisent leurs sorties et retrouvent la majorité des duplicats synthétiques désignés. Elles ne couvrent pas la grille complète de longueurs de segments, de faux positifs, d'erreur de bornes ou de concordance.</p>
      <figure role="img" aria-labelledby="fig3cap"><svg viewBox="0 0 680 210" width="100%" aria-label="Détection descriptive des duplicats désignés">
        <text x="8" y="18" font-size="13" font-weight="700">Proportion observée et IC95 de Wilson descriptif</text>
        <line x1="110" y1="165" x2="530" y2="165" stroke="#777"/><g font-size="9" fill="#555"><text x="105" y="182">0 %</text><text x="310" y="182">50 %</text><text x="515" y="182">100 %</text></g>
        <text x="15" y="70" font-size="11">Hap-IBD</text><line x1="{hap_x_low:.1f}" y1="66" x2="{hap_x_high:.1f}" y2="66" stroke="#1f4e79" stroke-width="4"/><circle cx="{hap_x:.1f}" cy="66" r="7" fill="#1f4e79"/><text x="545" y="70" font-size="10">23/25 · 92 %</text>
        <text x="15" y="122" font-size="11">Refined IBD</text><line x1="{refined_x_low:.1f}" y1="118" x2="{refined_x_high:.1f}" y2="118" stroke="#e85d38" stroke-width="4"/><circle cx="{refined_x:.1f}" cy="118" r="7" fill="#e85d38"/><text x="545" y="122" font-size="10">22/25 · 88 %</text>
      </svg><figcaption id="fig3cap">Figure 3 : Détection des duplicats désignés avec IC95 de Wilson bilatéral. Source : <code>explicit_ibd_r12.json</code> ; n = 25 par outil. Intervalles descriptifs, indépendance non démontrée.</figcaption></figure>
      <table><caption>Table 5 : Contrôles positifs IBD et limites. Source : évaluation IBD R12 close.</caption><thead><tr><th>Outil</th><th>Détection</th><th>IC95 descriptif</th><th>Segments / paires</th><th>Verdict R14b</th></tr></thead><tbody>
        <tr><td>Hap-IBD</td><td>23/25 (92 %)</td><td>{_fr(100 * hap_low, 1)}–{_fr(100 * hap_high, 1)} %</td><td>178 / 101</td><td>Intégration soutenue</td></tr>
        <tr><td>Refined IBD</td><td>22/25 (88 %)</td><td>{_fr(100 * refined_low, 1)}–{_fr(100 * refined_high, 1)} %</td><td>122 / 44</td><td>Intégration soutenue</td></tr>
      </tbody></table>
      <p>La différence apparente entre 92 % et 88 % n'est pas interprétée comme une supériorité d'un outil. Les mêmes unités synthétiques, la dépendance éventuelle et le faible effectif interdisent cette conclusion.</p>
    </section>

    <section id="s8"><h2>8. Ressources mesurées et faisabilité</h2>
      <p>R11b comprend 18 sondes IBD : deux outils, trois densités et trois répétitions par cellule. La Figure 4 montre la moyenne et l'étendue observée du temps mural. Les points ne sont pas des performances scientifiques ; ils décrivent uniquement le coût du harnais <code>RESOURCE_ONLY</code>.</p>
      <figure role="img" aria-labelledby="fig4cap"><svg viewBox="0 0 680 245" width="100%" aria-label="Temps mural IBD par outil et densité">
        <text x="8" y="18" font-size="13" font-weight="700">Temps mural IBD : moyenne et étendue de trois répétitions</text>{_ibd_resource_svg(summaries)}
        <g font-size="9" fill="#555"><line x1="180" y1="225" x2="540" y2="225" stroke="#777"/><text x="175" y="240">0 s</text><text x="525" y="240">{_fr(max(item['wall_max_seconds'] for item in summaries), 2)} s</text></g>
      </svg><figcaption id="fig4cap">Figure 4 : Temps mural moyen (point) et étendue min–max (trait), n = 3 par cellule. Source : 18 reçus R11b.</figcaption></figure>
      <table><caption>Table 6 : Résumé des 18 mesures IBD. Temps en secondes ; RSS en Mio. Source : <code>ibd_batch_measurements.json</code>.</caption><thead><tr><th>Densité</th><th>Outil</th><th>n</th><th>Temps moyen [min–max]</th><th>RSS moyen [min–max]</th></tr></thead><tbody>{ibd_rows}</tbody></table>
      <table><caption>Table 7 : Contrôles de faisabilité R11b. Les unités diffèrent selon la ligne et sont conservées telles que définies par le contrat.</caption><thead><tr><th>Contrôle</th><th>Projeté</th><th>Plafond</th><th>Statut</th></tr></thead><tbody>{resource_rows}</tbody></table>
      <p>Le calendrier projeté de R11b utilisait 24 workers, alors que l'hôte local n'offre que 8 cœurs logiques. Cette projection a servi à démontrer une faisabilité contractuelle de release, pas à annoncer la durée réelle sur cet ordinateur. L'appel SHAPEIT5 R13 observé en { _fr(wall_hours, 2) } heures est la mesure locale la plus pertinente.</p>
    </section>

    <section id="s9"><h2>9. Disposition des checkpoints et prévention des reprises accidentelles</h2>
      <div class="stat-grid">
        <div class="stat-card"><div class="num">1 180</div><div class="label">appels du plan historique</div></div>
        <div class="stat-card"><div class="num">1</div><div class="label">appel terminé et retenu</div></div>
        <div class="stat-card"><div class="num">1 179</div><div class="label">appels retirés non exécutés</div></div>
        <div class="stat-card"><div class="num">0</div><div class="label">nouvel appel autorisé</div></div>
      </div>
      <p>Le registre JSON historique conserve honnêtement 1 <code>COMPLETED</code> et 1 179 <code>PENDING</code>. L'overlay d'amendement donne la disposition décisionnelle : les lignes restantes sont retirées, et non déclarées réussies. Le chemin de réservation vérifie l'identité de l'amendement et l'empreinte du registre avant de refuser tout nouveau lot.</p>
      <pre><code>si overlay_R13b présent :
  vérifier identité, release, run et compteurs
  vérifier SHA-256 du registre historique
  exiger authorized_call_ids = []
  refuser claim_next_batch

résultat attendu : checkpoint_execution_closed_by_amendment</code></pre>
      <table><caption>Table 8 : Disposition amendée de la campagne R13. Source : reçu R13 d'intégration.</caption><thead><tr><th>Classe</th><th>Nombre</th><th>État interprétable</th></tr></thead><tbody>
        <tr><td>Retenu</td><td>1</td><td>Terminé, vérifié et utilisé comme preuve SHAPEIT5</td></tr>
        <tr><td>Retiré</td><td>1 179</td><td>Non exécuté ; aucune preuve de performance imputée</td></tr>
        <tr><td>Autorisé après amendement</td><td>0</td><td>Verrou d'exécution actif</td></tr>
      </tbody></table>
    </section>

    <section id="s10"><h2>10. Interprétation, risques et limites</h2>
      <table><caption>Table 9 : Ce qui est démontré et ce qui ne l'est pas. Source : synthèse contractuelle R12–R14b.</caption><thead><tr><th>Objet</th><th>État</th><th>Justification</th></tr></thead><tbody>
        <tr><td>Commande SHAPEIT5 conforme</td><td>Oui</td><td>Version, graine, thread, carte, référence, pedigree et sortie audités</td></tr>
        <tr><td>Sortie SHAPEIT5 exploitable techniquement</td><td>Oui</td><td>BCF et CSI valides, dimensions exactes, hashes vérifiés</td></tr>
        <tr><td>Adaptateurs IBD fonctionnels</td><td>Oui, avec portée limitée</td><td>Sorties non vides et contrôles positifs détectés</td></tr>
        <tr><td>Exactitude de phase complète</td><td>Non évaluée</td><td>Pas de SER aveugle sur la grille R07/R08</td></tr>
        <tr><td>Sensibilité/spécificité IBD complète</td><td>Non évaluée</td><td>Grille de longueurs et négatifs indépendants incomplète</td></tr>
        <tr><td>Généralisation aux données réelles</td><td>Non évaluée</td><td>Données synthétiques uniquement</td></tr>
      </tbody></table>
      <h3>Risques résiduels</h3><ul>
        <li>Un chromosome et une graine ne couvrent pas l'hétérogénéité des 22 autosomes.</li>
        <li>Les contrôles IBD R12 ne permettent pas d'estimer un faux positif rare à l'échelle prévue.</li>
        <li>Les cartes génétiques R13 sont synthétiques et ne représentent pas une carte humaine empirique.</li>
        <li>Les intervalles de Wilson de la Figure 3 supposent une indépendance non démontrée.</li>
        <li>La faisabilité à 24 workers ne décrit pas la durée sur l'hôte local à 8 cœurs.</li>
      </ul>
      <aside class="callout"><span class="label">Condition d'une future qualification complète</span>Une nouvelle campagne devrait être décidée prospectivement, financer les unités indépendantes manquantes, conserver l'aveuglement et produire les métriques R07/R08 avec leurs intervalles préspecifiés. R14b ne l'autorise pas implicitement.</aside>
    </section>

    <section id="s11"><h2>11. Décision finale R14b</h2>
      <table><caption>Table 10 : Décision finale par objet. Source : reçu R13 amendé et audit R14b.</caption><thead><tr><th>Objet</th><th>Décision</th><th>Conséquence</th></tr></thead><tbody>
        <tr><td>Intégration SHAPEIT5</td><td><strong>VÉRIFIÉE</strong></td><td>Le raccordement peut être conservé sous les versions et hashes liés</td></tr>
        <tr><td>Intégration Hap-IBD</td><td><strong>SOUTENUE</strong></td><td>Adaptateur fonctionnel ; performance exhaustive non qualifiée</td></tr>
        <tr><td>Intégration Refined IBD</td><td><strong>SOUTENUE</strong></td><td>Adaptateur fonctionnel ; performance exhaustive non qualifiée</td></tr>
        <tr><td>Campagne R13 exhaustive</td><td><strong>RETIRÉE</strong></td><td>1 179 appels non exécutés, zéro appel nouveau</td></tr>
        <tr><td>Qualification scientifique complète</td><td><strong>NON ACCORDÉE</strong></td><td>Performances R07/R08 non évaluées</td></tr>
        <tr><td>Effet fondateur</td><td><strong>NON DÉMONTRÉ</strong></td><td>Aucune conclusion biologique issue de R14b</td></tr>
        <tr><td>Promotion en production</td><td><strong>NON AUTORISÉE</strong></td><td>Les trois drapeaux restent faux</td></tr>
      </tbody></table>
      <p>La décision R14b est une clôture positive pour l'intégration et restrictive pour toute interprétation scientifique plus large. Elle remplace le prochain calcul lourd par une preuve proportionnée à la question réellement posée.</p>
    </section>

    <section id="s12"><h2>12. Annexes, références et traçabilité</h2>
      <h3>12.1 Références méthodologiques</h3><ol>
        <li>Hofmeister RJ et al. (2023). Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. <em>Nature Genetics</em>. DOI : <a href="https://doi.org/10.1038/s41588-023-01415-w">10.1038/s41588-023-01415-w</a>.</li>
        <li>Zhou Y, Browning SR, Browning BL (2020). A Fast and Simple Method for Detecting Identity-by-Descent Segments in Large-Scale Data. <em>American Journal of Human Genetics</em>. DOI : <a href="https://doi.org/10.1016/j.ajhg.2020.02.010">10.1016/j.ajhg.2020.02.010</a>.</li>
        <li>Browning BL, Browning SR (2013). Improving the Accuracy and Efficiency of Identity-by-Descent Detection in Population Data. <em>Genetics</em>. DOI : <a href="https://doi.org/10.1534/genetics.113.150029">10.1534/genetics.113.150029</a>.</li>
      </ol>
      <h3>12.2 Glossaire</h3>
      <table><caption>Table 11 : Glossaire des termes principaux.</caption><thead><tr><th>Terme</th><th>Définition dans ce rapport</th></tr></thead><tbody>
        <tr><td>BCF</td><td>Format binaire du Variant Call Format, ici sortie phasée de SHAPEIT5</td></tr>
        <tr><td>CSI</td><td>Index coordonné du BCF</td></tr>
        <tr><td>IBD</td><td>Identité par descendance</td></tr>
        <tr><td>RSS</td><td>Mémoire résidente maximale observée</td></tr>
        <tr><td>SER</td><td>Taux d'erreur de commutation de phase, non mesuré dans R14b</td></tr>
        <tr><td>Wilson</td><td>Intervalle de proportion binomiale utilisé ici à titre descriptif</td></tr>
      </tbody></table>
      <h3>12.3 Reproductibilité et impression</h3>
      <p>Le manifeste R14b lie le HTML, le PDF, la décision et l'inventaire par SHA-256. Le HTML contient ses styles, sa logique de pagination et Paged.js ; aucune feuille de style, police ou bibliothèque graphique distante n'est requise. Pour une impression manuelle : Chrome ou Firefox, format A4, marges par défaut et arrière-plans activés.</p>
      <h3>12.4 Contrôle qualité report-formater</h3>
      <p>Le contrôle vérifie l'absence de placeholder, la continuité de quatre figures et onze tables, les liens internes, les identifiants dupliqués, le pied de page, la pagination et les ressources externes. Les liens DOI sont des références cliquables ; ils ne sont pas nécessaires au rendu.</p>
    </section>
  </main>
  <footer class="page-footer"><span>Rapport final R14b · Effet_fondateur2 · {FIXED_DATE}</span></footer>
</body>"""


def _render_html(template: str, body: str) -> str:
    start = template.find("<body")
    end = template.rfind("</body>")
    if start < 0 or end < 0 or end <= start:
        _fail("r14b_template_body_missing")
    rendered = template[:start] + body + template[end + len("</body>"):]
    title_start = rendered.find("<title>")
    title_end = rendered.find("</title>", title_start)
    if title_start < 0 or title_end < 0:
        _fail("r14b_template_title_missing")
    return (
        rendered[: title_start + len("<title>")]
        + html.escape(REPORT_TITLE)
        + rendered[title_end:]
    )


def _quality_control(html_text: str, page_count: int) -> dict[str, Any]:
    identifiers = re.findall(r'\bid="([^"]+)"', html_text)
    internal_links = re.findall(r'href="#([^"]+)"', html_text)
    missing_links = sorted(set(internal_links) - set(identifiers))
    figure_numbers = [int(value) for value in re.findall(r"Figure (\d+) :", html_text)]
    table_numbers = [int(value) for value in re.findall(r"Table (\d+) :", html_text)]
    external_resources = re.findall(r'(?:src|href)="https?://[^"]+"', html_text)
    external_render_resources = [item for item in external_resources if item.startswith("src=")]
    return {
        "page_count": page_count,
        "figure_count": len(set(figure_numbers)),
        "table_count": len(set(table_numbers)),
        "figure_numbering_continuous": sorted(set(figure_numbers)) == list(range(1, 5)),
        "table_numbering_continuous": sorted(set(table_numbers)) == list(range(1, 12)),
        "duplicate_html_ids": len(identifiers) - len(set(identifiers)),
        "missing_internal_links": missing_links,
        "placeholder_count": html_text.count("[To complete:"),
        "external_render_resources": len(external_render_resources),
        "doi_link_count": len([item for item in external_resources if "doi.org" in item]),
        "footer_present": 'class="page-footer"' in html_text,
        "pagedjs_inlined": "Paged.js v0.4.3" in html_text,
    }


def build_r14b_report(
    repository_root: Path,
    r13_run_root: Path,
    r12_ibd_root: Path,
    r11b_root: Path,
    output_root: Path,
    template_path: Path,
    chrome_path: Path,
) -> dict[str, Any]:
    """Construit transactionnellement le rapport R14b et son PDF."""
    repository_root = repository_root.resolve(strict=True)
    r13_run_root = r13_run_root.resolve(strict=True)
    r12_ibd_root = r12_ibd_root.resolve(strict=True)
    r11b_root = r11b_root.resolve(strict=True)
    template_path = template_path.resolve(strict=True)
    chrome_path = chrome_path.resolve(strict=True)
    output_root = output_root.resolve(strict=False)
    if output_root.exists():
        _fail("r14b_output_already_exists", str(output_root))
    staging = output_root.with_name(f".{output_root.name}.building")
    if staging.exists():
        _fail("r14b_staging_already_exists", str(staging))
    data, sources = _validate_sources(
        repository_root,
        r13_run_root,
        r12_ibd_root,
        r11b_root,
        template_path,
    )
    amendment = data["amendment"]
    decision = {
        "schema_version": "1.0.0",
        "decision_id": REPORT_ID,
        "scope": REPORT_SCOPE,
        "release_id": amendment["release_id"],
        "closed_run_id": amendment["run_id"],
        "status": "COMPLETED_WITH_TOOL_INTEGRATION_VERIFIED",
        "shapeit5_integration": "VERIFIED",
        "hap_ibd_integration": "SUPPORTED",
        "refined_ibd_integration": "SUPPORTED",
        "full_r07_r08_performance": "NOT_EVALUATED",
        "scientific_qualification": "NOT_GRANTED",
        "founder_effect": "NOT_PROVEN",
        "production_promotion": "NOT_AUTHORIZED",
        "truth_opened": False,
        "historical_calls": 1_180,
        "retained_completed_calls": 1,
        "withdrawn_unexecuted_calls": 1_179,
        "newly_authorized_calls": 0,
        "qualification_flags": QUALIFICATION_FLAGS,
        "next_action": "HUMAN_DECISION_REQUIRED_BEFORE_ANY_FULL_PERFORMANCE_CAMPAIGN_OR_REAL_DATA_ANALYSIS",
    }
    try:
        staging.mkdir(parents=True)
        decision_path = staging / "R14B_DECISION.json"
        inventory_path = staging / "R14B_SOURCE_INVENTORY.json"
        html_path = staging / "rapport_final_r14b.html"
        pdf_path = staging / "rapport_final_r14b.pdf"
        _write_json(decision_path, decision)
        _write_json(
            inventory_path,
            {
                "schema_version": "1.0.0",
                "report_id": REPORT_ID,
                "source_authority": [
                    "CLOSED_RAW_RESULTS",
                    "TRACKED_GOVERNANCE_RECEIPTS",
                    "NORMATIVE_PLAN_AND_CONVENTIONS",
                    "METHOD_PUBLICATIONS",
                ],
                "template_path": template_path.as_posix(),
                "template_sha256": file_sha256(template_path),
                "sources": sources,
            },
        )
        body = _report_body(data, sources)
        html_path.write_text(
            _render_html(template_path.read_text(encoding="utf-8"), body),
            encoding="utf-8",
        )
        page_count = _run_chrome(chrome_path, html_path, pdf_path)
        qc = _quality_control(html_path.read_text(encoding="utf-8"), page_count)
        if (
            not qc["figure_numbering_continuous"]
            or not qc["table_numbering_continuous"]
            or qc["duplicate_html_ids"]
            or qc["missing_internal_links"]
            or qc["placeholder_count"]
            or qc["external_render_resources"]
            or not qc["footer_present"]
            or not qc["pagedjs_inlined"]
        ):
            _fail("r14b_quality_control_failed", json.dumps(qc, sort_keys=True))
        artifacts = [decision_path, inventory_path, html_path, pdf_path]
        manifest = {
            "schema_version": "1.0.0",
            "report_id": REPORT_ID,
            "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "scope": REPORT_SCOPE,
            "release_id": amendment["release_id"],
            "closed_run_id": amendment["run_id"],
            "rendering": {
                "engine": "REPORT_FORMATER_PAGEDJS_0.4.3_CHROME_DEVTOOLS_PRINT",
                "paper": "A4",
                "page_count": page_count,
                "self_contained_html": True,
            },
            "quality_control": qc,
            "artifacts": [
                {
                    "relative_path": path.name,
                    "sha256": file_sha256(path),
                    "size_bytes": path.stat().st_size,
                }
                for path in artifacts
            ],
            "qualification_flags": QUALIFICATION_FLAGS,
        }
        _write_json(staging / "R14B_REPORT_MANIFEST.json", manifest)
        staging.rename(output_root)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return validate_r14b_report(output_root)


def validate_r14b_report(output_root: Path) -> dict[str, Any]:
    """Valide les artefacts et la portée restrictive du rapport R14b."""
    output_root = output_root.resolve(strict=True)
    manifest = load_json(output_root / "R14B_REPORT_MANIFEST.json")
    decision = load_json(output_root / "R14B_DECISION.json")
    if manifest.get("report_id") != REPORT_ID or decision.get("decision_id") != REPORT_ID:
        _fail("r14b_identity_mismatch")
    if decision.get("shapeit5_integration") != "VERIFIED":
        _fail("r14b_shapeit_decision_invalid")
    if decision.get("scientific_qualification") != "NOT_GRANTED":
        _fail("r14b_scientific_qualification_invalid")
    if decision.get("qualification_flags") != QUALIFICATION_FLAGS:
        _fail("r14b_qualification_flags_invalid")
    if decision.get("newly_authorized_calls") != 0 or decision.get("truth_opened") is not False:
        _fail("r14b_scope_violation")
    for artifact in manifest.get("artifacts", []):
        path = output_root / artifact["relative_path"]
        if not path.is_file() or file_sha256(path) != artifact.get("sha256"):
            _fail("r14b_artifact_hash_mismatch", artifact.get("relative_path"))
        if path.stat().st_size != artifact.get("size_bytes"):
            _fail("r14b_artifact_size_mismatch", artifact.get("relative_path"))
    html_text = (output_root / "rapport_final_r14b.html").read_text(encoding="utf-8")
    for marker in (
        REPORT_TITLE,
        "Figure 1 :",
        "Figure 4 :",
        "Table 1 :",
        "Table 11 :",
        "IC95 de Wilson",
        "1 179 appels non exécutés",
        "performance complète non évaluée",
        "page-footer",
        "Paged.js v0.4.3",
    ):
        if marker not in html_text:
            _fail("r14b_html_marker_missing", marker)
    # Le chemin générique data/input/ peut apparaître dans la déclaration
    # d'absence de données réelles ; ce n'est ni un fichier ni une donnée.
    for marker in ("[To complete:", "calibrated=true"):
        if marker in html_text:
            _fail("r14b_html_forbidden_marker", marker)
    pdf_path = output_root / "rapport_final_r14b.pdf"
    if not pdf_path.read_bytes().startswith(b"%PDF-") or pdf_path.stat().st_size < 100_000:
        _fail("r14b_pdf_invalid")
    qc = manifest.get("quality_control", {})
    if qc.get("figure_count") != 4 or qc.get("table_count") != 11:
        _fail("r14b_visual_inventory_invalid")
    if manifest.get("rendering", {}).get("page_count", 0) < 12:
        _fail("r14b_page_count_invalid")
    return manifest


def readiness_receipt_document(output_root: Path) -> dict[str, Any]:
    """Construit le reçu suivi à partir d'un rendu R14b validé."""
    manifest = validate_r14b_report(output_root)
    decision = load_json(output_root.resolve(strict=True) / "R14B_DECISION.json")
    return {
        "schema_version": "1.0.0",
        "receipt_id": "R14B-INTEGRATION-CLOSURE-PASS",
        "report_id": REPORT_ID,
        "report_scope": REPORT_SCOPE,
        "accepted_output": output_root.as_posix(),
        "status": decision["status"],
        "shapeit5_integration": decision["shapeit5_integration"],
        "hap_ibd_integration": decision["hap_ibd_integration"],
        "refined_ibd_integration": decision["refined_ibd_integration"],
        "full_r07_r08_performance": decision["full_r07_r08_performance"],
        "scientific_qualification": decision["scientific_qualification"],
        "founder_effect": decision["founder_effect"],
        "production_promotion": decision["production_promotion"],
        "newly_authorized_calls": decision["newly_authorized_calls"],
        "qualification_flags": QUALIFICATION_FLAGS,
        "quality_control": manifest["quality_control"],
        "artifacts": manifest["artifacts"],
        "manifest_sha256": file_sha256(
            output_root.resolve(strict=True) / "R14B_REPORT_MANIFEST.json"
        ),
        "next_action": decision["next_action"],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "validate", "write-readiness"))
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--r13-run-root", type=Path)
    parser.add_argument("--r12-ibd-root", type=Path)
    parser.add_argument("--r11b-root", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--chrome", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.action == "build":
        required = (
            arguments.r13_run_root,
            arguments.r12_ibd_root,
            arguments.r11b_root,
            arguments.template,
            arguments.chrome,
        )
        if any(value is None for value in required):
            _fail("r14b_build_arguments_missing")
        result = build_r14b_report(
            arguments.repository_root,
            arguments.r13_run_root,
            arguments.r12_ibd_root,
            arguments.r11b_root,
            arguments.output_root,
            arguments.template,
            arguments.chrome,
        )
    elif arguments.action == "validate":
        result = validate_r14b_report(arguments.output_root)
    else:
        result = readiness_receipt_document(arguments.output_root)
        _write_json(arguments.repository_root / READINESS_PATH, result)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
