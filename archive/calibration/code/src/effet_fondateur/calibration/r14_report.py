"""Clôture R14 et rapport final du périmètre technique effectivement validé.

Le rapport ne recalcule aucune statistique. Il vérifie les reçus R11–R13,
reproduit leurs valeurs agrégées et rend une décision restrictive : intégrité
technique validée, qualification scientifique non accordée.
"""

from __future__ import annotations

import argparse
import asyncio
import base64
from datetime import datetime, timezone
from hashlib import sha256
import html
import json
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from typing import Any, Mapping
import urllib.request

import websockets

from effet_fondateur.calibration.contracts import file_sha256, load_json
from effet_fondateur.calibration.r13_evaluator import validate_evaluations
from effet_fondateur.calibration.r13_runner import DEFAULT_RUN_ID, validate_closed_run


REPORT_ID = "R14-TECHNICAL-CLOSURE-001"
REPORT_TITLE = "Rapport final R14 — validation technique et décision de non-qualification scientifique"
FIXED_DATE = "13 novembre 2026"


class R14ReportError(RuntimeError):
    """Erreur bloquante de construction ou d'audit du rapport R14."""

    def __init__(self, code: str, detail: str | None = None) -> None:
        self.code = code
        self.detail = detail
        super().__init__(code if detail is None else f"{code}: {detail}")


def _fail(code: str, detail: str | None = None) -> None:
    raise R14ReportError(code, detail)


def _write_json(path: Path, document: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _source_record(root: Path, path: Path, role: str) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    try:
        relative = resolved.relative_to(root).as_posix()
    except ValueError:
        relative = resolved.as_posix()
    return {
        "relative_path": relative,
        "role": role,
        "sha256": file_sha256(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def _validate_sources(
    repository_root: Path,
    candidate_root: Path,
    r13_workspace: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    run_manifest = validate_closed_run(
        repository_root,
        candidate_root,
        r13_workspace,
        DEFAULT_RUN_ID,
    )
    decision = validate_evaluations(repository_root, r13_workspace, DEFAULT_RUN_ID)
    receipt_path = repository_root / "calibration/readiness/R13_TECHNICAL_VALIDATION.json"
    receipt = load_json(receipt_path)
    projection_path = (
        repository_root
        / "calibration_workspace/resource_probes/"
        "r11_shapeit_sample_scaling_v1/projection.json"
    )
    projection = load_json(projection_path)
    if decision.get("status") != "PASS" or decision.get("stage_count") != 25:
        _fail("r13_technical_decision_not_pass")
    if decision.get("scientific_performance_evaluated") is not False:
        _fail("r13_scientific_scope_violation")
    if receipt.get("gate_g6_scientific_passed") is not False:
        _fail("scientific_gate_must_be_closed")
    if receipt.get("manifest_sha256") != run_manifest.get("manifest_sha256"):
        _fail("r13_receipt_manifest_mismatch")
    if projection.get("decision", {}).get("status") != "RESOURCE_BLOCKED":
        _fail("r11_resource_decision_mismatch")
    sources = [
        _source_record(
            repository_root,
            repository_root / "calibration/PLAN_DEVELOPPEMENT_PANEL_SIMULATION.md",
            "plan normatif R00–R14",
        ),
        _source_record(
            repository_root,
            repository_root / "calibration/CONVENTIONS_DEVELOPPEMENT.md",
            "contrat permanent des rapports",
        ),
        _source_record(repository_root, repository_root / "SESSION.md", "état de session"),
        _source_record(
            repository_root,
            repository_root / "calibration/readiness/R12_FROZEN_CANDIDATE.json",
            "rejet du candidat R12 initial",
        ),
        _source_record(
            repository_root,
            repository_root / "calibration/readiness/R12B_G5_READY.json",
            "autorisation G5 technique",
        ),
        _source_record(repository_root, receipt_path, "reçu R13 clos"),
        _source_record(repository_root, projection_path, "projection de ressources R11"),
        _source_record(
            repository_root,
            r13_workspace / f"runs/{DEFAULT_RUN_ID}/pipeline_run_manifest.json",
            "manifeste du run aveugle clos",
        ),
        _source_record(
            repository_root,
            r13_workspace / "evaluations/R13_TECHNICAL_DECISION.json",
            "décision technique R13",
        ),
        _source_record(
            repository_root,
            candidate_root / "governance/r12b_frozen_manifest.json",
            "gel R12b",
        ),
    ]
    return run_manifest, decision, projection, sources


def _source_table_rows(sources: list[dict[str, Any]]) -> str:
    rows = []
    for index, item in enumerate(sources, start=1):
        path = Path(item["relative_path"])
        display_path = "/".join(path.parts[-2:]) if len(path.parts) > 1 else path.name
        rows.append(
            "<tr>"
            f"<td>{index}</td>"
            f"<td><code>{html.escape(display_path)}</code></td>"
            f"<td>{html.escape(item['role'])}</td>"
            f"<td><code>{item['sha256'][:12]}…</code></td>"
            "</tr>"
        )
    return "".join(rows)


def _report_body(
    receipt: Mapping[str, Any],
    run_manifest: Mapping[str, Any],
    decision: Mapping[str, Any],
    projection: Mapping[str, Any],
    sources: list[dict[str, Any]],
) -> str:
    wall = projection["fits"]["wall_seconds"]["safe_projection_hours"]
    wall_ceiling = projection["budget_checks"][0]["ceiling"]
    cpu = projection["budget_checks"][1]["projected"]
    cpu_ceiling = projection["budget_checks"][1]["ceiling"]
    rss = projection["fits"]["maximum_rss_bytes"]["safe_projection_gibibytes"]
    run_hash = run_manifest["manifest_sha256"]
    source_rows = _source_table_rows(sources)
    return f"""
<body lang="fr">
  <section class="cover" id="cover">
    <div class="cover-logo"><span style="font-family:Georgia,serif;font-size:1.45em;font-weight:700;color:var(--accent);letter-spacing:.08em;">EFFET FONDATEUR 2</span></div>
    <div class="cover-accent-band"></div>
    <div class="cover-center">
      <h1>{REPORT_TITLE}</h1>
      <div class="cover-subtitle">Panel synthétique générique · clôture du périmètre technique R00–R14</div>
    </div>
    <div class="cover-meta">
      <table>
        <tr><td>Décision</td><td>Intégrité technique validée ; qualification scientifique non accordée</td></tr>
        <tr><td>Release</td><td><code>{receipt['release_id']}</code></td></tr>
        <tr><td>Date</td><td>{FIXED_DATE}</td></tr>
        <tr><td>Version</td><td>1.0.0 · {REPORT_ID}</td></tr>
      </table>
    </div>
  </section>

  <nav class="toc-page" id="toc">
    <h2>Sommaire</h2>
    <ol>
      <li><a class="toc-entry" href="#s1"><span class="toc-title">1. Résumé exécutif</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s2"><span class="toc-title">2. Contexte et question de décision</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s3"><span class="toc-title">3. Objectifs de R14</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s4"><span class="toc-title">4. Architecture et séparation</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s5"><span class="toc-title">5. Méthodologie et règles</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s6"><span class="toc-title">6. Résultats R11–R13</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s7"><span class="toc-title">7. Ressources produites</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s8"><span class="toc-title">8. Ressources de calcul</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s9"><span class="toc-title">9. Portée et impacts</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s10"><span class="toc-title">10. Risques et limites</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s11"><span class="toc-title">11. Décision R14</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
      <li><a class="toc-entry" href="#s12"><span class="toc-title">12. Annexes et traçabilité</span><span class="toc-dots"></span><span class="toc-page-num"></span></a></li>
    </ol>
  </nav>

  <main>
    <section id="s1">
      <div class="summary">
        <h2>1. Résumé exécutif</h2>
        <ul>
          <li><strong>R13 technique est clos avec 25/25 étapes PASS.</strong> Le runner n'a ouvert aucune vérité ; l'évaluateur l'a ouverte seulement après validation du run <code>CLOSED</code>.</li>
          <li><strong>Aucune performance scientifique n'a été évaluée.</strong> Aucun outil scientifique externe ni donnée réelle n'a été utilisé dans R12b–R13.</li>
          <li><strong>La qualification scientifique n'est pas accordée.</strong> R11 reste <code>RESOURCE_BLOCKED</code>, R07 reste <code>REVISION_REQUIRED</code> et G6 scientifique est fermée.</li>
          <li><strong>Aucun effet fondateur n'est démontré</strong> et aucune promotion en production n'est autorisée.</li>
          <li>Le présent document est le rendu HTML/PDF unique autorisé à R14 ; il reprend les sources closes sans recalcul statistique.</li>
        </ul>
      </div>
      <div class="stat-grid">
        <div class="stat-card"><div class="num">25 / 25</div><div class="label">étapes techniques PASS</div></div>
        <div class="stat-card"><div class="num">26</div><div class="label">sorties R13 empreintées</div></div>
        <div class="stat-card"><div class="num">0</div><div class="label">outil scientifique exécuté</div></div>
        <div class="stat-card"><div class="num">G6 fermé</div><div class="label">qualification scientifique</div></div>
      </div>
      <aside class="callout"><span class="label">Décision centrale</span>Le succès R13 prouve la traversabilité et l'isolement du dispositif technique. Il ne mesure ni sensibilité, ni spécificité, ni exactitude de phase, ni performance IBD, ni datation.</aside>
    </section>

    <section id="s2" style="break-before:page">
      <h2>2. Contexte et question de décision</h2>
      <p>Le dépôt met en œuvre un pipeline générique pour étudier un possible effet fondateur autour d'une variation autosomique configurée. Le plan unifié R00–R14 sépare la conception, le micro-panel technique, le pilote de calibration, le gel, la validation aveugle et la décision finale.</p>
      <p>Le pilote R11 a fourni des mesures réelles de ressources sur des données entièrement synthétiques. La projection de qualification dépasse deux plafonds R07. Le candidat R12 initial a en outre été rejeté faute de partition <code>VALIDATION_BLIND</code> et de composants R13 gelés. R12b a corrigé uniquement cette lacune technique, sans prétendre résoudre le dimensionnement scientifique.</p>
      <table>
        <caption>Table 1 : Question R14 et réponse fondée sur les portes closes. Source : plan R00–R14, reçus R12–R13.</caption>
        <thead><tr><th>Question</th><th>Preuve disponible</th><th>Décision</th></tr></thead>
        <tbody>
          <tr><td>Le chemin technique aveugle est-il fermé et vérifiable ?</td><td>R13 : 25/25 PASS, hashes et clôture vérifiés</td><td><strong>Oui</strong></td></tr>
          <tr><td>Les performances scientifiques sont-elles qualifiées ?</td><td>R11 bloqué ; R13 limité à l'intégrité technique</td><td><strong>Non</strong></td></tr>
          <tr><td>Un effet fondateur est-il démontré ?</td><td>Aucune métrique scientifique R13</td><td><strong>Non</strong></td></tr>
          <tr><td>La promotion en production est-elle autorisée ?</td><td>G6 scientifique fermée</td><td><strong>Non</strong></td></tr>
        </tbody>
      </table>
    </section>

    <section id="s3" style="break-before:page">
      <h2>3. Objectifs de R14</h2>
      <ol>
        <li>Assembler les preuves closes R11–R13 sans recalcul.</li>
        <li>Rendre visibles les effectifs, règles de décision, incertitudes et limites réellement disponibles.</li>
        <li>Décider séparément l'intégrité technique et la qualification scientifique.</li>
        <li>Fournir un rapport autonome, traçable et imprimable en PDF.</li>
      </ol>
      <aside class="callout"><span class="label">Hors périmètre</span>R14 ne complète pas les cellules scientifiques absentes, ne modifie pas R07/R08 et ne transforme aucun <code>NOT_EVALUABLE</code> ou blocage de ressources en succès.</aside>
    </section>

    <section id="s4" style="break-before:page">
      <h2>4. Architecture et séparation calibration–validation</h2>
      <figure role="img" aria-labelledby="fig1cap">
        <svg viewBox="0 0 680 145" width="100%" aria-label="Flux aveugle R12b à R14">
          <defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#1f4e79"/></marker></defs>
          <text x="8" y="18" font-size="13" font-weight="700">Séquence contrôlée et séparation des privilèges</text>
          <g font-size="11" text-anchor="middle">
            <rect x="8" y="42" width="116" height="58" rx="6" fill="#eaf0f6" stroke="#1f4e79"/><text x="66" y="64"><tspan x="66">R12b</tspan><tspan x="66" dy="16">vérité scellée</tspan></text>
            <rect x="145" y="42" width="116" height="58" rx="6" fill="#eaf0f6" stroke="#1f4e79"/><text x="203" y="64"><tspan x="203">Runner R13</tspan><tspan x="203" dy="16">vues aveugles</tspan></text>
            <rect x="282" y="42" width="116" height="58" rx="6" fill="#eaf0f6" stroke="#1f4e79"/><text x="340" y="64"><tspan x="340">CLOSE</tspan><tspan x="340" dy="16">26 sorties hashées</tspan></text>
            <rect x="419" y="42" width="116" height="58" rx="6" fill="#eaf0f6" stroke="#1f4e79"/><text x="477" y="64"><tspan x="477">Évaluateur</tspan><tspan x="477" dy="16">vérité minimale</tspan></text>
            <rect x="556" y="42" width="116" height="58" rx="6" fill="#dcefe2" stroke="#29633a"/><text x="614" y="64"><tspan x="614">R14</tspan><tspan x="614" dy="16">décision bornée</tspan></text>
          </g>
          <g stroke="#1f4e79" stroke-width="2" marker-end="url(#arrow)"><line x1="124" y1="71" x2="141" y2="71"/><line x1="261" y1="71" x2="278" y2="71"/><line x1="398" y1="71" x2="415" y2="71"/><line x1="535" y1="71" x2="552" y2="71"/></g>
          <text x="340" y="126" text-anchor="middle" font-size="11" fill="#555">La vérité n'est accessible qu'après CLOSE ; aucun chemin inverse n'est autorisé.</text>
        </svg>
        <figcaption id="fig1cap">Figure 1 : Flux d'aveuglement et de clôture. Source : manifeste R12b, manifeste du run R13 et décision R13.</figcaption>
      </figure>
      <table>
        <caption>Table 2 : État des phases déterminantes. Source : SESSION.md et reçus suivis.</caption>
        <thead><tr><th>Phase</th><th>État</th><th>Portée</th></tr></thead>
        <tbody>
          <tr><td>R10</td><td>Réalisé</td><td>Micro-panel d'intégrité technique</td></tr>
          <tr><td>R11</td><td><code>RESOURCE_BLOCKED</code></td><td>Pilote et projection de ressources</td></tr>
          <tr><td>R12</td><td>Rejeté</td><td>Candidat scientifique incomplet</td></tr>
          <tr><td>R12b / G5</td><td>Prêt</td><td>Validation aveugle technique seulement</td></tr>
          <tr><td>R13</td><td>25/25 <code>PASS</code></td><td>Intégrité technique seulement</td></tr>
          <tr><td>R14</td><td>Clôture restrictive</td><td>Non-qualification scientifique</td></tr>
        </tbody>
      </table>
    </section>

    <section id="s5" style="break-before:page">
      <h2>5. Méthodologie, équations et règles de décision</h2>
      <h3>5.1 Modèle de projection de ressources réellement utilisé</h3>
      <pre><code>log(y) = a + b × log(N)
y_sûr = U95_prédiction × r_allocation × 2

y : temps (s), RSS (octets) ou CPU (s)
N : sujets d'étude + références
a : intercept ; b : exposant estimé
r_allocation : ratio 4:3 mesuré vs allocation 1:1
2 : facteur de sécurité R07</code></pre>
      <p>Le modèle est une régression OLS log–log sur cinq effectifs égaux, avec borne prédictive unilatérale à 95 %, correction d'allocation et facteur de sécurité 2. La cible de 7 000 sujets est sept fois supérieure au plus grand total mesuré ; aucun crédit de parallélisation non mesuré n'est appliqué.</p>
      <h3>5.2 Équations de dimensionnement préspecifiées mais non évaluées par R13</h3>
      <pre><code>u₀(n) = 1 − 0,05^(1/n)
p_MC = (b + 1) / (N_DRAW + 1)

n : unités négatives indépendantes
b : tirages nuls au moins aussi extrêmes que l'observé
N_DRAW : nombre total de tirages nuls</code></pre>
      <p>Avec 3 000 unités négatives et aucun événement, le plan donne <code>u₀ = 0,000998</code>. Ces équations justifient R07/R08 ; elles ne produisent aucun résultat scientifique dans R13 technique.</p>
      <h3>5.3 Règle logique R13</h3>
      <pre><code>PASS_technique =
  run_CLOSED
  ∧ inventaire_aveugle_vérifié
  ∧ sorties_déclarées_et_hashées
  ∧ vérité_ouverte_après_clôture
  ∧ 25 évaluations techniques PASS</code></pre>
      <table>
        <caption>Table 3 : Unités, incertitudes et limites des mesures rapportées. Source : projection R11 et reçu R13.</caption>
        <thead><tr><th>Mesure</th><th>Unité indépendante</th><th>Incertitude</th><th>Usage</th></tr></thead>
        <tbody>
          <tr><td>Projection SHAPEIT5</td><td>5 effectifs synthétiques mesurés</td><td>Borne prédictive unilatérale 95 %</td><td>Décision de ressources</td></tr>
          <tr><td>Étapes R13</td><td>25 contrats d'étape</td><td>Contrôle déterministe, pas d'intervalle</td><td>Intégrité technique</td></tr>
          <tr><td>Sorties R13</td><td>26 artefacts déclarés</td><td>SHA-256 exact</td><td>Fermeture et non-altération</td></tr>
          <tr><td>Performance scientifique</td><td>NA</td><td>NA — non mesurée</td><td>Qualification interdite</td></tr>
        </tbody>
      </table>
    </section>

    <section id="s6" style="break-before:page">
      <h2>6. Validation et résultats R11–R13</h2>
      <h3>6.1 Résultat technique R13</h3>
      <table>
        <caption>Table 4 : Statut des 25 étapes R13. Source : décision technique R13 close.</caption>
        <thead><tr><th>Groupe 1</th><th>Groupe 2</th><th>Groupe 3</th><th>Groupe 4</th><th>Groupe 5</th></tr></thead>
        <tbody>
          <tr><td>00 — PASS</td><td>05 — PASS</td><td>10 — PASS</td><td>14 — PASS</td><td>16B — PASS</td></tr>
          <tr><td>01 — PASS</td><td>06 — PASS</td><td>10A — PASS</td><td>15 — PASS</td><td>16C — PASS</td></tr>
          <tr><td>02 — PASS</td><td>07 — PASS</td><td>11 — PASS</td><td>16 — PASS</td><td>17 — PASS</td></tr>
          <tr><td>03 — PASS</td><td>08 — PASS</td><td>12 — PASS</td><td>16A — PASS</td><td>18 — PASS</td></tr>
          <tr><td>04 — PASS</td><td>09 — PASS</td><td>12A — PASS</td><td>13 — PASS</td><td>19 — PASS</td></tr>
        </tbody>
      </table>
      <p>Le manifeste canonique du run porte l'empreinte <code>{run_hash}</code>. Le runner déclare <code>truth_opened=false</code>, <code>scientific_tools_executed=false</code> et <code>real_data_accessed=false</code>. L'évaluateur a ensuite routé les quatorze familles de vérité selon les consommateurs autorisés.</p>
      <h3>6.2 Portes et verdicts</h3>
      <figure role="img" aria-labelledby="fig2cap">
        <svg viewBox="0 0 680 255" width="100%" aria-label="État des portes R10 à R14">
          <text x="8" y="18" font-size="13" font-weight="700">État consolidé : succès technique, blocage scientifique maintenu</text>
          <g font-size="11">
            <text x="20" y="52">R10</text><rect x="90" y="35" width="500" height="24" rx="4" fill="#dcefe2"/><text x="600" y="52">réalisé</text>
            <text x="20" y="88">R11</text><rect x="90" y="71" width="500" height="24" rx="4" fill="#f7dfd8"/><text x="600" y="88">bloqué</text>
            <text x="20" y="124">R12</text><rect x="90" y="107" width="500" height="24" rx="4" fill="#f7dfd8"/><text x="600" y="124">rejeté</text>
            <text x="20" y="160">R12b</text><rect x="90" y="143" width="500" height="24" rx="4" fill="#dcefe2"/><text x="600" y="160">G5 tech.</text>
            <text x="20" y="196">R13</text><rect x="90" y="179" width="500" height="24" rx="4" fill="#dcefe2"/><text x="600" y="196">25 PASS</text>
            <text x="20" y="232">R14</text><rect x="90" y="215" width="500" height="24" rx="4" fill="#fff1c9"/><text x="600" y="232">non qualifié</text>
          </g>
        </svg>
        <figcaption id="fig2cap">Figure 2 : État des phases R10–R14. Vert : succès dans la portée annoncée ; rouge : blocage ; ambre : clôture restrictive. Source : SESSION.md et reçus R12–R13.</figcaption>
      </figure>
    </section>

    <section id="s7" style="break-before:page">
      <h2>7. Ressources produites et provenance</h2>
      <ul>
        <li>Release R12b non écrasable avec 11 vues <code>VALIDATION_BLIND</code> et 14 familles de vérité privées.</li>
        <li>Run clos <code>{DEFAULT_RUN_ID}</code> avec 26 sorties déclarées.</li>
        <li>25 manifestes d'évaluation et 25 tables de métriques techniques.</li>
        <li>Reçus G5 technique, R13 technique et décision R14 restrictive.</li>
        <li>Rapport HTML autonome et export PDF, sans donnée génétique individuelle.</li>
      </ul>
      <table>
        <caption>Table 5 : Preuves quantitatives R13. Source : reçu R13.</caption>
        <thead><tr><th>Élément</th><th>Valeur</th><th>Portée</th></tr></thead>
        <tbody>
          <tr><td>Étapes closes</td><td>25</td><td>Technique</td></tr>
          <tr><td>Sorties déclarées</td><td>26</td><td>Technique</td></tr>
          <tr><td>Évaluations / métriques</td><td>25 / 25</td><td>Technique</td></tr>
          <tr><td>Fichiers du workspace R13</td><td>266</td><td>Audit</td></tr>
          <tr><td>Volume logique</td><td>436 910 octets</td><td>Audit</td></tr>
          <tr><td>HTML/PDF avant R14</td><td>0</td><td>Conformité</td></tr>
        </tbody>
      </table>
    </section>

    <section id="s8" style="break-before:page">
      <h2>8. Ressources de calcul et budget</h2>
      <p>Aucun budget financier n'est défini dans les sources. La décision porte sur les plafonds de calcul R07. La projection sûre pour 4 000 sujets d'étude et 3 000 références est fondée sur un chromosome médian synthétique de 26 000 variants communs, 22 autosomes et 50 graines.</p>
      <figure role="img" aria-labelledby="fig3cap">
        <svg viewBox="0 0 680 250" width="100%" aria-label="Comparaison projections et plafonds R07">
          <text x="8" y="18" font-size="13" font-weight="700">Deux plafonds dépassés ; RSS dans la limite</text>
          <g font-size="11">
            <text x="20" y="53">Temps par appel</text><rect x="155" y="37" width="430" height="22" fill="#f7dfd8"/><rect x="155" y="37" width="244" height="22" fill="#1f4e79"/><text x="595" y="53">{wall:.2f} h / {wall_ceiling} h</text>
            <text x="20" y="103">CPU agrégé</text><rect x="155" y="87" width="430" height="22" fill="#f7dfd8"/><rect x="155" y="87" width="385" height="22" fill="#1f4e79"/><text x="595" y="103">{cpu:,.0f} / {cpu_ceiling:,} h-cœur</text>
            <text x="20" y="153">RSS par profil</text><rect x="155" y="137" width="430" height="22" fill="#dcefe2"/><rect x="155" y="137" width="15" height="22" fill="#1f4e79"/><text x="595" y="153">{rss:.2f} / 192 Gio</text>
          </g>
          <g font-size="10" fill="#555"><rect x="155" y="190" width="18" height="12" fill="#1f4e79"/><text x="180" y="201">plafond R07 dans l'échelle du projeté</text><rect x="410" y="190" width="18" height="12" fill="#f7dfd8"/><text x="435" y="201">dépassement</text></g>
          <text x="20" y="228" font-size="10" fill="#555">Les longueurs de barres sont normalisées séparément par ligne ; les valeurs imprimées font foi.</text>
        </svg>
        <figcaption id="fig3cap">Figure 3 : Projections de ressources et plafonds R07. Source : <code>projection.json</code>, modèle log–log avec borne prédictive 95 % et facteur de sécurité 2.</figcaption>
      </figure>
      <table>
        <caption>Table 6 : Décisions de ressources R11. Source : projection R11 close.</caption>
        <thead><tr><th>Contrôle</th><th>Projeté</th><th>Plafond</th><th>Statut</th></tr></thead>
        <tbody>
          <tr><td>Temps, un appel autosomal</td><td>{wall:.2f} h</td><td>{wall_ceiling} h</td><td><code>RESOURCE_BLOCKED</code></td></tr>
          <tr><td>CPU, 22 autosomes × 50 graines</td><td>{cpu:,.2f} h-cœur</td><td>{cpu_ceiling:,} h-cœur</td><td><code>RESOURCE_BLOCKED</code></td></tr>
          <tr><td>RSS, profil mesuré</td><td>{rss:.2f} Gio</td><td>192 Gio</td><td>Dans la limite</td></tr>
        </tbody>
      </table>
    </section>

    <section id="s9" style="break-before:page">
      <h2>9. Portée et impacts attendus</h2>
      <p>Le résultat réduit le risque de fuite de vérité et démontre qu'un run technique peut être fermé, empreinté puis évalué avec privilège minimal. Il fournit un socle reproductible pour une future validation scientifique correctement dimensionnée.</p>
      <figure role="img" aria-labelledby="fig4cap">
        <div class="stat-grid">
          <div class="stat-card"><div class="num">PASS</div><div class="label">contrats et aveuglement</div></div>
          <div class="stat-card"><div class="num">NA</div><div class="label">performance scientifique</div></div>
          <div class="stat-card"><div class="num">NON</div><div class="label">effet fondateur démontré</div></div>
          <div class="stat-card"><div class="num">NON</div><div class="label">promotion production</div></div>
        </div>
        <figcaption id="fig4cap">Figure 4 : Matrice synthétique de décision R14. Source : reçu R13 et règles G5–G7 du plan.</figcaption>
      </figure>
      <aside class="callout"><span class="label">Impact immédiat</span>Le pipeline de production n'est pas modifié. Le bénéfice porte sur la gouvernance de calibration, l'isolation des partitions et la traçabilité.</aside>
    </section>

    <section id="s10" style="break-before:page">
      <h2>10. Risques, limites et conditions de reprise</h2>
      <ul>
        <li>La cible de projection compte 7 000 sujets, soit sept fois le plus grand total mesuré.</li>
        <li>La mise à l'échelle des threads n'est pas mesurée ; aucun gain supposé ne peut être crédité.</li>
        <li>Les 134 000,44 heures-cœur n'incluent que le phasage commun ; les coûts IBD et les autres étapes ne sont pas inclus.</li>
        <li>R13 emploie un worker technique et n'exécute pas PLINK, KING, bcftools, SHAPEIT5, Hap-IBD ou Refined IBD.</li>
        <li>Les 25 PASS R13 ne doivent jamais être présentés comme 25 validations scientifiques.</li>
      </ul>
      <h3>Conditions minimales avant une qualification complète</h3>
      <ol>
        <li>Réviser prospectivement R07 sans réduire silencieusement les effectifs ou densités.</li>
        <li>Mesurer une stratégie de parallélisation et projeter les mondes IBD.</li>
        <li>Reprendre les étapes scientifiques affectées et obtenir une validation aveugle scientifique indépendante.</li>
        <li>Produire un nouveau reçu G6 avant toute décision de promotion.</li>
      </ol>
    </section>

    <section id="s11" style="break-before:page">
      <h2>11. Décision R14</h2>
      <table>
        <caption>Table 7 : Décision finale par objet. Source : synthèse contractuelle R11–R13.</caption>
        <thead><tr><th>Objet</th><th>Décision</th><th>Justification</th></tr></thead>
        <tbody>
          <tr><td>Intégrité technique R12b–R13</td><td><strong>VALIDÉE</strong></td><td>25/25 étapes PASS, clôture et hashes vérifiés</td></tr>
          <tr><td>Qualification scientifique</td><td><strong>NON ACCORDÉE</strong></td><td>R11 bloqué, R07 à réviser, G6 fermée</td></tr>
          <tr><td>Effet fondateur</td><td><strong>NON DÉMONTRÉ</strong></td><td>Aucune performance scientifique R13</td></tr>
          <tr><td>Promotion en production</td><td><strong>REFUSÉE À CE STADE</strong></td><td>Qualification scientifique absente</td></tr>
          <tr><td>Rapport R14</td><td><strong>PUBLIABLE AVEC RESTRICTIONS</strong></td><td>HTML/PDF technique, sources et limites visibles</td></tr>
        </tbody>
      </table>
      <aside class="callout"><span class="label">Verdict</span><strong>R14 est clos pour le périmètre technique.</strong> Les drapeaux restent <code>calibrated=false</code>, <code>founder_effect_proven=false</code> et <code>production_promotion=false</code>. La prochaine action scientifique est la révision R07.</aside>
    </section>

    <section id="s12" style="break-before:page">
      <h2>12. Annexes et traçabilité</h2>
      <h3>12.1 Inventaire des sources</h3>
      <table>
        <caption>Table 8 : Sources autoritatives utilisées par R14. Les empreintes complètes figurent dans le manifeste R14.</caption>
        <thead><tr><th>#</th><th>Chemin</th><th>Rôle</th><th>SHA-256</th></tr></thead>
        <tbody>{source_rows}</tbody>
      </table>
      <h3>12.2 Drapeaux et compatibilité</h3>
      <pre><code>calibrated = false
founder_effect_proven = false
production_promotion = false
r07_status = REVISION_REQUIRED
gate_g6_scientific_passed = false</code></pre>
      <p>Toute modification future de carte, panel, densité, outil, seuil, algorithme, format, contrat scientifique ou domaine d'usage impose une analyse d'invalidation et, selon son impact, un nouveau gel et une nouvelle validation.</p>
      <h3>12.3 Mode d'impression</h3>
      <p>Ouvrir le fichier HTML dans Chrome ou Firefox, utiliser <strong>Imprimer → Enregistrer au format PDF</strong>, conserver le format A4 et les marges par défaut, et activer les arrière-plans. Le PDF livré avec ce rapport a été généré depuis ce même HTML autonome.</p>
    </section>
  </main>

  <footer class="page-footer"><span>Rapport final R14 · Effet_fondateur2 · {FIXED_DATE}</span></footer>
</body>"""


def _render_html(template: str, body: str) -> str:
    start = template.find("<body")
    end = template.rfind("</body>")
    if start < 0 or end < 0 or end <= start:
        _fail("report_template_body_missing")
    rendered = template[:start] + body + template[end + len("</body>"):]
    title_start = rendered.find("<title>")
    title_end = rendered.find("</title>", title_start)
    if title_start < 0 or title_end < 0:
        _fail("report_template_title_missing")
    return (
        rendered[: title_start + len("<title>")]
        + html.escape(REPORT_TITLE)
        + rendered[title_end:]
    )


def _free_local_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as listener:
        listener.bind(("127.0.0.1", 0))
        return int(listener.getsockname()[1])


async def _print_after_pagination(websocket_url: str, pdf_path: Path) -> int:
    async with websockets.connect(websocket_url, max_size=20_000_000) as connection:
        next_identifier = 0

        async def call(method: str, parameters: Mapping[str, Any] | None = None) -> dict[str, Any]:
            nonlocal next_identifier
            next_identifier += 1
            identifier = next_identifier
            await connection.send(
                json.dumps(
                    {"id": identifier, "method": method, "params": dict(parameters or {})}
                )
            )
            while True:
                message = json.loads(await asyncio.wait_for(connection.recv(), timeout=15))
                if message.get("id") != identifier:
                    continue
                if "error" in message:
                    _fail("chrome_devtools_error", str(message["error"]))
                return dict(message.get("result", {}))

        await call("Page.enable")
        await call("Runtime.enable")
        previous = -1
        stable = 0
        page_count = 0
        for _ in range(120):
            result = await call(
                "Runtime.evaluate",
                {
                    "expression": "document.querySelectorAll('.pagedjs_page').length",
                    "returnByValue": True,
                },
            )
            page_count = int(result["result"].get("value", 0))
            stable = stable + 1 if page_count >= 12 and page_count == previous else 0
            previous = page_count
            if stable >= 4:
                break
            await asyncio.sleep(0.15)
        if page_count < 12 or stable < 4:
            _fail("pagedjs_pagination_incomplete", str(page_count))
        result = await call(
            "Page.printToPDF",
            {
                "printBackground": True,
                "preferCSSPageSize": True,
                "displayHeaderFooter": False,
                "transferMode": "ReturnAsBase64",
            },
        )
        pdf_path.write_bytes(base64.b64decode(result["data"]))
        return page_count


def _run_chrome(chrome_path: Path, html_path: Path, pdf_path: Path) -> int:
    profile_root = Path(tempfile.mkdtemp(prefix="r14-chrome-"))
    port = _free_local_port()
    process = subprocess.Popen(
        [
            str(chrome_path),
            "--headless=new",
            "--disable-gpu",
            f"--remote-debugging-port={port}",
            "--remote-allow-origins=*",
            f"--user-data-dir={profile_root}",
            html_path.as_uri(),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        target: Mapping[str, Any] | None = None
        for _ in range(100):
            try:
                with urllib.request.urlopen(
                    f"http://127.0.0.1:{port}/json",
                    timeout=1,
                ) as response:
                    targets = json.load(response)
                target = next(
                    (
                        item
                        for item in targets
                        if item.get("type") == "page"
                        and str(item.get("url", "")).startswith("file:")
                    ),
                    None,
                )
                if target is not None:
                    break
            except (OSError, ValueError):
                pass
            time.sleep(0.1)
        if target is None:
            _fail("chrome_devtools_target_missing")
        page_count = asyncio.run(
            _print_after_pagination(str(target["webSocketDebuggerUrl"]), pdf_path)
        )
    finally:
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
        shutil.rmtree(profile_root, ignore_errors=True)
    if not pdf_path.is_file() or not pdf_path.read_bytes().startswith(b"%PDF-"):
        _fail("pdf_output_invalid")
    return page_count


def build_r14_report(
    repository_root: Path,
    candidate_root: Path,
    r13_workspace: Path,
    output_root: Path,
    template_path: Path,
    chrome_path: Path,
) -> dict[str, Any]:
    """Construit transactionnellement la décision, le HTML et le PDF R14."""
    repository_root = repository_root.resolve(strict=True)
    candidate_root = candidate_root.resolve(strict=True)
    r13_workspace = r13_workspace.resolve(strict=True)
    template_path = template_path.resolve(strict=True)
    chrome_path = chrome_path.resolve(strict=True)
    output_root = output_root.resolve(strict=False)
    if output_root.exists():
        _fail("r14_output_already_exists", str(output_root))
    staging = output_root.with_name(f".{output_root.name}.building")
    if staging.exists():
        _fail("r14_staging_already_exists", str(staging))
    receipt = load_json(
        repository_root / "calibration/readiness/R13_TECHNICAL_VALIDATION.json"
    )
    run_manifest, decision, projection, sources = _validate_sources(
        repository_root,
        candidate_root,
        r13_workspace,
    )
    r14_decision = {
        "schema_version": "1.0.0",
        "decision_id": REPORT_ID,
        "release_id": receipt["release_id"],
        "closed_run_id": receipt["run_id"],
        "scope": "R14_TECHNICAL_CLOSURE_ONLY",
        "status": "COMPLETED_WITH_SCIENTIFIC_QUALIFICATION_NOT_GRANTED",
        "technical_integrity": "PASS",
        "scientific_qualification": "NOT_GRANTED",
        "founder_effect": "NOT_PROVEN",
        "production_promotion": "NOT_AUTHORIZED",
        "r07_status": "REVISION_REQUIRED",
        "r11_status": "RESOURCE_BLOCKED",
        "gate_g6_scientific_passed": False,
        "gate_g7_promotion_passed": False,
        "scientific_performance_evaluated": False,
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
        "next_required_action": "FORMAL_R07_REVISION_BEFORE_SCIENTIFIC_QUALIFICATION",
    }
    try:
        staging.mkdir(parents=True)
        decision_path = staging / "R14_DECISION.json"
        inventory_path = staging / "R14_SOURCE_INVENTORY.json"
        html_path = staging / "rapport_final_r14.html"
        pdf_path = staging / "rapport_final_r14.pdf"
        _write_json(decision_path, r14_decision)
        _write_json(
            inventory_path,
            {
                "schema_version": "1.0.0",
                "report_id": REPORT_ID,
                "template_path": template_path.as_posix(),
                "template_sha256": file_sha256(template_path),
                "sources": sources,
            },
        )
        body = _report_body(receipt, run_manifest, decision, projection, sources)
        html_path.write_text(
            _render_html(template_path.read_text(encoding="utf-8"), body),
            encoding="utf-8",
        )
        page_count = _run_chrome(chrome_path, html_path, pdf_path)
        manifest = {
            "schema_version": "1.0.0",
            "report_id": REPORT_ID,
            "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "release_id": receipt["release_id"],
            "closed_run_id": receipt["run_id"],
            "scope": "R14_TECHNICAL_CLOSURE_ONLY",
            "rendering": {
                "engine": "PAGEDJS_0.4.3_THEN_CHROME_DEVTOOLS_PRINT",
                "page_count": page_count,
                "paper": "A4",
            },
            "artifacts": [
                {"relative_path": path.name, "sha256": file_sha256(path), "size_bytes": path.stat().st_size}
                for path in (decision_path, inventory_path, html_path, pdf_path)
            ],
            "qualification_flags": r14_decision["qualification_flags"],
        }
        _write_json(staging / "R14_REPORT_MANIFEST.json", manifest)
        staging.rename(output_root)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return validate_r14_report(repository_root, output_root)


def validate_r14_report(repository_root: Path, output_root: Path) -> dict[str, Any]:
    """Vérifie les artefacts, leurs empreintes et les restrictions de décision."""
    repository_root = repository_root.resolve(strict=True)
    output_root = output_root.resolve(strict=True)
    manifest = load_json(output_root / "R14_REPORT_MANIFEST.json")
    decision = load_json(output_root / "R14_DECISION.json")
    if manifest.get("report_id") != REPORT_ID or decision.get("decision_id") != REPORT_ID:
        _fail("r14_report_identity_mismatch")
    if decision.get("scientific_qualification") != "NOT_GRANTED":
        _fail("r14_scientific_decision_mismatch")
    expected_flags = {
        "calibrated": False,
        "founder_effect_proven": False,
        "production_promotion": False,
    }
    if decision.get("qualification_flags") != expected_flags:
        _fail("r14_qualification_flags_mismatch")
    for artifact in manifest.get("artifacts", []):
        path = output_root / artifact["relative_path"]
        if not path.is_file() or file_sha256(path) != artifact.get("sha256"):
            _fail("r14_artifact_hash_mismatch", artifact.get("relative_path"))
        if path.stat().st_size != artifact.get("size_bytes"):
            _fail("r14_artifact_size_mismatch", artifact.get("relative_path"))
    html_text = (output_root / "rapport_final_r14.html").read_text(encoding="utf-8")
    required = (
        REPORT_TITLE,
        "25 / 25",
        "qualification scientifique non accordée",
        "Figure 1 :",
        "Figure 4 :",
        "Table 1 :",
        "Table 8 :",
        "page-footer",
        "Paged.js v0.4.3",
    )
    for marker in required:
        if marker not in html_text:
            _fail("r14_html_required_marker_missing", marker)
    forbidden = ("[To complete:", "data/input/", "calibrated=true")
    for marker in forbidden:
        if marker in html_text:
            _fail("r14_html_forbidden_marker", marker)
    pdf = output_root / "rapport_final_r14.pdf"
    if not pdf.read_bytes().startswith(b"%PDF-") or pdf.stat().st_size < 100_000:
        _fail("r14_pdf_invalid")
    if manifest.get("rendering", {}).get("page_count", 0) < 12:
        _fail("r14_pdf_page_count_invalid")
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--r13-workspace", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--chrome", type=Path, required=True)
    parser.add_argument("--validate-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.validate_only:
        manifest = validate_r14_report(arguments.repository_root, arguments.output_root)
    else:
        manifest = build_r14_report(
            arguments.repository_root,
            arguments.candidate_root,
            arguments.r13_workspace,
            arguments.output_root,
            arguments.template,
            arguments.chrome,
        )
    print(json.dumps(manifest, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
