"""Génère les figures ROH et âge utilisées dans le Résultat 7 du rapport."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import shutil
import tempfile
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TARGET_VARIANT_ID = "DOCK6_GRCh38_19_11237780_C_A"
TARGET_POSITION_BP = 11_237_780
FAMILY_COLORS = {"FAM001": "#2563eb", "FAM002": "#7c3aed", "FAM003": "#ea580c"}


class Result7VisualizationError(ValueError):
    """Signale une source manquante ou incohérente pour les figures."""


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="plot-result7-roh-age")
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    return parser


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise Result7VisualizationError(f"empty_source:{path}")
    return rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _common_overlap(rows: Sequence[dict[str, str]]) -> tuple[int, int]:
    """Retourne l'intersection physique de ROH couvrant la cible."""

    evaluated = [row for row in rows if row["TARGET_IN_ROH"] == "true"]
    if not evaluated:
        raise Result7VisualizationError("no_target_covering_roh")
    start = max(int(row["ROH_POS1_BP"]) for row in evaluated)
    end = min(int(row["ROH_POS2_BP"]) for row in evaluated)
    if start > TARGET_POSITION_BP or end < TARGET_POSITION_BP or start >= end:
        raise Result7VisualizationError("target_roh_overlap_missing")
    return start, end


def _quantiles(values: Sequence[float]) -> dict[str, float]:
    if not values:
        raise Result7VisualizationError("empty_quantile_input")
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(np.min(array)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.quantile(array, 0.5)),
        "q75": float(np.quantile(array, 0.75)),
        "maximum": float(np.max(array)),
    }


def _save(figure: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.svg"]
    figure.savefig(paths[0], dpi=190, bbox_inches="tight", facecolor="white")
    figure.savefig(paths[1], bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return paths


def _roh_figure(
    target_rows: Sequence[dict[str, str]], burden_rows: Sequence[dict[str, str]], output_dir: Path,
) -> tuple[list[Path], dict[str, object]]:
    homozygous_covering = [
        row for row in target_rows
        if row["TARGET_ZYGOSITY"] == "HOMOZYGOUS"
        and row["TARGET_GENOTYPE"] == "A/A"
        and row["TARGET_IN_ROH"] == "true"
    ]
    if len(homozygous_covering) != 5:
        raise Result7VisualizationError("unexpected_homozygous_target_roh_count")
    overlap_start, overlap_end = _common_overlap(homozygous_covering)
    independent_ids = {
        row["SAMPLE_ID"] for row in homozygous_covering
        if "target_carriers_independent" in row["ANALYTIC_ROLES"]
    }
    if len(independent_ids) != 3:
        raise Result7VisualizationError("unexpected_independent_target_roh_count")

    controls = [row for row in burden_rows if row["COHORT_ID"] == "controls_unrelated"]
    carriers = [row for row in burden_rows if row["COHORT_ID"] == "target_carriers_independent"]
    if len(controls) != 62 or len(carriers) != 3:
        raise Result7VisualizationError("unexpected_roh_burden_cohort_size")
    control_burden_mb = [float(row["TOTAL_ROH_KB"]) / 1000 for row in controls]
    carrier_burden_mb = [float(row["TOTAL_ROH_KB"]) / 1000 for row in carriers]
    carrier_by_id = {row["SAMPLE_ID"]: row for row in carriers}

    figure, axes = plt.subplots(1, 2, figsize=(15, 6.4), gridspec_kw={"width_ratios": [1.28, 1]})
    local_axis = axes[0]
    ordered = sorted(homozygous_covering, key=lambda row: (row["FID"], row["SAMPLE_ID"]))
    labels = []
    for index, row in enumerate(ordered):
        y = len(ordered) - index
        start_mb = int(row["ROH_POS1_BP"]) / 1_000_000
        end_mb = int(row["ROH_POS2_BP"]) / 1_000_000
        is_independent = row["SAMPLE_ID"] in independent_ids
        local_axis.plot(
            [start_mb, end_mb], [y, y], color=FAMILY_COLORS[row["FID"]],
            linewidth=10 if is_independent else 6, solid_capstyle="round", alpha=0.95 if is_independent else 0.55,
        )
        local_axis.scatter(
            [start_mb, end_mb], [y, y], s=34, color=FAMILY_COLORS[row["FID"]], zorder=3,
        )
        suffix = "indépendant" if is_independent else "autre A/A"
        labels.append(f"{row['FID']} — {suffix}")
    local_axis.axvspan(overlap_start / 1_000_000, overlap_end / 1_000_000,
                       color="#fbbf24", alpha=0.20, label="Recouvrement commun")
    local_axis.axvline(TARGET_POSITION_BP / 1_000_000, color="#dc2626", linewidth=2.2,
                       label="Variation DOCK6")
    local_axis.set_yticks(range(len(ordered), 0, -1), labels)
    local_axis.set_xlabel("Position sur le chromosome 19 (Mb, GRCh38)")
    local_axis.set_title("ROH couvrant la variation chez les cinq homozygotes A/A")
    local_axis.grid(axis="x", alpha=0.18)
    local_axis.legend(fontsize=8.5, loc="lower left")

    burden_axis = axes[1]
    rng = np.random.default_rng(20260817)
    control_y = 1 + rng.uniform(-0.08, 0.08, len(control_burden_mb))
    burden_axis.boxplot(
        [control_burden_mb], positions=[1], orientation="horizontal", widths=0.32, showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": "#dbeafe", "edgecolor": "#2563eb"},
        medianprops={"color": "#172554", "linewidth": 2},
        whiskerprops={"color": "#64748b"}, capprops={"color": "#64748b"},
    )
    burden_axis.scatter(control_burden_mb, control_y, s=22, alpha=0.45, color="#64748b", label="62 témoins")
    for row in carriers:
        x = float(row["TOTAL_ROH_KB"]) / 1000
        family = row["FID"]
        burden_axis.scatter(x, 2, s=95, marker="D", color=FAMILY_COLORS[family],
                            edgecolor="white", linewidth=0.8, zorder=4)
        burden_axis.annotate(family.replace("FAM00", "F"), (x, 2), xytext=(0, 10),
                             textcoords="offset points", ha="center", fontsize=9, fontweight="bold")
    burden_axis.set_yticks([1, 2], ["Témoins", "3 familles\nindépendantes"])
    burden_axis.set_xlabel("Longueur totale des ROH autosomiques (Mb)")
    burden_axis.set_title("Charge génomique en ROH")
    burden_axis.set_ylim(0.55, 2.45)
    burden_axis.grid(axis="x", alpha=0.18)
    burden_axis.legend(fontsize=9, loc="lower right")

    figure.suptitle("Régions d’homozygotie : signal local autour de DOCK6 et contexte génomique",
                    fontsize=16, fontweight="bold")
    figure.text(
        0.5, 0.01,
        "Les trois valeurs familiales sont descriptives : n=3 ne permet pas une comparaison statistique de charge avec les témoins.",
        ha="center", fontsize=9, color="#7c2d12",
    )
    figure.tight_layout(rect=(0, 0.045, 1, 0.94))
    metrics = {
        "homozygous_alt_target_in_roh_count": len(homozygous_covering),
        "independent_target_in_roh_count": len(independent_ids),
        "control_target_in_roh_count": sum(
            row["TARGET_IN_ROH"] == "true" and "controls_unrelated" in row["ANALYTIC_ROLES"]
            for row in target_rows
        ),
        "common_overlap_start_bp": overlap_start,
        "common_overlap_end_bp": overlap_end,
        "common_overlap_length_kb": (overlap_end - overlap_start + 1) / 1000,
        "control_total_roh_mb_quantiles": _quantiles(control_burden_mb),
        "independent_family_total_roh_mb": {
            family: float(carrier_by_id[sample_id]["TOTAL_ROH_KB"]) / 1000
            for family, sample_id in {
                "FAM001": "FAM001_ENFANT", "FAM002": "FAM002_MERE", "FAM003": "FAM003_ENFANT"
            }.items()
        },
    }
    return _save(figure, output_dir, "roh_local_and_genomewide"), metrics


def _age_figure(
    estimate_rows: Sequence[dict[str, str]], scenario_rows: Sequence[dict[str, str]], output_dir: Path,
) -> tuple[list[Path], dict[str, object]]:
    model_rows = {row["MODEL"]: row for row in estimate_rows}
    if set(model_rows) != {"CORRELATED", "INDEPENDENT"}:
        raise Result7VisualizationError("unexpected_age_models")
    generation_rows = [row for row in scenario_rows if row["SCENARIO_TYPE"] == "GENERATION_INTERVAL"]
    if [int(row["GENERATION_YEARS"]) for row in generation_rows] != [25, 28, 30]:
        raise Result7VisualizationError("unexpected_generation_intervals")

    figure, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), gridspec_kw={"width_ratios": [1, 1.25]})
    generation_axis = axes[0]
    models = [
        ("INDEPENDENT", "Modèle indépendant\n(sensibilité)", "#64748b"),
        ("CORRELATED", "Modèle corrélé\n(principal, exploratoire)", "#ea580c"),
    ]
    for y, (model, label, color) in enumerate(models, start=1):
        row = model_rows[model]
        estimate = float(row["ESTIMATE_GENERATIONS"])
        lower = float(row["CI_LOWER_GENERATIONS"])
        upper = float(row["CI_UPPER_GENERATIONS"])
        generation_axis.errorbar(
            estimate, y, xerr=[[estimate - lower], [upper - estimate]], fmt="o",
            markersize=9, color=color, ecolor=color, capsize=5, linewidth=2.2,
        )
        generation_axis.annotate(f"{estimate:.1f} générations", (estimate, y), xytext=(0, 11),
                                 textcoords="offset points", ha="center", fontsize=9, fontweight="bold")
    generation_axis.set_yticks([1, 2], [models[0][1], models[1][1]])
    generation_axis.set_xlabel("Âge estimé (générations)")
    generation_axis.set_title("Estimation Gamma et intervalle de confiance à 95 %")
    generation_axis.set_xlim(0, 105)
    generation_axis.set_ylim(0.65, 2.35)
    generation_axis.grid(axis="x", alpha=0.18)

    calendar_axis = axes[1]
    for y, row in enumerate(generation_rows, start=1):
        estimate = float(row["ESTIMATE_YEARS"])
        lower = float(row["CI_LOWER_YEARS"])
        upper = float(row["CI_UPPER_YEARS"])
        interval = int(row["GENERATION_YEARS"])
        calendar_axis.errorbar(
            estimate, y, xerr=[[estimate - lower], [upper - estimate]], fmt="o",
            markersize=9, color="#2563eb", ecolor="#60a5fa", capsize=5, linewidth=2.2,
        )
        calendar_axis.annotate(f"≈ {estimate:,.0f} ans".replace(",", " "), (estimate, y),
                               xytext=(0, 11), textcoords="offset points", ha="center",
                               fontsize=9, fontweight="bold")
        calendar_axis.text(25, y, f"{interval} ans/génération", va="center", fontsize=9, color="#334155")
    calendar_axis.set_yticks([1, 2, 3], ["25 ans", "28 ans", "30 ans"])
    calendar_axis.set_xlabel("Conversion calendaire (années avant le présent)")
    calendar_axis.set_title("Sensibilité à la durée d’une génération")
    calendar_axis.set_xlim(0, 3_050)
    calendar_axis.set_ylim(0.65, 3.35)
    calendar_axis.grid(axis="x", alpha=0.18)

    figure.suptitle("Âge exploratoire de la variation DOCK6", fontsize=16, fontweight="bold")
    figure.text(
        0.5, 0.01,
        "Trois familles indépendantes seulement : les intervalles sont larges et aucune datation historique précise n’est revendiquée.",
        ha="center", fontsize=9, color="#7c2d12",
    )
    figure.tight_layout(rect=(0, 0.05, 1, 0.93))
    correlated = model_rows["CORRELATED"]
    metrics = {
        "status": correlated["ANALYSIS_STATUS"],
        "included_independent_units": int(correlated["N_UNITS"]),
        "correlated_estimate_generations": float(correlated["ESTIMATE_GENERATIONS"]),
        "correlated_ci95_generations": [
            float(correlated["CI_LOWER_GENERATIONS"]), float(correlated["CI_UPPER_GENERATIONS"]),
        ],
        "independent_estimate_generations": float(model_rows["INDEPENDENT"]["ESTIMATE_GENERATIONS"]),
        "generation_interval_scenarios": {
            row["GENERATION_YEARS"]: {
                "estimate_years": float(row["ESTIMATE_YEARS"]),
                "ci95_years": [float(row["CI_LOWER_YEARS"]), float(row["CI_UPPER_YEARS"])],
            }
            for row in generation_rows
        },
    }
    return _save(figure, output_dir, "variant_age_forest"), metrics


def execute(run_dir: Path, output_dir: Path | None = None) -> Path:
    run_dir = run_dir.resolve()
    destination = (output_dir or run_dir / "derived/result7_visualizations").resolve()
    if destination.exists():
        raise Result7VisualizationError(f"output_dir_already_exists:{destination}")
    roh_dir = run_dir / "stages/16_analyze_roh/roh"
    age_dir = run_dir / "stages/14_estimate_variant_age/dating"
    sources = {
        "target_in_roh": roh_dir / "target_in_roh.tsv",
        "roh_burden": roh_dir / "roh_burden.tsv",
        "roh_summary": roh_dir / "roh_analysis_summary.json",
        "age_estimates": age_dir / "variant_age_estimates.tsv",
        "age_scenarios": age_dir / "variant_age_scenarios.tsv",
        "age_summary": age_dir / "variant_age_summary.json",
        "ld_v2_summary": run_dir / "derived/phase15_ld_v2/analysis_summary.json",
    }
    missing = [str(path) for path in sources.values() if not path.is_file()]
    if missing:
        raise Result7VisualizationError(f"missing_sources:{','.join(missing)}")
    target_rows = _read_tsv(sources["target_in_roh"])
    burden_rows = _read_tsv(sources["roh_burden"])
    estimate_rows = _read_tsv(sources["age_estimates"])
    scenario_rows = _read_tsv(sources["age_scenarios"])
    if {row["TARGET_VARIANT_ID"] for row in target_rows} != {TARGET_VARIANT_ID}:
        raise Result7VisualizationError("target_variant_mismatch")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        roh_paths, roh_metrics = _roh_figure(target_rows, burden_rows, temporary)
        age_paths, age_metrics = _age_figure(estimate_rows, scenario_rows, temporary)
        summary = {
            "schema_version": "1.0.0",
            "run_id": run_dir.name,
            "target_variant_id": TARGET_VARIANT_ID,
            "roh": roh_metrics,
            "age": age_metrics,
            "ld_v2": json.loads(sources["ld_v2_summary"].read_text(encoding="utf-8")),
            "source_sha256": {name: _sha256(path) for name, path in sources.items()},
            "interpretation_limits": [
                "ROH locale compatible avec une autozygotie autour de la variation, sans preuve autonome d'effet fondateur.",
                "Charge ROH des trois familles descriptive uniquement; aucune comparaison inférentielle.",
                "Datation exploratoire avec trois unités indépendantes et intervalle de confiance large.",
            ],
        }
        (temporary / "visualization_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8",
        )
        cards = "".join(
            f'<figure><img src="{html.escape(path.name)}" alt="{html.escape(path.stem)}">'
            f'<figcaption>{html.escape(path.stem)}</figcaption></figure>'
            for path in (*roh_paths, *age_paths) if path.suffix == ".png"
        )
        (temporary / "index.html").write_text(
            "<!doctype html><html lang=\"fr\"><meta charset=\"utf-8\"><title>Résultat 7 — ROH et âge</title>"
            "<style>body{font-family:system-ui;margin:2rem;background:#f8fafc;color:#172554}"
            "figure{background:white;padding:1rem;margin:1.5rem 0;border:1px solid #cbd5e1}"
            "img{max-width:100%;height:auto}figcaption{font-weight:700;margin-top:.5rem}</style>"
            "<h1>Résultat 7 — ROH et âge de la variation</h1>" + cards + "</html>\n",
            encoding="utf-8",
        )
        outputs = sorted(path for path in temporary.iterdir() if path.is_file())
        (temporary / "checksums.sha256").write_text(
            "".join(f"{_sha256(path)}  {path.name}\n" for path in outputs), encoding="utf-8",
        )
        temporary.rename(destination)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination


def main() -> int:
    args = _parser().parse_args()
    try:
        destination = execute(args.run_dir, args.output_dir)
    except (OSError, Result7VisualizationError) as error:
        raise SystemExit(f"result7_visualization_failed:{error}") from error
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
