"""Produit des vues exploratoires pseudonymisées des résultats de l'étape 16A."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from effet_fondateur.contracts import validate_json_document, validate_tsv_table


COLORS = {
    "AFR": "#7c3aed",
    "AMR": "#ea580c",
    "EAS": "#0891b2",
    "EUR": "#2563eb",
    "SAS": "#16a34a",
}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="plot-reference-ancestry")
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    return parser


def _save(figure: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.svg"]
    figure.savefig(paths[0], dpi=180, bbox_inches="tight", facecolor="white")
    figure.savefig(paths[1], bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return paths


def _scores(path: Path) -> list[dict[str, object]]:
    return validate_tsv_table(path, "ancestry_scores.schema.json").rows


def _xy(rows: list[dict[str, object]]) -> np.ndarray:
    return np.asarray([[float(row["PC1"]), float(row["PC2"])] for row in rows])


def _scatter_scope(
    rows: list[dict[str, object]], scope: str, output_dir: Path
) -> list[Path]:
    scoped = [row for row in rows if row["ANALYSIS_SCOPE"] == scope]
    references = [row for row in scoped if row["REFERENCE_INCLUDED"]]
    study = [row for row in scoped if row["PROJECTED"]]
    figure, axis = plt.subplots(figsize=(10.5, 7.5))
    for superpopulation in COLORS:
        selected = [row for row in references if row["SUPERPOPULATION"] == superpopulation]
        points = _xy(selected)
        axis.scatter(
            points[:, 0], points[:, 1], s=8, alpha=0.22,
            color=COLORS[superpopulation], label=f"1000G {superpopulation} (n={len(selected)})",
            linewidths=0,
        )
    if scope == "GLOBAL":
        points = _xy(study)
        axis.scatter(
            points[:, 0], points[:, 1], s=34, marker="x", color="#111827",
            linewidths=1.2, label=f"Étude projetée (n={len(study)})",
        )
        title = "PCA globale — projection de l’étude sur les références 1000G"
    else:
        noncarriers = [row for row in study if row["TARGET_COPY_STATUS"] == "NON_CARRIER_COPY"]
        carriers = [row for row in study if row["TARGET_COPY_STATUS"] == "CARRIER_COPY"]
        unreliable = [row for row in study if row["TARGET_COPY_STATUS"] == "CARRIER_COPY_UNRELIABLE"]
        points = _xy(noncarriers)
        axis.scatter(points[:, 0], points[:, 1], s=22, marker="x", color="#475569", alpha=0.65,
                     label=f"Copies non porteuses (n={len(noncarriers)})")
        points = _xy(carriers)
        axis.scatter(points[:, 0], points[:, 1], s=100, marker="*", color="#dc2626",
                     edgecolors="white", linewidths=0.5, label=f"Copies porteuses fiables (n={len(carriers)})")
        if unreliable:
            points = _xy(unreliable)
            axis.scatter(points[:, 0], points[:, 1], s=100, marker="*", color="#f59e0b",
                         edgecolors="#111827", linewidths=0.6,
                         label=f"Copie porteuse incertaine (n={len(unreliable)})")
        title = "PCA locale — haplotypes autour de la variation cible"
    axis.set_title(title, fontsize=15, fontweight="bold")
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    axis.grid(alpha=0.15)
    axis.legend(loc="best", fontsize=9, frameon=True)
    axis.text(
        0.01, -0.13,
        "Positionnement relatif uniquement : aucune attribution ethnique, preuve d’ascendance locale ou preuve IBD.",
        transform=axis.transAxes, fontsize=9, color="#7c2d12",
    )
    return _save(figure, output_dir, f"pca_{scope.lower()}_reference_cloud")


def _local_study_zoom(rows: list[dict[str, object]], output_dir: Path) -> list[Path]:
    study = [
        row for row in rows
        if row["ANALYSIS_SCOPE"] == "LOCAL" and row["PROJECTED"]
    ]
    figure, axis = plt.subplots(figsize=(9, 7))
    categories = (
        ("NON_CARRIER_COPY", "Copies non porteuses", "#64748b", "o", 28, 0.45),
        ("CARRIER_COPY", "Copies porteuses fiables", "#dc2626", "*", 150, 1.0),
        ("CARRIER_COPY_UNRELIABLE", "Copie porteuse incertaine", "#f59e0b", "*", 150, 1.0),
    )
    for status, label, color, marker, size, alpha in categories:
        selected = [row for row in study if row["TARGET_COPY_STATUS"] == status]
        if not selected:
            continue
        points = _xy(selected)
        axis.scatter(points[:, 0], points[:, 1], s=size, marker=marker, color=color,
                     alpha=alpha, edgecolors="white", linewidths=0.5,
                     label=f"{label} (n={len(selected)})")
    axis.set_title("Zoom sur les haplotypes de l’étude — PCA locale", fontsize=15, fontweight="bold")
    axis.set_xlabel("PC1 local")
    axis.set_ylabel("PC2 local")
    axis.grid(alpha=0.18)
    axis.legend(loc="best")
    axis.text(
        0.01, -0.14,
        "Un regroupement des copies porteuses est compatible avec un haplotype local partagé, sans démontrer une origine fondatrice.",
        transform=axis.transAxes, fontsize=9, color="#7c2d12",
    )
    return _save(figure, output_dir, "pca_local_study_carrier_zoom")


def _scree(eigenvalue_path: Path, output_dir: Path) -> list[Path]:
    rows = validate_tsv_table(eigenvalue_path, "ancestry_eigenvalues.schema.json").rows
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=False)
    for axis, scope in zip(axes, ("GLOBAL", "LOCAL")):
        selected = [row for row in rows if row["ANALYSIS_SCOPE"] == scope]
        ratios = [100 * float(row["EXPLAINED_VARIANCE_RATIO"]) for row in selected]
        labels = [str(row["COMPONENT"]) for row in selected]
        axis.bar(labels, ratios, color="#2563eb" if scope == "GLOBAL" else "#7c3aed", alpha=0.82)
        axis.set_title(f"PCA {scope.lower()}")
        axis.set_ylabel("Variance expliquée (%)")
        axis.tick_params(axis="x", rotation=45)
        axis.grid(axis="y", alpha=0.18)
        axis.text(0.02, 0.96, f"PC1+PC2 = {sum(ratios[:2]):.1f}%", transform=axis.transAxes,
                  va="top", fontweight="bold")
    figure.suptitle("Variance expliquée par les dix premières composantes", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return _save(figure, output_dir, "pca_variance_explained")


def _nearest_centroid_counts(
    rows: list[dict[str, object]], centroid_path: Path
) -> dict[str, dict[str, dict[str, int]]]:
    centroids = validate_tsv_table(
        centroid_path, "ancestry_population_centroids.schema.json"
    ).rows
    result: dict[str, dict[str, dict[str, int]]] = {}
    for scope in ("GLOBAL", "LOCAL"):
        centers = {
            str(row["GROUP_ID"]): np.asarray([float(row[f"PC{i}"]) for i in range(1, 5)])
            for row in centroids
            if row["ANALYSIS_SCOPE"] == scope and row["GROUP_LEVEL"] == "SUPERPOPULATION"
        }
        projected = [row for row in rows if row["ANALYSIS_SCOPE"] == scope and row["PROJECTED"]]
        counters: dict[str, Counter[str]] = {}
        for row in projected:
            values = np.asarray([float(row[f"PC{i}"]) for i in range(1, 5)])
            nearest = min(centers, key=lambda name: float(np.linalg.norm(values - centers[name])))
            counters.setdefault(str(row["TARGET_COPY_STATUS"]), Counter())[nearest] += 1
        result[scope] = {status: dict(sorted(counter.items())) for status, counter in counters.items()}
    return result


def _local_group_dispersion(rows: list[dict[str, object]]) -> dict[str, dict[str, float | int]]:
    """Décrit la compacité interne sans transformer celle-ci en test d'IBD."""

    local = [row for row in rows if row["ANALYSIS_SCOPE"] == "LOCAL" and row["PROJECTED"]]
    result: dict[str, dict[str, float | int]] = {}
    for status in ("CARRIER_COPY", "NON_CARRIER_COPY"):
        selected = [row for row in local if row["TARGET_COPY_STATUS"] == status]
        values = np.asarray([
            [float(row[f"PC{i}"]) for i in range(1, 11)] for row in selected
        ])
        distances = np.linalg.norm(values - values.mean(axis=0), axis=1)
        result[status] = {
            "entity_count": len(selected),
            "median_distance_to_group_centroid_pc1_pc10": float(np.median(distances)),
            "rms_distance_to_group_centroid_pc1_pc10": float(np.sqrt(np.mean(distances ** 2))),
        }
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def execute(run_dir: Path, output_dir: Path | None = None) -> Path:
    base = run_dir / "stages/16A_analyze_reference_ancestry/ancestry"
    scores_path = base / "ancestry_scores.tsv"
    eigenvalues_path = base / "ancestry_eigenvalues.tsv"
    centroids_path = base / "ancestry_population_centroids.tsv"
    summary_path = base / "reference_ancestry_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    validate_json_document(summary, "reference_ancestry_summary.schema.json")
    rows = _scores(scores_path)
    destination = output_dir or run_dir / "derived/reference_ancestry_visualizations"
    destination.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []
    generated += _scatter_scope(rows, "GLOBAL", destination)
    generated += _scatter_scope(rows, "LOCAL", destination)
    generated += _local_study_zoom(rows, destination)
    generated += _scree(eigenvalues_path, destination)
    metrics = {
        "schema_version": "1.0.0",
        "run_id": run_dir.name,
        "source_summary": summary,
        "nearest_superpopulation_centroid_descriptive_only": _nearest_centroid_counts(rows, centroids_path),
        "local_group_dispersion_descriptive_only": _local_group_dispersion(rows),
        "warnings": [
            "Les centroïdes 1000G sont des repères larges et non des attributions d'identité.",
            "Le regroupement local des copies porteuses ne prouve ni IBD ni effet fondateur.",
            "Les figures sont exploratoires et doivent être interprétées avec IBS, LD, ROH et datation séparément.",
        ],
        "sources": {
            path.name: _sha256(path)
            for path in (scores_path, eigenvalues_path, centroids_path, summary_path)
        },
        "figures": {path.name: _sha256(path) for path in generated},
    }
    metrics_path = destination / "visualization_summary.json"
    metrics_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    cards = "".join(
        f'<figure><img src="{html.escape(path.name)}" alt="{html.escape(path.stem)}"><figcaption>{html.escape(path.stem)}</figcaption></figure>'
        for path in generated if path.suffix == ".png"
    )
    gallery = destination / "index.html"
    gallery.write_text(
        "<!doctype html><html lang=\"fr\"><meta charset=\"utf-8\"><title>Visualisations 16A</title>"
        "<style>body{font-family:system-ui;margin:2rem;background:#f8fafc;color:#172554}"
        "figure{background:white;padding:1rem;margin:1.5rem 0;border:1px solid #cbd5e1}"
        "img{max-width:100%;height:auto}figcaption{font-weight:700;margin-top:.5rem}"
        ".warning{color:#7c2d12;background:#fff7ed;padding:1rem}</style>"
        f"<h1>Étape 16A — positionnement relatif aux références 1000G</h1>"
        f"<p>Global : {summary['global']['informative_variant_count']} variants ; "
        f"local : {summary['local']['informative_variant_count']} variants.</p>"
        "<p class=\"warning\">Aucune figure ne constitue une attribution ethnique, une preuve d’ascendance locale, d’IBD ou d’effet fondateur.</p>"
        f"{cards}</html>\n",
        encoding="utf-8",
    )
    return destination


def main() -> None:
    arguments = _parser().parse_args()
    print(execute(arguments.run_dir.resolve(), arguments.output_dir))


if __name__ == "__main__":
    main()
