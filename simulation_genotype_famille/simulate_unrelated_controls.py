"""Simulation reproductible de témoins indépendants compatibles avec un PED/MAP."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


LOGGER = logging.getLogger(__name__)
DEFAULT_RANDOM_SEED = 20260804
VALID_ALLELES = frozenset("ACGT")


class ControlSimulationError(RuntimeError):
    """Erreur empêchant une simulation de témoins fiable."""


@dataclass(frozen=True)
class SourceIndividual:
    """Individu lu dans le PED source, génotypes inclus."""

    family_id: str
    individual_id: str
    paternal_id: str
    maternal_id: str
    sex: str
    phenotype: str
    columns: tuple[str, ...]


@dataclass(frozen=True)
class MarkerFrequency:
    """Fréquence allélique et taux manquant estimés pour un marqueur."""

    chromosome: str
    marker_id: str
    genetic_distance_cm: str
    position_bp: str
    allele_1: str
    allele_2: str
    allele_1_frequency: float
    missing_rate: float
    called_genotype_count: int


@dataclass(frozen=True)
class SourceDataset:
    """PED/MAP validé et statistiques nécessaires à la simulation."""

    individuals: tuple[SourceIndividual, ...]
    markers: tuple[MarkerFrequency, ...]


def ensure_outputs_are_available(output_paths: Iterable[Path], force: bool) -> None:
    """Refuse l'écrasement d'une simulation précédente sans accord explicite."""
    existing_paths = [path for path in output_paths if path.exists()]
    if existing_paths and not force:
        existing = "\n- ".join(str(path) for path in existing_paths)
        raise ControlSimulationError(
            "Sorties déjà présentes. Utiliser --force pour les remplacer:\n- "
            + existing
        )


def write_lines_atomic(output_path: Path, lines: Iterable[str]) -> None:
    """Écrit une sortie complète puis la remplace atomiquement."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f".{output_path.name}.tmp")
    try:
        with temporary_path.open("w", encoding="utf-8", newline="\n") as output_file:
            for line in lines:
                output_file.write(line)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def read_map(map_path: Path) -> list[tuple[str, str, str, str]]:
    """Lit les quatre colonnes MAP et valide les identifiants de marqueurs."""
    markers: list[tuple[str, str, str, str]] = []
    with map_path.open("r", encoding="utf-8") as map_file:
        for line_number, line in enumerate(map_file, start=1):
            columns = line.split()
            if len(columns) != 4:
                raise ControlSimulationError(
                    f"Ligne MAP invalide {map_path}:{line_number}, 4 colonnes attendues"
                )
            markers.append(tuple(columns))
    if not markers:
        raise ControlSimulationError(f"Aucun marqueur dans {map_path}")
    duplicate_markers = [
        marker_id
        for marker_id, count in Counter(marker[1] for marker in markers).items()
        if count > 1
    ]
    if duplicate_markers:
        raise ControlSimulationError(
            "Identifiants MAP dupliqués: " + ", ".join(sorted(duplicate_markers))
        )
    return markers


def read_ped(ped_path: Path, marker_count: int) -> list[SourceIndividual]:
    """Lit un PED et vérifie le nombre de génotypes par individu."""
    expected_column_count = 6 + 2 * marker_count
    individuals: list[SourceIndividual] = []
    with ped_path.open("r", encoding="utf-8") as ped_file:
        for line_number, line in enumerate(ped_file, start=1):
            columns = tuple(line.split())
            if len(columns) != expected_column_count:
                raise ControlSimulationError(
                    f"Ligne PED invalide {ped_path}:{line_number}: "
                    f"{len(columns)} colonnes, {expected_column_count} attendues"
                )
            individuals.append(
                SourceIndividual(
                    family_id=columns[0],
                    individual_id=columns[1],
                    paternal_id=columns[2],
                    maternal_id=columns[3],
                    sex=columns[4],
                    phenotype=columns[5],
                    columns=columns,
                )
            )
    if not individuals:
        raise ControlSimulationError(f"Aucun individu dans {ped_path}")
    duplicate_individual_ids = [
        individual_id
        for individual_id, count in Counter(
            individual.individual_id for individual in individuals
        ).items()
        if count > 1
    ]
    if duplicate_individual_ids:
        raise ControlSimulationError(
            "IID source dupliqués entre familles: "
            + ", ".join(sorted(duplicate_individual_ids))
        )
    duplicate_keys = [
        key
        for key, count in Counter(
            (individual.family_id, individual.individual_id)
            for individual in individuals
        ).items()
        if count > 1
    ]
    if duplicate_keys:
        raise ControlSimulationError(
            "Individus PED dupliqués: "
            + ", ".join(f"{family_id}/{individual_id}" for family_id, individual_id in duplicate_keys)
        )
    return individuals


def estimate_marker_frequencies(
    map_rows: list[tuple[str, str, str, str]],
    individuals: list[SourceIndividual],
) -> list[MarkerFrequency]:
    """Estime les fréquences à partir des génotypes appelés du PED source."""
    frequencies: list[MarkerFrequency] = []
    for marker_index, map_row in enumerate(map_rows):
        allele_counts: Counter[str] = Counter()
        missing_genotype_count = 0
        genotype_offset = 6 + 2 * marker_index
        for individual in individuals:
            allele_1, allele_2 = individual.columns[
                genotype_offset : genotype_offset + 2
            ]
            if allele_1 == "0" or allele_2 == "0":
                missing_genotype_count += 1
                continue
            if allele_1 not in VALID_ALLELES or allele_2 not in VALID_ALLELES:
                raise ControlSimulationError(
                    f"Allèle invalide au marqueur {map_row[1]}: {allele_1}/{allele_2}"
                )
            allele_counts.update((allele_1, allele_2))

        observed_alleles = sorted(allele_counts)
        if len(observed_alleles) > 2:
            raise ControlSimulationError(
                f"Plus de deux allèles observés au marqueur {map_row[1]}: "
                + ",".join(observed_alleles)
            )
        if not observed_alleles:
            allele_1 = "0"
            allele_2 = "0"
            allele_1_frequency = 0.0
        elif len(observed_alleles) == 1:
            allele_1 = observed_alleles[0]
            allele_2 = observed_alleles[0]
            allele_1_frequency = 1.0
        else:
            allele_1, allele_2 = observed_alleles
            allele_1_frequency = allele_counts[allele_1] / sum(allele_counts.values())

        frequencies.append(
            MarkerFrequency(
                chromosome=map_row[0],
                marker_id=map_row[1],
                genetic_distance_cm=map_row[2],
                position_bp=map_row[3],
                allele_1=allele_1,
                allele_2=allele_2,
                allele_1_frequency=allele_1_frequency,
                missing_rate=missing_genotype_count / len(individuals),
                called_genotype_count=len(individuals) - missing_genotype_count,
            )
        )
    return frequencies


def load_source_dataset(source_prefix: Path) -> SourceDataset:
    """Charge un préfixe PED/MAP et calcule les fréquences empiriques."""
    ped_path = source_prefix.with_suffix(".ped")
    map_path = source_prefix.with_suffix(".map")
    missing_paths = [path for path in (ped_path, map_path) if not path.is_file()]
    if missing_paths:
        raise ControlSimulationError(
            "Fichiers source introuvables: "
            + ", ".join(str(path) for path in missing_paths)
        )
    map_rows = read_map(map_path)
    individuals = read_ped(ped_path, len(map_rows))
    markers = estimate_marker_frequencies(map_rows, individuals)
    return SourceDataset(individuals=tuple(individuals), markers=tuple(markers))


def read_source_groups(groups_path: Path | None) -> dict[str, str]:
    """Lit les groupes source au format IID puis groupe."""
    if groups_path is None:
        return {}
    if not groups_path.is_file():
        raise ControlSimulationError(f"Fichier de groupes introuvable: {groups_path}")
    groups: dict[str, str] = {}
    with groups_path.open("r", encoding="utf-8") as groups_file:
        for line_number, line in enumerate(groups_file, start=1):
            columns = line.split()
            if len(columns) != 2:
                raise ControlSimulationError(
                    f"Ligne de groupes invalide {groups_path}:{line_number}"
                )
            individual_id, group = columns
            if individual_id in groups:
                raise ControlSimulationError(
                    f"IID dupliqué dans {groups_path}: {individual_id}"
                )
            groups[individual_id] = group
    return groups


def sample_marker_genotype(
    marker: MarkerFrequency,
    random_generator: np.random.Generator,
    simulate_missingness: bool,
) -> tuple[str, str]:
    """Tire deux allèles indépendants selon HWE pour un marqueur."""
    if marker.allele_1 == "0":
        return "0", "0"
    if simulate_missingness and random_generator.random() < marker.missing_rate:
        return "0", "0"
    if marker.allele_1 == marker.allele_2:
        return marker.allele_1, marker.allele_1
    sampled_alleles = tuple(
        marker.allele_1
        if random_generator.random() < marker.allele_1_frequency
        else marker.allele_2
        for _ in range(2)
    )
    return sampled_alleles


def build_simulated_controls(
    markers: tuple[MarkerFrequency, ...],
    control_count: int,
    random_seed: int,
    identifier_prefix: str,
    simulate_missingness: bool,
) -> list[SourceIndividual]:
    """Génère des fondateurs indépendants avec FID et IID uniques."""
    if control_count <= 0:
        raise ControlSimulationError("Le nombre de témoins doit être strictement positif")
    if not identifier_prefix or any(character.isspace() for character in identifier_prefix):
        raise ControlSimulationError("Le préfixe des identifiants est vide ou contient un espace")

    random_generator = np.random.default_rng(random_seed)
    width = max(4, len(str(control_count)))
    controls: list[SourceIndividual] = []
    for control_index in range(1, control_count + 1):
        individual_id = f"{identifier_prefix}_{control_index:0{width}d}"
        genotype_columns: list[str] = []
        for marker in markers:
            genotype_columns.extend(
                sample_marker_genotype(marker, random_generator, simulate_missingness)
            )
        sex = str(1 + int(random_generator.integers(0, 2)))
        columns = tuple(
            [individual_id, individual_id, "0", "0", sex, "1"]
            + genotype_columns
        )
        controls.append(
            SourceIndividual(
                family_id=individual_id,
                individual_id=individual_id,
                paternal_id="0",
                maternal_id="0",
                sex=sex,
                phenotype="1",
                columns=columns,
            )
        )
    return controls


def simulate_unrelated_controls(
    source_prefix: Path,
    output_prefix: Path,
    control_count: int,
    random_seed: int = DEFAULT_RANDOM_SEED,
    identifier_prefix: str = "SIMCTRL",
    source_groups_path: Path | None = None,
    include_source: bool = True,
    simulate_missingness: bool = True,
    force: bool = False,
) -> dict[str, object]:
    """Produit un PED/MAP avec témoins simulés et un rapport de traçabilité."""
    source_prefix = source_prefix.resolve()
    output_prefix = output_prefix.resolve()
    if output_prefix.parent == source_prefix.parent:
        raise ControlSimulationError(
            "Le dossier de sortie doit être distinct du dossier source pour "
            "protéger les données existantes"
        )
    dataset = load_source_dataset(source_prefix)

    if source_groups_path is None:
        candidate_groups_path = source_prefix.parent / "groupes.txt"
        source_groups_path = candidate_groups_path if candidate_groups_path.is_file() else None
    source_groups = read_source_groups(source_groups_path)

    controls = build_simulated_controls(
        markers=dataset.markers,
        control_count=control_count,
        random_seed=random_seed,
        identifier_prefix=identifier_prefix,
        simulate_missingness=simulate_missingness,
    )
    source_individual_ids = {
        individual.individual_id for individual in dataset.individuals
    }
    identifier_collisions = sorted(
        control.individual_id
        for control in controls
        if control.individual_id in source_individual_ids
    )
    if identifier_collisions:
        raise ControlSimulationError(
            "IID simulés déjà présents dans la source: "
            + ", ".join(identifier_collisions)
        )
    output_directory = output_prefix.parent
    output_paths = {
        "ped": output_prefix.with_suffix(".ped"),
        "map": output_prefix.with_suffix(".map"),
        "groups": output_directory / "groupes.txt",
        "cases": output_directory / "cas.txt",
        "controls": output_directory / "temoins.txt",
        "frequencies": output_directory / "simulation_marker_frequencies.tsv",
        "report": output_directory / "control_simulation_report.json",
    }
    ensure_outputs_are_available(output_paths.values(), force)

    output_individuals = list(dataset.individuals) if include_source else []
    output_individuals.extend(controls)
    write_lines_atomic(
        output_paths["ped"],
        ("\t".join(individual.columns) + "\n" for individual in output_individuals),
    )
    write_lines_atomic(
        output_paths["map"],
        (
            "\t".join(
                [
                    marker.chromosome,
                    marker.marker_id,
                    marker.genetic_distance_cm,
                    marker.position_bp,
                ]
            )
            + "\n"
            for marker in dataset.markers
        ),
    )

    group_lines: list[str] = []
    if include_source:
        for individual in dataset.individuals:
            default_group = (
                "ATTEINT" if individual.phenotype == "2" else "SOURCE_NON_ATTEINT"
            )
            group_lines.append(
                f"{individual.individual_id}\t"
                f"{source_groups.get(individual.individual_id, default_group)}\n"
            )
    group_lines.extend(
        f"{control.individual_id}\tTEMOIN_SIMULE\n" for control in controls
    )
    write_lines_atomic(output_paths["groups"], group_lines)
    write_lines_atomic(
        output_paths["cases"],
        (
            f"{individual.family_id}\t{individual.individual_id}\n"
            for individual in dataset.individuals
            if include_source and individual.phenotype == "2"
        ),
    )
    write_lines_atomic(
        output_paths["controls"],
        (
            f"{control.family_id}\t{control.individual_id}\n"
            for control in controls
        ),
    )

    frequency_lines = [
        "chromosome\tmarker_id\tposition_bp\tallele_1\tallele_2"
        "\tallele_1_frequency\tmissing_rate\tcalled_genotype_count\n"
    ]
    frequency_lines.extend(
        "\t".join(
            [
                marker.chromosome,
                marker.marker_id,
                marker.position_bp,
                marker.allele_1,
                marker.allele_2,
                f"{marker.allele_1_frequency:.8f}",
                f"{marker.missing_rate:.8f}",
                str(marker.called_genotype_count),
            ]
        )
        + "\n"
        for marker in dataset.markers
    )
    write_lines_atomic(output_paths["frequencies"], frequency_lines)

    report: dict[str, object] = {
        "source_prefix": str(source_prefix),
        "output_prefix": str(output_prefix),
        "source_individual_count": len(dataset.individuals),
        "simulated_control_count": len(controls),
        "output_individual_count": len(output_individuals),
        "marker_count": len(dataset.markers),
        "random_seed": random_seed,
        "identifier_prefix": identifier_prefix,
        "include_source": include_source,
        "simulate_missingness": simulate_missingness,
        "frequency_source": "empirical_source_ped",
        "monomorphic_marker_count": sum(
            marker.allele_1 == marker.allele_2 and marker.allele_1 != "0"
            for marker in dataset.markers
        ),
        "all_missing_marker_count": sum(
            marker.allele_1 == "0" for marker in dataset.markers
        ),
        "limitations": [
            "Fréquences estimées sur une petite cohorte familiale enrichie en atteints.",
            "Marqueurs simulés indépendamment sous HWE: le LD et les haplotypes ne sont pas reproduits.",
            "Une fluctuation aléatoire peut placer certaines paires au-dessus d'un seuil KING faible.",
            "Témoins synthétiques réservés aux tests techniques et analyses exploratoires.",
            "Ne pas les utiliser comme preuve biologique d'un effet fondateur.",
        ],
        "outputs": {name: str(path) for name, path in output_paths.items()},
    }
    write_lines_atomic(
        output_paths["report"],
        [json.dumps(report, indent=2, ensure_ascii=False) + "\n"],
    )
    return report


def build_argument_parser() -> argparse.ArgumentParser:
    """Construit l'interface en ligne de commande du simulateur."""
    parser = argparse.ArgumentParser(
        description=(
            "Simule des témoins indépendants sous HWE à partir d'un PED/MAP. "
            "Ces témoins synthétiques ne reproduisent pas le LD."
        )
    )
    parser.add_argument("--source-prefix", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--n-controls", type=int)
    parser.add_argument("--seed", type=int, default=DEFAULT_RANDOM_SEED)
    parser.add_argument("--identifier-prefix", default="SIMCTRL")
    parser.add_argument("--source-groups", type=Path)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--no-simulate-missingness", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    """Demande le nombre de témoins si nécessaire puis lance la simulation."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    arguments = build_argument_parser().parse_args()
    try:
        control_count = arguments.n_controls
        if control_count is None:
            raw_value = input("Combien de témoins non apparentés simuler ? ").strip()
            try:
                control_count = int(raw_value)
            except ValueError as error:
                raise ControlSimulationError(
                    f"Nombre de témoins invalide: {raw_value}"
                ) from error
        report = simulate_unrelated_controls(
            source_prefix=arguments.source_prefix,
            output_prefix=arguments.output_prefix,
            control_count=control_count,
            random_seed=arguments.seed,
            identifier_prefix=arguments.identifier_prefix,
            source_groups_path=arguments.source_groups,
            include_source=not arguments.controls_only,
            simulate_missingness=not arguments.no_simulate_missingness,
            force=arguments.force,
        )
        LOGGER.warning(
            "Témoins synthétiques sans LD: ne pas les interpréter comme témoins biologiques."
        )
        LOGGER.info(
            "Simulation terminée: %s témoins, %s marqueurs, graine %s",
            report["simulated_control_count"],
            report["marker_count"],
            report["random_seed"],
        )
        return 0
    except ControlSimulationError as error:
        LOGGER.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
