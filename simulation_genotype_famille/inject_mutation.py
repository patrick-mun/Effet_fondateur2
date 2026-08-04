"""Injection traçable d'une mutation dans un jeu PLINK PED/MAP."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import shutil
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


LOGGER = logging.getLogger(__name__)
VALID_GENOTYPES = frozenset({"REF/REF", "REF/ALT", "ALT/ALT", "MISSING"})
VALID_BASES = frozenset("ACGT")


class MutationInjectionError(RuntimeError):
    """Erreur empêchant une injection de mutation sûre."""


@dataclass(frozen=True)
class MutationMetadata:
    """Métadonnées génomiques et cliniques associées à la mutation."""

    gene: str
    mutation_name: str
    chromosome: str
    position_bp: int
    reference_allele: str
    alternate_allele: str
    assembly: str = "GRCh38"
    marker_id: str = ""
    hgvs_c: str = ""
    hgvs_p: str = ""
    transcript: str = ""
    disease: str = ""
    inheritance: str = ""
    variant_type: str = ""
    dbsnp_id: str = ""
    clinvar_id: str = ""
    laboratory_method: str = ""
    source: str = ""
    notes: str = ""
    genetic_distance_cm: str = "0"


@dataclass(frozen=True)
class PedIndividual:
    """Ligne PED validée avant ajout du génotype de mutation."""

    family_id: str
    individual_id: str
    columns: tuple[str, ...]


def write_lines_atomic(output_path: Path, lines: Iterable[str]) -> None:
    """Écrit un fichier complet avant de remplacer atomiquement sa destination."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f".{output_path.name}.tmp")
    try:
        with temporary_path.open("w", encoding="utf-8", newline="\n") as output_file:
            output_file.writelines(lines)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def normalize_chromosome(value: object) -> str:
    """Normalise un chromosome PLINK en supprimant un éventuel préfixe ``chr``."""
    chromosome = str(value).strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    chromosome = chromosome.upper()
    if not chromosome or not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|XY|MT)", chromosome):
        raise MutationInjectionError(f"Chromosome invalide: {value!r}")
    return chromosome


def require_text(config: dict[str, object], field_name: str) -> str:
    """Retourne un champ texte obligatoire et non vide."""
    value = str(config.get(field_name, "")).strip()
    if not value:
        raise MutationInjectionError(
            f"Champ obligatoire absent du fichier mutation: {field_name}"
        )
    return value


def load_mutation_metadata(config_path: Path) -> MutationMetadata:
    """Charge et valide le fichier JSON décrivant la mutation."""
    if not config_path.is_file():
        raise MutationInjectionError(f"Configuration mutation introuvable: {config_path}")
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise MutationInjectionError(
            f"Configuration mutation JSON invalide: {config_path}: {error}"
        ) from error
    if not isinstance(config, dict):
        raise MutationInjectionError("La configuration mutation doit être un objet JSON")

    try:
        position_bp = int(config.get("position_bp", 0))
    except (TypeError, ValueError) as error:
        raise MutationInjectionError("position_bp doit être un entier") from error
    if position_bp <= 0:
        raise MutationInjectionError("position_bp doit être strictement positif")

    reference_allele = require_text(config, "reference_allele").upper()
    alternate_allele = require_text(config, "alternate_allele").upper()
    if reference_allele not in VALID_BASES or alternate_allele not in VALID_BASES:
        raise MutationInjectionError(
            "Les allèles doivent être des bases génomiques simples A, C, G ou T"
        )
    if reference_allele == alternate_allele:
        raise MutationInjectionError("Les allèles de référence et alternatif sont identiques")

    marker_id = str(config.get("marker_id", "")).strip() or f"rsMUT_{position_bp}"
    if not marker_id.startswith("rsMUT") or any(character.isspace() for character in marker_id):
        raise MutationInjectionError(
            "marker_id doit commencer par 'rsMUT' et ne contenir aucun espace"
        )

    optional_fields = {
        field_name: str(config.get(field_name, "")).strip()
        for field_name in (
            "hgvs_c",
            "hgvs_p",
            "transcript",
            "disease",
            "inheritance",
            "variant_type",
            "dbsnp_id",
            "clinvar_id",
            "laboratory_method",
            "source",
            "notes",
        )
    }
    return MutationMetadata(
        gene=require_text(config, "gene"),
        mutation_name=require_text(config, "mutation_name"),
        chromosome=normalize_chromosome(config.get("chromosome", "")),
        position_bp=position_bp,
        reference_allele=reference_allele,
        alternate_allele=alternate_allele,
        assembly=str(config.get("assembly", "GRCh38")).strip() or "GRCh38",
        marker_id=marker_id,
        genetic_distance_cm=str(config.get("genetic_distance_cm", "0")).strip() or "0",
        **optional_fields,
    )


def read_map(map_path: Path) -> list[tuple[str, str, str, int]]:
    """Lit un MAP mono-chromosome et vérifie son ordre physique."""
    if not map_path.is_file():
        raise MutationInjectionError(f"Fichier MAP source introuvable: {map_path}")
    rows: list[tuple[str, str, str, int]] = []
    with map_path.open("r", encoding="utf-8") as map_file:
        for line_number, line in enumerate(map_file, start=1):
            columns = line.split()
            if len(columns) != 4:
                raise MutationInjectionError(
                    f"Ligne MAP invalide {map_path}:{line_number}, 4 colonnes attendues"
                )
            try:
                position_bp = int(columns[3])
            except ValueError as error:
                raise MutationInjectionError(
                    f"Position MAP invalide {map_path}:{line_number}: {columns[3]}"
                ) from error
            rows.append(
                (normalize_chromosome(columns[0]), columns[1], columns[2], position_bp)
            )
    if not rows:
        raise MutationInjectionError(f"Aucun marqueur dans {map_path}")
    duplicate_ids = sorted(
        marker_id
        for marker_id, count in Counter(row[1] for row in rows).items()
        if count > 1
    )
    if duplicate_ids:
        raise MutationInjectionError(
            "Identifiants MAP dupliqués: " + ", ".join(duplicate_ids)
        )
    chromosomes = {row[0] for row in rows}
    if len(chromosomes) != 1:
        raise MutationInjectionError(
            "L'injecteur attend un MAP limité à un seul chromosome"
        )
    positions = [row[3] for row in rows]
    if positions != sorted(positions):
        raise MutationInjectionError("Les positions du MAP source ne sont pas triées")
    return rows


def read_ped(ped_path: Path, marker_count: int) -> list[PedIndividual]:
    """Lit le PED et contrôle le nombre de colonnes et l'unicité FID/IID."""
    if not ped_path.is_file():
        raise MutationInjectionError(f"Fichier PED source introuvable: {ped_path}")
    expected_columns = 6 + 2 * marker_count
    individuals: list[PedIndividual] = []
    with ped_path.open("r", encoding="utf-8") as ped_file:
        for line_number, line in enumerate(ped_file, start=1):
            columns = tuple(line.split())
            if len(columns) != expected_columns:
                raise MutationInjectionError(
                    f"Ligne PED invalide {ped_path}:{line_number}: "
                    f"{len(columns)} colonnes, {expected_columns} attendues"
                )
            individuals.append(PedIndividual(columns[0], columns[1], columns))
    if not individuals:
        raise MutationInjectionError(f"Aucun individu dans {ped_path}")
    duplicate_keys = sorted(
        key
        for key, count in Counter(
            (individual.family_id, individual.individual_id)
            for individual in individuals
        ).items()
        if count > 1
    )
    if duplicate_keys:
        formatted_keys = ", ".join(f"{fid}/{iid}" for fid, iid in duplicate_keys)
        raise MutationInjectionError(f"Individus PED dupliqués: {formatted_keys}")
    return individuals


def read_groups(groups_path: Path | None) -> dict[str, str]:
    """Lit les groupes au format ``IID groupe``."""
    if groups_path is None:
        return {}
    if not groups_path.is_file():
        raise MutationInjectionError(f"Fichier de groupes introuvable: {groups_path}")
    groups: dict[str, str] = {}
    with groups_path.open("r", encoding="utf-8") as groups_file:
        for line_number, line in enumerate(groups_file, start=1):
            columns = line.split()
            if len(columns) != 2:
                raise MutationInjectionError(
                    f"Ligne de groupes invalide {groups_path}:{line_number}"
                )
            individual_id, group = columns
            if individual_id in groups:
                raise MutationInjectionError(
                    f"IID dupliqué dans {groups_path}: {individual_id}"
                )
            groups[individual_id] = group
    return groups


def normalize_genotype(value: str) -> str:
    """Valide une notation abstraite REF/ALT indépendante du brin."""
    genotype = value.strip().upper()
    if genotype == "ALT/REF":
        genotype = "REF/ALT"
    if genotype not in VALID_GENOTYPES:
        raise MutationInjectionError(
            f"Génotype invalide {value!r}; valeurs: " + ", ".join(sorted(VALID_GENOTYPES))
        )
    return genotype


def parse_group_genotypes(values: Iterable[str]) -> dict[str, str]:
    """Transforme les options ``GROUPE=GENOTYPE`` en dictionnaire validé."""
    assignments: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise MutationInjectionError(
                f"Règle de groupe invalide {value!r}; format attendu GROUPE=GENOTYPE"
            )
        group, genotype = (part.strip() for part in value.split("=", maxsplit=1))
        if not group:
            raise MutationInjectionError(f"Groupe vide dans la règle {value!r}")
        if group in assignments:
            raise MutationInjectionError(f"Règle dupliquée pour le groupe {group}")
        assignments[group] = normalize_genotype(genotype)
    return assignments


def read_individual_genotypes(genotypes_path: Path | None) -> dict[tuple[str, str], str]:
    """Lit les génotypes individuels ``FID``, ``IID`` et ``GENOTYPE``."""
    if genotypes_path is None:
        return {}
    if not genotypes_path.is_file():
        raise MutationInjectionError(
            f"Fichier de génotypes individuels introuvable: {genotypes_path}"
        )
    assignments: dict[tuple[str, str], str] = {}
    with genotypes_path.open("r", encoding="utf-8", newline="") as genotype_file:
        reader = csv.DictReader(genotype_file, delimiter="\t")
        required_columns = {"FID", "IID", "GENOTYPE"}
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            raise MutationInjectionError(
                f"Colonnes obligatoires dans {genotypes_path}: FID, IID, GENOTYPE"
            )
        for line_number, row in enumerate(reader, start=2):
            key = (row["FID"].strip(), row["IID"].strip())
            if not all(key):
                raise MutationInjectionError(
                    f"FID ou IID vide dans {genotypes_path}:{line_number}"
                )
            if key in assignments:
                raise MutationInjectionError(
                    f"Individu dupliqué dans {genotypes_path}: {key[0]}/{key[1]}"
                )
            assignments[key] = normalize_genotype(row["GENOTYPE"])
    return assignments


def genotype_to_alleles(
    genotype: str, metadata: MutationMetadata
) -> tuple[str, str]:
    """Convertit REF/ALT vers les allèles génomiques écrits dans le PED."""
    return {
        "REF/REF": (metadata.reference_allele, metadata.reference_allele),
        "REF/ALT": (metadata.reference_allele, metadata.alternate_allele),
        "ALT/ALT": (metadata.alternate_allele, metadata.alternate_allele),
        "MISSING": ("0", "0"),
    }[genotype]


def sha256_file(path: Path) -> str:
    """Calcule l'empreinte SHA-256 d'un fichier sans le charger en mémoire."""
    digest = hashlib.sha256()
    with path.open("rb") as input_file:
        for block in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_outputs_are_available(output_paths: Iterable[Path], force: bool) -> None:
    """Refuse de remplacer une injection existante sans ``--force``."""
    existing_paths = [path for path in output_paths if path.exists()]
    if existing_paths and not force:
        formatted_paths = "\n- ".join(str(path) for path in existing_paths)
        raise MutationInjectionError(
            "Sorties déjà présentes. Utiliser --force pour les remplacer:\n- "
            + formatted_paths
        )


def build_assignments(
    individuals: list[PedIndividual],
    groups: dict[str, str],
    group_genotypes: dict[str, str],
    individual_genotypes: dict[tuple[str, str], str],
) -> list[tuple[PedIndividual, str, str, str]]:
    """Résout chaque génotype avec priorité aux affectations individuelles."""
    source_keys = {(individual.family_id, individual.individual_id) for individual in individuals}
    unknown_keys = sorted(set(individual_genotypes) - source_keys)
    if unknown_keys:
        formatted_keys = ", ".join(f"{fid}/{iid}" for fid, iid in unknown_keys)
        raise MutationInjectionError(
            "Génotypes fournis pour des individus absents du PED: " + formatted_keys
        )

    assignments: list[tuple[PedIndividual, str, str, str]] = []
    unresolved: list[str] = []
    for individual in individuals:
        key = (individual.family_id, individual.individual_id)
        group = groups.get(individual.individual_id, "")
        if key in individual_genotypes:
            genotype = individual_genotypes[key]
            assignment_source = "individual_override"
        elif group and group in group_genotypes:
            genotype = group_genotypes[group]
            assignment_source = f"group:{group}"
        else:
            unresolved.append(f"{individual.family_id}/{individual.individual_id}")
            continue
        assignments.append((individual, group, genotype, assignment_source))
    if unresolved:
        raise MutationInjectionError(
            "Aucun génotype explicite pour: " + ", ".join(unresolved)
        )
    return assignments


def build_project_info(metadata: MutationMetadata, marker_count: int) -> dict[str, object]:
    """Produit les libellés humains actuellement consommés par les rapports."""
    mutation_label_parts = [metadata.mutation_name]
    if metadata.hgvs_c:
        mutation_label_parts.append(metadata.hgvs_c)
    if metadata.hgvs_p:
        mutation_label_parts.append(metadata.hgvs_p)
    return {
        "Maladie": metadata.disease,
        "Gène": metadata.gene,
        "Mutation": " — ".join(mutation_label_parts),
        "Identifiant du marqueur": metadata.marker_id,
        "HGVS c.": metadata.hgvs_c,
        "HGVS p.": metadata.hgvs_p,
        "Transcrit": metadata.transcript,
        "Chromosome": metadata.chromosome,
        "Position génomique (bp)": metadata.position_bp,
        "Assemblage": metadata.assembly,
        "Allèle de référence": metadata.reference_allele,
        "Allèle alternatif": metadata.alternate_allele,
        "Type de variant": metadata.variant_type,
        "Transmission": metadata.inheritance,
        "dbSNP": metadata.dbsnp_id,
        "ClinVar": metadata.clinvar_id,
        "Méthode de laboratoire": metadata.laboratory_method,
        "Source": metadata.source,
        "Notes": metadata.notes,
        "Nombre de marqueurs": marker_count,
    }


def inject_mutation(
    source_prefix: Path,
    output_prefix: Path,
    mutation_config_path: Path,
    groups_path: Path | None = None,
    group_genotype_values: Iterable[str] = (),
    individual_genotypes_path: Path | None = None,
    force: bool = False,
) -> dict[str, object]:
    """Ajoute une mutation et ses génotypes explicites à une copie PED/MAP."""
    source_prefix = source_prefix.resolve()
    output_prefix = output_prefix.resolve()
    if source_prefix.parent == output_prefix.parent:
        raise MutationInjectionError(
            "Le dossier de sortie doit être distinct du dossier source"
        )

    source_ped_path = source_prefix.with_suffix(".ped")
    source_map_path = source_prefix.with_suffix(".map")
    metadata = load_mutation_metadata(mutation_config_path)
    map_rows = read_map(source_map_path)
    if map_rows[0][0] != metadata.chromosome:
        raise MutationInjectionError(
            f"La mutation est sur le chromosome {metadata.chromosome}, "
            f"mais le MAP contient le chromosome {map_rows[0][0]}"
        )
    marker_ids = {row[1] for row in map_rows}
    if metadata.marker_id in marker_ids:
        raise MutationInjectionError(
            f"Le marqueur {metadata.marker_id} existe déjà dans le MAP source"
        )
    individuals = read_ped(source_ped_path, len(map_rows))

    if groups_path is None:
        candidate_groups_path = source_prefix.parent / "groupes.txt"
        groups_path = candidate_groups_path if candidate_groups_path.is_file() else None
    groups = read_groups(groups_path)
    group_genotypes = parse_group_genotypes(group_genotype_values)
    individual_genotypes = read_individual_genotypes(individual_genotypes_path)
    assignments = build_assignments(
        individuals, groups, group_genotypes, individual_genotypes
    )

    insertion_index = next(
        (
            marker_index
            for marker_index, row in enumerate(map_rows)
            if row[3] > metadata.position_bp
        ),
        len(map_rows),
    )
    output_directory = output_prefix.parent
    output_paths = {
        "ped": output_prefix.with_suffix(".ped"),
        "map": output_prefix.with_suffix(".map"),
        "mutation_info": output_directory / "mutation_info.json",
        "project_info": output_directory / "project_info.json",
        "audit": output_directory / "mutation_genotype_audit.tsv",
        "report": output_directory / "mutation_injection_report.json",
    }
    source_sidecars = {
        name: source_prefix.parent / filename
        for name, filename in {
            "groups": "groupes.txt",
            "cases": "cas.txt",
            "controls": "temoins.txt",
        }.items()
        if (source_prefix.parent / filename).is_file()
    }
    output_sidecars = {
        name: output_directory / source_path.name
        for name, source_path in source_sidecars.items()
    }
    ensure_outputs_are_available(
        [*output_paths.values(), *output_sidecars.values()], force
    )

    mutation_row = (
        metadata.chromosome,
        metadata.marker_id,
        metadata.genetic_distance_cm,
        metadata.position_bp,
    )
    output_map_rows = [*map_rows[:insertion_index], mutation_row, *map_rows[insertion_index:]]
    write_lines_atomic(
        output_paths["map"],
        ("\t".join(map(str, row)) + "\n" for row in output_map_rows),
    )

    output_ped_lines: list[str] = []
    audit_lines = ["FID\tIID\tGROUP\tGENOTYPE\tALLELE_1\tALLELE_2\tASSIGNMENT_SOURCE\n"]
    genotype_counts: Counter[str] = Counter()
    for individual, group, genotype, assignment_source in assignments:
        allele_1, allele_2 = genotype_to_alleles(genotype, metadata)
        genotype_offset = 6 + 2 * insertion_index
        output_columns = [
            *individual.columns[:genotype_offset],
            allele_1,
            allele_2,
            *individual.columns[genotype_offset:],
        ]
        output_ped_lines.append("\t".join(output_columns) + "\n")
        audit_lines.append(
            "\t".join(
                [
                    individual.family_id,
                    individual.individual_id,
                    group,
                    genotype,
                    allele_1,
                    allele_2,
                    assignment_source,
                ]
            )
            + "\n"
        )
        genotype_counts[genotype] += 1
    write_lines_atomic(output_paths["ped"], output_ped_lines)
    write_lines_atomic(output_paths["audit"], audit_lines)

    for name, source_path in source_sidecars.items():
        output_path = output_sidecars[name]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, output_path)

    mutation_info = asdict(metadata)
    mutation_info["genotype_notation"] = {
        "REF": metadata.reference_allele,
        "ALT": metadata.alternate_allele,
        "MISSING": "0/0",
    }
    write_lines_atomic(
        output_paths["mutation_info"],
        [json.dumps(mutation_info, indent=2, ensure_ascii=False) + "\n"],
    )
    write_lines_atomic(
        output_paths["project_info"],
        [
            json.dumps(
                build_project_info(metadata, len(output_map_rows)),
                indent=2,
                ensure_ascii=False,
            )
            + "\n"
        ],
    )

    report: dict[str, object] = {
        "source_prefix": str(source_prefix),
        "output_prefix": str(output_prefix),
        "mutation_config": str(mutation_config_path.resolve()),
        "marker_id": metadata.marker_id,
        "chromosome": metadata.chromosome,
        "position_bp": metadata.position_bp,
        "insertion_index_zero_based": insertion_index,
        "source_marker_count": len(map_rows),
        "output_marker_count": len(output_map_rows),
        "individual_count": len(individuals),
        "genotype_counts": dict(sorted(genotype_counts.items())),
        "source_sha256": {
            "ped": sha256_file(source_ped_path),
            "map": sha256_file(source_map_path),
        },
        "output_sha256": {
            "ped": sha256_file(output_paths["ped"]),
            "map": sha256_file(output_paths["map"]),
        },
        "outputs": {
            **{name: str(path) for name, path in output_paths.items()},
            **{name: str(path) for name, path in output_sidecars.items()},
        },
        "warnings": [
            "Les génotypes injectés proviennent des règles explicites fournies; "
            "ils ne sont pas inférés des intensités ACPA.",
            "Les allèles REF/ALT doivent être renseignés sur le brin de référence "
            "de l'assemblage indiqué.",
            "Vérifier la cohérence mendélienne après injection avant toute interprétation.",
        ],
    }
    write_lines_atomic(
        output_paths["report"],
        [json.dumps(report, indent=2, ensure_ascii=False) + "\n"],
    )
    return report


def build_argument_parser() -> argparse.ArgumentParser:
    """Construit l'interface en ligne de commande de l'injecteur."""
    parser = argparse.ArgumentParser(
        description=(
            "Injecte une mutation documentée dans une copie PED/MAP. "
            "Chaque génotype doit être attribué explicitement."
        )
    )
    parser.add_argument("--source-prefix", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--mutation-config", type=Path, required=True)
    parser.add_argument("--groups", type=Path)
    parser.add_argument(
        "--group-genotype",
        action="append",
        default=[],
        metavar="GROUPE=GENOTYPE",
        help="Répétable; génotypes: REF/REF, REF/ALT, ALT/ALT ou MISSING",
    )
    parser.add_argument(
        "--individual-genotypes",
        type=Path,
        help="TSV FID/IID/GENOTYPE; prioritaire sur les règles de groupe",
    )
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    """Lance l'injection et présente une erreur concise à l'utilisateur."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    arguments = build_argument_parser().parse_args()
    try:
        report = inject_mutation(
            source_prefix=arguments.source_prefix,
            output_prefix=arguments.output_prefix,
            mutation_config_path=arguments.mutation_config,
            groups_path=arguments.groups,
            group_genotype_values=arguments.group_genotype,
            individual_genotypes_path=arguments.individual_genotypes,
            force=arguments.force,
        )
        LOGGER.info(
            "Mutation %s injectée pour %s individus; audit: %s",
            report["marker_id"],
            report["individual_count"],
            report["outputs"]["audit"],
        )
        return 0
    except MutationInjectionError as error:
        LOGGER.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
