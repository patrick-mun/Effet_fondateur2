"""Conversion multi-échantillons d'exports ChAS ACPA vers PLINK PED/MAP."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


LOGGER = logging.getLogger(__name__)

REQUIRED_ACPA_COLUMNS = {
    "Probe Set ID",
    "Forward Strand Base Calls",
    "dbSNP RS ID",
    "Chromosome",
    "Chromosomal Position",
}
REQUIRED_METADATA_COLUMNS = {
    "FILE",
    "FID",
    "IID",
    "PID",
    "MID",
    "SEX",
    "PHENOTYPE",
    "GROUP",
}
VALID_BASES = frozenset("ACGT")
VALID_SEX_VALUES = {"0", "1", "2"}
VALID_PHENOTYPE_VALUES = {"-9", "0", "1", "2"}


class ConversionError(RuntimeError):
    """Erreur de validation empêchant une conversion ACPA fiable."""


@dataclass(frozen=True)
class SampleMetadata:
    file_name: str
    file_path: Path
    family_id: str
    individual_id: str
    paternal_id: str
    maternal_id: str
    sex: str
    phenotype: str
    group: str


@dataclass(frozen=True)
class MarkerCall:
    probe_id: str
    rsid: str
    chromosome: str
    position_bp: int
    allele_1: str
    allele_2: str


@dataclass
class ParsedSample:
    metadata: SampleMetadata
    markers: dict[str, MarkerCall]
    source_metadata: dict[str, str]
    invalid_genotype_count: int
    skipped_row_count: int


@dataclass(frozen=True)
class OutputMarker:
    probe_id: str
    variant_id: str
    rsid: str
    chromosome: str
    position_bp: int


def normalize_chromosome(value: object) -> str:
    """Normalise un chromosome ChAS pour les comparaisons et le fichier MAP."""
    chromosome = str(value).strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    return chromosome.upper()


def parse_genotype(value: object) -> tuple[str, str]:
    """Convertit un appel forward-strand en deux allèles PLINK, ou 0/0."""
    genotype = "".join(str(value).strip().upper().split())
    if len(genotype) != 2 or any(base not in VALID_BASES for base in genotype):
        return "0", "0"
    return genotype[0], genotype[1]


def contains_whitespace(value: str) -> bool:
    """Indique si un identifiant contient un séparateur interdit par PLINK."""
    return any(character.isspace() for character in value)


def read_chas_metadata(file_path: Path) -> dict[str, str]:
    """Lit les métadonnées placées dans les lignes commentées ChAS."""
    metadata: dict[str, str] = {}
    with file_path.open("r", encoding="utf-8-sig", errors="replace") as input_file:
        for raw_line in input_file:
            line = raw_line.strip()
            if not line.startswith("#"):
                break
            key, separator, value = line[1:].partition(":")
            if separator:
                metadata[key.strip()] = value.strip()
    return metadata


def load_sample_metadata(metadata_path: Path, input_dir: Path) -> list[SampleMetadata]:
    """Charge et valide la table reliant fichiers ACPA et individus PLINK."""
    with metadata_path.open("r", encoding="utf-8-sig", newline="") as metadata_file:
        reader = csv.DictReader(metadata_file, delimiter="\t")
        fieldnames = {field.strip() for field in (reader.fieldnames or [])}
        missing_columns = REQUIRED_METADATA_COLUMNS - fieldnames
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ConversionError(f"Colonnes manquantes dans {metadata_path}: {missing}")

        samples: list[SampleMetadata] = []
        for row_number, raw_row in enumerate(reader, start=2):
            row = {str(key).strip(): str(value).strip() for key, value in raw_row.items()}
            file_name = row["FILE"]
            file_path = Path(file_name)
            if not file_path.is_absolute():
                file_path = input_dir / file_path
            file_path = file_path.resolve()

            if not file_path.is_file():
                raise ConversionError(
                    f"Fichier ACPA introuvable à la ligne {row_number}: {file_path}"
                )
            if not row["FID"] or not row["IID"]:
                raise ConversionError(f"FID/IID vide à la ligne {row_number}")
            pedigree_identifiers = {
                "FID": row["FID"],
                "IID": row["IID"],
                "PID": row["PID"] or "0",
                "MID": row["MID"] or "0",
            }
            invalid_identifiers = [
                name
                for name, value in pedigree_identifiers.items()
                if contains_whitespace(value)
            ]
            if invalid_identifiers:
                invalid = ", ".join(invalid_identifiers)
                raise ConversionError(
                    f"Identifiant PLINK contenant un espace à la ligne {row_number}: {invalid}"
                )
            if row["SEX"] not in VALID_SEX_VALUES:
                raise ConversionError(
                    f"SEX invalide à la ligne {row_number}: {row['SEX']}"
                )
            if row["PHENOTYPE"] not in VALID_PHENOTYPE_VALUES:
                raise ConversionError(
                    f"PHENOTYPE invalide à la ligne {row_number}: {row['PHENOTYPE']}"
                )
            if not row["GROUP"]:
                raise ConversionError(f"GROUP vide à la ligne {row_number}")

            samples.append(
                SampleMetadata(
                    file_name=file_name,
                    file_path=file_path,
                    family_id=row["FID"],
                    individual_id=row["IID"],
                    paternal_id=row["PID"] or "0",
                    maternal_id=row["MID"] or "0",
                    sex=row["SEX"],
                    phenotype=row["PHENOTYPE"],
                    group=row["GROUP"],
                )
            )

    if not samples:
        raise ConversionError(f"Aucun échantillon déclaré dans {metadata_path}")

    duplicate_iids = [
        individual_id
        for individual_id, count in Counter(
            sample.individual_id for sample in samples
        ).items()
        if count > 1
    ]
    if duplicate_iids:
        duplicates = ", ".join(sorted(duplicate_iids))
        raise ConversionError(f"IID dupliqué dans les métadonnées: {duplicates}")

    duplicate_files = [
        file_path
        for file_path, count in Counter(sample.file_path for sample in samples).items()
        if count > 1
    ]
    if duplicate_files:
        duplicates = ", ".join(str(path) for path in duplicate_files)
        raise ConversionError(f"Fichier ACPA déclaré plusieurs fois: {duplicates}")

    return samples


def parse_acpa_file(sample: SampleMetadata, chromosome: str) -> ParsedSample:
    """Lit un export ChAS et extrait les appels d'un chromosome."""
    source_metadata = read_chas_metadata(sample.file_path)
    genomic_build = source_metadata.get("UCSC Genomic Version", "")
    if genomic_build and genomic_build.lower() != "hg38":
        raise ConversionError(
            f"Build incompatible pour {sample.file_path.name}: {genomic_build}, hg38 attendu"
        )

    dataframe = pd.read_csv(
        sample.file_path,
        sep="\t",
        comment="#",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
    )
    dataframe.columns = [str(column).strip() for column in dataframe.columns]
    missing_columns = REQUIRED_ACPA_COLUMNS - set(dataframe.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ConversionError(
            f"Colonnes ACPA manquantes dans {sample.file_path.name}: {missing}"
        )

    target_chromosome = normalize_chromosome(chromosome)
    chromosome_mask = dataframe["Chromosome"].map(normalize_chromosome) == target_chromosome
    chromosome_data = dataframe.loc[chromosome_mask]

    markers: dict[str, MarkerCall] = {}
    invalid_genotype_count = 0
    skipped_row_count = 0

    for row_number, row in chromosome_data.iterrows():
        probe_id = str(row["Probe Set ID"]).strip()
        if not probe_id:
            skipped_row_count += 1
            continue
        if contains_whitespace(probe_id):
            raise ConversionError(
                f"Probe Set ID incompatible avec PLINK dans {sample.file_path.name}: {probe_id}"
            )
        if probe_id in markers:
            raise ConversionError(
                f"Probe Set ID dupliqué dans {sample.file_path.name}: {probe_id}"
            )

        try:
            position_bp = int(float(str(row["Chromosomal Position"]).strip()))
        except ValueError:
            LOGGER.warning(
                "Position invalide ignorée dans %s, ligne %s",
                sample.file_path.name,
                row_number + 2,
            )
            skipped_row_count += 1
            continue

        allele_1, allele_2 = parse_genotype(row["Forward Strand Base Calls"])
        if allele_1 == "0":
            invalid_genotype_count += 1

        rsid = str(row["dbSNP RS ID"]).strip()
        if rsid.lower() in {"", ".", "na", "nan", "none"}:
            rsid = ""

        markers[probe_id] = MarkerCall(
            probe_id=probe_id,
            rsid=rsid,
            chromosome=target_chromosome,
            position_bp=position_bp,
            allele_1=allele_1,
            allele_2=allele_2,
        )

    if not markers:
        raise ConversionError(
            f"Aucune sonde du chromosome {target_chromosome} dans {sample.file_path.name}"
        )

    return ParsedSample(
        metadata=sample,
        markers=markers,
        source_metadata=source_metadata,
        invalid_genotype_count=invalid_genotype_count,
        skipped_row_count=skipped_row_count,
    )


def select_candidate_probes(
    parsed_samples: list[ParsedSample], marker_mode: str
) -> set[str]:
    """Sélectionne l'intersection ou l'union des sondes disponibles."""
    marker_sets = [set(sample.markers) for sample in parsed_samples]
    if marker_mode == "intersection":
        return set.intersection(*marker_sets)
    return set.union(*marker_sets)


def build_output_markers(
    parsed_samples: list[ParsedSample],
    candidate_probes: Iterable[str],
    variant_id_source: str,
) -> tuple[list[OutputMarker], list[dict[str, str]]]:
    """Valide les marqueurs communs et construit leur ordre PLINK."""
    accepted: list[OutputMarker] = []
    excluded: list[dict[str, str]] = []

    for probe_id in candidate_probes:
        calls = [
            sample.markers[probe_id]
            for sample in parsed_samples
            if probe_id in sample.markers
        ]
        coordinates = {(call.chromosome, call.position_bp) for call in calls}
        if len(coordinates) != 1:
            excluded.append(
                {
                    "probe_id": probe_id,
                    "reason": "coordonnees_incoherentes",
                    "details": ";".join(
                        f"{chromosome}:{position}"
                        for chromosome, position in sorted(coordinates)
                    ),
                }
            )
            continue

        observed_alleles = {
            allele
            for call in calls
            for allele in (call.allele_1, call.allele_2)
            if allele != "0"
        }
        if len(observed_alleles) > 2:
            excluded.append(
                {
                    "probe_id": probe_id,
                    "reason": "plus_de_deux_alleles",
                    "details": ",".join(sorted(observed_alleles)),
                }
            )
            continue

        chromosome, position_bp = next(iter(coordinates))
        rsids = {call.rsid for call in calls if call.rsid}
        rsid = next(iter(rsids)) if len(rsids) == 1 else ""
        variant_id = rsid if variant_id_source == "rsid-preferred" and rsid else probe_id
        accepted.append(
            OutputMarker(
                probe_id=probe_id,
                variant_id=variant_id,
                rsid=rsid,
                chromosome=chromosome,
                position_bp=position_bp,
            )
        )

    duplicate_variant_ids = {
        variant_id
        for variant_id, count in Counter(marker.variant_id for marker in accepted).items()
        if count > 1
    }
    if duplicate_variant_ids:
        accepted = [
            OutputMarker(
                probe_id=marker.probe_id,
                variant_id=(
                    marker.probe_id
                    if marker.variant_id in duplicate_variant_ids
                    else marker.variant_id
                ),
                rsid=marker.rsid,
                chromosome=marker.chromosome,
                position_bp=marker.position_bp,
            )
            for marker in accepted
        ]

    accepted.sort(
        key=lambda marker: (
            int(marker.chromosome) if marker.chromosome.isdigit() else 10_000,
            marker.position_bp,
            marker.probe_id,
        )
    )
    return accepted, excluded


def ensure_outputs_are_available(output_paths: Iterable[Path], force: bool) -> None:
    """Empêche l'écrasement accidentel de données ou de résultats existants."""
    existing_paths = [path for path in output_paths if path.exists()]
    if existing_paths and not force:
        existing = "\n- ".join(str(path) for path in existing_paths)
        raise ConversionError(
            "Sorties déjà présentes. Utiliser --force pour les remplacer:\n- " + existing
        )


def write_lines_atomic(output_path: Path, lines: Iterable[str]) -> None:
    """Écrit un fichier complet puis le remplace atomiquement."""
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


def write_conversion_outputs(
    parsed_samples: list[ParsedSample],
    output_markers: list[OutputMarker],
    excluded_markers: list[dict[str, str]],
    output_prefix: Path,
    chromosome: str,
    marker_mode: str,
    variant_id_source: str,
    force: bool,
) -> dict[str, object]:
    """Génère PED, MAP, groupes, listes PLINK et rapports de conversion."""
    output_dir = output_prefix.parent
    output_paths = {
        "ped": output_prefix.with_suffix(".ped"),
        "map": output_prefix.with_suffix(".map"),
        "groups": output_dir / "groupes.txt",
        "cases": output_dir / "cas.txt",
        "controls": output_dir / "temoins.txt",
        "report": output_dir / "acpa_conversion_report.json",
        "excluded": output_dir / "acpa_excluded_markers.tsv",
        "marker_audit": output_dir / "acpa_marker_audit.tsv",
    }
    ensure_outputs_are_available(output_paths.values(), force)

    map_lines = (
        f"{marker.chromosome}\t{marker.variant_id}\t0\t{marker.position_bp}\n"
        for marker in output_markers
    )
    write_lines_atomic(output_paths["map"], map_lines)

    ped_lines: list[str] = []
    for sample in parsed_samples:
        pedigree_columns = [
            sample.metadata.family_id,
            sample.metadata.individual_id,
            sample.metadata.paternal_id,
            sample.metadata.maternal_id,
            sample.metadata.sex,
            sample.metadata.phenotype,
        ]
        genotype_columns: list[str] = []
        for marker in output_markers:
            call = sample.markers.get(marker.probe_id)
            if call is None:
                genotype_columns.extend(["0", "0"])
            else:
                genotype_columns.extend([call.allele_1, call.allele_2])
        ped_lines.append("\t".join(pedigree_columns + genotype_columns) + "\n")
    write_lines_atomic(output_paths["ped"], ped_lines)

    write_lines_atomic(
        output_paths["groups"],
        (
            f"{sample.metadata.individual_id}\t{sample.metadata.group}\n"
            for sample in parsed_samples
        ),
    )
    write_lines_atomic(
        output_paths["cases"],
        (
            f"{sample.metadata.family_id}\t{sample.metadata.individual_id}\n"
            for sample in parsed_samples
            if sample.metadata.group.upper() == "ATTEINT"
        ),
    )
    write_lines_atomic(
        output_paths["controls"],
        (
            f"{sample.metadata.family_id}\t{sample.metadata.individual_id}\n"
            for sample in parsed_samples
            if sample.metadata.group.upper() == "TEMOIN"
        ),
    )

    excluded_lines = ["probe_id\treason\tdetails\n"]
    excluded_lines.extend(
        f"{row['probe_id']}\t{row['reason']}\t{row['details']}\n"
        for row in excluded_markers
    )
    write_lines_atomic(output_paths["excluded"], excluded_lines)

    audit_lines = ["probe_id\tvariant_id\trsid\tchromosome\tposition_bp\n"]
    audit_lines.extend(
        "\t".join(
            [
                marker.probe_id,
                marker.variant_id,
                marker.rsid,
                marker.chromosome,
                str(marker.position_bp),
            ]
        )
        + "\n"
        for marker in output_markers
    )
    write_lines_atomic(output_paths["marker_audit"], audit_lines)

    builds = sorted(
        {
            sample.source_metadata.get("UCSC Genomic Version", "non_renseigne")
            for sample in parsed_samples
        }
    )
    sample_reports = [
        {
            "file": sample.metadata.file_name,
            "fid": sample.metadata.family_id,
            "iid": sample.metadata.individual_id,
            "group": sample.metadata.group,
            "marker_count": len(sample.markers),
            "invalid_genotype_count": sample.invalid_genotype_count,
            "skipped_row_count": sample.skipped_row_count,
        }
        for sample in parsed_samples
    ]
    report: dict[str, object] = {
        "chromosome": normalize_chromosome(chromosome),
        "genomic_builds": builds,
        "sample_count": len(parsed_samples),
        "marker_mode": marker_mode,
        "variant_id_source": variant_id_source,
        "output_marker_count": len(output_markers),
        "excluded_marker_count": len(excluded_markers),
        "samples": sample_reports,
        "outputs": {name: str(path) for name, path in output_paths.items()},
    }
    write_lines_atomic(
        output_paths["report"],
        [json.dumps(report, indent=2, ensure_ascii=False) + "\n"],
    )
    return report


def convert_acpa_to_plink(
    input_dir: Path,
    metadata_path: Path,
    output_prefix: Path,
    chromosome: str = "19",
    marker_mode: str = "intersection",
    variant_id_source: str = "probe",
    force: bool = False,
) -> dict[str, object]:
    """Convertit plusieurs exports ACPA en un jeu PLINK cohérent."""
    samples = load_sample_metadata(metadata_path.resolve(), input_dir.resolve())
    parsed_samples = [parse_acpa_file(sample, chromosome) for sample in samples]
    candidate_probes = select_candidate_probes(parsed_samples, marker_mode)
    output_markers, excluded_markers = build_output_markers(
        parsed_samples,
        candidate_probes,
        variant_id_source,
    )
    if not output_markers:
        raise ConversionError("Aucun marqueur valide disponible après les contrôles QC")

    return write_conversion_outputs(
        parsed_samples=parsed_samples,
        output_markers=output_markers,
        excluded_markers=excluded_markers,
        output_prefix=output_prefix.resolve(),
        chromosome=chromosome,
        marker_mode=marker_mode,
        variant_id_source=variant_id_source,
        force=force,
    )


def create_metadata_template(
    input_dir: Path, template_path: Path, force: bool = False
) -> None:
    """Crée une table à compléter à partir des noms de fichiers ACPA."""
    acpa_files = sorted(input_dir.glob("*.txt"))
    if not acpa_files:
        raise ConversionError(f"Aucun fichier .txt trouvé dans {input_dir}")
    ensure_outputs_are_available([template_path], force)

    lines = ["FILE\tFID\tIID\tPID\tMID\tSEX\tPHENOTYPE\tGROUP\n"]
    for acpa_file in acpa_files:
        individual_id = acpa_file.name.split("_(", maxsplit=1)[0]
        individual_id = re.sub(r"\s+", "_", individual_id).strip("_")
        lines.append(
            "\t".join(
                [
                    acpa_file.name,
                    individual_id,
                    individual_id,
                    "0",
                    "0",
                    "0",
                    "-9",
                    "NON_RENSEIGNE",
                ]
            )
            + "\n"
        )
    write_lines_atomic(template_path.resolve(), lines)


def build_argument_parser() -> argparse.ArgumentParser:
    """Construit l'interface en ligne de commande du convertisseur."""
    parser = argparse.ArgumentParser(
        description="Convertit des exports ChAS ACPA multi-échantillons en PED/MAP PLINK."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--output-prefix", type=Path)
    parser.add_argument("--chromosome", default="19")
    parser.add_argument(
        "--marker-mode",
        choices=["intersection", "union"],
        default="intersection",
    )
    parser.add_argument(
        "--variant-id-source",
        choices=["probe", "rsid-preferred"],
        default="probe",
    )
    parser.add_argument("--create-metadata-template", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    """Exécute la création du modèle ou la conversion demandée."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = build_argument_parser()
    arguments = parser.parse_args()

    try:
        if arguments.create_metadata_template:
            create_metadata_template(
                input_dir=arguments.input_dir,
                template_path=arguments.create_metadata_template,
                force=arguments.force,
            )
            LOGGER.info("Modèle créé: %s", arguments.create_metadata_template)
            return 0

        if arguments.metadata is None or arguments.output_prefix is None:
            parser.error(
                "--metadata et --output-prefix sont obligatoires pour une conversion"
            )

        report = convert_acpa_to_plink(
            input_dir=arguments.input_dir,
            metadata_path=arguments.metadata,
            output_prefix=arguments.output_prefix,
            chromosome=arguments.chromosome,
            marker_mode=arguments.marker_mode,
            variant_id_source=arguments.variant_id_source,
            force=arguments.force,
        )
        LOGGER.info(
            "Conversion terminée: %s individus, %s marqueurs",
            report["sample_count"],
            report["output_marker_count"],
        )
        return 0
    except ConversionError as error:
        LOGGER.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
