"""Normalisation contrôlée du paquet constructeur CytoScan D11b.

Le paquet officiel reste hors du dépôt et hors du cache du pipeline. Seule
sa base d'annotations NA36 est ouverte, en lecture seule, pour produire dans le
workspace privé une table de SNP autosomiques indépendante de tout échantillon.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile
import zipfile

from effet_fondateur.calibration.simulation.chip_markers import (
    CHIP_ASSEMBLY,
    CHIP_DESIGN_ID,
    CHIP_DESIGN_VERSION,
    ChipMarkerError,
)


ANNOTATION_DATABASE_MEMBER = "CytoScan750k_Accel_Array.na36.annot.db"
NORMALIZED_MARKER_FILENAME = "CytoScan750K_Accel_NA36.r1_GRCh38_SNP.tsv"
_EXPECTED_ARRAY_TYPE = "CytoScan750k_Accel"
_EXPECTED_NETAFFX_BUILD = "36"
_MAX_ANNOTATION_DATABASE_BYTES = 2_000_000_000
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_STABLE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
_BASES = frozenset("ACGT")
_FORBIDDEN_SOURCE_PARTS = {
    "dock6_real",
    "res_acpa",
    "temoins",
    "snp_acpa_accel",
}
_ANNOTATION_COLUMNS = (
    "ProbeSet_ID",
    "Chr_id",
    "Start",
    "Stop",
    "Strand",
    "dbSNP_RS_ID",
    "Strand_Vs_dbSNP",
    "Probe_Count",
    "Cytoband",
    "ChrX_PAR",
    "Flank",
    "Allele_A",
    "Allele_B",
    "Associated_Gene",
    "Genetic_Map",
    "Microsatellite",
    "Enzyme_Fragment",
    "CNV_Record",
    "SNP_Interference",
    "Percent_GC",
    "Heterozygous_Allele_Frequencies",
    "Allele_Frequency_Count",
    "Allele_Frequencies",
    "Minor_Allele",
    "Minor_Allele_Frequency",
    "OMIM",
    "ProbeSet_Type",
    "Probe_GC_A",
    "Process_Flag",
    "Fragment_GC",
)
_NORMALIZED_COLUMNS = (
    "MARKER_ID",
    "ARRAY_DESIGN",
    "DESIGN_VERSION",
    "ASSEMBLY",
    "CHROMOSOME",
    "POSITION_BP",
    "ALLELE_A",
    "ALLELE_B",
    "STRAND",
)


@dataclass(frozen=True)
class ChipPackageNormalization:
    """Résultat agrégé de la normalisation d'un paquet officiel."""

    package_sha256: str
    normalized_marker_sha256: str
    annotation_row_count: int
    snp_marker_count: int
    autosomal_marker_count: int
    non_autosomal_snp_count: int


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_package_path(path: Path, expected_sha256: str) -> str:
    if _SHA256.fullmatch(expected_sha256) is None:
        raise ChipMarkerError("chip_official_package_expected_sha256_invalid")
    if not path.is_file() or path.is_symlink() or path.suffix.casefold() != ".zip":
        raise ChipMarkerError("chip_official_package_invalid")
    resolved_parts = {part.casefold() for part in path.resolve().parts}
    if _FORBIDDEN_SOURCE_PARTS & resolved_parts or {"data", "input"} <= resolved_parts:
        raise ChipMarkerError("chip_official_package_source_forbidden")
    current = path.parent
    while current != current.parent:
        if current.is_symlink():
            raise ChipMarkerError("chip_official_package_symlink_forbidden")
        current = current.parent
    package_sha256 = _sha256_file(path)
    if package_sha256 != expected_sha256:
        raise ChipMarkerError("chip_official_package_sha256_mismatch")
    return package_sha256


def _annotation_member(archive: zipfile.ZipFile) -> zipfile.ZipInfo:
    matches = [
        member
        for member in archive.infolist()
        if member.filename == ANNOTATION_DATABASE_MEMBER
    ]
    if len(matches) != 1:
        raise ChipMarkerError("chip_official_annotation_database_missing")
    member = matches[0]
    if (
        member.is_dir()
        or member.file_size < 1
        or member.file_size > _MAX_ANNOTATION_DATABASE_BYTES
    ):
        raise ChipMarkerError("chip_official_annotation_database_invalid")
    return member


def _extract_annotation_database(
    package_path: Path,
    destination: Path,
) -> None:
    try:
        with zipfile.ZipFile(package_path) as archive:
            member = _annotation_member(archive)
            with archive.open(member) as source, destination.open("wb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
    except (OSError, zipfile.BadZipFile, RuntimeError) as error:
        raise ChipMarkerError("chip_official_package_unreadable") from error
    if destination.stat().st_size != member.file_size:
        raise ChipMarkerError("chip_official_annotation_database_truncated")


def _database_connection(path: Path) -> sqlite3.Connection:
    try:
        connection = sqlite3.connect(f"{path.as_uri()}?mode=ro&immutable=1", uri=True)
        connection.execute("PRAGMA query_only = ON")
        return connection
    except sqlite3.Error as error:
        raise ChipMarkerError("chip_official_annotation_database_unreadable") from error


def _validate_database(connection: sqlite3.Connection) -> tuple[int, int]:
    try:
        integrity = connection.execute("PRAGMA quick_check").fetchone()
        if integrity != ("ok",):
            raise ChipMarkerError("chip_official_annotation_database_integrity_failed")
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        if not {"Annotations", "Information"} <= tables:
            raise ChipMarkerError("chip_official_annotation_database_schema_invalid")
        columns = tuple(
            row[1] for row in connection.execute("PRAGMA table_info(Annotations)")
        )
        if columns != _ANNOTATION_COLUMNS:
            raise ChipMarkerError("chip_official_annotation_database_schema_invalid")
        information_rows = connection.execute(
            "SELECT key, value FROM Information"
        ).fetchall()
    except sqlite3.Error as error:
        raise ChipMarkerError("chip_official_annotation_database_unreadable") from error
    information: dict[str, str] = {}
    for key, value in information_rows:
        if not isinstance(key, str) or not isinstance(value, str) or key in information:
            raise ChipMarkerError("chip_official_annotation_metadata_invalid")
        information[key] = value
    if (
        information.get("array_type") != _EXPECTED_ARRAY_TYPE
        or information.get("genome-version") != "hg38"
        or information.get("genome-version-ncbi") != "38"
        or information.get("genome-version-ensembl-url-parameter") != CHIP_ASSEMBLY
        or information.get("netaffx-annotation-netaffx-build")
        != _EXPECTED_NETAFFX_BUILD
    ):
        raise ChipMarkerError("chip_official_annotation_metadata_mismatch")
    try:
        declared_annotation_count = int(information["annotation_row_count"])
    except (KeyError, TypeError, ValueError) as error:
        raise ChipMarkerError("chip_official_annotation_metadata_invalid") from error
    annotation_count = connection.execute(
        "SELECT COUNT(*) FROM Annotations"
    ).fetchone()[0]
    snp_count = connection.execute(
        "SELECT COUNT(*) FROM Annotations WHERE ProbeSet_Type = 'SNP'"
    ).fetchone()[0]
    if (
        annotation_count != declared_annotation_count
        or type(annotation_count) is not int
        or type(snp_count) is not int
        or snp_count < 1
    ):
        raise ChipMarkerError("chip_official_annotation_count_mismatch")
    return annotation_count, snp_count


def _write_normalized_markers(
    connection: sqlite3.Connection,
    output_path: Path,
) -> int:
    query = """
        SELECT ProbeSet_ID, Chr_id, Start, Stop, Allele_A, Allele_B, Strand
        FROM Annotations
        WHERE ProbeSet_Type = 'SNP' AND Chr_id BETWEEN 1 AND 22
        ORDER BY Chr_id, Start, ProbeSet_ID
    """
    marker_count = 0
    identifiers: set[str] = set()
    try:
        rows = connection.execute(query)
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(_NORMALIZED_COLUMNS)
            for marker_id, chromosome, start, stop, allele_a, allele_b, strand in rows:
                if (
                    not isinstance(marker_id, str)
                    or _STABLE_IDENTIFIER.fullmatch(marker_id) is None
                    or type(chromosome) is not int
                    or type(start) is not int
                    or type(stop) is not int
                    or not 1 <= chromosome <= 22
                    or start < 1
                    or start != stop
                    or allele_a not in _BASES
                    or allele_b not in _BASES
                    or allele_a == allele_b
                    or strand not in {"+", "-"}
                ):
                    raise ChipMarkerError("chip_official_autosomal_snp_invalid")
                identifier_key = marker_id.casefold()
                if identifier_key in identifiers:
                    raise ChipMarkerError("chip_official_marker_identifier_duplicate")
                identifiers.add(identifier_key)
                writer.writerow(
                    (
                        marker_id,
                        CHIP_DESIGN_ID,
                        CHIP_DESIGN_VERSION,
                        CHIP_ASSEMBLY,
                        chromosome,
                        start,
                        allele_a,
                        allele_b,
                        strand,
                    )
                )
                marker_count += 1
    except sqlite3.Error as error:
        raise ChipMarkerError("chip_official_annotation_database_unreadable") from error
    if marker_count < 1:
        raise ChipMarkerError("chip_official_autosomal_marker_empty")
    return marker_count


def normalize_official_chip_package(
    *,
    package_path: Path,
    expected_package_sha256: str,
    output_path: Path,
) -> ChipPackageNormalization:
    """Produit atomiquement le TSV SNP autosomique depuis le ZIP officiel.

    La liste `snplist.r3.txt` du paquet est réservée au calcul du SNPQC et ne
    définit donc pas l'univers observable. La table `Annotations`, spécifique
    au dessin Accel, est la source des coordonnées, allèles et brins.
    """
    package_sha256 = _validate_package_path(package_path, expected_package_sha256)
    if output_path.exists() or output_path.is_symlink():
        raise ChipMarkerError("chip_normalized_marker_output_exists")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".cytoscan_na36_normalization.",
        dir=output_path.parent,
    ) as temporary_name:
        temporary_root = Path(temporary_name)
        database_path = temporary_root / ANNOTATION_DATABASE_MEMBER
        normalized_path = temporary_root / NORMALIZED_MARKER_FILENAME
        _extract_annotation_database(package_path, database_path)
        connection = _database_connection(database_path)
        try:
            annotation_count, snp_count = _validate_database(connection)
            autosomal_count = _write_normalized_markers(connection, normalized_path)
        finally:
            connection.close()
        normalized_sha256 = _sha256_file(normalized_path)
        normalized_path.replace(output_path)
    return ChipPackageNormalization(
        package_sha256=package_sha256,
        normalized_marker_sha256=normalized_sha256,
        annotation_row_count=annotation_count,
        snp_marker_count=snp_count,
        autosomal_marker_count=autosomal_count,
        non_autosomal_snp_count=snp_count - autosomal_count,
    )
