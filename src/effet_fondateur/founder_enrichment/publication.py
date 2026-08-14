"""Helpers de publication sans logique scientifique supplémentaire."""

from __future__ import annotations

import csv
import gzip
import io
from collections.abc import Iterable, Sequence
from pathlib import Path

from .model import NullDraw


NULL_DRAW_COLUMNS = (
    "NULL_SOURCE", "STRATUM", "DRAW_INDEX", "ATTEMPT_INDEX", "UNIT_COUNT",
    "EVALUATION_STATUS", "LEFT_SHARED_CM", "RIGHT_SHARED_CM", "TOTAL_SHARED_CM",
    "LEFT_MARKER_COUNT", "RIGHT_MARKER_COUNT", "NON_EVALUABLE_REASON",
)


def _tsv_cell(value: object) -> object:
    """Sérialise les valeurs selon les conventions des contrats TSV V2."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    """Écrit une table déterministe avec cellules nulles vides."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = (
        io.TextIOWrapper(gzip.GzipFile(filename=path, mode="wb", mtime=0), encoding="utf-8", newline="")
        if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")
    )
    with handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _tsv_cell(row.get(column)) for column in columns})


def null_draw_rows(draws: Iterable[NullDraw]) -> Iterable[dict[str, object]]:
    """Agrège les métriques sans publier les identifiants des individus tirés."""
    for draw in draws:
        statistic = draw.statistic
        yield {
            "NULL_SOURCE": draw.source, "STRATUM": draw.stratum,
            "DRAW_INDEX": draw.draw_index, "ATTEMPT_INDEX": draw.attempt_index,
            "UNIT_COUNT": len(draw.individual_ids),
            "EVALUATION_STATUS": statistic.evaluation_status,
            "LEFT_SHARED_CM": statistic.left_shared_cm,
            "RIGHT_SHARED_CM": statistic.right_shared_cm,
            "TOTAL_SHARED_CM": statistic.total_shared_cm,
            "LEFT_MARKER_COUNT": statistic.left_marker_count,
            "RIGHT_MARKER_COUNT": statistic.right_marker_count,
            "NON_EVALUABLE_REASON": statistic.non_evaluable_reason,
        }
