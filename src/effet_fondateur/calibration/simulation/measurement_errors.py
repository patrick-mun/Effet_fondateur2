"""D17 : perturbations appariées de l'observation, jamais de la vérité biologique.

Les quotas sont arrondis à l'entier le plus proche (demi supérieur), puis les
positions sont tirées uniformément sans remise. Ce modèle de résistance MCAR
n'est pas une estimation des erreurs d'une puce. Une commutation inverse
l'orientation des hétérozygotes depuis une frontière jusqu'à la suivante ;
une erreur ponctuelle n'inverse qu'un site. Les deux peuvent se compenser.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, ROUND_HALF_UP
import hashlib
import random
from typing import Any, TextIO

import numpy as np
from numpy.typing import NDArray


METHOD_ID = "paired_measurement_perturbations_v1"
EVENT_TYPES = ("MARKER_DROP", "GENOTYPE", "PHASE_SWITCH_BOUNDARY",
               "PHASE_SWITCH_GT", "PHASE_POINT", "MISSING_GT")
GT_VALUES = ("0|0", "0|1", "1|0", "1|1", ".|.")


class MeasurementError(ValueError):
    """Violation D17, sans génotype ni identifiant individuel dans le message."""


@dataclass(frozen=True)
class MeasurementPanel:
    """Sondes ordonnées et copies H1/H2 ; -1 représente un appel absent complet."""

    samples: tuple[str, ...]
    variants: tuple[tuple[str, ...], ...]
    alleles: NDArray[np.int8]


@dataclass(frozen=True)
class PerturbedPanel:
    """Observation et registre chronologique privé ; aucun objet source muté."""

    alleles: NDArray[np.int8]
    retained: NDArray[np.bool_]
    events: tuple[dict[str, Any], ...]
    summary: dict[str, Any]


def quota(rate: float, population: int) -> int:
    """Nombre exact d'événements ; taux sans unité, population entière positive."""
    value = Decimal(str(rate))
    if not value.is_finite() or not 0 <= value <= 1 or population < 0:
        raise MeasurementError("d17_rate_or_population_invalid")
    return int((value * population).quantize(Decimal("1"), rounding=ROUND_HALF_UP))


def derive_seed(parent_seed: int, role: str, profile_id: str, operation: str) -> int:
    """Sous-flux SHA-256 distincts par rôle/profil/opération, pas par position de liste."""
    if type(parent_seed) is not int or not 0 <= parent_seed <= 4294967295:
        raise MeasurementError("d17_seed_invalid")
    if role not in {"CALIBRATION", "VALIDATION"}:
        raise MeasurementError("d17_panel_role_invalid")
    payload = f"{METHOD_ID}:D17:{parent_seed}:{role}:{profile_id}:{operation}"
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:4], "big")


def read_panel(lines: Iterable[str], *, phased: bool) -> MeasurementPanel:
    """Lit un VCF de sondes SNV bialléliques, GT seul ; ordre strict, sans doublon.

    Les génotypes partiels, multialléliques, champs annexes et identifiants
    dupliqués sont refusés. Seuls ././.|. complets sont admis comme manque.
    """
    samples: tuple[str, ...] = ()
    variants: list[tuple[str, ...]] = []
    values: list[list[tuple[int, int]]] = []
    previous = (0, 0)
    marker_ids: set[str] = set()
    separator = "|" if phased else "/"
    allowed = {f"{left}{separator}{right}": (left, right)
               for left in (0, 1) for right in (0, 1)}
    allowed[f".{separator}."] = (-1, -1)
    for line in lines:
        if line.startswith("##"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if line.startswith("#CHROM\t"):
            if samples or len(fields) < 10 or len(set(fields[9:])) != len(fields[9:]) or any(not name for name in fields[9:]):
                raise MeasurementError("d17_sample_header_invalid")
            samples = tuple(fields[9:])
            continue
        if not samples or len(fields) != 9 + len(samples) or fields[8] != "GT":
            raise MeasurementError("d17_vcf_columns_invalid")
        try:
            chromosome = int(fields[0].removeprefix("chr"))
            position = int(fields[1])
        except ValueError as error:
            raise MeasurementError("d17_vcf_coordinates_invalid") from error
        if not 1 <= chromosome <= 22 or position <= 0 or (chromosome, position) <= previous:
            raise MeasurementError("d17_vcf_order_invalid")
        if fields[3] not in "ACGT" or len(fields[3]) != 1 or fields[4] not in "ACGT" or len(fields[4]) != 1 or fields[3] == fields[4]:
            raise MeasurementError("d17_biallelic_snv_required")
        if not fields[2] or fields[2] == "." or fields[2] in marker_ids:
            raise MeasurementError("d17_marker_id_invalid")
        try:
            values.append([allowed[value] for value in fields[9:]])
        except KeyError as error:
            raise MeasurementError("d17_gt_invalid") from error
        previous = chromosome, position
        marker_ids.add(fields[2])
        variants.append(tuple(fields[:9]))
    if not variants:
        raise MeasurementError("d17_empty_panel")
    return MeasurementPanel(samples, tuple(variants), np.array(values, dtype=np.int8))


def project_dense_truth(
    dense_lines: Iterable[str], observed: MeasurementPanel,
) -> MeasurementPanel:
    """Projette la vérité dense sur les sondes sans perdre l'orientation H1/H2.

    Les appels non phasés D16 doivent correspondre exactement à leur dosage
    dense. Les autres variants sont ignorés après identification, jamais imputés.
    """
    wanted = {(row[0], row[1], row[3], row[4]): index
              for index, row in enumerate(observed.variants)}
    found: dict[int, list[str]] = {}
    header: str | None = None
    for line in dense_lines:
        if line.startswith("#CHROM\t"):
            if header is not None or tuple(line.rstrip().split("\t")[9:]) != observed.samples:
                raise MeasurementError("d17_dense_sample_order_mismatch")
            header = line
        elif not line.startswith("#"):
            fields = line.rstrip("\r\n").split("\t", 9)
            if header is None or len(fields) != 10:
                raise MeasurementError("d17_dense_header_or_columns_invalid")
            index = wanted.get((fields[0], fields[1], fields[3], fields[4]))
            if index is not None:
                if index in found or fields[8] != "GT":
                    raise MeasurementError("d17_dense_marker_duplicate_or_format_invalid")
                found[index] = fields[9].split("\t")
    if len(found) != len(wanted):
        raise MeasurementError("d17_dense_projection_incomplete")
    lines = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(observed.samples) + "\n"]
    lines.extend("\t".join(row) + "\t" + "\t".join(found[index]) + "\n"
                 for index, row in enumerate(observed.variants))
    truth = read_panel(lines, phased=True)
    if np.any(truth.alleles < 0) or not np.array_equal(np.sort(truth.alleles, axis=2), np.sort(observed.alleles, axis=2)):
        raise MeasurementError("d17_dense_observed_genotype_mismatch")
    return truth


def _gt(pair: Sequence[int]) -> str:
    return ".|." if pair[0] < 0 else f"{pair[0]}|{pair[1]}"


def inject_profile(
    panel: MeasurementPanel, profile: Mapping[str, Any], *, parent_seed: int,
    panel_role: str, target_chromosome: int, target_bp: int,
) -> PerturbedPanel:
    """Crée un profil indépendant depuis les copies vraies ; unités bp et taux.

    Ordre : réduction de densité, erreurs de dosage, commutations, erreurs de
    phase ponctuelles, manque. Les taux portent sur les cellules retenues,
    sauf la phase : hétérozygotes pour les points, frontières entre deux
    hétérozygotes successifs d'un même échantillon/chromosome pour les switches.
    Chaque changement effectif et chaque frontière sont enregistrés avant/après.
    """
    if np.any(panel.alleles < 0):
        raise MeasurementError("d17_complete_truth_required")
    profile_id = str(profile["profile_id"])
    seeds = {operation: derive_seed(parent_seed, panel_role, profile_id, operation)
             for operation in ("density", "genotype", "phase_switch", "phase_point", "missing")}
    generators = {name: random.Random(seed) for name, seed in seeds.items()}
    alleles = panel.alleles.copy()
    marker_count, sample_count, _ = alleles.shape
    retained = np.ones(marker_count, dtype=bool)
    events: list[dict[str, Any]] = []

    def record(kind: str, marker: int, sample: int | None, before: str, after: str) -> None:
        row = panel.variants[marker]
        field = "VARIANT_PRESENCE" if kind == "MARKER_DROP" else "PHASE_ORIENTATION" if kind == "PHASE_SWITCH_BOUNDARY" else "GT"
        events.append(dict(error_event_id=f"{profile_id}.D17_ERROR_{len(events) + 1:09d}", error_profile_id=profile_id,
                           sample_id=panel.samples[sample] if sample is not None else None,
                           error_type=kind, affected_field=field, original_value=before,
                           injected_value=after, chromosome=int(row[0].removeprefix("chr")),
                           position_bp=int(row[1]), ref=row[3], alt=row[4], marker_id=row[2], injected=True))

    protected = {index for index, row in enumerate(panel.variants)
                 if int(row[0].removeprefix("chr")) == target_chromosome and int(row[1]) == target_bp}
    keep_count = max(len(protected), quota(profile["marker_retention_fraction"], marker_count))
    if keep_count < 1:
        raise MeasurementError("d17_empty_retained_panel")
    candidates = sorted(set(range(marker_count)) - protected)
    for marker in sorted(generators["density"].sample(candidates, marker_count - keep_count)):
        retained[marker] = False
        record("MARKER_DROP", marker, None, "RETAINED", "REMOVED")
    cells = [int(marker) * sample_count + sample for marker in np.flatnonzero(retained)
             for sample in range(sample_count)]
    genotype_count = quota(profile["genotyping_error_rate"], len(cells))
    for cell in sorted(generators["genotype"].sample(cells, genotype_count)):
        marker, sample = divmod(cell, sample_count)
        pair = alleles[marker, sample]
        before = _gt(pair)
        dosage = generators["genotype"].choice([value for value in (0, 1, 2) if value != int(pair.sum())])
        pair[:] = (0, 0) if dosage == 0 else (1, 1) if dosage == 2 else generators["genotype"].choice(((0, 1), (1, 0)))
        record("GENOTYPE", marker, sample, before, _gt(pair))
    heterogeneous = retained[:, None] & (alleles[:, :, 0] != alleles[:, :, 1])
    tracks: list[tuple[int, list[int]]] = []
    boundaries: list[tuple[int, int]] = []
    chromosomes = tuple(dict.fromkeys(row[0] for row in panel.variants))
    for sample in range(sample_count):
        for chromosome in chromosomes:
            markers = [int(marker) for marker in np.flatnonzero(heterogeneous[:, sample])
                       if panel.variants[marker][0] == chromosome]
            tracks.append((sample, markers))
            boundaries.extend((sample, marker) for marker in markers[1:])
    switch_count = quota(profile["phase_switch_rate"], len(boundaries))
    selected_boundaries = set(generators["phase_switch"].sample(boundaries, switch_count))
    for sample, markers in tracks:
        state = 0
        for marker in markers:
            if (sample, marker) in selected_boundaries:
                record("PHASE_SWITCH_BOUNDARY", marker, sample, str(state), str(1 - state))
                state = 1 - state
            if state:
                before = _gt(alleles[marker, sample])
                alleles[marker, sample] = alleles[marker, sample, ::-1]
                record("PHASE_SWITCH_GT", marker, sample, before, _gt(alleles[marker, sample]))
    hetero_cells = np.flatnonzero(heterogeneous.ravel()).tolist()
    point_count = quota(profile["phase_point_rate"], len(hetero_cells))
    for cell in sorted(generators["phase_point"].sample(hetero_cells, point_count)):
        marker, sample = divmod(cell, sample_count)
        before = _gt(alleles[marker, sample])
        alleles[marker, sample] = alleles[marker, sample, ::-1]
        record("PHASE_POINT", marker, sample, before, _gt(alleles[marker, sample]))
    missing_count = quota(profile["missing_rate"], len(cells))
    for cell in sorted(generators["missing"].sample(cells, missing_count)):
        marker, sample = divmod(cell, sample_count)
        before = _gt(alleles[marker, sample])
        alleles[marker, sample] = (-1, -1)
        record("MISSING_GT", marker, sample, before, ".|.")
    counts = Counter(event["error_type"] for event in events)
    summary = dict(profile_id=profile_id, seeds=seeds, requested=dict(profile),
                   source_marker_count=marker_count, retained_marker_count=int(retained.sum()),
                   sample_count=sample_count, retained_cell_count=len(cells),
                   heterozygote_count=int(heterogeneous.sum()), eligible_switch_boundary_count=len(boundaries),
                   event_counts={kind: counts[kind] for kind in EVENT_TYPES}, event_count=len(events),
                   target_observed=bool(protected), target_retained=all(retained[index] for index in protected))
    return PerturbedPanel(alleles, retained, tuple(events), summary)


def replay_events(
    panel: MeasurementPanel, events: Iterable[Mapping[str, Any]], *, profile_id: str,
) -> tuple[NDArray[np.int8], NDArray[np.bool_]]:
    """Rejoue le registre seul : avant/après, type, ordre et frontières contrôlés.

    Cette seconde voie ne rappelle pas l'injecteur ni ses tirages. Une erreur
    d'export non enregistrée sera détectée en comparant ses sorties au rejeu.
    """
    alleles = panel.alleles.copy()
    retained = np.ones(len(panel.variants), dtype=bool)
    genotype_alleles = panel.alleles.copy()
    markers = {(int(row[0].removeprefix("chr")), int(row[1]), row[3], row[4], row[2]): index
               for index, row in enumerate(panel.variants)}
    samples = {sample: index for index, sample in enumerate(panel.samples)}
    orientation: dict[tuple[int, int], int] = {}
    seen: set[tuple[str, int, int | None]] = set()
    boundaries: set[tuple[int, int]] = set()
    switched: set[tuple[int, int]] = set()
    last_rank = -1
    for ordinal, event in enumerate(events, 1):
        kind = event["error_type"]
        if kind not in EVENT_TYPES or event["error_event_id"] != f"{profile_id}.D17_ERROR_{ordinal:09d}" or event["error_profile_id"] != profile_id or event["injected"] is not True:
            raise MeasurementError("d17_event_identity_or_type_invalid")
        rank = {"MARKER_DROP": 0, "GENOTYPE": 1, "PHASE_SWITCH_BOUNDARY": 2,
                "PHASE_SWITCH_GT": 2, "PHASE_POINT": 3, "MISSING_GT": 4}[kind]
        if rank < last_rank:
            raise MeasurementError("d17_event_order_invalid")
        last_rank = rank
        key = (event["chromosome"], event["position_bp"], event["ref"], event["alt"], event["marker_id"])
        if key not in markers or (event["sample_id"] is not None and event["sample_id"] not in samples):
            raise MeasurementError("d17_event_locus_or_sample_invalid")
        marker = markers[key]
        sample = samples.get(event["sample_id"])
        unique = (kind, marker, sample)
        if unique in seen:
            raise MeasurementError("d17_event_duplicate")
        seen.add(unique)
        before, after = event["original_value"], event["injected_value"]
        if kind == "MARKER_DROP":
            if sample is not None or before != "RETAINED" or after != "REMOVED" or event["affected_field"] != "VARIANT_PRESENCE":
                raise MeasurementError("d17_drop_event_invalid")
            retained[marker] = False
            continue
        if sample is None or not retained[marker]:
            raise MeasurementError("d17_event_on_removed_marker")
        if kind == "PHASE_SWITCH_BOUNDARY":
            state_key = sample, event["chromosome"]
            state = orientation.get(state_key, 0)
            if before != str(state) or after != str(1 - state) or event["affected_field"] != "PHASE_ORIENTATION":
                raise MeasurementError("d17_boundary_event_invalid")
            orientation[state_key] = 1 - state
            boundaries.add((marker, sample))
            continue
        current = _gt(alleles[marker, sample])
        if event["affected_field"] != "GT" or current != before or before == after or after not in GT_VALUES:
            raise MeasurementError("d17_event_before_after_mismatch")
        if kind.startswith("PHASE_") and (before not in {"0|1", "1|0"} or after != before[::-1]):
            raise MeasurementError("d17_phase_changed_dosage")
        if kind == "PHASE_SWITCH_GT" and not orientation.get((sample, event["chromosome"]), 0):
            raise MeasurementError("d17_switch_outside_block")
        if kind == "GENOTYPE" and (after == ".|." or before.count("1") == after.count("1")):
            raise MeasurementError("d17_genotype_dosage_unchanged")
        if kind == "MISSING_GT" and after != ".|.":
            raise MeasurementError("d17_missing_event_invalid")
        alleles[marker, sample] = (-1, -1) if after == ".|." else (int(after[0]), int(after[2]))
        if kind == "GENOTYPE":
            genotype_alleles[marker, sample] = alleles[marker, sample]
        elif kind == "PHASE_SWITCH_GT":
            switched.add((marker, sample))
    expected: set[tuple[int, int]] = set()
    eligible_boundaries: set[tuple[int, int]] = set()
    for sample in range(len(panel.samples)):
        for chromosome in dict.fromkeys(row[0] for row in panel.variants):
            markers_in_track = [int(marker) for marker in np.flatnonzero(
                retained & (genotype_alleles[:, sample, 0] != genotype_alleles[:, sample, 1]))
                if panel.variants[marker][0] == chromosome]
            eligible_boundaries.update((marker, sample) for marker in markers_in_track[1:])
            state = 0
            for marker in markers_in_track:
                if (marker, sample) in boundaries:
                    state = 1 - state
                if state:
                    expected.add((marker, sample))
    if boundaries - eligible_boundaries or switched != expected:
        raise MeasurementError("d17_switch_block_incomplete")
    return alleles, retained


def write_panel(
    handle: TextIO, panel: MeasurementPanel, alleles: NDArray[np.int8],
    retained: NDArray[np.bool_], *, phased: bool,
) -> None:
    """Écrit un VCF minimal sans profil, graine ou étiquette biologique cachée."""
    handle.write("##fileformat=VCFv4.2\n")
    for chromosome in dict.fromkeys(row[0] for row in panel.variants):
        handle.write(f"##contig=<ID={chromosome}>\n")
    handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(panel.samples) + "\n")
    for index in np.flatnonzero(retained):
        row = panel.variants[index]
        genotypes = [_gt(pair) if phased else "./." if pair[0] < 0 else "/".join(map(str, sorted(pair)))
                     for pair in alleles[index]]
        # INFO/QUAL ne transportent aucune annotation de simulation.
        handle.write("\t".join((*row[:5], ".", "PASS", ".", "GT", *genotypes)) + "\n")
