"""D19 : comparaisons techniques pures, sans appel d'outil ni inférence clinique.

La conversion est comparée à l'observation D17 remise au pipeline, pas à la
vérité biologique antérieure aux erreurs. Les manques comptent comme états
à conserver. BED variant-major : spécification primaire PLINK 1.9,
https://www.cog-genomics.org/plink/1.9/formats#bed (consultée 30 octobre 2026).
A1/A2 ne sont jamais assimilés à REF/ALT ; seuls les allèles décodés comptent.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from effet_fondateur.calibration.materialize_ibd import _safe_file
from effet_fondateur.calibration.simulation.measurement_errors import MeasurementPanel


METHOD_ID = "independent_technical_concordance_v1"
BASES = "ACGT"
STAGES = ("00", "01", "02", "03", "04")


class TechnicalEvaluationError(ValueError):
    """Erreur technique expurgée de toute donnée ou clé individuelle."""


@dataclass(frozen=True)
class GenotypeMatrix:
    """Axes ordonnés et génotypes nucléotidiques triés ; -1/-1 = absence.

    Les clés d'individu sont FID/IID ; celles des variants sont chr/bp/ID.
    Les deux indices alléliques vont de 0 à 3 pour A,C,G,T, pas REF/ALT.
    """

    samples: tuple[tuple[str, str], ...]
    variants: tuple[tuple[int, int, str], ...]
    alleles: NDArray[np.int8]


def result(stage: str, metrics: dict[str, int | float | None], reasons: Sequence[str] = (), *, missing: bool = False) -> dict[str, Any]:
    """Un défaut observé prime sur l'absence de preuve ; aucun manque n'est imputé."""
    return dict(report_id=f"CAL-REP-{stage}", stage_id=stage,
                status="FAIL" if reasons else "NOT_EVALUABLE" if missing else "PASS",
                reasons=list(reasons) or (["mandatory_evidence_missing"] if missing else []), metrics=metrics)


def unavailable(stage: str, reason: str) -> dict[str, Any]:
    """Absence d'exécution ou de vérité : métriques sans mesure explicitement nulles."""
    return dict(report_id=f"CAL-REP-{stage}", stage_id=stage, status="NOT_EVALUABLE",
                reasons=[reason], metrics={"evaluated_count": None})


def compare_rows(
    expected: Sequence[Mapping[str, Any]], actual: Sequence[Mapping[str, Any]],
    *, keys: Sequence[str], fields: Sequence[str],
) -> dict[str, int]:
    """Compte pertes, ajouts, doublons, désordres et substitutions séparément."""
    expected_keys = [tuple(row[field] for field in keys) for row in expected]
    actual_keys = [tuple(row[field] for field in keys) for row in actual]
    if len(set(expected_keys)) != len(expected_keys):
        raise TechnicalEvaluationError("d19_expected_duplicate_key")
    expected_by_key = dict(zip(expected_keys, expected, strict=True))
    counts = Counter(actual_keys)
    mismatches = sum(any(row[field] != expected_by_key[key][field] for field in fields)
                     for key, row in zip(actual_keys, actual, strict=True) if key in expected_by_key)
    return dict(lost_count=len(set(expected_keys) - set(actual_keys)),
                extra_count=len(set(actual_keys) - set(expected_keys)),
                duplicate_count=sum(count - 1 for count in counts.values()),
                order_mismatch_count=sum(left != right for left, right in zip(expected_keys, actual_keys)) + abs(len(expected_keys) - len(actual_keys)),
                value_mismatch_count=mismatches, compared_count=len(actual))


def row_result(stage: str, expected: Sequence[Mapping[str, Any]], actual: Sequence[Mapping[str, Any]], *, keys: Sequence[str], fields: Sequence[str]) -> dict[str, Any]:
    """Compare une table publiée à une référence indépendante, sans reporter les clés."""
    metrics = compare_rows(expected, actual, keys=keys, fields=fields)
    return result(stage, metrics, [name for name, value in metrics.items() if name != "compared_count" and value])


def observation_matrix(panel: MeasurementPanel, samples: Sequence[Mapping[str, Any]]) -> GenotypeMatrix:
    """Convertit un VCF observé en allèles physiques sans changer dosage ni manque."""
    if tuple(row["SAMPLE_ID"] for row in samples) != panel.samples:
        raise TechnicalEvaluationError("d19_expected_sample_order_mismatch")
    alleles = np.full(panel.alleles.shape, -1, dtype=np.int8)
    for index, variant in enumerate(panel.variants):
        for allele_index, base in enumerate((variant[3], variant[4])):
            alleles[index][panel.alleles[index] == allele_index] = BASES.index(base)
    return GenotypeMatrix(tuple((str(row["FID"]), str(row["IID"])) for row in samples),
                          tuple((int(row[0].removeprefix("chr")), int(row[1]), row[2]) for row in panel.variants),
                          np.sort(alleles, axis=2))


def read_plink(prefix: Path) -> GenotypeMatrix:
    """Décode BED/BIM/FAM autosomiques, y compris padding et allèles monomorphes.

    Refuse format sample-major, troncature, octets supplémentaires, doublons,
    padding non nul ou appel d'un allèle BIM absent (0). Aucun complément de
    brin, dosage de référence ou imputation n'est supposé.
    """
    for suffix in ("bed", "bim", "fam"):
        _safe_file(prefix.with_suffix("." + suffix))
    fam = [line.split() for line in prefix.with_suffix(".fam").read_text().splitlines() if line.strip()]
    bim = [line.split() for line in prefix.with_suffix(".bim").read_text().splitlines() if line.strip()]
    if not fam or not bim or any(len(row) != 6 for row in (*fam, *bim)):
        raise TechnicalEvaluationError("d19_plink_columns_invalid")
    samples = tuple((row[0], row[1]) for row in fam)
    try:
        variants = tuple((int(row[0]), int(row[3]), row[1]) for row in bim)
    except ValueError as error:
        raise TechnicalEvaluationError("d19_bim_coordinates_invalid") from error
    if len(set(samples)) != len(samples) or len({row[2] for row in variants}) != len(variants) or len({row[:2] for row in variants}) != len(variants):
        raise TechnicalEvaluationError("d19_plink_duplicate")
    if any(not 1 <= row[0] <= 22 or row[1] < 1 for row in variants):
        raise TechnicalEvaluationError("d19_bim_domain_invalid")
    block_size = (len(fam) + 3) // 4
    payload = prefix.with_suffix(".bed").read_bytes()
    if payload[:3] != b"\x6c\x1b\x01" or len(payload) != 3 + len(bim) * block_size:
        raise TechnicalEvaluationError("d19_bed_header_or_size_invalid")
    alleles = np.full((len(bim), len(fam), 2), -1, dtype=np.int8)
    for marker, row in enumerate(bim):
        bases = row[4:6]
        if any(base not in (*BASES, "0") for base in bases):
            raise TechnicalEvaluationError("d19_bim_allele_invalid")
        block = payload[3 + marker * block_size:3 + (marker + 1) * block_size]
        if len(fam) % 4 and block[-1] >> (2 * (len(fam) % 4)):
            raise TechnicalEvaluationError("d19_bed_padding_invalid")
        for sample in range(len(fam)):
            code = (block[sample // 4] >> (2 * (sample % 4))) & 3
            if code == 1:
                continue
            pair = (bases[0], bases[0]) if code == 0 else tuple(bases) if code == 2 else (bases[1], bases[1])
            if "0" in pair:
                raise TechnicalEvaluationError("d19_bed_called_absent_allele")
            alleles[marker, sample] = sorted(BASES.index(base) for base in pair)
    return GenotypeMatrix(samples, variants, alleles)


def compare_genotypes(expected: GenotypeMatrix, actual: GenotypeMatrix, *, stage: str = "03") -> dict[str, Any]:
    """Compare axes stricts et tous les génotypes alignables, manques inclus.

    Taux = nombre de cellules discordantes / cellules comparées ; dénominateur
    nul ou alignement impossible → taux nul au sens JSON (null), jamais 0.0.
    Une permutation correcte reste signalée comme désordre contractuel.
    """
    metrics: dict[str, int | float | None] = {}
    reasons: list[str] = []
    if len(set(expected.samples)) != len(expected.samples) or len(set(expected.variants)) != len(expected.variants):
        raise TechnicalEvaluationError("d19_expected_duplicate_key")
    for label, wanted, observed in (("sample", expected.samples, actual.samples), ("variant", expected.variants, actual.variants)):
        counts = Counter(observed)
        metrics.update({f"{label}_lost_count": len(set(wanted) - set(observed)),
                        f"{label}_extra_count": len(set(observed) - set(wanted)),
                        f"{label}_duplicate_count": sum(count - 1 for count in counts.values()),
                        f"{label}_order_mismatch_count": sum(left != right for left, right in zip(wanted, observed)) + abs(len(wanted) - len(observed))})
    reasons.extend(name for name, value in metrics.items() if value)
    if not reasons or all("order" in name for name in reasons):
        if expected.alleles.shape != (len(expected.variants), len(expected.samples), 2) or actual.alleles.shape != (len(actual.variants), len(actual.samples), 2):
            raise TechnicalEvaluationError("d19_matrix_shape_invalid")
        sample_indices = {key: index for index, key in enumerate(actual.samples)}
        variant_indices = {key: index for index, key in enumerate(actual.variants)}
        aligned = actual.alleles[[variant_indices[key] for key in expected.variants]][:, [sample_indices[key] for key in expected.samples]]
        mismatches = int(np.count_nonzero(np.any(expected.alleles != aligned, axis=2)))
        compared = len(expected.samples) * len(expected.variants)
        metrics.update(genotype_mismatch_count=mismatches, genotype_compared_count=compared,
                       genotype_discordance_rate=mismatches / compared if compared else None,
                       missing_call_mismatch_count=int(np.count_nonzero((expected.alleles[:, :, 0] < 0) != (aligned[:, :, 0] < 0))))
        if mismatches:
            reasons.append("genotype_mismatch")
        if compared == 0:
            return result(stage, metrics, reasons, missing=True)
    else:
        metrics.update(genotype_mismatch_count=None, genotype_compared_count=0,
                       genotype_discordance_rate=None, missing_call_mismatch_count=None)
    return result(stage, metrics, reasons)


def evaluate_target(expected_base: GenotypeMatrix, expected_target: GenotypeMatrix | None, actual: GenotypeMatrix | None) -> dict[str, Any]:
    """Évalue le jeu enrichi 04 : cible indépendante et fond préservé.

    La cible doit décrire un variant absent du jeu de base, comme l'exige 04.
    Elle provient d'une attente explicite fournie à l'évaluateur, jamais du
    groupe, du diagnostic ou de l'export observé employé comme sa propre vérité.
    """
    if actual is None or expected_target is None:
        return unavailable("04", "target_truth_or_output_missing")
    if expected_target.samples != expected_base.samples or len(expected_target.variants) != 1 or any(
        expected_target.variants[0][:2] == variant[:2] or expected_target.variants[0][2] == variant[2] for variant in expected_base.variants
    ):
        raise TechnicalEvaluationError("d19_target_expectation_invalid")
    wanted = GenotypeMatrix(expected_base.samples, expected_base.variants + expected_target.variants,
                            np.concatenate((expected_base.alleles, expected_target.alleles), axis=0))
    # PLINK peut placer la cible à sa coordonnée ; seul l'ordre canonique
    # chr/bp/ID est attendu dans le jeu enrichi, pas l'ordre de concaténation.
    order = sorted(range(len(wanted.variants)), key=lambda index: wanted.variants[index])
    wanted = GenotypeMatrix(wanted.samples, tuple(wanted.variants[index] for index in order), wanted.alleles[order])
    comparison = compare_genotypes(wanted, actual, stage="04")
    indices = [index for index, variant in enumerate(actual.variants) if variant == expected_target.variants[0]]
    observed_target = GenotypeMatrix(actual.samples, tuple(actual.variants[index] for index in indices), actual.alleles[indices])
    target_comparison = compare_genotypes(expected_target, observed_target, stage="04")
    comparison["metrics"].update({f"target_{name}": value for name, value in target_comparison["metrics"].items()})
    return comparison
