#!/usr/bin/env python3
"""Simulateur de commandes PLINK pour les seuls tests unitaires D19.

Ne calcule ni QC ni Mendel réel. Encode des BED complets pour que l'évaluateur
de concordance puisse contrôler les vrais adaptateurs 03 et 04 du pipeline.
"""

from pathlib import Path
import sys


def write_dataset(prefix: Path, families: list[list[str]], markers: list[list[str]], calls: list[list[tuple[str, str]]]) -> None:
    prefix.with_suffix(".fam").write_text("".join(" ".join(row) + "\n" for row in families))
    rows = []
    payload = bytearray.fromhex("6c1b01")
    for marker, genotypes in zip(markers, calls):
        bases = sorted({base for pair in genotypes for base in pair if base != "0"})
        bases = (bases + ["0", "0"])[:2]
        rows.append(" ".join([*marker[:4], *bases]) + "\n")
        block = bytearray((len(families) + 3) // 4)
        for index, pair in enumerate(genotypes):
            code = 1 if "0" in pair else 0 if pair == (bases[0], bases[0]) else 3 if pair == (bases[1], bases[1]) else 2
            block[index // 4] |= code << (2 * (index % 4))
        payload.extend(block)
    prefix.with_suffix(".bim").write_text("".join(rows))
    prefix.with_suffix(".bed").write_bytes(payload)


def read_dataset(prefix: Path) -> tuple[list[list[str]], list[list[str]], list[list[tuple[str, str]]]]:
    families = [row.split() for row in prefix.with_suffix(".fam").read_text().splitlines()]
    markers = [row.split() for row in prefix.with_suffix(".bim").read_text().splitlines()]
    payload = prefix.with_suffix(".bed").read_bytes()[3:]
    width = (len(families) + 3) // 4
    calls = []
    for marker_index, marker in enumerate(markers):
        genotypes = []
        for sample_index in range(len(families)):
            code = (payload[marker_index * width + sample_index // 4] >> (2 * (sample_index % 4))) & 3
            genotypes.append(("0", "0") if code == 1 else (marker[4], marker[4]) if code == 0 else (marker[5], marker[5]) if code == 3 else tuple(marker[4:6]))
        calls.append(genotypes)
    return families, markers, calls


def main() -> None:
    arguments = sys.argv
    if "--version" in arguments:
        print("PLINK fake D19 unit tests only")
        return
    output = Path(arguments[arguments.index("--out") + 1])
    if "--file" in arguments:
        source = Path(arguments[arguments.index("--file") + 1])
        pedigree = [row.split() for row in source.with_suffix(".ped").read_text().splitlines()]
        markers = [row.split() for row in source.with_suffix(".map").read_text().splitlines()]
        calls = [[tuple(row[6 + index * 2:8 + index * 2]) for row in pedigree] for index in range(len(markers))]
        write_dataset(output, [row[:6] for row in pedigree], markers, calls)
    elif "--bmerge" in arguments:
        base = Path(arguments[arguments.index("--bfile") + 1])
        additional = Path(arguments[arguments.index("--bmerge") + 1]).with_suffix("")
        families, markers, calls = read_dataset(base)
        other_families, other_markers, other_calls = read_dataset(additional)
        assert families == other_families
        markers += other_markers
        calls += other_calls
        order = sorted(range(len(markers)), key=lambda index: (int(markers[index][0]), int(markers[index][3]), markers[index][1]))
        write_dataset(output, families, [markers[index] for index in order], [calls[index] for index in order])
    elif "--mendel" in arguments:
        output.with_suffix(".mendel").write_text("FID KID CHR SNP CODE ERROR\n")
    else:
        raise SystemExit("unsupported fake PLINK command")


if __name__ == "__main__":
    main()
