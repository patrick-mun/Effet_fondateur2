import hashlib
import json
from pathlib import Path

import pytest

from effet_fondateur.ancestry import AncestryExtractCacheError, cache_reference_extract


def _selection(tmp_path: Path) -> tuple[Path, Path]:
    positions = tmp_path / "positions.tsv"
    samples = tmp_path / "samples.txt"
    positions.write_text("chr2\t100\t100\n", encoding="utf-8")
    samples.write_text("R1\nR2\n", encoding="utf-8")
    return positions, samples


def _call(tmp_path: Path, *, offline: bool, extractor):
    positions, samples = _selection(tmp_path)
    return cache_reference_extract(
        cache_root=tmp_path / "cache",
        panel_id="panel",
        assembly="GRCh38",
        chromosome=2,
        source_url="https://ftp.1000genomes.ebi.ac.uk/reference.chr2.vcf.gz",
        source_vcf_md5=hashlib.md5(b"source", usedforsecurity=False).hexdigest(),
        source_index_md5=hashlib.md5(b"index", usedforsecurity=False).hexdigest(),
        positions_path=positions,
        samples_path=samples,
        offline=offline,
        timeout_seconds=30,
        extractor=extractor,
    )


def test_reference_extract_is_populated_then_reused_offline(tmp_path: Path) -> None:
    calls = 0

    def extractor(source, positions, samples, vcf, index, timeout):
        nonlocal calls
        calls += 1
        vcf.write_bytes(b"synthetic bgzip extract")
        index.write_bytes(b"synthetic tabix index")

    first = _call(tmp_path, offline=False, extractor=extractor)
    second = _call(
        tmp_path,
        offline=True,
        extractor=lambda *_: pytest.fail("network extraction must not run on hit"),
    )

    assert first.status == "POPULATED"
    assert second.status == "HIT"
    assert calls == 1


def test_reference_extract_offline_miss_and_corruption_block(tmp_path: Path) -> None:
    with pytest.raises(AncestryExtractCacheError, match="offline_miss"):
        _call(tmp_path, offline=True, extractor=lambda *_: None)

    cached = _call(
        tmp_path,
        offline=False,
        extractor=lambda source, positions, samples, vcf, index, timeout: (
            vcf.write_bytes(b"vcf"), index.write_bytes(b"index")
        ),
    )
    cached.entry_dir.chmod(0o755)
    cached.vcf_path.chmod(0o644)
    cached.vcf_path.write_bytes(b"corrupt")
    cached.vcf_path.chmod(0o444)
    cached.entry_dir.chmod(0o555)
    with pytest.raises(AncestryExtractCacheError, match="corrupt"):
        _call(tmp_path, offline=True, extractor=lambda *_: None)


def test_reference_extract_accepts_only_checksum_verified_absolute_local_source(tmp_path: Path) -> None:
    positions, samples = _selection(tmp_path)
    source = (tmp_path / "reference.chr2.vcf.gz").resolve()
    index = Path(f"{source}.tbi")
    source.write_bytes(b"official local source")
    index.write_bytes(b"official local index")

    def call(vcf_md5: str):
        return cache_reference_extract(
            cache_root=tmp_path / "cache-local",
            panel_id="panel",
            assembly="GRCh38",
            chromosome=2,
            source_url=str(source),
            source_vcf_md5=vcf_md5,
            source_index_md5=hashlib.md5(index.read_bytes(), usedforsecurity=False).hexdigest(),
            positions_path=positions,
            samples_path=samples,
            offline=False,
            timeout_seconds=30,
            extractor=lambda source_url, positions_path, samples_path, vcf, tbi, timeout: (
                vcf.write_bytes(b"vcf"), tbi.write_bytes(b"index")
            ),
        )

    cached = call(hashlib.md5(source.read_bytes(), usedforsecurity=False).hexdigest())
    assert cached.status == "POPULATED"
    manifest = json.loads(cached.manifest_path.read_text(encoding="utf-8"))
    assert manifest["method_id"] == "bcftools_local_variant_extract_v1"
    with pytest.raises(AncestryExtractCacheError, match="checksum_mismatch"):
        call("0" * 32)
