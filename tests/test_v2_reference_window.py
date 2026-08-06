import hashlib
import json
import subprocess
from pathlib import Path
from typing import Sequence

import pytest

from effet_fondateur.references import (
    ReferenceAsset,
    ReferenceWindowError,
    ReferenceWindowIntegrityError,
    ResolvedReferencePanel,
    cache_reference_panel,
    extract_reference_window,
)


def _md5(content: bytes) -> str:
    return hashlib.md5(content, usedforsecurity=False).hexdigest()


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _cached_reference(tmp_path: Path):
    base_url = "https://ftp.1000genomes.ebi.ac.uk/synthetic"
    contents = {
        f"{base_url}/panel.vcf.gz": b"source vcf",
        f"{base_url}/panel.vcf.gz.tbi": b"source index",
        f"{base_url}/README": b"readme",
        f"{base_url}/manifest.txt": b"manifest",
    }
    resolved = ResolvedReferencePanel(
        catalog_id="synthetic_catalog_v1",
        catalog_sha256="a" * 64,
        panel_id="synthetic_panel_v1",
        provider="IGSR_1000_GENOMES",
        release_id="synthetic_release_v1",
        assembly="GRCh38",
        chromosome=19,
        sample_count=2,
        sample_scope="ALL_RELEASE_SAMPLES",
        population_scope="ALL_AVAILABLE_1KGP_POPULATIONS",
        readme_url=f"{base_url}/README",
        readme_sha256=_sha256(contents[f"{base_url}/README"]),
        manifest_url=f"{base_url}/manifest.txt",
        manifest_sha256=_sha256(contents[f"{base_url}/manifest.txt"]),
        vcf=ReferenceAsset(
            "panel.vcf.gz",
            f"{base_url}/panel.vcf.gz",
            _md5(contents[f"{base_url}/panel.vcf.gz"]),
        ),
        index=ReferenceAsset(
            "panel.vcf.gz.tbi",
            f"{base_url}/panel.vcf.gz.tbi",
            _md5(contents[f"{base_url}/panel.vcf.gz.tbi"]),
        ),
    )

    def downloader(url: str, destination: Path, timeout: int, chunk_size: int) -> None:
        destination.write_bytes(contents[url])

    cached = cache_reference_panel(
        resolved, tmp_path / "cache", downloader=downloader
    )
    return resolved, cached


class FakeBcftools:
    def __init__(
        self,
        *,
        source_samples: tuple[str, ...] = ("HG00001", "HG00002"),
        output_samples: tuple[str, ...] | None = None,
        contig_length: int = 58_617_616,
        unphased_positions: str = "",
        variant_count: int = 2,
    ) -> None:
        self.source_samples = source_samples
        self.output_samples = output_samples or source_samples
        self.contig_length = contig_length
        self.unphased_positions = unphased_positions
        self.variant_count = variant_count
        self.commands: list[list[str]] = []

    def __call__(
        self, command: Sequence[str], timeout_seconds: float
    ) -> subprocess.CompletedProcess[str]:
        arguments = list(command)
        self.commands.append(arguments)
        if arguments[1:] == ["--version"]:
            return subprocess.CompletedProcess(arguments, 0, "bcftools 1.21\n", "")
        if arguments[1:3] == ["view", "--header-only"]:
            header = (
                "##fileformat=VCFv4.2\n"
                f"##contig=<ID=chr19,length={self.contig_length}>\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">\n'
            )
            return subprocess.CompletedProcess(arguments, 0, header, "")
        if arguments[1:3] == ["query", "--list-samples"]:
            samples = (
                self.output_samples
                if Path(arguments[-1]).name.startswith("reference.chr")
                else self.source_samples
            )
            return subprocess.CompletedProcess(
                arguments, 0, "".join(f"{sample}\n" for sample in samples), ""
            )
        if arguments[1] == "view" and "--output" in arguments:
            output_path = Path(arguments[arguments.index("--output") + 1])
            output_path.write_bytes(b"bgzip reference window")
            return subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[1:3] == ["index", "--tbi"]:
            Path(f"{arguments[-1]}.tbi").write_bytes(b"tabix index")
            return subprocess.CompletedProcess(arguments, 0, "", "")
        if arguments[1:3] == ["index", "--nrecords"]:
            return subprocess.CompletedProcess(
                arguments, 0, f"{self.variant_count}\n", ""
            )
        if arguments[1] == "query" and "--include" in arguments:
            return subprocess.CompletedProcess(
                arguments, 0, self.unphased_positions, ""
            )
        raise AssertionError(f"unexpected command: {arguments}")


def test_extracts_all_samples_and_publishes_verified_manifest(tmp_path: Path) -> None:
    resolved, cached = _cached_reference(tmp_path)
    runner = FakeBcftools()

    extracted = extract_reference_window(
        resolved,
        cached,
        tmp_path / "window",
        start_bp=1_000,
        end_bp=2_000,
        bcftools_command="bcftools",
        threads=2,
        command_runner=runner,
    )

    assert extracted.sample_count == 2
    assert extracted.variant_count == 2
    assert extracted.vcf_path.is_file()
    assert extracted.index_path.is_file()
    manifest = json.loads(extracted.manifest_path.read_text(encoding="utf-8"))
    assert manifest["region"] == {"start_bp": 1_000, "end_bp": 2_000}
    assert manifest["source"]["cache_entry_key"] == cached.cache_entry_key
    extraction_command = next(
        command for command in runner.commands if "--regions" in command
    )
    assert extraction_command[extraction_command.index("--regions") + 1] == "chr19:1000-2000"
    assert "--samples" not in extraction_command
    assert "--samples-file" not in extraction_command


def test_blocks_sample_loss_or_reordering_without_publishing(tmp_path: Path) -> None:
    resolved, cached = _cached_reference(tmp_path)
    runner = FakeBcftools(output_samples=("HG00002", "HG00001"))
    output_dir = tmp_path / "window"

    with pytest.raises(
        ReferenceWindowIntegrityError, match="sample_set_or_order_mismatch"
    ):
        extract_reference_window(
            resolved,
            cached,
            output_dir,
            start_bp=1_000,
            end_bp=2_000,
            command_runner=runner,
        )

    assert not output_dir.exists()


def test_blocks_unphased_called_genotype_without_publishing(tmp_path: Path) -> None:
    resolved, cached = _cached_reference(tmp_path)
    runner = FakeBcftools(unphased_positions="chr19\t1500\n")
    output_dir = tmp_path / "window"

    with pytest.raises(ReferenceWindowIntegrityError, match="unphased_genotype"):
        extract_reference_window(
            resolved,
            cached,
            output_dir,
            start_bp=1_000,
            end_bp=2_000,
            command_runner=runner,
        )

    assert not output_dir.exists()


def test_blocks_non_grch38_contig_length(tmp_path: Path) -> None:
    resolved, cached = _cached_reference(tmp_path)

    with pytest.raises(ReferenceWindowIntegrityError, match="grch38_contig_mismatch"):
        extract_reference_window(
            resolved,
            cached,
            tmp_path / "window",
            start_bp=1_000,
            end_bp=2_000,
            command_runner=FakeBcftools(contig_length=59_128_983),
        )


def test_rejects_invalid_bounds_before_running_bcftools(tmp_path: Path) -> None:
    resolved, cached = _cached_reference(tmp_path)
    runner = FakeBcftools()

    with pytest.raises(ReferenceWindowError, match="invalid_reference_window_bounds"):
        extract_reference_window(
            resolved,
            cached,
            tmp_path / "window",
            start_bp=2_000,
            end_bp=1_000,
            command_runner=runner,
        )

    assert runner.commands == []


def test_detects_modified_cached_vcf_before_extraction(tmp_path: Path) -> None:
    resolved, cached = _cached_reference(tmp_path)
    cached.vcf_path.chmod(0o644)
    cached.vcf_path.write_bytes(b"modified")
    runner = FakeBcftools()

    with pytest.raises(ReferenceWindowIntegrityError, match="cache_size_mismatch"):
        extract_reference_window(
            resolved,
            cached,
            tmp_path / "window",
            start_bp=1_000,
            end_bp=2_000,
            command_runner=runner,
        )

    assert runner.commands == []
