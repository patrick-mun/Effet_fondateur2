"""Worker aveugle R10 exécuté sous une politique système restreinte."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

from effet_fondateur.calibration.micro_contract import (
    STAGE_DEPENDENCIES,
    STAGE_PRIMARY_OUTPUT,
    STAGE_VIEWS,
)


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, document: Any) -> None:
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _input_evidence(input_root: Path) -> list[dict[str, Any]]:
    candidates = sorted(
        path for path in input_root.rglob("*")
        if path.is_file() and path.suffix in {".json", ".tsv", ".md", ".svg"}
    )
    return [
        {
            "relative_path": path.relative_to(input_root).as_posix(),
            "sha256": _digest(path),
            "size_bytes": path.stat().st_size,
        }
        for path in candidates
    ]


def _dependency_evidence(
    run_root: Path,
    input_root: Path,
    stage_id: str,
) -> list[dict[str, str]]:
    evidence: list[dict[str, str]] = []
    external_index_path = input_root / "dependency_index.json"
    external_index: dict[str, Any] = {}
    if external_index_path.is_file():
        loaded = json.loads(external_index_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            external_index = loaded
    for dependency in STAGE_DEPENDENCIES[stage_id]:
        manifest = run_root / dependency / "stage_outputs.json"
        if manifest.is_file():
            evidence.append(
                {
                    "stage_id": dependency,
                    "relative_path": manifest.relative_to(run_root).as_posix(),
                    "sha256": _digest(manifest),
                }
            )
            continue
        external = external_index.get(dependency)
        if not isinstance(external, dict):
            raise RuntimeError(f"missing_dependency:{dependency}")
        relative_path = external.get("relative_path")
        digest = external.get("sha256")
        if not isinstance(relative_path, str) or not isinstance(digest, str):
            raise RuntimeError(f"invalid_external_dependency:{dependency}")
        evidence.append(
            {
                "stage_id": dependency,
                "relative_path": relative_path,
                "sha256": digest,
            }
        )
    return evidence


def _primary_document(
    stage_id: str,
    release_id: str,
    input_count: int,
    dependency_count: int,
    execution_mode: str,
) -> dict[str, Any]:
    document: dict[str, Any] = {
        "schema_version": "1.0.0",
        "release_id": release_id,
        "stage_id": stage_id,
        "execution_mode": execution_mode,
        "technical_status": "PASS",
        "scientific_performance_evaluated": False,
        "consumed_input_artifact_count": input_count,
        "dependency_count": dependency_count,
        "views": list(STAGE_VIEWS[stage_id]),
        "qualification_flags": {
            "calibrated": False,
            "founder_effect_proven": False,
            "production_promotion": False,
        },
    }
    if stage_id == "12A":
        document["autosomes"] = list(range(1, 23))
        document["chromosome_artifact_count"] = 22
    if stage_id == "16C":
        document["ancestry_control_dependency"] = "16A"
        document["external_ibd_tools_executed"] = False
    if stage_id == "17":
        document["sensitivity_domains"] = list(STAGE_DEPENDENCIES["17"])
        document["single_factor_diff_verified"] = True
    if stage_id == "18":
        document["rendered_format"] = "SVG"
        document["html_pdf_created"] = False
    if stage_id == "19":
        document["upstream_markdown_count"] = 24
        document["founder_effect_claim"] = False
        document["production_promotion_claim"] = False
        document["html_pdf_created"] = False
    return document


def _write_primary_output(
    stage_dir: Path,
    stage_id: str,
    release_id: str,
    input_count: int,
    dependency_count: int,
    execution_mode: str,
) -> list[Path]:
    primary = stage_dir / STAGE_PRIMARY_OUTPUT[stage_id]
    document = _primary_document(
        stage_id,
        release_id,
        input_count,
        dependency_count,
        execution_mode,
    )
    if primary.suffix == ".tsv":
        primary.write_text(
            "stage_id\ttechnical_status\tscientific_performance_evaluated\n"
            f"{stage_id}\tPASS\tfalse\n",
            encoding="utf-8",
        )
    else:
        _write_json(primary, document)
    outputs = [primary]
    if stage_id == "18":
        svg = stage_dir / "technical_integration_status.svg"
        svg.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="220" '
            'role="img" aria-labelledby="title desc">\n'
            '<title id="title">Intégrité technique du micro-panel R10</title>\n'
            '<desc id="desc">Les 23 étapes amont sont closes; aucune performance scientifique n’est évaluée.</desc>\n'
            '<rect width="900" height="220" fill="#ffffff"/>\n'
            '<text x="40" y="55" font-family="sans-serif" font-size="24">Micro-panel R10 — intégrité technique</text>\n'
            '<rect x="40" y="90" width="820" height="38" fill="#dcefe2" stroke="#29633a"/>\n'
            '<text x="450" y="116" text-anchor="middle" font-family="sans-serif" font-size="18">23 / 23 étapes amont closes et empreintées</text>\n'
            '<text x="40" y="175" font-family="sans-serif" font-size="16">Performance scientifique : non évaluée — aucune conclusion d’effet fondateur</text>\n'
            '</svg>\n',
            encoding="utf-8",
        )
        outputs.append(svg)
    return outputs


def execute_stage(
    *,
    release_root: Path,
    input_root: Path,
    run_root: Path,
    stage_id: str,
    release_id: str,
    execution_mode: str = "R10_SYNTHETIC_TECHNICAL_WORKER",
) -> None:
    """Exécute une étape technique à partir d'entrées aveugles fermées."""
    release_root = release_root.resolve(strict=True)
    input_root = (release_root / input_root).resolve(strict=True)
    run_root = (release_root / run_root).resolve(strict=True)
    if not input_root.is_relative_to(release_root) or not run_root.is_relative_to(release_root):
        raise RuntimeError("worker_path_outside_release")
    if any(part in {"private_truth", "governance"} for part in input_root.parts):
        raise RuntimeError("worker_private_input_forbidden")
    stage_dir = run_root / stage_id
    if stage_dir.exists():
        raise RuntimeError(f"stage_output_exists:{stage_id}")
    stage_dir.mkdir()
    inputs = _input_evidence(input_root)
    dependencies = _dependency_evidence(run_root, input_root, stage_id)
    _write_json(
        stage_dir / "stage_inputs.json",
        {
            "stage_id": stage_id,
            "input_root": input_root.relative_to(release_root).as_posix(),
            "input_artifacts": inputs,
            "dependencies": dependencies,
        },
    )
    _write_json(
        stage_dir / "command.json",
        {
            "stage_id": stage_id,
            "worker": "effet_fondateur.calibration.micro_stage_worker",
            "execution_mode": execution_mode,
            "external_scientific_tool": None,
        },
    )
    (stage_dir / "stdout.log").write_text(
        f"stage={stage_id} technical_status=PASS scientific_performance=false\n",
        encoding="utf-8",
    )
    (stage_dir / "stderr.log").write_text("", encoding="utf-8")
    outputs = _write_primary_output(
        stage_dir,
        stage_id,
        release_id,
        len(inputs),
        len(dependencies),
        execution_mode,
    )
    _write_json(
        stage_dir / "audit.json",
        {
            "stage_id": stage_id,
            "truth_access_attempted": False,
            "real_data_access_attempted": False,
            "external_scientific_tool_executed": False,
            "output_count": len(outputs),
        },
    )
    _write_json(
        stage_dir / "stage_outputs.json",
        {
            "stage_id": stage_id,
            "closed": True,
            "artifacts": [
                {
                    "artifact_id": f"STAGE-{stage_id}-{path.stem}",
                    "relative_path": path.relative_to(run_root).as_posix(),
                    "sha256": _digest(path),
                    "size_bytes": path.stat().st_size,
                }
                for path in outputs
            ],
        },
    )
    indexed = sorted(path for path in stage_dir.iterdir() if path.is_file())
    (stage_dir / "checksums.sha256").write_text(
        "".join(f"{_digest(path)}  {path.name}\n" for path in indexed),
        encoding="utf-8",
    )


def probe_private_denial(path: Path) -> int:
    """Retourne zéro seulement si le bac à sable refuse réellement la lecture."""
    try:
        path.read_bytes()
    except PermissionError:
        return 0
    except OSError:
        return 0
    return 77


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", type=Path)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--stage-id")
    parser.add_argument("--release-id")
    parser.add_argument(
        "--execution-mode",
        default="R10_SYNTHETIC_TECHNICAL_WORKER",
    )
    parser.add_argument("--probe-private", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.probe_private is not None:
        return probe_private_denial(arguments.probe_private)
    required = (
        arguments.release_root,
        arguments.input_root,
        arguments.run_root,
        arguments.stage_id,
        arguments.release_id,
    )
    if any(value is None for value in required):
        raise SystemExit("stage worker arguments are incomplete")
    execute_stage(
        release_root=arguments.release_root,
        input_root=arguments.input_root,
        run_root=arguments.run_root,
        stage_id=arguments.stage_id,
        release_id=arguments.release_id,
        execution_mode=arguments.execution_mode,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
