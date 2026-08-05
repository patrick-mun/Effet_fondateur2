"""Interface en ligne de commande du pipeline V2."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from effet_fondateur.contracts import ConfigurationError, load_pipeline_config
from effet_fondateur.orchestrator import PipelineError, resume_pipeline, run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="effet-fondateur")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser(
        "validate-config",
        help="Valider une configuration YAML sans lancer le pipeline.",
    )
    validate_parser.add_argument("config", type=Path)

    run_parser = subparsers.add_parser(
        "run",
        help="Créer et exécuter un run V2 dans un dossier immuable.",
    )
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--runs-dir", type=Path, default=Path("data/runs"))

    resume_parser = subparsers.add_parser(
        "resume",
        help="Reprendre un run V2 après validation de son état.",
    )
    resume_parser.add_argument("--run-dir", type=Path, required=True)
    return parser


def main(arguments: Sequence[str] | None = None) -> int:
    """Valide les arguments et exécute la commande V2 demandée."""
    parser = _build_parser()
    parsed_arguments = parser.parse_args(arguments)

    if parsed_arguments.command == "validate-config":
        try:
            load_pipeline_config(parsed_arguments.config)
        except ConfigurationError as error:
            parser.error(str(error))
        print(f"Configuration valide : {parsed_arguments.config}")
        return 0

    if parsed_arguments.command == "run":
        try:
            run_dir = run_pipeline(parsed_arguments.config, parsed_arguments.runs_dir)
        except (ConfigurationError, PipelineError, OSError, ValueError) as error:
            parser.error(str(error))
        print(f"Run V2 créé : {run_dir}")
        return 0

    if parsed_arguments.command == "resume":
        try:
            run_dir = resume_pipeline(parsed_arguments.run_dir)
        except (ConfigurationError, PipelineError, OSError, ValueError) as error:
            parser.error(str(error))
        print(f"Run V2 repris : {run_dir}")
        return 0

    parser.error(f"Commande inconnue : {parsed_arguments.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
