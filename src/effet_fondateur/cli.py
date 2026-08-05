"""Interface en ligne de commande du pipeline V2."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from effet_fondateur.contracts import ConfigurationError, load_pipeline_config


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="effet-fondateur")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser(
        "validate-config",
        help="Valider une configuration YAML sans lancer le pipeline.",
    )
    validate_parser.add_argument("config", type=Path)
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

    parser.error(f"Commande inconnue : {parsed_arguments.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
