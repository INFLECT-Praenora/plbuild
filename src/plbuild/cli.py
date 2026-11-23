"""Command-line entry point (spec section 13). Commands other than config checking are stubs until their milestone."""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .config import load_inputs
from .yamlload import ConfigError

COMMANDS = {
    "search": ("List candidate BioModels models; changes nothing", "L1"),
    "fetch": ("Download listed models into the cache", "L1"),
    "build": ("Extract, map, convert, merge; write curated.ant, vocabulary, report, manifest", "L5"),
    "validate": ("Run the validation checks on the current output", "L6"),
    "explain": ("Show how one reaction was handled, step by step", "L5"),
}


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="plbuild", description="Build the EnsGen process library from BioModels models.")
    p.add_argument("--version", action="version", version=f"plbuild {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>", required=True)
    for name, (help_text, _) in COMMANDS.items():
        sp = sub.add_parser(name, help=help_text, description=help_text)
        sp.add_argument("--config", default="plbuild.yaml", help="path to plbuild.yaml (default: %(default)s)")
        if name == "search":
            sp.add_argument("terms", nargs="+", help="search terms")
        elif name == "build":
            sp.add_argument("--no-validate", action="store_true", help="skip validation after the build")
            sp.add_argument("--strict", action="store_true", help="fail if any species is unresolved")
        elif name == "explain":
            sp.add_argument("model", help="BioModels id")
            sp.add_argument("reaction", help="reaction id in that model")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command != "search":
        try:
            load_inputs(args.config)
        except ConfigError as e:
            print(f"plbuild: error: {e}", file=sys.stderr)
            return 1
    milestone = COMMANDS[args.command][1]
    print(f"plbuild {args.command}: not implemented yet (milestone {milestone})", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
