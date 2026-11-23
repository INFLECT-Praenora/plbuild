"""Command-line entry point (spec section 13). Commands other than config checking are stubs until their milestone."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from . import build as build_pipeline
from .config import load_inputs
from .yamlload import ConfigError

COMMANDS = {
    "search": ("List candidate BioModels models; changes nothing", "L1"),
    "fetch": ("Download listed models into the cache", "L1"),
    "build": ("Extract, map, draw the library; write curated.ant, pool.ant, report, manifest", "L5"),
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


def _print_checks(checks) -> bool:
    ok = True
    for c in checks:
        print(f"  [{'pass' if c.ok else 'FAIL'}] {c.name}")
        for m in c.messages[:20]:
            print(f"         {m}")
        ok &= c.ok
    return ok


def _build(inputs, args) -> int:
    res = build_pipeline.build(inputs)
    for r in res.results:
        print(f"{r.model_id}: {r.n_reactions} reactions, {len(r.processes)} kept, {len(r.excluded)} excluded, "
              f"{len(r.unresolved)} unresolved names")
    gaps = json.loads(res.paths["report_json"].read_text(encoding="utf-8"))["totals"]
    print(f"library: {len(res.drawn)} of {len(res.pool)} candidates; {gaps['dropped_curated']} curated left out, "
          f"{gaps['gaps']} of them with no drawn replacement")
    print(f"wrote {res.paths['library']}")
    if args.strict and any(r.unresolved for r in res.results):
        print("plbuild: error: unresolved species (--strict)", file=sys.stderr)
        return 1
    if not args.no_validate:
        print("validate:")
        if not _print_checks(build_pipeline.validate(inputs)):
            return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    inputs = None
    if args.command != "search":
        try:
            inputs = load_inputs(args.config)
        except ConfigError as e:
            print(f"plbuild: error: {e}", file=sys.stderr)
            return 1
    try:
        if args.command == "build":
            return _build(inputs, args)
        if args.command == "validate":
            print("validate:")
            return 0 if _print_checks(build_pipeline.validate(inputs)) else 1
    except (build_pipeline.BuildError, ValueError, RuntimeError) as e:
        print(f"plbuild: error: {e}", file=sys.stderr)
        return 1
    milestone = COMMANDS[args.command][1]
    print(f"plbuild {args.command}: not implemented yet (milestone {milestone})", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
