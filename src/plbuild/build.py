"""The build and validate pipelines behind `plbuild build` and `plbuild validate`."""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from .config import Inputs
from .extract import ModelResult, extract_model
from .identity import assign_ids
from .process import Process
from .select import draw, find_gaps
from .templates import templates
from .validate import Check, check_antimony, check_bounds, check_vocabulary, read_library, reconstruct
from .writer import write_outputs


class BuildError(Exception):
    pass


@dataclass
class BuildResult:
    results: list[ModelResult]
    pool: list[Process]
    drawn: list[Process]
    paths: dict[str, Path]


def _hash_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sbml_path(entry) -> Path:
    if entry.sbml is None:
        raise BuildError(f"model {entry.id}: no local `sbml:` file given, and `plbuild fetch` is not implemented yet (L1)")
    if not entry.sbml.exists():
        raise BuildError(f"model {entry.id}: SBML file not found: {entry.sbml}")
    return entry.sbml


def build(inputs: Inputs) -> BuildResult:
    cfg = inputs.config
    results = [extract_model(e, _sbml_path(e), inputs) for e in cfg.models]
    curated: dict[str, Process] = {}
    for r in results:
        assign_ids(r.processes)
        for p in r.processes:
            if p.pid in curated:   # same process in two models: keep the first entry, list every source (L5 merges values)
                curated[p.pid].sources.extend(p.sources)
                p.pid = curated[p.pid].pid
            else:
                curated[p.pid] = p
    tmpl = templates(cfg.template_bounds)
    pool = list(curated.values()) + tmpl
    drawn = draw(pool, cfg.library_size, cfg.library_seed)
    gaps = find_gaps(list(curated.values()), drawn, tmpl)

    hashes = {"plbuild.yaml": _hash_file(cfg.path), "vocabulary.yaml": _hash_file(cfg.vocabulary_path)}
    if cfg.overrides_path:
        hashes["overrides.yaml"] = _hash_file(cfg.overrides_path)
    for r in results:
        hashes[f"sbml:{r.model_id}"] = _hash_file(r.sbml_path)
    paths = write_outputs(cfg.output_dir, results=results, pool=pool, drawn=drawn, gaps=gaps,
                          size=cfg.library_size, seed=cfg.library_seed, time_unit=cfg.time_unit, input_hashes=hashes)
    # keep what the reconstruction test needs next to the outputs
    recon = {r.model_id: {"sbml": str(r.sbml_path), "species_ids": r.species_ids, "initial": r.initial_conc,
                          "volumes": r.volumes, "complete": not r.excluded}
             for r in results}
    (cfg.output_dir / "reconstruction.json").write_text(json.dumps(recon, indent=2) + "\n", encoding="utf-8")
    return BuildResult(results, pool, drawn, paths)


def validate(inputs: Inputs) -> list[Check]:
    out = inputs.config.output_dir
    lib, pool_file, recon_file = out / "curated.ant", out / "pool.ant", out / "reconstruction.json"
    for f in (lib, pool_file, recon_file):
        if not f.exists():
            raise BuildError(f"{f} not found; run `plbuild build` first")
    library = read_library(lib)
    pool = read_library(pool_file)
    checks = [check_antimony(library), check_bounds(library), check_vocabulary(library, inputs)]
    recon = json.loads(recon_file.read_text(encoding="utf-8"))
    spans = {m.id: m.time_span for m in inputs.config.models}
    for model_id, info in recon.items():
        if not info["complete"]:
            checks.append(Check(f"reconstruction {model_id}", True,
                                ["skipped: some reactions were excluded (rate comparison at the initial state is not "
                                 "implemented yet)"]))
            continue
        checks.append(reconstruct(pool, model_id, Path(info["sbml"]), info["species_ids"], info["initial"],
                                  info["volumes"], spans.get(model_id, 100.0)))
    return checks
