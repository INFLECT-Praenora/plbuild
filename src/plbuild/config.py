"""Loading and checking the three curator files: plbuild.yaml, vocabulary.yaml and overrides.yaml (spec section 2).

Every problem is raised as ConfigError naming the file and line.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .templates import DEFAULT_BOUNDS as DEFAULT_TEMPLATE_BOUNDS
from .yamlload import ConfigError, LineDict, LineList, load_yaml

IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
BIOMODELS_ID = re.compile(r"^(BIOMD|MODEL)\d{10}$")
ENTITY_CLASSES = ("molecule", "complex", "cell", "other")


# --------------------------------------------------------------------------- checking helpers


class _Checker:
    def __init__(self, path: Path):
        self.path = path

    def fail(self, line: int | None, message: str):
        raise ConfigError(self.path, line, message)

    def mapping(self, value: Any, line: int | None, what: str) -> LineDict:
        if not isinstance(value, dict):
            self.fail(line, f"{what} must be a mapping")
        return value

    def seq(self, value: Any, line: int | None, what: str) -> LineList:
        if not isinstance(value, list):
            self.fail(line, f"{what} must be a list")
        return value

    def string(self, value: Any, line: int | None, what: str) -> str:
        if not isinstance(value, str) or not value:
            self.fail(line, f"{what} must be a non-empty string")
        return value

    def ident(self, value: Any, line: int | None, what: str) -> str:
        self.string(value, line, what)
        if not IDENT.match(value):
            self.fail(line, f"{what} '{value}' is not a valid Antimony identifier")
        return value

    def number(self, value: Any, line: int | None, what: str, positive: bool = False) -> float:
        if isinstance(value, str):   # PyYAML (YAML 1.1) reads 1e7 and 1.0e7 as strings; only 1.0e+7 is a float
            try:
                value = float(value)
            except ValueError:
                pass
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            self.fail(line, f"{what} must be a number")
        if positive and value <= 0:
            self.fail(line, f"{what} must be positive")
        return float(value)

    def boolean(self, value: Any, line: int | None, what: str) -> bool:
        if not isinstance(value, bool):
            self.fail(line, f"{what} must be true or false")
        return value

    def strings(self, value: Any, line: int | None, what: str) -> list[str]:
        items = self.seq(value, line, what)
        return [self.string(v, items.line_of(i), f"{what}[{i}]") for i, v in enumerate(items)]

    def keys(self, d: LineDict, allowed: set[str], what: str, required: set[str] = frozenset(),
             line: int | None = None):
        """Reject unknown keys; report a missing required key at `line` (the owning key's line) if given."""
        for k in d:
            if k not in allowed:
                self.fail(d.line_of(k), f"unknown key '{k}' in {what} (allowed: {', '.join(sorted(allowed))})")
        for k in sorted(required):
            if k not in d:
                self.fail(line or d.line, f"{what} is missing required key '{k}'")


# --------------------------------------------------------------------------- plbuild.yaml


@dataclass
class CompartmentSpec:
    name: str
    identifiers: list[str] = field(default_factory=list)
    synonyms: list[str] = field(default_factory=list)


@dataclass
class ModelEntry:
    id: str
    version: str | None = None
    allow_uncurated: bool = False
    default_compartment: str | None = None
    sbml: Path | None = None          # local SBML file; used instead of the cache until fetch (L1) exists
    time_span: float = 100.0          # simulation end time for the reconstruction test, in source time units


@dataclass
class BuildConfig:
    path: Path
    models: list[ModelEntry]
    compartments: dict[str, CompartmentSpec]
    time_unit: str
    bounds_factor: float
    template_bounds: tuple[float, float]  # (lower, upper) for template parameters other than Hill coefficients
    library_size: int | None            # number of entries drawn into the library; None = all
    library_seed: int                   # seed of that draw
    vocabulary_path: Path
    overrides_path: Path | None
    cache_dir: Path
    output_dir: Path


def _load_build_config(path: Path) -> BuildConfig:
    c = _Checker(path)
    top = load_yaml(path)
    top = c.mapping(top, 1, "plbuild.yaml")
    c.keys(top, {"models", "compartments", "library", "bounds", "vocabulary", "overrides", "cache", "output"},
           "plbuild.yaml", required={"models", "compartments", "vocabulary"})
    base = path.parent

    comps: dict[str, CompartmentSpec] = {}
    cd = c.mapping(top["compartments"], top.line_of("compartments"), "compartments")
    if not cd:
        c.fail(cd.line, "compartments must declare at least one compartment")
    for name, spec in cd.items():
        line = cd.line_of(name)
        c.ident(name, line, "compartment name")
        spec = LineDict() if spec is None else c.mapping(spec, line, f"compartment '{name}'")
        if spec:
            c.keys(spec, {"identifiers", "synonyms"}, f"compartment '{name}'")
        comps[name] = CompartmentSpec(
            name,
            c.strings(spec["identifiers"], spec.line_of("identifiers"), f"{name}.identifiers") if "identifiers" in spec else [],
            c.strings(spec["synonyms"], spec.line_of("synonyms"), f"{name}.synonyms") if "synonyms" in spec else [],
        )

    models: list[ModelEntry] = []
    seen: set[str] = set()
    ml = c.seq(top["models"], top.line_of("models"), "models")
    for i, entry in enumerate(ml):
        line = ml.line_of(i)
        if isinstance(entry, str):
            entry_d = LineDict({"id": entry})
            entry_d.line, entry_d.key_lines = line, {"id": line}
            entry = entry_d
        entry = c.mapping(entry, line, f"models[{i}]")
        c.keys(entry, {"id", "version", "allow_uncurated", "default_compartment", "sbml", "time_span"}, f"models[{i}]", required={"id"}, line=line)
        mid = c.string(entry["id"], entry.line_of("id"), "model id")
        if not BIOMODELS_ID.match(mid):
            c.fail(entry.line_of("id"), f"'{mid}' is not a BioModels id (BIOMD or MODEL followed by 10 digits)")
        if mid in seen:
            c.fail(entry.line_of("id"), f"model '{mid}' is listed twice")
        seen.add(mid)
        version = entry.get("version")
        if version is not None:
            if isinstance(version, bool) or not isinstance(version, (int, str)):
                c.fail(entry.line_of("version"), "version must be an integer or a string")
            version = str(version)
        allow = c.boolean(entry["allow_uncurated"], entry.line_of("allow_uncurated"), "allow_uncurated") \
            if "allow_uncurated" in entry else False
        dc = entry.get("default_compartment")
        if dc is not None and dc not in comps:
            c.fail(entry.line_of("default_compartment"),
                   f"default_compartment '{dc}' is not one of the declared compartments ({', '.join(comps)})")
        sbml = base / c.string(entry["sbml"], entry.line_of("sbml"), "sbml") if "sbml" in entry else None
        span = c.number(entry["time_span"], entry.line_of("time_span"), "time_span", positive=True)             if "time_span" in entry else 100.0
        models.append(ModelEntry(mid, version, allow, dc, sbml, span))

    time_unit, size, seed = "hour", None, 0
    if "library" in top:
        lib = c.mapping(top["library"], top.line_of("library"), "library")
        c.keys(lib, {"time_unit", "size", "seed"}, "library")
        if "time_unit" in lib:
            time_unit = c.string(lib["time_unit"], lib.line_of("time_unit"), "library.time_unit")
        if "size" in lib:
            if isinstance(lib["size"], bool) or not isinstance(lib["size"], int) or lib["size"] < 1:
                c.fail(lib.line_of("size"), "library.size must be a positive integer")
            size = lib["size"]
        if "seed" in lib:
            if isinstance(lib["seed"], bool) or not isinstance(lib["seed"], int) or lib["seed"] < 0:
                c.fail(lib.line_of("seed"), "library.seed must be a non-negative integer")
            seed = lib["seed"]

    factor, tbounds = 100.0, DEFAULT_TEMPLATE_BOUNDS
    if "bounds" in top:
        b = c.mapping(top["bounds"], top.line_of("bounds"), "bounds")
        c.keys(b, {"factor", "templates"}, "bounds")
        if "factor" in b:
            factor = c.number(b["factor"], b.line_of("factor"), "bounds.factor", positive=True)
            if factor <= 1:
                c.fail(b.line_of("factor"), "bounds.factor must be greater than 1")
        if "templates" in b:
            tl = b.line_of("templates")
            pair = c.seq(b["templates"], tl, "bounds.templates")
            if len(pair) != 2:
                c.fail(tl, "bounds.templates must be [lower, upper]")
            lo = c.number(pair[0], tl, "bounds.templates lower", positive=True)
            hi = c.number(pair[1], tl, "bounds.templates upper", positive=True)
            if not lo < 1.0 < hi:
                c.fail(tl, "bounds.templates must satisfy lower < 1 < upper (template parameters start at 1)")
            tbounds = (lo, hi)

    def rel(key: str, default: str | None) -> Path | None:
        if key not in top:
            return None if default is None else (base / default)
        return base / c.string(top[key], top.line_of(key), key)

    return BuildConfig(
        path=path, models=models, compartments=comps, time_unit=time_unit, bounds_factor=factor, template_bounds=tbounds,
        library_size=size, library_seed=seed,
        vocabulary_path=rel("vocabulary", None), overrides_path=rel("overrides", None),
        cache_dir=rel("cache", ".plbuild-cache"), output_dir=rel("output", "out"),
    )


# --------------------------------------------------------------------------- vocabulary.yaml


@dataclass
class Entity:
    name: str
    cls: str
    unit: str
    molar_mass: float | None = None
    identifiers: list[str] = field(default_factory=list)
    synonyms: list[str] = field(default_factory=list)


def _load_vocabulary(path: Path) -> dict[str, Entity]:
    c = _Checker(path)
    top = load_yaml(path)
    if top is None:
        c.fail(1, "vocabulary is empty")
    top = c.mapping(top, 1, "vocabulary.yaml")
    out: dict[str, Entity] = {}
    for name, spec in top.items():
        line = top.line_of(name)
        c.ident(name, line, "entity name")
        spec = c.mapping(spec, line, f"entity '{name}'")
        c.keys(spec, {"class", "unit", "molar_mass", "identifiers", "synonyms"}, f"entity '{name}'",
               required={"class", "unit"}, line=line)
        cls = c.string(spec["class"], spec.line_of("class"), f"{name}.class")
        if cls not in ENTITY_CLASSES:
            c.fail(spec.line_of("class"), f"{name}.class '{cls}' must be one of {', '.join(ENTITY_CLASSES)}")
        mm = None
        if "molar_mass" in spec:
            mm = c.number(spec["molar_mass"], spec.line_of("molar_mass"), f"{name}.molar_mass", positive=True)
        out[name] = Entity(
            name, cls, c.string(spec["unit"], spec.line_of("unit"), f"{name}.unit"), mm,
            c.strings(spec["identifiers"], spec.line_of("identifiers"), f"{name}.identifiers") if "identifiers" in spec else [],
            c.strings(spec["synonyms"], spec.line_of("synonyms"), f"{name}.synonyms") if "synonyms" in spec else [],
        )
    return out


# --------------------------------------------------------------------------- overrides.yaml


@dataclass
class ModelOverrides:
    model_id: str
    units: dict[str, str] = field(default_factory=dict)        # substance, volume, time
    species: dict[str, str] = field(default_factory=dict)      # SBML species id -> entity name
    compartments: dict[str, str] = field(default_factory=dict) # SBML compartment id -> EnsGen compartment
    exclude: list[str] = field(default_factory=list)           # reaction ids
    bounds: dict[str, tuple[float, float]] = field(default_factory=dict)
    hill: list[str] = field(default_factory=list)              # parameters to treat as Hill coefficients


def _load_overrides(path: Path, vocab: dict[str, Entity], comps: dict[str, CompartmentSpec]) -> dict[str, ModelOverrides]:
    c = _Checker(path)
    top = load_yaml(path)
    if top is None:
        return {}
    top = c.mapping(top, 1, "overrides.yaml")
    out: dict[str, ModelOverrides] = {}
    for mid, spec in top.items():
        line = top.line_of(mid)
        if not isinstance(mid, str) or not BIOMODELS_ID.match(mid):
            c.fail(line, f"'{mid}' is not a BioModels id")
        spec = c.mapping(spec, line, f"overrides for {mid}")
        c.keys(spec, {"units", "species", "compartments", "exclude", "bounds", "hill"}, f"overrides for {mid}")
        o = ModelOverrides(mid)
        if "units" in spec:
            u = c.mapping(spec["units"], spec.line_of("units"), f"{mid}.units")
            c.keys(u, {"substance", "volume", "time"}, f"{mid}.units")
            o.units = {k: c.string(v, u.line_of(k), f"{mid}.units.{k}") for k, v in u.items()}
        if "species" in spec:
            s = c.mapping(spec["species"], spec.line_of("species"), f"{mid}.species")
            for sid, ent in s.items():
                c.string(ent, s.line_of(sid), f"{mid}.species.{sid}")
                if ent not in vocab:
                    c.fail(s.line_of(sid), f"species '{sid}' maps to '{ent}', which is not in the vocabulary")
                o.species[str(sid)] = ent
        if "compartments" in spec:
            cm = c.mapping(spec["compartments"], spec.line_of("compartments"), f"{mid}.compartments")
            for cid, comp in cm.items():
                c.string(comp, cm.line_of(cid), f"{mid}.compartments.{cid}")
                if comp not in comps:
                    c.fail(cm.line_of(cid), f"compartment '{cid}' maps to '{comp}', which is not declared in plbuild.yaml")
                o.compartments[str(cid)] = comp
        if "exclude" in spec:
            o.exclude = c.strings(spec["exclude"], spec.line_of("exclude"), f"{mid}.exclude")
        if "bounds" in spec:
            b = c.mapping(spec["bounds"], spec.line_of("bounds"), f"{mid}.bounds")
            for p, pair in b.items():
                pl = b.line_of(p)
                pair = c.seq(pair, pl, f"{mid}.bounds.{p}")
                if len(pair) != 2:
                    c.fail(pl, f"{mid}.bounds.{p} must be [lower, upper]")
                lo = c.number(pair[0], pl, f"{mid}.bounds.{p} lower", positive=True)
                hi = c.number(pair[1], pl, f"{mid}.bounds.{p} upper", positive=True)
                if not lo < hi:
                    c.fail(pl, f"{mid}.bounds.{p}: lower must be less than upper")
                o.bounds[str(p)] = (lo, hi)
        if "hill" in spec:
            o.hill = c.strings(spec["hill"], spec.line_of("hill"), f"{mid}.hill")
        out[mid] = o
    return out


# --------------------------------------------------------------------------- everything together


@dataclass
class Inputs:
    config: BuildConfig
    vocabulary: dict[str, Entity]
    overrides: dict[str, ModelOverrides]


def load_inputs(config_path: Path | str) -> Inputs:
    """Load plbuild.yaml and the vocabulary and overrides files it points to."""
    config = _load_build_config(Path(config_path))
    vocab = _load_vocabulary(config.vocabulary_path)
    overrides = _load_overrides(config.overrides_path, vocab, config.compartments) if config.overrides_path else {}
    return Inputs(config, vocab, overrides)
