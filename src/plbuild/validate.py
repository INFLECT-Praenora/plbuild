"""Validation of the written library (spec section 12).

Checks run on the files as written, read back with a small reader of the restricted format:

- every entry loads in libAntimony on its own;
- every parameter has positive bounds with lower <= value <= upper (a Hill coefficient of 1 sits on its bound);
- every curated species id is a vocabulary entity plus a configured compartment;
- reconstruction: each source model, rebuilt from all of its entries in pool.ant, simulates like the original.

The EnsGen parser check (`ensgen_harness --check-library`) is added when EnsGen M2 exists.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import antimony
import numpy as np
import roadrunner

from .config import Inputs

_STATEMENT = re.compile(r"^(?P<lhs>[A-Za-z_]\w*)\s*=\s*(?P<val>\S+)\s*$")
_IN = re.compile(r"^(?P<s>[A-Za-z_]\w*)\s+in\s+(?P<c>[A-Za-z_]\w*)\s*$")
_RXN = re.compile(r"^J:\s*(?P<lhs>.*?)\s*=>\s*(?P<rhs>.*?)\s*;\s*(?P<rate>.+)$")


@dataclass
class Entry:
    pid: str
    kind: str = ""
    sources: list[tuple[str, str]] = field(default_factory=list)
    compartments: dict[str, str] = field(default_factory=dict)
    consumed: list[tuple[str, float]] = field(default_factory=list)
    produced: list[tuple[str, float]] = field(default_factory=list)
    rate: str = ""
    params: dict[str, float] = field(default_factory=dict)
    bounds: dict[str, tuple[float, float]] = field(default_factory=dict)
    text: str = ""
    line: int = 0


def _side(text: str) -> list[tuple[str, float]]:
    out = []
    for term in filter(None, (t.strip() for t in text.split("+"))):
        parts = term.split()
        out.append((parts[-1], float(parts[0]) if len(parts) == 2 else 1.0))
    return out


def read_library(path: Path) -> list[Entry]:
    """Read a library file back into entries. Raises ValueError with file:line on anything unexpected."""
    entries: list[Entry] = []
    cur: Entry | None = None
    for no, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if line.startswith("// @process"):
            cur = Entry(line.split()[2], line=no)
            entries.append(cur)
        if cur is not None:
            cur.text += raw + "\n"
        if not line or (line.startswith("//") and not line.startswith("// @")):
            continue
        if cur is None:
            raise ValueError(f"{path}:{no}: statement before the first // @process")
        if line.startswith("// @kind"):
            cur.kind = line.split()[2]
        elif line.startswith("// @source"):
            m, r = line.split()[2].split(":", 1)
            cur.sources.append((m, r))
        elif line.startswith("// @bounds"):
            _, _, name, lo, hi = line.split()
            cur.bounds[name] = (float(lo), float(hi))
        elif line.startswith("// @process"):
            pass
        elif (m := _IN.match(line)):
            cur.compartments[m["s"]] = m["c"]
        elif " identity " in line:
            pass
        elif (m := _RXN.match(line)):
            cur.consumed, cur.produced, cur.rate = _side(m["lhs"]), _side(m["rhs"]), m["rate"]
        elif (m := _STATEMENT.match(line)):
            cur.params[m["lhs"]] = float(m["val"])
        else:
            raise ValueError(f"{path}:{no}: not a statement of the library format: {line}")
    return entries


@dataclass
class Check:
    name: str
    ok: bool
    messages: list[str] = field(default_factory=list)


def check_antimony(entries: list[Entry]) -> Check:
    bad = []
    for e in entries:
        if antimony.loadAntimonyString(e.text) < 0:
            bad.append(f"{e.pid} (line {e.line}): {antimony.getLastError().strip()}")
    return Check("antimony", not bad, bad)


def check_bounds(entries: list[Entry]) -> Check:
    bad = []
    for e in entries:
        for name, value in e.params.items():
            if name not in e.bounds:
                bad.append(f"{e.pid}: {name} has no bounds")
                continue
            lo, hi = e.bounds[name]
            # spec 12 says lower < value < upper, but its own Hill bounds [1, 4] put n = 1 on the lower bound
            if not (0 < lo < hi) or not (lo <= value <= hi):
                bad.append(f"{e.pid}: {name} = {value} with bounds [{lo}, {hi}]")
    return Check("bounds", not bad, bad)


def check_vocabulary(entries: list[Entry], inputs: Inputs) -> Check:
    valid = {f"{ent}_{comp}" for ent in inputs.vocabulary for comp in inputs.config.compartments}
    bad = []
    for e in entries:
        if e.kind != "curated":
            continue
        for s, c in e.compartments.items():
            if s not in valid or c not in inputs.config.compartments:
                bad.append(f"{e.pid}: species '{s}' in '{c}' is not a vocabulary entity plus a configured compartment")
    return Check("vocabulary", not bad, bad)


def _rebuild_antimony(entries: list[Entry], initial: dict[str, float], volumes: dict[str, float]) -> str:
    """Spec 12, step 2: library rates are concentration per time, so each is multiplied by its compartment volume."""
    lines = ["model rebuilt"]
    for c, v in sorted(volumes.items()):
        lines.append(f"  compartment {c} = {v!r}")
    species: dict[str, str] = {}
    for k, e in enumerate(entries):
        species.update(e.compartments)
        rate = e.rate
        for p in sorted(e.params, key=len, reverse=True):
            rate = re.sub(rf"\b{re.escape(p)}\b", f"e{k}_{p}", rate)
        ref = (e.consumed or e.produced)[0][0]
        lhs = " + ".join(s if n == 1 else f"{n:g} {s}" for s, n in e.consumed)
        rhs = " + ".join(s if n == 1 else f"{n:g} {s}" for s, n in e.produced)
        lines.append(f"  R{k}: {lhs} => {rhs}; {e.compartments[ref]}*({rate})")
        for p, v in e.params.items():
            lines.append(f"  e{k}_{p} = {v!r}")
    for s, c in sorted(species.items()):
        lines.append(f"  species {s} in {c} = {initial.get(s, 0.0)!r}")
    lines.append("end")
    return "\n".join(lines)


def reconstruct(entries: list[Entry], model_id: str, sbml_path: Path, species_ids: dict[str, str],
                initial: dict[str, float], volumes: dict[str, float], t_end: float, n_points: int = 1001,
                rtol: float = 1e-6, floor: float = 1e-12) -> Check:
    """Spec 12, reconstruction: pass if every species agrees within rtol relative plus floor x the species maximum."""
    mine = [e for e in entries if any(m == model_id for m, _ in e.sources)]
    original = roadrunner.RoadRunner(str(sbml_path))
    if antimony.loadAntimonyString(_rebuild_antimony(mine, initial, volumes)) < 0:
        return Check(f"reconstruction {model_id}", False, [antimony.getLastError().strip()])
    rebuilt = roadrunner.RoadRunner(antimony.getSBMLString("rebuilt"))
    for rr in (original, rebuilt):
        rr.integrator.relative_tolerance = 1e-12
        rr.integrator.absolute_tolerance = 1e-12
        rr.integrator.maximum_num_steps = 100000
    src_ids = [s for s in species_ids if species_ids[s] in rebuilt.model.getFloatingSpeciesIds()
               + rebuilt.model.getBoundarySpeciesIds()]
    original.timeCourseSelections = ["time"] + [f"[{s}]" for s in src_ids]
    rebuilt.timeCourseSelections = ["time"] + [f"[{species_ids[s]}]" for s in src_ids]
    a = np.array(original.simulate(0, t_end, n_points))
    b = np.array(rebuilt.simulate(0, t_end, n_points))
    for j, s in enumerate(src_ids, start=1):
        ya, yb = a[:, j], b[:, j]
        tol = rtol * np.abs(ya) + floor * np.max(np.abs(ya))
        bad = np.nonzero(np.abs(ya - yb) > tol)[0]
        if bad.size:
            i = bad[0]
            return Check(f"reconstruction {model_id}", False,
                         [f"{s} differs first at t = {a[i, 0]:g}: original {ya[i]:.10g}, rebuilt {yb[i]:.10g}"])
    missing = sorted(set(species_ids) - set(src_ids))
    msg = [f"{len(src_ids)} species agree over 0..{t_end:g}"]
    if missing:
        msg.append(f"not in the rebuilt model (used by no kept reaction): {', '.join(missing)}")
    return Check(f"reconstruction {model_id}", True, msg)
