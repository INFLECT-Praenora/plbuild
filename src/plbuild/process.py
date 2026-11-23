"""The Process record shared by every stage, and how a process is written as a library entry."""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Param:
    name: str
    value: float
    lower: float
    upper: float
    meaning: str = ""          # one-line explanation for the comment block


@dataclass
class Process:
    pid: str                                  # P_<hash> for curated, A_*/M_* for templates
    kind: str                                 # curated | anonymous | modifier
    title: str                                # one line: what the process does
    consumed: list[tuple[str, float]]         # (species id, stoichiometry)
    produced: list[tuple[str, float]]
    influencing: list[str]
    expr: str                                 # rate law, concentration per time
    params: list[Param]
    compartments: dict[str, str] = field(default_factory=dict)   # species id -> compartment
    identities: dict[str, str] = field(default_factory=dict)     # species id -> identifiers.org URI
    sources: list[tuple[str, str]] = field(default_factory=list) # (model id, reaction id)
    original_rate: str = ""                   # the SBML kinetic law, for the comment
    notes: list[str] = field(default_factory=list)               # extra comment lines
    key: tuple = ()                           # identity key (identity.py)


def _term(species: str, stoich: float) -> str:
    if stoich == 1:
        return species
    s = f"{stoich:g}"
    return f"{s} {species}"


def reaction_line(p: Process) -> str:
    lhs = " + ".join(_term(s, n) for s, n in p.consumed)
    rhs = " + ".join(_term(s, n) for s, n in p.produced)
    return f"J: {lhs} => {rhs}; {p.expr}".replace("J:  =>", "J: =>").replace("=>  ;", "=> ;")


def fmt(x: float) -> str:
    """Shortest text that reads back as the same double, so rebuilds are byte-identical."""
    return repr(float(x)) if x != int(x) or abs(x) >= 1e16 else str(int(x))


def fmt_bound(x: float) -> str:
    """Bounds are soft limits, so 6 significant digits keep them readable."""
    return f"{x:.6g}"


def entry_text(p: Process) -> str:
    """The process as one commented library entry (EnsGen spec section 5 format)."""
    bar = "// " + "-" * 76
    out = [bar, f"// {p.title}"]
    for model, rxn in p.sources:
        out.append(f"// From {model}, reaction {rxn}")
    roles = []
    if p.consumed:
        roles.append("consumes " + ", ".join(_term(s, n) for s, n in p.consumed))
    if p.produced:
        roles.append("produces " + ", ".join(_term(s, n) for s, n in p.produced))
    if p.influencing:
        roles.append("influenced by " + ", ".join(p.influencing))
    if roles:
        out.append("// " + "; ".join(roles))
    if p.original_rate:
        out.append(f"// Original rate law: {p.original_rate}")
    for q in p.params:
        line = f"//   {q.name} = {fmt(q.value)}, allowed {fmt_bound(q.lower)} to {fmt_bound(q.upper)}"
        out.append(line + (f"  ({q.meaning})" if q.meaning else ""))
    out += [f"// {n}" for n in p.notes]
    out.append(bar)
    out.append(f"// @process {p.pid}")
    out.append(f"// @kind {p.kind}")
    out += [f"// @source {m}:{r}" for m, r in p.sources]
    for s in sorted(p.compartments):
        out.append(f"{s} in {p.compartments[s]}")
    for s in sorted(p.identities):
        out.append(f'{s} identity "{p.identities[s]}"')
    out.append(reaction_line(p))
    for q in p.params:
        out.append(f"{q.name} = {fmt(q.value)}")
    for q in p.params:
        out.append(f"// @bounds {q.name} {fmt_bound(q.lower)} {fmt_bound(q.upper)}")
    return "\n".join(out) + "\n"
