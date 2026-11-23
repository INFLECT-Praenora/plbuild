"""Drawing the library from the candidate pool, and finding which dropped curated processes no template can replace.

PI decision (2026-11-23): the first library holds `library.size` entries drawn at random (seeded) from the curated
processes plus the templates, so EnsGen does not simply rebuild the source model. The report lists every curated
process that was dropped and cannot be rebuilt from the templates that were drawn ("gaps").

A dropped process is replaceable when one drawn base template, with its placeholders bound to the process's own
species, times zero or more linear multipliers (M_LIN) on its influencing species, gives the same rate law up to
parameter names. Hill coefficients are fixed at their source values first, so a Hill law with n = 1 counts as
Michaelis-Menten. Only M_LIN is tried as a multiplier; M_ACT and M_INH are not yet.
"""
from __future__ import annotations

import itertools
import random
from dataclasses import dataclass

import sympy

from .identity import canonical_form, parse, role_tokens
from .process import Process


@dataclass
class Replacement:
    pid: str
    title: str
    source: str
    replaced_by: str | None        # e.g. "A_MM x M_LIN(Mek_active_cell)", None for a gap
    candidates: list[str]          # every template combination that would work, drawn or not


def draw(pool: list[Process], size: int | None, seed: int) -> list[Process]:
    """Seeded draw of `size` entries. The pool is sorted by id first, so the draw depends only on content and seed."""
    ordered = sorted(pool, key=lambda p: p.pid)
    if size is None or size >= len(ordered):
        return ordered
    picked = random.Random(seed).sample(ordered, size)
    return sorted(picked, key=lambda p: p.pid)


def _hill_fixed(p: Process) -> tuple[sympy.Expr, list[str]]:
    expr = parse(p.expr)
    hill = {sympy.Symbol(q.name): sympy.Float(q.value) if q.value != int(q.value) else sympy.Integer(int(q.value))
            for q in p.params if q.meaning == "Hill coefficient"}
    return expr.subs(hill), [q.name for q in p.params if q.meaning != "Hill coefficient"]


def _bindings(t: Process, p: Process):
    """Ways to bind template t's placeholders to p's species, matching role pattern and stoichiometry."""
    if sorted(n for _, n in t.consumed) != sorted(n for _, n in p.consumed):
        return
    if sorted(n for _, n in t.produced) != sorted(n for _, n in p.produced):
        return
    if len(t.influencing) > len(p.influencing):
        return
    if t.pid in ("A_TRANS", "A_EXCH"):
        return   # need the same entity in two compartments; never true within one curated reaction here
    for cons in itertools.permutations(p.consumed):
        if [n for _, n in cons] != [n for _, n in t.consumed]:
            continue
        for prod in itertools.permutations(p.produced):
            if [n for _, n in prod] != [n for _, n in t.produced]:
                continue
            for infl in itertools.permutations(p.influencing, len(t.influencing)):
                bind = {**{a: b for (a, _), (b, _) in zip(t.consumed, cons)},
                        **{a: b for (a, _), (b, _) in zip(t.produced, prod)},
                        **dict(zip(t.influencing, infl))}
                yield bind


def replacements(p: Process, templates: list[Process]) -> list[str]:
    """Template combinations that reproduce curated process p (see module docstring)."""
    target, target_params = _hill_fixed(p)
    tokens = role_tokens([s for s, _ in sorted(p.consumed)], [s for s, _ in sorted(p.produced)], sorted(p.influencing))
    want = canonical_form(target, tokens, target_params)
    have_lin = any(t.pid == "M_LIN" for t in templates)
    out: list[str] = []
    for t in templates:
        if t.kind != "anonymous":
            continue
        for bind in _bindings(t, p):
            base = parse(t.expr).subs({sympy.Symbol(a): sympy.Symbol(b) for a, b in bind.items()}, simultaneous=True)
            free = [s for s in p.influencing if s not in bind.values()]
            # try multiplying by each influencing species 0, 1 or 2 times (M_LIN applied repeatedly)
            for powers in itertools.product(range(3), repeat=len(free)):
                if any(powers) and not have_lin:
                    continue
                expr = base * sympy.Mul(*[sympy.Symbol(s) ** k for s, k in zip(free, powers)])
                used = [s for s, k in zip(free, powers) if k]
                if set(used) | set(bind.values()) != set(s for s, _ in p.consumed) | set(
                        s for s, _ in p.produced) | set(p.influencing):
                    continue   # every species of p must take part
                if canonical_form(expr, tokens, [q.name for q in t.params]) == want:
                    mult = "".join(f" x M_LIN({s})" * k for s, k in zip(free, powers))
                    out.append(f"{t.pid}{mult}")
    return out


def find_gaps(curated: list[Process], drawn: list[Process], templates: list[Process]) -> list[Replacement]:
    """For each curated process not drawn: which drawn templates replace it, if any."""
    drawn_ids = {p.pid for p in drawn}
    out = []
    for p in sorted(curated, key=lambda x: x.pid):
        if p.pid in drawn_ids:
            continue
        all_ways = replacements(p, templates)
        ok = [w for w in all_ways if all(part.split("(")[0].strip() in drawn_ids for part in w.split(" x "))]
        source = f"{p.sources[0][0]}:{p.sources[0][1]}" if p.sources else ""
        out.append(Replacement(p.pid, p.title, source, ok[0] if ok else None, all_ways))
    return out
