"""Identity keys and process ids (spec section 10), and the rate-law canonical form they rest on.

The canonical form replaces each species by a role token (C0.. consumed, P0.. produced, I0.. influencing) and each
parameter by a positional token, chosen so that renaming parameters never changes the form. `k*X` and `kdeg*X` share a
form; `Vm*X/(Km + X)` and `k*X` do not.
"""
from __future__ import annotations

import hashlib
import itertools
import re

import sympy

from .process import Process

_FUNCS = {"exp": sympy.exp, "ln": sympy.log, "log": sympy.log, "sqrt": sympy.sqrt, "abs": sympy.Abs,
          "pow": sympy.Pow, "min": sympy.Min, "max": sympy.Max}
_IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
MAX_PERMUTED = 6    # parameters up to this count are canonicalized exactly over all orderings


def parse(expr: str) -> sympy.Expr:
    """Parse an SBML L3 / Antimony rate expression; every identifier becomes a plain symbol."""
    names = {n for n in _IDENT.findall(expr) if n not in _FUNCS}
    local = {n: sympy.Symbol(n) for n in names}
    local.update(_FUNCS)
    return sympy.sympify(expr.replace("^", "**"), locals=local)


def _normal(e: sympy.Expr) -> sympy.Expr:
    """A normal form that renaming symbols cannot change: one fraction, numerator and denominator expanded."""
    try:
        return sympy.cancel(sympy.together(e))
    except (sympy.PolynomialError, TypeError):
        return sympy.simplify(e)


def role_tokens(consumed: list[str], produced: list[str], influencing: list[str]) -> dict[str, sympy.Symbol]:
    tok: dict[str, sympy.Symbol] = {}
    for prefix, species in (("C", consumed), ("P", produced), ("I", influencing)):
        for i, s in enumerate(species):
            tok.setdefault(s, sympy.Symbol(f"{prefix}{i}"))
    return tok


def canonical_form(expr: sympy.Expr, tokens: dict[str, sympy.Symbol], params: list[str]) -> str:
    """The expression with species as role tokens and parameters as positional tokens; invariant to parameter names."""
    e = expr.subs({sympy.Symbol(s): t for s, t in tokens.items()}, simultaneous=True)
    e = _normal(e)
    syms = [sympy.Symbol(p) for p in params if sympy.Symbol(p) in e.free_symbols]
    if not syms:
        return sympy.srepr(e)
    if len(syms) > MAX_PERMUTED:
        perms = [syms]
    else:
        perms = itertools.permutations(syms)
    best = None
    for perm in perms:
        cand = sympy.srepr(e.subs({s: sympy.Symbol(f"p{i}") for i, s in enumerate(perm)}, simultaneous=True))
        if best is None or cand < best:
            best = cand
    return best


def identity_key(p: Process) -> tuple:
    """Spec section 10: species with roles and stoichiometry, plus the canonical rate-law form."""
    consumed = sorted(p.consumed)
    produced = sorted(p.produced)
    influencing = sorted(p.influencing)
    tokens = role_tokens([s for s, _ in consumed], [s for s, _ in produced], influencing)
    form = canonical_form(parse(p.expr), tokens, [q.name for q in p.params])
    return (tuple(consumed), tuple(produced), tuple(influencing), form)


def process_id(key: tuple) -> str:
    return "P_" + hashlib.sha256(repr(key).encode("utf-8")).hexdigest()[:10]


def assign_ids(processes: list[Process]) -> None:
    for p in processes:
        p.key = identity_key(p)
        p.pid = process_id(p.key)
