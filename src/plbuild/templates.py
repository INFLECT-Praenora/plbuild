"""Generic template processes ("anonymous" in the EnsGen spec), defined here in code.

PI decision (2026-11-23): EnsGen has no built-in library. plbuild defines the templates and writes them into the same
library file as the curated processes. The set is the 14 from EnsGen spec section 5 plus 2X -> Y and X + Y -> Z.

Species names in a template are placeholders: EnsGen binds them to real species when it inserts the process.
Reactants become consumed roles, products produced roles, and any other species in the rate law an influencing role.
"""
from __future__ import annotations

from .process import Param, Process

LO, HI = 1e-3, 1e3   # default bounds for template parameters


def _p(name: str, meaning: str, value: float = 1.0, lower: float = LO, upper: float = HI) -> Param:
    return Param(name, value, lower, upper, meaning)


def _t(pid, kind, title, consumed, produced, influencing, expr, params, notes=()):
    return Process(pid=pid, kind=kind, title=title, consumed=consumed, produced=produced, influencing=influencing,
                   expr=expr, params=params, notes=list(notes))


def templates() -> list[Process]:
    """All template processes, in a fixed order."""
    return [
        _t("A_SYN", "anonymous", "Template: synthesis of X at a constant rate", [], [("X", 1)], [],
           "k0", [_p("k0", "synthesis rate")]),
        _t("A_DEG", "anonymous", "Template: first-order degradation or clearance of X", [("X", 1)], [], [],
           "kd*X", [_p("kd", "degradation rate constant")]),
        _t("A_CONV", "anonymous", "Template: first-order conversion of X into Y (e.g. activation, differentiation)",
           [("X", 1)], [("Y", 1)], [], "k*X", [_p("k", "conversion rate constant")]),
        _t("A_REVCONV", "anonymous", "Template: reversible first-order conversion between X and Y",
           [("X", 1)], [("Y", 1)], [], "k1*X - k2*Y",
           [_p("k1", "forward rate constant"), _p("k2", "backward rate constant")]),
        _t("A_MM", "anonymous", "Template: saturating (Michaelis-Menten) conversion of S into P",
           [("S", 1)], [("P", 1)], [], "Vm*S/(Km + S)",
           [_p("Vm", "maximal rate"), _p("Km", "half-saturation constant")]),
        _t("A_MMREV", "anonymous", "Template: reversible saturating conversion between S and P",
           [("S", 1)], [("P", 1)], [], "Vm*(S - P/Keq)/(Km*(1 + S/Km + P/Kmp))",
           [_p("Vm", "maximal rate"), _p("Km", "half-saturation constant for S"),
            _p("Keq", "equilibrium constant"), _p("Kmp", "half-saturation constant for P")]),
        _t("A_PROLIF", "anonymous", "Template: proliferation, X makes more X", [("X", 1)], [("X", 2)], [],
           "k*X", [_p("k", "proliferation rate constant")]),
        _t("A_LOGIST", "anonymous", "Template: logistic proliferation of X up to a capacity Xmax",
           [("X", 1)], [("X", 2)], [], "k*X*(1 - X/Xmax)",
           [_p("k", "proliferation rate constant"), _p("Xmax", "carrying capacity")]),
        _t("A_KILL", "anonymous", "Template: Y removes X (bimolecular), Y is not used up",
           [("X", 1)], [], ["Y"], "k*X*Y", [_p("k", "removal rate constant")],
           notes=["Written X => with Y influencing; the spec writes X + Y => Y, which is the same net change."]),
        _t("A_TRANS", "anonymous", "Template: directed transport of the same entity from compartment a to b",
           [("X_a", 1)], [("X_b", 1)], [], "k*X_a", [_p("k", "transport rate constant")],
           notes=["EnsGen binds X_a and X_b to the same entity in two different compartments."]),
        _t("A_EXCH", "anonymous", "Template: exchange of the same entity between compartments a and b",
           [("X_a", 1)], [("X_b", 1)], [], "k*(X_a - X_b)", [_p("k", "exchange rate constant")],
           notes=["EnsGen binds X_a and X_b to the same entity in two different compartments."]),
        _t("A_DIM", "anonymous", "Template: two X combine into Y (e.g. receptor dimerisation)",
           [("X", 2)], [("Y", 1)], [], "k*X^2", [_p("k", "association rate constant")]),
        _t("A_BIND", "anonymous", "Template: X and Y combine into Z (e.g. complex formation)",
           [("X", 1), ("Y", 1)], [("Z", 1)], [], "k*X*Y", [_p("k", "association rate constant")]),
        _t("M_LIN", "modifier", "Multiplier: the rate is proportional to A (e.g. A is an enzyme)",
           [], [], ["A"], "A", []),
        _t("M_ACT", "modifier", "Multiplier: A activates, with saturation (Hill form)",
           [], [], ["A"], "A^h/(K^h + A^h)",
           [_p("K", "half-activation constant"), _p("h", "Hill coefficient", 1.0, 1.0, 4.0)]),
        _t("M_INH", "modifier", "Multiplier: I inhibits (Hill form)",
           [], [], ["I"], "1/(1 + (I/K)^h)",
           [_p("K", "half-inhibition constant"), _p("h", "Hill coefficient", 1.0, 1.0, 4.0)]),
    ]
