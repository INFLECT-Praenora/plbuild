"""Source units and conversion factors (spec section 8).

For now only identity conversions are carried out: a reaction whose factors are not all 1 is excluded with `units`
and a note that general conversion arrives in milestone L4. The factors themselves are computed with pint.
"""
from __future__ import annotations

from functools import lru_cache

import libsbml
import pint


@lru_cache(maxsize=1)
def registry() -> pint.UnitRegistry:
    ureg = pint.UnitRegistry()
    ureg.define("item = [item]")          # SBML's 'item': a count of entities
    ureg.define("cell = [cell]")
    return ureg


def _unitdef_to_text(ud: libsbml.UnitDefinition) -> str:
    parts = []
    for u in ud.getListOfUnits():
        kind = libsbml.UnitKind_toString(u.getKind())
        factor = u.getMultiplier() * 10 ** u.getScale()
        base = kind if factor == 1 else f"({factor!r}*{kind})"
        exp = u.getExponentAsDouble()
        parts.append(base if exp == 1 else f"{base}**{exp:g}")
    return "*".join(parts) if parts else "dimensionless"


def source_units(model: libsbml.Model, overrides_units: dict[str, str]) -> tuple[dict[str, str], list[str]]:
    """Substance, volume and time units of a model, as pint-readable text. Returns (units, missing keys)."""
    out: dict[str, str] = {}
    attr = {"substance": model.getSubstanceUnits() if model.isSetSubstanceUnits() else "",
            "volume": model.getVolumeUnits() if model.isSetVolumeUnits() else "",
            "time": model.getTimeUnits() if model.isSetTimeUnits() else ""}
    for key in ("substance", "volume", "time"):
        if key in overrides_units:
            out[key] = overrides_units[key]
            continue
        uid = attr[key] or key                      # Level 2 uses unit definitions named substance/volume/time
        ud = model.getUnitDefinition(uid)
        if ud is not None:
            out[key] = _unitdef_to_text(ud)
        elif attr[key] and libsbml.UnitKind_isValidUnitKindString(attr[key], model.getLevel(), model.getVersion()):
            out[key] = attr[key]
    return out, [k for k in ("substance", "volume", "time") if k not in out]


def concentration_factor(substance: str, volume: str, entity_unit: str) -> float:
    """s_X: library concentration = s_X * source concentration. Raises pint errors if dimensions differ."""
    ureg = registry()
    q = (1 * ureg(substance) / ureg(volume)).to(ureg(entity_unit))
    return float(q.magnitude)


def time_factor(source_time: str, library_time: str) -> float:
    """k_t: number of source time units per library time unit."""
    ureg = registry()
    return float((1 * ureg(library_time)).to(ureg(source_time)).magnitude)
