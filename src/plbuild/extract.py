"""Reading SBML and turning each reaction into a candidate curated process (spec sections 5, 7, 8, 9).

Each reaction is either kept as a Process (with species already mapped to library ids such as ERK_cell) or excluded
with one of the reason codes of spec section 5.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import libsbml
import pint

from .config import Inputs, ModelEntry, ModelOverrides
from .mapping import Mapper, annotations
from .process import Param, Process
from . import units

ALLOWED_FUNCS = {
    libsbml.AST_PLUS, libsbml.AST_MINUS, libsbml.AST_TIMES, libsbml.AST_DIVIDE, libsbml.AST_POWER,
    libsbml.AST_FUNCTION_POWER, libsbml.AST_FUNCTION_EXP, libsbml.AST_FUNCTION_LN, libsbml.AST_FUNCTION_LOG,
    libsbml.AST_FUNCTION_ROOT, libsbml.AST_FUNCTION_ABS,
    libsbml.AST_INTEGER, libsbml.AST_REAL, libsbml.AST_REAL_E, libsbml.AST_RATIONAL, libsbml.AST_NAME,
}
for _n in ("AST_FUNCTION_MAX", "AST_FUNCTION_MIN"):
    if hasattr(libsbml, _n):
        ALLOWED_FUNCS.add(getattr(libsbml, _n))
TIME_TYPES = {libsbml.AST_NAME_TIME, libsbml.AST_FUNCTION_DELAY, libsbml.AST_FUNCTION_PIECEWISE}
HILL_NAME = re.compile(r"^(n|h|nH)(_.*)?$")


@dataclass
class Exclusion:
    reaction: str
    code: str
    detail: str


@dataclass
class ModelResult:
    model_id: str
    name: str
    sbml_path: Path
    n_reactions: int = 0
    processes: list[Process] = field(default_factory=list)
    excluded: list[Exclusion] = field(default_factory=list)
    unresolved: list[dict] = field(default_factory=list)       # species and compartments that did not map
    conversions: list[dict] = field(default_factory=list)      # notes such as compartment_divided
    species_ids: dict[str, str] = field(default_factory=dict)  # SBML species id -> library species id
    # for the reconstruction test (values in source units)
    initial_conc: dict[str, float] = field(default_factory=dict)
    volumes: dict[str, float] = field(default_factory=dict)    # library compartment -> source volume
    source_params: dict[str, dict[str, float]] = field(default_factory=dict)  # pid -> {param: value}


class _Excluded(Exception):
    def __init__(self, code: str, detail: str):
        self.code, self.detail = code, detail


def _walk(ast: libsbml.ASTNode):
    yield ast
    for i in range(ast.getNumChildren()):
        yield from _walk(ast.getChild(i))


def _names(ast: libsbml.ASTNode) -> list[str]:
    seen: list[str] = []
    for n in _walk(ast):
        if n.getType() == libsbml.AST_NAME and n.getName() not in seen:
            seen.append(n.getName())
    return seen


def _rename(ast: libsbml.ASTNode, mapping: dict[str, str]):
    """Rename identifiers in two phases, so a -> b and b -> c cannot chain."""
    tmp = {old: f"__plb{i}__" for i, old in enumerate(mapping)}
    for old, t in tmp.items():
        ast.renameSIdRefs(old, t)
    for old, t in tmp.items():
        ast.renameSIdRefs(t, mapping[old])


def _expand(doc: libsbml.SBMLDocument):
    """Inline function definitions; replace initial assignments by their values (spec section 5, steps 1 and 4)."""
    for option in ("expandFunctionDefinitions", "expandInitialAssignments"):
        props = libsbml.ConversionProperties()
        props.addOption(option, True)
        doc.convert(props)


def _initial_conc(sp: libsbml.Species, model: libsbml.Model) -> float:
    if sp.isSetInitialConcentration():
        return sp.getInitialConcentration()
    vol = model.getCompartment(sp.getCompartment()).getSize()
    return sp.getInitialAmount() / vol if vol else sp.getInitialAmount()


def extract_model(entry: ModelEntry, sbml_path: Path, inputs: Inputs) -> ModelResult:
    doc = libsbml.readSBMLFromFile(str(sbml_path))
    if doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL) or doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) or doc.getModel() is None:
        raise RuntimeError(f"{sbml_path}: cannot read SBML: {doc.getErrorLog().toString().strip()}")
    _expand(doc)
    model = doc.getModel()
    ov: ModelOverrides | None = inputs.overrides.get(entry.id)
    cfg = inputs.config
    mapper = Mapper(inputs.vocabulary, cfg.compartments)
    res = ModelResult(entry.id, model.getName() or model.getId(), sbml_path, model.getNumReactions())

    # --- things a rate law may not depend on
    dynamic: set[str] = set()
    for r in model.getListOfRules():
        if r.isRate() or r.isAlgebraic():
            if r.isRate():
                dynamic.add(r.getVariable())
            else:
                dynamic.update(_names(r.getMath()))
    for ev in model.getListOfEvents():
        for ea in ev.getListOfEventAssignments():
            dynamic.add(ea.getVariable())
    assign = {r.getVariable(): r.getMath() for r in model.getListOfRules() if r.isAssignment()}

    # --- species and compartments
    comp_map: dict[str, str | None] = {}
    for c in model.getListOfCompartments():
        name, how = mapper.compartment(c, ov, model.getNumCompartments(), entry.default_compartment)
        comp_map[c.getId()] = name
        if name is None:
            res.unresolved.append({"kind": "compartment", "id": c.getId(), "name": c.getName(), "reason": how})
        else:
            res.volumes[name] = c.getSize() if c.isSetSize() else 1.0
    ent_map: dict[str, str | None] = {}
    for sp in model.getListOfSpecies():
        ent, how = mapper.species(sp, ov)
        ent_map[sp.getId()] = ent
        comp = comp_map.get(sp.getCompartment())
        if ent is None:
            res.unresolved.append({"kind": "species", "id": sp.getId(), "name": sp.getName(), "reason": how,
                                   "annotations": annotations(sp)})
        elif comp is not None:
            lib_id = f"{ent}_{comp}"
            res.species_ids[sp.getId()] = lib_id
            res.initial_conc[lib_id] = _initial_conc(sp, model)

    # --- units
    src_units, missing = units.source_units(model, ov.units if ov else {})

    for rxn in model.getListOfReactions():
        try:
            if ov and rxn.getId() in ov.exclude:
                raise _Excluded("curator_excluded", "listed in overrides.yaml")
            p = _reaction(rxn, model, ov, cfg, res, ent_map, comp_map, assign, dynamic, src_units, missing, inputs)
            res.processes.append(p)
        except _Excluded as e:
            res.excluded.append(Exclusion(rxn.getId(), e.code, e.detail))
    return res


def _reaction(rxn, model, ov, cfg, res, ent_map, comp_map, assign, dynamic, src_units, missing, inputs) -> Process:
    rid = rxn.getId()
    kl = rxn.getKineticLaw()
    if kl is None or kl.getMath() is None:
        raise _Excluded("no_rate_law", "the reaction has no kinetic law")
    original = libsbml.formulaToL3String(kl.getMath())
    ast = kl.getMath().deepCopy()

    # local parameters keep their own names unless they clash with a global symbol
    local: dict[str, float] = {}
    lp_list = kl.getListOfLocalParameters() if model.getLevel() >= 3 else kl.getListOfParameters()
    rename: dict[str, str] = {}
    for lp in lp_list:
        name = lp.getId()
        if model.getElementBySId(name) is not None:
            rename[name] = f"{rid}_{name}"
        local[rename.get(name, name)] = lp.getValue()
    if rename:
        _rename(ast, rename)

    # inline assignment rules (recursively)
    for _ in range(50):
        hits = [n for n in _names(ast) if n in assign and n not in local]
        if not hits:
            break
        for n in hits:
            rule = assign[n]
            bad = [m for m in _names(rule) if m in dynamic] or [1 for x in _walk(rule) if x.getType() in TIME_TYPES]
            if bad:
                raise _Excluded("dynamic_rule", f"assignment rule for '{n}' depends on time or a changing quantity")
            ast.replaceArgument(n, rule)
    else:
        raise _Excluded("dynamic_rule", "assignment rules are circular")

    for node in _walk(ast):
        t = node.getType()
        if t in TIME_TYPES:
            raise _Excluded("time_dependent", f"rate law uses {node.getName() or libsbml.formulaToL3String(node)}")
        if t not in ALLOWED_FUNCS and not node.isNumber() and t not in (libsbml.AST_CONSTANT_E, libsbml.AST_CONSTANT_PI):
            raise _Excluded("unsupported_function", f"rate law uses '{libsbml.formulaToL3String(node)}'")
    names = _names(ast)
    for n in names:
        if n in dynamic and n not in local:
            raise _Excluded("dynamic_parameter", f"'{n}' is changed by an event or a rule")

    # stoichiometry must be fixed
    for ref in list(rxn.getListOfReactants()) + list(rxn.getListOfProducts()):
        if (model.getLevel() < 3 and ref.isSetStoichiometryMath()) or (ref.isSetId() and (
                ref.getId() in dynamic or ref.getId() in assign)):
            raise _Excluded("stoichiometry", f"stoichiometry of {ref.getSpecies()} is set by a rule or event")

    # map species; net stoichiometry decides the roles
    def lib(sid: str) -> str:
        if ent_map.get(sid) is None:
            raise _Excluded("unmapped_species", f"species '{sid}' is not in the vocabulary")
        if comp_map.get(model.getSpecies(sid).getCompartment()) is None:
            raise _Excluded("unmapped_compartment",
                            f"compartment '{model.getSpecies(sid).getCompartment()}' of '{sid}' is not configured")
        return res.species_ids[sid]

    net: dict[str, float] = {}
    order: list[str] = []
    fixed: list[str] = []
    for sign, refs in ((-1, rxn.getListOfReactants()), (1, rxn.getListOfProducts())):
        for ref in refs:
            sp = model.getSpecies(ref.getSpecies())
            if sp.getBoundaryCondition() or sp.getConstant():
                fixed.append(sp.getId())
                continue
            net[sp.getId()] = net.get(sp.getId(), 0.0) + sign * ref.getStoichiometry()
            if sp.getId() not in order:
                order.append(sp.getId())
    consumed = [(lib(s), -net[s]) for s in order if net[s] < 0]
    produced = [(lib(s), net[s]) for s in order if net[s] > 0]
    moving = {s for s in order if net[s] != 0}
    infl_src: list[str] = []
    for s in [*(m.getSpecies() for m in rxn.getListOfModifiers()), *fixed, *order, *names]:
        if model.getSpecies(s) is not None and s not in moving and s not in infl_src:
            infl_src.append(s)
    influencing = [lib(s) for s in infl_src]
    if not consumed and not produced:
        raise _Excluded("stoichiometry", "the reaction changes no species")

    # species-valued names in the rate law must use amounts only when volume is 1
    for s in names:
        sp = model.getSpecies(s)
        if sp is not None and sp.getHasOnlySubstanceUnits():
            vol = model.getCompartment(sp.getCompartment()).getSize()
            if vol != 1:
                raise _Excluded("units", f"'{s}' is an amount (hasOnlySubstanceUnits) in a compartment of volume "
                                         f"{vol}; converting this arrives in L4")

    # units: only identity conversion for now
    if missing:
        raise _Excluded("units", f"model declares no {', '.join(missing)} units and overrides.yaml gives none")
    # reference species: first consumed, else first produced (EnsGen rate convention)
    ref_sid = next((s for s in order if net[s] < 0), None) or next(s for s in order if net[s] > 0)
    try:
        kt = units.time_factor(src_units["time"], cfg.time_unit)
        factors = {s: units.concentration_factor(src_units["substance"], src_units["volume"],
                                                 inputs.vocabulary[ent_map[s]].unit) for s in moving}
    except (pint.errors.PintError, pint.errors.DimensionalityError, AttributeError) as e:
        raise _Excluded("units", f"cannot convert units: {e}")
    if not math.isclose(kt, 1.0) or any(not math.isclose(f, 1.0) for f in factors.values()):
        raise _Excluded("units", f"needs conversion (time factor {kt:g}, concentration factors "
                                 f"{ {s: round(f, 6) for s, f in factors.items()} }); conversion arrives in L4")

    # compartment factor (spec section 8, step 1)
    ref_comp_src = model.getSpecies(ref_sid).getCompartment()
    removed = _remove_factor(ast, ref_comp_src)
    if removed and ast.getType() == libsbml.AST_TIMES and ast.getNumChildren() == 1:
        ast = ast.getChild(0).deepCopy()
    if not removed:
        vol = model.getCompartment(ref_comp_src).getSize()
        if vol != 1:
            ast = libsbml.parseL3Formula(f"({libsbml.formulaToL3String(ast)})/{vol!r}")
        res.conversions.append({"reaction": rid, "note": "compartment_divided",
                                "detail": f"divided by volume {vol:g}" + (" (no change to the text)" if vol == 1 else "")})
    other_comps = [n for n in _names(ast) if model.getCompartment(n) is not None]
    for c in other_comps:
        ast.replaceArgument(c, libsbml.parseL3Formula(repr(model.getCompartment(c).getSize())))

    # parameters
    params: list[Param] = []
    for n in _names(ast):
        if model.getSpecies(n) is not None:
            continue
        if n in local:
            value = local[n]
        else:
            g = model.getParameter(n)
            if g is None:
                raise _Excluded("unsupported_function", f"unknown symbol '{n}' in rate law")
            value = g.getValue()
        params.append(_bounds(rid, n, value, ov, cfg.bounds_factor))

    _rename(ast, {s: res.species_ids[s] for s in _names(ast) if model.getSpecies(s) is not None})
    expr = libsbml.formulaToL3String(ast)

    title = rxn.getName() or rid
    p = Process(
        pid="", kind="curated", title=title,
        consumed=consumed, produced=produced, influencing=influencing, expr=expr, params=params,
        compartments={s: _comp_of(s, res) for s in
                      [*(x for x, _ in consumed), *(x for x, _ in produced), *influencing]},
        identities=_identities(inputs, [x for x, _ in consumed] + [x for x, _ in produced] + influencing),
        sources=[(res.model_id, rid)], original_rate=original,
    )
    return p


def _remove_factor(ast: libsbml.ASTNode, name: str) -> bool:
    """Remove `name` if it is a factor of the top-level product. Products may be nested ((V*k)*S), so descend
    through TIMES nodes only; a factor inside a sum or quotient is not top level."""
    if ast.getType() != libsbml.AST_TIMES:
        return False
    for i in range(ast.getNumChildren()):
        ch = ast.getChild(i)
        if ch.getType() == libsbml.AST_NAME and ch.getName() == name:
            ast.removeChild(i)
            return True
        if ch.getType() == libsbml.AST_TIMES and _remove_factor(ch, name):
            if ch.getNumChildren() == 1:
                ast.replaceChild(i, ch.getChild(0).deepCopy(), True)
            return True
    return False


def _comp_of(lib_id: str, res: ModelResult) -> str:
    for comp in sorted(res.volumes, key=len, reverse=True):
        if lib_id.endswith("_" + comp):
            return comp
    raise AssertionError(lib_id)


def _identities(inputs: Inputs, lib_ids: list[str]) -> dict[str, str]:
    out = {}
    for s in lib_ids:
        for comp in inputs.config.compartments:
            if s.endswith("_" + comp):
                ent = inputs.vocabulary.get(s[: -len(comp) - 1])
                if ent and ent.identifiers:
                    out[s] = ent.identifiers[0]
    return out


def _bounds(rid: str, name: str, value: float, ov: ModelOverrides | None, factor: float) -> Param:
    """Spec section 9."""
    if ov and name in ov.bounds:
        lo, hi = ov.bounds[name]
        return Param(name, value, lo, hi, "bounds from overrides.yaml")
    if (ov and name in ov.hill) or HILL_NAME.match(name):
        lo, hi = min(1.0, value), max(4.0, value)
        if value <= 0:
            raise _Excluded("parameter_value", f"Hill coefficient '{name}' = {value}")
        return Param(name, value, lo, hi, "Hill coefficient")
    if value <= 0:
        raise _Excluded("parameter_value", f"parameter '{name}' = {value}; give bounds in overrides.yaml")
    return Param(name, value, value / factor, value * factor)
