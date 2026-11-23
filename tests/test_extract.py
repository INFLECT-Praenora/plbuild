"""L2-style extraction: species roles, rate-law preparation, exclusion reason codes, compartment factor."""
import pytest

from plbuild.extract import extract_model
from fixtures import make_inputs


def extract(tmp_path, model, **kw):
    inputs, entry = make_inputs(tmp_path, model, **kw)
    return extract_model(entry, entry.sbml, inputs)


def only(res):
    assert not res.excluded, res.excluded
    assert len(res.processes) == 1
    return res.processes[0]


def code_of(res):
    assert len(res.excluded) == 1, (res.excluded, res.processes)
    return res.excluded[0].code


def test_reactant_product_and_rate_law_only_species(tmp_path):
    p = only(extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell, E in cell = 2\n"
                               "J1: A => B; k*A*E\nk = 0.5"))
    assert p.consumed == [("A_cell", 1)] and p.produced == [("B_cell", 1)] and p.influencing == ["E_cell"]
    assert p.expr == "k * A_cell * E_cell"
    assert [(q.name, q.value, q.lower, q.upper) for q in p.params] == [("k", 0.5, 0.005, 50.0)]


def test_enzyme_on_both_sides_is_influencing(tmp_path):
    p = only(extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell, E in cell = 2\n"
                               "J1: A + E => B + E; k*A*E\nk = 0.5"))
    assert p.consumed == [("A_cell", 1)] and p.produced == [("B_cell", 1)] and p.influencing == ["E_cell"]


def test_modifier_and_boundary_reactant_are_influencing(tmp_path):
    p = only(extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell, $E in cell = 2, C in cell = 1\n"
                               "J1: A + $E => B; k*A*E*C\nk = 0.5"))
    assert p.consumed == [("A_cell", 1)]
    assert sorted(p.influencing) == ["C_cell", "E_cell"]


def test_stoichiometry_is_kept(tmp_path):
    p = only(extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell\nJ1: 2 A => B; k*A^2\nk = 1"))
    assert p.consumed == [("A_cell", 2)]


def test_function_definitions_local_parameters_and_rules_are_resolved(tmp_path):
    model = """
        function mm(v, s, km)
          v*s/(km + s)
        end
        model m
          compartment cell = 1; species A in cell = 1, B in cell
          J1: A => B; mm(Vm, A, Km)*f
          Vm = 2; Km = 3
          f := g*2
          g := h + 1
          h = 4
        end
        """
    p = only(extract(tmp_path, model))
    assert "mm" not in p.expr and "f" not in p.expr.split()
    assert {q.name for q in p.params} == {"Vm", "Km", "h"}


def test_compartment_factor_is_removed(tmp_path):
    res = extract(tmp_path, "compartment cell = 2; species A in cell = 1, B in cell\nJ1: A => B; cell*k*A\nk = 1")
    assert only(res).expr == "k * A_cell"
    assert not res.conversions


def test_nested_product_compartment_factor_is_removed(tmp_path):
    res = extract(tmp_path, "compartment cell = 2; species A in cell = 1, B in cell, E in cell = 1\n"
                            "J1: A => B; k*cell*A*E\nk = 1")
    assert only(res).expr == "k * A_cell * E_cell"


def test_rate_without_compartment_factor_is_divided_and_reported(tmp_path):
    res = extract(tmp_path, "compartment cell = 2; species A in cell = 1, B in cell\nJ1: A => B; k*A\nk = 1")
    assert only(res).expr == "k * A_cell / 2"
    assert res.conversions[0]["note"] == "compartment_divided"


@pytest.mark.parametrize("model, code", [
    ("J1: A => B; k*A*time", "time_dependent"),
    ("J1: A => B; k*piecewise(1, A > 1, 0)", "time_dependent"),
    ("J1: A => B; k*sin(A)", "unsupported_function"),
    ("J1: A => B; k*A\nk' = 1", "dynamic_parameter"),
    ("J1: A => B; k*A\nat time > 1: k = 2", "dynamic_parameter"),
    ("J1: A => B; f*A\nf := k*time", "dynamic_rule"),
    ("J1: A => B; k*A\nk = 0", "parameter_value"),
    ("J1: A => Z; k*A\nspecies Z in cell", "unmapped_species"),
])
def test_exclusion_codes(tmp_path, model, code):
    base = "compartment cell = 1; species A in cell = 1, B in cell\nk = 1\n"
    assert code_of(extract(tmp_path, base + model)) == code


def test_unmapped_compartment(tmp_path):
    res = extract(tmp_path, "compartment cytoplasm = 1; compartment nucleus = 1\n"
                            "species A in cytoplasm = 1, B in cytoplasm\nJ1: A => B; k*A\nk = 1", default="cell")
    assert code_of(res) == "unmapped_compartment"
    assert {u["id"] for u in res.unresolved} == {"cytoplasm", "nucleus"}


def test_units_missing_is_excluded(tmp_path):
    res = extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell\nJ1: A => B; k*A\nk = 1",
                  overrides="")
    assert code_of(res) == "units"


def test_non_identity_units_wait_for_L4(tmp_path):
    res = extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell\nJ1: A => B; k*A\nk = 1",
                  overrides="BIOMD0000000001:\n  units: {substance: item, volume: litre, time: second}\n")
    assert code_of(res) == "units"
    assert "L4" in res.excluded[0].detail


def test_curator_excluded(tmp_path):
    res = extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell\nJ1: A => B; k*A\nk = 1",
                  overrides="BIOMD0000000001:\n  units: {substance: item, volume: litre, time: minute}\n"
                            "  exclude: [J1]\n")
    assert code_of(res) == "curator_excluded"


def test_hill_coefficient_bounds(tmp_path):
    p = only(extract(tmp_path, "compartment cell = 1; species A in cell = 1, B in cell\n"
                               "J1: A => B; V*A^n_A/(K^n_A + A^n_A)\nV = 1; K = 2; n_A = 6"))
    hill = next(q for q in p.params if q.name == "n_A")
    assert (hill.lower, hill.upper) == (1.0, 6.0)
