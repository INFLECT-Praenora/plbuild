"""Identity keys, the seeded draw and the replacement (gap) check."""
import antimony

from plbuild.identity import identity_key, process_id
from plbuild.process import Param, Process, entry_text
from plbuild.select import draw, find_gaps, replacements
from plbuild.templates import templates


def proc(expr, params, consumed=(("X_cell", 1),), produced=(), influencing=(), hill=()):
    return Process(pid="", kind="curated", title="t", consumed=list(consumed), produced=list(produced),
                   influencing=list(influencing), expr=expr,
                   params=[Param(n, v, v / 100, v * 100, "Hill coefficient" if n in hill else "") for n, v in params])


def test_same_form_same_key():
    assert identity_key(proc("k*X_cell", [("k", 1)])) == identity_key(proc("kdeg*X_cell", [("kdeg", 2)]))


def test_different_form_different_key():
    assert identity_key(proc("Vm*X_cell/(Km + X_cell)", [("Vm", 1), ("Km", 1)])) != \
        identity_key(proc("k*X_cell", [("k", 1)]))


def test_key_ignores_parameter_names_and_order():
    a = proc("Vm*X_cell/(Km + X_cell)", [("Vm", 1), ("Km", 1)])
    b = proc("V*X_cell/(X_cell + K)", [("K", 1), ("V", 1)])
    assert identity_key(a) == identity_key(b)
    assert process_id(identity_key(a)) == process_id(identity_key(b))


def test_templates_are_valid_entries():
    ts = templates()
    assert len(ts) == 16 and len({t.pid for t in ts}) == 16
    for t in ts:
        assert antimony.loadAntimonyString(entry_text(t)) >= 0, (t.pid, antimony.getLastError())
        for q in t.params:
            assert 0 < q.lower <= q.value <= q.upper


def test_draw_is_seeded_and_sized():
    pool = templates()
    a = draw(pool, 5, seed=3)
    assert [p.pid for p in a] == [p.pid for p in draw(list(reversed(pool)), 5, seed=3)]
    assert len(a) == 5
    assert [p.pid for p in draw(pool, None, 0)] == sorted(p.pid for p in pool)


def test_enzyme_hill_law_is_replaced_by_mm_times_linear():
    p = proc("E_cell*k*X_cell^n/(KM^n + X_cell^n)", [("k", 1), ("KM", 5), ("n", 1)], produced=[("Y_cell", 1)],
             influencing=["E_cell"], hill=("n",))
    assert replacements(p, templates()) == ["A_MM x M_LIN(E_cell)"]


def test_hill_law_with_n_2_has_no_template():
    p = proc("E_cell*k*X_cell^n/(KM^n + X_cell^n)", [("k", 1), ("KM", 5), ("n", 2)], produced=[("Y_cell", 1)],
             influencing=["E_cell"], hill=("n",))
    assert replacements(p, templates()) == []


def test_gap_when_the_needed_template_is_not_drawn():
    p = proc("kd*X_cell", [("kd", 1)])
    p.pid = "P_x"
    ts = templates()
    with_deg = find_gaps([p], [t for t in ts if t.pid == "A_DEG"], ts)
    assert with_deg[0].replaced_by == "A_DEG"
    without = find_gaps([p], [t for t in ts if t.pid == "A_SYN"], ts)
    assert without[0].replaced_by is None and without[0].candidates == ["A_DEG"]
