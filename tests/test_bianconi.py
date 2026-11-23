"""End to end on BIOMD0000000427 (Bianconi 2012): the first real library."""
import json
import shutil
from pathlib import Path

import pytest

from plbuild.build import build, validate
from plbuild.config import load_inputs
from plbuild.validate import read_library, reconstruct

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "bianconi2012"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("b427")
    for f in ("plbuild.yaml", "vocabulary.yaml", "overrides.yaml", "BIOMD0000000427.xml"):
        shutil.copy(EXAMPLE / f, d / f)
    inputs = load_inputs(d / "plbuild.yaml")
    return d, inputs, build(inputs)


def test_all_reactions_kept(built):
    _, _, res = built
    r = res.results[0]
    assert (r.n_reactions, len(r.processes), r.excluded, r.unresolved) == (22, 22, [], [])


def test_library_is_25_of_38(built):
    d, _, res = built
    assert len(res.pool) == 38 and len(res.drawn) == 25
    assert len(read_library(d / "out" / "curated.ant")) == 25
    assert len(read_library(d / "out" / "pool.ant")) == 38


def test_validation_passes(built):
    _, inputs, _ = built
    checks = validate(inputs)
    assert all(c.ok for c in checks), [(c.name, c.messages) for c in checks if not c.ok]
    assert any(c.name.startswith("reconstruction") for c in checks)


def test_report_lists_dropped_processes_and_gaps(built):
    d, _, res = built
    rep = json.loads((d / "out" / "report.json").read_text())
    dropped = rep["selection"]["dropped_curated"]
    assert len(dropped) == 22 - sum(p.kind == "curated" for p in res.drawn)
    for g in dropped:
        assert g["possible_replacements"], g           # every 427 process can be rebuilt from some template set
    assert rep["totals"]["gaps"] == sum(g["replaced_by"] is None for g in dropped)


def test_rebuild_is_byte_identical(built):
    d, inputs, _ = built
    first = (d / "out" / "curated.ant").read_bytes()
    build(inputs)
    assert (d / "out" / "curated.ant").read_bytes() == first


def test_broken_conversion_fails_reconstruction(built):
    d, inputs, _ = built
    pool = read_library(d / "out" / "pool.ant")
    victim = next(e for e in pool if e.kind == "curated" and "degradation" in e.text)
    k = next(iter(victim.params))
    victim.params[k] *= 10                             # a conversion factor off by 10
    info = json.loads((d / "out" / "reconstruction.json").read_text())["BIOMD0000000427"]
    check = reconstruct(pool, "BIOMD0000000427", Path(info["sbml"]), info["species_ids"], info["initial"],
                        info["volumes"], 2000)
    assert not check.ok
    assert "differs first at t =" in check.messages[0]
