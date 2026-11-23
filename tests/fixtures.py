"""Small hand-written models for the tests, written as Antimony and converted to SBML on the fly."""
from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import antimony

from plbuild.config import load_inputs

VOCAB = dedent("""\
    A: {class: molecule, unit: item/litre}
    B: {class: molecule, unit: item/litre}
    C: {class: molecule, unit: item/litre}
    E: {class: molecule, unit: item/litre}
    """)


def sbml_from_antimony(text: str, path: Path) -> Path:
    assert antimony.loadAntimonyString(text) >= 0, antimony.getLastError()
    path.write_text(antimony.getSBMLString(antimony.getMainModuleName()), encoding="utf-8")
    return path


def make_inputs(tmp_path: Path, model: str, *, overrides: str | None = None, compartments: str = "cell: {}",
                default: str = "cell", vocab: str = VOCAB, extra_models: str = "", library: str = ""):
    """One model, id BIOMD0000000001, in tmp_path. Returns (inputs, model entry)."""
    sbml_from_antimony(model, tmp_path / "m1.xml")
    if overrides is None:
        overrides = "BIOMD0000000001:\n  units: {substance: item, volume: litre, time: minute}\n"
    (tmp_path / "vocabulary.yaml").write_text(vocab)
    (tmp_path / "overrides.yaml").write_text(overrides)
    (tmp_path / "plbuild.yaml").write_text(dedent(f"""\
        vocabulary: vocabulary.yaml
        overrides: overrides.yaml
        output: out
        library:
          time_unit: minute
        {library}
        compartments:
          {compartments}
        models:
          - id: BIOMD0000000001
            sbml: m1.xml
            default_compartment: {default}
        {extra_models}
        """))
    inputs = load_inputs(tmp_path / "plbuild.yaml")
    return inputs, inputs.config.models[0]
