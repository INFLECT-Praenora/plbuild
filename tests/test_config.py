"""L0: malformed plbuild.yaml, vocabulary.yaml and overrides.yaml each fail with the file and line."""
from pathlib import Path
from textwrap import dedent

import pytest

from plbuild.config import load_inputs
from plbuild.yamlload import ConfigError

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "bianconi2012" / "plbuild.yaml"

GOOD_CONFIG = dedent("""\
    vocabulary: vocabulary.yaml
    overrides: overrides.yaml
    compartments:
      blood:
        synonyms: [plasma]
      tissue: {}
    models:
      - id: BIOMD0000000427
        version: 5
        default_compartment: tissue
      - BIOMD0000000001
    """)
GOOD_VOCAB = dedent("""\
    IL6:
      class: molecule
      unit: pg/mL
      molar_mass: 20900
      identifiers: [http://identifiers.org/uniprot/P05231]
      synonyms: [IL-6, interleukin 6]
    CD4_T:
      class: cell
      unit: cells/uL
    """)
GOOD_OVERRIDES = dedent("""\
    BIOMD0000000427:
      units: {substance: item, volume: litre, time: minute}
      species: {s1: IL6}
      compartments: {c1: blood}
      exclude: [J1]
      bounds: {k1: [0.1, 10]}
      hill: [n1]
    """)


def write(tmp_path, config=GOOD_CONFIG, vocab=GOOD_VOCAB, overrides=GOOD_OVERRIDES):
    (tmp_path / "plbuild.yaml").write_text(config)
    (tmp_path / "vocabulary.yaml").write_text(vocab)
    (tmp_path / "overrides.yaml").write_text(overrides)
    return tmp_path / "plbuild.yaml"


def expect_error(path, fname, line):
    with pytest.raises(ConfigError) as e:
        load_inputs(path)
    assert Path(e.value.path).name == fname, str(e.value)
    assert e.value.line == line, str(e.value)
    assert str(e.value).startswith(f"{e.value.path}:{line}:")
    return e.value


def test_good_inputs_load(tmp_path):
    inp = load_inputs(write(tmp_path))
    assert [m.id for m in inp.config.models] == ["BIOMD0000000427", "BIOMD0000000001"]
    assert inp.config.models[0].version == "5"
    assert inp.config.models[0].default_compartment == "tissue"
    assert inp.config.compartments["blood"].synonyms == ["plasma"]
    assert inp.config.time_unit == "hour" and inp.config.bounds_factor == 100
    assert inp.vocabulary["IL6"].molar_mass == 20900
    o = inp.overrides["BIOMD0000000427"]
    assert o.species == {"s1": "IL6"} and o.bounds == {"k1": (0.1, 10.0)} and o.units["time"] == "minute"


def test_example_library_config_loads():
    inp = load_inputs(EXAMPLE)
    assert inp.config.models[0].id == "BIOMD0000000427"
    assert len(inp.vocabulary) == 21


# --- plbuild.yaml


@pytest.mark.parametrize("config, line, text", [
    ("vocabulary: vocabulary.yaml\ncompartments: {blood: {}}\nmodels: [BIOMD0000000427\n", 4, "syntax"),
    (GOOD_CONFIG + "colour: red\n", 12, "unknown key 'colour'"),
    (GOOD_CONFIG.replace("  - BIOMD0000000001", "  - BIOMD1"), 11, "not a BioModels id"),
    (GOOD_CONFIG.replace("default_compartment: tissue", "default_compartment: cytoplasm"), 10, "cytoplasm"),
    (GOOD_CONFIG.replace("BIOMD0000000001", "BIOMD0000000427"), 11, "listed twice"),
    (GOOD_CONFIG.replace("  tissue: {}\n", "  2tissue: {}\n"), 6, "identifier"),
    (GOOD_CONFIG.replace("version: 5", "version: [5]"), 9, "version"),
    ("vocabulary: vocabulary.yaml\nmodels: []\n", 1, "compartments"),
])
def test_malformed_config(tmp_path, config, line, text):
    err = expect_error(write(tmp_path, config=config), "plbuild.yaml", line)
    assert text in str(err)


# --- vocabulary.yaml


@pytest.mark.parametrize("vocab, line, text", [
    ("IL6:\n  class: molecule\n  unit: pg/mL\n   bad: indent\n", 4, "syntax"),
    (GOOD_VOCAB.replace("class: cell", "class: organ"), 8, "must be one of"),
    (GOOD_VOCAB.replace("  unit: cells/uL\n", ""), 7, "missing required key 'unit'"),
    (GOOD_VOCAB.replace("molar_mass: 20900", "molar_mass: heavy"), 4, "number"),
    (GOOD_VOCAB.replace("IL6:", "IL-6:"), 1, "identifier"),
    (GOOD_VOCAB + "IL6:\n  class: molecule\n  unit: pg/mL\n", 10, "duplicate key"),
    (GOOD_VOCAB.replace("synonyms: [IL-6, interleukin 6]", "synonyms: IL-6"), 6, "list"),
])
def test_malformed_vocabulary(tmp_path, vocab, line, text):
    err = expect_error(write(tmp_path, vocab=vocab), "vocabulary.yaml", line)
    assert text in str(err)


# --- overrides.yaml


@pytest.mark.parametrize("overrides, line, text", [
    ("BIOMD0000000427:\n  units: {time: minute\n", 3, "syntax"),
    (GOOD_OVERRIDES.replace("species: {s1: IL6}", "species: {s1: IL7}"), 3, "not in the vocabulary"),
    (GOOD_OVERRIDES.replace("compartments: {c1: blood}", "compartments: {c1: nucleus}"), 4, "not declared"),
    (GOOD_OVERRIDES.replace("[0.1, 10]", "[10, 0.1]"), 6, "lower must be less"),
    (GOOD_OVERRIDES.replace("[0.1, 10]", "[0.1]"), 6, "[lower, upper]"),
    (GOOD_OVERRIDES.replace("time: minute", "speed: fast"), 2, "unknown key 'speed'"),
    (GOOD_OVERRIDES + "  excluded: [J2]\n", 8, "unknown key 'excluded'"),
    ("BIOMD12:\n  exclude: [J1]\n", 1, "not a BioModels id"),
])
def test_malformed_overrides(tmp_path, overrides, line, text):
    err = expect_error(write(tmp_path, overrides=overrides), "overrides.yaml", line)
    assert text in str(err)


def test_missing_vocabulary_file(tmp_path):
    cfg = write(tmp_path)
    (tmp_path / "vocabulary.yaml").unlink()
    with pytest.raises(ConfigError, match="file not found"):
        load_inputs(cfg)
