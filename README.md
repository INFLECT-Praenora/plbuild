# plbuild

`plbuild` builds the process library that [EnsGen](https://github.com/INFLECT-Praenora/EnsembleEvolutionary_Spec), the
evolutionary ensemble generator, loads. It reads a curated list of BioModels models and turns each reaction into a
process: species mapped to a shared vocabulary, compartments mapped to the configured set, rate laws converted to the
library's units, duplicates merged, provenance recorded. The output is `curated.ant`, a restricted Antimony file that
people can also edit by hand.

The specification is `docs/plbuild-spec.md` in
[INFLECT-Praenora/EnsembleEvolutionary_Spec](https://github.com/INFLECT-Praenora/EnsembleEvolutionary_Spec). The spec
changes as we learn what works.

## Status

| Milestone | Content | State |
| --- | --- | --- |
| L0 | Package skeleton, config loading, CLI stubs, CI | done |
| L1 | Fetch and cache, `search` | |
| L2 | Reaction extraction and exclusion rules | |
| L3 | Vocabulary, overrides, species and compartment mapping | |
| L4 | Units, compartment factor, stoichiometry scaling | |
| L5 | Identity keys, merging, writer, report, manifest | |
| L6 | Validation, reconstruction test, first real library | |

## Install and test

```
python -m pip install -e ".[test]"
python -m pytest
```

## Use

```
plbuild --help
plbuild build --config examples/bianconi2012/plbuild.yaml
```

All commands take `--config plbuild.yaml`, which names the vocabulary and overrides files. Every error in these files is
reported as `file:line: message`.

## First library

`examples/bianconi2012/` builds a library from BIOMD0000000427 (Bianconi 2012, EGFR and IGF1R pathway in lung cancer).
That model is also EnsGen's first ground-truth test: simulate it, give EnsGen the noise-free trajectories of all 21
species, and check that the ensemble reproduces them using only processes from this library.
