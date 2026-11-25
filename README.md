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
| L1 | Fetch and cache, `search` | not started; models are read from a local `sbml:` file |
| L2 | Reaction extraction and exclusion rules | done |
| L3 | Vocabulary, overrides, species and compartment mapping | done, except the vocabulary output file |
| L4 | Units, compartment factor, stoichiometry scaling | compartment factor done; only identity unit conversion |
| L5 | Identity keys, merging, writer, report, manifest | done, except merging parameter values across models |
| L6 | Validation, reconstruction test, first real library | first library from BIOMD0000000427; EnsGen parser check waits for EnsGen |

## What a build writes

| File | Content |
| --- | --- |
| `curated.ant` | The library EnsGen loads: `library.size` entries drawn (seeded) from the curated processes plus the templates. Every entry has a comment block explaining it, so people can read and edit the file. |
| `pool.ant` | Every candidate the library is drawn from. EnsGen does not load it. |
| `report.json`, `report.html` | Reactions kept and excluded (with reason), unresolved names, and every curated process left out of the library with the drawn templates that can stand in for it, or none (a "gap"). |
| `manifest.json` | Tool version, input and SBML hashes, output hashes, build date. |

## Templates

The generic processes (synthesis, degradation, conversion, Michaelis-Menten, transport, 2X -> Y, X + Y -> Z, and the
linear, activation and inhibition multipliers) are defined in `src/plbuild/templates.py` and written into the library
like any other entry. EnsGen has none built in.

Template parameters start at 1 and are bounded by `bounds.templates` in `plbuild.yaml` (default `[1.0e-6, 1.0e7]`);
Hill coefficients are always bounded 1 to 4. Curated parameters are bounded by `bounds.factor` around their source
value. The EnsGen spec's example range, 1e-3 to 1e3, is too narrow for models whose species are molecule counts.

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
species, and check that the ensemble reproduces them using only processes from this library. The library holds 25 of
the 38 candidates (22 curated + 16 templates), so EnsGen cannot simply rebuild the source model.
