# FCCgroup

[![PyPI](https://img.shields.io/pypi/v/fccgroup)](https://pypi.org/project/fccgroup/)
[![Python](https://img.shields.io/pypi/pyversions/fccgroup)](https://pypi.org/project/fccgroup/)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](LICENSE)
[![Issues](https://img.shields.io/github/issues/Food-Packaging-Forum/fccgroup)](https://github.com/Food-Packaging-Forum/fccgroup/issues)
[![Paper](https://img.shields.io/badge/DOI-10.1021%2Facs.est.5c15186-blue)](https://doi.org/10.1021/acs.est.5c15186)
[![DeepWiki](https://img.shields.io/badge/DeepWiki-FCCgroup-1f6feb)](https://deepwiki.com/Food-Packaging-Forum/FCCgroup)

FCCgroup is a modular Python package for transparent and reproducible chemical grouping. It
combines structure-based rules, curated chemical lists, and name- and formula-based patterns to
help prioritize and analyze food contact chemical datasets.

The three grouping methods can be selected independently or combined:

- Structural pattern matching with SMARTS fingerprints
- Functional list matching against packaged reference lists
- Name- and formula-based pattern matching with regular expressions

The package is developed under the organization of Food Packaging Forum.
Authored by Albert Anguera Sempere and Helene Wiesinger.

## Features

- Structural classification using SMARTS fingerprints
- Functional list matching from packaged assets
- Regex-based classification from names and formulas
- Automatic CompTox enrichment when selected methods require missing fields
- Flexible method selection through `GroupingConfig(methods=...)`
- Optional SMARTS fingerprint subsetting via `GroupingConfig(smarts_fingerprints=...)`
- Lazy loading of method-specific rules and data assets
- Optional CompTox enrichment for missing identifiers, structures, names, or formulas
- Package data bundled under `fccgroup/assets`

## How It Works

FCCgroup takes a pandas DataFrame and a `GroupingConfig`. `ColumnMapping` connects the input
columns to the chemical fields used by the selected methods. The package then:

1. Validates the configured columns and selected methods.
2. Resolves missing fields through the EPA CompTox service when enrichment is required.
3. Applies the selected SMARTS, list, and regex rules.
4. Returns the input identifiers and grouping results in a pandas DataFrame with two-level
  MultiIndex columns.

Resources are loaded only for the methods selected in the configuration. This makes it possible
to run a focused structural or list-based workflow without loading unrelated assets.

## Installation

Install from PyPI:

```bash
pip install fccgroup
```

Install from source:

```bash
git clone https://github.com/Food-Packaging-Forum/fccgroup.git
cd fccgroup
pip install -e .
```

Install development dependencies:

```bash
pip install -e .[dev]
```

## Quick Start

```python
import pandas as pd

from fccgroup import ChemicalGrouper, ColumnMapping, GroupingConfig, GroupingMethod

df = pd.DataFrame(
  {
    "CASRN": ["74-84-0"],
    "Structure": ["CC"],
    "Name": ["ethane"],
    "IUPAC": ["ethane"],
    "Formula": ["C2H6"],
  }
)

config = GroupingConfig(
  methods=[GroupingMethod.SMARTS, GroupingMethod.REGEX],
  column_mapping=ColumnMapping(
    cas="CASRN",
    smiles="Structure",
    name_columns=["Name", "IUPAC"],
    formula="Formula",
  ),
)

grouper = ChemicalGrouper(df=df, grouping_config=config)
results = grouper.group_chemicals(save=False)

# Columns are a MultiIndex: (group_label, column_name)
print(results.columns.tolist())
print(results.head())
```

## Selecting Grouping Methods

FCCgroup does not expose a `GroupingMode` enum. Method selection is configured with `GroupingMethod` values:

- `GroupingMethod.SMARTS`: structural pattern matching
- `GroupingMethod.LISTS`: functional list matching
- `GroupingMethod.REGEX`: regex-based grouping from names and formulas

Common configurations:

```python
GroupingConfig(methods=[GroupingMethod.SMARTS], column_mapping=...)
GroupingConfig(methods=[GroupingMethod.SMARTS, GroupingMethod.LISTS], column_mapping=...)
GroupingConfig(
  methods=[GroupingMethod.SMARTS, GroupingMethod.LISTS, GroupingMethod.REGEX],
  column_mapping=...,
)
```

### Filtering SMARTS fingerprints

To apply only a subset of the ~400 bundled SMARTS patterns, pass their names to `smarts_fingerprints`:

```python
GroupingConfig(
  methods=[GroupingMethod.SMARTS],
  column_mapping=...,
  smarts_fingerprints={"Alkanes", "PAH derivatives hydrocarbon"},
)
```

When `smarts_fingerprints` is `None` (default), all available patterns are applied.

### Custom assets path

By default `ChemicalGrouper` loads assets from the package installation directory. To point it at a different directory:

```python
ChemicalGrouper(df=df, grouping_config=config, assets_path="/path/to/custom/assets")
```

## Input Requirements

- `ChemicalGrouper` must be initialized with a non-empty pandas DataFrame.
- `ColumnMapping` must provide at least one of `cas` or `smiles` (the other may be `None`).
- `name_columns` and `formula` are optional at configuration time, but `REGEX` grouping may trigger
  CompTox enrichment when they are missing.
- Input column names can be custom; FCCgroup maps them to canonical internal fields.

CAS identifiers are the preferred input for list matching. SMILES are required for direct SMARTS
matching. When only one identifier type is supplied, CompTox can be used to resolve the fields
needed by the selected methods.

## Assets And External Services

- Packaged assets live under `fccgroup/assets`.
- `Mapping.xlsx` and the files in `fccgroup/assets/lists` are required for LISTS workflow.
- CompTox (EPA) is used only when the selected methods require fields that are not already
  available in the mapped input columns (e.g. SMILES needed for SMARTS but only CAS provided).
- CompTox enrichment requires a valid API key set in the `COMPTOX_API_KEY` environment variable.
- CompTox usage depends on network availability and the EPA CompTox service.

## Output

`group_chemicals(save=True)` returns a pandas DataFrame with a **MultiIndex** on columns. The first
level groups results by method; the second level is the column name. With the default `save=True`,
the same results are also written to `Grouping.xlsx` in the current working directory. Pass
`save=False` to keep the workflow in memory.

| Top-level label | Contents |
|---|---|
| `Identifier` | Internal identifier columns (`casId`, `SMILES`) |
| `Structural patterns` | `Chemical groups` and per-fingerprint columns (SMARTS method) |
| `Lists` | Per-list membership columns (LISTS method) |
| `Regex` | Pattern group columns (REGEX method) |

Example column access:

```python
# Access the SMILES identifier column
results[("Identifier", "SMILES")]

# Access the Chemical groups column
results[("Structural patterns", "Chemical groups")]
```

When `save=True` (default), results are also written to an Excel file in the current working directory.

## Runtime Dependencies

FCCgroup declares its runtime dependencies in [pyproject.toml](pyproject.toml), including pandas,
NumPy, RDKit, requests, joblib, openpyxl, and tqdm.

## Citation

If you use FCCgroup in your research, please cite the associated paper:

> Wiesinger, H., Parkinson, L. V., Geueke, B., Anguera Sempere, A., Boucher, J., Cabane, E.,
> Scheringer, M., Muncke, M. (2026). Prioritizing and Grouping Food Contact Chemicals – From
> Chaos to Clarity. *Environmental Science & Technology*. DOI:
> [10.1021/acs.est.5c15186](https://doi.org/10.1021/acs.est.5c15186).

The machine-readable citation is available in [citation.cff](citation.cff). The software can also
be cited as:

```text
@software{fccgroup,
  title={FCCgroup: A Modular Chemical Grouping Tool},
  author={Anguera Sempere, Albert and Wiesinger, Helene},
  organization={Food Packaging Forum},
  year={2026},
}
```

For project documentation and an interactive overview, see the [FCCgroup DeepWiki](https://deepwiki.com/Food-Packaging-Forum/FCCgroup).

## Contributing

Contributions are welcome through pull requests.

## Support

For issues, questions, or suggestions, open an issue at https://github.com/Food-Packaging-Forum/fccgroup/issues.

## License

Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0).
See LICENSE for details.

## Disclaimer

This software is provided "as is", without warranties of any kind, express or implied.
To the maximum extent permitted by applicable law, Food Packaging Forum and contributors
shall not be liable for any direct, indirect, incidental, special, exemplary, or
consequential damages arising from the use or misuse of this software.
