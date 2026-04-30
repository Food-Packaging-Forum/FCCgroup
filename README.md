# FCCgroup

[![PyPI](https://img.shields.io/pypi/v/fccgroup)](https://pypi.org/project/fccgroup/)
[![Python](https://img.shields.io/pypi/pyversions/fccgroup)](https://pypi.org/project/fccgroup/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Issues](https://img.shields.io/github/issues/Food-Packaging-Forum/fccgroup)](https://github.com/Food-Packaging-Forum/fccgroup/issues)

FCCgroup is a Python package for grouping chemicals with three complementary methods:

- Structural pattern matching with SMARTS fingerprints
- Functional list matching against packaged reference lists
- Regex-based grouping from chemical names and formulas

The package is developed under the organization of Food Packaging Forum.
Authored by Albert Anguera Sempere and Helene Wiesinger.

## Features

- Structural classification using SMARTS fingerprints
- Functional list matching from packaged assets
- Regex-based classification from names and formulas
- Automatic CompTox enrichment when selected methods require missing fields
- Flexible method selection through `GroupingConfig(methods=...)`
- Optional SMARTS fingerprint subsetting via `GroupingConfig(smarts_fingerprints=...)`
- Package data bundled under `fccgroup/assets`

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
- `name_columns` and `formula` are optional at configuration time, but `REGEX` grouping may trigger CompTox enrichment when they are missing.
- Input column names can be custom; FCCgroup maps them to canonical internal fields.

## Assets And External Services

- Packaged assets live under `fccgroup/assets`.
- `Mapping.xlsx` and the files in `fccgroup/assets/lists` are required for LISTS workflow.
- CompTox (EPA) is used only when the selected methods require fields that are not already available in the mapped input columns (e.g. SMILES needed for SMARTS but only CAS provided).
- CompTox enrichment requires a valid API key set in the `COMPTOX_API_KEY` environment variable.
- CompTox usage depends on network availability and the EPA CompTox service.

## Output

`group_chemicals(save=True)` returns a pandas DataFrame with a **MultiIndex** on columns. The first level groups results by method; the second level is the column name.

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

FCCgroup currently declares the runtime dependencies described in [requirements.txt](./requirements.txt)

## Citation

If you use FCCgroup in your research, please cite:

```text
@software{fccgroup,
  title={FCCgroup: Chemical Grouping and Classification Package},
  author={Anguera Sempere, Albert and Wiesinger, Helene},
  organization={Food Packaging Forum},
  year={2026},
}
```

## Contributing

Contributions are welcome through pull requests.

## Support

For issues, questions, or suggestions, open an issue at https://github.com/Food-Packaging-Forum/fccgroup/issues.

## License

MIT License. See LICENSE for details.
