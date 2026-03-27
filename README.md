# FCCgroup

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
- Automatic CIRpy enrichment when selected methods require missing fields
- Flexible method selection through `GroupingConfig(methods=...)`
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
results = grouper.group_chemicals()

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

## Input Requirements

- `ChemicalGrouper` must be initialized with a non-empty pandas DataFrame.
- `ColumnMapping` must provide at least one of `cas` or `smiles`.
- `name_columns` and `formula` are optional at configuration time, but `REGEX` grouping may trigger CIRpy enrichment when they are missing.
- Input column names can be custom; FCCgroup maps them to canonical internal fields.

## Assets And External Services

- Packaged assets live under `fccgroup/assets`.
- `Mapping.xlsx` and the files in `fccgroup/assets/lists` are required for LISTS workflow.
- CIRpy is used only when the selected methods require fields that are not already available in the mapped input columns.
- CIRpy usage depends on network availability and the external resolver service.

## Output

`group_chemicals()` returns a pandas DataFrame containing the normalized internal identifier columns plus the columns produced by the selected methods.

Typical outputs include:

- `SMILES` and/or `casId` internal identifier columns
- `Chemical groups` and SMARTS fingerprint columns when SMARTS is selected
- Functional list columns when LISTS is selected
- Regex-derived columns when REGEX is selected

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
