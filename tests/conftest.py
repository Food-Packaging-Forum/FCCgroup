"""Pytest fixtures and configuration."""

from pathlib import Path
import re

import pandas as pd
import pytest

from fccgroup import ChemicalGrouper, GroupingConfig, GroupingMethod, ColumnMapping
from fccgroup.constants import *
from fccgroup.patterns.library import fingerprints as pattern_fingerprints
from fccgroup.data.constants import regex_combination_dictionary as regex_combo_dict


@pytest.fixture
def grouper_smarts_only():
    """Create a ChemicalGrouper configured with the SMARTS method only."""
    fitted_df = pd.DataFrame({
        SMILES_COLUMN: ['CC'],
        CAS_COLUMN: ['74-84-0'],
        'Name': ['ethane'],
        'IUPAC': ['ethane'],
        'Formula': ['C2H6'],
    })
    config = GroupingConfig(
        methods=[GroupingMethod.SMARTS],
        column_mapping=ColumnMapping(
            cas=None,
            smiles=SMILES_COLUMN,
            name_columns=["Name", "IUPAC"],
            formula="Formula",
        ),
    )
    return ChemicalGrouper(df=fitted_df, grouping_config=config)


@pytest.fixture
def formaldehyde_df():
    """Create a test DataFrame for formaldehyde (CAS 50-00-0)."""
    return pd.DataFrame({
        SMILES_COLUMN: ['C=O'],
        'canonical_SMILES': ['C=O'],
        CAS_COLUMN: ['50-00-0'],
    })


@pytest.fixture
def ethane_df():
    """Create a test DataFrame for ethane (CAS 74-84-0)."""
    return pd.DataFrame({
        SMILES_COLUMN: ['CC'],
        'canonical_SMILES': ['CC'],
        CAS_COLUMN: ['74-84-0'],
    })

@pytest.fixture
def load_lookup_df():
    """Load the lookup DataFrame from the CSV file."""
    lookup_path = _assets_root() / SMILES_LOOKUP_FILE.name
    if not lookup_path.exists():
        pytest.skip("smiles_lookup.tsv not found in packaged assets")
    return pd.read_csv(lookup_path, sep='\t')


@pytest.fixture
def load_lookup_dict(load_lookup_df):
    """Build a canonical SMILES → row-index lookup dict from the lookup DataFrame."""
    from rdkit import Chem

    def canonical(s):
        try:
            mol = Chem.MolFromSmiles(s)
            return Chem.MolToSmiles(mol, canonical=True) if mol is not None else None
        except Exception:
            return None

    smiles_col = SMILES_COLUMN if SMILES_COLUMN in load_lookup_df.columns else load_lookup_df.columns[0]
    lookup = {}
    for idx, row in load_lookup_df.iterrows():
        key = canonical(str(row[smiles_col]))
        if key is not None:
            lookup[key] = idx
    return lookup


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _assets_root() -> Path:
    package_assets = _repo_root() / "fccgroup" / "assets"
    if package_assets.exists():
        return package_assets
    return _repo_root() / "assets"


@pytest.fixture
def pattern_df() -> pd.DataFrame:
    mapping_path = _assets_root() / MAPPING_FILE_NAME
    if not mapping_path.exists():
        pytest.skip("assets/Mapping.xlsx not found")

    raw = pd.read_excel(mapping_path, sheet_name=MAPPING_SHEET_SMARTS, header=1)
    if "Column_name" not in raw.columns or "SMARTS_pattern" not in raw.columns:
        pytest.skip("B - SMARTS sheet is missing expected columns")

    df = raw.rename(columns={"Column_name": "Group", "SMARTS_pattern": "SMARTS"})
    # Strip whitespace that may be present in the Excel source
    df["Group"] = df["Group"].str.strip()
    df["SMARTS"] = df["SMARTS"].apply(lambda v: v.strip() if isinstance(v, str) else v)
    # Only include active patterns (use=True) to avoid false positives from intentionally
    # disabled entries that have no corresponding fingerprint implementation
    if USE_COLUMN in df.columns:
        df = df[df[USE_COLUMN] == True].reset_index(drop=True)
    return df


@pytest.fixture
def fingerprints() -> dict:
    return pattern_fingerprints


@pytest.fixture
def regex_df() -> pd.DataFrame:
    mapping_path = _assets_root() / MAPPING_FILE_NAME
    if not mapping_path.exists():
        pytest.skip("assets/Mapping.xlsx not found")

    df = pd.read_excel(mapping_path, sheet_name=MAPPING_SHEET_KEYWORDS)
    if REGEX_COLUMN_NAME_COLUMN not in df.columns:
        pytest.skip("B - Keywords sheet is missing 'column_name'")
    return df


@pytest.fixture
def regex_combination_dictionary() -> dict:
    return regex_combo_dict


@pytest.fixture
def regex_columns(regex_df: pd.DataFrame, regex_combination_dictionary: dict) -> list[str]:
    base = [col for col in regex_df.get(REGEX_COLUMN_NAME_COLUMN, pd.Series([], dtype=str)).dropna().unique()]
    return base + ["Organic_C", "Inorganic_noC", "Metal_Metalloid"] + list(regex_combination_dictionary.keys())


@pytest.fixture
def universe(regex_columns: list[str]) -> pd.DataFrame:
    """Build a MultiIndex universe DataFrame from repository assets for comparison tests."""
    universe_path = _assets_root() / FCC_UNIVERSE_FILE.name.replace("_in", "_all")
    if not universe_path.exists():
        pytest.skip("assets/FCCuniverse_grouping_all.xlsx not found")

    df = pd.read_excel(universe_path, sheet_name="FCCs")
    if "casId_main" not in df.columns:
        pytest.skip("FCCuniverse_grouping_all.xlsx::FCCs missing 'casId_main'")

    identifier_cols = {
        "casId_main",
        SMILES_COLUMN,
        COMMON_NAME_COLUMN,
        IUPAC_NAME_COLUMN,
        FORMULA_COLUMN,
    }
    regex_cols = set(regex_columns)

    multi_cols = []
    for col in df.columns:
        if col in identifier_cols:
            multi_cols.append(("Identifier", col))
        elif col in regex_cols:
            multi_cols.append(("Regex", col))
        elif re.search(r"_G\d{2}$", str(col)):
            multi_cols.append(("Lists", col))
        else:
            multi_cols.append(("Other", str(col)))

    df = df.copy()
    df.columns = pd.MultiIndex.from_tuples(multi_cols)
    return df


@pytest.fixture
def lists_reference_df(universe: pd.DataFrame) -> pd.DataFrame:
    list_cols = [col for col in universe.columns if col[0] == "Lists"]
    if not list_cols:
        pytest.skip("No list columns found in universe fixture")

    out = universe[[('Identifier', 'casId_main'), *list_cols]].copy()
    out.columns = ["casId_main"] + [col[1] for col in list_cols]
    return out


@pytest.fixture
def regex_reference_df(universe: pd.DataFrame) -> pd.DataFrame:
    regex_cols = [col for col in universe.columns if col[0] == "Regex"]
    if not regex_cols:
        pytest.skip("No regex columns found in universe fixture")

    out = universe[[('Identifier', 'casId_main'), *regex_cols]].copy()
    out.columns = ["casId_main"] + [col[1] for col in regex_cols]
    return out


