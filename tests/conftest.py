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

def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent

@pytest.fixture
def pattern_df() -> pd.DataFrame:
    mapping_path = ASSETS_DIR / MAPPING_FILE_NAME
    if not mapping_path.exists():
        pytest.skip("assets/Mapping.xlsx not found")

    raw = pd.read_excel(mapping_path, sheet_name=MAPPING_SHEET_SMARTS)
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
    mapping_path = ASSETS_DIR / MAPPING_FILE_NAME
    if not mapping_path.exists():
        pytest.skip("assets/Mapping.xlsx not found")

    df = pd.read_excel(mapping_path, sheet_name=MAPPING_SHEET_KEYWORDS)
    if REGEX_COLUMN_NAME_COLUMN not in df.columns:
        pytest.skip("B - Keywords sheet is missing 'column_name'")
    df = df[df.use == True].reset_index(drop=True)
    return df

@pytest.fixture
def lists_df() -> pd.DataFrame:
    mapping_path = ASSETS_DIR / MAPPING_FILE_NAME
    if not mapping_path.exists():
        pytest.skip("assets/Mapping.xlsx not found")

    df = pd.read_excel(mapping_path, sheet_name=MAPPING_SHEET_LISTS)
    if GROUP_KEY_COLUMN not in df.columns:
        pytest.skip(f"B - Lists sheet is missing '{GROUP_KEY_COLUMN}'")
    if LIST_ID_COLUMN not in df.columns:
        pytest.skip(f"B - Lists sheet is missing '{LIST_ID_COLUMN}'")
    df = df[df.use == True].reset_index(drop=True)
    df[LIST_NAME_COLUMN] = df[GROUP_KEY_COLUMN] + "_" + df[LIST_ID_COLUMN]
    return df


@pytest.fixture
def regex_combination_dictionary() -> dict:
    return regex_combo_dict

@pytest.fixture
def universe() -> pd.DataFrame:
    """Build a MultiIndex universe DataFrame from repository assets for comparison tests."""
    universe_path = ASSETS_DIR / FCC_UNIVERSE_FILE
    if not universe_path.exists():
        pytest.skip("assets/grouped_chemicals.xlsx not found")

    df = pd.read_excel(universe_path, header=[0,1])
    return df


@pytest.fixture
def lists_reference_df(universe: pd.DataFrame) -> pd.DataFrame:
    list_cols = [col for col in universe.columns if col[0] == "Lists"]
    if not list_cols:
        pytest.skip("No list columns found in universe fixture")

    out = universe[[('Identifier', CAS_COLUMN)] + [("Lists", col[1]) for col in list_cols]].copy()
    out.columns = [CAS_COLUMN] + [col[1] for col in list_cols]
    return out


@pytest.fixture
def regex_reference_df(universe: pd.DataFrame) -> pd.DataFrame:
    regex_cols = [col for col in universe.columns if col[0] == "Regex"]
    if not regex_cols:
        pytest.skip("No regex columns found in universe fixture")

    out = universe[[('Identifier', CAS_COLUMN)] + [("Regex", col[1]) for col in regex_cols]].copy()
    out.columns = [CAS_COLUMN] + [col[1] for col in regex_cols]
    return out


