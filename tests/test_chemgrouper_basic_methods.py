"""Basic one-chemical tests for GroupingMethod combinations."""

from pathlib import Path

import pandas as pd
import pytest

from fccgroup import ChemicalGrouper, GroupingConfig, GroupingMethod, ColumnMapping
from fccgroup.constants import *
from fccgroup.constants import MULTIINDEX_IDENTIFIER_LABEL, MULTIINDEX_STRUCTURAL_LABEL


def _one_chemical_df() -> pd.DataFrame:
    """Single-chemical input with custom column names."""
    return pd.DataFrame(
        {
            "Structure": ["CC"],
            "CASRN": ["74-84-0"],
            "Name": ["ethane"],
            "IUPAC": ["ethane"],
            "Formula": ["C2H6"],
        }
    )


def _run_grouping(methods: list[GroupingMethod]) -> pd.DataFrame:
    """Run grouping for one chemical with selected methods."""
    config = GroupingConfig(
        methods=methods,
        column_mapping=ColumnMapping(
            cas="CASRN",
            smiles="Structure",
            name_columns=["Name", "IUPAC"],
            formula="Formula",
        ),
    )

    fitted_df = _one_chemical_df()
    grouper = ChemicalGrouper(df=fitted_df, grouping_config=config)
    return grouper.group_chemicals(save=False)


def _has_lists_and_regex_assets() -> bool:
    """Check whether optional assets for LISTS/REGEX are available."""

    mapping_file = ASSETS_DIR / MAPPING_FILE_NAME
    lists_dir = ASSETS_DIR / LISTS_DIR.name
    return mapping_file.exists() and lists_dir.exists()


def test_grouping_smarts_only_one_chemical() -> None:
    """SMARTS-only grouping should work for one chemical."""
    result = _run_grouping([GroupingMethod.SMARTS])
    
    print(result)
        
    assert len(result) == 1
    assert (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN) in result.columns
    assert result.loc[0, (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN)] == "CC"
    assert (MULTIINDEX_STRUCTURAL_LABEL, OUTPUT_COLUMN) in result.columns
    assert result.attrs.get("grouping_methods_applied") is None


def test_grouping_filters_out_unmapped_input_columns() -> None:
    """Only mapped input columns should be carried into the grouping workflow."""
    fitted_df = _one_chemical_df().copy()
    fitted_df["ExternalMetadata"] = ["should_not_be_used"]

    config = GroupingConfig(
        methods=[GroupingMethod.SMARTS],
        column_mapping=ColumnMapping(
            cas="CASRN",
            smiles="Structure",
            name_columns=["Name", "IUPAC"],
            formula="Formula",
        ),
    )

    result = ChemicalGrouper(df=fitted_df, grouping_config=config).group_chemicals(save=False)

    assert "ExternalMetadata" not in result.columns


@pytest.mark.skipif(not _has_lists_and_regex_assets(), reason="Mapping.xlsx/lists assets are not available")
def test_grouping_lists_only_one_chemical() -> None:
    """LISTS-only grouping should work for one chemical."""
    result = _run_grouping([GroupingMethod.LISTS])

    print(result)
    
    assert len(result) == 1
    assert (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN) in result.columns
    assert result.loc[0, (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN)] == "CC"
    assert result.attrs.get("grouping_methods_applied") is None


@pytest.mark.skipif(not _has_lists_and_regex_assets(), reason="Mapping.xlsx/lists assets are not available")
def test_grouping_regex_only_one_chemical() -> None:
    """REGEX-only grouping should work for one chemical."""
    result = _run_grouping([GroupingMethod.REGEX])

    print(result)
    
    assert len(result) == 1
    assert (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN) in result.columns
    assert result.loc[0, (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN)] == "CC"
    assert result.attrs.get("grouping_methods_applied") is None


@pytest.mark.skipif(not _has_lists_and_regex_assets(), reason="Mapping.xlsx/lists assets are not available")
def test_grouping_smarts_and_regex_one_chemical() -> None:
    """SMARTS + REGEX grouping should work for one chemical."""
    result = _run_grouping([GroupingMethod.SMARTS, GroupingMethod.REGEX])

    print(result)
    
    assert len(result) == 1
    assert (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN) in result.columns
    assert result.loc[0, (MULTIINDEX_IDENTIFIER_LABEL, SMILES_COLUMN)] == "CC"
    assert (MULTIINDEX_STRUCTURAL_LABEL, OUTPUT_COLUMN) in result.columns
    assert result.attrs.get("grouping_methods_applied") is None
