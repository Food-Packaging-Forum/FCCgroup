import numpy as np
import pandas as pd

from fccgroup.constants import *
from tests.conftest import _repo_root

def test_regex_columns_match(regex_df, regex_reference_df, regex_combination_dictionary):
    regex_columns = list(regex_df[REGEX_COLUMN_NAME_COLUMN].unique()) + ["Organic_C", "Inorganic_noC", "Metal_Metalloid"] + list(regex_combination_dictionary.keys())
    missing = [col for col in regex_columns if col not in regex_reference_df.columns]
    assert not missing, f"Some regex columns are missing in the output: {missing}"


def test_regex_value_comparison(universe, regex_df, regex_reference_df, regex_columns):
    regex_columns_test = [col for col in regex_columns if col in regex_reference_df.columns]
    for col in regex_columns_test:
        regex_reference_df[col] = regex_reference_df[col].astype('boolean').apply(lambda x: np.nan if type(x) is not bool else x) # Process column to avoid having empty strings
    regex_reference_df = regex_reference_df[["casId_main"] + regex_columns_test]
    regex_reference_df.columns = pd.MultiIndex.from_tuples([("Identifier", col) if col == "casId_main" else ("Regex", col) for col in regex_reference_df.columns])
    comparison_df = universe.merge(regex_reference_df, on=[("Identifier", "casId_main")], suffixes=("", "_x"))

    # Compare that values are the same in both datasets
    excel_filename = _repo_root() / "tests" / "Regex_comparison.xlsx"
    matching_values = True
    with pd.ExcelWriter(excel_filename) as writer:
        for parent, col in comparison_df.columns:
                if parent == "Regex_x" and ("Regex", col) in comparison_df.columns:
                    if not all(comparison_df[(parent, col)].fillna("") == comparison_df[("Regex", col)].fillna("")):
                        val = 'name'
                        try:
                            val = regex_df[regex_df[REGEX_COLUMN_NAME_COLUMN] == col][REGEX_KEYWORD_LOCATION_COLUMN].iloc[0]
                        except:
                            pass
                        if val == COMBINED_NAME_COLUMN:
                            val = ('Regex', val)
                        elif val == SMILES_COLUMN:
                            val = ('Identifier', SMILES_COLUMN)
                        else:
                            val = ('Identifier', val)
                        comparison_df[comparison_df[('Regex', col)].fillna("") != comparison_df[('Regex_x', col)].fillna("")][[val, ('Regex', col), ('Regex_x', col)]].to_excel(writer, sheet_name=col)
                        matching_values = False
    
    assert matching_values, f"Some regex columns do not match, check the file {excel_filename}"
