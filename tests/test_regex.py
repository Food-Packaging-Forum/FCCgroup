import numpy as np
import pandas as pd

from fccgroup.constants import *
from tests.conftest import _repo_root

def test_regex_columns_match(regex_df, regex_reference_df, regex_combination_dictionary):
    regex_columns = list(regex_df[REGEX_COLUMN_NAME_COLUMN].unique()) + ["Organic_C", "Inorganic_noC", "Metal_Metalloid"] + list(regex_combination_dictionary.keys())
    missing = [col for col in regex_columns if col not in regex_reference_df.columns]
    assert not missing, f"Some regex columns are missing in the output: {missing}"


