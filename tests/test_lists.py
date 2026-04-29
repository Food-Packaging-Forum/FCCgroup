import pandas as pd

from fccgroup.constants import CAS_COLUMN, LIST_NAME_COLUMN
from tests.conftest import _repo_root


def test_lists_columns_match(universe, lists_reference_df):
    lists_reference_df.columns = pd.MultiIndex.from_tuples([("Identifier", col) if col == CAS_COLUMN else ("Lists", col) for col in lists_reference_df.columns])
    missing = [col for col in lists_reference_df if col not in universe.columns]
    assert not missing, f"Some list columns are missing in the output: {missing}"
    