import pandas as pd

from tests.conftest import _repo_root


def test_lists_columns_match(universe, lists_reference_df):
    lists_reference_df.columns = pd.MultiIndex.from_tuples([("Identifier", col) if col == "casId_main" else ("Lists", col) for col in lists_reference_df.columns])
    missing = [col for col in lists_reference_df if col not in universe.columns]
    assert not missing, f"Some list columns are missing in the output: {missing}"


def test_lists_value_comparison(universe, lists_reference_df):
    idxs = lists_reference_df.columns.slice_locs("Alkylphenols_G01")
    lists_reference_df = lists_reference_df[["casId_main"] + list(lists_reference_df.columns[idxs[0]:idxs[1]])]
    lists_reference_df.columns = pd.MultiIndex.from_tuples([("Identifier", col) if col == "casId_main" else ("Lists", col) for col in lists_reference_df.columns])
    
    comparison_df = universe.merge(lists_reference_df, on=[("Identifier", "casId_main")], suffixes=("", "_x"))
    parent_col = "Lists"
    # Compare that values are the same in both datasets
    excel_filename = _repo_root() / "tests" / f"{parent_col}_comparison.xlsx"
    matching_values = True
    with pd.ExcelWriter(excel_filename) as writer:
        for parent, col in comparison_df.columns:
                if parent == f"{parent_col}_x" and (parent_col, col) in comparison_df.columns:
                    if not all(comparison_df[(parent, col)].fillna("") == comparison_df[("Lists", col)].fillna("")):
                        comparison_df[comparison_df[(parent_col, col)].fillna("") != comparison_df[(f'{parent_col}_x', col)].fillna("")][[(parent_col, col), (f'{parent_col}_x', col)]].to_excel(writer, sheet_name=col)
                        matching_values = False
    
    assert matching_values, f"Some list columns do not match, check the file {excel_filename}"