import inspect
import ast
import re
import textwrap
import pandas as pd

from tests.conftest import _repo_root

def test_pattern_names_in_fingerprints(pattern_df, fingerprints):
    columns_missmatch = []
    for column in pattern_df.Group:
        if column not in fingerprints.keys():
            if column in ['Molecular composition', 'Not groupable', 'SMILES check']:
                continue
            columns_missmatch.append(column)
    assert not columns_missmatch, f"Some columns in mapping are missing in fingerprints {columns_missmatch}"


def test_fingerprints_in_pattern_names(pattern_df, fingerprints):
    columns_missmatch = []
    for column in fingerprints.keys():
        if column not in pattern_df.Group.values:
            columns_missmatch.append(column)
    assert not columns_missmatch, f"Some columns in fingerprints are missing in mapping {columns_missmatch}"


def test_pattern_strings_in_fingerprints(pattern_df, fingerprints):
    unmatching_patterns = []
    filename = _repo_root() / "tests" / "Unmatching_patterns.xlsx"
    fingerprints_stringified = {key: (val if isinstance(val, str) else extract_pattern(val)) for key, val in fingerprints.items()}
    for _, row in pattern_df.iterrows():
        name = row["Group"]
        pattern = row["SMARTS"]
        if pd.isna(pattern) or pattern.strip() in ["", "Function"]:
            continue
        pattern_found = False
        for pattern_str in fingerprints_stringified.values():
            if pattern in str(pattern_str):
                pattern_found = True
                break
        if not pattern_found:
            actual_smarts = fingerprints_stringified.get(name, "")
            unmatching_patterns.append({"SMARTS": pattern, "Actual SMARTS": actual_smarts, "Check": ""})
    if unmatching_patterns:
        pd.DataFrame(unmatching_patterns).to_excel(filename, index=False)
    assert not unmatching_patterns, f"Pattern not found in fingerprints {unmatching_patterns}"

def extract_pattern(fn):
    """Extract first string literal from lambda body, including nested calls."""
    if not (callable(fn) and fn.__name__ == "<lambda>"):
        return None
    try:
        src = inspect.getsource(fn)
    except (OSError, TypeError):
        return None

    src = textwrap.dedent(src).strip()
    m = re.search(r"(lambda[^\n]*)", src)
    if not m:
        return None
    lambda_src = m.group(1)
    return lambda_src
