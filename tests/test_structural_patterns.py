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
            if column in ['Molecular composition', 'Not groupable']:
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
        if pattern not in fingerprints_stringified.values():
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
    try:
        tree = ast.parse(lambda_src)
    except SyntaxError:
        return None

    lam = next((n for n in ast.walk(tree) if isinstance(n, ast.Lambda)), None)
    if lam is None:
        return None

    # Recursively walk the lambda body to find string constants in calls
    for node in ast.walk(lam.body):
        if isinstance(node, ast.Call):
            for arg in node.args:
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    return arg.value
    return None

if __name__ == "__main__":
    import os,sys
    CODE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if CODE_PATH not in sys.path:
        sys.path.insert(0, CODE_PATH)
    # Import helper symbols from package modules in this repository
    from fccgroup.patterns.library import apply_pattern, cx_smarts_query
    from fccgroup.patterns.library import is_aromatic
    from fccgroup.data.constants import elements
    from rdkit import Chem
    filters = {
        "Alkane": lambda x: apply_pattern(Chem.AddHs(x), "CC"),
        "Carbonyl": lambda x: apply_pattern(Chem.AddHs(x), "C=O"),
        "Organo O (strict)": lambda x, row: {"C","O"}.issubset(row["Molecular Composition"].keys()),
        "Contains posttransition metal": lambda x, row: len(set(elements["posttransition_metals"]).intersection(set(row["Molecular Composition"]))),
        "Aromatics": is_aromatic,
        "Chlorophenols": lambda x, row: row["Contains C~Cl"] > 0 and cx_smarts_query(Chem.AddHs(x), 'c([#1,Cl])1c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1[#8][#6,#1]')
    }

    for name, fn in filters.items():
        pat = extract_pattern(fn)
        print(name, "→", pat or "None")
