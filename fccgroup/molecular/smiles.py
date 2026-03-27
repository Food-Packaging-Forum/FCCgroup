"""
SMILES string validation and canonical conversion utilities.

This module handles:
- SMILES validation and parsing
- Canonical SMILES generation
- SMILES sanitization
"""

from typing import Optional

import pandas as pd
from rdkit import Chem


def smiles_check(smiles: str) -> bool:
    """
    Checks whether a given SMILES string is valid.

    Parameters:
        smiles (str): The SMILES string to validate.

    Returns:
        bool: True if the SMILES string is valid and can be parsed into a molecule, False otherwise.

    Notes:
        - Returns False if the input is an empty string or NaN.
        - Uses RDKit's Chem.MolFromSmiles for validation.
    """
    if smiles == "" or pd.isna(smiles):
        return False
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception:
        return False


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """
    Canonicalize a SMILES string using RDKit.

    Converts a SMILES string to its canonical form, which provides a consistent
    representation for the same molecule. This is essential for comparing SMILES strings
    and ensuring consistent lookups in preprocessed databases.

    Args:
        smiles: SMILES string to canonicalize

    Returns:
        Canonicalized SMILES string, or None if the input is invalid
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None
