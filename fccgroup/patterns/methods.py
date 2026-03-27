from typing import Union, Callable, Dict

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolEnumerator

from ..constants import *
from ..molecular.composition import align_bundle_coords, molecule_composition
from .library import apply_pattern

def generate_fingerprints(
    index: str,
    smiles: str,
    is_valid_smiles: bool = True,
    fpps: Dict[str, Union[str, Callable]] = None
) -> Dict[str, Union[bool, str, float]]:
    """
    Generate a dictionary of fingerprints for a given SMILES string.
    
    Creates a copy of the fingerprints dictionary and applies all SMARTS patterns
    if the SMILES string is valid, otherwise returns NaN values.
    
    Args:
        index: Identifier for the molecule (in case SMILES is not unique)
        smiles: SMILES string representing the chemical structure
        is_valid_smiles: Whether the provided SMILES string is valid. Defaults to True
        fpps: Dictionary of patterns to apply. Defaults to module-level fingerprints
    
    Returns:
        Dictionary containing boolean fingerprints for each pattern, plus:
            - "ID": The provided index value
            - "Molecular Composition": For valid SMILES only
            - "Bond Composition": For valid SMILES only
    """
    if fpps is None:
        fpps = fingerprints

    FINGERPRINTS_DEFAULT: Dict[str, bool] = {k: False for k in fpps.keys()}
    FINGERPRINTS_INVALID_DEFAULT: Dict[str, float] = {k: np.nan for k in fpps.keys()}

    if is_valid_smiles:
        compound_fingerprint: Dict[str, Union[bool, str, float]] = FINGERPRINTS_DEFAULT.copy()
        compound_fingerprint[FINGERPRINT_ID_COLUMN] = index
        apply_all_patterns(smiles, fpps, compound_fingerprint)
    else:
        compound_fingerprint = FINGERPRINTS_INVALID_DEFAULT.copy()
        compound_fingerprint[FINGERPRINT_ID_COLUMN] = index
    
    return compound_fingerprint


def apply_all_patterns(
    smiles: str,
    fpps: Dict[str, Union[str, Callable]],
    compound_fingerprint: Dict[str, Union[bool, str, float]]
) -> None:
    """
    Apply SMARTS patterns and molecular analyses to a SMILES string.
    
    Modifies the compound_fingerprint dictionary in-place by:
    - Computing molecular composition and bond types
    - Applying each pattern in the fingerprints dictionary
    - For each pattern, storing the maximum value found
    
    Args:
        smiles: SMILES string representing the molecule to analyze
        fpps: Dictionary of fingerprint patterns to apply
        compound_fingerprint: Dictionary to update with results (modified in-place)
    
    Note:
        - Handles CX SMILES with multiple stereoisomers via rdMolEnumerator
        - Falls back to single molecule if enumeration fails
        - Skips recomputation if pattern result is already True
        
    Performance Notes:
        - O(400+) pattern iterations per molecule
        - Early exit optimization for boolean patterns
        - ~5-10ms typical execution time per molecule
        - inspect.signature() overhead (~0.01-0.05ms per callable) - optimization candidate
    """
    mol = Chem.MolFromSmiles(smiles)
    try:
        mols = [m for m in align_bundle_coords(rdMolEnumerator.Enumerate(mol))]
    except Exception:
        mols = [mol]
    
    for mol in mols:
        compound_fingerprint["Molecular Composition"] = molecule_composition(mol)
        for key, pattern in fpps.items():
            # Skip if already True (optimization for boolean fingerprints)
            if compound_fingerprint.get(key, False) is True:
                continue
            # Apply pattern and update with max value (handles value fingerprints)
            result = apply_pattern(mol, pattern, compound_fingerprint)
            compound_fingerprint[key] = max(result, compound_fingerprint.get(key, False))
