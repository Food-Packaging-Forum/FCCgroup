"""
Molecular composition analysis utilities.

Provides functions for analyzing the atomic composition and structural
properties of RDKit molecule objects.
"""

from typing import Dict, List

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdDepictor


def molecule_composition(mol: Chem.Mol) -> Dict[str, int]:
    """
    Calculates the atomic composition of a given RDKit molecule.

    Returns a dictionary with atomic symbols as keys and their counts as values.

    Args:
        mol (Chem.Mol): An RDKit molecule object.

    Returns:
        dict: A dictionary where keys are atomic symbols (str) and values are the counts (int)
            of each atom in the molecule, including explicit hydrogens.

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles('CCO')
        >>> molecule_composition(mol)
        {'C': 2, 'O': 1, 'H': 6}
    """
    composition = {}
    for atom in Chem.AddHs(mol).GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in composition:
            composition[symbol] = 0
        composition[symbol] += 1
    return composition


def get_bond_orders(mol: Chem.Mol) -> List[float]:
    """
    Returns a list of bond types for all bonds in the given molecule.

    Args:
        mol (Chem.Mol): An RDKit molecule object.

    Returns:
        List[float]: A list containing the bond type (as a float) for each bond
            in the molecule.

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles('C=C')
        >>> get_bond_orders(mol)
        [2.0]
    """
    return [bond.GetBondTypeAsDouble() for bond in mol.GetBonds()]


def align_bundle_coords(bndl: Chem.MolBundle) -> None:
    """
    Aligns the 2D coordinates of a bundle of molecules based on their Maximum Common Substructure (MCS).

    This function takes a list of RDKit molecule objects, sanitizes each molecule, finds the MCS among them,
    generates 2D coordinates for the MCS, and then aligns the 2D coordinates of each molecule in the bundle
    to match the MCS structure.

    Notes:
        - This function is used to handle CX SMILES that would break when applying a pattern. This occurs
          because molecule properties are not set correctly on enumeration.
        - Method adapted from: [Greg Landrum's post](https://greglandrum.github.io/rdkit-blog/posts/2021-05-13-intro-to-the-molecule-enumerator.html)

    Args:
        bndl (Chem.MolBundle): List of RDKit molecule objects to be aligned.

    Returns:
        None: The function modifies the input molecules in place.
    """
    for m in bndl:
        Chem.SanitizeMol(m)
    mcs = rdFMCS.FindMCS(bndl, completeRingsOnly=True)
    q = Chem.MolFromSmarts(mcs.smartsString)
    rdDepictor.Compute2DCoords(q)
    for m in bndl:
        rdDepictor.GenerateDepictionMatching2DStructure(m, q)
