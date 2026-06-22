"""
SMARTS patterns library for chemical fingerprinting.

This module re-exports the main fingerprints dictionary from smarts_methods.py
and provides utility functions for pattern matching and CX SMILES handling.

Contains ~400+ chemical substructure patterns covering:
- Elemental composition and inorganic compounds
- Hydrocarbons and polycyclic aromatics
- Oxygen, nitrogen, phosphorus, sulfur-containing groups
- Heterocyclic aromatics (triazoles, benzothiazoles, etc.)
- Halogenated compounds (PCBs, PBDEs, etc.)
- Organometallics and element-specific patterns
"""

import inspect

import pandas as pd
from typing import Optional, Union, Callable, Dict

from rdkit import Chem
from rdkit.Chem import rdGeneralizedSubstruct

from ..data.constants import elements


def cx_smarts_query(mol: Chem.Mol, pattern: Union[str, Callable]) -> bool:
    """
    Checks if a given molecule matches a SMARTS pattern using extended query matching.

    Args:
        mol (Chem.Mol): The molecule to be queried.
        pattern (Union[str, Callable]): The SMARTS pattern as a string or a callable.

    Returns:
        bool: True if the molecule matches the SMARTS pattern, False otherwise.
    """
    pattern_mol = Chem.MolFromSmarts(pattern)
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(pattern_mol)
    return rdGeneralizedSubstruct.MolHasSubstructMatch(mol, xqm)

def count_pattern_occurrences(mol: Chem.Mol, pattern: str) -> int:
    """
    Counts the number of occurrences of a specified SMARTS pattern in a given RDKit molecule.

    Args:
        mol (Chem.Mol): The RDKit molecule object in which to search for the pattern.
        pattern (str): The SMARTS pattern string to search for within the molecule.

    Returns:
        int: The number of times the SMARTS pattern is found in the molecule.
    """
    pattern = Chem.MolFromSmarts(pattern)
    return len(mol.GetSubstructMatches(pattern))


def apply_pattern(mol: Chem.Mol, pattern: Union[str, Callable], row: Optional[pd.Series] = None) -> Union[bool, int, float]:
    """
    Applies a SMARTS pattern or a callable function to a molecule to determine if it matches.

    Supports both SMARTS string patterns and callable functions. Callable functions can accept
    either one argument (mol) or two arguments (mol, row). This enables complex filtering logic
    that depends on computed features in the DataFrame row.

    Parameters:
        mol: The RDKit molecule object to be checked
        pattern: Either a SMARTS string or callable function
            - If a string: interpreted as a SMARTS pattern
            - If callable: should accept (mol) or (mol, row)
        row: Optional pandas Series with computed features for callable patterns

    Returns:
        bool or numeric value: True/False for boolean patterns, or numeric value for value fingerprints
    
    Note:
        Uses inspect.signature() to determine callables' parameter count. Callables with 2 parameters
        will be called with (mol, row).
        
    Performance Notes:
        - SMARTS string evaluation: ~1-5ms per pattern (RDKit optimized)
        - Callable overhead: Additional ~0.01-0.05ms for signature() inspection
        - Called ~400+ times per molecule in apply_all_patterns()
        - Optimization opportunity: Cache signatures for lambda patterns
    """
    if isinstance(pattern, str):
        pattern_mol = Chem.MolFromSmarts(pattern)
        return mol.HasSubstructMatch(pattern_mol)
    else:
        # Use inspect.signature to determine parameter count
        sig = inspect.signature(pattern)
        param_count = len(sig.parameters)
        
        if param_count == 1:
            return pattern(mol)
        else:
            return pattern(mol, row if row is not None else {})



def is_aromatic(mol: Chem.Mol) -> bool:
    """
    Determines if a given molecule, represented by its SMILES string, is aromatic.

    Args:
        mol (rdkit.Chem.Mol): The molecule to be checked, represented as an RDKit Mol object.

    Returns:
        bool: True if the molecule is aromatic, False otherwise.

    Notes:
    - A molecule is considered aromatic if all its atoms are aromatic.
    """
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic():
            return False
    return True

def is_salt(mol: Chem.Mol, avoid_single_molecules=False) -> bool:
    """
    Determines whether the given molecule contains both positive and negative ions.

    The function splits the molecule into fragments and checks if any fragment contains a positive ion ("+")
    and any fragment contains a negative ion ("-") in its SMILES representation. If both are present, it returns True.

    Args:
        mol (Chem.Mol): The RDKit molecule to analyze.
        avoid_single_molecules (bool, optional): If True, the function returns False for molecules with fewer than two fragments. Defaults to False.

    Returns:
        bool: True if both positive and negative ions are present in the molecule, False otherwise.
    """
    has_positive_ion = False
    has_negative_ion = False
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if avoid_single_molecules and len(fragments) < 2:
        return False
    for m in fragments:
        smiles = Chem.MolToSmiles(m)
        if not has_positive_ion and "+" in smiles:
            has_positive_ion = True
        elif not has_negative_ion and "-" in smiles:
            has_negative_ion = True
        if has_positive_ion and has_negative_ion:
            return True
    return False

def is_single_element(mol: Chem.Mol, allow_repeat: bool = False) -> bool:
    """
    Determines if a molecule consists of a single type of element.

    Args:
        mol (Chem.Mol): The molecule to analyze, represented as an RDKit Mol object.
        allow_repeat (bool, optional): If True, allows the same element to appear multiple times 
            in the molecule. Defaults to False.

    Returns:
        bool: True if the molecule contains only one type of element, False otherwise.
    """
    elements_set = set()
    for atom in Chem.AddHs(mol).GetAtoms():
        symbol = atom.GetSymbol()
        # If we don't allow repeat and atom is already in elements, return False
        if not allow_repeat and symbol in elements_set:
            return False
        # else, add the atom to the set of elements
        elements_set.add(symbol)
        # return False if we have more than one element
        if len(elements_set) > 1:
            return False
    return True

__all__ = ['fingerprints', 'cx_smarts_query', 'count_pattern_occurrences', 'apply_pattern']

# Main fingerprints dictionary: 400+ chemical substructure patterns
fingerprints: Dict[str, Union[str, Callable]] = {
    'Contains Two Charges': lambda x: is_salt(mol=x, avoid_single_molecules=False),
    'Salt': lambda x: is_salt(x, avoid_single_molecules=True),

    'Contains metal': lambda x, row: len(set(elements["metals"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains alkali metal': lambda x, row: len(set(elements["alkali_metals"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains alkaline metal': lambda x, row: len(set(elements["alkaline_metals"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains transition metal': lambda x, row: len(set(elements["transition_metals"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains posttransition metal': lambda x, row: len(set(elements["posttransition_metals"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains lanthanoid metal': lambda x, row: len(set(elements["lanthanoids"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains actinoid metal': lambda x, row: len(set(elements["actinoids"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains Cr': lambda x, row: row["Molecular Composition"].get("Cr", 0),
        'Contains Cr (VI)': lambda x: count_pattern_occurrences(x, '[#24+6]'),
        'Contains Ni': lambda x, row: row["Molecular Composition"].get("Ni", 0),
        'Contains Cd': lambda x, row: row["Molecular Composition"].get("Cd", 0),
        'Contains Hg': lambda x, row: row["Molecular Composition"].get("Hg", 0),
        'Contains Pb': lambda x, row: row["Molecular Composition"].get("Pb", 0),
        'Contains toxic heavy metal (As, Cd, Cr, Pb, Hg, Ni)': lambda x, row: any(row["Molecular Composition"].get(elem, 0) for elem in ["As", "Cd", "Cr", "Pb", "Hg", "Ni"]),

    'Contains noble gas': lambda x, row: len(set(elements["noble_gases"]).intersection(set(row["Molecular Composition"]))) > 0,
    
    'Contains metalloid': lambda x, row: len(set(elements["metalloids"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains As': lambda x, row: row["Molecular Composition"].get("As", 0),

    'Contains nonmetal': lambda x, row: len(set(elements["non_metals"]).intersection(set(row["Molecular Composition"]))) > 0,
        'Contains H': lambda x, row: row["Molecular Composition"].get("H", 0),
        'Contains halogen': lambda x, row: len(set(elements["halogens"]).intersection(set(row["Molecular Composition"]))) > 0,
    'Contains C': lambda x, row: row["Molecular Composition"].get("C", 0),
        'Contains C~H': lambda x: count_pattern_occurrences(Chem.AddHs(x), "[#6]~[#1]"),
        'Contains C~C or C~H': lambda x: count_pattern_occurrences(Chem.AddHs(x), "[#6]~[#1,#6]"),
        'Organic C': lambda x, row: row["Contains C"] != 0, 
        'Organic (CH bond)': lambda x, row: row["Contains C~H"] != 0, 
        'Organic (only nonmetals)': lambda x, row: row["Organic C"] and not row["Contains metal"] 
                                    and not row["Contains metalloid"],
        'Contains aliphatic C': lambda x: count_pattern_occurrences(x, "C"),
        'Contains aromatic C': lambda x: count_pattern_occurrences(x, "c"),
        'Contains aromatic': lambda x: count_pattern_occurrences(x, "a"),
        'Contains C-C': lambda x: count_pattern_occurrences(x, "C-C"),
        'Contains C=C': lambda x: count_pattern_occurrences(x, "C=C"),
        'Contains C#C': lambda x: count_pattern_occurrences(x, "C#C"),
        'Contains branched alkyl chain': lambda x: count_pattern_occurrences(x, "[#6]-C(-C)(-C)") > 0,
        'Contains only aromatic atoms': is_aromatic,
        'Contains only ring atoms': lambda x: not apply_pattern(x, "[*&R0]"),
        'Contains ring': "[*&R]",
        'Contains fused ring': "[*&R2]",
        'Contains benzene': lambda x: count_pattern_occurrences(x, 'c1ccccc1'),
        'Contains benzene (strict)': lambda x: count_pattern_occurrences(x, '[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1'),
        'Contains styrene': lambda x: count_pattern_occurrences(x, '[c&H1]1[c&H1][c&H1][c&H1][c&H1]c1-[C&H1]=[C&H2]'),
        'Contains styrene moiety': lambda x: count_pattern_occurrences(x, '[c&H1]1[c&H1][c&H1][c&H1][c&H1]c1-[C&H1]~[C&H2]'),
        'Contains naphthalene': lambda x: count_pattern_occurrences(x, 'c1c2ccccc2ccc1'),

    'Elements (loose)': lambda x, row: len(row["Molecular Composition"].keys()) == 1,
    'Elements (strict)': lambda x, row: len(row["Molecular Composition"].keys()) == 1 and row["Molecular Composition"][list(row["Molecular Composition"].keys())[0]] == 1,
    'Inorganic compound': lambda x, row: row["Contains C"] == 0,
        'Inorganic salt (no C)': lambda x, row: row["Inorganic compound"] and row["Salt"],
        'Inorganic other': lambda x, row: row["Inorganic compound"] and not row["Elements (loose)"] and not row["Inorganic salt (no C)"],
    'Contains Hydrocarbons': lambda x, row: {'C','H'}.issubset(row["Molecular Composition"]),
        'Hydrocarbons': lambda x, row: set(row["Molecular Composition"]) == {"C", "H"},
        'Aliphatic hydrocarbons': lambda x, row: row["Hydrocarbons"] and row["Contains aromatic C"] == 0,
        'Fully aromatic hydrocarbons': lambda x, row: row["Hydrocarbons"] and row["Contains only aromatic atoms"],
        'Aromatic hydrocarbons': lambda x, row: row["Hydrocarbons"] and row["Contains aromatic C"] != 0,
        'Alkynes': lambda x, row: row["Aliphatic hydrocarbons"] and row["Contains C#C"] != 0,
        'Alkenes': lambda x, row: row["Aliphatic hydrocarbons"] and row["Contains C#C"] == 0 and row["Contains C=C"] != 0,
        'Alkanes': lambda x, row: row["Aliphatic hydrocarbons"] and row["Contains C#C"] == 0 and row["Contains C=C"] == 0 and row["Contains aliphatic C"] != 0,
        'Branched alkanes': lambda x, row: row["Alkanes"] and row["Contains branched alkyl chain"],
        'Linear alkanes': lambda x, row: row["Alkanes"] and not row["Contains branched alkyl chain"],
        'Biphenyls/Terphenyls':  lambda x, row: row["Hydrocarbons"] and apply_pattern(x,'[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1-[c&R1]2[c&R1][c&R1][c&R1][c&R1][c&R1]2'),
        'Biphenyl/Terphenyl derivatives': '[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1-[c&R1]2[c&R1][c&R1][c&R1][c&R1][c&R1]2',
        'PAH derivatives': "c1[c&!R1&R](c)[c&!R1&R]([#6])ccc1",
        'PAH derivatives hydrocarbon': lambda x, row: row["Hydrocarbons"] and row['PAH derivatives'],
        'PAHs': lambda x, row: row["Hydrocarbons"] and row['PAH derivatives'] and row["Contains only ring atoms"],
        'Benzoids': lambda x, row: row["Hydrocarbons"] and not row["PAHs"] and row["Contains benzene"] > 0,
        'Abietic acid derivatives (tricyclic diterpenes backbone)': '[#6]12~[#6](~[#6]~[#6]~[#6]~[#6]~1)~[#6]~[#6]~[#6]3~[#6]~2~[#6]~[#6]~[#6]~[#6]~3',
        'Steroid-like compounds': '[#6]~12~[#6](~[#6]~[#6]~[#6]~[#6]~1)~[#6]~[#6]~[#6]~3~[#6]~2~[#6]~[#6]~[#6]4~[#6](~[#6]~[#6]~[#6]~4)~3',

    'Contains Organo O': lambda x, row: {'C','O'}.issubset(row["Molecular Composition"]),
        'Contains C~O': lambda x: count_pattern_occurrences(x, '[#8]~[#6]'),
        'Organo O (strict)': lambda x, row: {'C','O'}.issubset(row["Molecular Composition"].keys()) and
                                    set(row["Molecular Composition"].keys()).issubset({"C", "O", "H"}),
        'Contains O': lambda x, row: row["Molecular Composition"].get("O", 0),
        'Contains C=O': lambda x: count_pattern_occurrences(x, '[#8]=[#6]'),
        'Contains aldehyde': lambda x: count_pattern_occurrences(Chem.AddHs(x), '[#8]=[#6]([#1])([#6,#1])'),
        'Contains ketone': lambda x: count_pattern_occurrences(x, '[#8]=[#6]([#6])([#6])'),
        'Contains alcohol': lambda x: count_pattern_occurrences(x, '[#8&H1][#6]'),
        'Contains ether': lambda x: count_pattern_occurrences(x, '[#6][#8][#6]'),
        'Contains epoxide': lambda x: count_pattern_occurrences(x, '[#6]1[#8][#6]1'),
        'Contains carboxylic acid': lambda x: count_pattern_occurrences(x, '[#8&H1][#6]=[#8]'),
        'Contains carboxylic acid derivatives': lambda x: count_pattern_occurrences(x, '[#8][#6]=[#8]'),
        'Contains ester': lambda x: count_pattern_occurrences(x, '[#6][#8][#6]=[#8]'),
        'Contains benzoic acid': lambda x: count_pattern_occurrences(x, '[#8&H1][#6](=[#8])[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1'),
        'Contains benzoic acid derivatives': lambda x: count_pattern_occurrences(x, '[#8][#6](=[#8])c1ccccc1'),
        'Contains benzoic ester': lambda x: count_pattern_occurrences(x, '[#6][#8][#6](=[#8])[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1'),
        'Contains O-heterocycle': lambda x: count_pattern_occurrences(x, '[#8&R]'),
        'Contains phenol': lambda x: count_pattern_occurrences(x, 'c1cc([#8&H1])ccc1'),
        'Contains phenol derivatives': lambda x: count_pattern_occurrences(x, 'c1cc([#8])ccc1'),
        'Contains peroxide': lambda x: count_pattern_occurrences(x, '[#8][#8][#6]'),
        'Contains ethylene glycol': lambda x: count_pattern_occurrences(x, '[#8]-[C&H2]-[C&H2]-[#8]'),
        'Contains butylene glycol': lambda x: count_pattern_occurrences(x, '[#8]-[C&H2]-[C&H2]-[C&H2]-[C&H2]-[#8]'),
        'Contains adipic acid': lambda x: count_pattern_occurrences(x, '[#8]C(=[#8])[C&H2][C&H2][C&H2][C&H2]C(=[#8])[#8]'),
        'Contains lactic acid': lambda x: count_pattern_occurrences(x, '[#8]=[#6]([#8])C([C&H3])[#8]'),
        'Contains terephthalic acid': lambda x: count_pattern_occurrences(x, '[#8][#6](=[#8])c1[c&H1][c&H1]c([#6](=[#8])[#8])[c&H1][c&H1]1'),
        'Contains naphthalene-2,6-dicarboxylic acid': lambda x: count_pattern_occurrences(x, 'c1c2cc([#6](=[#8])[#8])ccc2cc([#6](=[#8])[#8])c1'),
        'Alkyl phenols': lambda x, row: cx_smarts_query(x, 'C*.[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1[#8&H1] |m:1:2.3.4|') and 
                row['Organo O (strict)'] and row['Contains O'] == 1 and row['Contains C=C'] == 0 and row['Contains C#C'] == 0 and row['Contains aromatic C'] == 6,
        'Alkyl phenol derivatives':  lambda x, row: cx_smarts_query(x, 'C*.[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1[#8] |m:1:2.3.4|') and 
                row['Organo O (strict)'] and row['Contains O'] == 1 and row['Contains C=C'] == 0 and row['Contains C#C'] == 0 and row['Contains aromatic C'] == 6,
        'Hindered phenols': 'c1([#8&H1])c(-[#6]([#6])([#6]))cc(*)cc1',
        'Hindered phenol derivatives': 'c1([#8])c(-[#6]([#6])([#6]))cc(*)cc1',
        'Hindered phenol backbone': '[#6]~1(~[#8])~[#6](~[#6](~[#6])(~[#6]))~[#6]~[#6](~*)~[#6]~[#6]1',
        'Hindered quinone derivatives': '[#6]1(=[#8])-[#6](-[#6]([#6])([#6]))=[#6]-[#6](~*)-[#6]=[#6]-1',
        'Bisphenol M': '[c&R1]1[c&R1]([#8&H1])[c&R1][c&R1][c&R1]([c&R1]1)-C-c1cccc(c1)-C-[c&R1]2[c&R1][c&R1][c&R1]([#8&H1])[c&R1][c&R1]2',
        'Bisphenol M derivatives': '[c&R1]1[c&R1]([#8])[c&R1][c&R1][c&R1]([c&R1]1)-C-c1cccc(c1)-C-[c&R1]2[c&R1][c&R1][c&R1]([#8])[c&R1][c&R1]2',
        'Bisphenol P': '[c&R1]1[c&R1]([#8&H1])[c&R1][c&R1][c&R1]([c&R1]1)-C-c1ccc(cc1)-C-[c&R1]2[c&R1][c&R1][c&R1]([#8&H1])[c&R1][c&R1]2',
        'Bisphenol P derivatives': '[c&R1]1[c&R1]([#8])[c&R1][c&R1][c&R1]([c&R1]1)-C-c1ccc(cc1)-C-[c&R1]2[c&R1][c&R1][c&R1]([#8])[c&R1][c&R1]2',
        'Bisphenol bridge': '[c&R1]1[c&R1]([#8&H1])[c&R1][c&R1][c&R1]([c&R1]1)-[*]-[c&R1]2[c&R1][c&R1][c&R1]([#8&H1])[c&R1][c&R1]2',
        'Bisphenol bridge derivatives': '[c&R1]1[c&R1]([#8])[c&R1][c&R1][c&R1]([c&R1]1)-[*]-[c&R1]2[c&R1][c&R1][c&R1]([#8])[c&R1][c&R1]2',
        'Bisphenol CS bridge': '[c&R1]1[c&R1]([#8&H1])[c&R1][c&R1][c&R1]([c&R1]1)-[C&R0,S&R0]-[c&R1]2[c&R1][c&R1][c&R1]([#8&H1])[c&R1][c&R1]2',
        'Bisphenol CS bridge derivatives': '[c&R1]1[c&R1]([#8])[c&R1][c&R1][c&R1]([c&R1]1)-[C&R0,S&R0]-[c&R1]2[c&R1][c&R1][c&R1]([#8])[c&R1][c&R1]2',
        'Bisphenols': lambda x, row: row['Bisphenol M'] or row['Bisphenol P'] or row['Bisphenol bridge'],
        'Bisphenol derivatives': lambda x, row: row['Bisphenol M derivatives'] or row['Bisphenol P derivatives'] or row['Bisphenol bridge derivatives'],
        'Benzylketone derivatives': 'c1ccccc1-[#6](=[#8])([#6])',
        '2-hydroxybenzophenones': 'c1ccccc1-[#6](=[#8])-c2c(~[#8])cccc2',
        'Acetophenone derivatives': '[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1C(=[#8])([C&H3])',
        'Benzoquinone 1,2 backbone': '[#6&R1]1=[#6&R1][#6&R1]=[#6&R1][#6&R1][#6&R1]1=[#8]',
        'Benzoquinone 1,4 backbone': '[#6&R1]1=[#6&R1][#6&R1][#6&R1]=[#6&R1][#6&R1]1=[#8]',
        'Benzoquinone backbone': lambda x, row: row['Benzoquinone 1,2 backbone'] or row['Benzoquinone 1,4 backbone'],
        'Benzophenone derivatives': '[c&R1]1[c&R1][c&R1][c&R1][c&R1][c&R1]1-[#6](=[#8])-[c&R1]2[c&R1][c&R1][c&R1][c&R1][c&R1]2',
        'Benzoquinone 1,2 derivatives': '[#6&R1]1=[#6&R1][#6&R1]=[#6&R1][#6&R1](=[#8])[#6&R1]1=[#8]',
        'Benzoquinone 1,4 derivatives': '[#6&R1]1=[#6&R1][#6&R1](=[#8])[#6&R1]=[#6&R1][#6&R1]1=[#8]',
        'Benzoquinone derivatives': lambda x, row: row['Benzoquinone 1,2 derivatives'] or row['Benzoquinone 1,4 derivatives'],
        'Anthraquinone derivatives': 'c1cccc2c1[#6](=[#8])c3ccccc3[#6](=[#8])2',
        'Xanthenes derivatives': 'c1cccc2c1[#6]c3ccccc3[#8]2',
        'Xanthones derivatives': 'c1cccc2c1[#6](=[#8])c3ccccc3[#8]2',
        'Adipic acid esters': '[#8]C(=[#8])[C&H2][C&H2][C&H2][C&H2]C(=[#8])[#8]',
        'Sebacic acid esters': '[#8]C(=[#8])[C&H2][C&H2][C&H2][C&H2][C&H2][C&H2]C(=[#8])[#8]',
        'Citric acid esters': '[#8]C(=[#8])[C&H2]C(-[#8])(C(=[#8])[#8])[C&H2]C(=[#8])[#8]',
        'Pentaerythritol derivatives': lambda x: apply_pattern(Chem.AddHs(x),'[C,#1][#8][C&H2]C([C&H2][#8][C,#1])([C&H2][#8][C,#1])([C&H2][#8][C,#1])'),
        'Acrylic acid esters': '[#6&H2]=[#6&H1]-[#6](=[#8])-[#8]',
        'Methacrylic acid esters': '[#6&H2]=[#6](-[#6&H3])-[#6](=[#8])-[#8]',
        'Glycerol derivatives': lambda x: apply_pattern(Chem.AddHs(x),'[C,#1][#8][C&H2][C&H1]([C&H2][#8][C,#1])[#8][C,#1]'),
        'Glycerol fatty acid esters': lambda x, row: row["Glycerol derivatives"] and row["Contains O"] <= 6 and row["Organo O (strict)"],
        'Maleic and fumaric acid esters': '[#8]C(=[#8])-[C&H1]=[C&H1]-C(=[#8])[#8]',
        'Stearic acid esters': '[C&H3][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H2]C(=[#8])[#8]',
        'Triglyceride': 'C(=[#8])[#8][C&H2][C&H1]([C&H2][#8]C(=[#8]))[#8]C(=[#8])',
        'Cinnamates': 'c1ccccc1-[#6]=[#6][#6](=[#8])[#8]',
        'Benzoates': '[#8][#6](=[#8])c1[c&H1][c&H1][c&H1][c&H1][c&H1]1',
        'Benzoates derivatives': '[#8]C(=[#8])c1ccccc1',
        'Parabens': '[#8&H1]c1[c;H1][c;H1]c(C(=[#8])[#8])[c;H1][c;H1]1',
        'Parabens derivatives': '[#8]c1[c;H1][c;H1]c(C(=[#8])[#8])[c;H1][c;H1]1',
        'Salicylates': '[#8&H1][#6](=[#8])c1[c&H1][c&H1][c&H1][c&H1]c1([#8])',
        'Salicylates derivatives': '[#8][#6](=[#8])c1[c&H1][c&H1][c&H1][c&H1]c1([#8])',
        'Ortho-phthalates': '[#8][#6](=[#8])c1[c&H1][c&H1][c&H1][c&H1]c1[#6](=[#8])[#8]',
        'Ortho-phthalates derivatives (no rings)': '[#8][#6](=[#8])c1ccccc1[#6](=[#8])[#8]',
        'Cyclohexane-1,2-dicarboxylates': '[#8][C](=[#8])[C&H1]1[C&H2][C&H2][C&H2][C&H2][C&H1]1[C](=[#8])[#8]',
        'Isophthalates': '[#8][#6](=[#8])c1[c&H1][c&H1][c&H1]c([#6](=[#8])[#8])[c&H1]1',
        'Isophthalates derivatives (no rings)': '[#8][#6](=[#8])c1cccc([#6](=[#8])[#8])c1',
        'Terephthalates': '[#8][#6](=[#8])c1[c&H1][c&H1]c([#6](=[#8])[#8])[c&H1][c&H1]1',
        'Terephthalates derivatives (no rings)': '[#8][#6](=[#8])c1ccc([#6](=[#8])[#8])cc1',
        'Trimellitates':  '[#8][#6](=[#8])c1[c&H1][c&H1]c([#6](=[#8])[#8])[c&H1]c1([#6](=[#8])[#8])',
        'Trimesitates': '[#8][#6](=[#8])c1[c&H1]c([#6](=[#8])[#8])[c&H1]c([#6](=[#8])[#8])[c&H1]1',
        'Hemimellitates':  '[#8][#6](=[#8])c1[c&H1][c&H1][c&H1]c([#6](=[#8])[#8])c1([#6](=[#8])[#8])',
        'Pyromellitates':  '[#8][#6](=[#8])c1[c&H1]c([#6](=[#8])[#8])c([#6](=[#8])[#8])[c&H1]c1([#6](=[#8])[#8])',
        'Benzofurans': '[c&R1]1[c&R1][c&R1][c&R1][c&R2]2[o&R1][c&R1,o&R1][c&R1,o&R1][c&R2]21',
        'Benzofuran backbone': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#6&R2]2~[#8&R1]~[#6&R1,#8&R1]~[#6&R1,#8&R1]~[#6&R2]~2~1',
        'Isobenzofuran backbone': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#6&R2]2~[#6&R1,#8&R1]~[#8&R1]~[#6&R1,#8&R1]~[#6&R2]~2~1',
        'Benzofuran derivatives': lambda x, row: row['Benzofuran backbone'] or row['Isobenzofuran backbone'],
        'Aromatic carboxylic acids': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains carboxylic acid"] != 0 and row["Contains O"] == 2,
        'Aromatic carboxylic acid esters and salts': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains carboxylic acid"] == 0 and row["Contains carboxylic acid derivatives"] != 0 and row["Contains O"] == 2,
        'Aromatic carbonyls (loose)': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    (row["Contains aldehyde"] + row["Contains ketone"]) == row["Contains O"],
        'Aromatic aldehydes': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains aldehyde"] != 0 and row["Contains O"] == 1,
        'Aromatic ketones': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains ketone"] != 0 and row["Contains O"] == 1,
        'Aromatic ethers or alcohols (loose)': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and
                                    (row["Contains ether"] + row["Contains alcohol"]) == row["Contains O"],
        'Aromatic alcohols': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains alcohol"] != 0 and row["Contains O"] == 1,
        'Aromatic ethers': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains ether"] != 0 and row["Contains O"] == 1 and row["Contains O-heterocycle"] == 0,
        'Aliphatic carboxylic acids esters and salts (loose)': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and
                                    2 *row["Contains carboxylic acid derivatives"] ==  row["Contains O"],
        'Aliphatic monocarboxylic acids': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains carboxylic acid"] != 0 and row["Contains O"] == 2,
        'Aliphatic monocarboxylic acid derivatives': lambda x, row: row["Organo O (strict)"] and not row["Contains aromatic C"] and row["Contains carboxylic acid derivatives"] and row["Contains O"] ==  2*row["Contains carboxylic acid derivatives"],
        'Aliphatic monocarboxylic acids esters and salts': lambda x, row: row["Organo O (strict)"] and not row["Contains aromatic C"] and row["Contains carboxylic acid derivatives"] == 1 and row["Contains O"] == 2,
        'Aliphatic dicarboxylic acids esters and salts': lambda x, row: row["Organo O (strict)"] and not row["Contains aromatic C"] and row["Contains carboxylic acid derivatives"] == 2 and row["Contains O"] == 4,
        'Aliphatic tricarboxylic acids esters and salts': lambda x, row: row["Organo O (strict)"] and not row["Contains aromatic C"] and row["Contains carboxylic acid derivatives"] == 3 and row["Contains O"] == 6,
        'Aliphatic tetracarboxylic acids esters and salts': lambda x, row: row["Organo O (strict)"] and not row["Contains aromatic C"] and row["Contains carboxylic acid derivatives"] == 4 and row["Contains O"] == 8,
        'Aliphatic carbonyls (loose)': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and
                                    (row["Contains aldehyde"] + row["Contains ketone"]) == row["Contains O"],
        'Aliphatic aldehydes': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains aldehyde"] != 0 and row["Contains O"] == 1,
        'Aliphatic ketones': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains ketone"] != 0 and row["Contains O"] == 1,
        'Aliphatic ethers or alcohols (loose)': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and
                                    (row["Contains ether"] + row["Contains alcohol"]) == row["Contains O"],
        'Aliphatic alcohols': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains alcohol"] != 0 and row["Contains O"] == 1,
        'Aliphatic ethers': lambda x, row: row["Organo O (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains ether"] != 0 and row["Contains O"] == 1,
    'Contains Organo N': lambda x, row: {'C','N'}.issubset(row["Molecular Composition"]),
        'Contains C~N': lambda x: count_pattern_occurrences(x, '[#7]~[#6]'),
        'Organo N (strict)': lambda x, row: {'C','N'}.issubset(row["Molecular Composition"].keys()) and
                                set(row["Molecular Composition"].keys()).issubset({"C", "N", "H"}),
        'Contains N':  lambda x, row: row["Molecular Composition"].get("N", 0),
        'Contains N-heterocycle': lambda x: count_pattern_occurrences(x, '[#7&R]'),
        'Contains amine (loose)': lambda x: count_pattern_occurrences(x, '[#7&X3]'),
        'Contains amine': lambda x: count_pattern_occurrences(Chem.AddHs(x), '[N&X3]([#6])([#1,#6])([#1,#6])'),
        'Contains primary amine': lambda x: count_pattern_occurrences(x, '[N&H2&X3]([#6])'),
        'Contains secondary amine': lambda x: count_pattern_occurrences(x, '[N&H1&X3]([#6])([#6])'),
        'Contains tertiary amine': lambda x: count_pattern_occurrences(x, '[N&X3]([#6])([#6])([#6])'),
        'Contains p-phenylenediamine': lambda x: count_pattern_occurrences(x, '[N&X3]c1ccc([N&X3])cc1'),
        'Contains hydroxylamine': lambda x: count_pattern_occurrences(x, '[N&X3](-[#6])(-[#6])(-[#8])'),
        'Contains quarternary ammonium ion': lambda x: apply_pattern(Chem.AddHs(x), '[#7+]([!#1])([!#1])([!#1])([!#1])'),
        'Contains hydrazine': lambda x: count_pattern_occurrences(x, '[#7&X3][#7&X3]'),
        'Contains imine': lambda x: count_pattern_occurrences(Chem.AddHs(x), '[#1,#6]-[#7]=[#6](-[#1,#6])(-[#6,#1])'),
        'Contains isocyanide': lambda x: count_pattern_occurrences(x, '[#6-]#[#7+&X2]'),
        'Contains cyanate': lambda x: count_pattern_occurrences(x, '[#8]-[#6]#[#7]'),
        'Contains isocyanate': lambda x: count_pattern_occurrences(x, '[#8]=[#6]=[#7&X2]'),
        'Contains nitrile': lambda x: count_pattern_occurrences(x, '[#6]#[#7&X1]'),
        'Contains amide': lambda x: count_pattern_occurrences(Chem.AddHs(x), '[#8]=[#6](-[#7&X3])(-[#1,#6])'),
        'Contains formamidine': lambda x: count_pattern_occurrences(x, '[#7&X2]=[#6]-[#7&X3]'),
        'Contains urea': lambda x: count_pattern_occurrences(x, '[#8]=[#6](-[#7&X3])(-[#7&X3])'),
        'Contains nitroso': lambda x: count_pattern_occurrences(x, '[#6][#7&X2]=[#8]'),
        'Contains nitrosamine': lambda x: count_pattern_occurrences(x, '[#6][#7][#7&X2]=[#8]'),
        'Contains nitro': lambda x: count_pattern_occurrences(x, '[#6][#7&X3](=[#8])([#8])'),
        'Contains azo': lambda x: count_pattern_occurrences(x, '[#6]-[#7]=[#7]-[#6]'),
        'Contains amino acid': lambda x: count_pattern_occurrences(x, '[#6](-[#8&H1])(=[#8])[#6]([#7&H1&R3])'),
        'Contains amino acid derivatives': lambda x: count_pattern_occurrences(x, '[#6](-[#8])(=[#8])-[#6](-[#7])'),
        'Contains caprolactam': lambda x: count_pattern_occurrences(x, 'O=C[C&H2][C&H2][C&H2][C&H2][C&H2]N'),
        'Contains hexamethylenediamine': lambda x: count_pattern_occurrences(x, 'N[C&H2][C&H2][C&H2][C&H2][C&H2][C&H2]N'),
        'Primary aromatic amines': 'c1ccccc1[N&H2&X3]',
        'Primary aromatic amines (any)': 'a[N&H2&X3]',
        'Secondary aromatic amines': 'c1ccccc1[N&H1&X3][#6]',
        'Secondary aromatic amines (any)': 'a[N&H1&X3][#6]',
        'Tertiary aromatic amines': 'c1ccccc1[N&X3]([#6])([#6])',
        'Tertiary aromatic amines (any)': 'a[N&X3]([#6])([#6])',
        'Imidazole': 'n1cncc1',
        'Imidazole derivatives': '[#7]1~[#6]~[#7]~[#6]~[#6]~1',
        'Pyrazoles': 'c1nncc1',
        'Pyrazole derivatives': '[#6]1~[#7]~[#7]~[#6]~[#6]~1',        
        'Triazines (1,2,3)': 'c1nnncc1',
        'Triazines (1,2,4)': 'c1nncnc1', 
        'Triazines (1,3,5)': 'c1ncncn1',
        'Triazines': lambda x, row: row["Triazines (1,2,3)"] or row["Triazines (1,2,4)"] or row["Triazines (1,3,5)"],
        'Triazines (1,2,3) derivatives': '[#6]1~[#7]~[#7]~[#7]~[#6]~[#6]~1',
        'Triazines (1,2,4) derivatives': '[#6]1~[#7]~[#7]~[#6]~[#7]~[#6]~1',
        'Triazines (1,3,5) derivatives': '[#6]1~[#7]~[#6]~[#7]~[#6]~[#7]~1',
        'Triazines derivatives': lambda x, row: row["Triazines (1,2,3) derivatives"] or row["Triazines (1,2,4) derivatives"] or row["Triazines (1,3,5) derivatives"],
        '2-hydroxyphenyl-diphenyl-triazines': 'n1c(-c2c(~[#8])cccc2)nc(-c3ccccc3)nc1-c4ccccc4',
        'Triazoles (1,2,3) derivatives': '[#6]1~[#7]~[#7]~[#7]~[#6]~1',
        'Triazoles (1,2,4) derivatives': '[#6]1~[#7]~[#7]~[#6]~[#7]~1',
        'Triazoles derivatives': lambda x, row: row["Triazoles (1,2,3) derivatives"] or row["Triazoles (1,2,4) derivatives"], 
        'Benzotriazines (1,2,3)': 'c1ccc2nnncc2c1',
        'Benzotriazines (1,2,4)': 'c1ccc2nncnc2c1',
        'Benzotriazines': lambda x, row: row["Benzotriazines (1,2,3)"] or row["Benzotriazines (1,2,4)"],
        'Benzotriazines (1,2,3) derivatives': '[#6]1~[#6]~[#6]~[#6]2~[#7]~[#7]~[#7]~[#6]~[#6]~2~[#6]~1',
        'Benzotriazines (1,2,4) derivatives': '[#6]1~[#6]~[#6]~[#6]2~[#7]~[#7]~[#6]~[#7]~[#6]~2~[#6]~1',
        'Benzotriazines derivatives': lambda x, row: row["Benzotriazines (1,2,3) derivatives"] or row["Benzotriazines (1,2,4) derivatives"],
        'Benzotriazoles': 'c1ccc2nnnc2c1',
        '2-hydroxyphenylbenzotriazoles': 'c1ccc2nn(-c3c(~[#8])cccc3)nc2c1',
        'Benzotriazoles derivatives': '[#6]1~[#6]~[#6]~[#6]2~[#7]~[#7]~[#7]~[#6]~2~[#6]~1',
        'Benzooxazoles': '[c&R1]1[c&R1][c&R1][c&R1][c&R2]2[n&R1][c&R1][o&R1][c&R2]21',
        'Benzooxazole derivatives': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#6&R2]2~[#7&R1]~[#6&R1]~[#8&R1]~[#6&R2]~2~1',
        'Benzoisooxazoles': '[c&R1]1[c&R1][c&R1][c&R1][c&R2]2[n&R1,o&R1][n&R1,o&R1][c&R1][c&R2]21',
        'Benzoisooxazole derivatives': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#6&R2]2~[#7&R1,#8&R1]~[#7&R1,#8&R1]~[#6&R1]~[#6&R2]~2~1',
        'Tetramethylpiperidine derivatives': '[#6]~1(~[#6])(~[#6])~[#7]~[#6](~[#6])(~[#6])~[#6]~[#6]~[#6]~1',
        'Diphenylformamidines': 'c1ccccc1-[#7&X2]=[#6]-[#7&X3]-c2ccccc2',
        'Oxalanilides': 'c1ccccc1-[#7][#6](=[#8])[#6](=[#8])[#7]-c2ccccc2',
        'Acridanone': 'c1cccc2c1[#7]c3ccccc3[#6](=[#8])2',
        'Acridanone derivatives': '[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]~1~[#7]~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6](=[#8])~2',
        'Aromatic amides (loose)': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "N", "O", "H"}) and 
                                    row["Contains aromatic C"] != 0 and row["Contains amide"] == row["Contains N"] == row["Contains O"],
        'Aromatic amides': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "N", "O", "H"})
                                    and row["Contains aromatic C"] != 0 and row["Contains amide"] != 0 
                                    and row["Contains O"] == 1 and row["Contains N"] == 1,
        'Aromatic nitriles': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains nitrile"] != 0 and row["Contains N"] == 1,
        'Aromatic amines (loose)': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] != 0 and 
                                    (row["Contains primary amine"] + row["Contains secondary amine"] + row["Contains tertiary amine"] == row["Contains N"]),
        'Aromatic amines': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains amine"] != 0 and row["Contains N"] == 1,
        'Aromatic imines': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains imine"] != 0 and row["Contains N"] == 1,
        'Aliphatic amides (loose)': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "N", "O", "H"}) and 
                                    {"C", "N", "O"}.issubset(row["Molecular Composition"].keys().difference({"H"})) and
                                    row["Contains aromatic C"] == 0 and  row["Contains amide"] == row["Contains N"] == row["Contains O"],
        'Aliphatic amides': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "N", "O", "H"}) and
                                    {"C", "N", "O"}.issubset(row["Molecular Composition"].keys().difference({"H"}))
                                    and row["Contains aromatic C"] == 0 and row["Contains amide"] != 0 
                                    and row["Contains O"] == 1 and row["Contains N"] == 1,
        'Aliphatic nitriles': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains nitrile"] != 0 and row["Contains N"] == 1,
        'Aliphatic amines (loose)': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] == 0 and 
                                    (row["Contains primary amine"] + row["Contains secondary amine"] + row["Contains tertiary amine"] == row["Contains N"]),
        'Aliphatic amines': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains amine"] != 0 and row["Contains N"] == 1,
        'Aliphatic imines': lambda x, row: row["Organo N (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains imine"] != 0 and row["Contains N"] == 1,

    'Contains Organo P': lambda x, row: {'C','P'}.issubset(row["Molecular Composition"]),
        'Contains C~P': lambda x: count_pattern_occurrences(x, '[#15]~[#6]'),
        'Organo P (strict)': lambda x, row: {'C','P'}.issubset(row["Molecular Composition"].keys()) and
                                set(row["Molecular Composition"].keys()).issubset({"C", "P", "H"}),
        'Contains P':  lambda x, row: row["Molecular Composition"].get("P", 0),
        'Contains P-heterocycle': lambda x: count_pattern_occurrences(x, '[#15&R]'),
        'Contains phosphates': lambda x: count_pattern_occurrences(x, '[#8]=[#15]([#8])([#8])([#8])'),
        'Contains phophonium': lambda x: count_pattern_occurrences(x, '[#15+]'),
        'Organophosphates': '[#8]=[#15]([#8][#6])([#8])([#8])',
        'Organophosphates (alkyl)': lambda x: apply_pattern(Chem.AddHs(x), '[#8]=[#15]([#8]C)([#8][#1,C])([#8][#1,C])'),
        'Organophosphates (aryl)': lambda x: apply_pattern(Chem.AddHs(x), '[#8]=[#15]([#8]c)([#8][#1,c])([#8][#1,c])'),
        'Organophosphites': '[#15&X3]([#8][#6])([#8])([#8])',
        'Organophosphites tautomers': '[#8]=[#15&X4&H1]([#8][#6])([#8])',
        'Organophosphonates': '[#8]=[#15]([#6])([#8])([#8])',
        'Organothiophosphates': '[#16]=[#15]([#8][#6])([#8])([#8])',
    
    'Contains Organo S': lambda x, row: {'C','S'}.issubset(row["Molecular Composition"]),
        'Contains C~S': lambda x: count_pattern_occurrences(x, '[#16]~[#6]'),
        'Organo S (strict)': lambda x, row: {'C','S'}.issubset(row["Molecular Composition"].keys()) and
                                set(row["Molecular Composition"].keys()).issubset({"C", "S", "H"}),
        'Contains S':  lambda x, row: row["Molecular Composition"].get("S", 0),
        'Contains thiol': lambda x: count_pattern_occurrences(x, '[#16&X2&H1]~[#6]'),
        'Contains sulfide': lambda x: count_pattern_occurrences(x, '[#6]-[#16&X2]-[#6]'),
        'Contains thioaldehyde': lambda x: count_pattern_occurrences(Chem.AddHs(x), '[#16&X1]=[#6]([#1])([#6,#1])'),
        'Contains thioketone': lambda x: count_pattern_occurrences(x, '[#16&X1]=[#6]([#6])([#6])'),
        'Contains disulfides': lambda x: count_pattern_occurrences(x, '[#6]-[#16&X2]-[#16&X2]-[#6]'),
        'Contains thioester': lambda x: count_pattern_occurrences(x, '[#6](=[#8])[#16&X2]'),
        'Contains thiosulfinate': lambda x: count_pattern_occurrences(x, '[#16&X3](=[#8])-[#16&X2]'),
        'Contains dithioester': lambda x: count_pattern_occurrences(x, '[#6](=[#16&X1])([#16&X2,#16-])'),
        'Contains dithiocarbamate': lambda x: count_pattern_occurrences(x, '[#7][#6](=[#16&X1])([#16&X2,#16-])'),
        'Contains thioamide': lambda x: count_pattern_occurrences(Chem.AddHs(x), '[#16&X1]=[#6](-[#7])(-[#1,#6])'),
        'Contains thiourea': lambda x: count_pattern_occurrences(x, '[#16&X1]=[#6](-[#7])(-[#7])'),
        'Contains sulfoxide': lambda x: count_pattern_occurrences(x, '[#16&X3](=[#8])([!#8])([!#8])'),
        'Contains sulfinate': lambda x: count_pattern_occurrences(x, '[#16&X3](=[#8])-[#8]'),
        'Contains sulfon': lambda x: count_pattern_occurrences(x, '[#16](=[#8])(=[#8])([!#8])([!#8])'),
        'Contains sulfonyl halide': lambda x: count_pattern_occurrences(x, '[#16](=[#8])(=[#8])([#9,#17,#35,#53])([!#8])'),
        'Contains sulfonamide': lambda x: count_pattern_occurrences(x, '[#16](=[#8])(=[#8])([#7&X3])([!#8])'),
        'Contains thiosulfonate': lambda x: count_pattern_occurrences(x, '[#16](=[#8])(=[#8])([#16&X2])([!#8])'),
        'Contains sulfate': lambda x: count_pattern_occurrences(x, '[#16](=[#8])(=[#8])([#8])([#8])'),
        'Contains sulfonate': lambda x: count_pattern_occurrences(x, '[#16](=[#8])(=[#8])([#8])([#6])'),
        'Contains sulfenes': lambda x: count_pattern_occurrences(x, '[#6]=[#16&X3](=[#8])(=[#8])'),
        'Contains thiocyanate': lambda x: count_pattern_occurrences(x, '[#16-,#16&X2]-[#6]#[#7]'),
        'Contains isothiocyanate': lambda x: count_pattern_occurrences(x, '[#16]=[#6]=[#7&X2]'),
        'Contains S-C-N': lambda x: count_pattern_occurrences(x, '[#16&X1,#16&X2]~[#6&X2]~[#7&X1,#7&X2]'),
        'Contains S-heterocycle': lambda x: count_pattern_occurrences(x, '[#16&R]'),
        'Organosulfate': '[#16](=[#8])([#8][#6])([#8])([#8])',
        'Thiodiproprionates': 'OC(=O)[C&H2][C&H2][S&X2][C&H2][C&H2]-C(=O)O',
        'Thiophenes': '[c&R1]1[c&R1][c&R1][c&R1][s&R1]1',
        'Thiophene derivatives': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#16&R1]~1',
        'Benzothiophenes (1)': '[c&R1]1[c&R1][c&R1][c&R1][c&R2]2[s&R1][c&R1,s&R1][c&R1,s&R1][c&R2]21',
        'Benzothiophenes (2)': '[c&R1]1[c&R1][c&R1][c&R1][c&R2]2[c&R1,s&R1][s&R1][c&R1,s&R1][c&R2]21',
        'Benzothiophenes': lambda x, row: row["Benzothiophenes (1)"] or row["Benzothiophenes (2)"],
        'Benzothiophenes (1) derivatives': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#6&R2]2~[#16&R1]~[#6&R1,#16&R1]~[#6&R1,#16&R1]~[#6&R2]~2~1',
        'Benzothiophenes (2) derivatives': '[#6&R1]1~[#6&R1]~[#6&R1]~[#6&R1]~[#6&R2]2~[#6&R1,#16&R1]~[#16&R1]~[#6&R1,#16&R1]~[#6&R2]~2~1',        
        'Benzothiophene derivatives': lambda x, row: row["Benzothiophenes (1) derivatives"] or row["Benzothiophenes (2) derivatives"],
        'Thiazoles': 'c1scnc1',
        'Thiazoles derivatives': '[#6]1[#16][#6][#7][#6]1',
        'Isothiazoles': 'c1sncc1',
        'Isothiazoles derivatives': '[#6]1[#16][#7][#6][#6]1',
        'Benzothiazoles': 'c1ccc2ncsc2c1',
        'Benzothiazoles derivatives': '[#6]1~[#6]~[#6]~[#6]2~[#7]~[#6]~[#16]~[#6]~2~[#6]~1',
        'Benzisothiazoles': 'c1ccc2sncc2c1',
        'Benzisothiazoles derivatives': '[#6]1[#6][#6][#6]2[#16][#7][#6][#6]2[#6]1',
        'Thioxantenes derivatives': 'c1cccc2c1[#6]c3ccccc3[#16]2',
        'Thioxantones derivatives': 'c1cccc2c1[#6](=[#8])c3ccccc3[#16]2',
        'Triarylsulfonium derivatives': 'c1ccccc1[#16](c2ccccc2)(c3ccccc3)',
        'Aromatic sulfonates': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "S", "O", "H"}) and 
                                    row["Contains aromatic C"] != 0 and row["Contains sulfonate"] != 0 and 
                                    row["Contains O"] == 3 and row["Contains S"] == 1,
        'Aromatic sulfons': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "S", "O", "H"}) and 
                                    row["Contains aromatic C"] != 0 and row["Contains sulfon"] != 0 and 
                                    row["Contains O"] == 2 and row["Contains S"] == 1,
        'Aromatic sulfinates': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "S", "O", "H"}) and 
                                    row["Contains aromatic C"] != 0 and row["Contains sulfinate"] != 0 and 
                                    row["Contains O"] == 2 and row["Contains S"] == 1,
        'Aromatic thioaldehydes': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains thioaldehyde"] != 0 and row["Contains S"] == 1,
        'Aromatic thioketones': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains thioketone"] != 0 and row["Contains S"] == 1,
        'Aromatic thiols': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains thiol"] != 0 and row["Contains S"] == 1,
        'Aromatic sulfides': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] != 0 and 
                                    row["Contains sulfide"] != 0 and row["Contains S"] == 1,
        'Aliphatic sulfonates': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "S", "O", "H"}) and 
                                    row["Contains aromatic C"] == 0 and row["Contains sulfonate"] != 0 and 
                                    row["Contains O"] == 3 and row["Contains S"] == 1,
        'Aliphatic sulfons': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "S", "O", "H"}) and 
                                    row["Contains aromatic C"] == 0 and row["Contains sulfon"] != 0 and 
                                    row["Contains O"] == 2 and row["Contains S"] == 1,
        'Aliphatic sulfinates': lambda x, row: set(row["Molecular Composition"].keys()).issubset({"C", "S", "O", "H"}) and 
                                    row["Contains aromatic C"] == 0 and row["Contains sulfinate"] != 0 and 
                                    row["Contains O"] == 2 and row["Contains S"] == 1,
        'Aliphatic thioaldehydes': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains thioaldehyde"] != 0 and row["Contains S"] == 1,
        'Aliphatic thioketones': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains thioketone"] != 0 and row["Contains S"] == 1,
        'Aliphatic thiols': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains thiol"] != 0 and row["Contains S"] == 1,
        'Aliphatic sulfides': lambda x, row: row["Organo S (strict)"] and row["Contains aromatic C"] == 0 and 
                                    row["Contains sulfide"] != 0 and row["Contains S"] == 1,
    
    
    'Contains Organo Se': lambda x, row: {'C','Se'}.issubset(row["Molecular Composition"]),
        'Contains Se':  lambda x, row: row["Molecular Composition"].get("Se", 0),
        'Organo Se (strict)': lambda x, row: {'C','Se'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "Se", "H"}),
        'Contains C~Se': lambda x: count_pattern_occurrences(x, '[#34]~[#6]'),

    'Contains Organo X': lambda x, row: row["Organic C"] and row["Contains halogen"],
        'Organo X (strict)': lambda x, row: row["Organic C"] and row["Contains halogen"] and
                                set(row["Molecular Composition"]).issubset({"C", "H", "F", "Cl", "Br", "I"}),
        'Contains C~X': lambda x: count_pattern_occurrences(x, '[#9,#17,#35,#53]~[#6]'),

    'Contains Organo F': lambda x, row: {'C','F'}.issubset(row["Molecular Composition"]),
        'Organo F (strict)': lambda x, row: {'C','F'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "F", "H"}),
        'Contains C~F': lambda x: count_pattern_occurrences(x, '[#9]~[#6]'),
        'Contains F':  lambda x, row: row["Molecular Composition"].get("F", 0),
        'PFAS': lambda x: apply_pattern(Chem.AddHs(x), '[C](F)(F)([!#1&!Cl&!Br&!I])([!#1&!Cl&!Br&!I])'),
        'Hydrofluorocarbons': lambda x, row: row["Organo F (strict)"] and row["Contains aromatic C"] == 0 and row["Contains C#C"] == 0 and row["Contains C=C"] == 0,

    'Contains Organo Cl': lambda x, row: {'C','Cl'}.issubset(row["Molecular Composition"]),
        'Organo Cl (strict)': lambda x, row: {'C','Cl'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "Cl", "H"}),
        'Contains C~Cl': lambda x: count_pattern_occurrences(x, '[#17]~[#6]'),
        'Contains Cl':  lambda x, row: row["Molecular Composition"].get("Cl", 0),
        'Contains aliphatic C~Cl': lambda x: count_pattern_occurrences(x, 'C~Cl'),
        'Contains aromatic c~Cl': lambda x: count_pattern_occurrences(x, 'c~Cl'),
        'Chlorinated alkynes': lambda x, row: row["Organo Cl (strict)"] and row["Contains aromatic C"] == 0
                                    and row["Contains C#C"] != 0,
        'Chlorinated alkenes': lambda x, row: row["Organo Cl (strict)"] and row["Contains aromatic C"] == 0
                                    and row["Contains C#C"] == 0 and row["Contains C=C"] != 0,
        'Chlorinated alkanes': lambda x, row: row["Organo Cl (strict)"] and row["Contains aromatic C"] == 0 
                                    and row["Contains C#C"] == 0 and row["Contains C=C"] == 0,
        'Chlorobenzenes': lambda x, row: row["Organo Cl (strict)"] and apply_pattern(Chem.AddHs(x), 'c1([Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])1'),
        'Contains vinyl chloride': lambda x: count_pattern_occurrences(x, '[C&H2]=[C&H1]-[#17]'),
        'Contains vinyl chloride moiety': lambda x: count_pattern_occurrences(x, '[C&H2]~[C&H1]-[#17]'),
        'Chlorophenols': lambda x, row: row["Contains aromatic c~Cl"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Cl])1c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1[#8&H1]'),
        'Chlorophenol derivatives': lambda x,row: row["Contains aromatic c~Cl"] > 0 and cx_smarts_query(x, '[#8]c1ccccc1.*Cl |m:7:2.3.4.5.6.7|'),
        'Chloroanisols': lambda x, row: row["Contains aromatic c~Cl"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Cl])1c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1[#8][C&H3]'),
        'Chloroanilin derivatives': lambda x, row: row["Contains aromatic c~Cl"] > 0 and cx_smarts_query(x, '[#6&H2]c1ccccc1.*Cl |m:7:2.3.4.5.6.7|'),
        'PCBs': lambda x, row: row["Contains C~Cl"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Cl])1c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1-c2c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])2'),
        'PCBs derivatives': lambda x: cx_smarts_query(x, 'c1ccccc1-c2ccccc2.*Cl |m:12:0.1.2.3.4.7.8.9.10.11|'),
        'PCDDs': lambda x, row: row["Contains C~Cl"] > 0 and apply_pattern(Chem.AddHs(x), "c12c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1[#8]c3c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c3[#8]2"),
        'PCDDs derivatives': lambda x: cx_smarts_query(x, 'c12ccccc1[#8]c3ccccc3[#8]2.*Cl |m:14:1.2.3.4.8.9.10.11|'),
        'PCDEs': lambda x, row: row["Contains C~Cl"] > 0 and apply_pattern(Chem.AddHs(x), "c([#1,Cl])1c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1[#8]c3c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])3"),
        'PCDEs derivatives': lambda x: cx_smarts_query(x, 'c1ccccc1[#8]c2ccccc2.*Cl |m:13:0.1.2.3.4.8.9.10.11.12|'),
        'PCDFs': lambda x, row: row["Contains C~Cl"] > 0 and apply_pattern(Chem.AddHs(x), "c12c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1[#8]c3c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c32"),
        'PCDFs derivatives': lambda x: cx_smarts_query(x, 'c12ccccc1[#8]c3ccccc32.*Cl |m:13:1.2.3.4.8.9.10.11|'),
        'DDT derivatives': lambda x: apply_pattern(Chem.AddHs(x), "c1([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c1~[#6](~[#6][Cl])~c2c([#1,Cl])c([#1,Cl])c([#1,Cl])c([#1,Cl])c2([#1,Cl])"),
        'Contains perchlorate': lambda x: count_pattern_occurrences(x, '[Cl]([#8])([#8])([#8])[#8]'),
        'Hydrochlorofluorocarbons': lambda x, row: row["Organo X (strict)"] and row["Contains F"]  > 0 and row["Contains Cl"] > 0 and (row["Contains Br"] + row["Contains I"] == 0) and 
                                                    row["Contains aromatic C"] == 0 and row["Contains C#C"] == 0 and row["Contains C=C"] == 0,
        'Chlorofluorocarbons': lambda x, row: row["Organo X (strict)"] and row["Contains F"]  > 0 and row["Contains Cl"] > 0 and (row["Contains H"] + row["Contains Br"] + row["Contains I"] == 0) 
                                                and (row["Contains aromatic C"] == 0) and (row["Contains C#C"] == 0) and (row["Contains C=C"] == 0),
        'Contains Organo Br': lambda x, row: {'C','Br'}.issubset(row["Molecular Composition"]),
        'Organo Br (strict)': lambda x, row: {'C','Br'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "Br", "H"}),
        'Contains C~Br': lambda x: count_pattern_occurrences(x, '[#35]~[#6]'),
        'Contains Br':  lambda x, row: row["Molecular Composition"].get("Br", 0),
        'Contains aliphatic C~Br': lambda x: count_pattern_occurrences(x, 'C~Br'),
        'Contains aromatic c~Br': lambda x: count_pattern_occurrences(x, 'c~Br'),
        'Bromophenols': lambda x, row: row["Contains aromatic c~Br"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Br])1c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c1[#8&H1]'),
        'Bromophenol derivatives': lambda x, row: row["Contains aromatic c~Br"] > 0 and cx_smarts_query(x, '[#8]c1ccccc1.*Br |m:7:2.3.4.5.6.7|'),
        'PBBs': lambda x, row: row["Contains C~Br"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Br])1c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c1-c2c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])2'),
        'PBBs derivatives': lambda x: cx_smarts_query(x, 'c1ccccc1-c2ccccc2.*Br |m:12:0.1.2.3.4.7.8.9.10.11|'),
        'PBDEs': lambda x, row: row["Contains C~Br"] > 0 and apply_pattern(Chem.AddHs(x), "c([#1,Br])1c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c1[#8]c2c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])2"),
        'PBDEs derivatives': lambda x: cx_smarts_query(x, 'c1ccccc1[#8]c2ccccc2.*Br |m:13:0.1.2.3.4.8.9.10.11.12|'),
        'PBDDs': lambda x, row: row["Contains C~Br"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Br])1c([#1,Br])c([#1,Br])c([#1,Br])c2c1Oc3c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c3[#8]2'),
        'PBDDs derivatives': lambda x: apply_pattern(x, 'c1cccc2c1[#8]c3ccccc3[#8]2.*Br |m:14:0.1.2.3.8.9.10.11|'),
        'PBDFs': lambda x, row: row["Contains C~Br"] > 0 and apply_pattern(Chem.AddHs(x), 'c([#1,Br])1c([#1,Br])c([#1,Br])c([#1,Br])c2c1-c3c([#1,Br])c([#1,Br])c([#1,Br])c([#1,Br])c3[#8]2'),
        'PBDFs derivatives': lambda x: apply_pattern(x, 'c1cccc2c1-c3ccccc3[#8]2.*Br |m:13:0.1.2.3.7.8.9.10|'),
    
    'Contains Organo I': lambda x, row: {'C','I'}.issubset(row["Molecular Composition"]),
        'Organo I (strict)': lambda x, row: {'C','I'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "I", "H"}),
        'Contains C~I': lambda x: count_pattern_occurrences(x, '[#53]~[#6]'),
        'Contains I':  lambda x, row: row["Molecular Composition"].get("I", 0),
        'Diaryliodonium derivatives': 'c1ccccc1-[I]-c2ccccc2',

    'Organic mixed': lambda x, row: row["Organic (only nonmetals)"] and not row["Organo O (strict)"] 
                            and not row["Organo N (strict)"] and not row["Organo P (strict)"] 
                            and not row["Organo S (strict)"] and not row["Organo Se (strict)"]
                            and not row["Organo X (strict)"] and not row["Hydrocarbons"],

    'Contains Organo B': lambda x, row: {'C','B'}.issubset(row["Molecular Composition"]),
        'Contains B':  lambda x, row: row["Molecular Composition"].get("B", 0),
        'Organo B (strict)': lambda x, row: {'C','B'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "B", "H"}),
        'Contains C~B': lambda x: count_pattern_occurrences(x, '[#5]~[#6]'),
        
    'Contains Organo Si': lambda x, row: {'C','Si'}.issubset(row["Molecular Composition"]),
        'Organo Si (strict)': lambda x, row: {'C','Si'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "Si", "H"}),
        'Contains C~Si': lambda x: count_pattern_occurrences(x, '[#14]~[#6]'),
        'Contains Si':  lambda x, row: row["Molecular Composition"].get("Si", 0),
        'Contains Si-heterocycle': lambda x: count_pattern_occurrences(x, '[#14&R]'),
        'Organosiloxanes':'[#14]~[#8]~[#14]~[#6]',
        'Siloxanes':'[#14][#8][#14]',
        'Organosilanes': lambda x: apply_pattern(Chem.AddHs(x), '[#6][#14]([#6,#1])([#6,#1])([#6,#1])'),
        'Contains Si~Cl':  lambda x: count_pattern_occurrences(x, '[#14]~[#17]'),
    
    'Contains Organo Sn': lambda x, row: {'C','Sn'}.issubset(row["Molecular Composition"]),
        'Contains Sn':  lambda x, row: row["Molecular Composition"].get("Sn", 0),
        'Organo Sn (strict)': lambda x, row: {'C','Sn'}.issubset(row["Molecular Composition"]) and
                                set(row["Molecular Composition"]).issubset({"C", "Sn", "H"}),
        'Contains C~Sn':lambda x: count_pattern_occurrences(x, '[#50]~[#6]'),

    'Contains Organometallics': lambda x, row: row["Organic C"] and (row["Contains metal"] or row["Contains metalloid"]),
        'Organometallics (strict)': lambda x, row: row["Organic C"] and (row["Contains metal"] or row["Contains metalloid"]) and
                                set(row["Molecular Composition"]).issubset(
                                    {"C", "H"}.union(set(elements["metals"])).union(set(elements["metalloids"]))
                                    ),
        'Organometallic salt': lambda x, row: row["Contains Organometallics"] and row["Salt"],
        'Organometallic salt (strict)': lambda x, row: row["Contains Organometallics"] and row["Salt"] and not (row["Contains Organo Si"] or row["Contains Organo B"] or row["Contains Organo Sn"]),
        'Organometallic other': lambda x, row: row["Contains Organometallics"] and not row["Organometallic salt"] 
                                                and not row['Contains Organo B'] and not row['Contains Organo Si'] 
                                                and not row['Contains Organo Sn'],
}
