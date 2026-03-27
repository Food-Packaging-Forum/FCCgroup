"""Shared runtime constants for FCCgroup.

This module centralizes column names, asset paths, mapping sheet names, and
list identifiers so the rest of the package can avoid duplicated literals.
"""

from pathlib import Path
from typing import Final


# ============================================================================
# CANONICAL CONFIGURATION KEYS
# ============================================================================

CANONICAL_CAS_KEY: Final[str] = "cas"
CANONICAL_SMILES_KEY: Final[str] = "smiles"
CANONICAL_NAME_COLUMNS_KEY: Final[str] = "name_columns"
CANONICAL_FORMULA_KEY: Final[str] = "formula"


# ============================================================================
# COLUMN NAMES
# ============================================================================

# Chemical data columns
CAS_COLUMN: Final[str] = "casId"
SMILES_COLUMN: Final[str] = "SMILES"
COMMON_NAME_COLUMN: Final[str] = "commonName"
IUPAC_NAME_COLUMN: Final[str] = "IUPAC_name"
FORMULA_COLUMN: Final[str] = "formula"

# Derived/internal columns
ENRICHED_NAME_COLUMNS_COLUMN: Final[str] = "column_names"
COMBINED_NAME_COLUMN: Final[str] = "name"
MOLECULAR_FORMULA_COLUMN: Final[str] = "MolFormula"
SMILES_CHECK_COLUMN: Final[str] = "SMILES_check"
STRUCTURAL_MATCHES_COUNT_COLUMN: Final[str] = "structural_matches_count"
FINGERPRINT_ID_COLUMN: Final[str] = "ID"

# Output/grouping columns
OUTPUT_COLUMN: Final[str] = "Chemical groups"
GROUPS_COLUMN: Final[str] = "groups"
FUNCTION_GROUP_COLUMN: Final[str] = "function_group_in"

# Metadata columns
GROUP_KEY_COLUMN: Final[str] = "Group_key"
LIST_ID_COLUMN: Final[str] = "list_id"
LIST_NAME_COLUMN: Final[str] = "list_name"
FILE_COLUMN: Final[str] = "file"
READ_PARAMETERS_COLUMN: Final[str] = "read_parameters"
USE_COLUMN: Final[str] = "use"
IS_COMPLEX_COLUMN: Final[str] = "is_complex"

# Regex mapping columns
REGEX_KEYWORD_COLUMN: Final[str] = "keyword"
REGEX_COLUMN_NAME_COLUMN: Final[str] = "column_name"
REGEX_KEYWORD_LOCATION_COLUMN: Final[str] = "keyword_loc"
REGEX_KEYWORD_TYPE_COLUMN: Final[str] = "keyword_type"
REGEX_GROUP_COLUMN: Final[str] = "group"
REGEX_SUPER_GROUP_COLUMN: Final[str] = "super_group"
REGEX_PATTERN_TYPE: Final[str] = "regex"

# CAS processing columns
CAS_VALID_COLUMN: Final[str] = "cas_valid"
CAS_VALID_FORMAT_COLUMN: Final[str] = "cas_valid_format"


# ============================================================================
# MAPPING FILE CONFIGURATION
# ============================================================================

MAPPING_FILE_NAME: Final[str] = "Mapping.xlsx"
MAPPING_SHEET_LISTS: Final[str] = "B - Lists"
MAPPING_SHEET_KEYWORDS: Final[str] = "B - Keywords"
MAPPING_SHEET_SMARTS: Final[str] = "B - SMARTS"

# Keyword to identify complex lists (containing sublists with groups and functions)
COMPLEX_LIST_KEYWORD: Final[str] = "*** Several ***"


# ============================================================================
# CAS NUMBER CONFIGURATION
# ============================================================================

# Possible CAS column names to check (case-insensitive matching)
CAS_COLUMN_ALIASES: Final[list[str]] = [
    "casid",
    "cas number",
    "casrn",
    "cas no",
    "cas no.",
    "casnumber",
    "cas rn",
    "cas\nnumber",
    "cas reg no (or other id)",
    "comptox_casrn",
    "cas_fixed",
    "cas",
    "cas_number",
    "casregistry",
    "cas_rn",
    "casrn_text",
]

# CAS number format validation regex (e.g., "123-45-67")
CAS_FORMAT_PATTERN: Final[str] = r"^\d{2,}-\d{2}-\d{1}$"
NO_CAS_TOKEN: Final[str] = "NOCAS"


# ============================================================================
# ASSET PATHS
# ============================================================================

# Asset root directory inside the installed package
ASSETS_DIR: Final[Path] = Path(__file__).parent / "assets"

# Asset subdirectories
LISTS_DIR: Final[Path] = ASSETS_DIR / "lists"
SMILES_LOOKUP_FILE: Final[Path] = ASSETS_DIR / "smiles_lookup.tsv"
FCC_UNIVERSE_FILE: Final[Path] = ASSETS_DIR / "FCCuniverse_grouping_in.xlsx"


# ============================================================================
# LIST IDENTIFIERS
# ============================================================================

LIST_FRAGRANCES_G05: Final[str] = "Fragrances_G05"
LIST_G04: Final[str] = "G04"
LIST_G13: Final[str] = "G13"
LIST_G23: Final[str] = "G23"
LIST_G24: Final[str] = "G24"
LIST_G25: Final[str] = "G25"
LIST_PESTICIDES_G17: Final[str] = "Pesticides_G17"
LIST_METABOLITES_G22: Final[str] = "Metabolites_G22"
LIST_METABOLITES_G28_01: Final[str] = "Metabolites_G28_01"
LIST_FOOD_RELATED_G28_02: Final[str] = "Food related compounds_G28_02"

PATTERN_COMPLEX_LIST_IDS: Final[tuple[str, ...]] = (LIST_G04, LIST_G13, LIST_G25)

SPECIAL_CAS_RENAME_COLUMNS: Final[dict[str, str]] = {
    LIST_FRAGRANCES_G05: "CAS No.",
    LIST_G13: "CASRN_text",
    LIST_G23: "cas_fixed",
    LIST_METABOLITES_G28_01: "CASRN",
    LIST_FOOD_RELATED_G28_02: "CASRN",
}

EXPANDABLE_CAS_LIST_SEPARATORS: Final[dict[str, list[str]]] = {
    LIST_FRAGRANCES_G05: [";"],
    LIST_G13: [",", "/", r"\(", "and", "&", "or"],
}

LIST_FUNCTION_COLUMNS: Final[dict[str, list[str]]] = {
    LIST_G04: ["Group name"],
    LIST_G13: ["Product Type"],
    LIST_G23: ["Harmonized_functions"],
    LIST_G24: ["function_name"],
    LIST_G25: ["function", "Other_function_in_plastic"],
}

PLASTICMAP_SHEET_CHEMICAL_FUNCTION: Final[str] = "ChemicalFunction"
PLASTICMAP_SHEET_CHEMICAL: Final[str] = "Chemical"
PLASTICMAP_SHEET_FUNCTION: Final[str] = "Function"
PLASTICMAP_SHEET_SOURCE: Final[str] = "Source"

PLASTCHEM_COLUMNS_KEY: Final[str] = "PlastChem_columns"
FUNCTION_SYNONYMS_G4_KEY: Final[str] = "function_synonyms_g4"
FUNCTION_SYNONYMS_REGEX_KEY: Final[str] = "function_synonyms_regex"
