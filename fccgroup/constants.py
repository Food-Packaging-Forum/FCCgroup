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
SMILES_CHECK_COLUMN: Final[str] = "SMILES check"
NOT_GROUPABLE_COLUMN: Final[str] = "Not groupable"
STRUCTURAL_MATCHES_COUNT_COLUMN: Final[str] = "structural_matches_count"
FINGERPRINT_ID_COLUMN: Final[str] = "ID"

# Output/grouping columns
OUTPUT_COLUMN: Final[str] = "Chemical groups"
GROUPS_CONCERN_COLUMN: Final[str] = "Groups of concern"
GROUPS_COLUMN: Final[str] = "groups"
FUNCTION_GROUP_COLUMN: Final[str] = "function_group_in"

# Metadata columns
GROUP_KEY_COLUMN: Final[str] = "Group_key"
LIST_ID_COLUMN: Final[str] = "list_id"
LIST_NAME_COLUMN: Final[str] = "list_name"
FILE_COLUMN: Final[str] = "file_location"
READ_PARAMETERS_COLUMN: Final[str] = "read_parameters"
USE_COLUMN: Final[str] = "use"
IS_COMPLEX_COLUMN: Final[str] = "is_complex"

# Regex mapping columns
REGEX_KEYWORD_COLUMN: Final[str] = "keyword"
REGEX_COLUMN_NAME_COLUMN: Final[str] = "column_name"
REGEX_KEYWORD_LOCATION_COLUMN: Final[str] = "keyword_loc"
REGEX_KEYWORD_TYPE_COLUMN: Final[str] = "keyword_type"
REGEX_GROUP_COLUMN: Final[str] = "Group_name"
REGEX_SUPER_GROUP_COLUMN: Final[str] = "Supergroup"
REGEX_PATTERN_TYPE: Final[str] = "regex"

# CAS processing columns
CAS_VALID_COLUMN: Final[str] = "cas_valid"
CAS_VALID_FORMAT_COLUMN: Final[str] = "cas_valid_format"

# ============================================================================
# GROUPS OF CONCERN
# ============================================================================

GROUPS_CONCERN_DICT = {
    'Inorganic compounds': {
            'Elements (loose)': 'Elements (loose)',
            'Inorganic salt (no C)': 'Inorganic salt (no C)',
            'Inorganic other': 'Inorganic other',
            'Contains toxic heavy metal (As, Cd, Cr, Pb, Hg, Ni)': 'Contains toxic heavy metal (As, Cd, Cr, Pb, Hg, Ni)',
            'Contains B': 'Contains B',
            'Contains perchlorate': 'Contains perchlorate',
    },
    'Organometallic compounds': {
            'Organosiloxanes': 'Organosiloxanes',
            'Contains Organo Sn': 'Contains Organo Sn',
            'Other organometallic compounds': 'Other organometallic compounds',
    },
    'Hydrocarbons': {
            'PAH derivatives hydrocarbon': 'PAH derivatives hydrocarbon',
            'Biphenyls/Terphenyls': 'Biphenyls/Terphenyls',
            'Benzoids': 'Benzoids',
            'Alkanes': 'Alkanes',
    },
    'Organooxygen compounds': {
            'Benzofuran derivatives': 'Benzofuran derivatives',
            'Alkyl phenol derivatives': 'Alkyl phenol derivatives',
            'Hindered phenol derivatives': 'Hindered phenol derivatives',
            'Bisphenol derivatives': 'Bisphenol derivatives',
            'Other phenol derivatives': 'Other phenol derivatives',
            'Aromatic ethers or alcohols (loose)': 'Aromatic ethers or alcohols (loose)',
            'Acetophenone derivatives': 'Acetophenone derivatives',
            'Benzophenone derivatives': 'Benzophenone derivatives',
            'Benzoquinone backbone': 'Benzoquinone backbone',
            'Other benzylketone derivatives': 'Other benzylketone derivatives',
            'Parabens derivatives': 'Parabens derivatives',
            'Salicylates derivatives': 'Salicylates derivatives',
            'Ortho-phthalates': 'Ortho-phthalates',
            'Terephthalates': 'Terephthalates',
            'Adipic acid esters': 'Adipic acid esters',
            'Citric acid esters': 'Citric acid esters',
            'Triglyceride': 'Triglyceride',
            'Aliphatic carbonyls (loose)': 'Aliphatic carbonyls (loose)',
    },
    'Organonitrogen compounds': {
            'Imidazole': 'Imidazole',
            'Benzotriazoles derivatives': 'Benzotriazoles derivatives',
            'Contains imine': 'Contains imine',
            'Contains isocyanate': 'Contains isocyanate',
            'Contains azo': 'Contains azo',
            'Contains nitrile': 'Contains nitrile',
            'Contains nitro': 'Contains nitro',
            'Contains nitrosamine': 'Contains nitrosamine',
            'Primary aromatic amines': 'Primary aromatic amines',
            'Secondary aromatic amines': 'Secondary aromatic amines',
            'Tertiary aromatic amines': 'Tertiary aromatic amines',
            'Aliphatic amines (loose)': 'Aliphatic amines (loose)',
            'Aliphatic amides (loose)': 'Aliphatic amides (loose)',
    },
    'Organophosphorus compounds': {
            'Organophosphates': 'Organophosphates',
            'Organophosphites': 'Organophosphites',
    },    
    'Organosulfur compounds': {
            'Benzothiazoles derivatives': 'Benzothiazoles derivatives',
            'Contains thiol': 'Contains thiol',
            'Contains sulfate': 'Contains sulfate',
            'Contains sulfinate': 'Contains sulfinate',
            'Contains dithiocarbamate': 'Contains dithiocarbamate',
            'Contains S-C-N': 'Contains S-C-N',
    },
    'Organohalogen compounds': {
            'PFAS': 'PFAS',
            'Chlorinated alkanes': 'Chlorinated alkanes',
            'Chlorinated alkenes': 'Chlorinated alkenes',
            'PCBs': 'PCBs',
            'Other organochlorine compounds': 'Other organochlorine compounds',
            'PBDEs': 'PBDEs',
            'PBBs': 'PBBs',
            'Other organobromine compounds': 'Other organobromine compounds',
            'Contains Organo I': 'Contains Organo I',
        },
}


GROUPS_CONCERN = []
RENAME_GROUPS_CONCERN = {}
for group_dict in GROUPS_CONCERN_DICT.values():
    for group, renamed_group in group_dict.items():
        GROUPS_CONCERN.append(group)
        RENAME_GROUPS_CONCERN[group] = renamed_group

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
FCC_UNIVERSE_FILE: Final[Path] = Path(__file__).parent / ".." / "tests" / "grouped_chemicals.xlsx"


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


# ============================================================================
# MULTIINDEX COLUMN GROUP LABELS
# ============================================================================

MULTIINDEX_IDENTIFIER_LABEL: Final[str] = "Identifier"
MULTIINDEX_STRUCTURAL_LABEL: Final[str] = "Structural patterns"
MULTIINDEX_LISTS_LABEL: Final[str] = "Lists"
MULTIINDEX_REGEX_LABEL: Final[str] = "Regex"


# ============================================================================
# COMPTOX CONFIGURATION
# ============================================================================

COMPTOX_SEARCH_EQUAL_URL: Final[str] = "https://comptox.epa.gov/ctx-api/chemical/search/equal/"
COMPTOX_DTXCID_DETAIL_URL: Final[str] = "https://comptox.epa.gov/ctx-api/chemical/detail/search/by-dtxcid/"
COMPTOX_API_KEY_ENV: Final[str] = "COMPTOX_API_KEY"
