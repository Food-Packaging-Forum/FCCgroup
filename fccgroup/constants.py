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
IMPORTANT_GROUPS_COLUMN: Final[str] = "Important groups"
PRIORITY_GROUPS_COLUMN: Final[str] = "Priority groups"
GROUPS_COLUMN: Final[str] = "groups"
FUNCTION_GROUP_COLUMN: Final[str] = "function_group_in"

# Metadata columns
GROUP_KEY_COLUMN: Final[str] = "Group_key"
LIST_ID_COLUMN: Final[str] = "list_id"
LIST_NAME_COLUMN: Final[str] = "list_name"
FILE_COLUMN: Final[str] = "list_file_location"
READ_PARAMETERS_COLUMN: Final[str] = "read_parameters"
USE_COLUMN: Final[str] = "use"
IS_COMPLEX_COLUMN: Final[str] = "is_complex"

# Regex mapping columns
REGEX_KEYWORD_COLUMN: Final[str] = "keyword"
REGEX_COLUMN_NAME_COLUMN: Final[str] = "Column_name"
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
            'Elements (loose)': {'rename_name': 'Elements (loose)', 'priority_group': True},
            'Inorganic salt (no C)': {'rename_name': 'Inorganic salt (no C)', 'priority_group': True},
            'Inorganic other': {'rename_name': 'Inorganic other', 'priority_group': True},
            'Contains toxic heavy metal (As, Cd, Cr, Pb, Hg, Ni)': {'rename_name': 'Contains toxic heavy metal (As, Cd, Cr, Pb, Hg, Ni)', 'priority_group': True},
            'Contains B': {'rename_name': 'Contains B', 'priority_group': True},
            'Contains perchlorate': {'rename_name': 'Contains perchlorate', 'priority_group': True},
    },
    'Organometallic compounds': {
            'Organosiloxanes': {'rename_name': 'Organosiloxanes', 'priority_group': False},
            'Contains Organo Sn': {'rename_name': 'Contains Organo Sn', 'priority_group': True},
            'Other organometallic compounds': {'rename_name': 'Other organometallic compounds', 'priority_group': False},
    },
    'Hydrocarbons': {
            'PAH derivatives hydrocarbon': {'rename_name': 'PAH derivatives hydrocarbon', 'priority_group': True},
            'Biphenyls/Terphenyls': {'rename_name': 'Biphenyls/Terphenyls', 'priority_group': True},
            'Benzoids': {'rename_name': 'Benzoids', 'priority_group': True},
            'Alkanes': {'rename_name': 'Alkanes', 'priority_group': False},
    },
    'Organooxygen compounds': {
            'Benzofuran derivatives': {'rename_name': 'Benzofuran derivatives', 'priority_group': True},
            'Alkyl phenol derivatives': {'rename_name': 'Alkyl phenol derivatives', 'priority_group': True},
            'Hindered phenol derivatives': {'rename_name': 'Hindered phenol derivatives', 'priority_group': False},
            'Bisphenol derivatives': {'rename_name': 'Bisphenol derivatives', 'priority_group': False},
            'Other phenol derivatives': {'rename_name': 'Other phenol derivatives', 'priority_group': True},
            'Aromatic ethers or alcohols (loose)': {'rename_name': 'Aromatic ethers or alcohols (loose)', 'priority_group': True},
            'Acetophenone derivatives': {'rename_name': 'Acetophenone derivatives', 'priority_group': True},
            'Benzophenone derivatives': {'rename_name': 'Benzophenone derivatives', 'priority_group': False},
            'Benzoquinone backbone': {'rename_name': 'Benzoquinone backbone', 'priority_group': False},
            'Other benzylketone derivatives': {'rename_name': 'Other benzylketone derivatives', 'priority_group': False},
            'Parabens derivatives': {'rename_name': 'Parabens derivatives', 'priority_group': False},
            'Salicylates derivatives': {'rename_name': 'Salicylates derivatives', 'priority_group': False},
            'Ortho-phthalates': {'rename_name': 'Ortho-phthalates', 'priority_group': True},
            'Terephthalates': {'rename_name': 'Terephthalates', 'priority_group': False},
            'Adipic acid esters': {'rename_name': 'Adipic acid esters', 'priority_group': False},
            'Citric acid esters': {'rename_name': 'Citric acid esters', 'priority_group': False},
            'Triglyceride': {'rename_name': 'Triglyceride', 'priority_group': False},
            'Aliphatic carbonyls (loose)': {'rename_name': 'Aliphatic carbonyls (loose)', 'priority_group': False},
    },
    'Organonitrogen compounds': {
            'Imidazole': {'rename_name': 'Imidazole', 'priority_group': True},
            'Benzotriazoles derivatives': {'rename_name': 'Benzotriazoles derivatives', 'priority_group': True},
            'Contains imine': {'rename_name': 'Contains imine', 'priority_group': True},
            'Contains isocyanate': {'rename_name': 'Contains isocyanate', 'priority_group': True},
            'Contains azo': {'rename_name': 'Contains azo', 'priority_group': False},
            'Contains nitrile': {'rename_name': 'Contains nitrile', 'priority_group': True},
            'Contains nitro': {'rename_name': 'Contains nitro', 'priority_group': True},
            'Contains nitrosamine': {'rename_name': 'Contains nitrosamine', 'priority_group': True},
            'Primary aromatic amines': {'rename_name': 'Primary aromatic amines', 'priority_group': True},
            'Secondary aromatic amines': {'rename_name': 'Secondary aromatic amines', 'priority_group': False},
            'Tertiary aromatic amines': {'rename_name': 'Tertiary aromatic amines', 'priority_group': False},
            'Aliphatic amines (loose)': {'rename_name': 'Aliphatic amines (loose)', 'priority_group': True},
            'Aliphatic amides (loose)': {'rename_name': 'Aliphatic amides (loose)', 'priority_group': True},
    },
    'Organophosphorus compounds': {
            'Organophosphates': {'rename_name': 'Organophosphates', 'priority_group': True},
            'Organophosphites': {'rename_name': 'Organophosphites', 'priority_group': False},
    },    
    'Organosulfur compounds': {
            'Benzothiazoles derivatives': {'rename_name': 'Benzothiazoles derivatives', 'priority_group': False},
            'Contains thiol': {'rename_name': 'Contains thiol', 'priority_group': True},
            'Contains sulfate': {'rename_name': 'Contains sulfate', 'priority_group': True},
            'Contains sulfinate': {'rename_name': 'Contains sulfinate', 'priority_group': True},
            'Contains dithiocarbamate': {'rename_name': 'Contains dithiocarbamate', 'priority_group': True},
            'Contains S-C-N': {'rename_name': 'Contains S-C-N', 'priority_group': True},
    },
    'Organohalogen compounds': {
            'PFAS': {'rename_name': 'PFAS', 'priority_group': True},
            'Chlorinated alkanes': {'rename_name': 'Chlorinated alkanes', 'priority_group': True},
            'Chlorinated alkenes': {'rename_name': 'Chlorinated alkenes', 'priority_group': True},
            'PCBs': {'rename_name': 'PCBs', 'priority_group': False},
            'Other organochlorine compounds': {'rename_name': 'Other organochlorine compounds', 'priority_group': True},
            'PBDEs': {'rename_name': 'PBDEs', 'priority_group': False},
            'PBBs': {'rename_name': 'PBBs', 'priority_group': False},
            'Other organobromine compounds': {'rename_name': 'Other organobromine compounds', 'priority_group': True},
            'Contains Organo I': {'rename_name': 'Contains Organo I', 'priority_group': True},
        },
}


IMPORTANT_GROUPS = []
RENAME_GROUPS_CONCERN = {}
PRIORITY_GROUPS_CONCERN = []
for group_dict in GROUPS_CONCERN_DICT.values():
    for group, info in group_dict.items():
        IMPORTANT_GROUPS.append(group)
        RENAME_GROUPS_CONCERN[group] = info['rename_name']
        if info['priority_group']:
            PRIORITY_GROUPS_CONCERN.append(group)

# ============================================================================
# MAPPING FILE CONFIGURATION
# ============================================================================

MAPPING_FILE_NAME: Final[str] = "Mapping.xlsx"
MAPPING_SHEET_LISTS: Final[str] = "B-Lists"
MAPPING_SHEET_KEYWORDS: Final[str] = "B-Keywords"
MAPPING_SHEET_SMARTS: Final[str] = "B-SMARTS"

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
