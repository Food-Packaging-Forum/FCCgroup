"""
Data processing and utilities subpackage.

Provides utilities for:
- CAS Registry Number validation and processing
- Chemical list and mapping file loading
- DataFrame manipulation and regex pattern application
"""

from .cas import (
    harmonize_cas_columns,
    unify_cas_columns,
    expand_cas_column,
    is_valid_cas,
    extract_and_reformat_date,
    clean_casrn,
    clean_cas_g13
)

from .lists import (
    load_mapping_file,
    load_lists,
    harmonize_function_columns,
)

from .frame import (
    apply_simple_regex,
    process_groups,
    combine_groups,
    generate_parent_groups,
    rename_columns_to_multiindex,
    column_contains_patterns,
    merge_dataframes,
)

__all__ = [
    # CAS processing
    "unify_cas_columns",
    "expand_cas_column",
    "is_valid_cas",
    "extract_and_reformat_date",
    "clean_casrn",
    "clean_cas_g13",
    # List loading
    "load_mapping_file",
    "load_lists",
    # DataFrame utilities
    "apply_simple_regex",
    "process_groups",
    "combine_groups",
    "generate_parent_groups",
    "rename_columns_to_multiindex",
    "column_contains_patterns",
    "merge_dataframes",
]
