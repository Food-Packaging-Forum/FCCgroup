"""
CAS Registry Number (CASRN) processing and validation utilities.

This module provides functions for standardizing, validating, and processing
CAS Registry Numbers across multiple DataFrames.
"""

import re
from typing import Dict, List, Optional

import pandas as pd

from ..constants import *


def harmonize_cas_columns(all_lists: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Harmonize CAS column names across all functional lists.
    
    Different chemical list files use different column naming conventions
    for CAS Registry Numbers (casid, CASRN, cas_number, etc.). This method
    standardizes all to the canonical 'casId' column name for consistent
    merging and lookups.

    Also handles:
    - Multi-CAS entries (splitting by separators like ";", "/", etc.)
    - CAS number format cleaning (removing invalid characters)
    - Validation using CAS checksum
    
    Args:
        all_lists: Dictionary mapping list names to their corresponding DataFrames
    Returns:
        Dictionary with standardized CAS column names and cleaned CAS numbers
    """
    all_lists = unify_cas_columns(lists=all_lists, cas_column_names=CAS_COLUMN_ALIASES)

    for list_name, col_name in SPECIAL_CAS_RENAME_COLUMNS.items():
        if list_name in all_lists:
            current_list = all_lists[list_name].rename(
                columns={col_name: CAS_COLUMN}
            )
            current_list[CAS_COLUMN] = current_list[CAS_COLUMN].fillna("").astype(str)
            all_lists[list_name] = current_list
    
    # Expand multi-CAS entries
    for curr_list, separators in EXPANDABLE_CAS_LIST_SEPARATORS.items():
        if curr_list in all_lists:
            curr_list_df = all_lists[curr_list]
            curr_list_df = expand_cas_column(
                curr_list_df,
                cas_column=CAS_COLUMN,
                separators=separators,
                verbose=False
            )
            all_lists[curr_list] = curr_list_df
    
    # Clean CAS numbers
    remove_leading = '0(['
    remove_leadtrail = '\'"^)]-'
    replace_dict = {
        'â€\x90': '-', 'â€Ž': '', 'â€"': '-', '┬á': '', '\xa0': '',
        'Â\xad': '', 'Â\x81': '', '_x0018_': '', '_x0004_': '', '_x0001_': '',
        ' ': '', '.': '-', '--': '-'
    }
    
    all_lists = clean_casrn(
        all_lists,
        cas_column=CAS_COLUMN,
        leading=remove_leading,
        lead_trailing=remove_leadtrail,
        replace_dict=replace_dict
    )
    
    # G13 specific cleaning
    if LIST_G13 in all_lists:
        all_lists[LIST_G13][CAS_COLUMN] = all_lists[LIST_G13][CAS_COLUMN].apply(
            clean_cas_g13
        )
    return all_lists

def unify_cas_columns(
    lists: Dict[str, pd.DataFrame],
    cas_column_names: List[str]
) -> Dict[str, pd.DataFrame]:
    """
    Standardizes CASRN columns across multiple DataFrames.
    
    For each DataFrame in the provided dictionary:
    - Checks for the presence of CASRN columns as specified in `cas_column_names`
    - Prints a message if no CASRN column is found
    - Prints a warning if multiple CASRN columns are found
    - Renames the matching CASRN column to 'casId'
    - Fills missing values with empty string and casts to string type
    
    Args:
        lists: Dictionary mapping list names to pandas DataFrames to process
        cas_column_names: List of possible column names that may contain CASRNs
    
    Returns:
        Dictionary with standardized CASRN columns named 'casId'
        
    Raises:
        ValueError: If lists dict is empty or cas_column_names list is empty
        
    Note:
        Returns a new dictionary; input dict is not modified in place
    """
    # Validation
    if not lists or len(lists) == 0:
        raise ValueError("lists dictionary cannot be empty")
    
    if not cas_column_names or len(cas_column_names) == 0:
        raise ValueError("cas_column_names list cannot be empty")
    
    new_lists: Dict[str, pd.DataFrame] = {}
    
    for curr_list, curr_list_df in lists.items():
        df_copy = curr_list_df.copy()

        # Identify if CASRNs columns are present
        # Convert columns to strings first to handle integer column names
        columns_lower = [str(col).lower() for col in df_copy.columns]
        overlap = set(cas_column_names) & set(columns_lower)
        overlap = [col for col in df_copy.columns if str(col).lower() in overlap]  # Preserve case-sensitivity
        if len(overlap) == 0: 
            print(f"\n{curr_list}: FALSE \n{list(df_copy.columns)}")
        
        if len(overlap) > 1: 
            print(f"\n{curr_list}: WARNING: More than one CASRN column found \n{overlap}")

        # Rename CASRN columns to the canonical CAS column
        if len(overlap) == 1:
            df_copy.rename(columns={overlap[0]: CAS_COLUMN}, inplace=True)
            df_copy[CAS_COLUMN] = df_copy[CAS_COLUMN].fillna('').astype(str)
        
        new_lists[curr_list] = df_copy
    
    return new_lists


def expand_cas_column(
    df: pd.DataFrame,
    cas_column: str,
    separators: Optional[List[str]] = None,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Splits rows with multiple CASRNs in the specified column into separate rows.
    
    Args:
        df: Input DataFrame
        cas_column: Name of the column containing CASRNs
        separators: List of separators used for splitting multiple CASRNs. Defaults to [';']
        verbose: If True, prints information about the operation
    
    Returns:
        DataFrame with expanded rows where each row contains a single CASRN
        
    Raises:
        ValueError: If df is empty or cas_column doesn't exist
        
    Note:
        Removes duplicate rows after expansion
    """
    # Validation
    if df is None or df.empty:
        raise ValueError("df DataFrame cannot be empty")
    
    if cas_column not in df.columns:
        raise ValueError(f"Column '{cas_column}' not found in DataFrame. Available columns: {list(df.columns)}")
    
    if separators is None:
        separators = [';']
    
    df_copy = df.copy()
    
    for separator in separators:
        df_copy[cas_column] = df_copy[cas_column].astype(str)
        number_separators = df_copy[cas_column].str.count(separator).sum()
        if verbose:
            print(f'Current length ({len(df_copy)}) + Number of separators ({number_separators}) = Final table ({len(df_copy)+number_separators})')

        if number_separators == 0:
            if verbose:
                print("No CASRNs to expand, returning original DataFrame.")
            continue
    
        # create expanded DataFrame with separate rows for each CASRN
        df_copy[cas_column] = df_copy[cas_column].str.split(separator)
        df_copy = df_copy.explode(cas_column).reset_index(drop=True)
    
    df_copy = df_copy.drop_duplicates().reset_index(drop=True)
    return df_copy


def is_valid_cas(cas: str, is_format: bool = False) -> bool:
    """
    Check if CASRN follows standard format and has a valid check digit.
    
    CASRN format: 123-45-67
    
    Check digit calculation:
        1. Remove hyphens: 123456 7
        2. Split CASRN into its digits: 1, 2, 3, 4, 5, 6, 7
        3. Last digit is the check digit: 7
        4. Other digits (1,2,3,4,5,6) are ordered based on their position (1-5 from right to left)
        5. Multiply each digit with its position and sum: 6*1 + 5*2 + 4*3 + 3*4 + 2*5 + 1*6 = ...
        6. Check digit is the last digit of the sum modulo 10
        7. Return True if last digit of the CASRN equals check digit
    
    Args:
        cas: CASRN string to validate
        is_format: If True, only check format (skip checksum validation)
    
    Returns:
        True if CASRN is valid (format and checksum if applicable), False otherwise
        
    Raises:
        TypeError: If cas is not a string
    """
    # Validation
    if not isinstance(cas, str):
        raise TypeError(f"cas must be a string, got {type(cas).__name__}")
    
    # Check if CASRN follows standard format
    if not re.search(CAS_FORMAT_PATTERN, cas): 
        return False
    
    if is_format:
        return True
    
    # Split CASRN into its digits
    cas_clean = cas.replace("-", "")
    check_digit = int(cas_clean[-1])
    digits = [int(d) for d in cas_clean[:-1]][::-1]
    positions = list(range(1, len(cas_clean)))

    # Multiply each digit with its position and compare to check digit
    sum_digit = sum(d * p for d, p in zip(digits, positions))
    return (sum_digit % 10 == check_digit)


def extract_and_reformat_date(casrn_date: str) -> str:
    """
    Check if text appears to be a date and convert it to CASRN format.
    
    Specifically:
        1. Searches for pattern "\\d{2,}-\\d{2}-\\d{2}" (2+ digits, hyphen, 2 digits, hyphen, 2 digits)
        2. If found, splits by '-' into year, month, day
        3. Converts last part to integer (removes leading zeros)
        4. Returns reformatted string or original text if no date pattern found
    
    Args:
        casrn_date: String that may be a date-format CASRN
    
    Returns:
        Original string if no date pattern found, otherwise reformatted as year-month-day
    """
    match = re.search(r'\d{2,}\-\d{2}\-\d{2}', casrn_date)
    if match:
        year, month, day = match.group().split('-')
        return f"{year}-{month}-{int(day)}"
    return casrn_date


def clean_casrn(
    all_lists: Dict[str, pd.DataFrame],
    cas_column: str,
    leading: str,
    lead_trailing: str,
    replace_dict: Dict[str, str],
    verbose: bool = False
) -> Dict[str, pd.DataFrame]:
    """
    Clean and validate CAS Registry Numbers in multiple DataFrames.
    
    For each DataFrame, this function:
        1. Removes specified leading and trailing characters
        2. Replaces special characters based on replace_dict
        3. Converts date-formatted CASRNs to standard CASRN format
        4. Validates each CASRN
        5. Removes duplicate rows
    
    Adds two new columns to each DataFrame:
        - 'cas_valid': Boolean indicating if CASRN is valid (format + checksum)
        - 'cas_valid_format': Boolean indicating if CASRN has valid format
    
    Args:
        all_lists: Dictionary mapping list names to DataFrames with CASRNs
        cas_column: Name of the column containing CAS values
        leading: Characters to remove from start of CASRN
        lead_trailing: Characters to remove from beginning and end of CASRN
        replace_dict: Dictionary mapping characters to their replacements
        verbose: If True, prints invalid CASRNs found
    
    Returns:
        Dictionary with cleaned DataFrames
        
    Raises:
        ValueError: If all_lists is empty
        Exception: If error occurs during processing for a specific list
        
    Performance Notes:
        - String operations: ~1-10ms per operation for 10K rows (vectorized)
        - .str methods: Use pandas C backend (highly optimized)
        - .apply() iteration: ~5-20ms per 10K rows
        - Groupby operation: ~5-20ms for deduplication
        - Total for 10K rows: ~50-100ms
    """
    # Validation
    if not all_lists or len(all_lists) == 0:
        raise ValueError("all_lists dictionary cannot be empty")
    
    new_lists: Dict[str, pd.DataFrame] = {}
    
    for curr_list, df in all_lists.items():
        try:
            df_copy = df.copy()
            df_copy[cas_column] = df_copy[cas_column].fillna('').astype(str)
            
            # Clean the CASRN column (vectorized string operations - fast)
            df_copy[cas_column] = df_copy[cas_column].str.strip()
            df_copy[cas_column] = df_copy[cas_column].str.strip(lead_trailing)
            df_copy[cas_column] = df_copy[cas_column].str.lstrip(leading)
            df_copy[cas_column] = df_copy[cas_column].replace(replace_dict, regex=False)
            
            # Convert date-formatted CASRNs (row iteration - slow, but small dataset typically)
            df_copy[cas_column] = df_copy[cas_column].apply(extract_and_reformat_date)
            
            # Validate CASRNs (vectorized - fast)
            df_copy[CAS_VALID_COLUMN] = df_copy[cas_column].apply(is_valid_cas)
            df_copy[CAS_VALID_FORMAT_COLUMN] = df_copy[cas_column].apply(lambda x: is_valid_cas(x, True))
            
            # Remove duplicates (groupby operation)
            df_copy = df_copy.groupby(cas_column, as_index=False).first()

            new_lists[curr_list] = df_copy
            
            # Report invalid CASRNs if verbose
            if verbose:
                curr_list_df_not_valid = df_copy.loc[
                    (~df_copy[CAS_VALID_COLUMN]) & ~(df_copy[cas_column].str.contains(NO_CAS_TOKEN))
                ]
                if len(curr_list_df_not_valid) > 0:
                    not_valid_cas = list(curr_list_df_not_valid[cas_column].unique())
                    print(f'{curr_list}: {len(not_valid_cas)} invalid CASRNs \n {not_valid_cas}')
        
        except Exception as e:
            raise Exception(
                f"Error processing list '{curr_list}' with columns {list(df.columns)}: {str(e)}"
            )
    
    return new_lists


def clean_cas_g13(cas_number: str) -> str:
    """
    Clean CAS numbers for G13 list (fragrance ingredients).
    
    Handles special formatting issues specific to the G13 fragrance list.
    
    Args:
        cas_number: CAS number string potentially needing correction
        
    Returns:
        Cleaned CAS number in standard format (XXX-XX-X)
    """
    if not re.search(r'^[\d\-]+$', cas_number):
        return cas_number
    if is_valid_cas(cas_number):
        return cas_number
    
    cas_nohyphen = cas_number.replace("-", "")
    cas_digits = [str(d) for d in cas_nohyphen]
    
    if re.search(r'\d{2}+\-\d{2}(\d{2})$', cas_number):
        if cas_digits[-2] == '0':
            cas_digits.pop(-2)
    
    part3 = cas_digits[-1]
    part2 = ''.join(cas_digits[-3:-1])
    part1 = ''.join(cas_digits[:-3])
    return f"{part1}-{part2}-{part3}"