"""
DataFrame manipulation and transformation utilities.

This module provides utilities for:
- Applying regex and substring patterns to DataFrames
- Processing and combining group columns
- Generating parent-child group relationships
- Creating MultiIndex column structures
"""

from typing import Dict, List, Union

import pandas as pd
import numpy as np
import re

from ..constants import *

def apply_simple_regex(df: pd.DataFrame, regex_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply regex or substring matching to DataFrame columns based on patterns.
    
    Iterates through pattern definitions and applies them to specified columns,
    creating new boolean columns for matches.
    
    Args:
        df: Input DataFrame containing data to be checked
        regex_df: DataFrame with pattern definitions containing columns:
            - 'keyword': Pattern or substring to search for
            - 'column_name': Name for the resulting boolean column
            - 'keyword_loc': Column in df to apply the pattern to
            - 'keyword_type': "regex" for regex pattern, any other value for substring
    
    Returns:
        Original DataFrame concatenated with new boolean columns for each pattern
        
    Raises:
        KeyError: If required columns are missing from regex_df or df
    """
    # Validate required columns
    required_cols = [
        REGEX_KEYWORD_COLUMN,
        REGEX_COLUMN_NAME_COLUMN,
        REGEX_KEYWORD_LOCATION_COLUMN,
        REGEX_KEYWORD_TYPE_COLUMN,
    ]
    if not all(col in regex_df.columns for col in required_cols):
        raise KeyError(f"regex_df must contain columns: {required_cols}")
    
    # Iterate through the regex patterns and apply them to the data
    data = {}
    for _, row in regex_df.iterrows():
        keyword = row[REGEX_KEYWORD_COLUMN]
        group_name = row[REGEX_COLUMN_NAME_COLUMN]
        column_to_check = row[REGEX_KEYWORD_LOCATION_COLUMN]
        is_regex = row[REGEX_KEYWORD_TYPE_COLUMN] == REGEX_PATTERN_TYPE
        
        if is_regex:
            # Compile pattern for performance
            pattern = re.compile(keyword, re.IGNORECASE)
            data[group_name] = df[column_to_check].apply(
                lambda x: bool(pattern.search(str(x))) if pd.notna(x) else False
            ).values
        else:
            # Simple substring matching (case-insensitive)
            keyword_lower = keyword.lower()
            data[group_name] = df[column_to_check].apply(
                lambda x: keyword_lower in str(x).lower() if pd.notna(x) else False
            ).values
    
    data = pd.DataFrame.from_dict(data)
    data.index = df.index
    df = pd.concat([df, data], axis=1).reset_index(drop=True)
    return df


def process_groups(row: pd.Series) -> Union[bool, float]:
    """
    Determine group membership based on regex matches in a row.
    
    Checks if at least one match was found for the various regex patterns
    applied to a row.
    
    Args:
        row: Series of boolean or numeric values indicating regex matches
    
    Returns:
        - True if at least one value indicates a match (sum > 0)
        - False if none match
        - NaN if all values are NaN
    """
    if row.isna().all():
        return np.nan
    return row.sum() > 0


def combine_groups(row: pd.Series, columns_in_group: List[str]) -> Union[bool, float]:
    """
    Determine if a new group can be generated based on column values.
    
    Checks if all subgroups in a specified list are fulfilled for a row.
    
    Args:
        row: Series representing a row
        columns_in_group: List of column names to check in the row
    
    Returns:
        - NaN if any value in the specified columns is NaN
        - True if all subgroups are fulfilled (count equals expected length)
        - False otherwise
    """
    if row[columns_in_group].isna().any():
        return np.nan
    return row[columns_in_group].sum() >= len(columns_in_group)


def generate_parent_groups(
    df: pd.DataFrame,
    group_dictionary: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Generate parent group columns based on child group matches.
    
    Creates new columns representing parent groups by aggregating
    matches from child group columns.
    
    Args:
        df: DataFrame containing boolean or numeric group match columns
        group_dictionary: Mapping of parent group names to lists of child column names.
            Example: {"Plasticizers": ["Phthalates", "Adipates"]}
    
    Returns:
        DataFrame with new parent group columns added
        
    Raises:
        KeyError: If child columns don't exist in df
    """
    new_cols = {}
    for parent, children in group_dictionary.items():
        if not parent:
            continue
        # Validate that child columns exist
        if not all(child in df.columns for child in children):
            raise KeyError(f"Child columns {children} not found in DataFrame")
        
        new_cols[parent] = df[children].apply(process_groups, axis=1)
    
    new_cols_df = pd.DataFrame(new_cols, index=df.index)
    df = pd.concat([df, new_cols_df], axis=1).reset_index(drop=True)
    return df


def rename_columns_to_multiindex(
    df: pd.DataFrame,
    column_names: List[str],
    column_indices: List[int]
) -> pd.DataFrame:
    """
    Rename DataFrame columns to a MultiIndex using group names and index ranges.
    
    Creates a two-level MultiIndex where the first level contains group names
    and the second level contains original column names.
    
    Args:
        df: DataFrame whose columns will be renamed
        column_names: List of group names for the MultiIndex's first level.
            Must have length one less than column_indices.
        column_indices: Indices indicating boundaries of each group in columns.
            Must have length one greater than column_names, where each pair
            (column_indices[i], column_indices[i+1]) defines the slice for
            group column_names[i]
    
    Returns:
        DataFrame with columns renamed to MultiIndex (group name, original column name)
    
    Raises:
        IndexError: If column_indices are out of bounds or inconsistent
        ValueError: If column_names and column_indices lengths are inconsistent
    """
    if isinstance(df.columns, pd.MultiIndex):
        return df
    
    # Validate input dimensions
    if len(column_names) != len(column_indices) - 1:
        raise ValueError(
            f"column_names must have length {len(column_indices) - 1}, "
            f"got {len(column_names)}"
        )
    
    # Validate indices are within bounds
    if column_indices[-1] > len(df.columns):
        raise IndexError(
            f"Last column_indices value {column_indices[-1]} exceeds "
            f"DataFrame columns length {len(df.columns)}"
        )
    
    new_column_names = []
    for i, name in enumerate(column_names):
        cols = df.columns[column_indices[i]:column_indices[i+1]]
        new_column_names += [(name, col) for col in cols]

    df.columns = pd.MultiIndex.from_tuples(new_column_names)
    return df


def column_contains_patterns(
    df: pd.DataFrame,
    column_name: str,
    patterns_dictionary: Dict[str, List[str]],
    is_regex: bool = False
) -> pd.DataFrame:
    """
    Add boolean columns indicating presence of patterns in a DataFrame column.
    
    For each key in patterns_dictionary, creates a new boolean column. Each row is set to True
    if any of the associated patterns are found in the specified column, otherwise False.
    
    Args:
        df: Input DataFrame to process
        column_name: Name of the column to search for patterns
        patterns_dictionary: Dict mapping new column names to lists of patterns to search for
        is_regex: If True, treats patterns as regular expressions; if False, as literal strings
    
    Returns:
        DataFrame with additional boolean columns for each key in patterns_dictionary
        
    Raises:
        ValueError: If df is empty, column_name doesn't exist, or patterns_dictionary is empty
    """
    # Validation
    if df is None or df.empty:
        raise ValueError("df DataFrame cannot be empty")
    
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame. Available columns: {list(df.columns)}")
    
    if not patterns_dictionary or len(patterns_dictionary) == 0:
        raise ValueError("patterns_dictionary cannot be empty")
    
    df_copy = df.copy()
    
    for curr_function, patterns_list in patterns_dictionary.items():
        # Initialize the new column to False
        df_copy[curr_function] = False
        
        # Set True when any pattern is found
        for curr_pattern in patterns_list:
            mask = df_copy[column_name].str.contains(curr_pattern, regex=is_regex, case=False, na=False)
            df_copy.loc[mask, curr_function] = True
    
    return df_copy


def merge_dataframes(
    target_dataframe: pd.DataFrame,
    source_dataframe: pd.DataFrame,
    target_id: str,
    source_id: str,
    rename_columns_dict: Dict[str, str]
) -> pd.DataFrame:
    """
    Merge two DataFrames on a specified key with optional column renaming.
    
    Args:
        target_dataframe: DataFrame to merge into
        source_dataframe: DataFrame to merge from
        target_id: Column name in target DataFrame to merge on
        source_id: Column name in source DataFrame to merge on
        rename_columns_dict: Dictionary mapping original column names (in source) to new names
    
    Returns:
        Merged DataFrame with:
            - Specified columns renamed as per rename_columns_dict
            - source_id column renamed to target_id (if different)
            - Only relevant columns kept
            - NaN values in renamed columns filled with False
            
    Raises:
        ValueError: If DataFrames are empty or required columns don't exist
    """
    # Validation
    if target_dataframe is None or target_dataframe.empty:
        raise ValueError("target_dataframe cannot be empty")
    
    if source_dataframe is None or source_dataframe.empty:
        raise ValueError("source_dataframe cannot be empty")
    
    if target_id not in target_dataframe.columns:
        raise ValueError(f"target_id '{target_id}' not found in target_dataframe")
    
    if source_id not in source_dataframe.columns:
        raise ValueError(f"source_id '{source_id}' not found in source_dataframe")
    
    # Check that all columns in rename_columns_dict exist in source_dataframe
    missing_cols = [col for col in rename_columns_dict.keys() if col not in source_dataframe.columns]
    if missing_cols:
        raise ValueError(f"Columns {missing_cols} not found in source_dataframe")
    
    df = source_dataframe.copy()
    df = df.rename(columns=rename_columns_dict)
    
    # Rename source_id to target_id if they differ
    if target_id != source_id:
        df.rename(columns={source_id: target_id}, inplace=True)
    
    # Keep only relevant columns before merging
    columns_to_keep = [target_id] + list(rename_columns_dict.values())
    df = df[columns_to_keep]
    
    # Merge DataFrames
    result = target_dataframe.merge(df, on=target_id, how='left')
    
    # Fill NaN values in merged columns with False
    for col in rename_columns_dict.values():
        if col in result.columns:
            result[col] = result[col].fillna(False)
    
    return result
