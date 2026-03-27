"""
Chemical list and mapping file loading utilities.

This module provides functions for loading and processing regulatory and
chemical inventory lists from Excel and other data sources.
"""

import json
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from ..constants import *


def load_mapping_file(
    path: str,
    keyword_complex: str = COMPLEX_LIST_KEYWORD,
    sheet_name: str = MAPPING_SHEET_LISTS,
) -> pd.DataFrame:
    """
    Load and process mapping configuration from Excel file.
    
    The mapping file defines which lists to load and their properties.
    Filters to only active lists (use == True) and creates list_name identifiers.
    
    Args:
        path: Path to the mapping Excel file
        keyword_complex: Keyword identifying complex lists with multiple subgroups
        sheet_name: Name of the sheet containing list mappings
    
    Returns:
        DataFrame with mapping configuration containing:
            - 'list_name': Identifier combining Group_key and list_id
            - 'is_complex': Boolean indicating if list contains sublists
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If required columns are missing
    """
    try:
        mapping_df = pd.read_excel(path, sheet_name=sheet_name)
    except FileNotFoundError:
        raise FileNotFoundError(f"Mapping file not found: {path}")
    except Exception as e:
        raise ValueError(f"Error loading mapping file: {str(e)}")
    
    # Filter to only active lists (use == True)
    if USE_COLUMN in mapping_df.columns:
        mapping_df = mapping_df.loc[mapping_df[USE_COLUMN] == True].reset_index(drop=True)
    
    # Create list_name column combining Group_key and list_id
    if GROUP_KEY_COLUMN in mapping_df.columns and LIST_ID_COLUMN in mapping_df.columns:
        mapping_df[LIST_NAME_COLUMN] = mapping_df[GROUP_KEY_COLUMN] + '_' + mapping_df[LIST_ID_COLUMN]
        mapping_df[LIST_NAME_COLUMN] = mapping_df[LIST_NAME_COLUMN].str.replace(keyword_complex + '_', '')
    
    # Flag complex lists (lists with sublists of groups/functions)
    if GROUP_KEY_COLUMN in mapping_df.columns:
        mapping_df[IS_COMPLEX_COLUMN] = mapping_df[GROUP_KEY_COLUMN] == keyword_complex
    else:
        mapping_df[IS_COMPLEX_COLUMN] = False
    
    return mapping_df


def load_lists(
    mapping_df: pd.DataFrame,
    data_path: str,
    verbose: bool = False
) -> Dict[str, pd.DataFrame]:
    """
    Load all chemical lists based on mapping configuration.
    
    Loads each list file specified in the mapping DataFrame from the data directory.
    
    Args:
        mapping_df: DataFrame with list mapping configuration (from load_mapping_file)
        data_path: Path to the directory containing list files
        verbose: If True, prints progress information
    
    Returns:
        Dictionary mapping list names to loaded DataFrames
        
    Raises:
        FileNotFoundError: If any list file cannot be found
        ValueError: If list files cannot be loaded
    """
    all_lists: Dict[str, pd.DataFrame] = {}
    
    for _, row in mapping_df.iterrows():
        list_name = row.get(LIST_NAME_COLUMN, 'Unknown')
        filename = row.get(FILE_COLUMN, '')
        read_parameters = row.get(READ_PARAMETERS_COLUMN, {})

        
        if LIST_G24 in list_name:
            if verbose:
                print(f"✓ Handling {list_name} separately")
            df = load_g24_list(row, data_path)
            all_lists[list_name] = df
            continue

        try:
            read_parameters = json.loads(read_parameters)
        except json.JSONDecodeError as e:
            raise json.JSONDecodeError(
                f"Invalid JSON in 'read_parameters' for {list_name}: {str(e)}",
                e.doc, e.pos
            )
        
        if not filename:
            continue
        
        try:
            file_path = Path(data_path) / filename
            
            # Determine file type
            if filename.endswith('.xlsx') or filename.endswith('.xls'):
                df = pd.read_excel(file_path, **read_parameters)
            elif filename.endswith('.tsv'):
                df = pd.read_csv(file_path, delimiter='\t', **read_parameters)
            elif filename.endswith('.csv'):
                df = pd.read_csv(file_path, **read_parameters)
            else:
                if verbose:
                    print(f"Skipping {filename}: Unsupported file format")
                continue
            
            all_lists[list_name] = df
            
            if verbose:
                print(f"✓ Loaded {list_name}: {len(df)} rows")
        
        except FileNotFoundError:
            if verbose:
                print(f"✗ Not found: {filename}")
        except Exception as e:
            if verbose:
                print(f"✗ Error loading {filename}: {str(e)}")
    
    return all_lists


def harmonize_function_columns(all_lists: Dict[str, pd.DataFrame], list_mapping: Optional[pd.DataFrame] = None) -> None:
    """
    Harmonize function/group columns for complex chemical lists.
    
    Different lists use different conventions for storing chemical functions
    and application groups. This method:
    
    1. Filters specific lists to relevant data
    2. Combines multiple function columns into a single 'function_group_in' column
    3. Creates boolean columns for special cases (endogenous compounds, etc.)
    
    Complex lists handled:
    - G04: Group name (food packaging functions)
    - G13: Product Type (fragrance applications)
    - G23: Harmonized functions (plasticizers)
    - G24: PlasticMAP functions
    - G25: Plastic functions
    """
    # Filter specific lists
    if LIST_PESTICIDES_G17 in all_lists:
        curr_list_df = all_lists[LIST_PESTICIDES_G17]
        all_lists[LIST_PESTICIDES_G17] = curr_list_df.loc[
            curr_list_df['Status under Reg. (EC) No 1107/2009'] == 'Approved'
        ]
    
    if LIST_METABOLITES_G22 in all_lists:
        G22_cols_remove = [
            'Subsystem: Drug metabolism', 'Subsystem: Alkaloid synthesis',
            'Subsystem: Miscellaneous', 'Subsystem: Stilbene, coumarine and lignin synthesis',
            'Subsystem: Selenoamino acid metabolism', 'Subsystem: Sphingolipid metabolism',
            'Subsystem: Limonene and pinene degradation', 'Subsystem: R group synthesis',
            'Subsystem: Xenobiotics metabolism'
        ]
        curr_list_df = all_lists[LIST_METABOLITES_G22]
        G22_cols = curr_list_df.columns
        G22_cols_updated = [
            col for col in G22_cols
            if 'Subsystem' in col and col not in G22_cols_remove
        ]
        curr_list_df['is_endogenous'] = curr_list_df[G22_cols_updated].any(axis=1)
        all_lists[LIST_METABOLITES_G22] = curr_list_df
    
    # Combine function columns for complex lists
    if list_mapping is not None:
        complex_lists = list(list_mapping.loc[list_mapping[IS_COMPLEX_COLUMN], LIST_NAME_COLUMN])
        for curr_list in complex_lists:
            if curr_list in all_lists:
                curr_list_df = all_lists[curr_list]
                function_cols = LIST_FUNCTION_COLUMNS.get(curr_list, [])
                
                if function_cols:
                    curr_list_df[FUNCTION_GROUP_COLUMN] = curr_list_df[function_cols].apply(
                        lambda x: ";".join(x.dropna().astype(str)), axis=1
                    )
                    all_lists[curr_list] = curr_list_df
    return all_lists

def load_g24_list(
        row: pd.Series,
        data_path: str
) -> pd.DataFrame:
    """
    Load the G24 list of chemicals with specific handling.
    
    This function is a placeholder for the actual loading logic for the G24 list,
    which may involve special processing steps or data sources.
    
    Args:
        row: The row from the mapping DataFrame corresponding to the G24 list
        data_path: Path to the directory containing list files

    Returns:
        DataFrame containing the G24 list of chemicals
    """
    filename = row.get(FILE_COLUMN, '')
    file_path = Path(data_path) / filename
    ## Load PlasticMAP database ('G24') sheets and merge them into one Dataframe
    G24_ChemicalFunction = pd.read_excel(file_path, sheet_name=PLASTICMAP_SHEET_CHEMICAL_FUNCTION)
    G24_Chemical = pd.read_excel(file_path, sheet_name=PLASTICMAP_SHEET_CHEMICAL)
    G24_Function = pd.read_excel(file_path, sheet_name=PLASTICMAP_SHEET_FUNCTION)
    G24_Source = pd.read_excel(file_path, sheet_name=PLASTICMAP_SHEET_SOURCE)

    ### Merge the sheets into one dataframe
    G24_df = pd.merge(G24_ChemicalFunction, G24_Chemical[['id','CAS','name']].rename(columns={'id':'Chemical_id'}), on='Chemical_id', how='left') # merge ChemicalFunction with Chemical
    G24_df = pd.merge(G24_df, G24_Function[['id','function_name']].rename(columns={'id':'Function_id'}), on='Function_id', how='left') # merge with Function
    G24_df = pd.merge(G24_df, G24_Source[['id','name']].rename(columns={'id':'Source_id','name':'Source_name'}), on='Source_id', how='left') # merge with Source
    return pd.DataFrame(G24_df)

