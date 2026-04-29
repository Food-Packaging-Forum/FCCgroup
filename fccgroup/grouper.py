"""
Chemical Grouping Module

This module provides functionality to group chemicals based on:
1. Structural patterns (SMARTS fingerprints)
2. Functional lists (regulatory databases, chemical inventories)
3. Regex patterns (name and formula matching)

Chemical information is retrieved dynamically using CompTox for user-provided CAS IDs.
"""

from typing import List, Dict, Optional, Union, Any, Set, Tuple
from pathlib import Path

import pandas as pd
import numpy as np
from joblib import Parallel, delayed

from tqdm.auto import tqdm

from .config import GroupingConfig, GroupingMethod, InputMode
from .constants import *
from .patterns.library import fingerprints
from .patterns.methods import apply_all_patterns, generate_fingerprints
from .molecular import smiles_check, molecule_composition
from .data import load_mapping_file, load_lists, harmonize_cas_columns, harmonize_function_columns
from .data import apply_simple_regex, process_groups, combine_groups, generate_parent_groups
from .data import column_contains_patterns, merge_dataframes
from .data.constants import individual_lists, synonym_lists, regex_combination_dictionary
from .comptox import fetch_chemical_info


class ChemicalGrouper:
    """
    Main class for chemical grouping operations.
    
    This class provides methods to classify chemicals based on their CAS numbers or SMILES
    through structural, functional, and nomenclature-based grouping. Supports flexible
    grouping method selection with lazy resource loading.

    Available methods are configured through GroupingConfig.methods:
    - SMARTS: Structural pattern matching
    - LISTS: Functional list matching
    - REGEX: Name- and formula-based regex grouping

    Resources are loaded on-demand only when needed by the selected methods, minimizing
    startup time and memory usage.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        assets_path: Optional[str] = None,
        grouping_config: Optional[GroupingConfig] = None,
    ) -> None:
        """
        Initialize the ChemicalGrouper with optional lazy loading based on config.
        
        Args:
            df: DataFrame this instance is fitted to. Method feasibility and
                enrichment requirements are resolved against this schema at init time.
            assets_path: Path to the assets folder. If None, uses default relative path.
            grouping_config: GroupingConfig object specifying which methods to apply.
        """
        if df is None or df.empty:
            raise ValueError("df DataFrame cannot be empty")

        if assets_path is None:
            self.assets_path: Path = ASSETS_DIR
        else:
            self.assets_path = Path(assets_path)
        
        self.lists_path = self.assets_path / LISTS_DIR.name
        self.df = df.copy()
        
        # Configuration
        if grouping_config is None:
            raise ValueError("grouping_config is required")
        self.config: GroupingConfig = grouping_config
        self.column_mapping = self.config.column_mapping.as_dict()
        self._filter_initial_dataframe_to_mapped_columns()
        
        self.selected_methods = self.config.resolved_methods()
        # If fingerprints are provided, make sure that they 
        self._smarts_fingerprints = self._resolve_smarts_fingerprints()
        # Validate that the provided columns are actual column names in the provided dataframe
        self._validate_declared_columns()
        # Validate method feasibility and whether CompTox enrichment is required.
        self._validate_method_requirements_and_set_comptox()
        
        self._lists_loaded: bool = False
        self._all_lists: Dict[str, pd.DataFrame] = {}
        self._list_mapping: Optional[pd.DataFrame] = None
        
        self._regex_loaded: bool = False
        self._regex_patterns: Optional[pd.DataFrame] = None

    @staticmethod
    def _is_provided(value: Optional[str]) -> bool:
        return value is not None and str(value).strip() != ""

    @staticmethod
    def _is_missing_value(value: Any) -> bool:
        if value is None:
            return True
        if isinstance(value, float) and pd.isna(value):
            return True
        return str(value).strip() == ""

    def _validate_declared_columns(self) -> None:
        """Fail fast only when mapped primary identifier columns are absent."""
        missing: List[str] = []

        for key in [CANONICAL_CAS_KEY, CANONICAL_SMILES_KEY]:
            value = getattr(self.config.column_mapping, key)
            if self._is_provided(value) and str(value) not in self.df.columns:
                missing.append(str(value))

        dedup_missing = []
        for col in missing:
            if col not in dedup_missing:
                dedup_missing.append(col)
        if dedup_missing:
            raise ValueError(
                "Mapped columns were declared but not found in initialization DataFrame: "
                + ", ".join(dedup_missing)
            )

    def _filter_initial_dataframe_to_mapped_columns(self) -> None:
        """Keep only mapped source columns from the initialization DataFrame."""
        mapped_columns: List[str] = []

        for key in [CANONICAL_CAS_KEY, CANONICAL_SMILES_KEY, CANONICAL_FORMULA_KEY]:
            value = self.column_mapping.get(key)
            if self._is_provided(value):
                mapped_columns.append(str(value).strip())

        mapped_columns.extend(self.config.resolved_name_columns())

        selected_columns: List[str] = []
        for col in mapped_columns:
            if col in self.df.columns and col not in selected_columns:
                selected_columns.append(col)

        self.df = self.df.loc[:, selected_columns].copy()

    def _validate_method_requirements_and_set_comptox(self) -> None:
        """Validate that the selected grouping methods are feasible with the provided column mapping. 
        If not, raise ValueError with clear message. Also determine whether CompTox enrichment will be needed."""

        has_cas = self._is_provided(self.column_mapping.get(CANONICAL_CAS_KEY))
        has_smiles = self._is_provided(self.column_mapping.get(CANONICAL_SMILES_KEY))
        if not (has_cas or has_smiles):
            raise ValueError("At least one of CAS or SMILES mapping must be provided")

        if has_cas and str(self.column_mapping.get(CANONICAL_CAS_KEY)) not in self.df.columns:
            raise ValueError(f"Mapped CAS column '{self.column_mapping.get(CANONICAL_CAS_KEY)}' not found in initialization DataFrame")
        if has_smiles and str(self.column_mapping.get(CANONICAL_SMILES_KEY)) not in self.df.columns:
            raise ValueError(f"Mapped SMILES column '{self.column_mapping.get(CANONICAL_SMILES_KEY)}' not found in initialization DataFrame")

        only_cas = has_cas and not has_smiles
        only_smiles = has_smiles and not has_cas
        methods = self.selected_methods

        needs_comptox = False
        if only_cas and methods != {GroupingMethod.LISTS}:
            needs_comptox = True
        if only_smiles and methods != {GroupingMethod.SMARTS}:
            needs_comptox = True

        regex_names = self.config.resolved_name_columns()
        has_formula = self._is_provided(self.column_mapping.get(CANONICAL_FORMULA_KEY))
        if GroupingMethod.REGEX in methods:
            if not (regex_names and has_formula and has_smiles):
                needs_comptox = True

        if GroupingMethod.SMARTS in methods and not has_smiles:
            self.column_mapping[CANONICAL_SMILES_KEY] = SMILES_COLUMN
            needs_comptox = True
        if GroupingMethod.LISTS in methods and not has_cas:
            self.column_mapping[CANONICAL_CAS_KEY] = CAS_COLUMN
            needs_comptox = True
        if GroupingMethod.REGEX in methods:
            if not regex_names:
                needs_comptox = True
            if not has_formula:
                self.column_mapping[CANONICAL_FORMULA_KEY] = FORMULA_COLUMN
                needs_comptox = True

        self._needs_comptox_enrichment = needs_comptox
    
    def _ensure_lists_loaded(self) -> None:
        """
        Load functional lists and related resources only when needed.
        
        This is called by group_chemicals() only if the GroupingConfig specifies
        that list-based grouping should be applied. Loading is deferred until
        necessary to minimize startup time and memory usage.
        
        Resources loaded:
        - Mapping file (Mapping.xlsx) for list configuration
        - All functional lists (70+ files from assets/lists/)
        - CAS column harmonization
        - Function column harmonization for complex lists
        
        Raises:
            FileNotFoundError: If required list files are missing
            ValueError: If list files cannot be loaded
        """
        if self._lists_loaded:
            return  # Already loaded, skip
        
        if not self.config.use_lists:
            return  # Not needed for this config
        
        print("  [LOAD] Loading functional lists (first use, ~1-3 seconds)...")
        
        # Load mapping file
        mapping_path = self.assets_path / MAPPING_FILE_NAME
        if not mapping_path.exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_path}")
        
        self._list_mapping = load_mapping_file(
            str(mapping_path),
            keyword_complex=COMPLEX_LIST_KEYWORD,
            sheet_name=MAPPING_SHEET_LISTS,
        )
        
        # Load all chemical lists
        if not self.lists_path.exists():
            raise FileNotFoundError(f"Lists directory not found: {self.lists_path}")
        
        self._all_lists = load_lists(
            mapping_df=self._list_mapping,
            data_path=str(self.lists_path),
            verbose=True
        )
        
        # Harmonize CAS and function columns
        self._all_lists = harmonize_cas_columns(self._all_lists)
        self._all_lists = harmonize_function_columns(self._all_lists, self._list_mapping)
        
        self._lists_loaded = True
        print(f"  [OK] Loaded {len(self._all_lists)} functional lists")
    
    def _ensure_regex_loaded(self) -> None:
        """
        Load regex patterns only when needed.
        
        Called by group_chemicals() only when REGEX grouping is selected.
        Loading is deferred to avoid overhead for users who do not need
        name/formula pattern matching.
        
        Resources loaded:
        - Regex patterns from Mapping.xlsx "B - Keywords" sheet
        
        Raises:
            FileNotFoundError: If Mapping.xlsx is missing
            ValueError: If regex patterns cannot be loaded
        """
        if self._regex_loaded:
            return  # Already loaded, skip
        
        if not self.config.use_regex:
            return  # Not needed for this config
        
        print("  [LOAD] Loading regex patterns (first use, ~0.5 seconds)...")
        
        mapping_path = self.assets_path / MAPPING_FILE_NAME
        if not mapping_path.exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_path}")
        
        self._regex_patterns = pd.read_excel(
            str(mapping_path),
            sheet_name=MAPPING_SHEET_KEYWORDS,
        )
        
        # Filter out disabled patterns
        self._regex_patterns = self._regex_patterns.loc[
            ~(self._regex_patterns[USE_COLUMN] == False)
        ]
        
        self._regex_loaded = True
        print(f"  [OK] Loaded {len(self._regex_patterns)} regex patterns")

    @staticmethod
    def _coerce_enriched_scalar(value: Any) -> Optional[str]:
        """Convert resolver outputs to a scalar string suitable for single-cell assignment."""
        if value is None:
            return None

        if isinstance(value, (list, tuple, set)):
            for item in value:
                if not ChemicalGrouper._is_missing_value(item):
                    text = str(item).strip()
                    if text:
                        return text
            return None

        if ChemicalGrouper._is_missing_value(value):
            return None

        text = str(value).strip()
        return text if text else None

    @staticmethod
    def _ensure_internal_columns(
        df: pd.DataFrame,
        column_mapping: Dict[str, str],
        input_mode: InputMode,
    ) -> Tuple[pd.DataFrame, str]:
        """Ensure canonical columns exist and return the canonical resolver column."""
        internal_mapping = {
            CAS_COLUMN: column_mapping[CANONICAL_CAS_KEY],
            SMILES_COLUMN: column_mapping[CANONICAL_SMILES_KEY],
            FORMULA_COLUMN: column_mapping[CANONICAL_FORMULA_KEY],
        }

        remap_columns = {}

        for canonical, mapped_col in internal_mapping.items():
            if canonical in df.columns:
                continue
            if mapped_col and mapped_col in df.columns:
                remap_columns[mapped_col] = canonical

        df = df.rename(columns=remap_columns)
        resolver_col = CAS_COLUMN if input_mode == InputMode.CAS_ID else SMILES_COLUMN
        return df, resolver_col

    def _missing_fields_for_row(self, row_data: Dict[str, Any]) -> Set[str]:
        """Compute missing canonical fields required by selected methods."""
        missing: Set[str] = set()

        if GroupingMethod.SMARTS in self.selected_methods and self._is_missing_value(row_data.get(SMILES_COLUMN)):
            missing.add(SMILES_COLUMN)

        if GroupingMethod.LISTS in self.selected_methods and self._is_missing_value(row_data.get(CAS_COLUMN)):
            missing.add(CAS_COLUMN)

        if GroupingMethod.REGEX in self.selected_methods:
            if not self.config.resolved_name_columns():
                missing.add(ENRICHED_NAME_COLUMNS_COLUMN)
            if self._is_missing_value(row_data.get(FORMULA_COLUMN)):
                missing.add(FORMULA_COLUMN)

        return missing

    @staticmethod
    def _build_comptox_batches(
        enrichment_queue: List[Tuple[Any, str]],
        max_payload_chars: int = 200,
    ) -> List[List[Tuple[Any, str]]]:
        """Create ordered batches where len("\\n".join(batch_identifiers)) is strictly below max."""
        batches: List[List[Tuple[Any, str]]] = []
        current_batch: List[Tuple[Any, str]] = []
        current_payload_len = 0

        for idx, identifier in enrichment_queue:
            identifier_len = len(identifier)
            if identifier_len >= max_payload_chars:
                print(
                    "  [WARN] Identifier too long for CompTox batch payload "
                    f"(len={identifier_len}) at row {idx}; skipping enrichment"
                )
                continue

            next_payload_len = identifier_len if not current_batch else current_payload_len + 1 + identifier_len

            if current_batch and next_payload_len >= max_payload_chars:
                batches.append(current_batch)
                current_batch = []
                current_payload_len = 0
                next_payload_len = identifier_len

            current_batch.append((idx, identifier))
            current_payload_len = next_payload_len

        if current_batch:
            batches.append(current_batch)

        return batches

    def _apply_enriched_row_values(
        self,
        df: pd.DataFrame,
        idx: Any,
        enriched_row: Dict[str, Any],
    ) -> None:
        """Assign all resolver-returned canonical fields for a single row."""
        for key in [CAS_COLUMN, SMILES_COLUMN, FORMULA_COLUMN, ENRICHED_NAME_COLUMNS_COLUMN]:
            value = enriched_row.get(key)
            if key == ENRICHED_NAME_COLUMNS_COLUMN:
                if isinstance(value, (list, tuple, set)):
                    unique_names = sorted(
                        {
                            str(v).strip()
                            for v in value
                            if not self._is_missing_value(v)
                        }
                    )
                    if unique_names:
                        df.at[idx, key] = '; '.join(unique_names)
                else:
                    normalized = self._coerce_enriched_scalar(value)
                    if normalized is not None:
                        df.at[idx, key] = normalized
                continue

            normalized = self._coerce_enriched_scalar(value)
            if normalized is not None:
                df.at[idx, key] = normalized

    def _resolve_smarts_fingerprints(self) -> Dict[str, Union[str, Any]]:
        """Resolve and validate which SMARTS fingerprints should be applied."""
        selected = self.config.smarts_fingerprints
        if selected is None:
            return fingerprints

        available = set(fingerprints.keys())
        unknown = sorted(selected - available)
        if unknown:
            raise ValueError(
                "Unknown SMARTS fingerprints in GroupingConfig.smarts_fingerprints: "
                + ", ".join(unknown)
            )

        filtered_fingerprints = {
            name: fingerprints[name]
            for name in selected
        }
        return filtered_fingerprints
    
    def group_chemicals(self, save=True) -> pd.DataFrame:
        """
        Group chemicals based on their CAS IDs or SMILES strings using configured methods.
        
        Uses the DataFrame bound at initialization time as the reference source.
        CompTox is queried only for fields that are required by selected methods and
        missing for each specific chemical.
        
        Args:
            save (bool): Whether to save the results to an Excel file. Input identifiers are read directly from the configured
            resolver column in the DataFrame provided at initialization.
        
        Returns:
            DataFrame with grouping results for each input. Columns vary based on methods used:
            - Internal identifier columns such as casId and/or SMILES when resolved
            - With SMARTS: Chemical groups
            - With Lists: List name columns
            - With Regex: Pattern group columns
            
        Raises:
            ValueError: If resolver column is missing or contains no valid values
            FileNotFoundError: If required resource files not found (for list/regex modes)
        """
        input_mode = self.config.resolver

        df = self.df.copy()
        resolver_col = self.column_mapping.get(input_mode.value)
        if not resolver_col or resolver_col not in df.columns:
            raise ValueError(
                f"Resolver column '{resolver_col}' for input mode '{input_mode.value}' "
                "was not found in initialization DataFrame"
            )

        if df.empty:
            raise ValueError(
                f"No valid identifiers found in resolver column '{resolver_col}'"
            )
        
        print(f"\n{'='*60}")
        print(f"Processing {len(df)} {input_mode.value} entries")
        print(f"Grouping mode: {self.config.description}")
        print(f"{'='*60}\n")
        
        # Pre-load resources if needed (lazy loading happens on first access)
        if self.config.use_lists:
            self._ensure_lists_loaded()
        if self.config.use_regex:
            self._ensure_regex_loaded()
        
        print(f"Step 0/3: Resolving {input_mode.value} input rows...")
        df, resolver_col = self._ensure_internal_columns(df, self.column_mapping, input_mode)

        if resolver_col not in df.columns:
            raise ValueError(
                f"Canonical resolver column '{resolver_col}' could not be resolved from mapping "
                f"for input mode '{input_mode.value}'"
            )

        comptox_used = 0
        if self._needs_comptox_enrichment:
            enrichment_queue: List[Tuple[Any, str]] = []

            with tqdm(df.iterrows(), total=len(df), desc="Step 0/3: Queueing CompTox lookups", unit="chem") as pbar:
                for idx, row in pbar:
                    row_data = row.to_dict()
                    identifier = row_data.get(resolver_col)
                    if self._is_missing_value(identifier):
                        pbar.write(
                            f"  [WARN] Missing identifier in resolver column for row {idx}, "
                            "skipping enrichment"
                        )
                        continue

                    enrichment_queue.append((idx, str(identifier)))

            comptox_batches = self._build_comptox_batches(enrichment_queue, max_payload_chars=200)
            comptox_used = len(enrichment_queue)

            with tqdm(comptox_batches, desc="Step 0/3: Resolving via CompTox", unit="batch") as pbar:
                for batch_entries in pbar:
                    batch_identifiers = [identifier for _, identifier in batch_entries]
                    try:
                        enriched_by_identifier = fetch_chemical_info(batch_identifiers)
                    except Exception as exc:
                        pbar.write(
                            "  [WARN] CompTox batch request failed for "
                            f"{len(batch_identifiers)} identifiers: {str(exc)}"
                        )
                        continue

                    for idx, identifier in batch_entries:
                        enriched_row = enriched_by_identifier.get(identifier)
                        if not isinstance(enriched_row, dict):
                            continue
                        self._apply_enriched_row_values(
                            df=df,
                            idx=idx,
                            enriched_row=enriched_row,
                        )

        if comptox_used:
            print(f"  [OK] CompTox enrichment used for {comptox_used}/{len(df)} entries")
        
        df.fillna('', inplace=True)

        # Snapshot column count by position before each method step.
        # Position-based tracking correctly labels columns even when names overlap across methods.
        n_identifier = len(df.columns)

        # Track which methods are applied
        methods_applied = []

        smarts_fingerprints = self._smarts_fingerprints

        # Step 1: Always apply SMARTS structural patterns (if config allows)
        if GroupingMethod.SMARTS in self.selected_methods:
            print("Step 1/3: Applying SMARTS structural patterns...")
            id_column = SMILES_COLUMN if input_mode == InputMode.SMILES else CAS_COLUMN
            df = self._apply_structural_patterns(
                df,
                id_column=id_column,
                smiles_column=SMILES_COLUMN,
                fingerprints_dict=smarts_fingerprints,
            )
            methods_applied.append("SMARTS")

        n_after_smarts = len(df.columns)

        # Step 2: Optionally apply functional lists
        if GroupingMethod.LISTS in self.selected_methods:
            print("Step 2/3: Matching against functional lists...")
            df = self._apply_functional_lists(df, cas_column=CAS_COLUMN)
            methods_applied.append("Functional Lists")

        n_after_lists = len(df.columns)

        # Step 3: Optionally apply regex patterns
        if GroupingMethod.REGEX in self.selected_methods:
            print("Step 3/3: Applying regex patterns...")
            regex_name_columns = [col for col in self.config.resolved_name_columns() if col in df.columns]
            if not regex_name_columns:
                regex_name_columns = [col for col in [ENRICHED_NAME_COLUMNS_COLUMN] if col in df.columns]

            regex_formula_column = FORMULA_COLUMN if FORMULA_COLUMN in df.columns else None
            regex_smiles_column = SMILES_COLUMN if SMILES_COLUMN in df.columns else None

            df = self._apply_regex_patterns(
                df,
                name_columns=regex_name_columns,
                formula_column=regex_formula_column,
                smiles_column=regex_smiles_column,
            )
            methods_applied.append("Regex Patterns")

        n_after_regex = len(df.columns)

        # Build MultiIndex labels by column position — each block is a contiguous slice
        # of columns appended by that method, so there is no ambiguity from name overlap.
        labels = (
            [MULTIINDEX_IDENTIFIER_LABEL] * n_identifier
            + [MULTIINDEX_STRUCTURAL_LABEL] * (n_after_smarts - n_identifier)
            + [MULTIINDEX_LISTS_LABEL] * (n_after_lists - n_after_smarts)
            + [MULTIINDEX_REGEX_LABEL] * (n_after_regex - n_after_lists)
        )
        df.columns = pd.MultiIndex.from_tuples(zip(labels, df.columns))

        print(f"\n{'='*60}")
        print(f"Grouping completed successfully!")
        print(f"Applied methods: {', '.join(methods_applied)}")
        print(f"{'='*60}\n")

        if save:
            excel_filename = f"Grouping.xlsx"
            with pd.ExcelWriter(excel_filename) as writer:    # Write DataFrames to the workbook
                df.to_excel(writer)
        return df
    
    def _apply_structural_patterns(
        self,
        df: pd.DataFrame,
        id_column: str,
        smiles_column: str,
        fingerprints_dict: Dict[str, Union[str, Any]],
    ) -> pd.DataFrame:
        """
        Apply SMARTS structural patterns to chemicals with valid SMILES.
        
        Processes each molecule's SMILES string through the complete set of fingerprints
        (structural patterns) in parallel. Results are merged back into the DataFrame with
        columns for each fingerprint result and a summary "Chemical groups" column.
        
        Args:
            df: DataFrame containing chemical data (must include SMILES column)
            id_column: Name of ID column to use ("casId" or "SMILES") as index for results
            
        Returns:
            DataFrame with original columns plus:
                - SMILES_check: Boolean indicating if SMILES is valid
                - Fingerprint columns for each pattern
                - structural_matches_count: Number of matching patterns
                - Chemical groups: Semicolon-separated list of matched group names
                
        Performance Notes:
            - Uses JobLib Parallel with n_jobs=-1 (auto CPU count)
            - Achieves ~85-95% CPU efficiency on 8-16 core systems
            - ~5-20ms per molecule processing time
            - Speedup: ~5-7x with 8 cores
            - Parallel overhead: ~100-500ms (negligible for 1000+ molecules)
            - Optimization candidate: Adaptive parallelization (disable for <1000 molecules)
        """
        if smiles_column not in df.columns:
            return df
        
        # Validate SMILES
        df[SMILES_CHECK_COLUMN] = df[smiles_column].apply(lambda x: smiles_check(x))
        
        # Apply fingerprints in parallel
        rows = list(df.iterrows())
        results = Parallel(n_jobs=-1, verbose=0)(
            delayed(generate_fingerprints)(
                row[id_column],
                row[smiles_column],
                row[SMILES_CHECK_COLUMN],
                fingerprints_dict
            )
            for _, row in tqdm(rows, desc="  Applying SMARTS fingerprints", unit="chem")
        )
        
        fpp = pd.DataFrame(results)
        fpp.rename(columns={FINGERPRINT_ID_COLUMN: id_column}, inplace=True)
        df = df.merge(fpp, on=id_column, how="left")

        fingerprint_columns = [name for name in fingerprints_dict if name in df.columns]
        if not fingerprint_columns:
            df[STRUCTURAL_MATCHES_COUNT_COLUMN] = df[SMILES_CHECK_COLUMN].apply(
                lambda is_valid: 0 if is_valid else "Invalid SMILES"
            )
            df[OUTPUT_COLUMN] = ""
            print(f"  [OK] Processed {len(df)} chemicals for structural patterns")
            return df
        
        # Calculate matches
        df[STRUCTURAL_MATCHES_COUNT_COLUMN] = (
            df[[SMILES_CHECK_COLUMN] + fingerprint_columns] > 0
        ).sum(axis=1) - 1
        df[STRUCTURAL_MATCHES_COUNT_COLUMN] = df[STRUCTURAL_MATCHES_COUNT_COLUMN].replace(
            -1, "Invalid SMILES"
        )
        df[OUTPUT_COLUMN] = df[fingerprint_columns].apply(
            lambda x: ",".join([name for name, val in x.items() if val and not pd.isna(val)]),
            axis=1
        )

        print(f"  [OK] Processed {len(df)} chemicals for structural patterns")
        return df
    
    def _apply_functional_lists(self, df: pd.DataFrame, cas_column: str) -> pd.DataFrame:
        """
        Match chemicals against functional regulatory and chemical inventory lists.
        
        Performs two types of matching:
        1. Simple lists: Direct CAS number matching (boolean presence/absence)
        2. Complex lists: Functional group matching with pattern recognition
           - G04: Food packaging functions (literal string matching)
           - G13: Fragrance applications (regex pattern matching)
           - G23: PlastChem categories (boolean column indicators)
           - G24: PlasticMAP functions (pivot table expansion)
           - G25: Plastic functions (regex pattern matching)
        
        Args:
            df: DataFrame containing chemical data (must include casId column)
            
        Returns:
            DataFrame with added list-based grouping columns
        """
        if cas_column not in df.columns or not self._all_lists:
            return df
        
        # Simple lists - direct CAS matching
        simple_lists = list(self._list_mapping.loc[~self._list_mapping[IS_COMPLEX_COLUMN], LIST_NAME_COLUMN]) \
            if self._list_mapping is not None else []
        
        for curr_list in tqdm(simple_lists, desc="  Matching simple lists", unit="list"):
            if curr_list in self._all_lists:
                curr_list_df = self._all_lists[curr_list]
                if CAS_COLUMN in curr_list_df.columns:
                    df[curr_list] = df[cas_column].isin(curr_list_df[CAS_COLUMN])
        
        # Complex lists - function-based matching
        complex_lists_handled = self._apply_complex_lists(df, cas_column=cas_column)
        
        total_lists = len(simple_lists) + complex_lists_handled
        print(f"  [OK] Matched against {len(simple_lists)} simple lists + {complex_lists_handled} complex lists = {total_lists} total")
        return df
    
    def _apply_complex_lists(self, df: pd.DataFrame, cas_column: str) -> int:
        """
        Apply complex list matching with function-based grouping.
        
        Handles special processing for complex chemical lists that contain
        functional subcategories:
        
        - G04 (Food Packaging): Uses literal string matching for function synonyms
        - G13 (Fragrances): Uses regex patterns for product type matching
        - G23 (PlastChem): Uses boolean indicator columns for chemical categories
        - G24 (PlasticMAP): Pivots function_group_in to create one column per function
        - G25 (Plastics): Uses regex patterns for function matching
        
        Each complex list is merged into the main DataFrame with renamed columns
        to include the list identifier (e.g., 'Antioxidant_G13').
        
        Args:
            df: DataFrame to add complex list columns to (modified in place)
            
        Returns:
            Number of complex lists successfully processed
        """
        if not self._all_lists:
            return 0
        
        data_column_name = FUNCTION_GROUP_COLUMN
        lists_processed = 0
        
        # G04, G13, G25: Pattern-based matching on function_group_in column
        complex_lists = PATTERN_COMPLEX_LIST_IDS
        pattern_dictionaries = [
            individual_lists.get("function_synonyms_g4", {}),
            synonym_lists.get("function_synonyms_regex", {}),
            synonym_lists.get("function_synonyms_regex", {}),
        ]
        regex_flags = [False, True, True]
        
        for dataset_name, patterns_dict, is_regex in zip(complex_lists, pattern_dictionaries, regex_flags):
            if dataset_name not in self._all_lists:
                continue
            
            if not patterns_dict:
                continue
            
            try:
                curr_list_df = self._all_lists[dataset_name]
                
                # Check if required column exists
                if data_column_name not in curr_list_df.columns:
                    continue
                
                # Apply pattern matching
                curr_list_df = column_contains_patterns(
                    curr_list_df,
                    column_name=data_column_name,
                    patterns_dictionary=patterns_dict,
                    is_regex=is_regex
                )
                
                # Rename columns to include list identifier
                rename_dict = {col: f"{col}_{dataset_name}" for col in patterns_dict}
                
                # Merge into main DataFrame
                df_merged = merge_dataframes(
                    target_dataframe=df,
                    source_dataframe=curr_list_df,
                    target_id=cas_column,
                    source_id=CAS_COLUMN,
                    rename_columns_dict=rename_dict
                )
                
                # Update df in place by adding the new columns
                for new_col in rename_dict.values():
                    if new_col in df_merged.columns:
                        df[new_col] = df_merged[new_col]
                
                lists_processed += 1
                
            except Exception as e:
                print(f"  [WARN] Could not process complex list {dataset_name}: {str(e)}")
        
        # G23: PlastChem boolean indicator columns
        dataset_name = LIST_G23
        if dataset_name in self._all_lists:
            try:
                platschem_columns = individual_lists.get(PLASTCHEM_COLUMNS_KEY, [])
                
                if platschem_columns:
                    curr_list_df = self._all_lists[dataset_name]
                    
                    # Keep only casId and PlastChem columns
                    available_cols = [CAS_COLUMN] + [col for col in platschem_columns if col in curr_list_df.columns]
                    curr_list_df = curr_list_df[available_cols].copy()
                    
                    # Convert to boolean: 1 -> True, otherwise False
                    curr_list_df[platschem_columns] = curr_list_df[platschem_columns] == 1
                    
                    # Rename columns to include list identifier
                    rename_dict = {col: f"{col}_{dataset_name}" for col in platschem_columns if col in curr_list_df.columns}
                    
                    # Merge into main DataFrame
                    df_merged = merge_dataframes(
                        target_dataframe=df,
                        source_dataframe=curr_list_df,
                        target_id=cas_column,
                        source_id=CAS_COLUMN,
                        rename_columns_dict=rename_dict
                    )
                    
                    # Update df in place
                    for new_col in rename_dict.values():
                        if new_col in df_merged.columns:
                            df[new_col] = df_merged[new_col]
                    
                    lists_processed += 1
                    
            except Exception as e:
                print(f"  [WARN] Could not process complex list {dataset_name}: {str(e)}")
        
        # G24: PlasticMAP - pivot function_group_in column
        dataset_name = LIST_G24
        if dataset_name in self._all_lists:
            try:
                curr_list_df = self._all_lists[dataset_name]
                
                if CAS_COLUMN in curr_list_df.columns and FUNCTION_GROUP_COLUMN in curr_list_df.columns:
                    # Create pivot table: one column per unique function
                    pivoted = curr_list_df[[CAS_COLUMN, FUNCTION_GROUP_COLUMN]].pivot_table(
                        index=CAS_COLUMN,
                        columns=FUNCTION_GROUP_COLUMN,
                        aggfunc=lambda x: True
                    ).reset_index()
                    
                    pivoted.fillna(False, inplace=True)
                    
                    # Get unique function names for renaming
                    unique_functions = curr_list_df[FUNCTION_GROUP_COLUMN].unique()
                    rename_dict = {col: f"{col}_{dataset_name}" for col in unique_functions if col in pivoted.columns}
                    
                    # Merge into main DataFrame
                    df_merged = merge_dataframes(
                        target_dataframe=df,
                        source_dataframe=pivoted,
                        target_id=cas_column,
                        source_id=CAS_COLUMN,
                        rename_columns_dict=rename_dict
                    )
                    
                    # Update df in place
                    for new_col in rename_dict.values():
                        if new_col in df_merged.columns:
                            df[new_col] = df_merged[new_col]
                    
                    lists_processed += 1
                    
            except Exception as e:
                print(f"  [WARN] Could not process complex list {dataset_name}: {str(e)}")
        
        return lists_processed
    
    def _apply_regex_patterns(
        self,
        df: pd.DataFrame,
        name_columns: List[str],
        formula_column: Optional[str],
        smiles_column: Optional[str]
    ) -> pd.DataFrame:
        """
        Apply regex patterns to chemical names and molecular formulas.
        
        Uses patterns from the Mapping.xlsx "B - Keywords" sheet to classify chemicals
        by their nomenclature and chemical composition. Generates hierarchical groupings
        through multiple stages:
        
        1. Pattern Matching: Applies regex patterns to name and formula columns
        2. Group Generation: Creates parent groups by aggregating child patterns
        3. Super Group Generation: Creates higher-level categories
        4. Derived Groups: Generates special categories (organic/inorganic, metals, UVCBs)
        5. Combination Groups: Applies logical AND operations for compound categories
        
        Args:
            df: DataFrame containing chemical data
            name_columns: List of column names to combine for name matching.
                All available columns will be combined with ' | ' separator and lowercased.
            formula_column: Column name containing molecular formula.
                Will be renamed to 'MolFormula' for pattern matching.
            smiles_column: Column name containing SMILES representation.
                Will be renamed to 'SMILES' for pattern matching.
        Returns:
            DataFrame with added regex-based grouping columns organized hierarchically
        """
        initial_ncols = len(df.columns)

        if self._regex_patterns is None or self._regex_patterns.empty:
            return df
        
        # Combine all available name columns with ' | ' separator
        if name_columns:
            df[COMBINED_NAME_COLUMN] = df[name_columns[0]].fillna('').astype(str)
            for col in name_columns[1:]:
                df[COMBINED_NAME_COLUMN] = df[COMBINED_NAME_COLUMN] + ' | ' + df[col].fillna('').astype(str)
            df[COMBINED_NAME_COLUMN] = df[COMBINED_NAME_COLUMN].str.lower()
        else:
            df[COMBINED_NAME_COLUMN] = ''

        # Ensure optional formula/smiles columns exist with fallback values
        if formula_column and formula_column in df.columns:
            df[MOLECULAR_FORMULA_COLUMN] = df[formula_column]
        elif MOLECULAR_FORMULA_COLUMN not in df.columns:
            df[MOLECULAR_FORMULA_COLUMN] = ''

        if smiles_column and smiles_column in df.columns:
            df[SMILES_COLUMN] = df[smiles_column]
        elif SMILES_COLUMN not in df.columns:
            df[SMILES_COLUMN] = ''
        
        # Stage 1: Apply all regex patterns from mapping file
        try:
            df = apply_simple_regex(df, self._regex_patterns)
            pattern_cols = self._regex_patterns[REGEX_COLUMN_NAME_COLUMN].unique()
            print(f"  [OK] Applied {len(self._regex_patterns)} regex patterns -> {len(pattern_cols)} pattern columns")
        except Exception as e:
            print(f"  [WARN] Could not apply regex patterns: {str(e)}")
            return df
        
        # Stage 2: Generate parent groups from pattern groupings
        try:
            group_dict = self._regex_patterns.groupby(REGEX_GROUP_COLUMN).agg(list)[REGEX_COLUMN_NAME_COLUMN].to_dict()
            df = generate_parent_groups(df, group_dict)
            print(f"  [OK] Generated {len(group_dict)} parent groups")
        except Exception as e:
            print(f"  [WARN] Could not generate parent groups: {str(e)}")
        
        # Stage 3: Generate super groups (highest level aggregation)
        try:
            super_group_dict = self._regex_patterns.groupby(REGEX_SUPER_GROUP_COLUMN).agg(list)[REGEX_COLUMN_NAME_COLUMN].to_dict()
            df = generate_parent_groups(df, super_group_dict)
            print(f"  [OK] Generated {len(super_group_dict)} super groups")
        except Exception as e:
            print(f"  [WARN] Could not generate super groups: {str(e)}")
        
        # Stage 4: Create derived groups based on existing patterns
        try:
            # Organic compounds (contains carbon)
            if 'Carbon' in df.columns:
                df['Organic_C'] = df['Carbon']
                df['Inorganic_noC'] = df['Carbon'].apply(
                    lambda x: not x if pd.notna(x) else np.nan
                )
            
            # Metal/Metalloid aggregation
            metal_metalloid_components = []
            if 'Metal' in df.columns:
                metal_metalloid_components.append('Metal')
            if 'Metalloid' in df.columns:
                metal_metalloid_components.append('Metalloid')
            
            if metal_metalloid_components:
                df = generate_parent_groups(df, {'Metal_Metalloid': metal_metalloid_components})
            
            # UVCB aggregation (Unknown or Variable composition, Complex reaction products, Biological materials)
            uvcb_components = []
            if 'UVCB - Biological origin' in df.columns:
                uvcb_components.append('UVCB - Biological origin')
            if 'UVCB - Process based' in df.columns:
                uvcb_components.append('UVCB - Process based')
            
            if uvcb_components:
                df = generate_parent_groups(df, {'UVCBs': uvcb_components})
            
            derived_count = sum([
                'Organic_C' in df.columns,
                'Inorganic_noC' in df.columns,
                'Metal_Metalloid' in df.columns,
                'UVCBs' in df.columns
            ])
            print(f"  [OK] Created {derived_count} derived groups (Organic/Inorganic, Metals, UVCBs)")
            
        except Exception as e:
            print(f"  [WARN] Could not create derived groups: {str(e)}")
        
        # Stage 5: Apply combination dictionary for compound classifications
        try:
            combinations_applied = 0
            for parent, children in regex_combination_dictionary.items():
                # Check if all required child columns exist
                if all(child in df.columns for child in children):
                    df[parent] = df[children].apply(
                        lambda x: combine_groups(x, children), axis=1
                    )
                    combinations_applied += 1
            
            if combinations_applied > 0:
                print(f"  [OK] Applied {combinations_applied} combination rules (OrganoMetallic, Salts)")
            
        except Exception as e:
            print(f"  [WARN] Could not apply combination groups: {str(e)}")
        
        total_regex_cols = len([col for col in df.columns if col])
        print(f"  [OK] Total regex-based columns added: {total_regex_cols - initial_ncols}")
        
        return df

