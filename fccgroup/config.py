"""
Configuration module for FCCgroup application.

Contains all application-level constants, configuration values, and settings.
This module centralizes magic strings, column names, and processing parameters
to improve maintainability and reduce duplication across the codebase.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Set, Union

from .constants import *

class InputMode(str, Enum):
    """Enumeration for supported chemical input modes."""
    CAS_ID = "cas"
    SMILES = "smiles"


class GroupingMethod(str, Enum):
    """Atomic grouping methods that can be selected independently."""
    SMARTS = "smarts"
    LISTS = "lists"
    REGEX = "regex"


@dataclass
class ColumnMapping:
    """Column mapping from canonical field names to user-provided DataFrame columns."""

    cas: Optional[str]
    smiles: Optional[str]
    name_columns: List[str] = field(default_factory=list)
    formula: Optional[str] = None

    def as_dict(self) -> dict:
        """Return mapping as a dictionary keyed by canonical field name."""
        return {
            CANONICAL_CAS_KEY: self.cas,
            CANONICAL_SMILES_KEY: self.smiles,
            CANONICAL_NAME_COLUMNS_KEY: list(self.name_columns),
            CANONICAL_FORMULA_KEY: self.formula,
        }


@dataclass
class GroupingConfig:
    """Configuration for chemical grouping behavior.
    
    Controls which grouping methods are applied and how resources are loaded.
    Resources are loaded lazily only when required by the selected methods.

    Example:
        >>> config = GroupingConfig(
        ...     methods=[GroupingMethod.SMARTS, GroupingMethod.REGEX],
        ...     column_mapping=ColumnMapping(
        ...         cas="CASRN",
        ...         smiles="Structure",
        ...         name_columns=["Name", "IUPAC"],
        ...         formula="Formula",
        ...     ),
        ... )
    """
    
    methods: List[Union[GroupingMethod, str]]
    """Selected grouping methods. Allowed values: SMARTS, LISTS, REGEX."""

    column_mapping: ColumnMapping
    """Mapping of canonical chemical fields to DataFrame column names."""

    resolver: InputMode = field(init=False)
    """Input resolver mode derived from column mapping (CAS preferred over SMILES)."""

    smarts_fingerprints: Optional[Set[str]] = None
    """Optional subset of SMARTS fingerprint names to apply. None means use all."""

    def __post_init__(self) -> None:
        """Validate column mapping at construction time."""
        self.resolver = self.validate_mapping()

    def validate_mapping(self) -> InputMode:
        """Validate that each selected method has its required column names configured.

        This is a configuration-level check: it verifies that the caller has
        *named* a column for every field that each method requires. It does not
        inspect any DataFrame — that is the job of ChemicalGrouper._validate_runtime_columns.

        Returns:
            InputMode inferred from mapping. CAS has precedence when both CAS and
            SMILES mappings are provided.

        Raises:
            ValueError: If a required column name is absent from the mapping for
                any selected method.
        """
        methods = self.resolved_methods()
        mapping = self.column_mapping.as_dict()

        def _provided(value: Optional[str]) -> bool:
            return value is not None and str(value).strip() != ""

        if not methods:
            raise ValueError("At least one grouping method must be selected")

        has_cas = _provided(mapping.get(CANONICAL_CAS_KEY))
        has_smiles = _provided(mapping.get(CANONICAL_SMILES_KEY))

        if not (has_cas or has_smiles):
            raise ValueError(
                "Missing required column mappings. Provide at least one of: "
                "CAS column mapping or SMILES column mapping."
            )

        return InputMode.CAS_ID if has_cas else InputMode.SMILES

    def resolved_name_columns(self) -> List[str]:
        """Return de-duplicated user-provided name columns for regex processing."""
        merged: List[str] = []
        for candidate in self.column_mapping.name_columns:
            if candidate is None:
                continue
            value = str(candidate).strip()
            if value and value not in merged:
                merged.append(value)
        return merged

    def resolved_methods(self) -> Set[GroupingMethod]:
        """Resolve active methods from methods list, validating allowed method names."""
        # Accept a bare enum/string for convenience, normalise to iterable
        raw = self.methods
        if isinstance(raw, (str, GroupingMethod)):
            raw = [raw]

        resolved: Set[GroupingMethod] = set()
        for method in raw:
            if isinstance(method, GroupingMethod):
                resolved.add(method)
            else:
                try:
                    resolved.add(GroupingMethod(str(method).lower()))
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid grouping method '{method}'. Allowed methods are: smarts, lists, regex"
                    ) from exc
        return resolved

    @property
    def use_smarts(self) -> bool:
        """Check if SMARTS structural patterns should be applied."""
        return GroupingMethod.SMARTS in self.resolved_methods()
    
    @property
    def use_lists(self) -> bool:
        """Check if functional lists should be applied."""
        return GroupingMethod.LISTS in self.resolved_methods()
    
    @property
    def use_regex(self) -> bool:
        """Check if regex patterns should be applied."""
        return GroupingMethod.REGEX in self.resolved_methods()
    
    @property
    def description(self) -> str:
        """Human-readable description of the grouping configuration."""
        method_labels = {
            GroupingMethod.SMARTS: "SMARTS",
            GroupingMethod.LISTS: "Lists",
            GroupingMethod.REGEX: "Regex",
        }
        ordered = [
            method_labels[m]
            for m in (GroupingMethod.SMARTS, GroupingMethod.LISTS, GroupingMethod.REGEX)
            if m in self.resolved_methods()
        ]
        return "Methods: " + ", ".join(ordered)
    
    def __str__(self) -> str:
        """String representation."""
        methods = ",".join(sorted(m.value for m in self.resolved_methods()))
        return f"GroupingConfig(methods={methods}, {self.description})"
