"""
FCCgroup - Chemical Grouping Package

A comprehensive package for grouping chemicals based on structural patterns,
functional lists, and nomenclature-based classification.

Main exports:
    - ChemicalGrouper: Main class for chemical grouping operations
    - GroupingMethod: Enumeration for atomic grouping methods
    - GroupingConfig: Configuration object for grouping behavior
    - ColumnMapping: Input column mapping configuration
"""

from .grouper import ChemicalGrouper
from .config import GroupingMethod, GroupingConfig, ColumnMapping

__version__ = "0.1.0"
__author__ = "Albert Anguera Sempere, Helene Wiesinger"
__organization__ = "Food Packaging Forum"

__all__ = [
    "ChemicalGrouper",
    "GroupingMethod",
    "GroupingConfig",
    "ColumnMapping",
]
