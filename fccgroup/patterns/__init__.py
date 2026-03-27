"""
Patterns and fingerprinting subpackage.

Provides SMARTS pattern libraries and fingerprinting utilities for
chemical classification based on molecular structure.
"""

from .library import (
    fingerprints,
    cx_smarts_query,
    count_pattern_occurrences,
    apply_pattern,
)

from .methods import generate_fingerprints, apply_all_patterns

__all__ = [
    "fingerprints",
    "cx_smarts_query",
    "count_pattern_occurrences",
    "apply_pattern",
    "generate_fingerprints",
    "apply_all_patterns",
]
