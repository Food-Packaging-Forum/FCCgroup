"""Tests that the groups of concern map onto patterns that actually exist."""

from fccgroup.constants import (
    GROUPS_CONCERN_DICT,
    IMPORTANT_GROUPS,
    PRIORITY_GROUPS_CONCERN,
    EXPOSURE_GROUPS_CONCERN,
    EXTERNAL_GROUPS_CONCERN,
)


def _all_groups():
    return [group for groups in GROUPS_CONCERN_DICT.values() for group in groups]


def test_groups_concern_not_empty():
    assert _all_groups()


def test_groups_concern_are_implemented_fingerprints(fingerprints):
    missing = sorted(set(_all_groups()) - set(fingerprints))
    assert not missing, f"Groups of concern without a fingerprint implementation: {missing}"


def test_groups_concern_are_active_in_smarts_sheet(pattern_df):
    missing = sorted(set(_all_groups()) - set(pattern_df["Group"]))
    assert not missing, f"Groups of concern not active (use=True) in B-SMARTS: {missing}"


def test_groups_concern_have_no_duplicates():
    groups = _all_groups()
    duplicated = sorted({g for g in groups if groups.count(g) > 1})
    assert not duplicated, f"Groups listed more than once: {duplicated}"


def test_groups_concern_flags_are_booleans():
    flags = {"priority_group", "exposure_group", "external_group"}
    for groups in GROUPS_CONCERN_DICT.values():
        for group, info in groups.items():
            assert set(info) == flags, group
            assert all(isinstance(v, bool) for v in info.values()), group


def test_derived_group_lists_are_consistent():
    assert IMPORTANT_GROUPS == _all_groups()
    for subset in (PRIORITY_GROUPS_CONCERN, EXPOSURE_GROUPS_CONCERN, EXTERNAL_GROUPS_CONCERN):
        assert set(subset) <= set(IMPORTANT_GROUPS)
