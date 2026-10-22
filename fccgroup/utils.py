"""Build the groups-of-concern dictionary from Mapping.xlsx.

Kept free of intra-package imports so ``constants`` can use it without
creating an import cycle.
"""

from pathlib import Path

import pandas as pd


def build_groups_concern_dict(
    mapping_path: Path,
    smarts_sheet_name: str,
    important_sheet_name: str,
    supergroup_column: str,
    column_name: str,
    priority_column: str,
    exposure_column: str,
    external_column: str,
) -> dict:
    """Build the groups-of-concern dictionary from Mapping.xlsx.

    Groups are the rows of the important-groups sheet, nested under their
    ``Supergroup`` as found in the SMARTS sheet (keyed by ``Column_name``)::

        {Supergroup: {Column_name: {'priority_group': bool,
                                    'exposure_group': bool,
                                    'external_group': bool}}}
    """
    smarts = pd.read_excel(
        mapping_path,
        sheet_name=smarts_sheet_name,
        header=1,
        usecols=[supergroup_column, column_name],
    )
    supergroups = smarts.drop_duplicates(column_name).set_index(column_name)[supergroup_column]

    important = pd.read_excel(
        mapping_path,
        sheet_name=important_sheet_name,
        header=1,
        usecols=[column_name, priority_column, exposure_column, external_column],
    ).dropna(subset=[column_name])

    unknown = sorted(set(important[column_name]) - set(supergroups.index))
    if unknown:
        raise ValueError(
            f"Groups in '{important_sheet_name}' missing from '{smarts_sheet_name}': {unknown}"
        )

    groups_dict: dict = {}
    for group, priority, exposure, external in important[
        [column_name, priority_column, exposure_column, external_column]
    ].itertuples(index=False):
        groups_dict.setdefault(supergroups[group], {})[group] = {
            'priority_group': bool(priority),
            'exposure_group': bool(exposure),
            'external_group': bool(external),
        }
    return groups_dict
