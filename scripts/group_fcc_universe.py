import os
import sys
sys.path

if "..\\fccgroup" not in sys.path:
    sys.path.insert(0, "..\\fccgroup")

from fccgroup import ChemicalGrouper, GroupingMethod, GroupingConfig, ColumnMapping
import pandas as pd

df = pd.read_excel("tests/FCCuniverse.xlsx")

config = GroupingConfig(
    methods=[GroupingMethod.SMARTS, GroupingMethod.LISTS, GroupingMethod.REGEX],
    column_mapping=ColumnMapping(**{
        "cas": "casId",
        "smiles": "SMILES",
        "name_columns": ["Name"],
        "formula": "Mol Formula",
    })
)

grouper = ChemicalGrouper(df, assets_path="fccgroup/assets", grouping_config=config)

grouper.group_chemicals()
