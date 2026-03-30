import pandas as pd

from fccgroup.constants import FILE_COLUMN, LIST_NAME_COLUMN, READ_PARAMETERS_COLUMN
from fccgroup.data.lists import load_lists


def test_load_lists_skips_missing_plasticmap_file(tmp_path):
    mapping_df = pd.DataFrame(
        [
            {
                LIST_NAME_COLUMN: "PlasticMAP_G24",
                FILE_COLUMN: "plasticmap_missing.xlsx",
                READ_PARAMETERS_COLUMN: "{}",
            }
        ]
    )

    loaded = load_lists(mapping_df=mapping_df, data_path=str(tmp_path), verbose=False)

    assert loaded == {}
