"""CompTox API integration helpers."""

from typing import Any, Dict, List, Optional, Set, Union
import os

import requests

from .config import InputMode
from .constants import (
    CAS_COLUMN,
    COMPTOX_API_KEY_ENV,
    COMPTOX_DTXCID_DETAIL_URL,
    COMPTOX_SEARCH_EQUAL_URL,
    ENRICHED_NAME_COLUMNS_COLUMN,
    FORMULA_COLUMN,
    SMILES_COLUMN,
)


NO_DATA = "No data"
API_KEY = os.getenv(COMPTOX_API_KEY_ENV, "")
HEADERS = {
    "accept": "application/json",
    "content-type": "application/json",
    "x-api-key": API_KEY,
}

def _is_missing_value(value: Any) -> bool:
    if value is None:
        return True
    return str(value).strip() == ""

# TODO: handle request errors

def _fetch_dtxcid_records(identifiers: List[str]) -> List[str]:
    response = requests.post(
        COMPTOX_SEARCH_EQUAL_URL,
        headers=HEADERS,
        data="\n".join(identifiers).strip(),
        timeout=30,
    )
    response.raise_for_status()
    dtxcids = []
    for entry in response.json():
        dtxcids.append(entry["dtxcid"])
    return dtxcids


def _fetch_detailed_record(dtxcids: List[str], identifiers: List[str]) -> Dict[str, Dict[str, Any]]:
    response = requests.post(
        COMPTOX_DTXCID_DETAIL_URL,
        headers=HEADERS,
        json=dtxcids,
        timeout=30,
    )
    response.raise_for_status()
    
    data: Dict[str, Dict[str, Any]] = {}

    records = response.json()
    records_by_dtxcid: Dict[str, Dict[str, Any]] = {}
    for record in records:
        dtxcid = str(record.get("dtxcid", "")).strip()
        if dtxcid:
            records_by_dtxcid[dtxcid] = record

    for dtxcid, identifier in zip(dtxcids, identifiers):
        record = records_by_dtxcid.get(str(dtxcid).strip())
        if not record:
            continue

        info = {}
        info[CAS_COLUMN] = record.get("casrn", "")
        info[SMILES_COLUMN] = record.get("smiles", "")
        info[FORMULA_COLUMN] = record.get("molFormula", "")
        info[ENRICHED_NAME_COLUMNS_COLUMN] = [
            name
            for name in [record.get("preferredName", ""), record.get("iupacName", "")]
            if not _is_missing_value(name)
        ]
        data[identifier] = info
    return data


def fetch_chemical_info(
    identifiers: List[str]
) -> Dict[str, Optional[Union[str, List[str]]]]:
    """Fetch missing chemical information from EPA CompTox using DTXCID detail lookups."""

    if not API_KEY:
        raise ValueError("CompTox API key is not set in environment variable 'COMPTOX_API_KEY'")

    try:
        dtxcids = _fetch_dtxcid_records(identifiers)

        results = {identifier: {
            CAS_COLUMN: None,
            SMILES_COLUMN: None,
            ENRICHED_NAME_COLUMNS_COLUMN: None,
            FORMULA_COLUMN: None,
        } for identifier in identifiers}

        if not dtxcids:
            return results
        
        found_dtxcids = []
        found_identifiers = []
        for i, dtxcid in enumerate(dtxcids):
            if not dtxcid:
                continue
            found_dtxcids.append(dtxcid)
            found_identifiers.append(identifiers[i])

        if not found_dtxcids:
            return results

        detailed_records = _fetch_detailed_record(found_dtxcids, found_identifiers)
        results.update(detailed_records)

    except Exception as e:
        print(f"  [WARN] Could not fetch data for {identifiers}: {str(e)}")

    return results
