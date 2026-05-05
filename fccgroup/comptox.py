"""CompTox API integration helpers."""

import time
from typing import Any, Dict, List, Optional, Union
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

_RATE_LIMIT_INTERVAL = 1.0 / 5  # 5 requests per second
_last_request_time: float = 0.0


def _is_missing_value(value: Any) -> bool:
    if value is None:
        return True
    return str(value).strip() == ""


def _rate_limit() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < _RATE_LIMIT_INTERVAL:
        time.sleep(_RATE_LIMIT_INTERVAL - elapsed)
    _last_request_time = time.monotonic()


def _post_with_retry(url: str, max_retries: int = 1, **kwargs) -> requests.Response:
    for attempt in range(max_retries + 1):
        _rate_limit()
        try:
            response = requests.post(url, headers=HEADERS, **kwargs)
            response.raise_for_status()
            return response
        except Exception:
            if attempt >= max_retries:
                raise
            time.sleep(2.0)
    raise RuntimeError("unreachable")


def _fetch_dtxcid_records(identifiers: List[str]) -> List[str]:
    response = _post_with_retry(
        COMPTOX_SEARCH_EQUAL_URL,
        data="\n".join(identifiers).strip(),
        timeout=30,
    )
    dtxcids = []
    for entry in response.json():
        dtxcids.append(entry["dtxcid"])
    return dtxcids


def _fetch_detailed_record(dtxcids: List[str], identifiers: List[str]) -> Dict[str, Dict[str, Any]]:
    response = _post_with_retry(
        COMPTOX_DTXCID_DETAIL_URL,
        json=dtxcids,
        timeout=30,
    )

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

    identifiers = [id for id in identifiers if not _is_missing_value(id)]
    if not identifiers:
        return {}

    results: Dict[str, Any] = {
        identifier: {
            CAS_COLUMN: None,
            SMILES_COLUMN: None,
            ENRICHED_NAME_COLUMNS_COLUMN: None,
            FORMULA_COLUMN: None,
        }
        for identifier in identifiers
    }

    try:
        dtxcids = _fetch_dtxcid_records(identifiers)

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
