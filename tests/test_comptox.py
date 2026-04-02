"""Unit tests for CompTox helpers."""

from fccgroup.comptox import _fetch_detailed_record
from fccgroup.constants import (
    CAS_COLUMN,
    ENRICHED_NAME_COLUMNS_COLUMN,
    FORMULA_COLUMN,
    SMILES_COLUMN,
)


class _DummyResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def test_fetch_detailed_record_reorders_using_dtxcid(monkeypatch):
    """Detail records should map to identifiers by DTXCID, not API response order."""

    # Response order is intentionally swapped relative to request dtxcids.
    payload = [
        {
            "dtxcid": "DTXCID2",
            "casrn": "222-22-2",
            "smiles": "CC",
            "molFormula": "C2H6",
            "preferredName": "Ethane",
            "iupacName": "ethane",
        },
        {
            "dtxcid": "DTXCID1",
            "casrn": "111-11-1",
            "smiles": "C=O",
            "molFormula": "CH2O",
            "preferredName": "Formaldehyde",
            "iupacName": "methanal",
        },
    ]

    def _fake_post(url, headers, json, timeout):
        return _DummyResponse(payload)

    monkeypatch.setattr("fccgroup.comptox.requests.post", _fake_post)

    result = _fetch_detailed_record(
        dtxcids=["DTXCID1", "DTXCID2"],
        identifiers=["id-A", "id-B"],
    )

    assert result["id-A"][CAS_COLUMN] == "111-11-1"
    assert result["id-A"][SMILES_COLUMN] == "C=O"
    assert result["id-A"][FORMULA_COLUMN] == "CH2O"
    assert result["id-A"][ENRICHED_NAME_COLUMNS_COLUMN] == ["Formaldehyde", "methanal"]

    assert result["id-B"][CAS_COLUMN] == "222-22-2"
    assert result["id-B"][SMILES_COLUMN] == "CC"
    assert result["id-B"][FORMULA_COLUMN] == "C2H6"
    assert result["id-B"][ENRICHED_NAME_COLUMNS_COLUMN] == ["Ethane", "ethane"]
