"""Tests for the SMARTS fingerprint dependency graph and subset selection."""

import pytest
import pandas as pd
from rdkit import Chem

from fccgroup import ChemicalGrouper, GroupingConfig, GroupingMethod, ColumnMapping
from fccgroup.constants import MULTIINDEX_STRUCTURAL_LABEL, OUTPUT_COLUMN
from fccgroup.molecular.composition import molecule_composition
from fccgroup.patterns.library import fingerprints
from fccgroup.patterns.methods import generate_fingerprints
from fccgroup.patterns._dependency_graph import FINGERPRINT_DEPENDENCIES
from fccgroup.patterns.dependencies import (
    dependency_order,
    extract_dependencies,
    find_cycle,
    render_graph_module,
    resolve_dependencies,
)

SAMPLE_SMILES = ["CC", "CC(C)C", "C=CC", "C#C", "c1ccccc1", "Oc1ccccc1", "C=O", "FC(F)(Cl)Cl",
                 "[Na+].[Cl-]", "CC(C)(C)c1ccc(O)cc1"]


def _structural(smiles_list, **config_kwargs):
    df = pd.DataFrame({"Structure": smiles_list})
    config = GroupingConfig(
        methods=[GroupingMethod.SMARTS],
        column_mapping=ColumnMapping(cas=None, smiles="Structure"),
        **config_kwargs,
    )
    return ChemicalGrouper(df=df, grouping_config=config).group_chemicals(save=False)


class TestDependencyGraph:
    def test_committed_graph_is_up_to_date(self):
        """Regenerate with scripts/build_fingerprint_graph.py if this fails."""
        from fccgroup.patterns import _dependency_graph
        import inspect

        graph = extract_dependencies(fingerprints)
        assert inspect.getsource(_dependency_graph).replace("\r\n", "\n") == render_graph_module(graph)

    def test_graph_has_no_cycles(self):
        assert find_cycle(extract_dependencies(fingerprints)) == []

    def test_known_dependencies_are_detected(self):
        assert "Aliphatic hydrocarbons" in FINGERPRINT_DEPENDENCIES["Alkanes"]
        assert "Alkanes" in FINGERPRINT_DEPENDENCIES["Branched alkanes"]

    def test_closure_is_transitive(self):
        closure = resolve_dependencies({"Branched alkanes"}, FINGERPRINT_DEPENDENCIES)
        assert {"Alkanes", "Aliphatic hydrocarbons", "Hydrocarbons"} <= closure

    def test_multiline_lambda_is_detected(self):
        fps = {
            "A": "C",
            "B": lambda x, row: row["A"] or
                row["A"],
        }
        assert extract_dependencies(fps)["B"] == {"A"}

    def test_nested_comprehension_is_detected(self):
        fps = {"A": "C", "B": lambda x, row: any(row[k] for k in ["A"])}
        assert extract_dependencies(fps)["B"] == {"A"}

    def test_dependency_order_puts_dependencies_first(self):
        ordered = dependency_order(set(fingerprints), FINGERPRINT_DEPENDENCIES, fingerprints)
        position = {name: i for i, name in enumerate(ordered)}
        assert len(ordered) == len(fingerprints)
        for name, deps in FINGERPRINT_DEPENDENCIES.items():
            assert all(position[d] < position[name] for d in deps), name

    def test_dependency_order_keeps_definition_order_when_valid(self):
        closure = resolve_dependencies({"Alkanes"}, FINGERPRINT_DEPENDENCIES)
        ordered = dependency_order(closure, FINGERPRINT_DEPENDENCIES, fingerprints)
        assert ordered == [n for n in fingerprints if n in closure]

    def test_static_graph_covers_runtime_lookups(self):
        """Every row[...] key actually read at runtime must be in the static graph."""

        class RecordingRow(dict):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.read = set()

            def __getitem__(self, key):
                self.read.add(key)
                return super().__getitem__(key) if key in self else 0

        for smiles in SAMPLE_SMILES:
            mol = Chem.MolFromSmiles(smiles)
            for name, pattern in fingerprints.items():
                if isinstance(pattern, str):
                    continue
                row = RecordingRow({"Molecular Composition": molecule_composition(mol)})
                try:
                    if pattern.__code__.co_argcount == 1:
                        continue
                    pattern(mol, row)
                except Exception:
                    pass
                missed = (row.read & set(fingerprints)) - set(FINGERPRINT_DEPENDENCIES.get(name, ()))
                assert not missed, f"{name} reads {missed} but the graph does not list it"


REFERENCE_SMILES = ["CC(C)C", "C=CC", "Oc1ccccc1", "FC(F)(Cl)Cl"]


@pytest.fixture(scope="module")
def full_reference():
    """Full fingerprint run (dependency order) for each reference molecule, computed once."""
    ordered = {
        n: fingerprints[n]
        for n in dependency_order(set(fingerprints), FINGERPRINT_DEPENDENCIES, fingerprints)
    }
    return {s: generate_fingerprints("id", s, True, ordered) for s in REFERENCE_SMILES}


class TestSubsetEvaluation:
    @pytest.mark.parametrize("name", sorted(fingerprints))
    def test_every_fingerprint_can_be_computed_alone(self, name, full_reference):
        """Selecting one fingerprint plus its dependencies reproduces the full result."""
        required = resolve_dependencies({name}, FINGERPRINT_DEPENDENCIES)
        subset = {n: fingerprints[n] for n in dependency_order(required, FINGERPRINT_DEPENDENCIES, fingerprints)}

        for smiles in REFERENCE_SMILES:
            got = generate_fingerprints("id", smiles, True, subset)[name]
            assert got == full_reference[smiles][name], f"{name} differs for {smiles}"

    def test_dependencies_are_added_automatically(self):
        results = _structural(["CC", "CC(C)C", "C=CC"], smarts_fingerprints={"Branched alkanes"})
        col = lambda n: (MULTIINDEX_STRUCTURAL_LABEL, n)

        assert list(results[col("Branched alkanes")]) == [False, True, False]
        assert list(results[col(OUTPUT_COLUMN)]) == ["", "Branched alkanes", ""]

    def test_dependency_columns_hidden_by_default(self):
        results = _structural(["CC(C)C"], smarts_fingerprints={"Branched alkanes"})
        columns = {c[1] for c in results.columns}

        assert "Branched alkanes" in columns
        assert not {"Alkanes", "Aliphatic hydrocarbons", "Hydrocarbons"} & columns

    def test_dependency_columns_kept_when_requested(self):
        results = _structural(
            ["CC(C)C"],
            smarts_fingerprints={"Branched alkanes"},
            include_fingerprint_dependencies=True,
        )
        columns = {c[1] for c in results.columns}

        assert {"Branched alkanes", "Alkanes", "Aliphatic hydrocarbons", "Hydrocarbons"} <= columns

    def test_subset_matches_full_run(self):
        smiles = ["CC", "CC(C)C", "C=CC", "c1ccccc1", "Oc1ccccc1"]
        selected = {"Alkanes", "Alkenes", "Aromatic hydrocarbons"}
        subset = _structural(smiles, smarts_fingerprints=selected)
        full = _structural(smiles)

        for name in selected:
            col = (MULTIINDEX_STRUCTURAL_LABEL, name)
            assert list(subset[col]) == list(full[col]), name
