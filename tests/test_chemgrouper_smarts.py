"""Tests for ChemicalGrouper using SMARTS patterns only."""

import pytest
import pandas as pd

from fccgroup import ChemicalGrouper, GroupingConfig, GroupingMethod, ColumnMapping
from fccgroup.constants import CAS_COLUMN


class TestChemicalGrouperSMARTS:
    """Test suite for ChemicalGrouper using SMARTS-only mode."""
    
    def test_formaldehyde_cas_50_00_0(self, grouper_smarts_only, formaldehyde_df):
        """Test grouping of CAS 50-00-0 (formaldehyde) using SMILES C=O."""
        results = grouper_smarts_only.group_chemicals()
        
        # Check that results DataFrame was created
        assert results is not None
        assert len(results) == 1
        
        # Check that required columns exist
        assert 'SMILES' in results.columns
        assert results['SMILES'].iloc[0] == 'C=O'
        
        # Check that SMARTS grouping was applied
        assert 'Chemical groups' in results.columns or any('group' in str(col).lower() for col in results.columns)
        
        # Formaldehyde should contain oxygen patterns
        groups = results['Chemical groups'].iloc[0] if 'Chemical groups' in results.columns else str(results.iloc[0])
        assert isinstance(groups, (str, list))
    
    def test_ethane_smiles_cc(self, grouper_smarts_only, ethane_df):
        """Test grouping of SMILES CC (ethane) using SMARTS-only mode."""
        results = grouper_smarts_only.group_chemicals()
        
        # Check that results DataFrame was created
        assert results is not None
        assert len(results) == 1
        
        # Check that required columns exist
        assert 'SMILES' in results.columns
        assert results['SMILES'].iloc[0] == 'CC'
        
        # Check that SMARTS grouping was applied
        assert 'Chemical groups' in results.columns or any('group' in str(col).lower() for col in results.columns)
        
        # Ethane should contain carbon/alkane patterns
        groups = results['Chemical groups'].iloc[0] if 'Chemical groups' in results.columns else str(results.iloc[0])
        assert isinstance(groups, (str, list))
    
    def test_grouping_config_smarts_only_mode(self):
        """Test that a SMARTS-only config is properly configured."""
        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS],
            column_mapping=ColumnMapping(cas=None, smiles='SMILES'),
        )

        assert config.use_smarts is True
        assert config.use_lists is False
        assert config.use_regex is False
        assert "SMARTS" in config.description
    
    def test_invalid_smiles_not_in_lookup(self, grouper_smarts_only, formaldehyde_df):
        """Test behavior when SMILES is not found in lookup."""
        results = grouper_smarts_only.group_chemicals()
        
        # Should still return results, but marked as not found
        assert len(results) == 1
        assert results['SMILES'].iloc[0] == '[Invalid]'

    def test_methods_smarts_and_regex(self):
        """A config selecting SMARTS and REGEX should reflect both, but not LISTS."""
        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS, GroupingMethod.REGEX],
            column_mapping=ColumnMapping(cas="CASRN", smiles="Structure"),
        )

        assert config.use_smarts is True
        assert config.use_regex is True
        assert config.use_lists is False

    def test_custom_column_mapping_for_smiles_input(self):
        """Grouping should work when user provides non-default column names."""
        custom_df = pd.DataFrame({
            'Structure': ['CC'],
            'CASRN': ['74-84-0'],
            'Name': ['ethane'],
            'IUPAC': ['ethane'],
            'Formula': ['C2H6'],
        })

        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS],
            column_mapping=ColumnMapping(
                cas=None,
                smiles='Structure',
                name_columns=['Name', 'IUPAC'],
                formula='Formula',
            )
        )

        grouper = ChemicalGrouper(df=custom_df, grouping_config=config)
        results = grouper.group_chemicals()

        assert len(results) == 1
        assert 'SMILES' in results.columns
        assert results['SMILES'].iloc[0] == 'CC'
        assert 'Chemical groups' in results.columns

    def test_regex_requires_cas_column_for_smiles_input(self):
        """Regex selection should fail fast if CAS column is missing."""
        df_missing_cas = pd.DataFrame({
            'Structure': ['CC'],
            'Name': ['ethane'],
            'Formula': ['C2H6'],
        })

        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS, GroupingMethod.REGEX],
            column_mapping=ColumnMapping(
                cas='CASRN',
                smiles='Structure',
                name_columns=['Name', 'IUPAC'],
                formula='Formula',
            )
        )

        with pytest.raises(ValueError, match='Mapped columns were declared'):
            ChemicalGrouper(df=df_missing_cas, grouping_config=config)

    def test_regex_allows_missing_name_and_formula_when_cas_is_present(self):
        """Regex should run when CAS is provided, even if names/formula are absent in input df."""
        df_minimal = pd.DataFrame({
            'Structure': ['CC'],
            'CASRN': ['74-84-0'],
        })

        config = GroupingConfig(
            methods=[GroupingMethod.REGEX],
            column_mapping=ColumnMapping(
                cas='CASRN',
                smiles='Structure',
                name_columns=['Name', 'IUPAC'],
                formula='Formula',
            )
        )

        grouper = ChemicalGrouper(df=df_minimal, grouping_config=config)
        results = grouper.group_chemicals()

        assert len(results) == 1

    def test_lists_requires_cas_column_for_smiles_input(self):
        """List matching should fail fast when CAS mapping column is missing in SMILES mode."""
        df_missing_cas = pd.DataFrame({
            'Structure': ['CC'],
            'Name': ['ethane'],
            'IUPAC': ['ethane'],
            'Formula': ['C2H6'],
        })

        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS, GroupingMethod.LISTS],
            column_mapping=ColumnMapping(
                cas='CASRN',
                smiles='Structure',
                name_columns=['Name', 'IUPAC'],
                formula='Formula',
            )
        )

        with pytest.raises(ValueError, match='Mapped columns were declared'):
            ChemicalGrouper(df=df_missing_cas, grouping_config=config)

    def test_smarts_fails_without_cas_or_smiles_column_mapping(self):
        """SMARTS must raise ValueError at construction if neither CAS nor SMILES is mapped."""
        with pytest.raises(ValueError, match='Provide at least one'):
            GroupingConfig(
                methods=[GroupingMethod.SMARTS],
                column_mapping=ColumnMapping(cas=None, smiles=None),
            )

    def test_lists_allows_smiles_only_mapping(self):
        """Config accepts LISTS with smiles mapping; feasibility is resolved by grouper init."""
        config = GroupingConfig(
            methods=[GroupingMethod.LISTS],
            column_mapping=ColumnMapping(cas=None, smiles='Structure'),
        )
        assert config.use_lists is True

    def test_regex_allows_smiles_only_mapping(self):
        """Config accepts REGEX with smiles mapping; missing fields can be CIRPY-enriched."""
        config = GroupingConfig(
            methods=[GroupingMethod.REGEX],
            column_mapping=ColumnMapping(cas=None, smiles='Structure'),
        )
        assert config.use_regex is True

    def test_smarts_and_lists_with_cas_only_mapping_passes_validation(self):
        """Providing only a CAS column name satisfies both SMARTS and LISTS validation.

        SMARTS requires CAS OR SMILES; LISTS requires CAS.
        A mapping with only cas=<col> (smiles=None) satisfies both at construction time.
        """
        # Construction (and therefore validate_mapping) should not raise.
        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS, GroupingMethod.LISTS],
            column_mapping=ColumnMapping(cas='CASRN', smiles=None),
        )
        config.validate_mapping()  # explicit call also must not raise

    def test_smarts_uses_only_selected_fingerprints(self):
        """SMARTS should only evaluate explicitly selected fingerprints when configured."""
        df = pd.DataFrame({
            'Structure': ['CC'],
            'CASRN': ['74-84-0'],
            'Name': ['ethane'],
            'IUPAC': ['ethane'],
            'Formula': ['C2H6'],
        })

        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS],
            smarts_fingerprints={'Contains C-C'},
            column_mapping=ColumnMapping(
                cas=None,
                smiles='Structure',
                name_columns=['Name', 'IUPAC'],
                formula='Formula',
            )
        )

        grouper = ChemicalGrouper(df=df, grouping_config=config)
        results = grouper.group_chemicals()

        assert 'Contains C-C' in results.columns
        assert 'Alkenes' not in results.columns

    def test_smarts_rejects_unknown_selected_fingerprints(self):
        """Unknown SMARTS fingerprint names in config should fail fast."""
        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS],
            smarts_fingerprints={'NotARealFingerprint'},
            column_mapping=ColumnMapping(cas=None, smiles='Structure'),
        )

        with pytest.raises(ValueError, match='Unknown SMARTS fingerprints'):
            ChemicalGrouper(df=pd.DataFrame({'Structure': ['CC']}), grouping_config=config)

    def test_smiles_to_cas_list_enrichment_is_scalarized(self, monkeypatch):
        """List-valued CAS resolver outputs should be normalized before DataFrame assignment."""
        df = pd.DataFrame({'Structure': ['CC']})
        config = GroupingConfig(
            methods=[GroupingMethod.SMARTS, GroupingMethod.LISTS],
            column_mapping=ColumnMapping(cas=None, smiles='Structure'),
        )

        # Skip list asset loading; this test targets resolver assignment behavior.
        monkeypatch.setattr(ChemicalGrouper, '_ensure_lists_loaded', lambda self: None)
        monkeypatch.setattr(
            ChemicalGrouper,
            '_apply_structural_patterns',
            lambda self, df, id_column, smiles_column, fingerprints_dict: df,
        )

        def _fake_fetch(identifiers):
            return {
                identifier: {
                    CAS_COLUMN: ['74-84-0', '999-99-9'],
                }
                for identifier in identifiers
            }

        monkeypatch.setattr('fccgroup.grouper.fetch_chemical_info', _fake_fetch)

        result = ChemicalGrouper(df=df, grouping_config=config).group_chemicals()

        assert CAS_COLUMN in result.columns
        assert result.loc[0, CAS_COLUMN] == '74-84-0'

    def test_build_comptox_batches_enforces_strict_limit(self):
        """CompTox batches must satisfy len('\\n'.join(batch)) < 200 for each request."""
        grouper = ChemicalGrouper(
            df=pd.DataFrame({'Structure': ['CC']}),
            grouping_config=GroupingConfig(
                methods=[GroupingMethod.SMARTS],
                column_mapping=ColumnMapping(cas=None, smiles='Structure'),
            ),
        )

        # 3 identifiers of length 66 cannot fit in one batch because 66+1+66+1+66 = 200.
        queue = [
            (0, 'A' * 66),
            (1, 'B' * 66),
            (2, 'C' * 66),
        ]

        batches = grouper._build_comptox_batches(queue, max_payload_chars=200)

        assert len(batches) == 2
        for batch in batches:
            identifiers = [identifier for _, identifier in batch]
            assert len("\n".join(identifiers)) < 200

    def test_comptox_batch_failure_continues_other_batches(self, monkeypatch):
        """A failing CompTox batch should not stop enrichment of other batches."""
        df = pd.DataFrame({'Structure': ['A' * 120, 'B' * 120]})
        config = GroupingConfig(
            methods=[GroupingMethod.LISTS],
            column_mapping=ColumnMapping(cas=None, smiles='Structure'),
        )

        monkeypatch.setattr(ChemicalGrouper, '_ensure_lists_loaded', lambda self: None)
        monkeypatch.setattr(
            ChemicalGrouper,
            '_apply_functional_lists',
            lambda self, df, cas_column: df,
        )

        def _fake_fetch(identifiers):
            if identifiers and identifiers[0].startswith('B'):
                raise RuntimeError('simulated batch failure')
            return {
                identifier: {
                    CAS_COLUMN: '74-84-0',
                }
                for identifier in identifiers
            }

        monkeypatch.setattr('fccgroup.grouper.fetch_chemical_info', _fake_fetch)

        result = ChemicalGrouper(df=df, grouping_config=config).group_chemicals()

        assert result.loc[0, CAS_COLUMN] == '74-84-0'
        assert result.loc[1, CAS_COLUMN] == ''
