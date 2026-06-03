"""Tests for the aerosysprop.fingerprints module."""

import numpy as np
import pandas as pd
import pytest

from aerosysprop.fingerprints import (
    smarts_hits, smarts_hits_all, available_methods, _SMARTS_PROVIDERS,
)


def test_available_methods_returns_known_set():
    methods = available_methods()
    # core ones must be present
    assert 'aiomfac' in methods
    assert 'simpol' in methods
    assert 'nannoolal_primary' in methods
    assert 'evaporation' in methods


def test_aiomfac_hits_basic():
    """A minimal three-SMILES batch should produce sensible AIOMFAC hits."""
    df = smarts_hits(
        ['O', 'CCO', 'Cc1ccccc1'],
        method='aiomfac', drop_zero_columns=True,
    )
    assert isinstance(df, pd.DataFrame)
    assert list(df.index) == ['O', 'CCO', 'Cc1ccccc1']
    # Water should match subgroup 16 (H2O) exactly once
    assert df.loc['O', 16] == 1
    # Ethanol should match subgroup 14 (OH), 1 (CH3), 2 (CH2) once each
    assert df.loc['CCO', 14] == 1
    assert df.loc['CCO', 1] == 1
    assert df.loc['CCO', 2] == 1
    # Toluene should match ACH (9) five times and ACCH3 (11) once
    assert df.loc['Cc1ccccc1', 9] == 5
    assert df.loc['Cc1ccccc1', 11] == 1


def test_simpol_hits():
    df = smarts_hits(['O', 'CCO'], method='simpol', drop_zero_columns=True)
    assert isinstance(df, pd.DataFrame)
    # ethanol should have some non-zero count
    assert (df.loc['CCO'] > 0).any()


def test_invalid_smiles_yields_zero_row():
    df = smarts_hits(['not-a-smiles', 'CCO'], method='aiomfac',
                     drop_zero_columns=True)
    # The invalid row should be all zeros
    assert (df.loc['not-a-smiles'] == 0).all()


def test_drop_zero_columns_actually_drops():
    df_full = smarts_hits(['O'], method='aiomfac', drop_zero_columns=False)
    df_drop = smarts_hits(['O'], method='aiomfac', drop_zero_columns=True)
    assert df_drop.shape[1] < df_full.shape[1]
    # And the dropped frame has only non-zero columns
    assert (df_drop.values.sum(axis=0) > 0).all()


def test_as_dataframe_false_returns_tuple():
    res = smarts_hits(['O', 'CCO'], method='aiomfac', as_dataframe=False)
    assert isinstance(res, tuple) and len(res) == 3
    matrix, rows, cols = res
    assert isinstance(matrix, np.ndarray)
    assert matrix.shape[0] == 2
    assert rows == ['O', 'CCO']


def test_unknown_method_raises():
    with pytest.raises(ValueError, match='Unknown method'):
        smarts_hits(['O'], method='not-a-method')


def test_smarts_hits_all_multiindex_columns():
    df = smarts_hits_all(
        ['O', 'CCO'],
        methods=['aiomfac', 'simpol'],
        drop_zero_columns=True,
    )
    # Columns should be a MultiIndex with method + id levels
    assert isinstance(df.columns, pd.MultiIndex)
    assert df.columns.names == ['method', 'id']
    assert set(df.columns.get_level_values('method')) == {'aiomfac', 'simpol'}


def test_smarts_hits_all_defaults_to_all_methods():
    df = smarts_hits_all(['O'], drop_zero_columns=False)
    # Every registered method should appear
    methods = set(df.columns.get_level_values('method'))
    assert 'aiomfac' in methods
    assert 'simpol' in methods


def test_column_order_is_sorted_ids():
    df = smarts_hits(['O', 'CCO'], method='aiomfac', drop_zero_columns=False)
    cols = list(df.columns)
    assert cols == sorted(cols), 'columns should be in ascending ID order'


def test_index_label_customisable():
    df = smarts_hits(['O', 'CCO'], method='aiomfac',
                     index_label='my_label')
    assert df.index.name == 'my_label'


def test_columns_name_is_method():
    df = smarts_hits(['O'], method='simpol')
    assert df.columns.name == 'simpol'
