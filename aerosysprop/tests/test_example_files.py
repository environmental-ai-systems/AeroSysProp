"""The downloadable example files offered in the app must load cleanly
through the same loaders the upload tab uses."""

import os
import tempfile

from rdkit import Chem

import pytest
from aerosysprop._ui import (
    _example_smiles_txt, _example_named_txt, _example_named_csv,
    _example_sdf, _EXAMPLE_COMPOUNDS,
)
from aerosysprop.io import read_smiles_txt, read_smiles_sdf


def _write(data, suffix):
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.write(fd, data)
    os.close(fd)
    return path


def test_smiles_list_loads_with_default_settings():
    path = _write(_example_smiles_txt(), '.txt')
    try:
        mols = read_smiles_txt(path, smiles_col=0)   # app upload default
    finally:
        os.unlink(path)
    assert len(mols) == len(_EXAMPLE_COMPOUNDS)
    assert all(m is not None for _, m in mols)


def test_named_list_loads_with_documented_columns():
    path = _write(_example_named_txt(), '.txt')
    try:
        mols = read_smiles_txt(path, smiles_col=1, name_col=0)
    finally:
        os.unlink(path)
    names = [n for n, _ in mols]
    assert names == [n for n, _ in _EXAMPLE_COMPOUNDS]
    assert all(m is not None for _, m in mols)


def test_sdf_loads_and_preserves_names():
    path = _write(_example_sdf(), '.sdf')
    try:
        mols = read_smiles_sdf(path)
    finally:
        os.unlink(path)
    names = [n for n, _ in mols]
    assert names == [n for n, _ in _EXAMPLE_COMPOUNDS]
    assert all(m is not None for _, m in mols)


def test_named_csv_loads():
    path = _write(_example_named_csv(), '.csv')
    try:
        mols = read_smiles_txt(path, smiles_col=1, name_col=0)
    finally:
        os.unlink(path)
    names = [n for n, _ in mols]
    assert names == [n for n, _ in _EXAMPLE_COMPOUNDS]


def test_examples_are_valid_smiles():
    for _, smi in _EXAMPLE_COMPOUNDS:
        assert Chem.MolFromSmiles(smi) is not None


# ---- delimiter auto-detection in read_smiles_txt --------------------------

@pytest.mark.parametrize('sep', ['\t', ',', ' '])
def test_auto_detect_separator(sep):
    content = f'water{sep}O\nethanol{sep}CCO\n'
    path = _write(content.encode(), '.txt')
    try:
        mols = read_smiles_txt(path, smiles_col=1, name_col=0)
    finally:
        os.unlink(path)
    assert [n for n, _ in mols] == ['water', 'ethanol']
    assert Chem.MolToSmiles(mols[1][1]) == 'CCO'


def test_comma_with_trailing_space_is_stripped():
    path = _write(b'water, O\nethanol, CCO\n', '.csv')
    try:
        mols = read_smiles_txt(path, smiles_col=1, name_col=0)
    finally:
        os.unlink(path)
    assert Chem.MolToSmiles(mols[0][1]) == 'O'
    assert Chem.MolToSmiles(mols[1][1]) == 'CCO'


def test_single_column_still_works():
    path = _write(b'CCO\nCC(=O)O\n', '.txt')
    try:
        mols = read_smiles_txt(path, smiles_col=0)
    finally:
        os.unlink(path)
    assert len(mols) == 2


def test_explicit_delimiter_overrides_autodetect():
    path = _write(b'name|CCO\n', '.txt')
    try:
        mols = read_smiles_txt(path, smiles_col=1, name_col=0, delimiter='|')
    finally:
        os.unlink(path)
    assert Chem.MolToSmiles(mols[0][1]) == 'CCO'
