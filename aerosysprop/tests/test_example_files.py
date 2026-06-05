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
from aerosysprop.io import read_smiles_txt, read_smiles_sdf, detect_columns


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


# ---- column auto-detection (detect_columns) -------------------------------

def test_detect_single_smiles_column():
    assert detect_columns(['O', 'CCO', 'CC(=O)O']) == (0, None, 1)


def test_detect_name_then_smiles():
    lines = ['water,O', 'ethanol,CCO', 'acetic_acid,CC(=O)O']
    assert detect_columns(lines) == (1, 0, 2)


def test_detect_smiles_then_name():
    lines = ['O,water', 'CCO,ethanol', 'CC(=O)O,acetic_acid']
    assert detect_columns(lines) == (0, 1, 2)


def test_detect_smiles_in_third_column():
    lines = ['1,water,O', '2,ethanol,CCO', '3,acetic_acid,CC(=O)O']
    s, n, ncols = detect_columns(lines)
    assert s == 2 and ncols == 3


def test_detect_skips_comment_and_blank_lines():
    lines = ['# header comment', '', 'water\tO', 'ethanol\tCCO']
    assert detect_columns(lines) == (1, 0, 2)


def test_detect_empty_input():
    assert detect_columns(['# only a comment', '']) == (0, None, 0)


def test_detect_then_read_roundtrip():
    """The detected columns must actually load the file."""
    content = 'glycerol,OCC(O)CO\nlevoglucosan,OC[C@H]1OC2OC1C(O)C2O\n'
    path = _write(content.encode(), '.csv')
    try:
        s, n, _ = detect_columns(content.splitlines())
        kwargs = {'smiles_col': s}
        if n is not None:
            kwargs['name_col'] = n
        mols = read_smiles_txt(path, **kwargs)
    finally:
        os.unlink(path)
    assert [nm for nm, _ in mols] == ['glycerol', 'levoglucosan']
