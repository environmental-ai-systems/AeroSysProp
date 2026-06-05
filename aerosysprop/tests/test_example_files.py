"""The downloadable example files offered in the app must load cleanly
through the same loaders the upload tab uses."""

import os
import tempfile

from rdkit import Chem

from aerosysprop._ui import (
    _example_smiles_txt, _example_named_txt, _example_sdf, _EXAMPLE_COMPOUNDS,
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


def test_examples_are_valid_smiles():
    for _, smi in _EXAMPLE_COMPOUNDS:
        assert Chem.MolFromSmiles(smi) is not None
