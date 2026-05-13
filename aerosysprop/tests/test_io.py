"""
Tests for aerosysprop.io — SMILES loading from text, SQLite, and SDF files.
"""
import os
import sqlite3
import tempfile

import pytest
from rdkit import Chem
from rdkit.Chem import SDWriter

from aerosysprop.io import (
    create_smiles_db,
    read_smiles_db,
    read_smiles_sdf,
    read_smiles_txt,
)

FIXTURES = os.path.join(os.path.dirname(__file__), 'io_fixtures')

# Four compounds used across all fixture files
EXPECTED = {
    'methyl_hydroperoxide': 'COO',
    'ethylene_glycol':      'C(O)CO',
    'acetone':              'CC(C)=O',
    'acetic_acid':          'CC(=O)O',
}
EXPECTED_SMILES = set(EXPECTED.values())


def _canonical(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


def _check_results(results, expected_names=None):
    """Common assertions for all loader outputs."""
    assert len(results) == len(EXPECTED)
    for name, mol in results:
        assert mol is not None, 'mol is None for %r' % name
        assert isinstance(mol, Chem.rdchem.Mol)
    if expected_names:
        assert set(name for name, _ in results) == set(expected_names)


# ---------------------------------------------------------------------------
# read_smiles_txt
# ---------------------------------------------------------------------------

class TestReadSmilesTxt:

    def test_single_column(self):
        """File with one SMILES per line (no name column)."""
        path = os.path.join(FIXTURES, 'single_col.txt')
        results = read_smiles_txt(path, smiles_col=0)
        _check_results(results)
        # when there is no name column, the SMILES is used as the name
        canonical_names = {_canonical(name) for name, _ in results}
        canonical_expected = {_canonical(s) for s in EXPECTED_SMILES}
        assert canonical_names == canonical_expected

    def test_two_column(self):
        """File with name\\tSMILES columns."""
        path = os.path.join(FIXTURES, 'two_col.txt')
        results = read_smiles_txt(path, smiles_col=1, name_col=0)
        _check_results(results, expected_names=EXPECTED.keys())

    def test_comment_lines_skipped(self):
        """Lines starting with # should be ignored."""
        path = os.path.join(FIXTURES, 'single_col.txt')
        results = read_smiles_txt(path, smiles_col=0)
        # single_col.txt has one comment line; should still get 4 compounds
        assert len(results) == 4

    def test_invalid_smiles_skipped(self, tmp_path):
        """Invalid SMILES should be skipped with a warning, not raise."""
        p = tmp_path / 'bad.txt'
        p.write_text('COO\nNOT_A_SMILES\nCC(C)=O\n')
        with pytest.warns(UserWarning, match='invalid SMILES'):
            results = read_smiles_txt(str(p), smiles_col=0)
        assert len(results) == 2

    def test_empty_lines_skipped(self, tmp_path):
        p = tmp_path / 'gaps.txt'
        p.write_text('COO\n\nCC(C)=O\n\n')
        results = read_smiles_txt(str(p), smiles_col=0)
        assert len(results) == 2

    def test_mols_usable_with_methods(self):
        """Molecules returned should work with the property methods."""
        from aerosysprop.methods import boiling_points
        path = os.path.join(FIXTURES, 'single_col.txt')
        results = read_smiles_txt(path, smiles_col=0)
        for _, mol in results:
            bp = boiling_points.nannoolal(mol)
            assert isinstance(bp, float)
            assert bp > 0


# ---------------------------------------------------------------------------
# read_smiles_db
# ---------------------------------------------------------------------------

class TestReadSmilesDb:

    def test_reads_all_compounds(self):
        path = os.path.join(FIXTURES, 'compounds.db')
        results = read_smiles_db(path)
        _check_results(results, expected_names=EXPECTED.keys())

    def test_custom_table_and_columns(self, tmp_path):
        """Loader should respect non-default table/column names."""
        db = str(tmp_path / 'custom.db')
        conn = sqlite3.connect(db)
        conn.execute('CREATE TABLE mols (cid TEXT, smi TEXT NOT NULL)')
        conn.executemany('INSERT INTO mols VALUES (?, ?)',
                         [('A', 'COO'), ('B', 'CC(C)=O')])
        conn.commit()
        conn.close()

        results = read_smiles_db(db, table='mols', smiles_col='smi', name_col='cid')
        assert len(results) == 2
        assert {n for n, _ in results} == {'A', 'B'}

    def test_missing_name_column_falls_back_to_smiles(self, tmp_path):
        """If the name column does not exist, SMILES is used as the name."""
        db = str(tmp_path / 'noname.db')
        conn = sqlite3.connect(db)
        conn.execute('CREATE TABLE compounds (smiles TEXT NOT NULL)')
        conn.execute("INSERT INTO compounds VALUES ('COO')")
        conn.commit()
        conn.close()

        results = read_smiles_db(db, name_col='name')
        assert len(results) == 1
        name, mol = results[0]
        assert _canonical(name) == _canonical('COO')

    def test_invalid_smiles_skipped(self, tmp_path):
        db = str(tmp_path / 'bad.db')
        create_smiles_db(db, [('good', 'COO'), ('bad', 'NOT_SMILES')])
        with pytest.warns(UserWarning, match='invalid SMILES'):
            results = read_smiles_db(db)
        assert len(results) == 1


# ---------------------------------------------------------------------------
# read_smiles_sdf
# ---------------------------------------------------------------------------

class TestReadSmilesSdf:

    def test_reads_all_compounds(self):
        path = os.path.join(FIXTURES, 'compounds.sdf')
        results = read_smiles_sdf(path)
        _check_results(results, expected_names=EXPECTED.keys())

    def test_name_from_sd_property(self):
        """name_property kwarg should pull from an SD tag."""
        path = os.path.join(FIXTURES, 'compounds.sdf')
        results = read_smiles_sdf(path, name_property='compound_id')
        names = {n for n, _ in results}
        assert names == {k.upper() for k in EXPECTED.keys()}

    def test_fallback_name_when_no_title(self, tmp_path):
        """Molecules without a title line should get a sequential name."""
        sdf = str(tmp_path / 'notitle.sdf')
        writer = SDWriter(sdf)
        for smi in ['COO', 'CC(C)=O']:
            mol = Chem.MolFromSmiles(smi)
            writer.write(mol)
        writer.close()

        results = read_smiles_sdf(sdf)
        assert len(results) == 2
        assert {n for n, _ in results} == {'mol_1', 'mol_2'}

    def test_invalid_record_skipped(self, tmp_path, recwarn):
        """A corrupt SDF record should be skipped with a warning."""
        sdf = str(tmp_path / 'corrupt.sdf')
        # write one valid molecule then append a truncated record
        writer = SDWriter(sdf)
        writer.write(Chem.MolFromSmiles('COO'))
        writer.close()
        with open(sdf, 'a') as f:
            f.write('\n$$$$\n')   # empty record — RDKit returns None
        results = read_smiles_sdf(sdf)
        assert len(results) == 1


# ---------------------------------------------------------------------------
# create_smiles_db (round-trip)
# ---------------------------------------------------------------------------

class TestCreateSmilesDb:

    def test_roundtrip(self, tmp_path):
        """create_smiles_db then read_smiles_db should return the same data."""
        db = str(tmp_path / 'rt.db')
        compounds = [('ethanol', 'CCO'), ('propane', 'CCC')]
        create_smiles_db(db, compounds)
        results = read_smiles_db(db)
        assert len(results) == 2
        assert {n for n, _ in results} == {'ethanol', 'propane'}

    def test_append_to_existing(self, tmp_path):
        """Calling create_smiles_db twice should append, not overwrite."""
        db = str(tmp_path / 'append.db')
        create_smiles_db(db, [('ethanol', 'CCO')])
        create_smiles_db(db, [('propane', 'CCC')])
        results = read_smiles_db(db)
        assert len(results) == 2
