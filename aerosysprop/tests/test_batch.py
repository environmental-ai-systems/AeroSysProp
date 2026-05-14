"""
Tests for aerosysprop.batch — parallel property prediction.
"""
import csv
import os
import tempfile

import pytest
from rdkit import Chem

from aerosysprop.batch import ALL_METHODS, _OUTPUT_KEYS, resolve_methods, run_all
from aerosysprop import db

SMILES = [
    ('methyl_hydroperoxide', 'COO'),
    ('ethylene_glycol',      'C(O)CO'),
    ('acetone',              'CC(C)=O'),
    ('acetic_acid',          'CC(=O)O'),
]


def _mols():
    return [(name, Chem.MolFromSmiles(smi)) for name, smi in SMILES]


# ---------------------------------------------------------------------------
# resolve_methods
# ---------------------------------------------------------------------------

class TestResolveMethods:

    def test_no_deps_unchanged(self):
        assert resolve_methods(['vp_simpol']) == {'vp_simpol'}

    def test_direct_dep_added(self):
        resolved = resolve_methods(['vp_nannoolal'])
        assert 'bp_nannoolal' in resolved

    def test_transitive_deps(self):
        resolved = resolve_methods(['density_schroeder'])
        assert 'cp_nannoolal' in resolved
        assert 'bp_nannoolal' in resolved

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match='Unknown method IDs'):
            run_all(_mols(), methods=['not_a_method'], workers=1)

    def test_invalid_bp_source_raises(self):
        with pytest.raises(ValueError, match='bp_source must be one of'):
            resolve_methods(['vp_nannoolal'], bp_source='not_a_bp')

    def test_bp_source_joback_resolves_correctly(self):
        resolved = resolve_methods(['vp_nannoolal'], bp_source='bp_joback')
        assert 'bp_joback' in resolved
        assert 'bp_nannoolal' not in resolved


# ---------------------------------------------------------------------------
# run_all — method selection
# ---------------------------------------------------------------------------

class TestMethodSelection:

    def test_only_requested_keys_returned(self):
        results = run_all(_mols(), methods=['bp_nannoolal'], workers=1, use_threads=True)
        for r in results:
            assert 'bp_nannoolal_K' in r
            assert 'bp_joback_K' not in r
            assert 'vp_simpol_log10atm' not in r

    def test_dependency_computed_silently(self):
        """vp_nannoolal needs bp_nannoolal — both should appear in results."""
        results = run_all(_mols(), methods=['vp_nannoolal'], workers=1, use_threads=True)
        for r in results:
            assert 'vp_nannoolal_log10atm' in r
            assert 'bp_nannoolal_K' in r

    def test_all_methods_default(self):
        results = run_all(_mols(), methods=None, workers=1, use_threads=True)
        all_keys = {col for m in ALL_METHODS for col in _OUTPUT_KEYS[m]}
        for r in results:
            assert all_keys <= r.keys()

    def test_multiple_methods(self):
        results = run_all(_mols(), methods=['bp_nannoolal', 'vp_simpol'],
                          workers=1, use_threads=True)
        for r in results:
            assert 'bp_nannoolal_K' in r
            assert 'vp_simpol_log10atm' in r
            assert 'bp_joback_K' not in r


# ---------------------------------------------------------------------------
# run_all — general behaviour
# ---------------------------------------------------------------------------

class TestRunAll:

    def test_returns_one_result_per_molecule(self):
        results = run_all(_mols(), workers=2, use_threads=True)
        assert len(results) == len(SMILES)

    def test_result_order_matches_input(self):
        mols = _mols()
        results = run_all(mols, workers=2, use_threads=True)
        for (name, _), result in zip(mols, results):
            assert result['name'] == name

    def test_numeric_predictions_are_floats(self):
        results = run_all(_mols(), methods=['bp_nannoolal', 'vp_simpol'],
                          workers=1, use_threads=True)
        for r in results:
            assert isinstance(r['bp_nannoolal_K'], float)
            assert r['bp_nannoolal_K'] > 0

    def test_single_worker_matches_multi(self):
        mols = _mols()
        r1 = run_all(mols, methods=['bp_nannoolal'], workers=1, use_threads=True)
        r2 = run_all(mols, methods=['bp_nannoolal'], workers=2, use_threads=True)
        for a, b in zip(r1, r2):
            assert a['bp_nannoolal_K'] == pytest.approx(b['bp_nannoolal_K'])

    def test_bp_source_affects_results(self):
        """Different bp_source choices produce different VP values."""
        mols = _mols()
        r_nannoolal = run_all(mols, methods=['vp_nannoolal'], workers=1,
                              use_threads=True, bp_source='bp_nannoolal')
        r_joback = run_all(mols, methods=['vp_nannoolal'], workers=1,
                           use_threads=True, bp_source='bp_joback')
        vp_n = [r['vp_nannoolal_log10atm'] for r in r_nannoolal if r]
        vp_j = [r['vp_nannoolal_log10atm'] for r in r_joback if r]
        assert vp_n != vp_j, 'Different BP inputs should yield different VP outputs'

    def test_progress_callback_called(self):
        calls = []
        run_all(_mols(), methods=['bp_nannoolal'], workers=1,
                use_threads=True, progress_callback=lambda d, t: calls.append((d, t)))
        assert len(calls) == len(SMILES)
        assert calls[-1] == (len(SMILES), len(SMILES))

    def test_csv_output(self, tmp_path):
        csv_path = str(tmp_path / 'out.csv')
        run_all(_mols(), methods=['bp_nannoolal', 'vp_simpol'],
                workers=1, use_threads=True, output_csv=csv_path)
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == len(SMILES)
        assert 'bp_nannoolal_K' in rows[0]
        assert 'bp_joback_K' not in rows[0]

    def test_invalid_smiles_returns_none(self):
        from aerosysprop.batch import _run_one
        name, result = _run_one(('bad', 'NOT_A_SMILES', 298.15, frozenset(['bp_nannoolal']), 'bp_nannoolal'))
        assert result is None


# ---------------------------------------------------------------------------
# db round-trip
# ---------------------------------------------------------------------------

class TestDbRoundTrip:

    def test_save_and_reload(self, tmp_path):
        db_path = str(tmp_path / 'session.db')
        results = run_all(_mols(), methods=['bp_nannoolal', 'vp_simpol'],
                          workers=1, use_threads=True)
        db.save_results(db_path, results, temperature=298.15)
        loaded = db.load_results(db_path)
        assert len(loaded) == len(SMILES)
        names = {r['name'] for r in loaded}
        assert names == {name for name, _ in SMILES}

    def test_list_compounds(self, tmp_path):
        db_path = str(tmp_path / 'session.db')
        results = run_all(_mols(), methods=['bp_nannoolal'],
                          workers=1, use_threads=True)
        db.save_results(db_path, results, temperature=298.15)
        compounds = db.list_compounds(db_path)
        assert len(compounds) == len(SMILES)

    def test_export_csv(self, tmp_path):
        db_path  = str(tmp_path / 'session.db')
        csv_path = str(tmp_path / 'export.csv')
        results = run_all(_mols(), methods=['bp_nannoolal'],
                          workers=1, use_threads=True)
        db.save_results(db_path, results, temperature=298.15)
        db.export_csv(db_path, csv_path)
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == len(SMILES)
