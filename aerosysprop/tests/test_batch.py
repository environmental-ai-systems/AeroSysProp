"""
Tests for aerosysprop.batch — parallel property prediction.
"""
import os

import pytest
from rdkit import Chem

from aerosysprop.batch import run_all

FIXTURES = os.path.join(os.path.dirname(__file__), 'io_fixtures')
SMILES = [
    ('methyl_hydroperoxide', 'COO'),
    ('ethylene_glycol',      'C(O)CO'),
    ('acetone',              'CC(C)=O'),
    ('acetic_acid',          'CC(=O)O'),
]


def _mols():
    return [(name, Chem.MolFromSmiles(smi)) for name, smi in SMILES]


class TestRunAll:

    def test_returns_one_result_per_molecule(self):
        results = run_all(_mols(), workers=2)
        assert len(results) == len(SMILES)

    def test_result_order_matches_input(self):
        mols = _mols()
        results = run_all(mols, workers=2)
        for (name, _), result in zip(mols, results):
            assert result['name'] == name

    def test_expected_keys_present(self):
        results = run_all(_mols(), workers=2)
        expected_keys = {
            'name', 'smiles',
            'bp_joback_K', 'bp_stein_K', 'bp_nannoolal_K',
            'vp_nannoolal_log10atm', 'vp_myrdal_log10atm',
            'vp_evaporation_log10atm', 'vp_evaporation2_log10atm',
            'vp_simpol_log10atm',
            'tc_nannoolal_K', 'pc_nannoolal_bar',
            'density_girolami_gcc', 'density_schroeder_gcc',
            'density_le_bas_gcc', 'density_tyn_calus_gcc',
        }
        for result in results:
            assert expected_keys <= result.keys()

    def test_numeric_predictions_are_floats(self):
        results = run_all(_mols(), workers=2)
        for result in results:
            assert isinstance(result['bp_nannoolal_K'], float)
            assert result['bp_nannoolal_K'] > 0
            assert isinstance(result['vp_simpol_log10atm'], float)

    def test_single_worker(self):
        """workers=1 should produce identical results to multi-worker."""
        mols = _mols()
        r1 = run_all(mols, workers=1)
        r2 = run_all(mols, workers=2)
        for a, b in zip(r1, r2):
            assert a['bp_nannoolal_K'] == pytest.approx(b['bp_nannoolal_K'])

    def test_csv_output(self, tmp_path):
        """output_csv should produce a CSV with a header and one row per mol."""
        import csv
        csv_path = str(tmp_path / 'out.csv')
        run_all(_mols(), workers=2, output_csv=csv_path)
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == len(SMILES)
        assert 'bp_nannoolal_K' in rows[0]

    def test_invalid_smiles_returns_none(self):
        """A molecule that fails to parse should appear as None in results."""
        mols = [('bad', Chem.MolFromSmiles('C')), ('good', Chem.MolFromSmiles('CCO'))]
        # Corrupt the first mol object by replacing it with a valid one then
        # testing that None SMILES mol is handled — simulate via direct call
        from aerosysprop.batch import _run_one
        name, result = _run_one(('bad', 'NOT_A_SMILES', 298.15))
        assert result is None
