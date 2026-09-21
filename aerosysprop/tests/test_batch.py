"""
Tests for aerosysprop.batch — parallel property prediction.
"""
import csv
import os
import tempfile

import pytest
from rdkit import Chem

from aerosysprop.batch import (
    ALL_METHODS, _OUTPUT_KEYS, resolve_methods, run_all,
    resolved_output_columns,
    PROPERTY_METHODS, PROPERTY_METHOD_BY_ID, PropertyMethod,
    METHOD_META, METHOD_GROUPS, BP_METHODS,
)
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
# PropertyMethod registry — derived tables must stay consistent with it
# ---------------------------------------------------------------------------

class TestPropertyRegistry:

    def test_registry_nonempty(self):
        assert PROPERTY_METHODS
        assert all(isinstance(m, PropertyMethod) for m in PROPERTY_METHODS)

    def test_all_methods_matches_registry_order(self):
        assert ALL_METHODS == [m.id for m in PROPERTY_METHODS]

    def test_method_meta_derived(self):
        for m in PROPERTY_METHODS:
            assert METHOD_META[m.id] == (m.label, m.unit, m.group, m.deps)

    def test_output_keys_derived(self):
        for m in PROPERTY_METHODS:
            assert _OUTPUT_KEYS[m.id] == m.output_keys

    def test_method_groups_derived_and_ordered(self):
        # every method appears under its declared group, in registry order
        for m in PROPERTY_METHODS:
            assert m.id in METHOD_GROUPS[m.group]
        # group order = first appearance in the registry
        seen = []
        for m in PROPERTY_METHODS:
            if m.group not in seen:
                seen.append(m.group)
        assert list(METHOD_GROUPS) == seen

    def test_bp_methods_are_the_boiling_point_group(self):
        assert BP_METHODS == frozenset(
            m.id for m in PROPERTY_METHODS if m.group == 'Boiling Point'
        )

    def test_lookup_by_id(self):
        assert PROPERTY_METHOD_BY_ID['bp_nannoolal'].label == 'Nannoolal et al. (2004)'

    def test_bp_dependent_flag_consistency(self):
        # density methods depend on cp_nannoolal, which is bp_dependent ->
        # they must be flagged bp_dependent too
        for mid in ('vp_nannoolal', 'vp_myrdal', 'cp_nannoolal',
                    'density_schroeder', 'density_le_bas', 'density_tyn_calus'):
            assert PROPERTY_METHOD_BY_ID[mid].bp_dependent
        for mid in ('bp_joback', 'vp_simpol', 'density_girolami'):
            assert not PROPERTY_METHOD_BY_ID[mid].bp_dependent


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
        all_keys = set(resolved_output_columns(set(ALL_METHODS), ['bp_nannoolal']))
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

    def test_single_bp_source_column_names_unsuffixed(self):
        """A single bp source must keep the original (unsuffixed) column names
        for backward compatibility."""
        results = run_all(_mols(), methods=['vp_nannoolal'], workers=1,
                          use_threads=True, bp_sources=['bp_stein'])
        for r in results:
            if r:
                assert 'vp_nannoolal_log10atm' in r
                assert not any('__bp' in k for k in r)

    def test_multi_bp_sources_suffix_columns(self):
        """Multiple bp sources produce one suffixed column per source."""
        results = run_all(
            _mols(), methods=['vp_nannoolal', 'cp_nannoolal'], workers=1,
            use_threads=True, bp_sources=['bp_nannoolal', 'bp_stein'],
        )
        r = next(x for x in results if x)
        for key in ('vp_nannoolal_log10atm', 'tc_nannoolal_K'):
            assert f'{key}__bpNannoolal' in r
            assert f'{key}__bpStein' in r
        # Unsuffixed key must NOT be present in the multi-source case
        assert 'vp_nannoolal_log10atm' not in r

    def test_multi_bp_sources_values_differ(self):
        """The per-source VP values should differ when the BP methods differ."""
        results = run_all(
            _mols(), methods=['vp_nannoolal'], workers=1, use_threads=True,
            bp_sources=['bp_nannoolal', 'bp_joback'],
        )
        diffs = [
            r['vp_nannoolal_log10atm__bpNannoolal']
            != r['vp_nannoolal_log10atm__bpJoback']
            for r in results
            if r and r.get('vp_nannoolal_log10atm__bpNannoolal') is not None
            and r.get('vp_nannoolal_log10atm__bpJoback') is not None
        ]
        assert diffs and any(diffs)

    def test_resolved_output_columns_single_vs_multi(self):
        mset = resolve_methods(['vp_nannoolal'], bp_sources=['bp_nannoolal'])
        single = resolved_output_columns(mset, ['bp_nannoolal'])
        assert 'vp_nannoolal_log10atm' in single
        assert not any('__bp' in c for c in single)

        mset2 = resolve_methods(['vp_nannoolal'],
                                bp_sources=['bp_nannoolal', 'bp_stein'])
        multi = resolved_output_columns(mset2, ['bp_nannoolal', 'bp_stein'])
        assert 'vp_nannoolal_log10atm__bpNannoolal' in multi
        assert 'vp_nannoolal_log10atm__bpStein' in multi

    def test_invalid_bp_source_in_list_raises(self):
        with pytest.raises(ValueError, match='bp_source must be one of'):
            run_all(_mols(), methods=['vp_nannoolal'], workers=1,
                    use_threads=True, bp_sources=['bp_nannoolal', 'not_a_bp'])

    def test_multi_bp_csv_header(self, tmp_path):
        """CSV output must include the suffixed columns for multiple sources."""
        csv_path = str(tmp_path / 'multi.csv')
        run_all(_mols(), methods=['vp_nannoolal'], workers=1, use_threads=True,
                bp_sources=['bp_nannoolal', 'bp_stein'], output_csv=csv_path)
        with open(csv_path) as f:
            header = f.readline()
        assert 'vp_nannoolal_log10atm__bpNannoolal' in header
        assert 'vp_nannoolal_log10atm__bpStein' in header

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


class TestParameterSets:
    """vp_evaporation has two parameter sets; set 2 is the default."""

    def test_default_is_set_2(self):
        r = run_all(_mols(), methods=['vp_evaporation'], workers=1, use_threads=True)[0]
        assert 'vp_evaporation_log10atm__set2' in r
        assert 'vp_evaporation_log10atm__set1' not in r
        assert 'vp_evaporation_log10atm' not in r

    def test_both_sets(self):
        r = run_all(_mols(), methods=['vp_evaporation'], workers=1, use_threads=True,
                    parameter_sets={'vp_evaporation': [1, 2]})[0]
        assert {'vp_evaporation_log10atm__set1', 'vp_evaporation_log10atm__set2'} <= r.keys()

    def test_sets_match_the_underlying_functions(self):
        from aerosysprop.methods import vapour_pressures
        mol = Chem.MolFromSmiles('OC=CC(=O)O')  # an enol: the sets treat it differently
        r = run_all([('enol', mol)], methods=['vp_evaporation'], workers=1, use_threads=True,
                    parameter_sets={'vp_evaporation': [1, 2]})[0]
        assert r['vp_evaporation_log10atm__set1'] == pytest.approx(
            vapour_pressures.evaporation(mol, 298.15))
        assert r['vp_evaporation_log10atm__set2'] == pytest.approx(
            vapour_pressures.evaporation2(mol, 298.15))

    def test_deprecated_id_still_works(self):
        with pytest.warns(FutureWarning):
            r = run_all(_mols(), methods=['vp_evaporation2'], workers=1, use_threads=True)[0]
        assert 'vp_evaporation_log10atm__set2' in r

    def test_unknown_set_is_rejected(self):
        with pytest.raises(ValueError):
            run_all(_mols(), methods=['vp_evaporation'], workers=1, use_threads=True,
                    parameter_sets={'vp_evaporation': [3]})
