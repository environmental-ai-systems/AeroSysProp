"""EVAPORATION parameter sets 1 and 2 (method vp_evaporation).

The two sets share equations and coefficients and differ in three SMARTS
patterns (804, 805 and 1902). These tests freeze where they differ on the
MCM, what they give there, and which groups cause it, so any change to
either set shows up here. Regenerate evaporation_set_fixtures.json with
generate_evaporation_set_fixtures.py only for an approved change.

Two kinds of molecule differ:

* enols (78 MCM species): the intended change in set 2, which treats the
  enol OH as an alkyl alcohol, following the published method;
* alpha-hydroxy carboxylic acids (13 species, e.g. glycolic acid): set 2
  widened the diol pattern 1902 to [CX4,CX3] to catch enols, and it now also
  matches the carboxylic acid carbon, adding a diol correction. This looks
  unintended and is recorded in the method card for the core team; these
  tests describe current behaviour, not endorse it.
"""
import json
import warnings
from pathlib import Path

import pytest
from rdkit import Chem

from aerosysprop import groups
from aerosysprop.batch import run_all
from aerosysprop.methods import vapour_pressures as vp

FIX = json.loads((Path(__file__).parent / 'evaporation_set_fixtures.json').read_text())
T = FIX['temperature_K']
DIFFERING = FIX['differing']
IDS = [d['name'] for d in DIFFERING]


def _nonzero(d):
    return {k: v for k, v in d.items() if v}


def test_fixture_summary():
    cats = [d['category'] for d in DIFFERING]
    assert len(DIFFERING) == 91
    assert cats.count('enol') == 78 and cats.count('alpha_hydroxy_acid') == 13


@pytest.mark.parametrize('case', DIFFERING, ids=IDS)
def test_set_values_are_unchanged(case):
    mol = Chem.MolFromSmiles(case['smiles'])
    assert vp.evaporation(mol, T) == pytest.approx(case['set1_log10atm'], rel=1e-9)
    assert vp.evaporation2(mol, T) == pytest.approx(case['set2_log10atm'], rel=1e-9)


@pytest.mark.parametrize('case', DIFFERING, ids=IDS)
def test_group_differences_are_unchanged(case):
    mol = Chem.MolFromSmiles(case['smiles'])
    g1, g2 = _nonzero(groups.evaporation(mol)), _nonzero(groups.evaporation2(mol))
    got = {k: [g1.get(k, 0), g2.get(k, 0)] for k in set(g1) | set(g2)
           if g1.get(k, 0) != g2.get(k, 0)}
    assert got == case['group_differences']


def test_alpha_hydroxy_acids_differ_only_by_the_diol_correction():
    for case in DIFFERING:
        if case['category'] == 'alpha_hydroxy_acid':
            (key, (n1, n2)), = case['group_differences'].items()
            assert key == '19' and n2 - n1 == 2, case['name']


def test_both_sets_through_run_all():
    mols = [(d['name'], Chem.MolFromSmiles(d['smiles'])) for d in DIFFERING]
    rows = run_all(mols, methods=['vp_evaporation'], temperature=T, workers=1,
                   use_threads=True, parameter_sets={'vp_evaporation': [1, 2]},
                   progress_callback=lambda *_: None)
    for row, case in zip(rows, DIFFERING):
        assert row['vp_evaporation_log10atm__set1'] == pytest.approx(case['set1_log10atm'], rel=1e-9)
        assert row['vp_evaporation_log10atm__set2'] == pytest.approx(case['set2_log10atm'], rel=1e-9)


@pytest.mark.parametrize('control', FIX['controls'], ids=[c['smiles'] for c in FIX['controls']])
def test_sets_agree_on_controls(control):
    mol = Chem.MolFromSmiles(control['smiles'])
    assert vp.evaporation(mol, T) == pytest.approx(control['log10atm'], rel=1e-9)
    assert vp.evaporation2(mol, T) == pytest.approx(control['log10atm'], rel=1e-9)


def test_sets_differ_on_no_other_mcm_species():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        '_gen_evap', Path(__file__).parent / 'generate_evaporation_set_fixtures.py')
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    mcm_species = gen.mcm_species
    known = {d['smiles'] for d in DIFFERING}
    extra = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for name, smi in mcm_species():
            if smi in known:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                if abs(vp.evaporation(mol, T) - vp.evaporation2(mol, T)) > 1e-9:
                    extra.append(name)
            except Exception:
                continue
    assert extra == []
