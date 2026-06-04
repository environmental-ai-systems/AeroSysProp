"""Tests for AIOMFAC mixture persistence in aerosysprop.db."""

import os
import tempfile

import pytest

from aerosysprop import db


@pytest.fixture
def db_path():
    fd, path = tempfile.mkstemp(suffix='.db')
    os.close(fd)
    yield path
    if os.path.exists(path):
        os.unlink(path)


def _example_rows():
    return [
        {'species': 'water',   'smiles': 'O',   'kind': 'neutral',
         'fraction': 0.85, 'gamma': 1.0006},
        {'species': 'ethanol', 'smiles': 'CCO', 'kind': 'neutral',
         'fraction': 0.10, 'gamma': 2.7},
        {'species': 'cation 202', 'smiles': None, 'kind': 'cation',
         'fraction': None, 'gamma': 0.74},
        {'species': 'anion 242',  'smiles': None, 'kind': 'anion',
         'fraction': None, 'gamma': 0.73},
        {'species': 'NaCl (mean ionic)', 'smiles': None, 'kind': 'salt_mean',
         'fraction': None, 'gamma': 0.735},
    ]


def test_save_and_list_mixture(db_path):
    mid = db.save_mixture_result(
        db_path, name='test_mix', species_rows=_example_rows(),
        temperature=298.15, basis='mass', salts=['NaCl'],
    )
    assert isinstance(mid, int) and mid > 0

    mixes = db.list_mixtures(db_path)
    assert len(mixes) == 1
    m = mixes[0]
    assert m['name'] == 'test_mix'
    assert m['temperature'] == 298.15
    assert m['basis'] == 'mass'
    assert m['salts'] == 'NaCl'


def test_load_mixture_species_roundtrip(db_path):
    rows = _example_rows()
    mid = db.save_mixture_result(
        db_path, name='roundtrip', species_rows=rows,
        temperature=300.0, salts=['NaCl'],
    )
    loaded = db.load_mixture_species(db_path, mid)
    assert len(loaded) == len(rows)
    # neutral water row preserved
    water = next(r for r in loaded if r['species'] == 'water')
    assert water['kind'] == 'neutral'
    assert water['smiles'] == 'O'
    assert water['fraction'] == pytest.approx(0.85)
    assert water['gamma'] == pytest.approx(1.0006)
    # ion row has null fraction
    cation = next(r for r in loaded if r['kind'] == 'cation')
    assert cation['fraction'] is None
    assert cation['gamma'] == pytest.approx(0.74)


def test_multiple_mixtures_newest_first(db_path):
    db.save_mixture_result(db_path, name='first',  species_rows=_example_rows(),
                           temperature=298.15, salts=['NaCl'])
    db.save_mixture_result(db_path, name='second', species_rows=_example_rows(),
                           temperature=310.0, salts=['KCl'])
    mixes = db.list_mixtures(db_path)
    assert len(mixes) == 2
    # newest first
    assert mixes[0]['name'] == 'second'
    assert mixes[1]['name'] == 'first'


def test_no_salts_mixture(db_path):
    rows = [r for r in _example_rows() if r['kind'] == 'neutral']
    mid = db.save_mixture_result(
        db_path, name='sr_only', species_rows=rows,
        temperature=298.15, salts=None,
    )
    mixes = db.list_mixtures(db_path)
    assert mixes[0]['salts'] == ''
    loaded = db.load_mixture_species(db_path, mid)
    assert all(r['kind'] == 'neutral' for r in loaded)


def test_mixture_tables_coexist_with_compound_results(db_path):
    """The new mixture tables must not interfere with the existing
    compound/results persistence."""
    db.save_results(
        db_path,
        [{'name': 'ethanol', 'smiles': 'CCO', 'bp_nannoolal_K': 351.0}],
        temperature=298.15,
    )
    db.save_mixture_result(
        db_path, name='mix', species_rows=_example_rows(),
        temperature=298.15, salts=['NaCl'],
    )
    # Both stores readable
    assert len(db.list_compounds(db_path)) == 1
    assert len(db.list_mixtures(db_path)) == 1
