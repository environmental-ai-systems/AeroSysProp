"""
Regression tests for all property-prediction methods.

Reference values are stored in fixtures.json, generated once from the
current implementation via generate_fixtures.py.  Re-run the generator
only when you intentionally change a method.
"""
import json
import math
import os

import pytest
from rdkit import Chem

from aerosysprop.methods import (
    boiling_points,
    critical_properties,
    liquid_densities,
    vapour_pressures,
)

TEMPERATURE = 298.15
_FIXTURE_FILE = os.path.join(os.path.dirname(__file__), 'fixtures.json')

with open(_FIXTURE_FILE) as _f:
    _FIXTURES = json.load(_f)

_SMILES = list(_FIXTURES.keys())


def _mol(smiles):
    return Chem.MolFromSmiles(smiles)


def _approx(expected):
    """pytest.approx wrapper that passes None (NaN/inf) through unchanged."""
    if expected is None:
        return expected
    return pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Boiling points
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('smiles', _SMILES)
def test_bp_joback_and_reid(smiles):
    expected = _FIXTURES[smiles]['boiling_point']['joback_and_reid']
    assert boiling_points.joback_and_reid(_mol(smiles)) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_bp_stein_and_brown(smiles):
    expected = _FIXTURES[smiles]['boiling_point']['stein_and_brown']
    assert boiling_points.stein_and_brown(_mol(smiles)) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_bp_nannoolal(smiles):
    expected = _FIXTURES[smiles]['boiling_point']['nannoolal']
    assert boiling_points.nannoolal(_mol(smiles)) == _approx(expected)


# ---------------------------------------------------------------------------
# Vapour pressures  (use fixture boiling points to avoid propagating BP drift)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_nannoolal_bp_nannoolal(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['nannoolal']
    expected = _FIXTURES[smiles]['vapour_pressure']['nannoolal_bp_nannoolal']
    assert vapour_pressures.nannoolal(mol, TEMPERATURE, bp) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_nannoolal_bp_stein_and_brown(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['stein_and_brown']
    expected = _FIXTURES[smiles]['vapour_pressure']['nannoolal_bp_stein_and_brown']
    assert vapour_pressures.nannoolal(mol, TEMPERATURE, bp) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_nannoolal_bp_joback_and_reid(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['joback_and_reid']
    expected = _FIXTURES[smiles]['vapour_pressure']['nannoolal_bp_joback_and_reid']
    assert vapour_pressures.nannoolal(mol, TEMPERATURE, bp) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_myrdal_yalkowsky_bp_nannoolal(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['nannoolal']
    expected = _FIXTURES[smiles]['vapour_pressure']['myrdal_yalkowsky_bp_nannoolal']
    assert vapour_pressures.myrdal_and_yalkowsky(mol, TEMPERATURE, bp) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_myrdal_yalkowsky_bp_stein_and_brown(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['stein_and_brown']
    expected = _FIXTURES[smiles]['vapour_pressure']['myrdal_yalkowsky_bp_stein_and_brown']
    assert vapour_pressures.myrdal_and_yalkowsky(mol, TEMPERATURE, bp) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_myrdal_yalkowsky_bp_joback_and_reid(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['joback_and_reid']
    expected = _FIXTURES[smiles]['vapour_pressure']['myrdal_yalkowsky_bp_joback_and_reid']
    assert vapour_pressures.myrdal_and_yalkowsky(mol, TEMPERATURE, bp) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_evaporation(smiles):
    expected = _FIXTURES[smiles]['vapour_pressure']['evaporation']
    assert vapour_pressures.evaporation(_mol(smiles), TEMPERATURE) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_evaporation2(smiles):
    expected = _FIXTURES[smiles]['vapour_pressure']['evaporation2']
    assert vapour_pressures.evaporation2(_mol(smiles), TEMPERATURE) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_vp_simpol(smiles):
    expected = _FIXTURES[smiles]['vapour_pressure']['simpol']
    assert vapour_pressures.simpol(_mol(smiles), TEMPERATURE) == _approx(expected)


# ---------------------------------------------------------------------------
# Liquid densities
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('smiles', _SMILES)
def test_density_girolami(smiles):
    expected = _FIXTURES[smiles]['liquid_density']['girolami']
    result = liquid_densities.girolami(_mol(smiles))
    assert float(result.real) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_density_schroeder_cp_nannoolal(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['nannoolal']
    cp = critical_properties.nannoolal(mol, bp)
    expected = _FIXTURES[smiles]['liquid_density']['schroeder_cp_nannoolal']
    assert float(liquid_densities.schroeder(mol, TEMPERATURE, cp).real) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_density_schroeder_cp_joback(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['joback_and_reid']
    cp = critical_properties.joback_and_reid(mol, bp)
    expected = _FIXTURES[smiles]['liquid_density']['schroeder_cp_joback']
    assert float(liquid_densities.schroeder(mol, TEMPERATURE, cp).real) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_density_le_bas_cp_nannoolal(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['nannoolal']
    cp = critical_properties.nannoolal(mol, bp)
    expected = _FIXTURES[smiles]['liquid_density']['le_bas_cp_nannoolal']
    assert float(liquid_densities.le_bas(mol, TEMPERATURE, cp).real) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_density_le_bas_cp_joback(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['joback_and_reid']
    cp = critical_properties.joback_and_reid(mol, bp)
    expected = _FIXTURES[smiles]['liquid_density']['le_bas_cp_joback']
    assert float(liquid_densities.le_bas(mol, TEMPERATURE, cp).real) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_density_tyn_calus_cp_nannoolal(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['nannoolal']
    cp = critical_properties.nannoolal(mol, bp)
    expected = _FIXTURES[smiles]['liquid_density']['tyn_calus_cp_nannoolal']
    assert float(liquid_densities.tyn_and_calus(mol, TEMPERATURE, cp).real) == _approx(expected)


@pytest.mark.parametrize('smiles', _SMILES)
def test_density_tyn_calus_cp_joback(smiles):
    mol = _mol(smiles)
    bp = _FIXTURES[smiles]['boiling_point']['joback_and_reid']
    cp = critical_properties.joback_and_reid(mol, bp)
    expected = _FIXTURES[smiles]['liquid_density']['tyn_calus_cp_joback']
    assert float(liquid_densities.tyn_and_calus(mol, TEMPERATURE, cp).real) == _approx(expected)
