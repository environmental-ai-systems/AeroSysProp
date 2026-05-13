#!/usr/bin/env python
"""
Generate regression fixture data for unit tests.

Run once from the repo root:

    python -m aerosysprop.tests.generate_fixtures

or directly:

    python aerosysprop/tests/generate_fixtures.py

Only re-run when you intentionally change a method's implementation.
The output (fixtures.json) should be committed to the repo so that the
test suite can be run without executing this script first.
"""
import json
import math
import os
import sys

from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..'))

from aerosysprop.methods import (
    boiling_points,
    critical_properties,
    liquid_densities,
    vapour_pressures,
)

TEMPERATURE = 298.15
SMILES_FILE = os.path.join(HERE, 'SMILES', 'short_list.txt')
FIXTURE_FILE = os.path.join(HERE, 'fixtures.json')


def _read_smiles(path):
    """Read SMILES from a tab-separated file: <index>\t<SMILES>\t[value]"""
    smiles = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                smiles.append(parts[1])
    return smiles


def _safe(val):
    """Return None for NaN/inf so JSON stays valid."""
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return None
    return val


def _compute(mol):
    bp_joback    = boiling_points.joback_and_reid(mol)
    bp_stein     = boiling_points.stein_and_brown(mol)
    bp_nannoolal = boiling_points.nannoolal(mol)

    cp_joback    = critical_properties.joback_and_reid(mol, bp_joback)
    cp_nannoolal = critical_properties.nannoolal(mol, bp_nannoolal)

    return {
        'boiling_point': {
            'joback_and_reid': _safe(bp_joback),
            'stein_and_brown': _safe(bp_stein),
            'nannoolal':       _safe(bp_nannoolal),
        },
        'vapour_pressure': {
            'nannoolal_bp_nannoolal':              _safe(vapour_pressures.nannoolal(mol, TEMPERATURE, bp_nannoolal)),
            'nannoolal_bp_stein_and_brown':        _safe(vapour_pressures.nannoolal(mol, TEMPERATURE, bp_stein)),
            'nannoolal_bp_joback_and_reid':        _safe(vapour_pressures.nannoolal(mol, TEMPERATURE, bp_joback)),
            'myrdal_yalkowsky_bp_nannoolal':       _safe(vapour_pressures.myrdal_and_yalkowsky(mol, TEMPERATURE, bp_nannoolal)),
            'myrdal_yalkowsky_bp_stein_and_brown': _safe(vapour_pressures.myrdal_and_yalkowsky(mol, TEMPERATURE, bp_stein)),
            'myrdal_yalkowsky_bp_joback_and_reid': _safe(vapour_pressures.myrdal_and_yalkowsky(mol, TEMPERATURE, bp_joback)),
            'evaporation':  _safe(vapour_pressures.evaporation(mol, TEMPERATURE)),
            'evaporation2': _safe(vapour_pressures.evaporation2(mol, TEMPERATURE)),
            'simpol':       _safe(vapour_pressures.simpol(mol, TEMPERATURE)),
        },
        'liquid_density': {
            'girolami':               _safe(liquid_densities.girolami(mol)),
            'schroeder_cp_nannoolal': _safe(liquid_densities.schroeder(mol, TEMPERATURE, cp_nannoolal)),
            'schroeder_cp_joback':    _safe(liquid_densities.schroeder(mol, TEMPERATURE, cp_joback)),
            'le_bas_cp_nannoolal':    _safe(liquid_densities.le_bas(mol, TEMPERATURE, cp_nannoolal)),
            'le_bas_cp_joback':       _safe(liquid_densities.le_bas(mol, TEMPERATURE, cp_joback)),
            'tyn_calus_cp_nannoolal': _safe(liquid_densities.tyn_and_calus(mol, TEMPERATURE, cp_nannoolal)),
            'tyn_calus_cp_joback':    _safe(liquid_densities.tyn_and_calus(mol, TEMPERATURE, cp_joback)),
        },
    }


if __name__ == '__main__':
    smiles_list = _read_smiles(SMILES_FILE)
    fixtures = {}
    errors = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print('SKIP  %s  (invalid SMILES)' % smi)
            continue
        try:
            fixtures[smi] = _compute(mol)
            print('OK    %s' % smi)
        except Exception as e:
            errors.append((smi, str(e)))
            print('ERR   %s: %s' % (smi, e))

    with open(FIXTURE_FILE, 'w') as f:
        json.dump(fixtures, f, indent=2)

    print('\nWrote %d entries to %s' % (len(fixtures), FIXTURE_FILE))
    if errors:
        print('%d failed:' % len(errors))
        for smi, err in errors:
            print('  %s: %s' % (smi, err))
