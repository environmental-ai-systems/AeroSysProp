"""Tests for the bisulfate dissociation equilibrium HSO4- <-> H+ + SO4--
(AIOMFAC Phase 3), added as an option to the electrolyte model.

The equilibrium constant and condition follow Fortran AIOMFAC
(SubModDissociationEquil.f90). A Fortran *reference* fixture for H2SO4 is still
to be generated; these tests lock the equilibrium-constant value, the
conservation laws, the self-consistency of the solve, the physical dilution
trend, the on/off option, and a regression snapshot of the Python port.
"""
import numpy as np
import pytest

from aerosysprop.methods.activity_coefficients import (
    AIOMFAC, aiomfac_from_smiles, _component_from_smiles, _lnK_bisulfate,
    _ION_H, _ION_HSO4, _ION_SO4,
)


def _water_h2so4(bisulfate=True):
    return AIOMFAC([_component_from_smiles('O')], salts=['H2SO4'],
                   dissociate_bisulfate=bisulfate)


# ---- equilibrium constant -------------------------------------------------

def test_lnK_at_298_matches_fortran():
    # SubModDissociationEquil.f90 reduces to this constant at the reference T.
    assert _lnK_bisulfate(298.15) == pytest.approx(-4.54916799587, abs=1e-11)
    assert np.exp(_lnK_bisulfate(298.15)) == pytest.approx(0.01058, abs=1e-4)


def test_lnK_temperature_dependence_and_clamping():
    # K decreases with T over ambient range (bisulfate associates more when warm).
    assert _lnK_bisulfate(310.0) < _lnK_bisulfate(298.15) < _lnK_bisulfate(280.0)
    # Out-of-range temperatures clamp rather than extrapolate / error.
    assert _lnK_bisulfate(100.0) == _lnK_bisulfate(180.0)
    assert _lnK_bisulfate(900.0) == _lnK_bisulfate(473.0)


# ---- structure: HSO4- enters the system -----------------------------------

def test_h2so4_system_includes_hso4_anion():
    s = _water_h2so4()
    assert s.dissociate_bisulfate is True
    assert _ION_HSO4 in s.anion_ids          # 248 added
    assert _ION_SO4 in s.anion_ids           # 261 retained
    res = s.gamma([0.9, 0.1], 298.15)
    assert _ION_HSO4 in res.anion_ids
    assert np.all(np.isfinite(res.gamma_anion)) and np.all(res.gamma_anion > 0)
    assert np.all(np.isfinite(res.gamma_cation)) and np.all(res.gamma_cation > 0)


# ---- conservation + self-consistency --------------------------------------

@pytest.mark.parametrize('wf', [0.01, 0.05, 0.10, 0.30, 0.50])
def test_conservation_and_equilibrium_residual(wf):
    s = _water_h2so4()
    iH, iS, iB = s._idx_H, s._idx_SO4, s._idx_HSO4
    comp = s._composition([1.0 - wf, wf])
    H_tot, SO4_tot = comp['SMC'][iH], comp['SMA'][iS]   # reference has HSO4-=0
    c2 = s._solve_bisulfate(comp, 298.15)
    mH, mSO4, mHSO4 = c2['SMC'][iH], c2['SMA'][iS], c2['SMA'][iB]
    # Proton balance and total-sulfate conservation.
    assert mH + mHSO4 == pytest.approx(H_tot, rel=1e-9)
    assert mSO4 + mHSO4 == pytest.approx(SO4_tot, rel=1e-9)
    assert mH > 0 and mSO4 > 0 and mHSO4 > 0
    # Equilibrium is actually satisfied: ln(mH*mSO4/mHSO4)+Δlnγ == lnK.
    _, lc, la = s._ln_gamma_all(c2, 298.15)
    resid = (lc[iH] + la[iS] - la[iB]) + np.log(mH * mSO4 / mHSO4) \
        - _lnK_bisulfate(298.15)
    assert abs(resid) < 1e-6


# ---- robustness with organics present (no oscillation / stall) ------------

@pytest.mark.parametrize('smiles,wtf', [
    (['O', 'CCO'],        [0.6, 0.3, 0.1]),          # ethanol
    (['O', 'OCC(O)CO'],   [0.5, 0.4, 0.1]),          # glycerol
    (['O', 'CC(=O)O'],    [0.6, 0.3, 0.1]),          # acetic acid
    (['O', 'OCC(O)CO'],   [0.2, 0.7, 0.1]),          # organic-rich
    (['O', 'CCO', 'CC(=O)O'], [0.4, 0.3, 0.2, 0.1]),  # two organics
])
def test_converges_with_organics(smiles, wtf):
    # The bisulfate solve is a bracketed bisection (not a fixed-point iteration),
    # so organic interactions shift the result but cannot cause it to oscillate or
    # fail to converge: the equilibrium residual must still hit ~0.
    s = AIOMFAC([_component_from_smiles(x) for x in smiles], salts=['H2SO4'])
    iH, iS, iB = s._idx_H, s._idx_SO4, s._idx_HSO4
    c2 = s._solve_bisulfate(s._composition(wtf), 298.15)
    mH, mS, mB = c2['SMC'][iH], c2['SMA'][iS], c2['SMA'][iB]
    assert mH > 0 and mS > 0 and mB > 0
    _, lc, la = s._ln_gamma_all(c2, 298.15)
    resid = (lc[iH] + la[iS] - la[iB]) + np.log(mH * mS / mB) \
        - _lnK_bisulfate(298.15)
    assert abs(resid) < 1e-6


# ---- physical dilution trend ----------------------------------------------

def test_dissociation_increases_on_dilution():
    s = _water_h2so4()
    iH, iS, iB = s._idx_H, s._idx_SO4, s._idx_HSO4

    def alpha(wf):
        comp = s._composition([1.0 - wf, wf])
        denom = min(comp['SMC'][iH], comp['SMA'][iS])
        c2 = s._solve_bisulfate(comp, 298.15)
        return 1.0 - c2['SMA'][iB] / denom

    a_dilute, a_conc = alpha(1e-4), alpha(1e-2)
    assert a_dilute > a_conc            # more dilute -> more dissociated
    assert a_dilute > 0.8               # approaches full dissociation when dilute


# ---- the option: on vs off -------------------------------------------------

def test_option_off_is_full_dissociation():
    s = _water_h2so4(bisulfate=False)
    assert s.dissociate_bisulfate is False
    assert _ION_HSO4 not in s.anion_ids
    comp = s._composition([0.9, 0.1])
    # Fully dissociated: m_H == 2 * m_SO4 (H2SO4 -> 2 H+ + SO4--).
    cidx = s._cation_idx[_ION_H]; aidx = s._anion_idx[_ION_SO4]
    assert comp['SMC'][cidx] == pytest.approx(2.0 * comp['SMA'][aidx], rel=1e-9)


def test_option_changes_results():
    on = aiomfac_from_smiles(['O'], ['H2SO4'], [0.9, 0.1], 298.15,
                             strict='silent', bisulfate=True)
    off = aiomfac_from_smiles(['O'], ['H2SO4'], [0.9, 0.1], 298.15,
                              strict='silent', bisulfate=False)
    assert _ION_HSO4 in on.anion_ids and _ION_HSO4 not in off.anion_ids
    assert not np.isclose(on.gamma_cation[0], off.gamma_cation[0])


def test_non_sulfuric_salt_unaffected():
    # NaCl has neither H+ nor SO4-- -> bisulfate handling stays off.
    s = AIOMFAC([_component_from_smiles('O')], salts=['NaCl'])
    assert s.dissociate_bisulfate is False


# ---- regression snapshot (Python port; Fortran reference pending) ----------

def test_regression_snapshot_water_h2so4():
    """Locks the current Python-port values for water + 10 wt% H2SO4 at 298.15 K.
    Update intentionally if the model changes; replace with a Fortran-reference
    assertion once that fixture is generated."""
    res = aiomfac_from_smiles(['O'], ['H2SO4'], [0.90, 0.10], 298.15,
                              strict='silent')
    # anions reported as [SO4--, HSO4-]
    assert res.anion_ids == [_ION_SO4, _ION_HSO4]
    assert res.gamma_anion[0] == pytest.approx(0.030764, rel=2e-4)   # SO4--
    assert res.gamma_anion[1] == pytest.approx(1.311333, rel=2e-4)   # HSO4-
    assert res.gamma_cation[0] == pytest.approx(0.635590, rel=2e-4)  # H+
    assert res.gamma_neutral[0] == pytest.approx(0.999138, rel=2e-4)  # water
