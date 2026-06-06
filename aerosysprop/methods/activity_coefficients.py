# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2026 David Topping.
# All Rights Reserved.
# This file is part of aerosysprop.
#
# aerosysprop is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# aerosysprop is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# aerosysprop.  If not, see <http://www.gnu.org/licenses/>.

"""Activity coefficients via AIOMFAC (short-range part).

This module implements the **short-range (SR)** part of the AIOMFAC model
(Zuend et al. 2008, 2011) — a modified UNIFAC formulation that predicts
mole-fraction-based activity coefficients gamma_i for non-electrolyte
mixtures (organic compounds + water).

The SR part alone is a complete predictor for non-electrolyte systems; the
middle-range (MR, organic-ion and ion-ion) and long-range (LR, Debye-Hueckel)
contributions of full AIOMFAC are zero in the absence of ions and are not
implemented here. See ``activity_coefficients.md`` in this directory for full
theory, validation methodology, compound coverage (1520 regression fixtures),
known limitations, and citations.

Public API
----------
- ``AIOMFAC_SR(components)`` -- cached system for repeated evaluation
- ``aiomfac_sr_from_smiles(smiles, x, T_K)`` -- one-shot SMILES API

Example
-------
    >>> from aerosysprop.methods.activity_coefficients import (
    ...     aiomfac_sr_from_smiles,
    ... )
    >>> # gamma_water, gamma_ethanol at x = 50:50, T = 25 C
    >>> aiomfac_sr_from_smiles(['O', 'CCO'], [0.5, 0.5], 298.15)
    array([1.26196..., 1.18692...])

Equations and parameter conventions follow the AIOMFAC Fortran reference
(``ModSRunifac.f90`` / ``ModSRparam.f90`` in the public AIOMFAC repo:
https://github.com/andizuend/AIOMFAC). The SR formulation:

    ln gamma_i^SR = ln gamma_i^C + ln gamma_i^R - ln gamma_i^R,ref

with Voutsas-Tassios (1997) combinatorial term and standard UNIFAC residual
term; the temperature law is the AIOMFAC-default 1-parameter form
Psi(m,k) = exp(-a_mk / T). The pure-component reference subtraction enforces
gamma_i -> 1 as x_i -> 1.

Citations
---------
Zuend, A., Marcolli, C., Luo, B. P., and Peter, T. (2008). ACP 8, 4559-4593.
Zuend, A., et al. (2011). ACP 11, 9155-9206.
Fredenslund, A., et al. (1975). AIChE J. 21, 1086-1099.
Voutsas, E. C., and Tassios, D. P. (1997). Ind. Eng. Chem. Res. 36, 4965-4972.
"""

from collections.abc import Mapping, Sequence

import numpy as np
from rdkit import Chem

from aerosysprop import data, groups


_Z_COORDINATION_HALF = 5.0  # z/2 with z = 10 (hardcoded in AIOMFAC SRgcomb)


class AIOMFAC_SR:
    """Cached system for short-range AIOMFAC activity-coefficient evaluation.

    A component is a mapping ``{subgroup_id: count}``; counts may be int or
    float (group SMARTS matching can yield fractional counts for ring-shared
    atoms). Subgroup IDs follow the AIOMFAC numbering convention (1-200 for
    neutral organic subgroups); use ``aerosysprop.groups.aiomfac(rdkit_mol)``
    to produce these dicts from SMILES.

    Build the system once for a fixed list of components, then call
    ``ln_gamma(x, T_K)`` or ``gamma(x, T_K)`` for any composition. The
    per-system setup (subgroup-index map, nu matrix, R/Q/main-group lookups,
    pure-component reference state) is computed at construction; only the
    composition-dependent and temperature-dependent terms are evaluated
    per call. Pure-component reference Gamma values are cached and reused
    until ``T_K`` changes between calls.

    Parameters
    ----------
    components : sequence of mapping ``{int: number}``
        One dict per component, in the order the mole-fraction vector
        will use. Counts of 0 are ignored.

    Examples
    --------
    >>> from aerosysprop import groups
    >>> from rdkit import Chem
    >>> water = groups.aiomfac(Chem.MolFromSmiles('O'))
    >>> ethanol = groups.aiomfac(Chem.MolFromSmiles('CCO'))
    >>> sys = AIOMFAC_SR([water, ethanol])
    >>> sys.gamma([0.5, 0.5], 298.15)
    array([1.26196..., 1.18692...])

    Notes
    -----
    Internally uses the AIOMFAC-default 1-parameter Psi temperature law
    (Psi(m,k) = exp(-a_mk / T)) and Voutsas-Tassios (1997) combinatorial
    form (coordination number z = 10 hardcoded, matching the Fortran).
    The pure-component reference state (gamma_i -> 1 as x_i -> 1) is
    enforced by subtracting Gamma_k^(i) computed at x_j = delta_ij.
    """

    def __init__(self, components, n_neutral=None):
        comps = [dict(c) for c in components]
        # Sorted union of subgroup IDs actually present anywhere in the mixture
        subgroup_ids = sorted({
            g for c in comps for g, n in c.items() if n != 0
        })
        if not subgroup_ids:
            raise ValueError('At least one component must contain subgroups')

        ng = len(subgroup_ids)
        nk = len(comps)
        # n_neutral is the count of leading neutral components. Anything past
        # that index in `components` is treated as an ion (single ion subgroup
        # >= 201) and gets the infinite-dilution combinatorial correction the
        # Fortran applies in SRgcomb (lnGaCinf, reference is pure water).
        if n_neutral is None:
            n_neutral = nk
        self.n_neutral = n_neutral
        self.components = comps
        self.subgroup_ids = subgroup_ids
        self._sg_index = {g: i for i, g in enumerate(subgroup_ids)}

        # nu[i, k] = count of subgroup k in component i
        nu = np.zeros((nk, ng), dtype=np.float64)
        for i, c in enumerate(comps):
            for g, n in c.items():
                if n != 0:
                    nu[i, self._sg_index[g]] = float(n)
        self.nu = nu

        # Per-subgroup R_k, Q_k and main group IDs
        self.R = np.array(
            [data.AIOMFAC_RI[g] for g in subgroup_ids], dtype=np.float64
        )
        self.Q = np.array(
            [data.AIOMFAC_QI[g] for g in subgroup_ids], dtype=np.float64
        )
        self.main_group = np.array(
            [data.AIOMFAC_MAIN_GROUP[g] for g in subgroup_ids], dtype=np.intp
        )

        # Per-component aggregated van der Waals volume r_i and surface q_i
        self.r = nu @ self.R
        self.q = nu @ self.Q

        # Subgroup-pair main-group interaction parameter a[k,l] (asymmetric).
        # AIOMFAC convention: interactions between same main group are 0;
        # this is what the reference Fortran asserts in SRsystm.
        a = np.zeros((ng, ng), dtype=np.float64)
        for k in range(ng):
            mk = int(self.main_group[k])
            row = data.AIOMFAC_SR_INTERACTIONS[mk]
            for l in range(ng):
                ml = int(self.main_group[l])
                a[k, l] = 0.0 if mk == ml else row[ml]
        self._a = a

        # Cache of pure-component reference ln gamma^R (T-dependent => filled
        # lazily on each ln_gamma call when temperature changes).
        self._cached_T = None
        self._cached_Psi = None
        self._cached_ln_gamma_R_ref = None

        # Per-ion infinite-dilution combinatorial correction (lnGaCinf in
        # SRgcomb). Reference solvent is pure water: R_w = 0.92, Q_w = 1.40,
        # XL_w = -2.32. This is composition-independent so we compute once.
        self._lnGaCinf = np.zeros(nk, dtype=np.float64)
        if nk > n_neutral:
            r_ion = self.r[n_neutral:]
            q_ion = self.q[n_neutral:]
            XL_ion = 5.0 * (r_ion - q_ion) - r_ion + 1.0
            self._lnGaCinf[n_neutral:] = (
                np.log(r_ion / 0.92)
                + 5.0 * q_ion * np.log((q_ion / 1.40) * (0.92 / r_ion))
                + XL_ion
                - (r_ion / 0.92) * (-2.32)
            )

    # ---- composition-independent, T-dependent pieces ------------------------
    def _psi(self, T_K):
        return np.exp(-self._a / T_K)

    def _ln_gamma_R_at(self, x_groups_weighted, Psi):
        """Compute group residual ln Gamma_k given the raw group counts in the
        mixture (S1_k in the Fortran). Returns array of Gamma_k (length NG)."""
        S2 = x_groups_weighted.sum()
        # group mole fraction in mixture
        Xg = x_groups_weighted / S2
        S3 = (self.Q * Xg).sum()
        theta = self.Q * Xg / S3  # area fraction per group
        S4 = theta @ Psi  # S4_k = sum_m theta_m * Psi[m, k]
        # term: sum_m theta_m * Psi[k, m] / S4_m  -> shape (NG,)
        outgoing = (Psi / S4) @ theta
        return self.Q * (1.0 - np.log(S4) - outgoing)

    def _ensure_temperature_cache(self, T_K):
        if self._cached_T == T_K:
            return
        Psi = self._psi(T_K)
        # Reference Gamma_k^(i) per component, evaluated at pure x_i = 1.
        ln_gamma_R_ref = np.zeros(len(self.components), dtype=np.float64)
        for i in range(len(self.components)):
            nu_i = self.nu[i]
            if (nu_i > 0).sum() <= 1:
                # single subgroup => trivially zero residual contribution
                continue
            Gamma_pure = self._ln_gamma_R_at(nu_i, Psi)
            ln_gamma_R_ref[i] = (nu_i * Gamma_pure).sum()
        self._cached_T = T_K
        self._cached_Psi = Psi
        self._cached_ln_gamma_R_ref = ln_gamma_R_ref

    # ---- public API ---------------------------------------------------------
    def ln_gamma(self, mole_fractions, T_K):
        """Return ln(gamma_i^SR) for the given composition and temperature.

        Computes the AIOMFAC SR activity coefficient as the sum of the
        combinatorial (Voutsas-Tassios form) and residual (UNIFAC form)
        contributions, with pure-component reference state subtraction.

        Parameters
        ----------
        mole_fractions : sequence of float
            Mole fractions of the configured components, in the same order
            as passed at construction. Does not need to sum exactly to 1
            (UNIFAC is a function of x, not normalised x).
        T_K : float
            Temperature in Kelvin. Setting Psi(m,k) = exp(-a_mk / T).

        Returns
        -------
        numpy.ndarray
            Natural logarithm of the activity coefficients, shape
            ``(n_components,)``.

        Raises
        ------
        ValueError
            If ``len(mole_fractions)`` does not equal the number of
            components configured at construction.
        """
        x = np.asarray(mole_fractions, dtype=np.float64)
        if x.shape != (len(self.components),):
            raise ValueError(
                f'mole_fractions length {x.shape} != number of components '
                f'{len(self.components)}'
            )
        self._ensure_temperature_cache(T_K)
        Psi = self._cached_Psi

        # Combinatorial part (Voutsas-Tassios form; SRgcomb in Fortran)
        QSS = (self.q * x).sum()
        RSS = (self.r * x).sum()
        ViV = self.r / RSS
        VQR = ViV * QSS / self.q
        ln_gamma_C = (
            np.log(ViV) + 1.0 - ViV
            - _Z_COORDINATION_HALF * self.q * (np.log(VQR) + 1.0 - VQR)
        )

        # Residual part at mixture composition (SRgres)
        # S1_k = sum_i nu[i,k] * x_i
        S1 = x @ self.nu
        Gamma_mix = self._ln_gamma_R_at(S1, Psi)
        ln_gamma_R = self.nu @ Gamma_mix

        # ln gamma^SR = ln_C + ln_R - ln_R_ref. For ions (idx >= n_neutral),
        # subtract the infinite-dilution combinatorial correction lnGaCinf
        # (Fortran SRgcomb lines 514-520; reference solvent is pure water).
        ln_gamma_SR = ln_gamma_C + ln_gamma_R - self._cached_ln_gamma_R_ref
        ln_gamma_SR -= self._lnGaCinf
        return ln_gamma_SR

    def gamma(self, mole_fractions, T_K):
        """Return gamma_i^SR (linear, not log) for the given composition / T.

        Equivalent to ``numpy.exp(self.ln_gamma(mole_fractions, T_K))``;
        see :meth:`ln_gamma` for parameter semantics.
        """
        return np.exp(self.ln_gamma(mole_fractions, T_K))


def _component_from_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smi!r}')
    return {k: v for k, v in groups.aiomfac(mol).items() if v != 0}


# ----------------------------------------------------------------------------
# System caches keyed by (smiles tuple, salts tuple). When `aiomfac_from_smiles`
# or `aiomfac_sr_from_smiles` is called repeatedly with the same SMILES list
# (the typical pattern for numerical simulations -- partitioning solvers,
# Monte Carlo, equilibrium iterations -- composition varies but compounds are
# fixed), the system construction is amortised to a single one-time cost.
# Each cached system holds its own internal temperature cache and per-system
# parameter arrays, so threads / callers sharing a cached system don't fight.
# ----------------------------------------------------------------------------
import functools


@functools.lru_cache(maxsize=128)
def _cached_aiomfac_system(smiles_tuple, salts_tuple, dissociate_bisulfate=True):
    neutrals = [_component_from_smiles(s) for s in smiles_tuple]
    return AIOMFAC(neutrals, salts=list(salts_tuple) if salts_tuple else None,
                   dissociate_bisulfate=dissociate_bisulfate)


@functools.lru_cache(maxsize=128)
def _cached_aiomfac_sr_system(smiles_tuple):
    components = [_component_from_smiles(s) for s in smiles_tuple]
    return AIOMFAC_SR(components)


def clear_system_cache():
    """Drop the cached AIOMFAC and AIOMFAC_SR systems built from
    ``aiomfac_from_smiles`` / ``aiomfac_sr_from_smiles`` calls. Useful in
    long-running workloads that build systems from an unbounded number of
    distinct SMILES lists."""
    _cached_aiomfac_system.cache_clear()
    _cached_aiomfac_sr_system.cache_clear()


def system_cache_info():
    """Return ``(aiomfac, aiomfac_sr)`` :class:`functools._CacheInfo`
    tuples for the two cached system constructors."""
    return (
        _cached_aiomfac_system.cache_info(),
        _cached_aiomfac_sr_system.cache_info(),
    )


def aiomfac_sr_from_smiles(smiles, mole_fractions, T_K, strict='warn',
                           return_diagnostics=False):
    """SMILES-friendly entry point for the short-range AIOMFAC calculation.

    Decomposes each SMILES via :func:`aerosysprop.groups.aiomfac` (which now
    correctly emits subgroup 153 for the alcohol-OH convention), builds a
    fresh :class:`AIOMFAC_SR` system, and returns the activity coefficients
    at the requested composition and temperature. For repeated evaluation
    at many compositions, instantiate :class:`AIOMFAC_SR` once and reuse it
    rather than calling this function in a loop.

    Parameters
    ----------
    smiles : sequence of str
        SMILES strings for each non-electrolyte component, in mole-fraction
        order. The first component is conventionally water (``'O'``) but
        any neutral organic-soluble mixture works.
    mole_fractions : sequence of float
        Mole fractions, same order and length as ``smiles``. Does not need
        to sum exactly to 1 (UNIFAC is a function of x, not normalised x).
    T_K : float
        Temperature in Kelvin.

    Returns
    -------
    numpy.ndarray
        Activity coefficients gamma_i on a mole-fraction basis with
        pure-component reference state, shape ``(len(smiles),)``.

    Raises
    ------
    ValueError
        If any SMILES is malformed (RDKit can't parse it).

    Notes
    -----
    Some compounds with heteroatom groups outside the AIOMFAC SMARTS set
    (e.g., -NO2, -SH, -ONO2, hydroperoxide -OOH, peroxynitrate) are
    silently decomposed only for the parts the matcher recognises. The
    resulting gamma reflects only the recognised functional groups.
    See ``activity_coefficients.md`` (Known limitations) for the full list
    of compounds affected.

    Examples
    --------
    >>> # water + ethanol at 50:50, 25 C
    >>> aiomfac_sr_from_smiles(['O', 'CCO'], [0.5, 0.5], 298.15)
    array([1.26196..., 1.18692...])

    >>> # Ternary water + glycerol + glucose at a humid SOA-relevant point
    >>> aiomfac_sr_from_smiles(
    ...     ['O', 'OCC(O)CO', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'],
    ...     [0.7, 0.2, 0.1], 298.15,
    ... )  # doctest: +SKIP

    The ``strict`` argument controls how composition / matcher issues are
    surfaced (see :func:`aerosysprop.diagnostics.diagnose_aiomfac_sr`):

    - ``'warn'`` (default): emits ``warnings.warn`` for any issue, then
      proceeds with the calculation.
    - ``'raise'``: raises ``ValueError`` if any issue is found.
    - ``'silent'``: no diagnostics, original behaviour.

    Set ``return_diagnostics=True`` to receive a ``(gamma, report)``
    tuple instead of the bare ``gamma`` ndarray; the report is the same
    :class:`aerosysprop.diagnostics.DiagnosticReport` that drives the
    ``strict`` behaviour and can be inspected after the call.
    """
    report = None
    if strict != 'silent' or return_diagnostics:
        from aerosysprop.diagnostics import diagnose_aiomfac_sr, apply_strict
        report = diagnose_aiomfac_sr(smiles, mole_fractions, T_K)
        if strict != 'silent':
            apply_strict(report, strict)
    system = _cached_aiomfac_sr_system(tuple(smiles))
    gamma = system.gamma(mole_fractions, T_K)
    if return_diagnostics:
        return gamma, report
    return gamma


# ============================================================================
# Phase 2: Electrolyte support (MR + LR Debye-Hueckel) for salt-containing
# mixtures. Implementation follows ModMRpart.f90:LR_MR_activity in the public
# AIOMFAC repo; see activity_coefficients.md for theory and citations.
# ============================================================================

from collections import namedtuple

# Hardcoded T-dependent Debye-Hueckel reference constants for water solvent.
# Same as Fortran ModMRpart.f90 lines 109-114.
_M_WATER     = 0.01801528           # [kg/mol] molar mass of water
_DH_DENS_W   = 997.0                # [kg/m^3] density of water at 298.15 K
_DH_DIEL_W   = 78.54                # static dielectric constant at 298.15 K
_DH_A_COEFF  = 1.327757e5           # for A_DebyeHw
_DH_B_COEFF  = 6.359696             # for B_DebyeHw

# Aggregated activity-coefficient result for a full AIOMFAC calculation.
# The `diagnostics` field is populated by `aiomfac_from_smiles` when a
# pre-flight diagnostic report was generated (i.e. strict != 'silent');
# it is None when the calculation was driven through `AIOMFAC.gamma`
# directly (where SMILES decomposition didn't happen) or when the user
# opted out via strict='silent'. The field defaults to None for
# backward compatibility with code that constructs AIOMFACResult
# positionally.
AIOMFACResult = namedtuple(
    'AIOMFACResult',
    ['gamma_neutral', 'gamma_cation', 'gamma_anion', 'gamma_mean',
     'cation_ids', 'anion_ids', 'salt_names', 'diagnostics'],
    defaults=[None],
)


def _lookup_salt(name):
    """Look up a salt by name and return (salt_id, cation_id, n_cation, anion_id,
    n_anion, molar_mass_kg_per_mol). Raises KeyError if not found."""
    matches = [
        sid for sid, sname in data.AIOMFAC_SALT_NAME.items() if sname == name
    ]
    if not matches:
        raise KeyError(
            f'salt {name!r} not in AIOMFAC catalog; known salts include '
            + ', '.join(sorted(set(data.AIOMFAC_SALT_NAME.values()))[:8]) + ', ...'
        )
    sid = matches[0]
    return (
        sid,
        int(data.AIOMFAC_SALT_CATION_GROUP[sid]),
        float(data.AIOMFAC_SALT_CATION_STOICH[sid]),
        int(data.AIOMFAC_SALT_ANION_GROUP[sid]),
        float(data.AIOMFAC_SALT_ANION_STOICH[sid]),
        float(data.AIOMFAC_SALT_MASS[sid]) * 1.0e-3,   # g/mol -> kg/mol
    )


def _component_molar_mass(subgroups):
    """Molar mass [kg/mol] of a neutral component from its subgroup decomposition."""
    return sum(count * data.AIOMFAC_MASS[sg] for sg, count in subgroups.items())


# ---------------------------------------------------------------------------
# Bisulfate dissociation equilibrium  HSO4- <-> H+ + SO4--
# ---------------------------------------------------------------------------
# Temperature-dependent dissociation constant K_HSO4(T) from the Knopf et al.
# (2003) parameterization, exactly as in Fortran AIOMFAC
# (SubModDissociationEquil.f90). Valid ~179-474 K.
_BS_T0    = 298.15
_BS_DH0   = -1.8554748598e4
_BS_CP0   = -2.05594443486e2
_BS_DCPDT = -9.94992864240e-1
_BS_FE1   = 1.0 / 8.314472
_BS_FE2   = _BS_DH0 - _BS_CP0 * _BS_T0 + 0.5 * _BS_DCPDT * _BS_T0 * _BS_T0
_BS_FE3   = 1.0 / _BS_T0
_BS_FE4   = _BS_CP0 - _BS_DCPDT * _BS_T0
# AIOMFAC ion ids involved in the equilibrium.
_ION_H, _ION_HSO4, _ION_SO4 = 205, 248, 261


def _lnK_bisulfate(T_K):
    """``ln`` of the HSO4- <-> H+ + SO4-- dissociation constant at ``T_K``
    (Knopf et al. 2003; matches Fortran AIOMFAC). Clamped to the valid T range.
    At 298.15 K this returns -4.54916799587 (K ~= 0.01057)."""
    T = min(max(float(T_K), 180.0), 473.0)
    return (
        -4.54916799587
        - _BS_FE1 * (
            _BS_FE2 * (1.0 / T - _BS_FE3)
            - _BS_FE4 * np.log(T / _BS_T0)
            - 0.5 * _BS_DCPDT * (T - _BS_T0)
        )
    )


class AIOMFAC:
    """Full AIOMFAC activity-coefficient calculation: SR + MR + LR.

    Handles neutral organic + water + salt mixtures. For non-electrolyte
    mixtures, MR and LR are identically zero and this reduces to the SR-only
    case (see :class:`AIOMFAC_SR`).

    Parameters
    ----------
    neutral_components : sequence of mapping ``{int: number}``
        One subgroup-count dict per neutral component, in the order the
        mass-fraction vector will use. Water should be the first (canonical
        Fortran convention).
    salts : sequence of str, optional
        Salt names as listed in ``data.AIOMFAC_SALT_NAME`` (e.g., 'NaCl',
        'Na2SO4', 'NH4NO3'). Looked up against ``aiomfac_salts.data``.

    Example
    -------
    >>> sys = AIOMFAC(
    ...     neutral_components=[{16: 1}],         # water
    ...     salts=['NaCl'],
    ... )
    >>> result = sys.gamma(mass_fractions=[0.95, 0.05], T_K=298.15)
    >>> result.gamma_neutral   # gamma for water
    array([...])
    >>> result.gamma_mean      # mean ionic gamma_+- for NaCl
    array([...])
    """

    def __init__(self, neutral_components, salts=None, dissociate_bisulfate=True):
        self.neutral_components = [dict(c) for c in neutral_components]
        self.n_neutral = len(self.neutral_components)
        self.salt_names = list(salts) if salts else []

        # Look up each salt -> (cation_id, n_cation, anion_id, n_anion, M_salt)
        self.salt_records = [_lookup_salt(name) for name in self.salt_names]

        # Distinct cation/anion IDs in first-encountered order across the
        # user-declared salt list. This matches the Fortran output convention
        # (ions reported in salt-declaration order) so result indices stay
        # aligned with the user's salt input order.
        self.cation_ids, self.anion_ids = [], []
        _seen_c, _seen_a = set(), set()
        for rec in self.salt_records:
            cid, aid = rec[1], rec[3]
            if cid not in _seen_c:
                self.cation_ids.append(cid); _seen_c.add(cid)
            if aid not in _seen_a:
                self.anion_ids.append(aid); _seen_a.add(aid)
        self.n_cation = len(self.cation_ids)
        self.n_anion  = len(self.anion_ids)
        self._cation_idx = {cid: i for i, cid in enumerate(self.cation_ids)}
        self._anion_idx  = {aid: i for i, aid in enumerate(self.anion_ids)}

        # Bisulfate equilibrium: when H+ (205) and SO4-- (261) coexist, the model
        # must also carry HSO4- (248) so the equilibrium HSO4- <-> H+ + SO4-- can
        # be solved (HSO4- is already fully parameterized in the ion data). Add it
        # to the anion set; the per-trial molality split is solved in gamma().
        self.dissociate_bisulfate = bool(
            dissociate_bisulfate
            and _ION_H in self._cation_idx
            and _ION_SO4 in self._anion_idx
        )
        if self.dissociate_bisulfate and _ION_HSO4 not in self._anion_idx:
            self.anion_ids.append(_ION_HSO4)
            self._anion_idx[_ION_HSO4] = len(self.anion_ids) - 1
            self.n_anion = len(self.anion_ids)
        if self.dissociate_bisulfate:
            self._idx_H    = self._cation_idx[_ION_H]
            self._idx_HSO4 = self._anion_idx[_ION_HSO4]
            self._idx_SO4  = self._anion_idx[_ION_SO4]

        # Charges and absolute charges, per cation/anion index
        self.ZC = np.array(
            [data.AIOMFAC_ION_CHARGE[cid] for cid in self.cation_ids],
            dtype=np.float64,
        )
        self.ZA = np.array(
            [abs(data.AIOMFAC_ION_CHARGE[aid]) for aid in self.anion_ids],
            dtype=np.float64,
        )

        # Molar masses [kg/mol] in component order: neutrals then salts
        self.Mmass = np.array(
            [_component_molar_mass(c) for c in self.neutral_components]
            + [rec[5] for rec in self.salt_records],
            dtype=np.float64,
        )

        # Build the SR sub-system, treating each ion as its own "component"
        # with a single subgroup. SR's n_neutral is the count of leading
        # neutrals so the ion combinatorial correction kicks in.
        sr_components = (
            list(self.neutral_components)
            + [{cid: 1} for cid in self.cation_ids]
            + [{aid: 1} for aid in self.anion_ids]
        )
        self.sr = AIOMFAC_SR(sr_components, n_neutral=self.n_neutral)

        # Pre-fetch MR interaction matrices indexed by main_group or ion_id
        self._build_mr_tables()

    def _build_mr_tables(self):
        """Build dense numpy arrays for the MR ion-ion and org-ion params,
        sized to the cations/anions present in this system."""
        nc = self.n_cation
        na = self.n_anion

        # Ion-ion binary: bac, cac, omega; cn1, cn2, omega2 (shape nc x na)
        self.bac    = np.zeros((nc, na), dtype=np.float64)
        self.cac    = np.zeros((nc, na), dtype=np.float64)
        self.omega  = np.full((nc, na), 0.8, dtype=np.float64)
        self.cnac1  = np.zeros((nc, na), dtype=np.float64)
        self.cnac2  = np.zeros((nc, na), dtype=np.float64)
        self.omega2 = np.full((nc, na), 0.6, dtype=np.float64)
        for i, cid in enumerate(self.cation_ids):
            for j, aid in enumerate(self.anion_ids):
                self.bac[i, j]    = data.AIOMFAC_MR_ION_ION_INTERACTIONS_b1[cid].get(aid, 0.0)
                self.cac[i, j]    = data.AIOMFAC_MR_ION_ION_INTERACTIONS_b2[cid].get(aid, 0.0)
                self.cnac1[i, j]  = data.AIOMFAC_MR_ION_ION_INTERACTIONS_b3[cid].get(aid, 0.0)
                self.cnac2[i, j]  = data.AIOMFAC_MR_ION_ION_INTERACTIONS_c2[cid].get(aid, 0.0)
                # c1 = omegaTAB (default 0.8); c3 = omega2TAB (default 0.6)
                self.omega[i, j]  = data.AIOMFAC_MR_ION_ION_INTERACTIONS_c1[cid].get(aid, 0.8)
                self.omega2[i, j] = data.AIOMFAC_MR_ION_ION_INTERACTIONS_c3[cid].get(aid, 0.6)

        # Org-ion: bnc, cnc indexed by main_group present in the neutrals
        # We collect the unique main groups for the neutral components.
        neutral_main_groups = []
        # Reuse the SR sub-system's subgroup -> main-group mapping for the
        # leading neutral components.
        for i in range(self.n_neutral):
            for sg, ct in self.neutral_components[i].items():
                if ct == 0:
                    continue
                mg = int(data.AIOMFAC_MAIN_GROUP[sg])
                if mg not in neutral_main_groups:
                    neutral_main_groups.append(mg)
        self.main_groups = neutral_main_groups
        nmg = len(self.main_groups)

        self.bnc = np.zeros((nmg, nc), dtype=np.float64)
        self.cnc = np.zeros((nmg, nc), dtype=np.float64)
        self.bna = np.zeros((nmg, na), dtype=np.float64)
        self.cna = np.zeros((nmg, na), dtype=np.float64)
        # Detect unparameterised pairs (Fortran's -7.7777e4 sentinel) and
        # raise a clear error rather than silently using 0. Matches the
        # Fortran AIOMFAC error 1/2 behaviour in MRinteractcoeff.
        SENTINEL = -1000.0  # Fortran flags anything < -1000 as "not estimated"
        missing = []
        for k, mg in enumerate(self.main_groups):
            row_b = data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b1.get(mg, {})
            row_c = data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b2.get(mg, {})
            for i, cid in enumerate(self.cation_ids):
                b = row_b.get(cid, 0.0); c = row_c.get(cid, 0.0)
                if b < SENTINEL or c < SENTINEL:
                    missing.append((mg, cid, 'cation'))
                self.bnc[k, i] = b
                self.cnc[k, i] = c
            for j, aid in enumerate(self.anion_ids):
                b = row_b.get(aid, 0.0); c = row_c.get(aid, 0.0)
                if b < SENTINEL or c < SENTINEL:
                    missing.append((mg, aid, 'anion'))
                self.bna[k, j] = b
                self.cna[k, j] = c
        if missing:
            details = ', '.join(
                f'main_group {mg} <-> ion {ion} ({kind})'
                for mg, ion, kind in missing
            )
            raise ValueError(
                'AIOMFAC: required organic-ion MR interaction parameter '
                f'not estimated in current AIOMFAC parameter set: {details}. '
                'These (main group, ion) pairs are flagged with the '
                '"-7.7777e4" sentinel in the Fortran source as well — Fortran '
                'AIOMFAC also rejects such mixtures (error code 1 or 2). '
                'Consider substituting an organic that uses parameterised '
                'main groups, or wait for an AIOMFAC parameter-set update.'
            )

        # Subgroup mole-fraction -> main-group mole-fraction projection:
        # for each neutral component, build vector of (main_group_count) so
        # later we can compute mgroupxsf = sum_i ITAB_mg(i) * XN(i).
        # ITAB_mg[i, k] = count of main_group k atoms in neutral component i.
        self._ITAB_mg = np.zeros((self.n_neutral, nmg), dtype=np.float64)
        # SubGroupMW per (neutral_i, main_group_k): mass-weighted sum
        self._subgroup_mass_mg = np.zeros((self.n_neutral, nmg), dtype=np.float64)
        # We also need the molar mass of each main group as it appears in this
        # mixture (Fortran mgroupMK, computed from subgroup composition).
        mg_to_subgroups = {mg: [] for mg in self.main_groups}
        for i, c in enumerate(self.neutral_components):
            for sg, ct in c.items():
                if ct == 0:
                    continue
                mg = int(data.AIOMFAC_MAIN_GROUP[sg])
                k = self.main_groups.index(mg)
                self._ITAB_mg[i, k] += ct
                self._subgroup_mass_mg[i, k] += ct * data.AIOMFAC_MASS[sg]
                if sg not in mg_to_subgroups[mg]:
                    mg_to_subgroups[mg].append(sg)

    # ---- composition & ionic-strength helpers --------------------------------
    def _composition(self, mass_fractions):
        """Return all derived composition quantities needed by SR + MR + LR.

        Returns dict with keys: XN, mNeutral, SMC, SMA, X_sr, meanSolventMW,
        SumIonMolalities, SI, Z, mgroupxsf, mgroupMK, avgmaingrMW.
        """
        wtf = np.asarray(mass_fractions, dtype=np.float64)
        if wtf.shape != (self.n_neutral + len(self.salt_records),):
            raise ValueError(
                f'mass_fractions length {wtf.shape} != n_neutral '
                f'({self.n_neutral}) + n_salts ({len(self.salt_records)})'
            )

        # Salt-free neutral mole fractions XN (Fortran ModCalcActCoeff.f90:283-285)
        wtfbyMW = wtf[: self.n_neutral] / self.Mmass[: self.n_neutral]
        sumXN = wtfbyMW.sum()
        XN = wtfbyMW / sumXN if sumXN > 0 else np.zeros(self.n_neutral)
        meanSolventMW = float((self.Mmass[: self.n_neutral] * XN).sum())
        mNeutral = XN / meanSolventMW if meanSolventMW > 0 else XN

        # Ion molalities (Fortran MassFrac2IonMolalities)
        SMC = np.zeros(self.n_cation, dtype=np.float64)
        SMA = np.zeros(self.n_anion, dtype=np.float64)
        sumWN = wtf[: self.n_neutral].sum()
        for k, rec in enumerate(self.salt_records):
            _, cid, n_c, aid, n_a, M_salt = rec
            comp_idx = self.n_neutral + k
            mol_per_kg = wtf[comp_idx] / M_salt / sumWN if sumWN > 0 else 0.0
            SMC[self._cation_idx[cid]] += mol_per_kg * n_c
            SMA[self._anion_idx[aid]]  += mol_per_kg * n_a
        SumIonMolalities = SMC.sum() + SMA.sum()

        # Ionic strength and net charge density
        SI = 0.5 * (SMA * self.ZA ** 2).sum() + 0.5 * (SMC * self.ZC ** 2).sum()
        Z  = (SMA * self.ZA).sum() + (SMC * self.ZC).sum()

        # Dissociated mole fractions X for SR (neutrals first, then cations,
        # then anions) — Fortran ModCalcActCoeff.f90:300-303
        sum_molalities = mNeutral.sum() + SumIonMolalities
        X_sr = np.concatenate([mNeutral, SMC, SMA]) / sum_molalities

        # Main-group mole fractions on a salt-free basis (Fortran lines 121-148)
        # mgroupxsf[k] = sum_i ITAB_mg[i, k] * XN[i]  (then normalised)
        if self.main_groups:
            mg_xsf_raw = self._ITAB_mg.T @ XN     # [n_mg]
            sumxg = mg_xsf_raw.sum()
            mgroupxsf = mg_xsf_raw / sumxg if sumxg > 0 else mg_xsf_raw
            # mgroupMK[k] = sum_i (Mass_mg[i,k] / mg_xsf_raw[k]) * XN[i]
            # but Fortran normalises differently; we follow Fortran's formula:
            # mgroupMK[k] = sum_i SubGroupMW[k] * subgroupxsf[i,k] / mgroupxsf[k]
            # Equivalent for us: mgroupMK = mass-weighted average per main group.
            with np.errstate(divide='ignore', invalid='ignore'):
                mgroupMK = np.where(
                    mgroupxsf > 0,
                    (self._subgroup_mass_mg.T @ XN) / np.where(mg_xsf_raw > 0, mg_xsf_raw, 1.0),
                    0.0,
                )
            avgmaingrMW = float((mgroupMK * mgroupxsf).sum())
        else:
            mgroupxsf = np.zeros(0)
            mgroupMK  = np.zeros(0)
            avgmaingrMW = 0.0

        return {
            'XN': XN, 'mNeutral': mNeutral, 'SMC': SMC, 'SMA': SMA,
            'X_sr': X_sr, 'meanSolventMW': meanSolventMW,
            'SumIonMolalities': SumIonMolalities, 'SI': SI, 'Z': Z,
            'mgroupxsf': mgroupxsf, 'mgroupMK': mgroupMK,
            'avgmaingrMW': avgmaingrMW,
        }

    # ---- LR (Debye-Hueckel) -------------------------------------------------
    def _ln_gamma_LR(self, T_K, SI):
        """Long-range Debye-Hueckel contributions, mole-fraction basis.

        Returns (gnlrln, gclrln, galrln) arrays for neutrals, cations, anions.
        Follows Fortran ModMRpart.f90:186-195.
        """
        sqrtDens = _DH_DENS_W ** 0.5
        diel = _DH_DIEL_W
        A = _DH_A_COEFF * sqrtDens / (diel * T_K) ** 1.5
        B = _DH_B_COEFF * sqrtDens / (diel * T_K) ** 0.5
        SI2 = SI ** 0.5
        bb  = 1.0 + B * SI2
        TermI = bb - 1.0 / bb - 2.0 * np.log(bb) if SI > 0 else 0.0

        # Neutrals: gnlrln = Mmass * (2 A / B^3) * TermI
        coeff_neutral = (2.0 * A / B ** 3) * TermI
        gnlrln = self.Mmass[: self.n_neutral] * coeff_neutral

        # Ions: gclrln = -Z_c^2 * A * sqrt(I) / (1 + B sqrt(I)); same for anions
        if SI > 0:
            lrw1 = A * SI2 / bb
        else:
            lrw1 = 0.0
        gclrln = -self.ZC ** 2 * lrw1
        galrln = -self.ZA ** 2 * lrw1
        return gnlrln, gclrln, galrln

    # ---- MR (organic-ion + ion-ion) ------------------------------------------
    def _ln_gamma_MR(self, comp):
        """Middle-range contributions, mole-fraction basis. Returns (gnmrln,
        gcmrln, gamrln). Follows Fortran ModMRpart.f90:200-372 (excluding the
        triple-ion Qcca and cation-cation Rcc terms)."""
        SI = comp['SI']
        if SI <= 0 or not self.salt_records:
            return (
                np.zeros(self.n_neutral),
                np.zeros(self.n_cation),
                np.zeros(self.n_anion),
            )

        SI2 = SI ** 0.5
        SMC = comp['SMC']
        SMA = comp['SMA']
        Z = comp['Z']
        mgroupxsf = comp['mgroupxsf']
        mgroupMK  = comp['mgroupMK']
        avgmaingrMW = comp['avgmaingrMW']

        # Ion-ion B and Cn (with ionic-strength scaling, Fortran lines 205-217)
        Bca  = self.bac   + self.cac   * np.exp(-self.omega  * SI2)
        Cnca = self.cnac1 + self.cnac2 * np.exp(-self.omega2 * SI2)
        # Derivatives w.r.t. ionic strength (lines 218-220):
        # BSca = (-0.5/SI2) * omega  * (Bca  - bac)
        # CSca = (-0.5/SI2) * omega2 * (Cnca - cnac1)
        lrw1 = -0.5 / SI2 if SI2 > 0 else 0.0
        BSca = lrw1 * self.omega  * (Bca  - self.bac)
        CSca = lrw1 * self.omega2 * (Cnca - self.cnac1)

        # Org-ion B (lines 223-234)
        Bkionc = self.bnc + self.cnc * np.exp(-1.2 * SI2)
        Bkiona = self.bna + self.cna * np.exp(-1.2 * SI2)
        BSkionc = lrw1 * 1.2 * (Bkionc - self.bnc)
        BSkiona = lrw1 * 1.2 * (Bkiona - self.bna)

        # ---- Neutral (main-group level then back to molecules) --------------
        # SumBm[k] = sum_j Bkiona[k,j] * SMA[j] + Bkionc[k,j] * SMC[j]
        SumBm = Bkiona @ SMA + Bkionc @ SMC

        # Sumkion1 = sum over ions and main groups of
        #   (B + SI*BS) * mgroupxsf * SM
        if mgroupxsf.size:
            Sumkion1 = (
                ((Bkiona + SI * BSkiona) * mgroupxsf[:, None] * SMA[None, :]).sum()
                + ((Bkionc + SI * BSkionc) * mgroupxsf[:, None] * SMC[None, :]).sum()
            )
        else:
            Sumkion1 = 0.0

        # SumCA1 = sum_{i,j} (Bca + SI*BSca + 2Z*Cnca + Z*SI*Csca) * SMC[i]*SMA[j]
        twoZ = 2.0 * Z
        ZSI = Z * SI
        SumCA1 = (
            (Bca + SI * BSca + twoZ * Cnca + ZSI * CSca)
            * SMC[:, None] * SMA[None, :]
        ).sum()

        # gkmrln per main group, then mapped back to neutral component
        if mgroupxsf.size:
            gkmrln = SumBm - mgroupMK / avgmaingrMW * Sumkion1 - mgroupMK * SumCA1
        else:
            gkmrln = np.zeros(0)
        # gnmrln[i] = sum_k ITAB_mg[i, k] * gkmrln[k]
        gnmrln = self._ITAB_mg @ gkmrln if gkmrln.size else np.zeros(self.n_neutral)

        # ---- Cations (lines 297-322) ----------------------------------------
        if mgroupxsf.size:
            SumBcx = Bkionc.T @ mgroupxsf
            SumBax = Bkiona.T @ mgroupxsf
        else:
            SumBcx = np.zeros(self.n_cation)
            SumBax = np.zeros(self.n_anion)

        # SumCBca[i] = sum_j (Bca[i,j] + Z*Cnca[i,j]) * SMA[j]
        SumCBca = ((Bca + Z * Cnca) * SMA[None, :]).sum(axis=1)
        # SumABca[i] = sum_j (Bca[j,i] + Z*Cnca[j,i]) * SMC[j]
        SumABca = ((Bca + Z * Cnca) * SMC[:, None]).sum(axis=0)

        smcsma = SMC[:, None] * SMA[None, :]
        SumBsCA = (BSca * smcsma).sum()
        SumCnca = (Cnca * smcsma).sum()
        SumCsca = (CSca * smcsma).sum()

        if mgroupxsf.size:
            Sumkion2 = (
                (BSkiona * mgroupxsf[:, None] * SMA[None, :]).sum()
                + (BSkionc * mgroupxsf[:, None] * SMC[None, :]).sum()
            )
        else:
            Sumkion2 = 0.0

        if avgmaingrMW > 0:
            gcmrln = (
                SumBcx / avgmaingrMW
                + self.ZC ** 2 / (2.0 * avgmaingrMW) * Sumkion2
                + SumCBca
                + 0.5 * self.ZC ** 2 * SumBsCA
                + self.ZC * SumCnca
                + 0.5 * self.ZC ** 2 * Z * SumCsca
            )
            gamrln = (
                SumBax / avgmaingrMW
                + self.ZA ** 2 / (2.0 * avgmaingrMW) * Sumkion2
                + SumABca
                + 0.5 * self.ZA ** 2 * SumBsCA
                + self.ZA * SumCnca
                + 0.5 * self.ZA ** 2 * Z * SumCsca
            )
        else:
            gcmrln = SumCBca + 0.5 * self.ZC ** 2 * SumBsCA + self.ZC * SumCnca
            gamrln = SumABca + 0.5 * self.ZA ** 2 * SumBsCA + self.ZA * SumCnca

        return gnmrln, gcmrln, gamrln

    # ---- combined contributions & bisulfate equilibrium ---------------------
    def _ln_gamma_all(self, comp, T_K):
        """Combined SR+LR+MR ln activity coefficients for ``comp``: returns
        (neutral [mole-fraction basis], cation, anion [both molality basis])."""
        ln_gamma_SR = self.sr.ln_gamma(comp['X_sr'], T_K)
        gn_sr = ln_gamma_SR[: self.n_neutral]
        gc_sr = ln_gamma_SR[self.n_neutral : self.n_neutral + self.n_cation]
        ga_sr = ln_gamma_SR[self.n_neutral + self.n_cation :]
        gn_lr, gc_lr, ga_lr = self._ln_gamma_LR(T_K, comp['SI'])
        gn_mr, gc_mr, ga_mr = self._ln_gamma_MR(comp)
        ln_gamma_neutral = gn_sr + gn_lr + gn_mr
        # Ions: subtract Tmolal to convert mole-fraction basis -> molality.
        Tmolal = np.log(
            _M_WATER / comp['meanSolventMW']
            + _M_WATER * comp['SumIonMolalities']
        )
        ln_gamma_cation = gc_sr + gc_lr + gc_mr - Tmolal
        ln_gamma_anion  = ga_sr + ga_lr + ga_mr - Tmolal
        return ln_gamma_neutral, ln_gamma_cation, ln_gamma_anion

    def _update_ionic(self, comp, SMC, SMA):
        """Return a composition dict with new ion molalities and the derived
        ionic-strength / SR mole-fraction quantities recomputed (the neutral-only
        quantities are unchanged by the bisulfate split)."""
        c2 = dict(comp)
        c2['SMC'] = SMC
        c2['SMA'] = SMA
        c2['SI'] = 0.5 * (SMA * self.ZA ** 2).sum() + 0.5 * (SMC * self.ZC ** 2).sum()
        c2['Z']  = (SMA * self.ZA).sum() + (SMC * self.ZC).sum()
        c2['SumIonMolalities'] = SMC.sum() + SMA.sum()
        sm = comp['mNeutral'].sum() + c2['SumIonMolalities']
        c2['X_sr'] = np.concatenate([comp['mNeutral'], SMC, SMA]) / sm
        return c2

    def _solve_bisulfate(self, comp, T_K):
        """Partition sulfate between SO4-- and HSO4- by solving the equilibrium
        HSO4- <-> H+ + SO4-- self-consistently with the activity coefficients.

        ``comp`` enters as the fully-dissociated reference (HSO4- molality 0).
        Returns a new composition dict at the equilibrium ion molalities, found by
        driving ``ln(m_H m_SO4 / m_HSO4) + (lnγ_H + lnγ_SO4 - lnγ_HSO4) = lnK``.
        Conservation: m_H = H_tot - x, m_SO4 = SO4_tot - x, m_HSO4 = x.
        """
        iH, iHSO4, iSO4 = self._idx_H, self._idx_HSO4, self._idx_SO4
        SMC0, SMA0 = comp['SMC'], comp['SMA']
        H_tot   = SMC0[iH]   + SMA0[iHSO4]      # proton balance
        SO4_tot = SMA0[iSO4] + SMA0[iHSO4]      # total sulfate
        if H_tot <= 0.0 or SO4_tot <= 0.0:
            return comp
        lnK = _lnK_bisulfate(T_K)
        xmax = min(H_tot, SO4_tot)

        def diffK(x):
            SMC = SMC0.copy()
            SMA = SMA0.copy()
            SMC[iH]    = H_tot - x
            SMA[iSO4]  = SO4_tot - x
            SMA[iHSO4] = x
            _, lc, la = self._ln_gamma_all(self._update_ionic(comp, SMC, SMA), T_K)
            return (
                (lc[iH] + la[iSO4] - la[iHSO4])
                + np.log((H_tot - x) * (SO4_tot - x) / x)
                - lnK
            )

        # As x->0+ the log term -> +inf (f>0); as x->xmax- -> -inf (f<0): a root
        # is bracketed. Bisection (no scipy dependency).
        lo, hi = xmax * 1e-12, xmax * (1.0 - 1e-12)
        flo, fhi = diffK(lo), diffK(hi)
        if not (np.isfinite(flo) and np.isfinite(fhi)) or flo * fhi > 0.0:
            # No bracketed sign change (or non-finite activity coefficients from a
            # degenerate composition) -> keep the fully-dissociated reference
            # rather than risk a meaningless root. Bisection below is otherwise
            # unconditionally convergent.
            return comp
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            fm = diffK(mid)
            if fm == 0.0 or (hi - lo) <= xmax * 1e-13:
                break
            if (fm > 0.0) == (flo > 0.0):
                lo, flo = mid, fm
            else:
                hi = mid
        x = 0.5 * (lo + hi)
        SMC = SMC0.copy()
        SMA = SMA0.copy()
        SMC[iH]    = H_tot - x
        SMA[iSO4]  = SO4_tot - x
        SMA[iHSO4] = x
        return self._update_ionic(comp, SMC, SMA)

    # ---- public API ----------------------------------------------------------
    def gamma(self, mass_fractions, T_K):
        """Return activity coefficients for the configured system.

        Parameters
        ----------
        mass_fractions : sequence of float
            Mass fractions of each component, in the order
            [neutral_1, neutral_2, ..., salt_1, salt_2, ...]. Should sum to 1.
        T_K : float
            Temperature in Kelvin.

        Returns
        -------
        AIOMFACResult
            Named tuple with fields:
              - gamma_neutral: mole-fraction-basis gamma for each neutral
              - gamma_cation:  molality-basis gamma for each distinct cation
              - gamma_anion:   molality-basis gamma for each distinct anion
              - gamma_mean:    mean ionic gamma_+- for each salt
              - cation_ids, anion_ids, salt_names: parallel lookup lists
        """
        comp = self._composition(mass_fractions)
        T_K = float(T_K)

        # Solve the bisulfate equilibrium (updating the ion molalities) before
        # evaluating the final activity coefficients, when applicable.
        if self.dissociate_bisulfate:
            comp = self._solve_bisulfate(comp, T_K)

        # Combined SR + LR + MR (Fortran ModCalcActCoeff.f90:168-176).
        ln_gamma_neutral, ln_gamma_cation, ln_gamma_anion = \
            self._ln_gamma_all(comp, T_K)

        gamma_n = np.exp(ln_gamma_neutral)
        gamma_c = np.exp(ln_gamma_cation)
        gamma_a = np.exp(ln_gamma_anion)

        # Mean ionic activity coefficient per salt:
        # gamma_+- = (gamma_c^nu_c * gamma_a^nu_a) ^ (1 / (nu_c + nu_a))
        gamma_mean = np.zeros(len(self.salt_records))
        for k, (sid, cid, n_c, aid, n_a, _M) in enumerate(self.salt_records):
            gc = gamma_c[self._cation_idx[cid]]
            ga = gamma_a[self._anion_idx[aid]]
            gamma_mean[k] = (gc ** n_c * ga ** n_a) ** (1.0 / (n_c + n_a))

        return AIOMFACResult(
            gamma_neutral=gamma_n,
            gamma_cation=gamma_c,
            gamma_anion=gamma_a,
            gamma_mean=gamma_mean,
            cation_ids=list(self.cation_ids),
            anion_ids=list(self.anion_ids),
            salt_names=list(self.salt_names),
        )


def aiomfac_from_smiles(smiles, salts, mass_fractions, T_K, strict='warn',
                        bisulfate=True):
    """Full AIOMFAC entry point taking SMILES + salt names + mass fractions.

    Wraps :class:`AIOMFAC` with SMILES-based neutral-component decomposition.

    Parameters
    ----------
    smiles : sequence of str
        SMILES strings for each neutral component, in mass-fraction order.
    salts : sequence of str
        Salt names from the AIOMFAC catalog (e.g., 'NaCl', 'Na2SO4').
    mass_fractions : sequence of float
        Mass fractions of all components: [neutral_1, ..., salt_1, ...].
        Should sum to 1.
    T_K : float
        Temperature in Kelvin.

    strict : {'warn', 'raise', 'silent'}, default 'warn'
        How to surface composition / matcher issues. ``'warn'`` emits
        ``warnings.warn`` for each issue and proceeds; ``'raise'`` raises
        ``ValueError`` if any issue is detected; ``'silent'`` skips the
        pre-flight check. See :mod:`aerosysprop.diagnostics`.
    bisulfate : bool, default True
        When True and both H+ and SO4-- are present (e.g. H2SO4), solve the
        bisulfate dissociation equilibrium HSO4- <-> H+ + SO4-- so the results
        include HSO4- with its equilibrium molality (matching Fortran AIOMFAC's
        treatment). Set False to force full dissociation to 2 H+ + SO4--.

    Returns
    -------
    AIOMFACResult
        See :class:`AIOMFAC.gamma` for fields.
    """
    report = None
    if strict != 'silent':
        from aerosysprop.diagnostics import diagnose_aiomfac, apply_strict
        report = diagnose_aiomfac(smiles, salts, mass_fractions, T_K)
        apply_strict(report, strict)
    system = _cached_aiomfac_system(tuple(smiles), tuple(salts) if salts else (),
                                    bool(bisulfate))
    result = system.gamma(mass_fractions, T_K)
    # Attach the pre-flight diagnostic report so users can introspect issues
    # alongside the numeric outputs (e.g. ``result.diagnostics.warnings``).
    return result._replace(diagnostics=report)
