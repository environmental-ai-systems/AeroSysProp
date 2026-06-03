# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2026 David Topping.
# All Rights Reserved.
# This file is part of aerosysprop.
"""Diagnostics framework for AeroSysProp predictive methods.

This module lets you screen SMILES strings and mixture compositions for
known issues **before** calling any aerosysprop method, OR have those
issues surfaced as warnings during a normal call. It catches silent
failure modes like:

  - Heavy atoms in a SMILES that the matcher doesn't have SMARTS for
    (silent functional-group truncation -- e.g. nitromethane decomposed
    as if it were CH3).
  - AIOMFAC subgroups 281-285 (aerosysprop's "Mark Barley" extensions)
    that the Fortran reference silently drops.
  - Subgroup IDs outside the Fortran AIOMFAC range (1..265).
  - Fractional or negative subgroup counts produced by matcher edge
    cases.
  - Missing AIOMFAC MR organic-ion parameter for a (main group, ion)
    pair present in the mixture.
  - Mass / mole fractions that don't sum to 1 or contain out-of-range
    values.
  - Electrolyte mixtures without water present.
  - Salts whose mass fraction exceeds the Fortran-tested upper bound.
  - Bisulfate / Qcca cases that the current port doesn't yet support.

Public API
----------
- :class:`Issue` -- a single diagnostic finding
- :class:`DiagnosticReport` -- container; ``is_ok``, ``has_errors``,
  ``errors``, ``warnings``, ``summary()``, iterable
- :func:`diagnose` -- dispatch a SMILES (and optionally a composition)
  through a method's pre-flight checks; returns a report without
  running the method
- :func:`diagnose_smiles` -- generic per-molecule check
  (validity + atom coverage against a given SMARTS set)
- :func:`diagnose_aiomfac_sr` -- SR-only AIOMFAC pre-flight
- :func:`diagnose_aiomfac` -- full SR+MR+LR AIOMFAC pre-flight
- :data:`STRICT_VALUES` -- the three allowed strict-mode values

Example
-------
>>> from aerosysprop.diagnostics import diagnose_aiomfac
>>> report = diagnose_aiomfac(
...     smiles=['O', 'Cc1ccccc1'],         # water + toluene
...     salts=['NH42SO4'],
...     mass_fractions=[0.85, 0.10, 0.05],
...     T_K=298.15,
... )
>>> report.is_ok
False
>>> for issue in report.errors:
...     print(issue.code, issue.message)
MR_PARAM_MISSING main_group 4 (ACCHn) has no parameter against anion 261 (SO4--) ...
"""

from dataclasses import dataclass, field
from typing import Optional, Any, Callable, Sequence, Mapping
import functools
import warnings

from rdkit import Chem

from aerosysprop import data, groups


# ---------------------------------------------------------------------------
# Issue codes (standardised identifiers used in DiagnosticReport)
# ---------------------------------------------------------------------------

# SMILES / matcher
INVALID_SMILES        = 'INVALID_SMILES'
EMPTY_DECOMPOSITION   = 'EMPTY_DECOMPOSITION'
UNMATCHED_ATOMS       = 'UNMATCHED_ATOMS'
EXTENSION_SUBGROUP    = 'EXTENSION_SUBGROUP'      # AIOMFAC 281..285
SUBGROUP_OUT_OF_RANGE = 'SUBGROUP_OUT_OF_RANGE'
FRACTIONAL_COUNT      = 'FRACTIONAL_COUNT'
NEGATIVE_COUNT        = 'NEGATIVE_COUNT'

# AIOMFAC composition / parameters
MR_PARAM_MISSING      = 'MR_PARAM_MISSING'
COMPOSITION_SUM       = 'COMPOSITION_SUM'
COMPOSITION_VALUE     = 'COMPOSITION_VALUE'
NO_WATER              = 'NO_WATER'
HIGH_CONCENTRATION    = 'HIGH_CONCENTRATION'
BISULFATE_NEEDED      = 'BISULFATE_NEEDED'
QCCA_NEEDED           = 'QCCA_NEEDED'
UNKNOWN_SALT          = 'UNKNOWN_SALT'

# Strict-mode arg values
STRICT_VALUES = ('silent', 'warn', 'raise')


# ---------------------------------------------------------------------------
# Core types
# ---------------------------------------------------------------------------

@dataclass
class Issue:
    """A single diagnostic finding."""
    severity: str       # 'error' | 'warning' | 'info'
    code: str           # one of the constants above
    message: str        # human-readable
    details: dict = field(default_factory=dict)

    def __str__(self):
        return f'[{self.severity.upper():7s}] {self.code:22s} {self.message}'


@dataclass
class DiagnosticReport:
    """Aggregate of issues found while screening a SMILES or composition."""
    issues: list = field(default_factory=list)
    context: dict = field(default_factory=dict)

    def add(self, severity, code, message, **details):
        self.issues.append(Issue(severity, code, message, dict(details)))

    @property
    def errors(self):
        return [i for i in self.issues if i.severity == 'error']

    @property
    def warnings(self):
        return [i for i in self.issues if i.severity == 'warning']

    @property
    def has_errors(self):
        return any(i.severity == 'error' for i in self.issues)

    @property
    def has_warnings(self):
        return any(i.severity == 'warning' for i in self.issues)

    @property
    def is_ok(self):
        return not self.issues

    def summary(self):
        if not self.issues:
            return 'No issues found.'
        return '\n'.join(str(i) for i in self.issues)

    def __iter__(self):
        return iter(self.issues)

    def __len__(self):
        return len(self.issues)


# ---------------------------------------------------------------------------
# Method registry: method -> list of SMARTS-collection variable names
# ---------------------------------------------------------------------------

def _smarts_for_aiomfac():
    return [data.AIOMFAC_UNIFAC_SMARTS, data.AIOMFAC_HCTAIL_SMARTS]


def _smarts_for_nannoolal():
    return [data.NANNOOLAL_SMARTS_PRIMARY, data.NANNOOLAL_SMARTS_SECONDARY]


def _smarts_for_stein_brown():
    return [data.STEIN_AND_BROWN_SMARTS]


def _smarts_for_simpol():
    return [data.SIMPOL_SMARTS]


def _smarts_for_evaporation():
    return [data.EVAPORATION_SMARTS]


def _smarts_for_evaporation2():
    return [data.EVAPORATION2_SMARTS]


def _smarts_for_girolami():
    return [data.GIROLAMI_SMARTS]


def _smarts_for_le_bas():
    return [data.LE_BAS_SMARTS]


def _smarts_for_schroeder():
    return [data.SCHROEDER_SMARTS]


# Registry built lazily so method modules don't have to import on package
# load. Keys are dotted-path strings; values are no-arg callables returning
# a list of SMARTS dicts to apply for coverage.
_METHOD_REGISTRY = {
    'aerosysprop.methods.boiling_points.nannoolal':       _smarts_for_nannoolal,
    'aerosysprop.methods.boiling_points.stein_and_brown': _smarts_for_stein_brown,
    'aerosysprop.methods.boiling_points.joback_and_reid': _smarts_for_stein_brown,
    'aerosysprop.methods.vapour_pressures.nannoolal':           _smarts_for_nannoolal,
    'aerosysprop.methods.vapour_pressures.myrdal_and_yalkowsky': _smarts_for_nannoolal,
    'aerosysprop.methods.vapour_pressures.simpol':              _smarts_for_simpol,
    'aerosysprop.methods.vapour_pressures.evaporation':         _smarts_for_evaporation,
    'aerosysprop.methods.vapour_pressures.evaporation2':        _smarts_for_evaporation2,
    'aerosysprop.methods.critical_properties.nannoolal':       _smarts_for_nannoolal,
    'aerosysprop.methods.critical_properties.joback_and_reid': _smarts_for_stein_brown,
    'aerosysprop.methods.liquid_densities.girolami':   _smarts_for_girolami,
    'aerosysprop.methods.liquid_densities.le_bas':     _smarts_for_le_bas,
    'aerosysprop.methods.liquid_densities.schroeder':  _smarts_for_schroeder,
    # AIOMFAC entries register the matcher set; the AIOMFAC-specific
    # composition checks live in diagnose_aiomfac() / diagnose_aiomfac_sr().
    'aerosysprop.methods.activity_coefficients.aiomfac_sr_from_smiles': _smarts_for_aiomfac,
    'aerosysprop.methods.activity_coefficients.aiomfac_from_smiles':    _smarts_for_aiomfac,
}


def _method_key(method):
    """Resolve a callable to its registry key."""
    if isinstance(method, str):
        return method
    return f'{method.__module__}.{method.__name__}'


# ---------------------------------------------------------------------------
# Per-molecule checks
# ---------------------------------------------------------------------------

def _heavy_atom_indices(mol):
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]


def _atoms_covered_by(mol, smarts_dicts):
    """Run every SMARTS pattern in `smarts_dicts` over the molecule and
    return the set of atom indices that matched at least one pattern."""
    covered = set()
    for smarts_dict in smarts_dicts:
        for sg_id, smarts_obj in smarts_dict.items():
            for match in smarts_obj.findall(mol):
                covered.update(match)
    return covered


def diagnose_smiles(smiles, smarts_dicts=None, report=None):
    """Screen a single SMILES string against an optional list of SMARTS
    dicts (e.g., the ones used by a particular predictive method).

    Returns a :class:`DiagnosticReport`. If ``smarts_dicts`` is None, only
    SMILES-validity checks are performed.
    """
    if report is None:
        report = DiagnosticReport()
    report.context['smiles'] = smiles

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        report.add('error', INVALID_SMILES,
                   f'RDKit could not parse SMILES {smiles!r}',
                   smiles=smiles)
        return report

    if smarts_dicts is not None:
        heavy = set(_heavy_atom_indices(mol))
        covered = _atoms_covered_by(mol, smarts_dicts)
        unmatched = heavy - covered
        if unmatched:
            unmatched_atoms = [
                (i, mol.GetAtomWithIdx(i).GetSymbol())
                for i in sorted(unmatched)
            ]
            report.add(
                'warning', UNMATCHED_ATOMS,
                f'{len(unmatched)} heavy atom(s) in {smiles!r} are not matched '
                'by any SMARTS pattern; the method will silently behave as '
                'if these atoms (and their attached functional groups) were '
                'absent.',
                unmatched=unmatched_atoms,
                smiles=smiles,
            )
    return report


@functools.lru_cache(maxsize=4096)
def _cached_aiomfac_smiles_analysis(smiles):
    """Memoised per-SMILES analysis: returns ``(decomp, frozen_issues)`` where
    ``decomp`` is the subgroup dict (or None if invalid) and
    ``frozen_issues`` is a tuple of (severity, code, message, frozen_details)
    quadruples that callers can replay into a DiagnosticReport.

    Repeated calls with the same SMILES are essentially free — only the
    composition-dependent portion of `diagnose_aiomfac` needs to re-run on
    every method call. Critical for hot loops with thousands of repeated
    AIOMFAC evaluations.
    """
    issues_out = []
    def emit(severity, code, message, **details):
        issues_out.append((severity, code, message,
                           tuple(sorted(details.items()))))

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        emit('error', INVALID_SMILES,
             f'RDKit could not parse SMILES {smiles!r}', smiles=smiles)
        return None, tuple(issues_out)
    decomp = {k: v for k, v in groups.aiomfac(mol).items() if v != 0}
    if not decomp:
        emit('error', EMPTY_DECOMPOSITION,
             f'AIOMFAC matcher returned no subgroups for {smiles!r}',
             smiles=smiles)
        return decomp, tuple(issues_out)

    ext = [k for k in decomp if 281 <= k <= 285]
    if ext:
        emit('warning', EXTENSION_SUBGROUP,
             f'{smiles!r} decomposes to extension subgroup(s) {ext}: aerosysprop '
             'inherits these from a 2012 "Mark Barley" extension and they are not '
             'part of the Fortran AIOMFAC parameter set (topsubno=265). Fortran '
             'AIOMFAC would silently drop them; gamma values for this compound '
             'may be unreliable.',
             smiles=smiles, subgroups=tuple(ext))

    oor = [k for k in decomp if k < 1 or k > 265]
    if oor:
        emit('error', SUBGROUP_OUT_OF_RANGE,
             f'{smiles!r} produces subgroup IDs outside Fortran AIOMFAC range '
             f'(1..265): {oor}',
             smiles=smiles, subgroups=tuple(oor))

    neg = [(k, v) for k, v in decomp.items() if v < 0]
    if neg:
        emit('error', NEGATIVE_COUNT,
             f'{smiles!r} has negative subgroup count(s) {neg} — likely a '
             'matcher bug.',
             smiles=smiles, subgroups=tuple(neg))

    frac = [(k, v) for k, v in decomp.items() if abs(v - round(v)) > 1e-9]
    if frac:
        emit('warning', FRACTIONAL_COUNT,
             f'{smiles!r} has fractional subgroup count(s) {frac}; can arise '
             'from ring-shared SMARTS matches and may signal an incomplete '
             'decomposition.',
             smiles=smiles, subgroups=tuple(frac))

    # Atom-coverage check (silent matcher truncation)
    heavy = set(_heavy_atom_indices(mol))
    covered = _atoms_covered_by(mol, _smarts_for_aiomfac())
    unmatched = heavy - covered
    if unmatched:
        unmatched_atoms = tuple(
            (i, mol.GetAtomWithIdx(i).GetSymbol())
            for i in sorted(unmatched)
        )
        emit('warning', UNMATCHED_ATOMS,
             f'{len(unmatched)} heavy atom(s) in {smiles!r} are not matched '
             'by any SMARTS pattern; the method will silently behave as if '
             'these atoms (and their attached functional groups) were absent.',
             unmatched=unmatched_atoms, smiles=smiles)

    return decomp, tuple(issues_out)


def _check_aiomfac_decomposition(smiles, report):
    """AIOMFAC-specific subgroup-level checks. Uses
    :func:`_cached_aiomfac_smiles_analysis` so repeated calls with the same
    SMILES amortise to near-zero overhead. Returns the decomposition dict
    (or None if the SMILES is invalid)."""
    decomp, cached_issues = _cached_aiomfac_smiles_analysis(smiles)
    for severity, code, message, frozen_details in cached_issues:
        report.add(severity, code, message, **dict(frozen_details))
    return decomp


def clear_diagnostic_cache():
    """Drop the memoised per-SMILES analysis cache. Useful if the AIOMFAC
    matcher SMARTS data is reloaded at runtime (rare)."""
    _cached_aiomfac_smiles_analysis.cache_clear()


def diagnostic_cache_info():
    """Return the lru_cache info tuple for the per-SMILES analysis cache."""
    return _cached_aiomfac_smiles_analysis.cache_info()


def diagnose_aiomfac_sr(smiles, mole_fractions=None, T_K=None):
    """Pre-flight check for `aiomfac_sr_from_smiles`. Returns a report
    covering SMILES validity, matcher coverage, AIOMFAC subgroup
    well-formedness, and composition-vector sanity.

    Pass ``mole_fractions`` and ``T_K`` to also check composition values.
    """
    report = DiagnosticReport()
    report.context['method'] = 'aiomfac_sr_from_smiles'
    for smi in smiles:
        _check_aiomfac_decomposition(smi, report)

    if mole_fractions is not None:
        if len(mole_fractions) != len(smiles):
            report.add(
                'error', COMPOSITION_VALUE,
                f'mole_fractions length {len(mole_fractions)} does not match '
                f'number of SMILES ({len(smiles)})',
            )
        else:
            if any(x < 0.0 or x > 1.0 for x in mole_fractions):
                report.add(
                    'error', COMPOSITION_VALUE,
                    f'mole_fractions has values outside [0, 1]: {mole_fractions}',
                )
            total = sum(mole_fractions)
            if abs(total - 1.0) > 1e-2:
                report.add(
                    'warning', COMPOSITION_SUM,
                    f'mole_fractions sum to {total:.6f}, not 1.0 (off by '
                    f'{abs(total-1):.2e}); UNIFAC operates on absolute mole '
                    'numbers so the result is well-defined but you may have '
                    'intended a normalised composition.',
                    total=total,
                )
    if T_K is not None and (T_K < 100.0 or T_K > 1000.0):
        report.add(
            'warning', COMPOSITION_VALUE,
            f'T_K = {T_K} is outside the physically reasonable range '
            '[100, 1000] K',
        )
    return report


def _lookup_salt_record(name):
    matches = [
        sid for sid, sname in data.AIOMFAC_SALT_NAME.items() if sname == name
    ]
    if not matches:
        return None
    sid = matches[0]
    return (
        sid,
        int(data.AIOMFAC_SALT_CATION_GROUP[sid]),
        int(data.AIOMFAC_SALT_CATION_STOICH[sid]),
        int(data.AIOMFAC_SALT_ANION_GROUP[sid]),
        int(data.AIOMFAC_SALT_ANION_STOICH[sid]),
    )


def diagnose_aiomfac(smiles, salts, mass_fractions, T_K=None):
    """Pre-flight check for `aiomfac_from_smiles` (full SR+MR+LR).

    Catches every issue the diagnostic framework knows about: SMILES
    validity and matcher coverage on each neutral, AIOMFAC subgroup
    well-formedness, unknown salt names, mass-fraction sums and value
    bounds, missing organic-ion MR parameters for any (main group, ion)
    pair present in the mixture, charge balance, water presence,
    bisulfate / Qcca caveats.
    """
    report = DiagnosticReport()
    report.context['method'] = 'aiomfac_from_smiles'

    # ---- Neutral SMILES checks ----
    neutral_decomps = []
    for smi in smiles:
        decomp = _check_aiomfac_decomposition(smi, report)
        neutral_decomps.append(decomp)

    # ---- Salts ----
    salt_records = []
    for name in salts:
        rec = _lookup_salt_record(name)
        if rec is None:
            report.add(
                'error', UNKNOWN_SALT,
                f'salt {name!r} not in AIOMFAC catalog; known names include '
                + ', '.join(sorted(set(data.AIOMFAC_SALT_NAME.values()))[:10])
                + ', ...',
                salt=name,
            )
        else:
            salt_records.append((name, *rec[1:]))

    # ---- Composition sanity ----
    n_expected = len(smiles) + len(salts)
    if len(mass_fractions) != n_expected:
        report.add(
            'error', COMPOSITION_VALUE,
            f'mass_fractions length {len(mass_fractions)} does not match '
            f'n_neutral ({len(smiles)}) + n_salts ({len(salts)})',
        )
        return report
    if any(w < 0.0 or w > 1.0 for w in mass_fractions):
        report.add(
            'error', COMPOSITION_VALUE,
            f'mass_fractions has values outside [0, 1]: {mass_fractions}',
        )
    total = sum(mass_fractions)
    if abs(total - 1.0) > 1e-3:
        report.add(
            'warning', COMPOSITION_SUM,
            f'mass_fractions sum to {total:.6f}, not 1.0; AIOMFAC requires '
            'normalised compositions.',
            total=total,
        )

    # Water presence (Fortran emits a warning for electrolyte mixes without
    # water as the canonical first solvent).
    if salts and not any(s in ('O', '[H]O[H]') for s in smiles):
        report.add(
            'warning', NO_WATER,
            'Electrolyte mixture without water (SMILES "O") as the first '
            'neutral; Fortran AIOMFAC assumes water is the reference solvent '
            'in the LR Debye-Hueckel and the ion molality convention.',
        )

    # ---- High concentration check (per-salt) ----
    # Fortran ships per-(NC,NA) maximum-validated mass fractions in
    # TABhighestwtf. We don't load it currently; warn for any salt > 30%
    # mass fraction as a coarse heuristic.
    for k, name in enumerate(salts):
        idx = len(smiles) + k
        if idx < len(mass_fractions) and mass_fractions[idx] > 0.30:
            report.add(
                'warning', HIGH_CONCENTRATION,
                f'{name} at mass fraction {mass_fractions[idx]:.3f} is above '
                '~30% — the AIOMFAC parameterisation may not have been fitted '
                'to this concentration range; gamma values become less reliable '
                'near saturation.',
                salt=name, mass_fraction=mass_fractions[idx],
            )

    # ---- Bisulfate / Qcca caveats ----
    salt_names_set = set(salts)
    has_H2SO4 = 'H2SO4' in salt_names_set
    has_NH4 = any(n.startswith('NH4') or n.startswith('NH42') for n in salt_names_set)
    has_H = any(n.startswith('H') and n not in ('Hg', 'He') for n in salt_names_set)
    if has_H2SO4:
        report.add(
            'warning', BISULFATE_NEEDED,
            'H2SO4 in mixture — Fortran AIOMFAC solves the bisulfate '
            'equilibrium HSO4- <-> H+ + SO4-- iteratively; the Python port '
            'treats H2SO4 as fully dissociated to 2 H+ + SO4--, so gammas '
            'will not match Fortran.',
        )
    if has_NH4 and has_H:
        report.add(
            'warning', QCCA_NEEDED,
            'NH4+ and H+ both present — Fortran AIOMFAC applies extra MR '
            'triple-ion (Qcca) and cation-cation (Rcc) terms which are not '
            'implemented in the Python port; gammas will not match Fortran '
            'for this combination.',
        )

    # ---- MR parameter availability ----
    if salt_records and all(d for d in neutral_decomps):
        SENTINEL = -1000.0
        # main groups present in the neutrals
        present_mgs = set()
        for decomp in neutral_decomps:
            for sg in decomp:
                mg = data.AIOMFAC_MAIN_GROUP.get(sg)
                if mg is not None:
                    present_mgs.add(int(mg))
        # ion IDs in the mix
        cation_ids = sorted({r[1] for r in salt_records})
        anion_ids  = sorted({r[3] for r in salt_records})
        # check each (main_group, ion) pair
        for mg in sorted(present_mgs):
            row_b = data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b1.get(mg, {})
            row_c = data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b2.get(mg, {})
            for ion_id in cation_ids + anion_ids:
                b = row_b.get(ion_id, 0.0)
                c = row_c.get(ion_id, 0.0)
                if b < SENTINEL or c < SENTINEL:
                    kind = 'cation' if ion_id < 240 else 'anion'
                    report.add(
                        'error', MR_PARAM_MISSING,
                        f'main_group {mg} has no MR parameter against {kind} '
                        f'{ion_id} (Fortran sentinel -7.7777e4). The Fortran '
                        'AIOMFAC reference would reject this mixture (error 1 '
                        'or 2); the Python port also raises ValueError when '
                        'this combination is used.',
                        main_group=mg, ion=ion_id, kind=kind,
                    )

    return report


# ---------------------------------------------------------------------------
# Top-level dispatcher
# ---------------------------------------------------------------------------

def diagnose(method, *args, **kwargs):
    """Dispatch a SMILES (or mixture) through a method's pre-flight checks.

    Parameters
    ----------
    method : callable or str
        An aerosysprop method (e.g., ``boiling_points.nannoolal``) or its
        fully-qualified dotted path (``'aerosysprop.methods.boiling_points.nannoolal'``).
    *args, **kwargs
        For AIOMFAC entry points, pass the same arguments you would pass to
        the method itself (smiles list, salts, mass_fractions, T_K). For
        property methods (boiling points, vapour pressures, etc.), the
        first positional argument is the molecule (SMILES str or RDKit
        Mol); other arguments are ignored by the diagnostic.

    Returns
    -------
    DiagnosticReport
    """
    key = _method_key(method)
    name = key.rsplit('.', 1)[1] if '.' in key else key

    # AIOMFAC dispatch
    if name == 'aiomfac_from_smiles':
        smiles    = kwargs.get('smiles', args[0] if len(args) > 0 else None)
        salts     = kwargs.get('salts',  args[1] if len(args) > 1 else [])
        mass_frac = kwargs.get('mass_fractions', args[2] if len(args) > 2 else None)
        T_K       = kwargs.get('T_K',    args[3] if len(args) > 3 else None)
        return diagnose_aiomfac(smiles, salts, mass_frac, T_K)
    if name == 'aiomfac_sr_from_smiles':
        smiles = kwargs.get('smiles', args[0] if len(args) > 0 else None)
        mole_f = kwargs.get('mole_fractions', args[1] if len(args) > 1 else None)
        T_K    = kwargs.get('T_K',   args[2] if len(args) > 2 else None)
        return diagnose_aiomfac_sr(smiles, mole_f, T_K)

    # Property-method dispatch: first arg is a SMILES or RDKit Mol
    smarts_provider = _METHOD_REGISTRY.get(key)
    if smarts_provider is None:
        report = DiagnosticReport()
        report.add('warning', 'NO_DIAGNOSTIC',
                   f'No registered diagnostic for {key}; only generic '
                   'SMILES validity is checked.')
    mol_arg = args[0] if args else kwargs.get('mol') or kwargs.get('compound')
    if isinstance(mol_arg, str):
        smi = mol_arg
    elif mol_arg is not None and hasattr(mol_arg, 'GetAtoms'):
        smi = Chem.MolToSmiles(mol_arg)
    else:
        report = DiagnosticReport()
        report.add('error', INVALID_SMILES,
                   'Could not interpret first argument as SMILES or RDKit Mol')
        return report
    smarts_dicts = smarts_provider() if smarts_provider else None
    return diagnose_smiles(smi, smarts_dicts=smarts_dicts)


# ---------------------------------------------------------------------------
# strict= integration helper for method entry points
# ---------------------------------------------------------------------------

def with_diagnostics(func):
    """Decorator: add a ``strict='warn'|'raise'|'silent'`` argument to a
    property method that runs the registered pre-flight diagnostic over
    its first positional argument (assumed to be a SMILES string or an
    RDKit Mol object). Use on ``methods.*`` functions; AIOMFAC entry
    points integrate diagnostics inline instead.

    Example::

        @with_diagnostics
        def nannoolal(compound, temperature, boiling_point):
            ...
    """
    import functools

    @functools.wraps(func)
    def wrapper(compound, *args, strict='warn', **kwargs):
        if strict != 'silent':
            report = diagnose(func, compound)
            apply_strict(report, strict)
        return func(compound, *args, **kwargs)
    return wrapper


def apply_strict(report, strict):
    """Apply a strict-mode argument to a DiagnosticReport.

    Used by method entry points to translate the strict argument into
    warnings or exceptions.
    """
    if strict not in STRICT_VALUES:
        raise ValueError(
            f'strict must be one of {STRICT_VALUES}, got {strict!r}'
        )
    if strict == 'silent':
        return
    if strict == 'raise':
        if report.has_errors or report.has_warnings:
            raise ValueError(
                'aerosysprop diagnostics flagged issues:\n' + report.summary()
            )
        return
    if strict == 'warn':
        for issue in report.issues:
            warnings.warn(str(issue), UserWarning, stacklevel=3)
