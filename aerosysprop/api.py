# This file is part of AeroSysProp.
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
"""
aerosysprop.api
===============
Agent-facing API: every function takes plain values (strings, numbers, lists)
and returns plain JSON-compatible dicts with explicit units, so the same layer
can back the MCP server, a REST endpoint or ordinary scripts.

Design rules
------------
* SMILES strings in, never RDKit ``Mol`` objects.
* Outputs are labelled dicts, never positional tuples. NaN/inf become ``None``.
* Every numeric result is returned alongside the diagnostics that apply to it,
  so a caller cannot see a value without also seeing its caveats.
* Anticipated problems raise :class:`ApiError`, which carries a stable
  ``code``, a message and a ``hint`` telling the caller how to fix the call.
* Nothing in this module writes to disk. Session databases are opened
  read-only, and file access is confined to the ``AEROSYSPROP_DATA_DIR``
  directory.
"""
from __future__ import annotations

import functools
import json
import math
import os
import re
import sqlite3
import urllib.error
import urllib.request
import warnings
from pathlib import Path
from typing import Any, Optional

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

from aerosysprop import __version__, batch, data, db
from aerosysprop import io as asp_io
from aerosysprop.diagnostics import (
    DiagnosticReport,
    diagnose,
    diagnose_aiomfac,
    diagnose_aiomfac_sr,
)
from aerosysprop.mixture_methods import MIXTURE_METHODS

__all__ = [
    'ApiError', 'list_methods', 'list_salts', 'predict_properties',
    'diagnose_smiles', 'activity_coefficients', 'list_data_files',
    'load_smiles_file', 'read_session_database', 'opera_predict',
    'DATA_DIR_ENV', 'MAX_MOLECULES',
]

DATA_DIR_ENV = 'AEROSYSPROP_DATA_DIR'
MAX_MOLECULES = int(os.environ.get('AEROSYSPROP_MAX_MOLECULES', '1000'))
MAX_OPERA_MOLECULES = 100
SMILES_FILE_SUFFIXES = {'.txt', '.prop', '.smi', '.csv', '.tsv'}
DB_SUFFIXES = {'.db', '.sqlite', '.sqlite3'}
SDF_SUFFIXES = {'.sdf'}


# ---------------------------------------------------------------------------
# Errors and JSON helpers
# ---------------------------------------------------------------------------

class ApiError(ValueError):
    """An anticipated problem with a request, with a hint on how to fix it."""

    def __init__(self, code: str, message: str, hint: str = '', details=None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.hint = hint
        self.details = details

    def to_dict(self) -> dict:
        out = {'code': self.code, 'message': self.message}
        if self.hint:
            out['hint'] = self.hint
        if self.details is not None:
            out['details'] = json_safe(self.details)
        return out

    def as_text(self) -> str:
        return json.dumps({'error': self.to_dict()}, ensure_ascii=False)


def json_safe(value: Any) -> Any:
    """Recursively convert a value into something ``json.dumps`` accepts,
    turning NaN/inf into ``None``, numpy scalars/arrays into Python values,
    complex numbers into their real part and sets/tuples into lists."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (complex, np.complexfloating)):
        value = value.real
    if isinstance(value, (float, np.floating)):
        f = float(value)
        return f if math.isfinite(f) else None
    if isinstance(value, np.ndarray):
        return [json_safe(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return [json_safe(v) for v in sorted(value, key=str)]
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return str(value)


def _issues(report: Optional[DiagnosticReport]) -> list:
    if report is None:
        return []
    out = []
    for i in report.issues:
        d = {'severity': i.severity, 'code': i.code, 'message': i.message}
        if i.details:
            d['details'] = json_safe(i.details)
        out.append(d)
    return out


def _parse_smiles(smiles: str):
    if not smiles:
        return None
    with rdBase.BlockLogs():
        return Chem.MolFromSmiles(smiles)


def _as_str_list(value, argname: str) -> list:
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, (list, tuple)) or not value:
        raise ApiError('INVALID_ARGUMENT', f'{argname} must be a non-empty list of strings.')
    out = []
    for v in value:
        if not isinstance(v, str):
            raise ApiError('INVALID_ARGUMENT',
                           f'{argname} must contain only strings; got {v!r}.')
        out.append(v.strip())
    return out


def _as_float_list(value, argname: str, expected_len: int) -> list:
    if not isinstance(value, (list, tuple)):
        raise ApiError('INVALID_ARGUMENT', f'{argname} must be a list of numbers.')
    try:
        out = [float(v) for v in value]
    except (TypeError, ValueError):
        raise ApiError('INVALID_ARGUMENT', f'{argname} must contain only numbers.')
    if len(out) != expected_len:
        raise ApiError(
            'INVALID_ARGUMENT',
            f'{argname} has {len(out)} values but {expected_len} were expected.',
            hint='Give one fraction per neutral component, followed by one per salt, '
                 'in the same order as the smiles and salts lists.')
    return out


def _check_count(n: int, limit: int = MAX_MOLECULES):
    if n > limit:
        raise ApiError(
            'TOO_MANY_MOLECULES',
            f'{n} molecules requested; the limit per call is {limit}.',
            hint='Split the list into smaller batches.')


def _check_temperature(T) -> tuple:
    try:
        T = float(T)
    except (TypeError, ValueError):
        raise ApiError('INVALID_ARGUMENT', 'temperature_K must be a number in kelvin.')
    if not math.isfinite(T) or T <= 0:
        raise ApiError('INVALID_ARGUMENT',
                       f'temperature_K must be a positive number of kelvin; got {T}.',
                       hint='Convert Celsius to kelvin by adding 273.15.')
    notes = []
    if T < 100.0 or T > 1000.0:
        notes.append({'severity': 'warning', 'code': 'TEMPERATURE_RANGE',
                      'message': f'temperature_K = {T} is outside the physically '
                                 'reasonable range [100, 1000] K.'})
    return T, notes


# ---------------------------------------------------------------------------
# Method catalogue
# ---------------------------------------------------------------------------

# Method ID -> dotted path of the underlying function in the diagnostics
# registry. Methods absent here have no matcher-coverage diagnostic.
_DIAG_PATH = {
    'bp_joback':         'aerosysprop.methods.boiling_points.joback_and_reid',
    'bp_stein':          'aerosysprop.methods.boiling_points.stein_and_brown',
    'bp_nannoolal':      'aerosysprop.methods.boiling_points.nannoolal',
    'vp_nannoolal':      'aerosysprop.methods.vapour_pressures.nannoolal',
    'vp_myrdal':         'aerosysprop.methods.vapour_pressures.myrdal_and_yalkowsky',
    'vp_evaporation':    'aerosysprop.methods.vapour_pressures.evaporation',
    'vp_evaporation2':   'aerosysprop.methods.vapour_pressures.evaporation2',
    'vp_simpol':         'aerosysprop.methods.vapour_pressures.simpol',
    'cp_nannoolal':      'aerosysprop.methods.critical_properties.nannoolal',
    'density_girolami':  'aerosysprop.methods.liquid_densities.girolami',
    'density_schroeder': 'aerosysprop.methods.liquid_densities.schroeder',
    'density_le_bas':    'aerosysprop.methods.liquid_densities.le_bas',
}


def _slug(text: str) -> str:
    return re.sub(r'[^a-z]+', '_', text.lower()).strip('_')


_GROUP_ALIASES = {_slug(m.group): m.group for m in batch.PROPERTY_METHODS}
_GROUP_ALIASES.update({
    'vapor_pressure': 'Vapour Pressure', 'vapour_pressures': 'Vapour Pressure',
    'vapor_pressures': 'Vapour Pressure', 'boiling_points': 'Boiling Point',
    'critical_property': 'Critical Properties', 'density': 'Liquid Density',
    'liquid_densities': 'Liquid Density', 'densities': 'Liquid Density',
})

_UNIT_SUFFIX = {'_K': 'K', '_bar': 'bar', '_log10atm': 'log10(atm)', '_gcc': 'g/cm3'}


def _unit_for(column: str) -> str:
    base = column.split('__bp', 1)[0]
    for suffix, unit in _UNIT_SUFFIX.items():
        if base.endswith(suffix):
            return unit
    return ''


def _normalise_methods(methods) -> Optional[list]:
    if methods is None or methods == 'all' or (isinstance(methods, (list, tuple)) and not methods):
        return None
    out = []
    for m in _as_str_list(methods, 'methods'):
        if m in batch.PROPERTY_METHOD_BY_ID:
            if m not in out:
                out.append(m)
            continue
        if _slug(m) in _GROUP_ALIASES:
            group = _GROUP_ALIASES[_slug(m)]
            out.extend(i for i in batch.METHOD_GROUPS[group] if i not in out)
            continue
        raise ApiError(
            'UNKNOWN_METHOD', f'Unknown method {m!r}.',
            hint='Call list_methods for valid IDs. Property names such as '
                 '"vapour_pressure" or "boiling_point" select every method in that group.',
            details={'valid_method_ids': batch.ALL_METHODS,
                     'valid_groups': sorted(set(_GROUP_ALIASES))})
    return out


def _normalise_bp_sources(bp_sources) -> list:
    if bp_sources is None or (isinstance(bp_sources, (list, tuple)) and not bp_sources):
        return ['bp_nannoolal']
    out = []
    for b in _as_str_list(bp_sources, 'bp_sources'):
        key = b.lower()
        key = key if key.startswith('bp_') else 'bp_' + key
        key = {'bp_stein_and_brown': 'bp_stein',
               'bp_joback_and_reid': 'bp_joback'}.get(key, key)
        if key not in batch.BP_METHODS:
            raise ApiError('UNKNOWN_METHOD', f'Unknown boiling point source {b!r}.',
                           details={'valid_bp_sources': sorted(batch.BP_METHODS)})
        if key not in out:
            out.append(key)
    return out


def list_methods() -> dict:
    """Describe every available property method, its unit and dependencies."""
    single = []
    for m in batch.PROPERTY_METHODS:
        single.append({
            'id': m.id,
            'label': m.label,
            'property': m.group,
            'output_columns': {c: _unit_for(c) for c in m.output_keys},
            'depends_on': sorted('a boiling point source' if d == batch.BP_ANY else d
                                 for d in m.deps),
            'uses_boiling_point_source': m.bp_dependent,
            'has_coverage_diagnostic': m.id in _DIAG_PATH,
        })
    mixture = [{
        'id': m.id, 'name': m.name, 'status': m.status, 'summary': m.summary,
        'when_to_use': m.when_to_use, 'diagnostic_notes': m.diagnostics_notes,
    } for m in MIXTURE_METHODS]
    return {
        'aerosysprop_version': __version__,
        'single_component_methods': single,
        'property_groups': list(batch.METHOD_GROUPS),
        'boiling_point_sources': sorted(batch.BP_METHODS),
        'default_boiling_point_source': 'bp_nannoolal',
        'mixture_methods': mixture,
        'notes': [
            'All values are group-contribution estimates, not measurements.',
            'Vapour pressures are log10 of the pure-component saturation '
            'vapour pressure in atmospheres at the requested temperature.',
            'Methods that need a boiling point use the chosen boiling point '
            'source; results can differ substantially between sources.',
            'Temperatures are in kelvin.',
        ],
    }


# ---------------------------------------------------------------------------
# Salts and ions (AIOMFAC)
# ---------------------------------------------------------------------------

def _salt_key(name: str) -> str:
    return re.sub(r'[\s()]', '', name).lower()


@functools.lru_cache(maxsize=1)
def _ion_names() -> dict:
    """Map AIOMFAC ion group IDs to readable names, taken from the comments
    in ``aiomfac_ions.data``."""
    names = {}
    path = Path(data.__file__).parent / 'aiomfac_ions.data'
    for line in path.read_text().splitlines():
        m = re.match(r'\s*(\d+)\s+[-+\d.]+\s*#\s*(\S+)', line)
        if m:
            names[int(m.group(1))] = m.group(2)
    return names


def _ion_label(ion_id) -> str:
    return _ion_names().get(int(ion_id), f'ion_{int(ion_id)}')


@functools.lru_cache(maxsize=1)
def _salt_lookup() -> dict:
    return {_salt_key(name): name for name in data.AIOMFAC_SALT_NAME.values()}


def _normalise_salts(salts) -> list:
    if salts is None or (isinstance(salts, (list, tuple)) and not salts):
        return []
    out = []
    for s in _as_str_list(salts, 'salts'):
        name = _salt_lookup().get(_salt_key(s))
        if name is None:
            raise ApiError(
                'UNKNOWN_SALT', f'Salt {s!r} is not in the AIOMFAC catalogue.',
                hint='Call list_salts for the accepted names.',
                details={'known_salts': sorted(set(data.AIOMFAC_SALT_NAME.values()))})
        out.append(name)
    return out


def list_salts() -> dict:
    """List the salts AIOMFAC supports, with their constituent ions."""
    salts = []
    for sid, name in sorted(data.AIOMFAC_SALT_NAME.items()):
        also = []
        plain = re.sub(r'[\s()]', '', name)
        if plain != name:
            also.append(plain)
        if name == 'NH42SO4':
            also.append('(NH4)2SO4')
        salts.append({
            'name': name,
            'also_accepted': also,
            'cation': _ion_label(data.AIOMFAC_SALT_CATION_GROUP[sid]),
            'cation_stoichiometry': json_safe(data.AIOMFAC_SALT_CATION_STOICH[sid]),
            'anion': _ion_label(data.AIOMFAC_SALT_ANION_GROUP[sid]),
            'anion_stoichiometry': json_safe(data.AIOMFAC_SALT_ANION_STOICH[sid]),
            'molar_mass_g_per_mol': json_safe(data.AIOMFAC_SALT_MASS[sid]),
        })
    return {'salts': salts,
            'note': 'Salt names are matched ignoring brackets and case, so '
                    '"(NH4)2SO4" and "CaCl2" are accepted.'}


# ---------------------------------------------------------------------------
# Single-component properties
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=16384)
def _coverage_issues(method_id: str, smiles: str) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        report = diagnose(_DIAG_PATH[method_id], smiles)
    return json.dumps(_issues(report))


def _method_diagnostics(smiles: str, method_ids) -> dict:
    out = {}
    for mid in method_ids:
        if mid not in _DIAG_PATH:
            continue
        issues = json.loads(_coverage_issues(mid, smiles))
        if issues:
            out[mid] = issues
    return out


def _molecule_summary(mol) -> dict:
    return {
        'canonical_smiles': Chem.MolToSmiles(mol),
        'formula': rdMolDescriptors.CalcMolFormula(mol),
        'heavy_atoms': mol.GetNumHeavyAtoms(),
        'formal_charge': Chem.GetFormalCharge(mol),
    }


def predict_properties(smiles, names=None, methods=None, temperature_K=298.15,
                       bp_sources=None, include_diagnostics=True) -> dict:
    """Predict single-component properties for one or more molecules.

    Parameters mirror :func:`aerosysprop.batch.run_all`, but take SMILES
    strings, accept property-group names (e.g. ``'vapour_pressure'``) as well
    as method IDs, and never write output files.
    """
    smiles = _as_str_list(smiles, 'smiles')
    _check_count(len(smiles))
    if names is not None:
        names = _as_str_list(names, 'names')
        if len(names) != len(smiles):
            raise ApiError('INVALID_ARGUMENT', 'names must be the same length as smiles.')
    else:
        names = list(smiles)
    method_ids = _normalise_methods(methods)
    bp = _normalise_bp_sources(bp_sources)
    T, notes = _check_temperature(temperature_K)
    resolved = [m for m in batch.ALL_METHODS
                if m in batch.resolve_methods(method_ids or batch.ALL_METHODS, bp_sources=bp)]

    records, mols = [], []
    for i, (smi, name) in enumerate(zip(smiles, names)):
        rec = {'index': i, 'name': name, 'input_smiles': smi}
        mol = _parse_smiles(smi)
        if mol is None:
            rec['error'] = {'code': 'INVALID_SMILES',
                            'message': f'RDKit could not parse SMILES {smi!r}.'}
        else:
            rec.update(_molecule_summary(mol))
            # Unique internal keys: run_all indexes results by name.
            mols.append((f'__{i}', mol))
        records.append(rec)

    by_key = {}
    if mols:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rows = batch.run_all(
                mols, methods=method_ids, temperature=T,
                workers=min(4, len(mols)), use_threads=True,
                progress_callback=lambda *_: None, bp_sources=bp)
        by_key = {key: row for (key, _), row in zip(mols, rows)}

    columns = {}
    for rec in records:
        if 'error' in rec:
            continue
        row = by_key.get(f'__{rec["index"]}') or {}
        values = {k: json_safe(v) for k, v in row.items() if k not in ('name', 'smiles')}
        for k in values:
            columns.setdefault(k, _unit_for(k))
        rec['values'] = values
        failed = sorted(k for k, v in values.items() if v is None)
        if failed:
            rec['not_computed'] = failed
        if include_diagnostics:
            diag = _method_diagnostics(rec['canonical_smiles'], resolved)
            if diag:
                rec['diagnostics'] = diag

    return {
        'temperature_K': T,
        'boiling_point_sources': bp,
        'methods_run': resolved,
        'columns': columns,
        'results': records,
        'summary': {
            'molecules': len(records),
            'invalid_smiles': sum(1 for r in records if 'error' in r),
            'with_diagnostic_warnings': sum(1 for r in records if r.get('diagnostics')),
        },
        'warnings': notes,
        'interpretation': (
            'Values are group-contribution estimates. Where a molecule lists '
            'diagnostics for a method, that method does not recognise some of '
            'its atoms or groups, so treat that value with caution. '
            '"not_computed" lists outputs the method could not produce.'),
    }


def diagnose_smiles(smiles, methods=None, include_aiomfac=True) -> dict:
    """Pre-flight check of SMILES strings against each method's group set."""
    smiles = _as_str_list(smiles, 'smiles')
    _check_count(len(smiles))
    method_ids = _normalise_methods(methods) or list(batch.ALL_METHODS)
    checked = [m for m in method_ids if m in _DIAG_PATH]
    if include_aiomfac:
        checked.append('aiomfac')
    results = []
    for i, smi in enumerate(smiles):
        rec = {'index': i, 'input_smiles': smi}
        mol = _parse_smiles(smi)
        rec['valid'] = mol is not None
        if mol is None:
            rec['error'] = {'code': 'INVALID_SMILES',
                            'message': f'RDKit could not parse SMILES {smi!r}.'}
            results.append(rec)
            continue
        rec.update(_molecule_summary(mol))
        diag = _method_diagnostics(rec['canonical_smiles'], method_ids)
        if include_aiomfac:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                aio = _issues(diagnose_aiomfac_sr((rec['canonical_smiles'],)))
            if aio:
                diag['aiomfac'] = aio
        rec['issues_by_method'] = diag
        rec['methods_without_issues'] = [m for m in checked if m not in diag]
        results.append(rec)
    return {'results': results,
            'methods_checked': checked,
            'methods_without_coverage_check': [m for m in method_ids if m not in _DIAG_PATH]}


# ---------------------------------------------------------------------------
# Mixtures: AIOMFAC activity coefficients
# ---------------------------------------------------------------------------

def activity_coefficients(smiles, temperature_K=298.15, mole_fractions=None,
                          mass_fractions=None, salts=None, bisulfate=True) -> dict:
    """AIOMFAC activity coefficients for an organic/water(/salt) mixture.

    * Without salts, give ``mole_fractions`` (short-range AIOMFAC only) or
      ``mass_fractions`` (full model, no ions).
    * With salts, give ``mass_fractions`` for every neutral component
      followed by every salt.
    """
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_from_smiles,
        aiomfac_sr_from_smiles,
    )
    smiles = _as_str_list(smiles, 'smiles')
    salts = _normalise_salts(salts)
    T, notes = _check_temperature(temperature_K)
    if (mole_fractions is None) == (mass_fractions is None):
        raise ApiError('INVALID_ARGUMENT',
                       'Give exactly one of mole_fractions or mass_fractions.',
                       hint='Use mass_fractions when salts are present.')
    if salts and mole_fractions is not None:
        raise ApiError('INVALID_ARGUMENT',
                       'Mixtures containing salts must be specified with mass_fractions.',
                       hint='List mass fractions for each neutral component, then each salt.')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if mole_fractions is not None:
            x = _as_float_list(mole_fractions, 'mole_fractions', len(smiles))
            report = diagnose_aiomfac_sr(tuple(smiles), tuple(x), T)
            if report.has_errors:
                raise ApiError('DIAGNOSTIC_ERRORS', 'The mixture failed pre-flight checks.',
                               details=_issues(report))
            gamma = np.asarray(aiomfac_sr_from_smiles(smiles, x, T, strict='silent'))
            model = 'AIOMFAC short-range (modified UNIFAC)'
            neutral = [{'smiles': s, 'mole_fraction': xi,
                        'activity_coefficient': json_safe(g),
                        'activity': json_safe(g * xi)}
                       for s, xi, g in zip(smiles, x, gamma)]
            ions, salt_rows, basis = [], [], 'mole_fraction'
        else:
            w = _as_float_list(mass_fractions, 'mass_fractions', len(smiles) + len(salts))
            report = diagnose_aiomfac(smiles, salts, w, T)
            if report.has_errors:
                raise ApiError('DIAGNOSTIC_ERRORS', 'The mixture failed pre-flight checks.',
                               details=_issues(report))
            res = aiomfac_from_smiles(smiles, salts, w, T, strict='silent',
                                      bisulfate=bool(bisulfate))
            model = 'AIOMFAC (short-range + middle-range + long-range Debye-Hueckel)'
            neutral = [{'smiles': s, 'mass_fraction': wi,
                        'activity_coefficient': json_safe(g)}
                       for s, wi, g in zip(smiles, w, np.atleast_1d(res.gamma_neutral))]
            ions = [{'ion': _ion_label(i), 'kind': 'cation', 'activity_coefficient': json_safe(g)}
                    for i, g in zip(res.cation_ids, np.atleast_1d(res.gamma_cation))]
            ions += [{'ion': _ion_label(i), 'kind': 'anion', 'activity_coefficient': json_safe(g)}
                     for i, g in zip(res.anion_ids, np.atleast_1d(res.gamma_anion))]
            salt_rows = [{'salt': n, 'mass_fraction': wi,
                          'mean_ionic_activity_coefficient': json_safe(g)}
                         for n, wi, g in zip(res.salt_names, w[len(smiles):],
                                             np.atleast_1d(res.gamma_mean))]
            basis = 'mass_fraction'

    return {
        'model': model,
        'temperature_K': T,
        'composition_basis': basis,
        'bisulfate_equilibrium': bool(bisulfate) if salts else None,
        'neutral_components': neutral,
        'ions': ions,
        'salts': salt_rows,
        'reference_states': {
            'neutral_components': 'mole-fraction basis, pure-component (Raoult) reference',
            'ions': 'molality basis, infinite dilution in water reference',
        },
        'diagnostics': _issues(report),
        'warnings': notes,
    }


# ---------------------------------------------------------------------------
# Read-only file and database access, confined to AEROSYSPROP_DATA_DIR
# ---------------------------------------------------------------------------

def _data_root() -> Path:
    root = os.environ.get(DATA_DIR_ENV)
    if not root:
        raise ApiError(
            'FILE_ACCESS_DISABLED',
            'File access is disabled because no data directory is configured.',
            hint=f'Start the server with the {DATA_DIR_ENV} environment variable '
                 'set to the folder that holds your SMILES files and session databases.')
    path = Path(root).expanduser().resolve()
    if not path.is_dir():
        raise ApiError('FILE_ACCESS_DISABLED', f'{DATA_DIR_ENV} does not point to a directory.')
    return path


def _resolve_in_root(relative: str) -> Path:
    if not isinstance(relative, str) or not relative.strip():
        raise ApiError('INVALID_ARGUMENT', 'path must be a non-empty string.')
    root = _data_root()
    candidate = (root / relative).resolve()
    if candidate != root and root not in candidate.parents:
        raise ApiError('PATH_NOT_ALLOWED',
                       f'{relative!r} is outside the permitted data directory.',
                       hint='Use a path relative to the data directory, as returned by list_data_files.')
    if not candidate.is_file():
        raise ApiError('FILE_NOT_FOUND', f'No file named {relative!r} in the data directory.',
                       hint='Call list_data_files to see what is available.')
    return candidate


def list_data_files() -> dict:
    """List SMILES files and session databases inside the data directory."""
    root = _data_root()
    kinds = {**{s: 'smiles_text' for s in SMILES_FILE_SUFFIXES},
             **{s: 'sqlite' for s in DB_SUFFIXES},
             **{s: 'sdf' for s in SDF_SUFFIXES}}
    files = []
    for p in sorted(root.rglob('*')):
        rel = p.relative_to(root)
        if any(part.startswith('.') for part in rel.parts):
            continue
        if p.is_file() and p.suffix.lower() in kinds:
            files.append({'path': rel.as_posix(), 'kind': kinds[p.suffix.lower()],
                          'size_bytes': p.stat().st_size})
            if len(files) >= 500:
                break
    return {'data_directory': str(root), 'files': files}


def load_smiles_file(path: str, smiles_column: Optional[int] = None,
                     name_column: Optional[int] = None, table: str = 'compounds',
                     name_property: Optional[str] = None,
                     max_molecules: int = MAX_MOLECULES) -> dict:
    """Read molecules from a text, SQLite or SDF file in the data directory.

    For text files the SMILES column is detected automatically unless
    ``smiles_column`` is given (0-based).
    """
    p = _resolve_in_root(path)
    suffix = p.suffix.lower()
    with rdBase.BlockLogs(), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if suffix in SMILES_FILE_SUFFIXES:
            if smiles_column is None:
                with open(p) as fh:
                    smiles_column, detected_name, _ = asp_io.detect_columns(fh)
                if name_column is None:
                    name_column = detected_name if detected_name is not None else 0
            mols = asp_io.read_smiles_txt(str(p), smiles_col=int(smiles_column),
                                          name_col=int(name_column or 0))
        elif suffix in DB_SUFFIXES:
            try:
                mols = asp_io.read_smiles_db(str(p), table=table)
            except sqlite3.Error as e:
                raise ApiError('READ_FAILED', f'Could not read table {table!r}: {e}')
        elif suffix in SDF_SUFFIXES:
            mols = asp_io.read_smiles_sdf(str(p), name_property=name_property)
        else:
            raise ApiError('UNSUPPORTED_FORMAT', f'Unsupported file type {suffix!r}.',
                           hint='Supported: ' + ', '.join(sorted(
                               SMILES_FILE_SUFFIXES | DB_SUFFIXES | SDF_SUFFIXES)))
    molecules = [{'name': n, 'smiles': Chem.MolToSmiles(m)}
                 for n, m in mols[:max(0, int(max_molecules))]]
    return {'path': path, 'count': len(molecules), 'total_in_file': len(mols),
            'truncated': len(mols) > len(molecules), 'molecules': molecules,
            'note': 'Rows whose SMILES could not be parsed are skipped by the loader.'}


def read_session_database(path: str, what: str = 'results', methods=None,
                          temperature_K: Optional[float] = None,
                          mixture_id: Optional[int] = None, limit: int = 500) -> dict:
    """Read saved results from an AeroSysProp session database, read-only.

    ``what`` is one of ``'results'``, ``'compounds'``, ``'mixtures'`` or
    ``'mixture_species'`` (which needs ``mixture_id``).
    """
    p = _resolve_in_root(path)
    try:
        if what == 'results':
            cols = _as_str_list(methods, 'methods') if methods else None
            rows = db.load_results(str(p), methods=cols, temperature=temperature_K,
                                   read_only=True)
        elif what == 'compounds':
            rows = [{'name': n, 'smiles': s}
                    for n, s in db.list_compounds(str(p), read_only=True)]
        elif what == 'mixtures':
            rows = db.list_mixtures(str(p), read_only=True)
        elif what == 'mixture_species':
            if mixture_id is None:
                raise ApiError('INVALID_ARGUMENT', 'mixture_species needs a mixture_id.',
                               hint="Call with what='mixtures' first to list the IDs.")
            rows = db.load_mixture_species(str(p), int(mixture_id), read_only=True)
        else:
            raise ApiError('INVALID_ARGUMENT', f'Unknown value for what: {what!r}.',
                           details={'valid': ['results', 'compounds', 'mixtures',
                                              'mixture_species']})
    except sqlite3.DatabaseError as e:
        raise ApiError('READ_FAILED', f'Could not read {path!r}: {e}',
                       hint='The file may not be an AeroSysProp session database, '
                            'or may predate the table requested.')
    limit = max(0, int(limit))
    out = {'path': path, 'what': what, 'total_rows': len(rows),
           'truncated': len(rows) > limit, 'rows': json_safe(rows[:limit])}
    if what == 'results':
        cols = sorted({k for r in rows[:limit] for k in r} - {'name', 'smiles'})
        out['columns'] = {c: _unit_for(c) for c in cols}
    return out


# ---------------------------------------------------------------------------
# OPERA external service
# ---------------------------------------------------------------------------

def _opera_url() -> str:
    from aerosysprop.services import SERVICE_BY_ID
    return SERVICE_BY_ID['opera'].service_url().rstrip('/')


def _http_json(url: str, body=None, timeout: float = 5.0):
    req = urllib.request.Request(
        url, data=None if body is None else json.dumps(body).encode(),
        headers={'Content-Type': 'application/json'},
        method='GET' if body is None else 'POST')
    # The OPERA service is local; bypass any HTTP proxy set in the environment.
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def opera_predict(smiles, names=None, endpoints=None, timeout_s: float = 600) -> dict:
    """Run the OPERA QSAR models through the optional OPERA container."""
    from aerosysprop.services import _OPERA_ENDPOINTS
    smiles = _as_str_list(smiles, 'smiles')
    _check_count(len(smiles), MAX_OPERA_MOLECULES)
    if names is not None:
        names = _as_str_list(names, 'names')
        if len(names) != len(smiles):
            raise ApiError('INVALID_ARGUMENT', 'names must be the same length as smiles.')
    codes = {c.lower(): c for _, c, _ in _OPERA_ENDPOINTS}
    if endpoints:
        chosen = []
        for e in _as_str_list(endpoints, 'endpoints'):
            if e.lower() not in codes:
                raise ApiError('UNKNOWN_ENDPOINT', f'Unknown OPERA endpoint {e!r}.',
                               details={'valid_endpoints': {c: d for _, c, d in _OPERA_ENDPOINTS}})
            chosen.append(codes[e.lower()])
        endpoints = chosen
    else:
        endpoints = None
    base = _opera_url()
    unavailable = ApiError(
        'SERVICE_UNAVAILABLE', f'The OPERA service is not reachable at {base}.',
        hint='Start it with `docker compose up opera` (or set OPERA_URL), then retry. '
             "AeroSysProp's own methods do not need OPERA.")
    try:
        health = _http_json(base + '/health', timeout=3)
    except (urllib.error.URLError, OSError, ValueError):
        raise unavailable
    if not health.get('opera_available', True):
        raise ApiError('SERVICE_UNAVAILABLE',
                       'The OPERA container is running but the OPERA program is not '
                       'available inside it.')
    body = {'smiles': smiles, 'names': names, 'endpoints': endpoints, 'format': 'json'}
    try:
        res = _http_json(base + '/predict', body, timeout=timeout_s)
    except urllib.error.HTTPError as e:
        raise ApiError('SERVICE_ERROR',
                       f'OPERA returned HTTP {e.code}: {e.read().decode(errors="replace")[:500]}')
    except (urllib.error.URLError, OSError):
        raise unavailable
    return {'service': 'OPERA', 'service_url': base,
            'endpoints': endpoints or 'all',
            'n': res.get('n'), 'columns': res.get('columns'),
            'records': json_safe(res.get('records')),
            'note': 'OPERA values are QSAR model predictions with their own units '
                    '(see list of endpoints); they are independent of '
                    "AeroSysProp's group-contribution estimates."}
