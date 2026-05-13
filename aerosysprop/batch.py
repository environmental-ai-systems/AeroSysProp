"""
aerosysprop.batch
=================
Parallel batch processing of large compound sets.

All prediction methods are unchanged — this module only controls *how* they
are dispatched across molecules, not what they compute.

Typical usage
-------------
>>> from aerosysprop.io import read_smiles_txt
>>> from aerosysprop.batch import run_all, ALL_METHODS
>>>
>>> mols = read_smiles_txt('compounds.txt', smiles_col=0)
>>> results = run_all(mols, methods=['bp_nannoolal', 'vp_simpol'],
...                   temperature=298.15, output_csv='results.csv')

Each element of *results* is a dict keyed by output column name.
Failed predictions for individual methods are stored as ``None``.
"""
import csv
import math
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from rdkit import Chem
from tqdm import tqdm

_MP_CONTEXT = multiprocessing.get_context(
    'fork' if sys.platform != 'win32' else 'spawn'
)

# ---------------------------------------------------------------------------
# Method registry
# ---------------------------------------------------------------------------

# (display label, unit string, property group, direct dependencies)
METHOD_META = {
    'bp_joback':         ('Joback & Reid (1987)',            'K',           'Boiling Point',       frozenset()),
    'bp_stein':          ('Stein & Brown (1994)',             'K',           'Boiling Point',       frozenset()),
    'bp_nannoolal':      ('Nannoolal et al. (2004)',          'K',           'Boiling Point',       frozenset()),
    'vp_nannoolal':      ('Nannoolal et al. (2008)',          'log₁₀(atm)', 'Vapour Pressure',     frozenset({'bp_nannoolal'})),
    'vp_myrdal':         ('Myrdal & Yalkowsky (1997)',        'log₁₀(atm)', 'Vapour Pressure',     frozenset({'bp_nannoolal'})),
    'vp_evaporation':    ('EVAPORATION (Compernolle 2011)',   'log₁₀(atm)', 'Vapour Pressure',     frozenset()),
    'vp_evaporation2':   ('EVAPORATION revised',             'log₁₀(atm)', 'Vapour Pressure',     frozenset()),
    'vp_simpol':         ('SIMPOL.1 (Pankow & Asher 2008)',  'log₁₀(atm)', 'Vapour Pressure',     frozenset()),
    'cp_nannoolal':      ('Nannoolal et al. (2004)',          'K / bar',     'Critical Properties', frozenset({'bp_nannoolal'})),
    'density_girolami':  ('Girolami (1994)',                  'g/cc',        'Liquid Density',      frozenset()),
    'density_schroeder': ('Schroeder',                        'g/cc',        'Liquid Density',      frozenset({'cp_nannoolal'})),
    'density_le_bas':    ('Le Bas',                           'g/cc',        'Liquid Density',      frozenset({'cp_nannoolal'})),
    'density_tyn_calus': ('Tyn & Calus',                     'g/cc',        'Liquid Density',      frozenset({'cp_nannoolal'})),
}

# Canonical ordering for display and CSV columns
ALL_METHODS = [
    'bp_joback', 'bp_stein', 'bp_nannoolal',
    'vp_nannoolal', 'vp_myrdal', 'vp_evaporation', 'vp_evaporation2', 'vp_simpol',
    'cp_nannoolal',
    'density_girolami', 'density_schroeder', 'density_le_bas', 'density_tyn_calus',
]

# Methods grouped for UI display
METHOD_GROUPS = {
    'Boiling Point':       ['bp_joback', 'bp_stein', 'bp_nannoolal'],
    'Vapour Pressure':     ['vp_nannoolal', 'vp_myrdal', 'vp_evaporation', 'vp_evaporation2', 'vp_simpol'],
    'Critical Properties': ['cp_nannoolal'],
    'Liquid Density':      ['density_girolami', 'density_schroeder', 'density_le_bas', 'density_tyn_calus'],
}

# Output column names produced by each method
_OUTPUT_KEYS = {
    'bp_joback':         ('bp_joback_K',),
    'bp_stein':          ('bp_stein_K',),
    'bp_nannoolal':      ('bp_nannoolal_K',),
    'vp_nannoolal':      ('vp_nannoolal_log10atm',),
    'vp_myrdal':         ('vp_myrdal_log10atm',),
    'vp_evaporation':    ('vp_evaporation_log10atm',),
    'vp_evaporation2':   ('vp_evaporation2_log10atm',),
    'vp_simpol':         ('vp_simpol_log10atm',),
    'cp_nannoolal':      ('tc_nannoolal_K', 'pc_nannoolal_bar'),
    'density_girolami':  ('density_girolami_gcc',),
    'density_schroeder': ('density_schroeder_gcc',),
    'density_le_bas':    ('density_le_bas_gcc',),
    'density_tyn_calus': ('density_tyn_calus_gcc',),
}

_CSV_FIELDS = ['name', 'smiles'] + [
    col for m in ALL_METHODS for col in _OUTPUT_KEYS[m]
]


def resolve_methods(requested):
    """
    Expand a list of method IDs to include all transitive dependencies.

    Dependencies are computed internally but included in the returned set so
    callers know what outputs to expect.

    Parameters
    ----------
    requested : iterable of str
        Method IDs from ``ALL_METHODS``.

    Returns
    -------
    set of str
        The requested methods plus any dependencies needed to compute them.
    """
    needed = set(requested)
    changed = True
    while changed:
        changed = False
        for m in list(needed):
            for dep in METHOD_META[m][3]:
                if dep not in needed:
                    needed.add(dep)
                    changed = True
    return needed


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def _safe(val):
    if isinstance(val, complex):
        val = val.real
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return None
    return val


def _run_one(args):
    """
    Compute selected property predictions for a single molecule.

    Args are (name, smiles, temperature, frozenset_of_method_ids) — plain
    types so they pickle safely across process boundaries.
    """
    name, smiles, temperature, methods = args

    from aerosysprop.methods import (
        boiling_points,
        critical_properties,
        liquid_densities,
        vapour_pressures,
    )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return name, None

    out = {'name': name, 'smiles': smiles}

    def _call(fn, *a):
        try:
            return _safe(fn(*a))
        except Exception:
            return None

    # Boiling points
    bp_j = _call(boiling_points.joback_and_reid, mol) if 'bp_joback' in methods else None
    bp_s = _call(boiling_points.stein_and_brown, mol)  if 'bp_stein' in methods else None
    bp_n = _call(boiling_points.nannoolal, mol)         if 'bp_nannoolal' in methods else None

    if 'bp_joback'    in methods: out['bp_joback_K']  = bp_j
    if 'bp_stein'     in methods: out['bp_stein_K']   = bp_s
    if 'bp_nannoolal' in methods: out['bp_nannoolal_K'] = bp_n

    # Vapour pressures
    if 'vp_nannoolal' in methods:
        out['vp_nannoolal_log10atm'] = _call(vapour_pressures.nannoolal, mol, temperature, bp_n) if bp_n else None
    if 'vp_myrdal' in methods:
        out['vp_myrdal_log10atm'] = _call(vapour_pressures.myrdal_and_yalkowsky, mol, temperature, bp_n) if bp_n else None
    if 'vp_evaporation'  in methods: out['vp_evaporation_log10atm']  = _call(vapour_pressures.evaporation,  mol, temperature)
    if 'vp_evaporation2' in methods: out['vp_evaporation2_log10atm'] = _call(vapour_pressures.evaporation2, mol, temperature)
    if 'vp_simpol'       in methods: out['vp_simpol_log10atm']       = _call(vapour_pressures.simpol,       mol, temperature)

    # Critical properties (Nannoolal only — depends on bp_nannoolal)
    cp_n = None
    if 'cp_nannoolal' in methods:
        if bp_n is not None:
            try:
                cp_n = critical_properties.nannoolal(mol, bp_n)
                out['tc_nannoolal_K']   = _safe(cp_n.temperature)
                out['pc_nannoolal_bar'] = _safe(cp_n.pressure)
            except Exception:
                out['tc_nannoolal_K']   = None
                out['pc_nannoolal_bar'] = None
        else:
            out['tc_nannoolal_K']   = None
            out['pc_nannoolal_bar'] = None

    # Liquid densities
    if 'density_girolami' in methods:
        try:
            out['density_girolami_gcc'] = _safe(float(liquid_densities.girolami(mol).real))
        except Exception:
            out['density_girolami_gcc'] = None

    # Densities that need critical properties — resolve cp_n if not already done
    needs_cp = any(m in methods for m in ('density_schroeder', 'density_le_bas', 'density_tyn_calus'))
    if needs_cp and cp_n is None and bp_n is not None and 'cp_nannoolal' not in methods:
        try:
            cp_n = critical_properties.nannoolal(mol, bp_n)
        except Exception:
            cp_n = None

    for method_id, fn, key in (
        ('density_schroeder', liquid_densities.schroeder, 'density_schroeder_gcc'),
        ('density_le_bas',    liquid_densities.le_bas,    'density_le_bas_gcc'),
        ('density_tyn_calus', liquid_densities.tyn_and_calus, 'density_tyn_calus_gcc'),
    ):
        if method_id in methods:
            if cp_n is not None:
                try:
                    out[key] = _safe(float(fn(mol, temperature, cp_n).real))
                except Exception:
                    out[key] = None
            else:
                out[key] = None

    return name, out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_all(mols, methods=None, temperature=298.15, workers=None,
            output_csv=None, progress_callback=None, use_threads=False):
    """
    Run property predictions on a list of molecules in parallel.

    Parameters
    ----------
    mols : list of (str, rdkit.Chem.Mol)
        (name, mol) pairs — the direct output of any ``aerosysprop.io`` loader.
    methods : list of str or None
        Method IDs to compute (see ``ALL_METHODS``).  Dependencies are
        resolved automatically — e.g. requesting ``'vp_nannoolal'`` will
        also compute ``'bp_nannoolal'`` internally.  ``None`` runs all methods.
    temperature : float
        Temperature in K for vapour pressure and density calculations.
        Default 298.15 K.
    workers : int or None
        Number of parallel workers.  ``None`` uses all available CPU cores.
    output_csv : str or None
        If provided, results are written to this CSV file incrementally.
    progress_callback : callable or None
        Called as ``progress_callback(done, total)`` after each molecule
        completes.  Use this to update a UI progress bar.
    use_threads : bool
        Use ``ThreadPoolExecutor`` instead of ``ProcessPoolExecutor``.
        Recommended when calling from a multi-threaded host (e.g. Streamlit).

    Returns
    -------
    list of dict
        One dict per molecule in input order.  Keys are output column names
        (e.g. ``'bp_nannoolal_K'``).  Failed predictions are ``None``.

    Examples
    --------
    >>> results = run_all(mols, methods=['bp_nannoolal', 'vp_simpol'])
    >>> results = run_all(mols, output_csv='out.csv', progress_callback=print)
    """
    if methods is None:
        methods_to_run = frozenset(ALL_METHODS)
    else:
        unknown = set(methods) - set(ALL_METHODS)
        if unknown:
            raise ValueError('Unknown method IDs: %s' % sorted(unknown))
        methods_to_run = frozenset(resolve_methods(methods))

    if workers is None:
        workers = os.cpu_count() or 1

    args = [
        (name, Chem.MolToSmiles(mol), temperature, methods_to_run)
        for name, mol in mols
    ]

    results_by_name = {}
    csv_file = None
    writer = None

    if output_csv:
        requested_fields = ['name', 'smiles'] + [
            col for m in ALL_METHODS
            if m in methods_to_run
            for col in _OUTPUT_KEYS[m]
        ]
        csv_file = open(output_csv, 'w', newline='')
        writer = csv.DictWriter(csv_file, fieldnames=requested_fields,
                                extrasaction='ignore')
        writer.writeheader()

    Executor = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    pool_kwargs = {} if use_threads else {'mp_context': _MP_CONTEXT}

    try:
        with Executor(max_workers=workers, **pool_kwargs) as pool:
            futures = {pool.submit(_run_one, a): a[0] for a in args}
            total = len(futures)
            done_count = 0
            with tqdm(total=total, desc='Processing', unit='mol',
                      disable=progress_callback is not None) as bar:
                for future in as_completed(futures):
                    name, result = future.result()
                    results_by_name[name] = result
                    done_count += 1
                    if writer and result:
                        writer.writerow(result)
                        csv_file.flush()
                    bar.update(1)
                    if progress_callback:
                        progress_callback(done_count, total)
    finally:
        if csv_file:
            csv_file.close()

    return [results_by_name.get(name) for name, _ in mols]
