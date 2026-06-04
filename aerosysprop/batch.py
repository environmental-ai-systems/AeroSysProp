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

# Sentinel used in METHOD_META deps to mean "whichever boiling point method
# the user has chosen" — resolved to the actual method ID at runtime.
BP_ANY = 'bp_any'

# The three boiling point methods, for validation
BP_METHODS = frozenset({'bp_joback', 'bp_stein', 'bp_nannoolal'})

# (display label, unit string, property group, direct dependencies)
# Use BP_ANY for any method that needs *a* boiling point — not a specific one.
METHOD_META = {
    'bp_joback':         ('Joback & Reid (1987)',            'K',           'Boiling Point',       frozenset()),
    'bp_stein':          ('Stein & Brown (1994)',             'K',           'Boiling Point',       frozenset()),
    'bp_nannoolal':      ('Nannoolal et al. (2004)',          'K',           'Boiling Point',       frozenset()),
    'vp_nannoolal':      ('Nannoolal et al. (2008)',          'log₁₀(atm)', 'Vapour Pressure',     frozenset({BP_ANY})),
    'vp_myrdal':         ('Myrdal & Yalkowsky (1997)',        'log₁₀(atm)', 'Vapour Pressure',     frozenset({BP_ANY})),
    'vp_evaporation':    ('EVAPORATION (Compernolle 2011)',   'log₁₀(atm)', 'Vapour Pressure',     frozenset()),
    'vp_evaporation2':   ('EVAPORATION revised',             'log₁₀(atm)', 'Vapour Pressure',     frozenset()),
    'vp_simpol':         ('SIMPOL.1 (Pankow & Asher 2008)',  'log₁₀(atm)', 'Vapour Pressure',     frozenset()),
    'cp_nannoolal':      ('Nannoolal et al. (2004)',          'K / bar',     'Critical Properties', frozenset({BP_ANY})),
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

# Short labels used to suffix output columns when a calculation is run across
# more than one boiling-point source.
_BP_SHORT = {
    'bp_joback': 'Joback',
    'bp_stein': 'Stein',
    'bp_nannoolal': 'Nannoolal',
}

# Methods whose output depends on the chosen boiling-point source. When more
# than one source is selected, their columns are produced once per source with
# a '__bp<Source>' suffix; with a single source the column names are unchanged.
_BP_DEPENDENT_KEYS = {
    'vp_nannoolal':      ('vp_nannoolal_log10atm',),
    'vp_myrdal':         ('vp_myrdal_log10atm',),
    'cp_nannoolal':      ('tc_nannoolal_K', 'pc_nannoolal_bar'),
    'density_schroeder': ('density_schroeder_gcc',),
    'density_le_bas':    ('density_le_bas_gcc',),
    'density_tyn_calus': ('density_tyn_calus_gcc',),
}


def resolved_output_columns(methods_to_run, bp_sources):
    """Return the ordered output-column names for a run, accounting for the
    boiling-point-source suffix when multiple sources are selected."""
    multi = len(bp_sources) > 1
    cols = []
    for m in ALL_METHODS:
        if m not in methods_to_run:
            continue
        if m in _BP_DEPENDENT_KEYS:
            for key in _BP_DEPENDENT_KEYS[m]:
                if multi:
                    cols.extend(key + '__bp' + _BP_SHORT[b] for b in bp_sources)
                else:
                    cols.append(key)
        else:
            cols.extend(_OUTPUT_KEYS[m])
    return cols


def resolve_methods(requested, bp_source='bp_nannoolal', bp_sources=None):
    """
    Expand a list of method IDs to include all transitive dependencies.

    Parameters
    ----------
    requested : iterable of str
        Method IDs from ``ALL_METHODS``.
    bp_source : str
        Single boiling-point method used as input for vapour pressure,
        critical property, and density methods that need one. Must be one of
        ``'bp_joback'``, ``'bp_stein'``, ``'bp_nannoolal'``. Default
        ``'bp_nannoolal'``. Ignored if ``bp_sources`` is given.
    bp_sources : sequence of str or None
        Multiple boiling-point sources to run the dependent methods across.
        When given, each ``BP_ANY`` dependency resolves to all of these.
        ``None`` falls back to ``[bp_source]``.

    Returns
    -------
    set of str
        The requested methods plus any dependencies needed to compute them.
    """
    if bp_sources is None:
        bp_sources = [bp_source]
    for b in bp_sources:
        if b not in BP_METHODS:
            raise ValueError('bp_source must be one of %s' % sorted(BP_METHODS))
    needed = set(requested)
    changed = True
    while changed:
        changed = False
        for m in list(needed):
            for dep in METHOD_META[m][3]:
                actuals = bp_sources if dep == BP_ANY else [dep]
                for actual in actuals:
                    if actual not in needed:
                        needed.add(actual)
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

    Args are (name, smiles, temperature, frozenset_of_method_ids,
    tuple_of_bp_sources) — plain types so they pickle safely across process
    boundaries. Methods that depend on a boiling point are computed once per
    boiling-point source; when more than one source is given the resulting
    columns get a '__bp<Source>' suffix.
    """
    name, smiles, temperature, methods, bp_sources = args

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

    # --- Boiling points (each computed once; needed as output and/or input) ---
    bp_vals = {}
    for bp_id, fn, key in (
        ('bp_joback',    boiling_points.joback_and_reid, 'bp_joback_K'),
        ('bp_stein',     boiling_points.stein_and_brown, 'bp_stein_K'),
        ('bp_nannoolal', boiling_points.nannoolal,       'bp_nannoolal_K'),
    ):
        if bp_id in methods or bp_id in bp_sources:
            bp_vals[bp_id] = _call(fn, mol)
        if bp_id in methods:
            out[key] = bp_vals.get(bp_id)

    # --- Boiling-point-independent vapour pressures / density ---
    if 'vp_evaporation'  in methods: out['vp_evaporation_log10atm']  = _call(vapour_pressures.evaporation,  mol, temperature)
    if 'vp_evaporation2' in methods: out['vp_evaporation2_log10atm'] = _call(vapour_pressures.evaporation2, mol, temperature)
    if 'vp_simpol'       in methods: out['vp_simpol_log10atm']       = _call(vapour_pressures.simpol,       mol, temperature)
    if 'density_girolami' in methods:
        try:
            out['density_girolami_gcc'] = _safe(float(liquid_densities.girolami(mol).real))
        except Exception:
            out['density_girolami_gcc'] = None

    # --- Boiling-point-dependent block: once per boiling-point source ---
    multi = len(bp_sources) > 1
    needs_cp = any(m in methods for m in
                   ('cp_nannoolal', 'density_schroeder', 'density_le_bas', 'density_tyn_calus'))
    for bp_src in bp_sources:
        suffix = ('__bp' + _BP_SHORT[bp_src]) if multi else ''
        bp_for = bp_vals.get(bp_src)

        if 'vp_nannoolal' in methods:
            out['vp_nannoolal_log10atm' + suffix] = (
                _call(vapour_pressures.nannoolal, mol, temperature, bp_for)
                if bp_for else None
            )
        if 'vp_myrdal' in methods:
            out['vp_myrdal_log10atm' + suffix] = (
                _call(vapour_pressures.myrdal_and_yalkowsky, mol, temperature, bp_for)
                if bp_for else None
            )

        # Critical properties from THIS boiling-point source (also used for
        # the density methods below, for consistency).
        cp_n = None
        if needs_cp and bp_for is not None:
            try:
                cp_n = critical_properties.nannoolal(mol, bp_for)
            except Exception:
                cp_n = None

        if 'cp_nannoolal' in methods:
            out['tc_nannoolal_K'   + suffix] = _safe(cp_n.temperature) if cp_n else None
            out['pc_nannoolal_bar' + suffix] = _safe(cp_n.pressure)    if cp_n else None

        for method_id, fn, key in (
            ('density_schroeder', liquid_densities.schroeder,     'density_schroeder_gcc'),
            ('density_le_bas',    liquid_densities.le_bas,        'density_le_bas_gcc'),
            ('density_tyn_calus', liquid_densities.tyn_and_calus, 'density_tyn_calus_gcc'),
        ):
            if method_id in methods:
                if cp_n is not None:
                    try:
                        out[key + suffix] = _safe(float(fn(mol, temperature, cp_n).real))
                    except Exception:
                        out[key + suffix] = None
                else:
                    out[key + suffix] = None

    return name, out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_all(mols, methods=None, temperature=298.15, workers=None,
            output_csv=None, progress_callback=None, use_threads=False,
            bp_source='bp_nannoolal', bp_sources=None):
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
    bp_source : str
        Single boiling-point method used as input for vapour pressure,
        critical property, and density methods that need one. One of
        ``'bp_joback'``, ``'bp_stein'``, ``'bp_nannoolal'``. Default
        ``'bp_nannoolal'``. Ignored if ``bp_sources`` is given.
    bp_sources : sequence of str or None
        Multiple boiling-point sources to run the dependent methods across.
        With one source the output column names are unchanged; with several,
        each dependent column is produced once per source with a
        ``'__bp<Source>'`` suffix (e.g. ``'vp_nannoolal_log10atm__bpStein'``).
        ``None`` falls back to ``[bp_source]``.

    Returns
    -------
    list of dict
        One dict per molecule in input order.  Keys are output column names
        (e.g. ``'bp_nannoolal_K'``).  Failed predictions are ``None``.

    Examples
    --------
    >>> results = run_all(mols, methods=['bp_nannoolal', 'vp_simpol'])
    >>> results = run_all(mols, methods=['vp_nannoolal'],
    ...                   bp_sources=['bp_nannoolal', 'bp_stein'])
    """
    # Normalise the boiling-point source list (dedupe, preserve order, validate)
    if bp_sources is None:
        bp_sources = [bp_source]
    bp_sources = list(dict.fromkeys(bp_sources))
    for b in bp_sources:
        if b not in BP_METHODS:
            raise ValueError('bp_source must be one of %s' % sorted(BP_METHODS))

    if methods is None:
        methods_to_run = frozenset(resolve_methods(ALL_METHODS, bp_sources=bp_sources))
    else:
        unknown = set(methods) - set(ALL_METHODS)
        if unknown:
            raise ValueError('Unknown method IDs: %s' % sorted(unknown))
        methods_to_run = frozenset(resolve_methods(methods, bp_sources=bp_sources))

    if workers is None:
        workers = os.cpu_count() or 1

    bp_sources_t = tuple(bp_sources)
    args = [
        (name, Chem.MolToSmiles(mol), temperature, methods_to_run, bp_sources_t)
        for name, mol in mols
    ]

    results_by_name = {}
    csv_file = None
    writer = None

    if output_csv:
        requested_fields = ['name', 'smiles'] + resolved_output_columns(
            methods_to_run, bp_sources)
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
