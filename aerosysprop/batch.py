"""
aerosysprop.batch
=================
Parallel batch processing of large compound sets.

All prediction methods are unchanged — this module only controls *how* they
are dispatched across molecules, not what they compute.

Typical usage
-------------
>>> from aerosysprop.io import read_smiles_txt
>>> from aerosysprop.batch import run_all
>>>
>>> mols = read_smiles_txt('compounds.txt', smiles_col=0)
>>> results = run_all(mols, temperature=298.15, output_csv='results.csv')

Each element of *results* is a dict with keys like ``bp_nannoolal_K``,
``vp_simpol_log10atm``, etc.  Failed predictions for individual methods are
stored as ``None`` rather than raising.
"""
import csv
import math
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

from rdkit import Chem
from tqdm import tqdm

# macOS and Windows default to 'spawn', which re-imports the calling script and
# breaks when run outside of an if __name__ == '__main__' guard.  'fork' is
# safe on Linux/macOS for this workload and sidesteps the bootstrapping issue.
# On Windows 'fork' is unavailable so we fall back to 'spawn'.
_MP_CONTEXT = multiprocessing.get_context(
    'fork' if sys.platform != 'win32' else 'spawn'
)

# ---------------------------------------------------------------------------
# Worker (runs in a subprocess — must be importable at module level)
# ---------------------------------------------------------------------------

def _safe(val):
    """Return None for NaN/inf/complex imaginary part so CSV stays clean."""
    if isinstance(val, complex):
        val = val.real
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return None
    return val


def _run_one(args):
    """
    Compute all property predictions for a single molecule.

    Receives a (name, smiles, temperature) tuple so the arguments are
    plain strings/floats — safe to pickle across process boundaries.
    Returns (name, result_dict) where failed methods appear as None.
    """
    name, smiles, temperature = args

    # Import inside the worker so each subprocess initialises its own copy
    # of the aerosysprop data tables (avoids shared-memory issues).
    from aerosysprop.methods import (
        boiling_points,
        critical_properties,
        liquid_densities,
        vapour_pressures,
    )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return name, None

    def _call(fn, *a):
        try:
            return _safe(fn(*a))
        except Exception:
            return None

    # --- Boiling points ---
    bp_j = _call(boiling_points.joback_and_reid, mol)
    bp_s = _call(boiling_points.stein_and_brown, mol)
    bp_n = _call(boiling_points.nannoolal, mol)

    out = {
        'name':              name,
        'smiles':            smiles,
        'bp_joback_K':       bp_j,
        'bp_stein_K':        bp_s,
        'bp_nannoolal_K':    bp_n,
    }

    # --- Vapour pressures that need a boiling point (use Nannoolal BP) ---
    out['vp_nannoolal_log10atm']  = _call(vapour_pressures.nannoolal,           mol, temperature, bp_n) if bp_n else None
    out['vp_myrdal_log10atm']     = _call(vapour_pressures.myrdal_and_yalkowsky, mol, temperature, bp_n) if bp_n else None

    # --- Vapour pressures that don't need a boiling point ---
    out['vp_evaporation_log10atm']  = _call(vapour_pressures.evaporation,  mol, temperature)
    out['vp_evaporation2_log10atm'] = _call(vapour_pressures.evaporation2, mol, temperature)
    out['vp_simpol_log10atm']       = _call(vapour_pressures.simpol,       mol, temperature)

    # --- Critical properties (Nannoolal, requires Nannoolal BP) ---
    cp_n = None
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

    # --- Liquid densities ---
    try:
        out['density_girolami_gcc'] = _safe(float(liquid_densities.girolami(mol).real))
    except Exception:
        out['density_girolami_gcc'] = None

    if cp_n is not None:
        try:
            out['density_schroeder_gcc'] = _safe(float(liquid_densities.schroeder(mol, temperature, cp_n).real))
        except Exception:
            out['density_schroeder_gcc'] = None
        try:
            out['density_le_bas_gcc'] = _safe(float(liquid_densities.le_bas(mol, temperature, cp_n).real))
        except Exception:
            out['density_le_bas_gcc'] = None
        try:
            out['density_tyn_calus_gcc'] = _safe(float(liquid_densities.tyn_and_calus(mol, temperature, cp_n).real))
        except Exception:
            out['density_tyn_calus_gcc'] = None
    else:
        out['density_schroeder_gcc'] = None
        out['density_le_bas_gcc']    = None
        out['density_tyn_calus_gcc'] = None

    return name, out


# Column order for CSV output
_CSV_FIELDS = [
    'name', 'smiles',
    'bp_joback_K', 'bp_stein_K', 'bp_nannoolal_K',
    'vp_nannoolal_log10atm', 'vp_myrdal_log10atm',
    'vp_evaporation_log10atm', 'vp_evaporation2_log10atm', 'vp_simpol_log10atm',
    'tc_nannoolal_K', 'pc_nannoolal_bar',
    'density_girolami_gcc', 'density_schroeder_gcc',
    'density_le_bas_gcc', 'density_tyn_calus_gcc',
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_all(mols, temperature=298.15, workers=None, output_csv=None):
    """
    Run all property predictions on a list of molecules in parallel.

    Parameters
    ----------
    mols : list of (str, rdkit.Chem.Mol)
        (name, mol) pairs — the direct output of any ``aerosysprop.io`` loader.
    temperature : float
        Temperature in K for vapour pressure and density calculations.
        Default 298.15 K.
    workers : int or None
        Number of parallel worker processes.  ``None`` uses all available CPU
        cores (``os.cpu_count()``).
    output_csv : str or None
        If provided, results are written to this CSV file as they complete.

    Returns
    -------
    list of dict
        One dict per molecule, in the same order as *mols*.  Keys include
        ``name``, ``smiles``, and one entry per predicted property.
        Failed predictions are stored as ``None``.

    Examples
    --------
    >>> from aerosysprop.io import read_smiles_txt
    >>> from aerosysprop.batch import run_all
    >>> mols = read_smiles_txt('compounds.txt', smiles_col=0)
    >>> results = run_all(mols, output_csv='results.csv')
    """
    if workers is None:
        workers = os.cpu_count() or 1

    # Convert mol objects to SMILES for pickling across process boundaries
    args = [
        (name, Chem.MolToSmiles(mol), temperature)
        for name, mol in mols
    ]

    results_by_name = {}
    csv_file = None
    writer = None

    if output_csv:
        csv_file = open(output_csv, 'w', newline='')
        writer = csv.DictWriter(csv_file, fieldnames=_CSV_FIELDS, extrasaction='ignore')
        writer.writeheader()

    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=_MP_CONTEXT) as pool:
            futures = {pool.submit(_run_one, a): a[0] for a in args}
            with tqdm(total=len(futures), desc='Processing', unit='mol') as bar:
                for future in as_completed(futures):
                    name, result = future.result()
                    results_by_name[name] = result
                    if writer and result:
                        writer.writerow(result)
                        csv_file.flush()
                    bar.update(1)
    finally:
        if csv_file:
            csv_file.close()

    # Return in original input order
    return [results_by_name.get(name) for name, _ in mols]
