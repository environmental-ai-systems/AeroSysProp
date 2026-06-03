# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2026 David Topping.
# All Rights Reserved.
# This file is part of aerosysprop.
"""Batch SMARTS-hit fingerprints for SMILES inspection.

Given a list of SMILES strings, returns a matrix of how many times each
SMARTS pattern in a chosen method's pattern set matches each molecule.
Useful for:

  - Sanity-checking matcher decompositions across a batch of compounds
  - Quick functional-group "fingerprint" comparison
  - Pre-filtering compounds by group presence (e.g. "which of these
    have any nitrate-related pattern?")
  - Diagnosing matcher coverage on a custom SMILES list

Public API
----------
- :func:`smarts_hits` -- hits matrix for one SMARTS-method name
- :func:`smarts_hits_all` -- side-by-side hits across many methods
- :func:`available_methods` -- registered method-name keys

Example
-------
>>> from aerosysprop.fingerprints import smarts_hits, available_methods
>>> available_methods()
['aiomfac', 'aiomfac_hctail', 'aiomfac_ions', 'composition', ...]

>>> df = smarts_hits(['O', 'CCO', 'Cc1ccccc1'], method='aiomfac')
>>> df.loc['CCO']        # rows are the SMILES
1      1
2      1
...
14     1
"""

import numpy as np
from rdkit import Chem

from aerosysprop import data


# Method name -> attribute on aerosysprop.data holding a {id: smarts}
# dictionary. Each entry is a self-contained SMARTS family that produces
# a sensible hit-count matrix.
_SMARTS_PROVIDERS = {
    # AIOMFAC sets
    'aiomfac':              'AIOMFAC_UNIFAC_SMARTS',
    'aiomfac_hctail':       'AIOMFAC_HCTAIL_SMARTS',
    'aiomfac_ions':         'AIOMFAC_ION_SMARTS',
    # Boiling-point / VP / Tc methods
    'nannoolal_primary':    'NANNOOLAL_SMARTS_PRIMARY',
    'nannoolal_secondary':  'NANNOOLAL_SMARTS_SECONDARY',
    'stein_and_brown':      'STEIN_AND_BROWN_SMARTS',
    # Vapour-pressure-specific
    'simpol':               'SIMPOL_SMARTS',
    'evaporation':          'EVAPORATION_SMARTS',
    'evaporation2':         'EVAPORATION2_SMARTS',
    # Density methods
    'girolami':             'GIROLAMI_SMARTS',
    'le_bas':               'LE_BAS_SMARTS',
    'schroeder':            'SCHROEDER_SMARTS',
    # General
    'composition':          'COMPOSITION_SMARTS',
    'moller':               'MOLLER_SMARTS',
    # Fingerprint sets
    'fp4':                  'FP4_SMARTS',
    'maccs':                'MACCS_SMARTS',
    # Functional-group SMARTS used by specific methods that we expose for
    # introspection
    'delplatt':             None,   # placeholder if not loaded
}


def available_methods():
    """Return the list of registered SMARTS-method names."""
    return sorted(
        name for name, attr in _SMARTS_PROVIDERS.items()
        if attr is not None and hasattr(data, attr)
    )


def _resolve_smarts(method):
    if method not in _SMARTS_PROVIDERS:
        raise ValueError(
            f'Unknown method {method!r}. Available: {available_methods()}'
        )
    attr = _SMARTS_PROVIDERS[method]
    if attr is None or not hasattr(data, attr):
        raise ValueError(
            f'Method {method!r} is registered but its SMARTS data '
            f'attribute {attr!r} is not loaded.'
        )
    return getattr(data, attr)


def _hits_for_one(mol, smarts_dict):
    """Run every SMARTS in `smarts_dict` over `mol`, return dict {id: count}."""
    if mol is None:
        return {}
    return {sid: len(smarts_obj.findall(mol))
            for sid, smarts_obj in smarts_dict.items()}


def smarts_hits(smiles, method='aiomfac', as_dataframe=True,
                drop_zero_columns=False, index_label='SMILES'):
    """Return a hit-count matrix for a batch of SMILES against one method's
    SMARTS set.

    Parameters
    ----------
    smiles : sequence of str
        SMILES strings to evaluate. Invalid SMILES yield an all-zero row.
    method : str, default 'aiomfac'
        Name of the SMARTS-method to run. See :func:`available_methods`.
    as_dataframe : bool, default True
        If True, return a pandas ``DataFrame`` (rows = smiles, cols =
        SMARTS IDs sorted ascending, values = hit counts). If False,
        return a ``(numpy.ndarray, list_of_smiles, list_of_ids)`` tuple.
    drop_zero_columns : bool, default False
        If True, drop SMARTS IDs that are zero across every SMILES in the
        batch — useful for large pattern sets where most patterns won't
        match the batch at all.
    index_label : str, default 'SMILES'
        Name for the row index (when as_dataframe=True).

    Returns
    -------
    pandas.DataFrame or tuple
        See ``as_dataframe``.
    """
    smarts_dict = _resolve_smarts(method)
    sorted_ids = sorted(smarts_dict.keys())
    n_smi = len(smiles)
    n_pat = len(sorted_ids)
    matrix = np.zeros((n_smi, n_pat), dtype=int)
    col_idx = {sid: i for i, sid in enumerate(sorted_ids)}

    for r, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        for sid, smarts_obj in smarts_dict.items():
            n = len(smarts_obj.findall(mol))
            if n > 0:
                matrix[r, col_idx[sid]] = n

    if drop_zero_columns:
        mask = (matrix.sum(axis=0) > 0)
        matrix = matrix[:, mask]
        sorted_ids = [sid for sid, keep in zip(sorted_ids, mask) if keep]

    if as_dataframe:
        import pandas as pd
        df = pd.DataFrame(matrix, index=list(smiles), columns=sorted_ids)
        df.index.name = index_label
        df.columns.name = method
        return df
    return matrix, list(smiles), sorted_ids


def smarts_hits_all(smiles, methods=None, as_dataframe=True,
                    drop_zero_columns=False):
    """Run :func:`smarts_hits` across multiple methods and concatenate the
    results side-by-side with a ``MultiIndex(method, id)`` on the columns.

    Parameters
    ----------
    smiles : sequence of str
    methods : sequence of str, optional
        Methods to run; defaults to all in :func:`available_methods`.
    as_dataframe : bool, default True
    drop_zero_columns : bool, default False
        If True, drop columns that are zero across every SMILES in the
        batch (applied per-method).

    Returns
    -------
    pandas.DataFrame
        Columns are a ``MultiIndex`` with levels ``('method', 'id')``.
    """
    import pandas as pd
    methods = list(methods) if methods is not None else available_methods()
    frames = []
    for m in methods:
        df = smarts_hits(
            smiles, method=m,
            as_dataframe=True,
            drop_zero_columns=drop_zero_columns,
        )
        df.columns = pd.MultiIndex.from_product(
            [[m], df.columns], names=['method', 'id'],
        )
        frames.append(df)
    return pd.concat(frames, axis=1)
