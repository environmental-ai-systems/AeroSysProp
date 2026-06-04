"""
aerosysprop.db
==============
Session database for persisting compounds and prediction results.

All data lives in a single SQLite file.  Two tables:

- ``compounds``  — SMILES strings with a name/identifier
- ``results``    — long-form predictions (one row per compound × method)

This schema makes it easy to query subsets by method, add new methods to an
existing session without schema changes, and export filtered views to CSV.
"""
import csv
import sqlite3
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_DDL = """
CREATE TABLE IF NOT EXISTS compounds (
    id       INTEGER PRIMARY KEY AUTOINCREMENT,
    name     TEXT,
    smiles   TEXT NOT NULL,
    added_at TEXT DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))
);

CREATE UNIQUE INDEX IF NOT EXISTS idx_compounds_smiles ON compounds(smiles);

CREATE TABLE IF NOT EXISTS results (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    compound_id INTEGER NOT NULL REFERENCES compounds(id) ON DELETE CASCADE,
    method      TEXT NOT NULL,
    temperature REAL NOT NULL,
    value       REAL,
    run_at      TEXT DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))
);

CREATE INDEX IF NOT EXISTS idx_results_compound ON results(compound_id);
CREATE INDEX IF NOT EXISTS idx_results_method   ON results(method);

-- Mixture (AIOMFAC activity-coefficient) calculations. A mixture is a set of
-- neutral components plus optional salts at a given composition / temperature;
-- each species (neutral, cation, anion, or per-salt mean ionic) gets a row in
-- mixture_species with its computed activity coefficient.
CREATE TABLE IF NOT EXISTS mixtures (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    name        TEXT,
    temperature REAL NOT NULL,
    basis       TEXT NOT NULL DEFAULT 'mass',
    salts       TEXT,
    created_at  TEXT DEFAULT (strftime('%Y-%m-%dT%H:%M:%SZ', 'now'))
);

CREATE TABLE IF NOT EXISTS mixture_species (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    mixture_id INTEGER NOT NULL REFERENCES mixtures(id) ON DELETE CASCADE,
    species    TEXT NOT NULL,
    smiles     TEXT,
    kind       TEXT NOT NULL,          -- 'neutral' | 'cation' | 'anion' | 'salt_mean'
    fraction   REAL,                   -- mass fraction (NULL for ions / mean ionic)
    gamma      REAL
);

CREATE INDEX IF NOT EXISTS idx_mixspecies_mixture ON mixture_species(mixture_id);
"""


def _connect(path):
    conn = sqlite3.connect(path)
    conn.execute('PRAGMA foreign_keys = ON')
    conn.executescript(_DDL)
    return conn


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def save_results(path, results, temperature):
    """
    Persist batch results to a session database.

    Compounds are upserted by SMILES (existing rows are not duplicated).
    Result rows are inserted fresh each run so the history is preserved.

    Parameters
    ----------
    path : str
        Path to the ``.db`` session file (created if it does not exist).
    results : list of dict
        The output of ``batch.run_all`` — one dict per molecule.
    temperature : float
        The temperature (K) used for the predictions.
    """
    conn = _connect(path)
    run_at = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    try:
        for row in results:
            if row is None:
                continue
            conn.execute(
                'INSERT OR IGNORE INTO compounds (name, smiles) VALUES (?, ?)',
                (row['name'], row['smiles'])
            )
            (cid,) = conn.execute(
                'SELECT id FROM compounds WHERE smiles = ?', (row['smiles'],)
            ).fetchone()

            for key, value in row.items():
                if key in ('name', 'smiles') or value is None:
                    continue
                conn.execute(
                    'INSERT INTO results (compound_id, method, temperature, value, run_at) '
                    'VALUES (?, ?, ?, ?, ?)',
                    (cid, key, temperature, value, run_at)
                )
        conn.commit()
    finally:
        conn.close()


def load_results(path, methods=None, temperature=None):
    """
    Load results from a session database as a wide-form list of dicts.

    Parameters
    ----------
    path : str
        Path to the ``.db`` session file.
    methods : list of str or None
        Restrict to these output column names (e.g. ``['bp_nannoolal_K']``).
        ``None`` returns all stored methods.
    temperature : float or None
        Restrict to results computed at this temperature.  ``None`` returns
        all temperatures.

    Returns
    -------
    list of dict
        One dict per compound.  Keys: ``name``, ``smiles``, plus one key per
        method with the most recent value for that compound.
    """
    conn = _connect(path)
    try:
        where_clauses = []
        params = []
        if methods:
            placeholders = ','.join('?' * len(methods))
            where_clauses.append('r.method IN (%s)' % placeholders)
            params.extend(methods)
        if temperature is not None:
            where_clauses.append('r.temperature = ?')
            params.append(temperature)

        where = ('WHERE ' + ' AND '.join(where_clauses)) if where_clauses else ''

        # Fetch the most recent value per (compound, method)
        rows = conn.execute(
            """
            SELECT c.id, c.name, c.smiles, r.method, r.value
            FROM compounds c
            JOIN results r ON r.compound_id = c.id
            {where}
            ORDER BY c.id, r.method, r.run_at DESC
            """.format(where=where),
            params
        ).fetchall()

        compounds = {}
        seen = set()
        for cid, name, smiles, method, value in rows:
            if cid not in compounds:
                compounds[cid] = {'name': name, 'smiles': smiles}
            key = (cid, method)
            if key not in seen:
                compounds[cid][method] = value
                seen.add(key)

        return list(compounds.values())
    finally:
        conn.close()


def list_compounds(path):
    """Return all (name, smiles) pairs stored in the session database."""
    conn = _connect(path)
    try:
        return conn.execute(
            'SELECT name, smiles FROM compounds ORDER BY id'
        ).fetchall()
    finally:
        conn.close()


def save_mixture_result(path, name, species_rows, temperature,
                        basis='mass', salts=None):
    """Persist an AIOMFAC mixture activity-coefficient calculation.

    Parameters
    ----------
    path : str
        Path to the ``.db`` session file (created if it does not exist).
    name : str
        A label for this mixture calculation.
    species_rows : list of dict
        One dict per species with keys ``species``, ``smiles`` (optional),
        ``kind`` ('neutral' | 'cation' | 'anion' | 'salt_mean'),
        ``fraction`` (optional), and ``gamma``.
    temperature : float
        Temperature (K) used.
    basis : str
        Composition basis ('mass' or 'mole').
    salts : sequence of str or None
        Salt names included in the mixture.

    Returns
    -------
    int
        The new mixture row id.
    """
    conn = _connect(path)
    try:
        cur = conn.execute(
            'INSERT INTO mixtures (name, temperature, basis, salts) '
            'VALUES (?, ?, ?, ?)',
            (name, temperature, basis,
             ','.join(salts) if salts else ''),
        )
        mixture_id = cur.lastrowid
        for r in species_rows:
            conn.execute(
                'INSERT INTO mixture_species '
                '(mixture_id, species, smiles, kind, fraction, gamma) '
                'VALUES (?, ?, ?, ?, ?, ?)',
                (mixture_id, r['species'], r.get('smiles'), r['kind'],
                 r.get('fraction'), r.get('gamma')),
            )
        conn.commit()
        return mixture_id
    finally:
        conn.close()


def list_mixtures(path):
    """Return stored mixtures as a list of dicts (id, name, temperature,
    basis, salts, created_at), newest first."""
    conn = _connect(path)
    try:
        rows = conn.execute(
            'SELECT id, name, temperature, basis, salts, created_at '
            'FROM mixtures ORDER BY id DESC'
        ).fetchall()
        return [
            {'id': r[0], 'name': r[1], 'temperature': r[2],
             'basis': r[3], 'salts': r[4], 'created_at': r[5]}
            for r in rows
        ]
    finally:
        conn.close()


def load_mixture_species(path, mixture_id):
    """Return the species rows for one stored mixture as a list of dicts."""
    conn = _connect(path)
    try:
        rows = conn.execute(
            'SELECT species, smiles, kind, fraction, gamma '
            'FROM mixture_species WHERE mixture_id = ? ORDER BY id',
            (mixture_id,),
        ).fetchall()
        return [
            {'species': r[0], 'smiles': r[1], 'kind': r[2],
             'fraction': r[3], 'gamma': r[4]}
            for r in rows
        ]
    finally:
        conn.close()


def export_csv(path, output_path, methods=None, temperature=None):
    """
    Export session results to a CSV file.

    Parameters
    ----------
    path : str
        Path to the ``.db`` session file.
    output_path : str
        Destination CSV path.
    methods : list of str or None
        Restrict columns.  ``None`` exports all stored methods.
    temperature : float or None
        Restrict to a specific temperature.
    """
    rows = load_results(path, methods=methods, temperature=temperature)
    if not rows:
        return
    fields = list(rows[0].keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
