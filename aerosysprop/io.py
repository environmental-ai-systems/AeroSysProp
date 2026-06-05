"""
aerosysprop.io
==============
Utilities for loading SMILES strings into RDKit molecule objects from
common file formats, ready to pass directly to the aerosysprop methods.

Supported sources
-----------------
- Plain text files  (.txt / .prop / any extension)
- SQLite databases  (.db)   — see the schema note below
- SDF files         (.sdf)

All loader functions return a list of (name, mol) tuples where *name* is a
string identifier and *mol* is an RDKit Mol object.  Invalid SMILES are
skipped with a warning rather than raising an exception, so a partial result
is always returned.

SQLite schema
-------------
The expected schema (written by the companion storage tool) is::

    CREATE TABLE compounds (
        id      INTEGER PRIMARY KEY AUTOINCREMENT,
        name    TEXT,
        smiles  TEXT NOT NULL
    );

The table name and column names are configurable via keyword arguments so
that existing databases with different schemas can be used without
modification.
"""
import sqlite3
import warnings

from rdkit import Chem
from rdkit.Chem import SDMolSupplier


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _mol_from_smiles(smiles, name):
    """Return (name, mol) or None if the SMILES is invalid."""
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        warnings.warn('Skipping invalid SMILES for %r: %s' % (name, smiles))
        return None
    return (name, mol)


def _split_line(line, delimiter=None):
    """Split one line into columns. With delimiter=None, auto-detect:
    tab > comma > whitespace (SMILES tokens contain none of these)."""
    if delimiter is None:
        if '\t' in line:
            return [p.strip() for p in line.split('\t')]
        if ',' in line:
            return [p.strip() for p in line.split(',')]
        return line.split()
    return [p.strip() for p in line.split(delimiter)]


def detect_columns(lines, delimiter=None, comment='#', sample=50):
    """Inspect delimited text and decide which column holds the SMILES.

    ``lines`` is any iterable of raw text lines (e.g. a file object or the
    splitlines() of uploaded content). The column whose cells most often
    parse as valid SMILES is taken as the SMILES column; the first other
    column becomes the name column.

    Returns ``(smiles_col, name_col, ncols)`` — both indices 0-based,
    ``name_col`` is ``None`` for single-column input, ``ncols`` is the widest
    row seen (0 if the input had no data lines).
    """
    rows = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(comment):
            continue
        rows.append(_split_line(line, delimiter))
        if len(rows) >= sample:
            break
    if not rows:
        return 0, None, 0
    ncols = max(len(r) for r in rows)
    if ncols == 1:
        return 0, None, 1
    # Probing each cell as a SMILES is expected to fail for name columns;
    # silence RDKit's parse-error logging for the duration so it doesn't spam
    # the console, then restore it.
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.error')
    try:
        valid = [0] * ncols
        for r in rows:
            for c in range(min(len(r), ncols)):
                if r[c] and Chem.MolFromSmiles(r[c]) is not None:
                    valid[c] += 1
    finally:
        RDLogger.EnableLog('rdApp.error')
    smiles_col = max(range(ncols), key=lambda c: valid[c])
    name_col = next((c for c in range(ncols) if c != smiles_col), None)
    return smiles_col, name_col, ncols


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def read_smiles_txt(path, smiles_col=1, name_col=0, delimiter=None, comment='#'):
    """
    Load SMILES from a plain-text file.

    Each non-empty, non-comment line is split on *delimiter* (any whitespace
    by default).  *smiles_col* and *name_col* are zero-based column indices.

    If the line has only one column it is assumed to be the SMILES itself and
    *name_col* is ignored (the SMILES string is used as the name).

    Parameters
    ----------
    path : str
        Path to the text file.
    smiles_col : int
        Column index (0-based) containing the SMILES string.  Default 1,
        matching the tab-separated ``index\\tSMILES`` format used in the
        existing AeroSysProp example files.
    name_col : int
        Column index (0-based) for the compound name/ID.  Default 0.
    delimiter : str or None
        Column separator. ``None`` (default) auto-detects per line: a tab if
        present, else a comma if present, else any whitespace. (SMILES tokens
        never contain commas, tabs or spaces, so this is unambiguous.) Pass an
        explicit string to force a specific separator.
    comment : str
        Lines whose first non-whitespace character matches this string are
        ignored.

    Returns
    -------
    list of (str, rdkit.Chem.Mol)
        (name, mol) pairs for every valid SMILES in the file.

    Examples
    --------
    >>> # File with one SMILES per line
    >>> mols = read_smiles_txt('compounds.txt', smiles_col=0)

    >>> # Tab-separated index\\tSMILES (AeroSysProp default format)
    >>> mols = read_smiles_txt('compounds.prop')
    """
    results = []
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith(comment):
                continue
            parts = _split_line(line, delimiter)
            if len(parts) == 1:
                smiles = parts[0]
                name = smiles
            elif max(smiles_col, name_col) < len(parts):
                smiles = parts[smiles_col]
                name = parts[name_col]
            else:
                warnings.warn('Skipping line %d in %s: not enough columns' % (lineno, path))
                continue
            entry = _mol_from_smiles(smiles, name)
            if entry:
                results.append(entry)
    return results


def read_smiles_db(path, table='compounds', smiles_col='smiles', name_col='name'):
    """
    Load SMILES from a SQLite database file.

    Parameters
    ----------
    path : str
        Path to the ``.db`` SQLite file.
    table : str
        Table name to query.  Default ``'compounds'``.
    smiles_col : str
        Column containing the SMILES strings.  Default ``'smiles'``.
    name_col : str
        Column used as the compound name/ID.  If the column does not exist
        in the table the SMILES string is used as the name instead.

    Returns
    -------
    list of (str, rdkit.Chem.Mol)
        (name, mol) pairs for every valid SMILES in the database.

    Examples
    --------
    >>> mols = read_smiles_db('compounds.db')
    >>> mols = read_smiles_db('mydb.db', table='molecules', smiles_col='smi', name_col='compound_id')
    """
    conn = sqlite3.connect(path)
    try:
        cursor = conn.execute('SELECT * FROM %s LIMIT 0' % table)
        col_names = [d[0] for d in cursor.description]

        has_name = name_col in col_names
        if has_name:
            query = 'SELECT %s, %s FROM %s' % (name_col, smiles_col, table)
        else:
            query = 'SELECT %s FROM %s' % (smiles_col, table)

        results = []
        for row in conn.execute(query):
            if has_name:
                name, smiles = str(row[0]), row[1]
            else:
                smiles = row[0]
                name = smiles
            entry = _mol_from_smiles(smiles, name)
            if entry:
                results.append(entry)
    finally:
        conn.close()
    return results


def read_smiles_sdf(path, name_property=None):
    """
    Load molecules from an SDF file.

    The compound name is taken from (in order of preference):

    1. The SD property named by *name_property* (if supplied).
    2. The molecule's title line (the first line of each record, accessible
       via ``mol.GetProp('_Name')``).
    3. A sequential index (``mol_1``, ``mol_2``, …) as a fallback.

    Parameters
    ----------
    path : str
        Path to the ``.sdf`` file.
    name_property : str or None
        Name of an SD tag to use as the compound identifier, e.g.
        ``'CAS'`` or ``'compound_id'``.

    Returns
    -------
    list of (str, rdkit.Chem.Mol)
        (name, mol) pairs for every valid structure in the file.

    Examples
    --------
    >>> mols = read_smiles_sdf('compounds.sdf')
    >>> mols = read_smiles_sdf('pubchem_download.sdf', name_property='PUBCHEM_COMPOUND_CID')
    """
    supplier = SDMolSupplier(path, removeHs=True, sanitize=True)
    results = []
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            warnings.warn('Skipping unreadable structure %d in %s' % (idx, path))
            continue
        if name_property and mol.HasProp(name_property):
            name = mol.GetProp(name_property)
        elif mol.HasProp('_Name') and mol.GetProp('_Name').strip():
            name = mol.GetProp('_Name').strip()
        else:
            name = 'mol_%d' % idx
        results.append((name, mol))
    return results


def create_smiles_db(path, compounds, table='compounds'):
    """
    Create (or extend) a SQLite database with a compounds table.

    This is a convenience function for the companion storage tool.  Calling
    it on an existing database appends rows without duplicating the table.

    Parameters
    ----------
    path : str
        Path to the ``.db`` file to create or update.
    compounds : iterable of (str, str)
        Sequence of ``(name, smiles)`` pairs to insert.
    table : str
        Table name.  Default ``'compounds'``.

    Examples
    --------
    >>> create_smiles_db('mydata.db', [('ethanol', 'CCO'), ('acetone', 'CC(C)=O')])
    """
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            'CREATE TABLE IF NOT EXISTS %s ('
            '  id     INTEGER PRIMARY KEY AUTOINCREMENT,'
            '  name   TEXT,'
            '  smiles TEXT NOT NULL'
            ')' % table
        )
        conn.executemany(
            'INSERT INTO %s (name, smiles) VALUES (?, ?)' % table,
            compounds
        )
        conn.commit()
    finally:
        conn.close()
