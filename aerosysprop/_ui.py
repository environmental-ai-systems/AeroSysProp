"""
aerosysprop._ui
===============
Shared Streamlit UI helpers used by the multipage app (``app.py`` and the
modules under ``aerosysprop/pages/``). Keeping these in one place avoids
duplicating the compound-loading, results-table, and session-database
widgets across the single-component and mixture pages.

These functions assume an active Streamlit script context (they call
``st.*``). They are not part of the headless prediction API.
"""
import io
import os
import tempfile

import pandas as pd
import streamlit as st
from rdkit import Chem

from aerosysprop import db
from aerosysprop.io import (
    read_smiles_db,
    read_smiles_sdf,
    read_smiles_txt,
)


# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------

def ensure_session_state(extra_defaults=None):
    """Initialise common session-state keys shared across pages.

    Pages may pass ``extra_defaults`` (a dict) for page-specific keys.
    """
    defaults = {
        'mols': [],                 # single-component compound set
        'results': [],              # single-component results
        'mix_components': [],       # mixture neutral components
        'mix_results': None,        # last AIOMFACResult-derived dict
        'session_db_path': None,
    }
    if extra_defaults:
        defaults.update(extra_defaults)
    for key, default in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = default


# ---------------------------------------------------------------------------
# File-format guidance
# ---------------------------------------------------------------------------

# A small atmospherically-relevant set used to generate the downloadable
# example files offered in the format-help expander.
_EXAMPLE_COMPOUNDS = [
    ('water',         'O'),
    ('ethanol',       'CCO'),
    ('acetic_acid',   'CC(=O)O'),
    ('glycerol',      'OCC(O)CO'),
    ('levoglucosan',  'OC[C@H]1OC2OC1C(O)C2O'),
    ('pinonic_acid',  'CC(=O)C1CCC(C1(C)C)CC(=O)O'),
]


def _example_smiles_txt():
    """One SMILES per line — loads with the default upload settings."""
    lines = ['# One SMILES per line. Upload with "SMILES column" = 0 (default).']
    lines += [smi for _, smi in _EXAMPLE_COMPOUNDS]
    return ('\n'.join(lines) + '\n').encode()


def _example_named_txt():
    """name <TAB> SMILES — set "SMILES column" = 1, "Name column" = 0."""
    lines = ['# name <TAB> SMILES. Upload with "SMILES column" = 1, "Name column" = 0.']
    lines += [f'{name}\t{smi}' for name, smi in _EXAMPLE_COMPOUNDS]
    return ('\n'.join(lines) + '\n').encode()


def _example_sdf():
    """A small SDF (with 2D coordinates and molecule names)."""
    from rdkit.Chem import AllChem
    blocks = []
    for name, smi in _EXAMPLE_COMPOUNDS:
        mol = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(mol)
        mol.SetProp('_Name', name)
        blocks.append(Chem.MolToMolBlock(mol) + '$$$$\n')
    return ''.join(blocks).encode()


def format_help():
    with st.expander('📋 Supported file formats', expanded=False):
        st.markdown("""
| Format | Extensions | Structure |
|---|---|---|
| **Plain text** | `.txt` `.csv` `.tsv` `.prop` | One SMILES per line, or two columns: `name\\tSMILES` |
| **SDF** | `.sdf` | Standard MDL Structure-Data File — exported from ChemDraw, MarvinSketch, PubChem, etc. |
| **Session database** | `.db` | SQLite file previously saved from AeroSysProp |

**Plain text examples**

Single-column (SMILES only):
```
COO
CC(C)=O
CC(=O)O
```

Two-column, tab-separated (name then SMILES):
```
methyl_hydroperoxide    COO
acetone                 CC(C)=O
acetic_acid             CC(=O)O
```

Lines beginning with `#` are treated as comments and skipped.

**SDF files**

The compound name is read from the molecule title line. If no title is present,
molecules are named `mol_1`, `mol_2`, etc.

**Session database**

A `.db` file previously created by clicking **Save session**. Loading one in
the sidebar reloads all previously stored compounds and results.
""")

        st.markdown('**Download example files** to try the loaders:')
        c1, c2, c3 = st.columns(3)
        c1.download_button(
            'SMILES list (.txt)',
            data=_example_smiles_txt(),
            file_name='example_smiles.txt', mime='text/plain',
            width='stretch', key='ex_smiles_txt',
        )
        c2.download_button(
            'Named list (.txt)',
            data=_example_named_txt(),
            file_name='example_named.txt', mime='text/plain',
            width='stretch', key='ex_named_txt',
        )
        c3.download_button(
            'SDF (.sdf)',
            data=_example_sdf(),
            file_name='example_compounds.sdf',
            mime='chemical/x-mdl-sdfile',
            width='stretch', key='ex_sdf',
        )
        st.caption(
            'The **SMILES list** loads with the default upload settings. The '
            '**Named list** is tab-separated `name⇥SMILES` — set "SMILES '
            'column" = 1 and "Name column" = 0 when uploading it. The **SDF** '
            'loads directly; names come from each molecule title.'
        )


# ---------------------------------------------------------------------------
# Sidebar: session database loader (shared)
# ---------------------------------------------------------------------------

def sidebar_session_loader():
    """Render the 'Load existing session (.db)' uploader in the sidebar and
    store the resulting path in session state. Safe to call on every page."""
    st.subheader('Session database')
    session_db_file = st.file_uploader(
        'Load existing session (.db)',
        type=['db'],
        key='session_upload',
    )
    if session_db_file:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.db')
        tmp.write(session_db_file.read())
        tmp.flush()
        st.session_state.session_db_path = tmp.name
        st.success(f'Session loaded: {session_db_file.name}')


# ---------------------------------------------------------------------------
# Compound-loading input panel (shared by single-component and mixture pages)
# ---------------------------------------------------------------------------

def load_compounds_panel(state_key='mols', key_prefix='', help_text=None):
    """Render the compound-loading panel (SMILES text / file upload / session
    load) and store the loaded ``(name, mol)`` list in
    ``st.session_state[state_key]``.

    Parameters
    ----------
    state_key : str
        The session-state key under which to store the loaded molecules. Use
        a distinct key per page (e.g. 'mols' vs 'mix_components') so the
        single-component and mixture compound sets don't collide.
    key_prefix : str
        Prefix for widget keys, to keep them unique across pages.
    help_text : str, optional
        Extra markdown shown above the SMILES box.
    """
    format_help()

    tab_smiles, tab_upload, tab_session = st.tabs([
        'Enter SMILES', 'Upload file', 'Load from session',
    ])

    with tab_smiles:
        if help_text:
            st.markdown(help_text)
        else:
            st.markdown(
                'Enter one or more SMILES strings, one per line. Optionally '
                'prefix with a name: `name SMILES`'
            )
        smiles_input = st.text_area(
            'SMILES',
            placeholder='COO\nacetone CC(C)=O\nCC(=O)O',
            height=160,
            label_visibility='collapsed',
            key=f'{key_prefix}smiles_text',
        )
        if st.button('Load SMILES', key=f'{key_prefix}load_smiles'):
            mols, errors = [], []
            for i, line in enumerate(smiles_input.strip().splitlines(), start=1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) == 1:
                    smi, name = parts[0], parts[0]
                else:
                    name, smi = parts[0], parts[1]
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    errors.append(f'Line {i}: invalid SMILES — {smi!r}')
                else:
                    mols.append((name, mol))
            st.session_state[state_key] = mols
            if errors:
                st.warning('\n'.join(errors))
            if mols:
                st.success(f'{len(mols)} molecule(s) loaded.')

    with tab_upload:
        uploaded = st.file_uploader(
            'Upload a compound file',
            type=['txt', 'csv', 'tsv', 'prop', 'sdf', 'db'],
            key=f'{key_prefix}file_upload',
        )
        col_smiles, col_name = st.columns(2)
        smiles_col_idx = col_smiles.number_input(
            'SMILES column (0-based, text files only)', value=0, min_value=0,
            key=f'{key_prefix}smiles_col',
        )
        name_col_idx = col_name.number_input(
            'Name column (-1 = none, text files only)', value=-1, min_value=-1,
            key=f'{key_prefix}name_col',
        )
        if uploaded and st.button('Load file', key=f'{key_prefix}load_file'):
            suffix = '.' + uploaded.name.rsplit('.', 1)[-1].lower()
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(uploaded.read())
                tmp_path = tmp.name
            try:
                if suffix == '.sdf':
                    mols = read_smiles_sdf(tmp_path)
                elif suffix == '.db':
                    mols = read_smiles_db(tmp_path)
                else:
                    kwargs = {'smiles_col': int(smiles_col_idx)}
                    if name_col_idx >= 0:
                        kwargs['name_col'] = int(name_col_idx)
                    mols = read_smiles_txt(tmp_path, **kwargs)
                st.session_state[state_key] = mols
                st.success(f'{len(mols)} molecule(s) loaded from {uploaded.name}.')
            except Exception as e:
                st.error(f'Failed to read file: {e}')
            finally:
                os.unlink(tmp_path)

    with tab_session:
        if st.session_state.session_db_path:
            compounds = db.list_compounds(st.session_state.session_db_path)
            st.info(f'{len(compounds)} compounds stored in session.')
            if st.button('Reload compounds from session',
                         key=f'{key_prefix}reload_session'):
                mols = []
                for name, smi in compounds:
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        mols.append((name, mol))
                st.session_state[state_key] = mols
                st.success(f'{len(mols)} molecules loaded.')
        else:
            st.info('Load a session database in the sidebar to use this tab.')

    mols = st.session_state.get(state_key, [])
    if mols:
        st.caption(f'**{len(mols)} molecule(s) ready.**')
    return mols


# ---------------------------------------------------------------------------
# Results table + CSV export (shared)
# ---------------------------------------------------------------------------

def results_table(df, key_prefix='', table_threshold=500,
                  csv_filename='aerosysprop_results.csv'):
    """Render summary statistics, a paginated table (for small sets), and a
    CSV download button for a results DataFrame."""
    st.subheader('Summary statistics')
    numeric_cols = df.select_dtypes('number').columns.tolist()
    if numeric_cols:
        st.dataframe(df[numeric_cols].describe().round(4),
                     use_container_width=True)

    if len(df) <= table_threshold:
        st.subheader('Results table')
        page_size = 50
        n_pages = max(1, (len(df) + page_size - 1) // page_size)
        page = st.number_input('Page', min_value=1, max_value=n_pages,
                               value=1, step=1, key=f'{key_prefix}page') - 1
        start, end = page * page_size, (page + 1) * page_size
        st.dataframe(df.iloc[start:end], use_container_width=True)
        st.caption(f'Showing rows {start + 1}–{min(end, len(df))} of {len(df)}')
    else:
        st.info(f'Table hidden ({len(df)} rows). Download the full results below.')

    csv_buf = io.StringIO()
    df.to_csv(csv_buf, index=False)
    st.download_button(
        '⬇ Download CSV',
        data=csv_buf.getvalue().encode(),
        file_name=csv_filename,
        mime='text/csv',
        use_container_width=True,
        key=f'{key_prefix}csv_download',
    )
