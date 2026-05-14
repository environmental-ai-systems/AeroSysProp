"""
AeroSysProp — Streamlit web interface.

Launch with:
    streamlit run aerosysprop/app.py
"""
import io
import os
import tempfile

import pandas as pd
import streamlit as st
from rdkit import Chem

from aerosysprop import db
from aerosysprop.services import REGISTRY
from aerosysprop.batch import (
    ALL_METHODS,
    BP_METHODS,
    METHOD_GROUPS,
    METHOD_META,
    _OUTPUT_KEYS,
    run_all,
)
from aerosysprop.io import (
    create_smiles_db,
    read_smiles_db,
    read_smiles_sdf,
    read_smiles_txt,
)

st.set_page_config(
    page_title='AeroSysProp',
    page_icon='🧪',
    layout='wide',
)

TABLE_THRESHOLD = 500   # show full table below this, summary stats above

# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------

for key, default in [
    ('mols', []),
    ('results', []),
    ('session_db_path', None),
    ('opera_exe', ''),
]:
    if key not in st.session_state:
        st.session_state[key] = default


# ---------------------------------------------------------------------------
# Sidebar — method selection and settings
# ---------------------------------------------------------------------------

with st.sidebar:
    st.title('AeroSysProp')
    st.caption('Thermophysical property prediction for atmospheric aerosol research.')
    st.divider()

    st.subheader('Methods')
    selected_methods = []

    for group, group_methods in METHOD_GROUPS.items():
        st.markdown(f'**{group}**')
        col1, col2 = st.columns([1, 1])
        if col1.button('All', key=f'all_{group}', use_container_width=True):
            for m in group_methods:
                st.session_state[f'method_{m}'] = True
        if col2.button('None', key=f'none_{group}', use_container_width=True):
            for m in group_methods:
                st.session_state[f'method_{m}'] = False

        for m in group_methods:
            label, unit, _, deps = METHOD_META[m]
            dep_note = ''
            if deps:
                dep_note = ' *(needs a boiling point)*'
            checked = st.checkbox(
                f'{label} [{unit}]{dep_note}',
                key=f'method_{m}',
                value=m in ('bp_nannoolal', 'vp_simpol', 'vp_evaporation'),
            )
            if checked:
                selected_methods.append(m)

    st.divider()
    st.subheader('Settings')

    bp_source = st.selectbox(
        'Boiling point source for VP / CP / density',
        options=['bp_nannoolal', 'bp_stein', 'bp_joback'],
        format_func=lambda m: METHOD_META[m][0],
        help=(
            'Vapour pressure, critical property, and liquid density methods that '
            'require a boiling point will use this method as their input.'
        ),
    )

    temperature = st.number_input(
        'Temperature (K)',
        min_value=200.0, max_value=400.0,
        value=298.15, step=1.0, format='%.2f',
    )
    workers = st.number_input(
        'Parallel workers',
        min_value=1, max_value=os.cpu_count() or 1,
        value=min(4, os.cpu_count() or 1),
        help='Number of CPU cores to use. More workers = faster for large sets.',
    )

    st.divider()
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
# File format guidance
# ---------------------------------------------------------------------------

def _format_help():
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
molecules are named `mol_1`, `mol_2`, etc. You can also specify an SD property
tag to use as the name (available via the API — see `read_smiles_sdf`).

**Session database**

A `.db` file previously created by clicking **Save session** below.
Loading one in the sidebar reloads all previously stored compounds and results.
""")


# ---------------------------------------------------------------------------
# Input section
# ---------------------------------------------------------------------------

st.header('Input')
_format_help()

tab_smiles, tab_upload, tab_session, tab_pubchem = st.tabs([
    'Enter SMILES', 'Upload file', 'Load from session', '🔍 Search by name (coming soon)'
])

with tab_smiles:
    st.markdown('Enter one or more SMILES strings, one per line. Optionally prefix with a name: `name SMILES`')
    smiles_input = st.text_area(
        'SMILES',
        placeholder='COO\nacetone CC(C)=O\nCC(=O)O',
        height=160,
        label_visibility='collapsed',
    )
    if st.button('Load SMILES', key='load_smiles'):
        mols = []
        errors = []
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
        st.session_state.mols = mols
        if errors:
            st.warning('\n'.join(errors))
        if mols:
            st.success(f'{len(mols)} molecule(s) loaded.')

with tab_upload:
    uploaded = st.file_uploader(
        'Upload a compound file',
        type=['txt', 'csv', 'tsv', 'prop', 'sdf', 'db'],
    )
    col_smiles, col_name = st.columns(2)
    smiles_col_idx = col_smiles.number_input('SMILES column (0-based, text files only)', value=0, min_value=0)
    name_col_idx   = col_name.number_input('Name column (-1 = none, text files only)', value=-1, min_value=-1)

    if uploaded and st.button('Load file', key='load_file'):
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
            st.session_state.mols = mols
            st.success(f'{len(mols)} molecule(s) loaded from {uploaded.name}.')
        except Exception as e:
            st.error(f'Failed to read file: {e}')
        finally:
            os.unlink(tmp_path)

with tab_session:
    if st.session_state.session_db_path:
        compounds = db.list_compounds(st.session_state.session_db_path)
        st.info(f'{len(compounds)} compounds stored in session.')
        if st.button('Reload compounds from session'):
            mols = []
            for name, smi in compounds:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    mols.append((name, mol))
            st.session_state.mols = mols
            st.success(f'{len(mols)} molecules loaded.')
    else:
        st.info('Load a session database in the sidebar to use this tab.')

with tab_pubchem:
    st.info(
        '**Coming soon** — search by compound name, CAS number, or InChIKey '
        'to retrieve SMILES directly from PubChem.',
        icon='🔍',
    )

# Show what's currently loaded
if st.session_state.mols:
    st.caption(f'**{len(st.session_state.mols)} molecule(s) ready to process.**')


# ---------------------------------------------------------------------------
# Run predictions
# ---------------------------------------------------------------------------

st.divider()
st.header('Predictions')

if not selected_methods:
    st.warning('Select at least one method in the sidebar.')
elif not st.session_state.mols:
    st.info('Load compounds above, then click Run.')
else:
    n = len(st.session_state.mols)
    st.write(f'Will run **{len(selected_methods)} method(s)** on **{n} compound(s)** at **{temperature} K**.')

    if n > TABLE_THRESHOLD:
        st.info(
            f'For sets larger than {TABLE_THRESHOLD} compounds the results table is hidden '
            'to keep the UI responsive. Use **Export CSV** or **Save session** to access full results.',
            icon='ℹ️',
        )

    if st.button('▶ Run', type='primary', use_container_width=True):
        progress_bar = st.progress(0.0)
        status_text  = st.empty()

        def _update(done, total):
            progress_bar.progress(done / total)
            status_text.caption(f'Processing {done} / {total} molecules…')

        with st.spinner('Running predictions…'):
            results = run_all(
                st.session_state.mols,
                methods=selected_methods,
                temperature=temperature,
                workers=int(workers),
                progress_callback=_update,
                use_threads=True,
                bp_source=bp_source,
            )

        progress_bar.progress(1.0)
        status_text.caption(f'Done — {len(results)} molecules processed.')
        st.session_state.results = results
        st.success('Predictions complete.')


# ---------------------------------------------------------------------------
# Results display
# ---------------------------------------------------------------------------

if st.session_state.results:
    results = st.session_state.results
    valid   = [r for r in results if r is not None]
    failed  = len(results) - len(valid)

    if failed:
        st.warning(f'{failed} molecule(s) could not be processed (invalid SMILES or method failure).')

    if valid:
        df = pd.DataFrame(valid)

        # --- Summary statistics (always shown) ---
        st.subheader('Summary statistics')
        numeric_cols = df.select_dtypes('number').columns.tolist()
        if numeric_cols:
            st.dataframe(
                df[numeric_cols].describe().round(4),
                use_container_width=True,
            )

        # --- Full table (only for small sets) ---
        if len(valid) <= TABLE_THRESHOLD:
            st.subheader('Results table')
            page_size = 50
            n_pages   = max(1, (len(valid) + page_size - 1) // page_size)
            page      = st.number_input('Page', min_value=1, max_value=n_pages,
                                        value=1, step=1) - 1
            start, end = page * page_size, (page + 1) * page_size
            st.dataframe(df.iloc[start:end], use_container_width=True)
            st.caption(f'Showing rows {start + 1}–{min(end, len(valid))} of {len(valid)}')
        else:
            st.info(f'Table hidden ({len(valid)} rows). Download the full results below.')

        # --- Export ---
        st.subheader('Export')
        col_csv, col_db = st.columns(2)

        csv_buf = io.StringIO()
        df.to_csv(csv_buf, index=False)
        col_csv.download_button(
            '⬇ Download CSV',
            data=csv_buf.getvalue().encode(),
            file_name='aerosysprop_results.csv',
            mime='text/csv',
            use_container_width=True,
        )

        if col_db.button('💾 Save to session database', use_container_width=True):
            with tempfile.NamedTemporaryFile(delete=False, suffix='.db',
                                             prefix='aerosysprop_session_') as tmp:
                session_path = tmp.name
            db.save_results(session_path, valid, temperature)
            st.session_state.session_db_path = session_path

            # Offer the .db file for download
            with open(session_path, 'rb') as f:
                st.download_button(
                    '⬇ Download session (.db)',
                    data=f.read(),
                    file_name='aerosysprop_session.db',
                    mime='application/octet-stream',
                    use_container_width=True,
                )
            st.success('Session saved. You can reload it next time via the sidebar.')


# ---------------------------------------------------------------------------
# External Services
# ---------------------------------------------------------------------------

st.divider()
st.header('External Services')
st.caption(
    'Send your loaded compounds to third-party tools that accept SMILES input. '
    'Results from these services can be compared with AeroSysProp predictions.'
)

if not st.session_state.mols:
    st.info('Load compounds in the Input section above to enable exports.')
else:
    service_tabs = st.tabs([s.name for s in REGISTRY])

    for tab, service in zip(service_tabs, REGISTRY):
        with tab:
            # Header row: name + homepage link
            col_title, col_link = st.columns([3, 1])
            col_title.subheader(service.name)
            col_link.markdown(
                f'[Documentation / Download ↗]({service.homepage})',
                unsafe_allow_html=False,
            )

            st.markdown(f'*{service.tagline}*')
            st.markdown(service.description)
            st.divider()

            if service.status == 'coming_soon':
                st.info('Support for this service is coming soon.', icon='🔧')

            elif service.status == 'file_export':
                # --- Export file ---
                st.markdown('#### Export for ' + service.name)
                st.markdown(service.export_format_note)

                # Build the .smi file content
                lines = [
                    f'{smi} {name}'
                    for name, mol in st.session_state.mols
                    for smi in [Chem.MolToSmiles(mol)]
                ]
                smi_content = '\n'.join(lines) + '\n'

                st.download_button(
                    f'⬇ Download {service.export_filename}',
                    data=smi_content.encode(),
                    file_name=service.export_filename,
                    mime='text/plain',
                    use_container_width=True,
                )

                # --- Installation / usage instructions ---
                with st.expander('Installation and usage instructions'):
                    st.markdown(service.install_notes)

                # --- Direct invocation (subprocess) ---
                st.markdown('#### Run locally from this UI')
                st.info(
                    'Direct invocation requires verifying the CLI flags against a local '
                    f'{service.name} installation. Once confirmed, set the executable path '
                    'below and the Run button will become active.',
                    icon='⚙️',
                )

                exe_key = f'{service.id}_exe'
                if exe_key not in st.session_state:
                    st.session_state[exe_key] = ''

                exe_path = st.text_input(
                    f'{service.name} executable path',
                    value=st.session_state[exe_key],
                    placeholder='/path/to/OPERA2.x_CL',
                    key=f'exe_input_{service.id}',
                )
                st.session_state[exe_key] = exe_path

                run_disabled = not (exe_path and os.path.isfile(exe_path))
                if st.button(
                    f'▶ Run {service.name}',
                    disabled=run_disabled,
                    key=f'run_{service.id}',
                    use_container_width=True,
                ):
                    # TODO: implement subprocess call once CLI flags are verified.
                    # Expected pattern:
                    #   import subprocess, tempfile, pathlib
                    #   with tempfile.TemporaryDirectory() as out_dir:
                    #       smi_tmp = pathlib.Path(out_dir) / 'input.smi'
                    #       smi_tmp.write_text(smi_content)
                    #       argv = service.subprocess_cmd(str(smi_tmp), out_dir)
                    #       subprocess.run(argv, check=True)
                    #       result_csv = next(pathlib.Path(out_dir).glob('*.csv'))
                    #       opera_df = pd.read_csv(result_csv)
                    #   st.dataframe(opera_df, use_container_width=True)
                    st.warning(
                        'Subprocess invocation not yet implemented — '
                        'export the file above and run manually for now.'
                    )

            elif service.status == 'api':
                # Placeholder for future API-based services
                st.info('API-based integration — implementation in progress.', icon='🌐')
