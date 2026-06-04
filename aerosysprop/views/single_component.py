"""
Single-component properties page.

Predicts per-compound boiling points, vapour pressures, critical properties,
and liquid densities for a loaded set of molecules. This is the original
AeroSysProp prediction flow, now one page of the multipage app.
"""
import os
import tempfile

import pandas as pd
import streamlit as st

from aerosysprop import db
from aerosysprop._ui import (
    ensure_session_state,
    load_compounds_panel,
    results_table,
    sidebar_session_loader,
)
from aerosysprop.batch import (
    METHOD_GROUPS,
    METHOD_META,
    run_all,
)

ensure_session_state()

TABLE_THRESHOLD = 500


# ---------------------------------------------------------------------------
# Sidebar — method selection and settings
# ---------------------------------------------------------------------------

with st.sidebar:
    st.title('Single-component')
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
            dep_note = ' *(needs a boiling point)*' if deps else ''
            checked = st.checkbox(
                f'{label} [{unit}]{dep_note}',
                key=f'method_{m}',
                value=m in ('bp_nannoolal', 'vp_simpol', 'vp_evaporation'),
            )
            if checked:
                selected_methods.append(m)

    # Boiling-point input — sits with the methods (not under generic settings)
    # because it is an input dependency of the VP / CP / density methods above.
    # Selecting more than one source runs each dependent method once per source.
    needs_bp = any(METHOD_META[m][3] for m in selected_methods)
    st.markdown('**Boiling-point input**')
    st.caption(
        'Source(s) feeding the vapour-pressure, critical-property and '
        'density methods that need a boiling point. Pick several to compare '
        'the sensitivity to the boiling-point choice — each dependent result '
        'is then reported once per source.'
    )
    bp_sources = st.multiselect(
        'Boiling-point source(s)',
        options=['bp_nannoolal', 'bp_stein', 'bp_joback'],
        default=['bp_nannoolal'],
        format_func=lambda m: METHOD_META[m][0],
        label_visibility='collapsed',
    )
    if needs_bp and not bp_sources:
        st.warning('Select at least one boiling-point source for the '
                   'boiling-point-dependent methods.')

    st.divider()
    st.subheader('Settings')
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
    sidebar_session_loader()


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

st.header('🧪 Single-component properties')
st.subheader('Input')
mols = load_compounds_panel(state_key='mols', key_prefix='sc_')


# ---------------------------------------------------------------------------
# Run predictions
# ---------------------------------------------------------------------------

st.divider()
st.header('Predictions')

if not selected_methods:
    st.warning('Select at least one method in the sidebar.')
elif not mols:
    st.info('Load compounds above, then click Run.')
else:
    n = len(mols)
    st.write(
        f'Will run **{len(selected_methods)} method(s)** on **{n} compound(s)** '
        f'at **{temperature} K**.'
    )
    if n > TABLE_THRESHOLD:
        st.info(
            f'For sets larger than {TABLE_THRESHOLD} compounds the results '
            'table is hidden to keep the UI responsive. Use Download CSV or '
            'Save session to access full results.',
            icon='ℹ️',
        )
    if st.button('▶ Run', type='primary', use_container_width=True):
        progress_bar = st.progress(0.0)
        status_text = st.empty()

        def _update(done, total):
            progress_bar.progress(done / total)
            status_text.caption(f'Processing {done} / {total} molecules…')

        with st.spinner('Running predictions…'):
            results = run_all(
                mols,
                methods=selected_methods,
                temperature=temperature,
                workers=int(workers),
                progress_callback=_update,
                use_threads=True,
                bp_sources=bp_sources or ['bp_nannoolal'],
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
    valid = [r for r in results if r is not None]
    failed = len(results) - len(valid)

    if failed:
        st.warning(
            f'{failed} molecule(s) could not be processed '
            '(invalid SMILES or method failure).'
        )
    if valid:
        df = pd.DataFrame(valid)
        results_table(df, key_prefix='sc_', table_threshold=TABLE_THRESHOLD,
                      csv_filename='aerosysprop_results.csv')

        st.subheader('Save to session database')
        if st.button('💾 Save to session database', use_container_width=True):
            with tempfile.NamedTemporaryFile(
                delete=False, suffix='.db',
                prefix='aerosysprop_session_',
            ) as tmp:
                session_path = tmp.name
            db.save_results(session_path, valid, temperature)
            st.session_state.session_db_path = session_path
            with open(session_path, 'rb') as f:
                st.download_button(
                    '⬇ Download session (.db)',
                    data=f.read(),
                    file_name='aerosysprop_session.db',
                    mime='application/octet-stream',
                    use_container_width=True,
                )
            st.success('Session saved. Reload it next time via the sidebar.')
