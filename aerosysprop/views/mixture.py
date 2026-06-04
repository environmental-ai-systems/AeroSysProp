"""
Mixture activity coefficients page.

Computes AIOMFAC activity coefficients (short-range modified UNIFAC +
middle-range organic-ion / ion-ion + long-range Debye-Hückel) for a mixture
of neutral organic components, water, and optional salts at a user-specified
composition and temperature.

Unlike the single-component page (which predicts each compound independently),
this page treats the loaded compounds as ONE coupled mixture: every species'
activity coefficient depends on the whole composition.
"""
import tempfile

import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem

from aerosysprop import data, db
from aerosysprop._ui import (
    ensure_session_state,
    load_compounds_panel,
    sidebar_session_loader,
)
from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles

ensure_session_state()

# Salt catalog (names sorted for the multiselect)
_SALT_NAMES = sorted(set(data.AIOMFAC_SALT_NAME.values()))


with st.sidebar:
    st.title('Mixture properties')
    st.divider()
    st.subheader('Method')
    method = st.selectbox(
        'Mixture-property method',
        options=['AIOMFAC activity coefficients'],
        help=(
            'AIOMFAC: short-range modified UNIFAC + middle-range organic-ion / '
            'ion-ion + long-range Debye–Hückel. Further mixture-property '
            'methods can be added here.'
        ),
    )
    st.divider()
    st.subheader('Settings')
    temperature = st.number_input(
        'Temperature (K)',
        min_value=200.0, max_value=400.0,
        value=298.15, step=1.0, format='%.2f',
        key='mix_temperature',
    )
    st.divider()
    sidebar_session_loader()


st.header('Mixture properties')
st.markdown(
    'Predict properties of a **mixture** rather than each compound '
    'independently. The active method is **AIOMFAC activity coefficients** '
    '(short-range modified UNIFAC + middle-range + long-range Debye–Hückel): '
    'every species\' activity coefficient depends on the whole composition.'
)

# ---------------------------------------------------------------------------
# Step 1 — neutral components
# ---------------------------------------------------------------------------

st.subheader('1. Neutral components')
st.markdown(
    'Load the neutral organic components (and water) as SMILES. **Water is '
    '`O`** and, by AIOMFAC convention, is normally the first component. You '
    'can type SMILES, upload a file, or reload from a saved session — the '
    'same loaders as the single-component page.'
)

components = load_compounds_panel(
    state_key='mix_components',
    key_prefix='mix_',
    help_text=(
        'Enter one SMILES per line (optionally `name SMILES`). Include water '
        'as `O` for aqueous mixtures. Example:\n\n```\nwater O\nethanol CCO\n'
        'glycerol OCC(O)CO\n```'
    ),
)

# ---------------------------------------------------------------------------
# Step 2 — salts
# ---------------------------------------------------------------------------

st.subheader('2. Salts (optional)')
salts = st.multiselect(
    'Select inorganic salts from the AIOMFAC catalog',
    options=_SALT_NAMES,
    help=(
        'Leave empty for a non-electrolyte (organic + water) mixture. Salts '
        'dissociate into their ions; activity coefficients are returned per '
        'ion plus a mean ionic value per salt.'
    ),
)

# ---------------------------------------------------------------------------
# Step 3 — composition
# ---------------------------------------------------------------------------

st.subheader('3. Composition (mass fractions)')

if not components:
    st.info('Load neutral components above to set the composition.')
else:
    component_names = [name for name, _ in components]
    all_species = component_names + list(salts)

    # Build / re-use an editable composition table
    default_frac = round(1.0 / max(len(all_species), 1), 4)
    comp_df = pd.DataFrame({
        'component': all_species,
        'mass_fraction': [default_frac] * len(all_species),
    })
    edited = st.data_editor(
        comp_df,
        key='mix_composition_editor',
        use_container_width=True,
        hide_index=True,
        column_config={
            'component': st.column_config.TextColumn('Component', disabled=True),
            'mass_fraction': st.column_config.NumberColumn(
                'Mass fraction', min_value=0.0, max_value=1.0,
                step=0.01, format='%.4f',
            ),
        },
    )
    total = float(edited['mass_fraction'].sum())
    if abs(total - 1.0) > 1e-3:
        st.warning(
            f'Mass fractions sum to {total:.4f}, not 1.0. They will be used '
            'as-is; tick the box below to auto-normalise.'
        )
    normalise = st.checkbox('Auto-normalise mass fractions to sum to 1',
                            value=True)

    # ---------------------------------------------------------------------
    # Step 4 — run
    # ---------------------------------------------------------------------
    st.subheader('4. Run')
    if st.button('▶ Compute activity coefficients', type='primary',
                 use_container_width=True):
        smiles = [Chem.MolToSmiles(mol) for _, mol in components]
        fracs = edited['mass_fraction'].to_numpy(dtype=float)
        if normalise and fracs.sum() > 0:
            fracs = fracs / fracs.sum()

        try:
            result = aiomfac_from_smiles(
                smiles=smiles,
                salts=list(salts),
                mass_fractions=list(fracs),
                T_K=float(temperature),
                strict='warn',
            )
        except Exception as e:
            st.error(f'AIOMFAC calculation failed: {e}')
            result = None

        if result is not None:
            # Assemble the per-species result rows
            rows = []
            for i, (name, smi) in enumerate(zip(component_names, smiles)):
                rows.append({
                    'species': name, 'smiles': smi, 'kind': 'neutral',
                    'fraction': float(fracs[i]),
                    'gamma': float(result.gamma_neutral[i]),
                })
            for cid, g in zip(result.cation_ids, result.gamma_cation):
                rows.append({
                    'species': f'cation {cid}', 'smiles': None,
                    'kind': 'cation', 'fraction': None, 'gamma': float(g),
                })
            for aid, g in zip(result.anion_ids, result.gamma_anion):
                rows.append({
                    'species': f'anion {aid}', 'smiles': None,
                    'kind': 'anion', 'fraction': None, 'gamma': float(g),
                })
            for sname, g in zip(result.salt_names, result.gamma_mean):
                rows.append({
                    'species': f'{sname} (mean ionic)', 'smiles': None,
                    'kind': 'salt_mean', 'fraction': None, 'gamma': float(g),
                })

            st.session_state.mix_results = {
                'rows': rows,
                'temperature': float(temperature),
                'salts': list(salts),
                'diagnostics': result.diagnostics,
            }
            st.success('Activity coefficients computed.')


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

if st.session_state.get('mix_results'):
    res = st.session_state.mix_results
    rows = res['rows']
    df = pd.DataFrame(rows)

    st.divider()
    st.header('Results')
    st.caption(
        f'AIOMFAC activity coefficients at {res["temperature"]} K'
        + (f' with salts: {", ".join(res["salts"])}' if res['salts'] else
           ' (non-electrolyte mixture)')
    )

    # Neutral, ion, and mean-ionic sub-tables
    neutrals = df[df['kind'] == 'neutral']
    if not neutrals.empty:
        st.subheader('Neutral components (mole-fraction basis)')
        st.dataframe(
            neutrals[['species', 'fraction', 'gamma']].rename(
                columns={'fraction': 'mass_fraction', 'gamma': 'gamma_x'}
            ),
            use_container_width=True, hide_index=True,
        )

    ions = df[df['kind'].isin(['cation', 'anion'])]
    if not ions.empty:
        st.subheader('Ions (molality basis)')
        st.dataframe(
            ions[['species', 'gamma']].rename(columns={'gamma': 'gamma_m'}),
            use_container_width=True, hide_index=True,
        )

    means = df[df['kind'] == 'salt_mean']
    if not means.empty:
        st.subheader('Mean ionic activity coefficients')
        st.dataframe(
            means[['species', 'gamma']].rename(columns={'gamma': 'gamma_±'}),
            use_container_width=True, hide_index=True,
        )

    # Diagnostics
    report = res.get('diagnostics')
    if report is not None and len(report) > 0:
        sev = 'error' if report.has_errors else 'warning'
        with st.expander(
            f'⚠️ {len(report)} input diagnostic(s) — click to view',
            expanded=report.has_errors,
        ):
            for issue in report:
                if issue.severity == 'error':
                    st.error(f'**{issue.code}** — {issue.message}')
                else:
                    st.warning(f'**{issue.code}** — {issue.message}')

    # Export
    st.subheader('Export')
    col_csv, col_db = st.columns(2)
    csv = df.to_csv(index=False).encode()
    col_csv.download_button(
        '⬇ Download CSV', data=csv,
        file_name='aiomfac_mixture.csv', mime='text/csv',
        use_container_width=True,
    )

    mix_name = st.text_input('Mixture name (for the session database)',
                             value='mixture_1')
    if col_db.button('💾 Save to session database', use_container_width=True):
        path = st.session_state.session_db_path
        if path is None:
            with tempfile.NamedTemporaryFile(
                delete=False, suffix='.db',
                prefix='aerosysprop_session_',
            ) as tmp:
                path = tmp.name
            st.session_state.session_db_path = path
        db.save_mixture_result(
            path, name=mix_name, species_rows=rows,
            temperature=res['temperature'], basis='mass', salts=res['salts'],
        )
        with open(path, 'rb') as f:
            st.download_button(
                '⬇ Download session (.db)',
                data=f.read(),
                file_name='aerosysprop_session.db',
                mime='application/octet-stream',
                use_container_width=True,
            )
        st.success(f'Mixture "{mix_name}" saved to the session database.')
