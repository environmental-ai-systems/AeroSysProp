"""Homepage / landing view for the AeroSysProp multipage app."""
import streamlit as st

from aerosysprop._ui import sidebar_session_loader

with st.sidebar:
    st.title('AeroSysProp')
    st.caption('Thermophysical property prediction for atmospheric aerosol research.')
    st.divider()
    sidebar_session_loader()

st.title('🧪 AeroSysProp')
st.markdown(
    'Thermophysical property prediction for atmospheric aerosol research, '
    'built on group-contribution methods. Choose a workflow from the sidebar, '
    'or read the summaries below.'
)

st.divider()

col1, col2, col3 = st.columns(3)

with col1:
    st.subheader('🧪 Single-component properties')
    st.markdown(
        """
Predict **per-compound** thermophysical properties for a list of molecules:

- **Boiling points** — Joback & Reid, Stein & Brown, Nannoolal
- **Vapour pressures** — Nannoolal, Myrdal & Yalkowsky, EVAPORATION, SIMPOL
- **Critical properties** — Nannoolal
- **Liquid densities** — Girolami, Schroeder, Le Bas, Tyn & Calus

Load SMILES by typing, file upload, or from a saved session; run in
parallel; export to CSV or a session database.
"""
    )

with col2:
    st.subheader('⚗️ Mixture properties')
    st.markdown(
        """
Predict properties of a **mixture** rather than each compound
independently. Currently available:

- **AIOMFAC activity coefficients** — short-range modified UNIFAC +
  middle-range organic-ion / ion-ion + long-range Debye–Hückel,
  validated against the Fortran AIOMFAC reference. Handles organic +
  water mixtures and electrolytes with common atmospheric salts
  (NaCl, (NH₄)₂SO₄, NH₄NO₃, …), returning per-species γ and mean
  ionic γ± from a composition you specify.

Further mixture-property methods can be added alongside AIOMFAC here.
"""
    )

with col3:
    st.subheader('🔗 External services')
    st.markdown(
        """
Send the compounds you've loaded to **third-party tools** that accept
SMILES input (e.g. the OPERA property suite), to compare against
AeroSysProp predictions. Export the right file format and, where
supported, run the tool from this interface.
"""
    )

st.divider()

st.info(
    'Use the **sidebar** to navigate between workflows. A session database '
    '(.db) loaded in the sidebar is shared across all pages — compounds and '
    'results persist as you move between them.',
    icon='👈',
)
