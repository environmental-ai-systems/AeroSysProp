"""
AeroSysProp — Streamlit web interface (navigation router).

Launch with:
    streamlit run aerosysprop/app.py

This is a multipage app built with st.navigation. The workflows live under
aerosysprop/views/ and appear in the sidebar:

  • Homepage                       — overview / landing
  • Single-component properties    — per-compound boiling points, vapour
                                     pressures, critical properties, densities
  • Mixture properties             — activity coefficients for a mixture
                                     (currently AIOMFAC SR + MR + LR)
  • External services              — send loaded compounds to third-party tools
"""
import sys as _sys
import pathlib as _pathlib

# Make the aerosysprop package importable when launched via
# `streamlit run aerosysprop/app.py` (Streamlit puts the script's own
# directory on sys.path, not the repo root). Walk up until we find the
# directory that contains the package, and add it to sys.path. (A
# `pip install -e .` makes this redundant, but the bootstrap keeps the app
# runnable straight from a checkout.)
for _parent in _pathlib.Path(__file__).resolve().parents:
    if (_parent / 'aerosysprop' / '__init__.py').exists():
        if str(_parent) not in _sys.path:
            _sys.path.insert(0, str(_parent))
        break

import time

import streamlit as st

from aerosysprop._ui import ensure_session_state

st.set_page_config(
    page_title='AeroSysProp',
    page_icon='🧪',
    layout='wide',
)

ensure_session_state()

# Logo assets bundled with the package.
_ASSETS = _pathlib.Path(__file__).parent / 'assets'
_LOGO_GIF = _ASSETS / 'aerosysprop_logo.gif'
_LOGO_STATIC = _ASSETS / 'aerosysprop_logo_static.png'

# Length of the opening animation; the splash plays it through once on the
# first load of a session, then the app shows the static logo top-left.
_SPLASH_SECONDS = 7.3


# ---------------------------------------------------------------------------
# Opening page — play the animated logo once per session, then continue.
# ---------------------------------------------------------------------------
if not st.session_state.get('_intro_played', False):
    if _LOGO_GIF.exists():
        # Centre the animated logo; it streams to the browser and plays while
        # the script waits, after which we flip the flag and rerun into the app.
        left, mid, right = st.columns([1, 3, 1])
        with mid:
            st.image(str(_LOGO_GIF), width=720)
            st.caption('Loading AeroSysProp…')
        time.sleep(_SPLASH_SECONDS)
    st.session_state['_intro_played'] = True
    st.rerun()

# Persistent static logo, shown top-left above the sidebar navigation on
# every page once the intro has played.
if _LOGO_STATIC.exists():
    st.logo(str(_LOGO_STATIC), size='large')

# Page script paths are given relative to this entry script's directory
# (Streamlit resolves st.Page paths against the main script's location).
pages = [
    st.Page('views/home.py',
            title='Homepage', icon='🏠', default=True),
    st.Page('views/single_component.py',
            title='Single-component properties', icon='🧪'),
    st.Page('views/mixture.py',
            title='Mixture properties', icon='⚗️'),
    st.Page('views/external_services.py',
            title='External services', icon='🔗'),
]

pg = st.navigation(pages)
pg.run()
