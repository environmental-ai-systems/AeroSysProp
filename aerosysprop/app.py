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

import base64
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


def _data_uri(path, mime):
    b64 = base64.b64encode(path.read_bytes()).decode('ascii')
    return f'data:{mime};base64,{b64}'


# Enlarge the persistent top-left logo (st.logo's built-in sizes top out at
# 'large', which is still fairly small for a wide wordmark).
st.markdown(
    '<style>[data-testid="stLogo"]{height:3.2rem;width:auto;}</style>',
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------------------
# Opening page — play the animated logo once per session, then continue.
# ---------------------------------------------------------------------------
if not st.session_state.get('_intro_played', False):
    if _LOGO_GIF.exists():
        # Embed the GIF as a base64 <img> so the browser plays the animation.
        # (st.image renders only the first frame of an animated GIF.) The
        # element streams to the browser and animates while the script waits,
        # then we flip the flag and rerun into the app.
        st.markdown(
            f'''
            <div style="display:flex;flex-direction:column;align-items:center;
                        justify-content:center;min-height:60vh;">
                <img src="{_data_uri(_LOGO_GIF, "image/gif")}"
                     style="width:min(720px,80vw);height:auto;" />
                <p style="margin-top:1rem;color:#888;font-size:0.9rem;">
                    Loading AeroSysProp…</p>
            </div>
            ''',
            unsafe_allow_html=True,
        )
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
