"""
External services page.

Send loaded compounds to third-party tools that accept SMILES input (e.g. the
OPERA property suite) and, where supported, run them from this interface.
Results can be compared against AeroSysProp predictions.
"""
import os

import streamlit as st
from rdkit import Chem

from aerosysprop._ui import ensure_session_state, sidebar_session_loader
from aerosysprop.services import REGISTRY

ensure_session_state()

with st.sidebar:
    st.title('External services')
    st.divider()
    sidebar_session_loader()

st.header('External services')
st.caption(
    'Send your loaded compounds to third-party tools that accept SMILES input. '
    'Results from these services can be compared with AeroSysProp predictions. '
    'Compounds are taken from whatever you last loaded on the '
    '**Single-component properties** page.'
)

mols = st.session_state.get('mols', [])

if not mols:
    st.info(
        'Load compounds on the **Single-component properties** page first to '
        'enable exports.'
    )
else:
    st.caption(f'**{len(mols)} compound(s)** available to export.')
    service_tabs = st.tabs([s.name for s in REGISTRY])

    for tab, service in zip(service_tabs, REGISTRY):
        with tab:
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
                st.markdown('#### Export for ' + service.name)
                st.markdown(service.export_format_note)

                lines = [
                    f'{smi} {name}'
                    for name, mol in mols
                    for smi in [Chem.MolToSmiles(mol)]
                ]
                smi_content = '\n'.join(lines) + '\n'

                st.download_button(
                    f'⬇ Download {service.export_filename}',
                    data=smi_content.encode(),
                    file_name=service.export_filename,
                    mime='text/plain',
                    use_container_width=True,
                    key=f'dl_{service.id}',
                )

                with st.expander('Installation and usage instructions'):
                    st.markdown(service.install_notes)

                st.markdown('#### Run locally from this UI')
                st.info(
                    'Direct invocation requires verifying the CLI flags against '
                    f'a local {service.name} installation. Once confirmed, set '
                    'the executable path below and the Run button will become '
                    'active.',
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
                    st.warning(
                        'Subprocess invocation not yet implemented — '
                        'export the file above and run manually for now.'
                    )

            elif service.status == 'api':
                st.info('API-based integration — implementation in progress.',
                        icon='🌐')
