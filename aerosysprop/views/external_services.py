"""
External services page.

Send loaded compounds to third-party tools that accept SMILES input (e.g. the
OPERA property suite) and, where supported, run them from this interface.
Results can be compared against AeroSysProp predictions.
"""
import time

import streamlit as st
from rdkit import Chem

from aerosysprop._ui import ensure_session_state, sidebar_session_loader
from aerosysprop.services import REGISTRY

ensure_session_state()


@st.cache_data(ttl=10, show_spinner=False)
def _service_health(url):
    """Probe an external HTTP service's ``/health`` endpoint.

    Cached briefly so Streamlit reruns (which fire on every interaction) don't
    re-hit the network each time; the 'Re-check' button clears the cache.
    Returns a dict: ``reachable``, ``ready`` (opera_available), ``cmd``,
    ``error``.
    """
    import requests
    try:
        r = requests.get(f'{url}/health', timeout=1.5)
        r.raise_for_status()
        info = r.json()
        return {'reachable': True, 'ready': bool(info.get('opera_available')),
                'cmd': info.get('opera_cmd'), 'error': None}
    except Exception as e:
        return {'reachable': False, 'ready': False, 'cmd': None, 'error': str(e)}

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

                # Methodology / interpretation note for QSAR services.
                if service.endpoint_options:
                    st.info(
                        f'{service.name} predictions are **QSAR-model estimates**. '
                        f'Consult the [{service.name} documentation]'
                        f'({service.homepage}) for the model methodologies, '
                        'training data and references, and check each endpoint\'s '
                        'applicability-domain (`AD`, `AD_index`) and confidence '
                        '(`Conf_index`) columns before relying on a value.',
                        icon='📚',
                    )

                # ---- Run via an HTTP service (e.g. the OPERA container) ----
                svc_url = service.service_url()
                if svc_url:
                    import requests  # bundled with streamlit

                    st.markdown(f'#### Run via the {service.name} service')

                    # Automatic availability check (cached; see _service_health).
                    health = _service_health(svc_url)
                    sc1, sc2 = st.columns([5, 1])
                    if health['reachable'] and health['ready']:
                        sc1.success(
                            f'✅ {service.name} service detected and ready at '
                            f'`{svc_url}`.'
                        )
                    elif health['reachable']:
                        sc1.warning(
                            f'⚠️ A service is running at `{svc_url}` but '
                            f'{service.name} is not configured inside it '
                            f'(`opera_available=false`; cmd `{health["cmd"]}`).'
                        )
                    else:
                        sc1.info(
                            f'ℹ️ {service.name} is **optional** and is not '
                            f'running. The app looked for it at `{svc_url}` '
                            f'(set by `{service.service_url_env}`).'
                        )
                    if sc2.button('↻ Re-check', key=f'recheck_{service.id}',
                                  use_container_width=True):
                        _service_health.clear()
                        st.rerun()

                    run_ready = health['reachable'] and health['ready']

                    if not run_ready:
                        st.markdown(
                            f'**Option A — run our containerised {service.name} '
                            '(recommended).** From the repository root:'
                        )
                        st.code(
                            'docker compose --profile opera up --build opera',
                            language='bash',
                        )
                        st.caption(
                            f'Builds a self-contained {service.name} service '
                            '(MATLAB Runtime bundled). On Apple Silicon, Docker + '
                            'Rosetta run the x86 image automatically — nothing to '
                            'configure. Details in `opera/EXTERNAL_TOOLS.md`.'
                        )
                        with st.expander(
                            "Option B — I don't want Docker. What do I provide?"
                        ):
                            st.markdown(
                                f'Manage your own {service.name} install and '
                                'point the app at it:\n\n'
                                f'1. **Install {service.name}** on an x86 OS '
                                '(Intel Mac, Linux, Windows, or a VM) — download '
                                f'from [its releases page]({service.homepage}). '
                                'It is an x86_64 + MATLAB-Runtime tool, so it '
                                'cannot run natively on Apple Silicon.\n'
                                '2. **Then pick how to use it:**\n'
                                '   - *Manual:* download the `.smi` above and run '
                                'it yourself —\n'
                                '     ```\n'
                                '     run_OPERA.sh <mcr_dir> -s '
                                'aerosysprop_for_opera.smi -o predictions.csv\n'
                                '     ```\n'
                                '   - *In-app button:* wrap your install with the '
                                'bundled HTTP service (`opera/server.py`) and set '
                                f'the **`{service.service_url_env}`** environment '
                                'variable to its address before launching the app '
                                f'(currently `{svc_url}`), then **Re-check**.\n\n'
                                'The app talks to the service over the network, so '
                                f'your {service.name} can run on another machine '
                                'or architecture.'
                            )

                    # Reminder: the service is a separate, independent container.
                    if run_ready:
                        st.caption(
                            f'ℹ️ The {service.name} service runs as a **separate '
                            'container** and keeps running after you close the '
                            'app — stop it with `docker compose --profile opera '
                            'down` (or `docker compose stop opera`).'
                        )

                    # Optional endpoint subset — fewer models = faster big runs.
                    selected_codes = []
                    if run_ready and service.endpoint_options:
                        _label_to_code = {lbl: code for lbl, code, _ in
                                          service.endpoint_options}
                        _picked = st.multiselect(
                            'Endpoints to predict (leave empty = all)',
                            options=[lbl for lbl, _, _ in service.endpoint_options],
                            key=f'eps_{service.id}',
                            help='Pick a subset to speed up large runs — OPERA '
                                 'loads fewer models. Empty = all endpoints.',
                        )
                        selected_codes = [_label_to_code[l] for l in _picked]
                        if selected_codes:
                            st.caption('Will request `-e '
                                       + ' '.join(selected_codes) + '`')
                        with st.expander('ⓘ What these endpoints mean'):
                            st.markdown('\n'.join(
                                f'- **{lbl}** — {desc}'
                                for lbl, _, desc in service.endpoint_options))

                    # --- Run as an asynchronous job with live progress ---
                    # Large inputs can take many minutes in OPERA, so we submit a
                    # background job and poll it rather than holding one long
                    # request open (which would time out).
                    job_key = f'{service.id}_job_id'
                    active_job = st.session_state.get(job_key)

                    if st.button(f'▶ Run {service.name}', type='primary',
                                 disabled=not run_ready or bool(active_job),
                                 key=f'runsvc_{service.id}',
                                 use_container_width=True):
                        names = [name for name, _ in mols]
                        smis = [Chem.MolToSmiles(mol) for _, mol in mols]
                        body = {'smiles': smis, 'names': names}
                        if selected_codes:
                            body['endpoints'] = selected_codes
                        try:
                            r = requests.post(f'{svc_url}/jobs', json=body,
                                              timeout=10)
                            r.raise_for_status()
                            st.session_state[job_key] = r.json()['id']
                            st.rerun()
                        except Exception as e:
                            st.error(f'Could not submit {service.name} job: {e}')

                    active_job = st.session_state.get(job_key)
                    if active_job:
                        try:
                            jr = requests.get(
                                f'{svc_url}/jobs/{active_job}',
                                params={'results': 'false'}, timeout=5,
                            )
                            jr.raise_for_status()
                            job = jr.json()
                        except Exception as e:
                            st.warning(f'Lost contact with the job ({e}); it may '
                                       'still be running on the service.')
                            job = None

                        if job is not None:
                            status = job.get('status')
                            n_total = max(int(job.get('n_total') or 1), 1)
                            n_done = int(job.get('n_done') or 0)
                            if status in ('queued', 'running'):
                                frac = min(max(n_done / n_total, 0.0), 1.0)
                                st.progress(frac, text=(
                                    f'{job.get("message", "Running")} — '
                                    f'{n_done}/{n_total} descriptors · '
                                    f'{job.get("elapsed", 0):.0f}s'
                                ))
                                if st.button('✖ Cancel',
                                             key=f'cancel_{service.id}'):
                                    try:
                                        requests.delete(
                                            f'{svc_url}/jobs/{active_job}',
                                            timeout=5)
                                    except Exception:
                                        pass
                                    st.session_state.pop(job_key, None)
                                    st.rerun()
                                time.sleep(2)
                                st.rerun()
                            elif status == 'done':
                                import pandas as pd
                                try:
                                    fr = requests.get(
                                        f'{svc_url}/jobs/{active_job}',
                                        timeout=120)
                                    fr.raise_for_status()
                                    full = fr.json()
                                    df = pd.DataFrame(full.get('records', []))
                                    st.success(
                                        f'{service.name} returned '
                                        f'{full.get("n", len(df))} row(s).')
                                    st.dataframe(df, use_container_width=True,
                                                 hide_index=True)
                                    st.download_button(
                                        f'⬇ Download {service.id}_results.csv',
                                        data=df.to_csv(index=False).encode(),
                                        file_name=f'{service.id}_results.csv',
                                        mime='text/csv',
                                        key=f'dlres_{service.id}',
                                        use_container_width=True)
                                except Exception as e:
                                    st.error(f'Could not fetch results: {e}')
                                if st.button('Clear result',
                                             key=f'clear_{service.id}'):
                                    st.session_state.pop(job_key, None)
                                    st.rerun()
                            else:  # error or cancelled
                                st.error(
                                    f'{service.name} job {status}: '
                                    f'{job.get("error") or job.get("message")}')
                                if st.button('Dismiss',
                                             key=f'dismiss_{service.id}'):
                                    st.session_state.pop(job_key, None)
                                    st.rerun()
                else:
                    st.markdown('#### Run from this UI')
                    st.info(
                        'Direct invocation is not configured for this service. '
                        'Export the file above and run the tool manually, or '
                        'expose it as an HTTP service (see the Docker docs).',
                        icon='⚙️',
                    )

            elif service.status == 'api':
                st.info('API-based integration — implementation in progress.',
                        icon='🌐')
