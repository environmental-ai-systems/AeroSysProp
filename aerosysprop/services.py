"""
aerosysprop.services
====================
Registry of external property-prediction services that can consume the
same SMILES inputs as AeroSysProp.

Adding a new service
--------------------
Append an ``ExternalService`` entry to ``REGISTRY``.  The UI in ``app.py``
reads this list and renders a tab per service automatically — no UI code
changes needed.

Service status values
---------------------
``'file_export'``
    The service is a locally-installed CLI tool.  AeroSysProp generates a
    correctly-formatted input file for the user to pass to the tool.
    Once the CLI interface is verified, add a ``subprocess_cmd`` function
    to enable direct invocation from the UI.

``'api'``
    The service exposes a public REST API.  AeroSysProp calls it directly
    and displays results alongside its own predictions.

``'coming_soon'``
    Placeholder — shown in the UI but not yet functional.
"""
import os
from dataclasses import dataclass, field
from typing import Callable, List, Optional


@dataclass
class ExternalService:
    id: str
    name: str
    tagline: str
    description: str
    homepage: str
    status: str                        # 'file_export' | 'api' | 'coming_soon'
    export_filename: str = ''          # suggested download filename
    export_format_note: str = ''       # shown to user above the download button
    install_notes: str = ''            # shown when status == 'file_export'
    # Optional HTTP service (e.g. an OPERA container exposing /predict). When
    # ``service_url_env`` is set the UI offers a "Run via service" panel that
    # POSTs SMILES to ``<url>/predict``. The URL is read from the named
    # environment variable, falling back to ``default_service_url`` — so the
    # same code works whether the app runs natively (localhost) or in a
    # container on the compose network (service name).
    service_url_env: str = ''
    default_service_url: str = ''
    # Optional endpoint catalogue for services that can compute a subset of
    # outputs. Each item is ``(label, code, description)``; ``code`` is passed
    # through to the service (for OPERA, mapped to its ``-e`` flag). Empty
    # selection => all.
    endpoint_options: List[tuple] = field(default_factory=list)
    # Filled in once the CLI interface has been verified against a real install:
    # subprocess_cmd: Optional[Callable[[str, str], List[str]]] = None
    # Arguments: (input_smi_path, output_dir) -> argv list

    def service_url(self) -> str:
        """Resolved base URL of the backing HTTP service, or '' if none."""
        if not self.service_url_env:
            return ''
        return os.environ.get(self.service_url_env, self.default_service_url)


# OPERA endpoint catalogue: ``(label, OPERA -e code, one-line description)``.
# Codes are case-insensitive (see OPERA's help.txt). Selecting a subset speeds up
# large runs because OPERA loads fewer models. Individual endpoints first, then
# convenience groups.
_OPERA_ENDPOINTS = [
    ('Boiling point (BP)', 'BP', 'Normal boiling point (°C).'),
    ('Melting point (MP)', 'MP', 'Melting point (°C).'),
    ('Vapour pressure (logVP)', 'VP', 'log₁₀ vapour pressure (log mmHg, 25 °C).'),
    ('Water solubility (logWS)', 'WS', 'log₁₀ water solubility (log mol/L, 25 °C).'),
    ("Henry's law (logHL)", 'HL', "log₁₀ Henry's law constant (log atm·m³/mol)."),
    ('Octanol–water (LogP)', 'logP', 'Octanol–water partition coefficient (log).'),
    ('LogD (pH 5.5/7.4)', 'logD', 'Octanol–water distribution coeff. at pH 5.5 and 7.4.'),
    ('pKa', 'pKa', 'Strongest acidic / basic dissociation constants.'),
    ('Octanol–air (LogKOA)', 'KOA', 'Octanol–air partition coefficient (log).'),
    ('Atmospheric OH rate (logOH)', 'AOH', 'OH-radical reaction rate (log cm³/molecule·s).'),
    ('Bioconcentration (logBCF)', 'BCF', 'Fish bioconcentration factor (log L/kg).'),
    ('Soil adsorption (logKoc)', 'Koc', 'Organic-carbon adsorption coefficient (log L/kg).'),
    ('Biodegradation half-life', 'BioDeg', 'Biotransformation half-life (log days).'),
    ('Ready biodegradability', 'RB', 'Ready-biodegradability class (yes/no).'),
    ('HPLC retention time (RT)', 'RT', 'Reversed-phase HPLC retention time (min).'),
    ('▸ Group: all physicochemical', 'PC', 'Every physicochemical endpoint above in one pass.'),
    ('▸ Group: all environmental fate', 'EF', 'AOH, BCF, biodegradation, Koc and KM together.'),
    ('▸ Group: all toxicity', 'Tox', 'Estrogen/androgen receptor (CERAPP/CoMPARA) + acute oral toxicity (CATMoS).'),
    ('▸ Group: all ADME', 'ADME', 'Fraction unbound, intrinsic clearance and Caco-2 permeability.'),
]


REGISTRY: List[ExternalService] = [

    ExternalService(
        id='opera',
        name='OPERA',
        tagline='QSAR models for physicochemical and environmental fate properties',
        description=(
            'OPERA (OPEn (Q)SAR Application) is an open-source tool developed by '
            'NIEHS and the US EPA that provides predictions for 30+ endpoints including '
            'vapour pressure, boiling point, LogP, water solubility, biodegradation, '
            'and toxicological properties.  It uses machine-learning models trained on '
            'curated experimental data and complements AeroSysProp\'s group-contribution '
            'methods well — particularly for cross-validation of boiling point and '
            'vapour pressure estimates.'
        ),
        homepage='https://github.com/NIEHS/OPERA',
        status='file_export',
        # Optional OPERA HTTP service (opera/ container). When it is running,
        # the External Services page can POST SMILES and show results inline.
        # Native app -> localhost; containerised app -> the `opera` compose
        # service (overridden via the OPERA_URL env var in docker-compose.yml).
        service_url_env='OPERA_URL',
        default_service_url='http://localhost:8601',
        endpoint_options=_OPERA_ENDPOINTS,
        export_filename='aerosysprop_for_opera.smi',
        export_format_note=(
            'OPERA accepts a plain SMILES file — one compound per line, '
            'space-separated SMILES and name.  Download the file below and pass '
            'it to the OPERA CLI as your input.'
        ),
        install_notes=(
            'Download the latest OPERA release for your platform from the link above.  '
            'Once installed, run:\n\n'
            '```\n'
            'OPERA2.x_CL -s aerosysprop_for_opera.smi -o opera_results/\n'
            '```\n\n'
            'The output directory will contain a CSV with predictions for all compounds.  '
            'Direct invocation from this UI will be enabled once the exact CLI flags '
            'have been verified against a local OPERA installation.'
        ),
    ),

    # ---------------------------------------------------------------------------
    # Template for the next service — copy, fill in, and uncomment:
    # ---------------------------------------------------------------------------
    # ExternalService(
    #     id='my_service',
    #     name='My Service',
    #     tagline='One-line description',
    #     description='Longer description shown in the UI.',
    #     homepage='https://...',
    #     status='coming_soon',       # change to 'file_export' or 'api' when ready
    # ),

]

# Convenience lookup
SERVICE_BY_ID = {s.id: s for s in REGISTRY}
