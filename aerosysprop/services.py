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
    # Filled in once the CLI interface has been verified against a real install:
    # subprocess_cmd: Optional[Callable[[str, str], List[str]]] = None
    # Arguments: (input_smi_path, output_dir) -> argv list


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
