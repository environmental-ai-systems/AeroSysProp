# This file is part of AeroSysProp.
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
"""
aerosysprop.mcp_server
======================
Model Context Protocol (MCP) server exposing AeroSysProp to AI agents.

Run it with::

    aerosysprop-mcp                      # stdio, for Claude Desktop / Claude Code / local agents
    aerosysprop-mcp --http --port 8765   # Streamable HTTP, for networked clients

Every tool is a thin wrapper over :mod:`aerosysprop.api` or
:mod:`aerosysprop.runs`. Only batch runs write anything: each run creates a new
folder inside ``AEROSYSPROP_OUTPUT_DIR`` (disabled if unset), and nothing is
ever overwritten, deleted, or written to a session database. File reading is
limited to ``AEROSYSPROP_DATA_DIR`` (disabled if unset).

Nothing may be printed to stdout in this module: with the stdio transport,
stdout carries the protocol messages.
"""
from __future__ import annotations

import argparse
import functools
import logging
import os
import warnings
from pathlib import Path
from typing import Annotated, Literal, Optional

from pydantic import Field

try:
    from mcp.server.mcpserver import MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
    from mcp.types import ToolAnnotations
except ImportError as exc:  # pragma: no cover - exercised only without the extra
    raise ImportError(
        'The AeroSysProp MCP server needs the optional MCP dependency. '
        'Install it with:  pip install "aerosysprop[mcp]"'
    ) from exc

from aerosysprop import __version__, api, runs

_PKG = Path(__file__).parent
_GUIDE = _PKG / 'agent_guide.md'
_AIOMFAC_DOC = _PKG / 'methods' / 'activity_coefficients.md'

INSTRUCTIONS = """\
AeroSysProp estimates thermophysical properties of organic molecules relevant \
to atmospheric aerosol (boiling points, vapour pressures, critical properties, \
liquid densities) from SMILES using published group-contribution methods, and \
AIOMFAC activity coefficients for organic/water/salt mixtures.

How to use it well:
1. Call list_methods once to see method IDs, units and dependencies.
2. Give molecules as SMILES strings. If unsure a SMILES is right, call \
diagnose_smiles first.
3. Always read the diagnostics returned with each value. A coverage warning \
means the method ignored part of the molecule, so its value may be unreliable; \
say so when reporting it.
4. Vapour pressures are log10(atm). Temperatures are kelvin.
5. For mixtures with salts, call list_salts for accepted names and give \
mass fractions (neutrals first, then salts).
6. For more than about 50 molecules, use start_batch_run instead of \
predict_properties: it writes a run folder (CSV plus a provenance manifest) \
and returns a summary. Poll batch_run_status, then use read_batch_results \
to look at flagged molecules or where methods disagree most. Follow the \
summary's reporting_requirements: give the run folder, quote provenance \
labels, state the coverage with any aggregate, and look up specific \
molecules by name or SMILES rather than inferring them.
Values are estimates, not measurements. When methods disagree, report the \
spread rather than picking one silently. The resource aerosysprop://guide has \
fuller guidance."""

server = MCPServer(
    name='aerosysprop',
    title='AeroSysProp',
    description='Thermophysical property and AIOMFAC activity coefficient '
                'estimation for atmospheric organic compounds.',
    instructions=INSTRUCTIONS,
    version=__version__,
)

_READ_ONLY = ToolAnnotations(read_only_hint=True, destructive_hint=False,
                             idempotent_hint=True, open_world_hint=False)


def _wrap(fn):
    """Translate ApiError into ToolError so the agent sees the code, the
    message and the hint, and can correct its call."""
    @functools.wraps(fn)
    def inner(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except api.ApiError as e:
            raise ToolError(e.as_text()) from e
    return inner


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

SmilesList = Annotated[list[str], Field(
    description='SMILES strings, one per molecule, e.g. ["CCO", "OC(=O)CCC(=O)O"].')]
Temperature = Annotated[float, Field(
    description='Temperature in kelvin (25 °C = 298.15 K).', gt=0)]


@server.tool(title='List property methods', annotations=_READ_ONLY)
@_wrap
def list_methods() -> dict:
    """List every property method AeroSysProp offers: IDs, units, output
    columns, dependencies, valid boiling point sources and mixture methods.
    Call this first to learn the method IDs used by predict_properties."""
    return api.list_methods()


@server.tool(title='List AIOMFAC salts', annotations=_READ_ONLY)
@_wrap
def list_salts() -> dict:
    """List the inorganic salts AIOMFAC supports, with their cation, anion,
    stoichiometry and molar mass. Use these names in activity_coefficients."""
    return api.list_salts()


@server.tool(title='Predict pure-component properties', annotations=_READ_ONLY)
@_wrap
def predict_properties(
    smiles: SmilesList,
    methods: Annotated[Optional[list[str]], Field(
        description='Method IDs from list_methods (e.g. "vp_simpol", "bp_nannoolal") '
                    'or property names ("boiling_point", "vapour_pressure", '
                    '"critical_properties", "liquid_density"). Omit to run all methods. '
                    'Dependencies such as the boiling point are added automatically.')] = None,
    temperature_K: Temperature = 298.15,
    bp_sources: Annotated[Optional[list[str]], Field(
        description='Boiling point method(s) feeding the vapour pressure, critical '
                    'property and density methods that need one: "bp_nannoolal" '
                    '(default), "bp_stein", "bp_joback". Give several to compare; '
                    'dependent columns then get a "__bp<Source>" suffix.')] = None,
    names: Annotated[Optional[list[str]], Field(
        description='Optional labels for the molecules, same length as smiles.')] = None,
    include_diagnostics: Annotated[bool, Field(
        description='Attach matcher-coverage diagnostics to each molecule. Leave on: '
                    'with it off, a molecule showing no diagnostics has simply not been '
                    'checked, and must not be reported as unproblematic.')] = True,
) -> dict:
    """Estimate pure-component properties (boiling point K, vapour pressure
    log10(atm), critical temperature K and pressure bar, liquid density g/cm3)
    for up to 1000 molecules using group-contribution methods.

    Each result carries canonical SMILES, formula, the values, any outputs
    that could not be computed, and per-method diagnostics. A diagnostic of
    code UNMATCHED_ATOMS means that method ignored part of the molecule."""
    return api.predict_properties(smiles, names=names, methods=methods,
                                  temperature_K=temperature_K, bp_sources=bp_sources,
                                  include_diagnostics=include_diagnostics)


@server.tool(title='Check SMILES before predicting', annotations=_READ_ONLY)
@_wrap
def diagnose_smiles(
    smiles: SmilesList,
    methods: Annotated[Optional[list[str]], Field(
        description='Limit the check to these method IDs or property names. '
                    'Omit to check every method.')] = None,
    include_aiomfac: Annotated[bool, Field(
        description='Also check decomposition into AIOMFAC groups.')] = True,
) -> dict:
    """Check SMILES validity and, for each method, whether every heavy atom is
    recognised by that method's functional groups. Use it to choose methods
    suited to a molecule, or to explain a suspicious value. No properties are
    computed."""
    return api.diagnose_smiles(smiles, methods=methods, include_aiomfac=include_aiomfac)


@server.tool(title='AIOMFAC activity coefficients', annotations=_READ_ONLY)
@_wrap
def activity_coefficients(
    smiles: Annotated[list[str], Field(
        description='SMILES of each neutral component, e.g. ["O", "CCO"] for water '
                    'and ethanol. Water is "O".')],
    temperature_K: Temperature = 298.15,
    mole_fractions: Annotated[Optional[list[float]], Field(
        description='Mole fraction of each neutral component, same order as smiles. '
                    'Only for salt-free mixtures; runs short-range AIOMFAC.')] = None,
    mass_fractions: Annotated[Optional[list[float]], Field(
        description='Mass fractions: one per neutral component, then one per salt, '
                    'summing to 1. Required when salts are present.')] = None,
    salts: Annotated[Optional[list[str]], Field(
        description='Salt names from list_salts, e.g. ["NaCl", "(NH4)2SO4"].')] = None,
    bisulfate: Annotated[bool, Field(
        description='Solve the HSO4- <-> H+ + SO4-- equilibrium when H+ and SO4-- '
                    'are both present. Default true.')] = True,
) -> dict:
    """Compute AIOMFAC activity coefficients for a liquid mixture of water,
    organics and optionally inorganic salts. Give exactly one of
    mole_fractions (no salts) or mass_fractions.

    Returns activity coefficients for each neutral component (mole-fraction
    basis), each ion (molality basis) and the mean ionic coefficient per salt,
    with the pre-flight diagnostics. Calls with diagnostic errors are refused
    with the reasons listed."""
    return api.activity_coefficients(smiles, temperature_K=temperature_K,
                                     mole_fractions=mole_fractions,
                                     mass_fractions=mass_fractions,
                                     salts=salts, bisulfate=bisulfate)


@server.tool(title='List data files', annotations=_READ_ONLY)
@_wrap
def list_data_files() -> dict:
    """List SMILES files (.txt .prop .smi .csv .tsv), SDF files and session
    databases (.db) in the configured data directory. Paths returned here are
    what load_smiles_file and read_session_database expect."""
    return api.list_data_files()


@server.tool(title='Load molecules from a file', annotations=_READ_ONLY)
@_wrap
def load_smiles_file(
    path: Annotated[str, Field(description='Path relative to the data directory, '
                                           'as returned by list_data_files.')],
    smiles_column: Annotated[Optional[int], Field(
        description='0-based SMILES column for text files. Omit to auto-detect.')] = None,
    name_column: Annotated[Optional[int], Field(
        description='0-based name column for text files. Omit to auto-detect.')] = None,
    table: Annotated[str, Field(description='Table name for SQLite files.')] = 'compounds',
    max_molecules: Annotated[int, Field(description='Maximum molecules to return.',
                                        ge=1, le=api.MAX_MOLECULES)] = 200,
) -> dict:
    """Read molecule names and SMILES from a file in the data directory, ready
    to pass to predict_properties. Nothing is computed or written."""
    return api.load_smiles_file(path, smiles_column=smiles_column,
                                name_column=name_column, table=table,
                                max_molecules=max_molecules)


@server.tool(title='Read saved session results', annotations=_READ_ONLY)
@_wrap
def read_session_database(
    path: Annotated[str, Field(description='Session database path relative to the '
                                           'data directory.')],
    what: Annotated[Literal['results', 'compounds', 'mixtures', 'mixture_species'], Field(
        description='"results" for saved property predictions, "compounds" for the '
                    'molecule list, "mixtures" for saved AIOMFAC runs, '
                    '"mixture_species" for one run\'s species (needs mixture_id).')] = 'results',
    methods: Annotated[Optional[list[str]], Field(
        description='For "results": restrict to these output columns, e.g. '
                    '["bp_nannoolal_K"].')] = None,
    temperature_K: Annotated[Optional[float], Field(
        description='For "results": only results computed at this temperature.')] = None,
    mixture_id: Annotated[Optional[int], Field(
        description='For "mixture_species": the mixture ID.')] = None,
    limit: Annotated[int, Field(description='Maximum rows to return.', ge=1, le=5000)] = 500,
) -> dict:
    """Read results previously saved by the AeroSysProp app. The database is
    opened read-only; this server cannot add, change or delete saved results."""
    return api.read_session_database(path, what=what, methods=methods,
                                     temperature_K=temperature_K,
                                     mixture_id=mixture_id, limit=limit)


@server.tool(title='OPERA QSAR predictions',
             annotations=ToolAnnotations(read_only_hint=True, destructive_hint=False,
                                         idempotent_hint=True, open_world_hint=True))
@_wrap
def opera_predict(
    smiles: Annotated[list[str], Field(description='SMILES strings, up to 100.')],
    endpoints: Annotated[Optional[list[str]], Field(
        description='OPERA endpoint codes, e.g. ["BP", "VP", "logP", "HL"]. Group '
                    'codes: "PC" physicochemical, "EF" environmental fate, "Tox", '
                    '"ADME". Omit for all (slow).')] = None,
    names: Annotated[Optional[list[str]], Field(
        description='Optional labels, same length as smiles.')] = None,
) -> dict:
    """Get independent QSAR predictions from the US EPA / NIEHS OPERA models
    through the optional OPERA container, for example to cross-check
    AeroSysProp boiling points and vapour pressures. Returns a
    SERVICE_UNAVAILABLE error if the container is not running. Note OPERA's
    units differ (boiling point in °C, vapour pressure in log10 mmHg)."""
    return api.opera_predict(smiles, names=names, endpoints=endpoints)


_BATCH_WRITE = ToolAnnotations(read_only_hint=False, destructive_hint=False,
                               idempotent_hint=False, open_world_hint=False)


@server.tool(title='Start a batch run', annotations=_BATCH_WRITE)
@_wrap
def start_batch_run(
    input_path: Annotated[Optional[str], Field(
        description='A file in the data directory (from list_data_files) holding the '
                    'molecules. Give this or smiles, not both.')] = None,
    smiles: Annotated[Optional[list[str]], Field(
        description='SMILES strings, when not reading from a file.')] = None,
    names: Annotated[Optional[list[str]], Field(
        description='Optional labels, same length as smiles.')] = None,
    methods: Annotated[Optional[list[str]], Field(
        description='Method IDs or property names, as for predict_properties. Omit for '
                    'all methods.')] = None,
    temperature_K: Temperature = 298.15,
    bp_sources: Annotated[Optional[list[str]], Field(
        description='Boiling point source(s), as for predict_properties.')] = None,
    label: Annotated[Optional[str], Field(
        description='Short label for the run folder name, e.g. "mcm vapour pressure".')] = None,
    force_new: Annotated[bool, Field(
        description='Run again even if an identical completed run exists. Leave false: '
                    'an identical request returns the existing run, avoiding duplicates.')] = False,
) -> dict:
    """Start a background batch run for a large set of molecules (up to 100,000).
    Results are written to a new run folder (results.csv, diagnostics.csv and a
    run.json provenance manifest) in the output folder; nothing is overwritten.
    Returns a run_id at once; poll batch_run_status until complete. The full
    MCM (about 5,800 species) takes around 20 s on 4 cores."""
    return runs.start_batch_run(input_path=input_path, smiles=smiles, names=names,
                                methods=methods, temperature_K=temperature_K,
                                bp_sources=bp_sources, label=label, force_new=force_new,
                                requested_by='aerosysprop-mcp')


@server.tool(title='Batch run status and summary', annotations=_READ_ONLY)
@_wrap
def batch_run_status(
    run_id: Annotated[str, Field(description='The run_id from start_batch_run.')],
) -> dict:
    """Progress of a batch run, or once complete its summary: a coverage statement,
    caveats and reporting requirements, the range of every output column for clean
    molecules (closed-shell, unflagged) and for all molecules, and where methods
    disagree most, with provenance labels for every method. Follow the reporting
    requirements, and pass the coverage statement and caveats on with any
    aggregate figure."""
    return runs.batch_run_status(run_id)


@server.tool(title='Read batch run results', annotations=_READ_ONLY)
@_wrap
def read_batch_results(
    run_id: Annotated[str, Field(description='The run_id of a completed run.')],
    view: Annotated[Literal['rows', 'flagged', 'largest_spread'], Field(
        description='"rows" for all molecules in order, "flagged" for molecules with '
                    'diagnostic warnings, "largest_spread" for molecules where methods '
                    'disagree most within a group.')] = 'rows',
    group: Annotated[Optional[str], Field(
        description='For largest_spread: a group name exactly as listed under '
                    'method_spread in the run summary, e.g. "Vapour Pressure (log10(atm))".')] = None,
    method: Annotated[Optional[str], Field(
        description='For flagged: only molecules flagged by this method ID.')] = None,
    columns: Annotated[Optional[list[str]], Field(
        description='Only return these value columns (plus name and SMILES).')] = None,
    names: Annotated[Optional[list[str]], Field(
        description='Exact molecule names to look up. Use this (or smiles) whenever '
                    'reporting a specific molecule.')] = None,
    smiles: Annotated[Optional[list[str]], Field(
        description='SMILES to look up; matched after canonicalisation.')] = None,
    name_contains: Annotated[Optional[str], Field(
        description='Only molecules whose name contains this text.')] = None,
    limit: Annotated[int, Field(description='Rows to return (at most 200).', ge=1, le=200)] = 50,
    offset: Annotated[int, Field(description='Rows to skip, for paging.', ge=0)] = 0,
) -> dict:
    """Read a slice of a completed batch run without loading the whole table.
    Look up specific molecules by exact name or SMILES; never infer a molecule's
    value from the summary. Also reports whether results.csv still matches the
    fingerprint recorded when the run wrote it."""
    return runs.read_batch_results(run_id, view=view, columns=columns, method=method,
                                   name_contains=name_contains, group=group,
                                   limit=limit, offset=offset, names=names, smiles=smiles)


@server.tool(title='List batch runs', annotations=_READ_ONLY)
@_wrap
def list_batch_runs(
    limit: Annotated[int, Field(description='How many runs to list, newest first.',
                                ge=1, le=200)] = 20,
) -> dict:
    """List previous batch runs in the output folder with their status, methods,
    temperature and size on disk, plus the folder's total size."""
    return runs.list_batch_runs(limit=limit)


# ---------------------------------------------------------------------------
# Resources: reference documents the agent can read on demand
# ---------------------------------------------------------------------------

@server.resource('aerosysprop://guide', name='agent_guide', title='Guide for agents',
                 description='How to choose methods and interpret AeroSysProp results.',
                 mime_type='text/markdown')
def agent_guide() -> str:
    return _GUIDE.read_text(encoding='utf-8')


@server.resource('aerosysprop://docs/aiomfac', name='aiomfac_reference',
                 title='AIOMFAC implementation notes',
                 description='Theory, validation and known limitations of the AIOMFAC '
                             'implementation.',
                 mime_type='text/markdown')
def aiomfac_reference() -> str:
    return _AIOMFAC_DOC.read_text(encoding='utf-8')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='aerosysprop-mcp',
        description='Run the AeroSysProp MCP server (stdio by default).')
    parser.add_argument('--http', action='store_true',
                        help='Serve over Streamable HTTP instead of stdio.')
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--port', type=int, default=8765)
    args = parser.parse_args(argv)
    # Logs go to stderr; stdout is reserved for the protocol on stdio.
    logging.basicConfig(level=logging.WARNING)
    # Diagnostics are returned in results and written to diagnostics.csv, so the
    # per-molecule UserWarnings the methods also raise would only flood the log.
    # Setting it in the environment reaches batch-run worker processes too.
    os.environ.setdefault('PYTHONWARNINGS', 'ignore::UserWarning')
    warnings.simplefilter('ignore', UserWarning)
    if args.http:
        server.run('streamable-http', host=args.host, port=args.port)
    else:
        server.run('stdio')


if __name__ == '__main__':
    main()
