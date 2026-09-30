"""The guidance level changes descriptive text only; values are identical.

The agent evaluation compares answers given with full and bare guidance, so
anything other than documentation differing between the levels would
confound it. These tests compare the levels byte for byte.
"""
import json
import subprocess
import sys
import textwrap
import time

import pytest

from aerosysprop import api, guidance, runs

MOLS = ['CCO', 'OCC(=O)O', 'CC[N+](=O)[O-]', 'OC(=O)CCC(=O)O', 'CC(C)=O', '[CH2]CO']


def _canon(obj) -> str:
    return json.dumps(obj, sort_keys=True)


def _at(level, monkeypatch, fn, *args, **kwargs):
    monkeypatch.setenv(guidance.GUIDANCE_ENV, level)
    return fn(*args, **kwargs)


# ---------------------------------------------------------------------------
# The level itself
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('value,expected', [(None, 'full'), ('', 'full'), ('full', 'full'),
                                            ('bare', 'bare'), (' BARE ', 'bare')])
def test_level_from_environment(monkeypatch, value, expected):
    if value is None:
        monkeypatch.delenv(guidance.GUIDANCE_ENV, raising=False)
    else:
        monkeypatch.setenv(guidance.GUIDANCE_ENV, value)
    assert guidance.level() == expected


def test_unknown_level_is_refused(monkeypatch):
    monkeypatch.setenv(guidance.GUIDANCE_ENV, 'minimal')
    with pytest.raises(ValueError, match='AEROSYSPROP_GUIDANCE'):
        guidance.level()


# ---------------------------------------------------------------------------
# API: values identical, only descriptive text removed
# ---------------------------------------------------------------------------

def test_predict_properties_identical_apart_from_interpretation(monkeypatch):
    kw = dict(methods=['vapour_pressure', 'vp_evaporation@1', 'vp_evaporation@2'])
    full = _at('full', monkeypatch, api.predict_properties, MOLS, **kw)
    bare = _at('bare', monkeypatch, api.predict_properties, MOLS, **kw)
    assert 'interpretation' in full and 'interpretation' not in bare
    del full['interpretation']
    assert _canon(full) == _canon(bare)


@pytest.mark.parametrize('call', [
    lambda: api.diagnose_smiles(MOLS),
    lambda: api.activity_coefficients(['O', 'CCO'], mole_fractions=[0.8, 0.2]),
    lambda: api.activity_coefficients(['O', 'CCO'], salts=['NaCl'],
                                      mass_fractions=[0.7, 0.2, 0.1]),
    lambda: api.list_salts(),
])
def test_other_results_are_byte_identical(monkeypatch, call):
    full = _at('full', monkeypatch, call)
    bare = _at('bare', monkeypatch, call)
    assert _canon(full) == _canon(bare)


def test_errors_are_identical(monkeypatch):
    def err(level):
        monkeypatch.setenv(guidance.GUIDANCE_ENV, level)
        with pytest.raises(api.ApiError) as e:
            api.predict_properties(['CCO'], methods=['nonsense'])
        return e.value.to_dict()
    assert err('full') == err('bare')


def test_bare_list_methods_keeps_ids_columns_and_units_only(monkeypatch):
    full = _at('full', monkeypatch, api.list_methods)
    bare = _at('bare', monkeypatch, api.list_methods)
    assert set(bare) == {'aerosysprop_version', 'single_component_methods', 'property_groups',
                         'boiling_point_sources', 'mixture_methods'}
    assert [m['id'] for m in bare['single_component_methods']] == \
           [m['id'] for m in full['single_component_methods']]
    for f, b in zip(full['single_component_methods'], bare['single_component_methods']):
        assert set(b) == {'id', 'output_columns'}
        assert b['output_columns'] == f['output_columns']
    # Parameter sets remain visible only through the column names.
    evap = next(b for b in bare['single_component_methods'] if b['id'] == 'vp_evaporation')
    assert {'vp_evaporation_log10atm__set1', 'vp_evaporation_log10atm__set2'} <= set(evap['output_columns'])
    assert bare['mixture_methods'] == [{'id': m['id']} for m in full['mixture_methods']]
    assert 'card' in full['single_component_methods'][0] and 'notes' in full


# ---------------------------------------------------------------------------
# Batch runs: the run folder is the same; only the served summary is trimmed
# ---------------------------------------------------------------------------

@pytest.fixture
def out_dir(tmp_path, monkeypatch):
    monkeypatch.setenv(runs.OUTPUT_DIR_ENV, str(tmp_path))
    yield tmp_path
    runs._jobs.clear()


def _finished(run_id):
    for _ in range(600):
        st = runs.batch_run_status(run_id)
        if st['status'] != 'running':
            return st
        time.sleep(0.1)
    raise AssertionError('run did not finish')


def test_batch_summary_trimmed_when_served_but_run_folder_unchanged(out_dir, monkeypatch):
    monkeypatch.setenv(guidance.GUIDANCE_ENV, 'bare')
    started = runs.start_batch_run(smiles=MOLS * 3, methods=['vapour_pressure'], label='g')
    bare = _finished(started['run_id'])
    assert bare['status'] == 'complete'
    manifest = out_dir / started['run_id'] / runs.MANIFEST_FILE
    written = manifest.read_bytes()

    full = _at('full', monkeypatch, runs.batch_run_status, started['run_id'])
    assert manifest.read_bytes() == written  # serving never rewrites run.json

    # run.json keeps the full summary even for a run started at the bare level.
    stored = json.loads(written)['summary']
    assert stored['caveats'] and stored['reporting_requirements']
    assert _canon(full['summary']) == _canon(stored)

    fs, bs = full['summary'], bare['summary']
    assert 'caveats' not in bs and 'reporting_requirements' not in bs
    assert 'statement' in fs['coverage'] and 'statement' not in bs['coverage']
    fs = dict(fs, coverage={k: v for k, v in fs['coverage'].items() if k != 'statement'})
    del fs['caveats'], fs['reporting_requirements']
    assert _canon(fs) == _canon(bs)
    for k in ('provenance', 'outputs', 'warnings', 'folder'):
        assert _canon(full[k]) == _canon(bare[k])


def test_identical_request_reuses_the_run_across_levels(out_dir, monkeypatch):
    first = _at('full', monkeypatch, runs.start_batch_run, smiles=MOLS, methods=['vp_simpol'])
    _finished(first['run_id'])
    again = _at('bare', monkeypatch, runs.start_batch_run, smiles=MOLS, methods=['vp_simpol'])
    assert again['reused'] and again['run_id'] == first['run_id']


# ---------------------------------------------------------------------------
# MCP server, started as it is in use: a fresh process per level
# ---------------------------------------------------------------------------

_PROBE = textwrap.dedent('''
    import asyncio, json, sys
    from mcp import Client
    from aerosysprop.mcp_server import server
    async def main():
        async with Client(server) as c:
            init = c.initialize_result if hasattr(c, 'initialize_result') else None
            tools = (await c.list_tools()).tools
            resources = (await c.list_resources()).resources
            res = await c.call_tool('predict_properties', {
                'smiles': ['OCC(=O)O', 'CC[N+](=O)[O-]'],
                'methods': ['vp_evaporation@1', 'vp_evaporation@2', 'vp_simpol']})
            sc = res.structured_content
            result = sc.get('result', sc) if isinstance(sc, dict) else json.loads(res.content[0].text)
            print(json.dumps({
                'instructions': server.instructions,
                'tools': {t.name: {'title': t.title, 'description': t.description,
                                   'schema': t.input_schema} for t in tools},
                'resources': sorted(str(r.uri) for r in resources),
                'result': result}))
    asyncio.run(main())
''')


def _probe(level):
    env = {**__import__('os').environ, guidance.GUIDANCE_ENV: level}
    p = subprocess.run([sys.executable, '-c', _PROBE], env=env, capture_output=True,
                       text=True, timeout=300)
    return p


@pytest.fixture(scope='module')
def served():
    pytest.importorskip('mcp')
    out = {}
    for level in guidance.LEVELS:
        p = _probe(level)
        assert p.returncode == 0, p.stderr[-2000:]
        out[level] = json.loads(p.stdout.strip().splitlines()[-1])
    return out


def _without_descriptions(schema):
    if isinstance(schema, dict):
        return {k: _without_descriptions(v) for k, v in schema.items() if k != 'description'}
    if isinstance(schema, list):
        return [_without_descriptions(v) for v in schema]
    return schema


def test_bare_server_serves_one_line_tools_and_no_resources(served):
    from aerosysprop.mcp_server import BARE_DESCRIPTIONS, BARE_INSTRUCTIONS, INSTRUCTIONS
    full, bare = served['full'], served['bare']
    assert set(full['tools']) == set(bare['tools']) == set(BARE_DESCRIPTIONS)
    for name, t in bare['tools'].items():
        assert t['description'] == BARE_DESCRIPTIONS[name]
        assert t['title'] is None
        assert '"description"' not in json.dumps(t['schema'])
        # Same arguments, types, defaults and bounds: validation is identical.
        assert t['schema'] == _without_descriptions(full['tools'][name]['schema'])
    assert bare['instructions'] == BARE_INSTRUCTIONS and full['instructions'] == INSTRUCTIONS
    assert bare['resources'] == [] and 'aerosysprop://guide' in full['resources']


def test_values_over_mcp_are_identical(served):
    full, bare = served['full']['result'], served['bare']['result']
    del full['interpretation']
    assert _canon(full) == _canon(bare)


def test_server_refuses_to_start_with_unknown_level():
    p = _probe('minimal')
    assert p.returncode != 0 and 'AEROSYSPROP_GUIDANCE' in p.stderr
