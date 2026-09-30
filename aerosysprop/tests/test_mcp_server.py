"""Tests for the MCP server, driven through an in-process MCP client."""
import asyncio
import json

import pytest

pytest.importorskip('mcp')

from mcp import Client  # noqa: E402

from aerosysprop.mcp_server import server  # noqa: E402

EXPECTED_TOOLS = {
    'list_methods', 'list_salts', 'predict_properties', 'diagnose_smiles',
    'activity_coefficients', 'list_data_files', 'load_smiles_file',
    'read_session_database', 'opera_predict',
    'start_batch_run', 'batch_run_status', 'read_batch_results', 'list_batch_runs',
}
WRITES = {'start_batch_run'}  # writes only new run folders; never destructive


def _run(coro_fn):
    async def go():
        async with Client(server) as client:
            return await coro_fn(client)
    return asyncio.run(go())


def _payload(result):
    """Structured content if present, else parse the first text block."""
    if getattr(result, 'structured_content', None):
        sc = result.structured_content
        return sc.get('result', sc) if isinstance(sc, dict) else sc
    return json.loads(result.content[0].text)


def test_tool_set_and_write_annotations():
    tools = _run(lambda c: c.list_tools()).tools
    names = {t.name for t in tools}
    assert names == EXPECTED_TOOLS
    for t in tools:
        assert t.annotations is not None, t.name
        assert t.annotations.read_only_hint is (t.name not in WRITES), t.name
        assert t.annotations.destructive_hint is False, t.name
        assert t.description and len(t.description) > 40, t.name


def test_predict_properties_over_mcp():
    async def call(c):
        return await c.call_tool('predict_properties',
                                 {'smiles': ['CCO'], 'methods': ['bp_nannoolal']})
    res = _run(call)
    assert not res.is_error
    out = _payload(res)
    assert out['results'][0]['values']['bp_nannoolal_K'] > 300


def test_activity_coefficients_over_mcp():
    async def call(c):
        return await c.call_tool('activity_coefficients', {
            'smiles': ['O', 'CCO'], 'salts': ['NaCl'], 'mass_fractions': [0.7, 0.2, 0.1]})
    out = _payload(_run(call))
    assert out['salts'][0]['salt'] == 'NaCl'


def test_api_error_reaches_agent_as_structured_tool_error():
    async def call(c):
        return await c.call_tool('predict_properties',
                                 {'smiles': ['CCO'], 'methods': ['nonsense']})
    res = _run(call)
    assert res.is_error
    text = res.content[0].text
    assert 'UNKNOWN_METHOD' in text and 'list_methods' in text


def test_guide_resource_is_served():
    async def call(c):
        return await c.read_resource('aerosysprop://guide')
    res = _run(call)
    assert 'log10' in res.contents[0].text


def test_batch_run_over_mcp(tmp_path, monkeypatch):
    monkeypatch.setenv('AEROSYSPROP_OUTPUT_DIR', str(tmp_path))

    async def go(c):
        started = _payload(await c.call_tool('start_batch_run', {
            'smiles': ['CCO', 'CC[N+](=O)[O-]'], 'methods': ['vp_simpol'], 'label': 'mcp test'}))
        for _ in range(300):
            st = _payload(await c.call_tool('batch_run_status', {'run_id': started['run_id']}))
            if st['status'] != 'running':
                break
            await asyncio.sleep(0.1)
        rows = _payload(await c.call_tool('read_batch_results', {
            'run_id': started['run_id'], 'view': 'flagged'}))
        return st, rows
    st, rows = _run(go)
    assert st['status'] == 'complete'
    assert st['summary']['coverage']['flagged_fraction_by_method'] == {'vp_simpol': 0.5}
    assert rows['total_matching'] == 1
