"""Batch runs: run folders, provenance manifest, summaries and write rules."""
import csv
import hashlib
import json
import time

import pytest

from aerosysprop import api, runs

MOLS = ['CCO', 'CCCO', 'CCCCO', 'CC(=O)O', 'OC(=O)CCC(=O)O', 'CC[N+](=O)[O-]',
        'CC(C)=O', 'c1ccccc1O', 'OCC(O)CO', 'CCOCC']


def _wait(run_id, timeout=120):
    t = time.time()
    while time.time() - t < timeout:
        st = runs.batch_run_status(run_id)
        if st['status'] != 'running':
            return st
        time.sleep(0.2)
    raise AssertionError('run did not finish')


@pytest.fixture
def dirs(tmp_path, monkeypatch):
    data, out = tmp_path / 'data', tmp_path / 'out'
    data.mkdir()
    out.mkdir()
    (data / 'mols.txt').write_text(''.join(f'm{i}\t{s}\n' for i, s in enumerate(MOLS * 6)))
    monkeypatch.setenv(api.DATA_DIR_ENV, str(data))
    monkeypatch.setenv(runs.OUTPUT_DIR_ENV, str(out))
    yield data, out
    runs._jobs.clear()


def test_disabled_without_output_folder(monkeypatch):
    monkeypatch.delenv(runs.OUTPUT_DIR_ENV, raising=False)
    with pytest.raises(api.ApiError) as e:
        runs.start_batch_run(smiles=['CCO'])
    assert e.value.code == 'BATCH_RUNS_DISABLED'


def test_file_run_writes_folder_with_provenance(dirs):
    data, out = dirs
    started = runs.start_batch_run(input_path='mols.txt', methods=['vapour_pressure'],
                                   label='Test VP run')
    assert started['status'] == 'running' and started['molecules'] == 60
    st = _wait(started['run_id'])
    assert st['status'] == 'complete', st
    d = out / started['run_id']
    assert {p.name for p in d.iterdir()} == {'results.csv', 'diagnostics.csv', 'run.json'}
    m = json.loads((d / 'run.json').read_text())
    # the manifest fingerprints the input and every output
    assert m['request']['input']['sha256'] == hashlib.sha256((data / 'mols.txt').read_bytes()).hexdigest()
    for f in ('results.csv', 'diagnostics.csv'):
        assert m['outputs'][f]['sha256'] == hashlib.sha256((d / f).read_bytes()).hexdigest()
    assert m['provenance']['methods']['vp_simpol']['label'].startswith('vp_simpol@1 (')
    assert m['provenance']['software']['aerosysprop_version']
    rows = list(csv.DictReader(open(d / 'results.csv')))
    assert len(rows) == 60 and 'vp_simpol_log10atm' in rows[0]
    # nitroethane is flagged, and the summary reports method spread
    s = st['summary']
    assert s['molecules'] == 60 and s['flags_by_method'].get('vp_simpol') == 6
    assert 'Vapour Pressure (log10(atm))' in s['method_spread']


def test_results_match_run_all(dirs):
    from rdkit import Chem
    from aerosysprop import batch
    st = _wait(runs.start_batch_run(smiles=MOLS, methods=['bp_nannoolal'])['run_id'])
    got = runs.read_batch_results(st['run_id'], limit=200)['rows']
    ref = batch.run_all([(s, Chem.MolFromSmiles(s)) for s in MOLS], methods=['bp_nannoolal'],
                        use_threads=True, progress_callback=lambda *_: None)
    assert [r['bp_nannoolal_K'] for r in got] == pytest.approx([r['bp_nannoolal_K'] for r in ref])


def test_views_and_tamper_detection(dirs):
    _, out = dirs
    rid = _wait(runs.start_batch_run(input_path='mols.txt', methods=['vapour_pressure'])['run_id'])['run_id']
    flagged = runs.read_batch_results(rid, view='flagged', method='vp_simpol')
    assert flagged['total_matching'] == 6
    spread = runs.read_batch_results(rid, view='largest_spread',
                                     group='Vapour Pressure (log10(atm))', limit=3)
    s = [r['spread'] for r in spread['rows']]
    assert s == sorted(s, reverse=True)
    page = runs.read_batch_results(rid, columns=['vp_simpol_log10atm'], limit=5, offset=5)
    assert page['returned'] == 5 and set(page['rows'][0]) <= {
        'name', 'canonical_smiles', 'vp_simpol_log10atm', 'flagged_methods'}
    assert page['results_match_manifest'] is True
    f = out / rid / 'results.csv'
    f.write_text(f.read_text().replace('m0', 'edited', 1))
    assert runs.read_batch_results(rid)['results_match_manifest'] is False


def test_every_run_gets_a_new_folder_and_earlier_runs_are_untouched(dirs):
    _, out = dirs
    a = _wait(runs.start_batch_run(smiles=['CCO'])['run_id'])['run_id']
    before = {p.name: p.read_bytes() for p in (out / a).iterdir()}
    b = _wait(runs.start_batch_run(smiles=['CCO'])['run_id'])['run_id']
    assert a != b
    assert {p.name: p.read_bytes() for p in (out / a).iterdir()} == before
    assert [r['run_id'] for r in runs.list_batch_runs()['runs']] == sorted([a, b], reverse=True)


def test_only_one_run_at_a_time(dirs, monkeypatch):
    real = api._compute

    def slow(*a, **k):
        time.sleep(1.5)
        return real(*a, **k)
    monkeypatch.setattr(api, '_compute', slow)
    first = runs.start_batch_run(smiles=['CCO'])['run_id']
    with pytest.raises(api.ApiError) as e:
        runs.start_batch_run(smiles=['CCO'])
    assert e.value.code == 'BUSY'
    _wait(first)


def test_failed_run_is_recorded(dirs, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError('simulated failure')
    monkeypatch.setattr(api, '_compute', boom)
    st = _wait(runs.start_batch_run(smiles=['CCO'])['run_id'])
    assert st['status'] == 'failed' and 'simulated failure' in st['error']


def test_interrupted_run_is_reported(dirs):
    _, out = dirs
    rid = _wait(runs.start_batch_run(smiles=['CCO'])['run_id'])['run_id']
    m = json.loads((out / rid / 'run.json').read_text())
    m.update(status='running', pid=2 ** 22 + 12345)  # a process that does not exist
    (out / rid / 'run.json').write_text(json.dumps(m))
    runs._jobs.pop(rid)
    assert runs.batch_run_status(rid)['status'] == 'interrupted'


@pytest.mark.parametrize('bad', ['../etc', '2026-01-01_000000__zzzzzz', '', 'x/y'])
def test_run_ids_are_validated(dirs, bad):
    with pytest.raises(api.ApiError) as e:
        runs.batch_run_status(bad)
    assert e.value.code == 'UNKNOWN_RUN'


def test_reading_an_incomplete_run_is_refused(dirs, monkeypatch):
    real = api._compute
    monkeypatch.setattr(api, '_compute', lambda *a, **k: (time.sleep(1.0), real(*a, **k))[1])
    rid = runs.start_batch_run(smiles=['CCO'])['run_id']
    with pytest.raises(api.ApiError) as e:
        runs.read_batch_results(rid)
    assert e.value.code == 'RUN_NOT_COMPLETE'
    _wait(rid)


def test_radicals_are_counted(dirs):
    st = _wait(runs.start_batch_run(smiles=['CCO', 'CC[O]', 'CCO[O]'], methods=['vp_simpol'])['run_id'])
    assert st['summary']['radical_species'] == 2
    rows = runs.read_batch_results(st['run_id'])['rows']
    assert [r['radical_electrons'] for r in rows] == ['0', '1', '1']
