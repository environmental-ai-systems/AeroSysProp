"""Tests for the agent-facing API (aerosysprop.api)."""
import hashlib
import json
import sqlite3

import pytest
from rdkit import Chem

from aerosysprop import api, batch, db


def _strict_json(obj):
    """Fail if the object is not plain JSON (no NaN, numpy, tuples...)."""
    return json.loads(json.dumps(obj, allow_nan=False))


# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------

def test_list_methods_matches_registry_and_is_json():
    out = _strict_json(api.list_methods())
    ids = [m['id'] for m in out['single_component_methods']]
    assert ids == batch.ALL_METHODS
    vp = next(m for m in out['single_component_methods'] if m['id'] == 'vp_nannoolal')
    assert vp['output_columns'] == {'vp_nannoolal_log10atm': 'log10(atm)'}
    assert vp['depends_on'] == ['a boiling point source']
    assert out['mixture_methods'][0]['id'] == 'aiomfac'


def test_list_salts_has_ion_names_and_aliases():
    out = _strict_json(api.list_salts())
    by_name = {s['name']: s for s in out['salts']}
    assert by_name['NaCl']['cation'] == 'Na+'
    assert by_name['NaCl']['anion'] == 'Cl-'
    assert '(NH4)2SO4' in by_name['NH42SO4']['also_accepted']


# ---------------------------------------------------------------------------
# predict_properties
# ---------------------------------------------------------------------------

def test_predict_matches_run_all():
    smiles = ['CCO', 'OC(=O)CCC(=O)O']
    out = _strict_json(api.predict_properties(smiles, methods=['bp_nannoolal', 'vp_simpol']))
    ref = batch.run_all([(s, Chem.MolFromSmiles(s)) for s in smiles],
                        methods=['bp_nannoolal', 'vp_simpol'], use_threads=True,
                        progress_callback=lambda *_: None)
    for rec, row in zip(out['results'], ref):
        assert rec['values']['bp_nannoolal_K'] == pytest.approx(row['bp_nannoolal_K'])
        assert rec['values']['vp_simpol_log10atm'] == pytest.approx(row['vp_simpol_log10atm'])
    assert out['columns']['bp_nannoolal_K'] == 'K'


def test_predict_handles_invalid_smiles_without_failing_the_batch():
    out = _strict_json(api.predict_properties(['CCO', 'not_a_smiles'], methods=['bp_joback']))
    assert out['results'][0]['values']['bp_joback_K'] > 0
    assert out['results'][1]['error']['code'] == 'INVALID_SMILES'
    assert out['summary']['invalid_smiles'] == 1


def test_predict_duplicate_names_do_not_collide():
    out = api.predict_properties(['CCO', 'CCCO'], names=['x', 'x'], methods=['bp_joback'])
    a, b = (r['values']['bp_joback_K'] for r in out['results'])
    assert a != b


def test_predict_accepts_property_group_names():
    out = api.predict_properties(['CCO'], methods=['vapour_pressure'])
    for mid in batch.METHOD_GROUPS['Vapour Pressure']:
        assert mid in out['methods_run']


def test_predict_multiple_bp_sources_suffixes_columns():
    out = api.predict_properties(['CCCCO'], methods=['vp_nannoolal'],
                                 bp_sources=['nannoolal', 'stein'])
    assert 'vp_nannoolal_log10atm__bpNannoolal' in out['columns']
    assert 'vp_nannoolal_log10atm__bpStein' in out['columns']


def test_predict_reports_coverage_diagnostics():
    # A nitro group is outside several methods' group sets.
    out = _strict_json(api.predict_properties(['CC[N+](=O)[O-]'], methods=['vapour_pressure']))
    diag = out['results'][0]['diagnostics']
    assert {'vp_evaporation', 'vp_simpol'} <= set(diag)
    assert diag['vp_simpol'][0]['code'] == 'UNMATCHED_ATOMS'
    assert out['summary']['with_diagnostic_warnings'] == 1


@pytest.mark.parametrize('kwargs, code', [
    ({'methods': ['not_a_method']}, 'UNKNOWN_METHOD'),
    ({'bp_sources': ['bp_magic']}, 'UNKNOWN_METHOD'),
    ({'temperature_K': -5}, 'INVALID_ARGUMENT'),
    ({'names': ['a', 'b']}, 'INVALID_ARGUMENT'),
])
def test_predict_errors_are_structured(kwargs, code):
    with pytest.raises(api.ApiError) as e:
        api.predict_properties(['CCO'], **kwargs)
    assert e.value.code == code
    _strict_json(json.loads(e.value.as_text()))


def test_predict_limit(monkeypatch):
    monkeypatch.setattr(api, 'MAX_MOLECULES', 2)
    with pytest.raises(api.ApiError) as e:
        api._check_count(3, api.MAX_MOLECULES)
    assert e.value.code == 'TOO_MANY_MOLECULES'


def test_predict_never_writes_csv(monkeypatch):
    seen = {}
    real = batch.run_all

    def spy(*a, **k):
        seen.update(k)
        return real(*a, **k)
    monkeypatch.setattr(batch, 'run_all', spy)
    api.predict_properties(['CCO'], methods=['bp_joback'])
    assert seen.get('output_csv') is None


def test_diagnose_smiles():
    out = _strict_json(api.diagnose_smiles(['CCO', 'C1CC']))
    assert out['results'][0]['valid'] is True
    assert out['results'][0]['formula'] == 'C2H6O'
    assert out['results'][1]['valid'] is False
    assert 'density_tyn_calus' in out['methods_without_coverage_check']


# ---------------------------------------------------------------------------
# activity_coefficients
# ---------------------------------------------------------------------------

def test_aiomfac_sr_matches_library():
    from aerosysprop.methods.activity_coefficients import aiomfac_sr_from_smiles
    out = _strict_json(api.activity_coefficients(['O', 'CCO'], mole_fractions=[0.5, 0.5]))
    ref = aiomfac_sr_from_smiles(['O', 'CCO'], [0.5, 0.5], 298.15, strict='silent')
    got = [c['activity_coefficient'] for c in out['neutral_components']]
    assert got == pytest.approx(list(ref))
    assert out['composition_basis'] == 'mole_fraction'


def test_aiomfac_with_salts_and_alias():
    from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles
    out = _strict_json(api.activity_coefficients(
        ['O', 'CCO'], salts=['NaCl', '(NH4)2SO4'], mass_fractions=[0.6, 0.2, 0.1, 0.1]))
    ref = aiomfac_from_smiles(['O', 'CCO'], ['NaCl', 'NH42SO4'], [0.6, 0.2, 0.1, 0.1],
                              298.15, strict='silent')
    assert [s['salt'] for s in out['salts']] == ['NaCl', 'NH42SO4']
    assert [s['mean_ionic_activity_coefficient'] for s in out['salts']] == \
        pytest.approx(list(ref.gamma_mean))
    assert {i['ion'] for i in out['ions']} >= {'Na+', 'Cl-', 'NH4+'}


@pytest.mark.parametrize('kwargs, code', [
    ({}, 'INVALID_ARGUMENT'),
    ({'mole_fractions': [0.5, 0.5], 'mass_fractions': [0.5, 0.5]}, 'INVALID_ARGUMENT'),
    ({'salts': ['NaCl'], 'mole_fractions': [0.5, 0.5]}, 'INVALID_ARGUMENT'),
    ({'salts': ['Unobtainium'], 'mass_fractions': [0.5, 0.4, 0.1]}, 'UNKNOWN_SALT'),
    ({'mole_fractions': [1.0]}, 'INVALID_ARGUMENT'),
])
def test_aiomfac_argument_errors(kwargs, code):
    with pytest.raises(api.ApiError) as e:
        api.activity_coefficients(['O', 'CCO'], **kwargs)
    assert e.value.code == code


# ---------------------------------------------------------------------------
# File and database access: confined and read-only
# ---------------------------------------------------------------------------

def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def data_dir(tmp_path, monkeypatch):
    (tmp_path / 'mols.txt').write_text('ethanol\tCCO\nbutanol\tCCCCO\n')
    results = batch.run_all([('ethanol', Chem.MolFromSmiles('CCO'))],
                            methods=['bp_joback'], use_threads=True,
                            progress_callback=lambda *_: None)
    db.save_results(str(tmp_path / 'session.db'), results, 298.15)
    monkeypatch.setenv(api.DATA_DIR_ENV, str(tmp_path))
    return tmp_path


def test_file_access_disabled_without_env(monkeypatch):
    monkeypatch.delenv(api.DATA_DIR_ENV, raising=False)
    with pytest.raises(api.ApiError) as e:
        api.list_data_files()
    assert e.value.code == 'FILE_ACCESS_DISABLED'


def test_list_and_load_files(data_dir):
    files = {f['path']: f['kind'] for f in api.list_data_files()['files']}
    assert files == {'mols.txt': 'smiles_text', 'session.db': 'sqlite'}
    out = _strict_json(api.load_smiles_file('mols.txt'))
    assert [m['name'] for m in out['molecules']] == ['ethanol', 'butanol']


@pytest.mark.parametrize('bad', ['../outside.txt', '/etc/passwd', 'sub/../../x.db'])
def test_path_traversal_is_refused(data_dir, bad):
    with pytest.raises(api.ApiError) as e:
        api.load_smiles_file(bad)
    assert e.value.code in {'PATH_NOT_ALLOWED', 'FILE_NOT_FOUND'}


def test_session_database_read_is_read_only(data_dir):
    dbfile = data_dir / 'session.db'
    before = _sha(dbfile)
    out = _strict_json(api.read_session_database('session.db'))
    assert out['rows'][0]['bp_joback_K'] > 0
    assert out['columns'] == {'bp_joback_K': 'K'}
    api.read_session_database('session.db', what='compounds')
    assert _sha(dbfile) == before


def test_read_only_connection_rejects_writes(data_dir):
    conn = db._connect(str(data_dir / 'session.db'), read_only=True)
    with pytest.raises(sqlite3.OperationalError):
        conn.execute("INSERT INTO compounds (name, smiles) VALUES ('x', 'C')")
    conn.close()


def test_read_only_does_not_create_missing_database(tmp_path):
    missing = tmp_path / 'nope.db'
    with pytest.raises(FileNotFoundError):
        db.load_results(str(missing), read_only=True)
    assert not missing.exists()


def test_old_database_without_mixture_tables(data_dir):
    legacy = data_dir / 'legacy.db'
    conn = sqlite3.connect(legacy)
    conn.execute('CREATE TABLE compounds (id INTEGER PRIMARY KEY, name TEXT, smiles TEXT)')
    conn.commit(); conn.close()
    with pytest.raises(api.ApiError) as e:
        api.read_session_database('legacy.db', what='mixtures')
    assert e.value.code == 'READ_FAILED'


def test_api_exposes_no_write_functions():
    for name in api.__all__:
        assert not any(w in name for w in ('save', 'write', 'delete', 'create', 'export'))


# ---------------------------------------------------------------------------
# OPERA
# ---------------------------------------------------------------------------

def test_opera_unavailable_is_a_clear_error(monkeypatch):
    monkeypatch.setenv('OPERA_URL', 'http://127.0.0.1:9')  # nothing listens here
    with pytest.raises(api.ApiError) as e:
        api.opera_predict(['CCO'], endpoints=['BP'])
    assert e.value.code == 'SERVICE_UNAVAILABLE'


def test_opera_unknown_endpoint(monkeypatch):
    with pytest.raises(api.ApiError) as e:
        api.opera_predict(['CCO'], endpoints=['NOPE'])
    assert e.value.code == 'UNKNOWN_ENDPOINT'
