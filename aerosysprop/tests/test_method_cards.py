"""Method cards: every method is documented and its parameter files match
the fingerprints recorded when its parameter set was approved."""
import pytest

from aerosysprop import api, batch, method_cards


def test_every_method_has_a_card():
    cards = method_cards.load_cards()
    missing = [m for m in batch.ALL_METHODS + ['aiomfac'] if m not in cards]
    assert not missing, f'methods without a card: {missing}'


@pytest.mark.parametrize('method_id', sorted(method_cards.load_cards()))
def test_card_is_sound_and_parameters_unchanged(method_id):
    problems = method_cards.check_card(method_cards.card_for(method_id))
    assert not problems, '\n'.join(problems)


def test_same_model_links_point_to_real_cards():
    cards = method_cards.load_cards()
    for cid, card in cards.items():
        if 'same_model_as' in card:
            assert card['same_model_as'] in cards, cid


def test_fingerprint_detects_a_one_character_change(tmp_path, monkeypatch):
    f = tmp_path / 'data' / 'x.data'
    f.parent.mkdir()
    f.write_text('1 0.5\n')
    monkeypatch.setattr(method_cards, 'PACKAGE_DIR', tmp_path)
    method_cards.file_sha256.cache_clear()
    before = method_cards.file_sha256('data/x.data')
    f.write_text('1 0.6\n')
    method_cards.file_sha256.cache_clear()
    after = method_cards.file_sha256('data/x.data')
    method_cards.file_sha256.cache_clear()
    assert before != after and len(before) == 64


def test_changed_file_is_reported(tmp_path, monkeypatch):
    f = tmp_path / 'data' / 'x.data'
    f.parent.mkdir()
    f.write_text('1 0.5\n')
    monkeypatch.setattr(method_cards, 'PACKAGE_DIR', tmp_path)
    method_cards.file_sha256.cache_clear()
    card = {'id': 'x', 'name': 'x', 'property': 'p', 'status': 'in-house', 'model': {},
            'review': {'state': 'draft'},
            'parameter_set': {'version': 1,
                              'files': {'data/x.data': method_cards.file_sha256('data/x.data')}}}
    assert method_cards.check_card(card) == []
    f.write_text('1 0.6\n')
    method_cards.file_sha256.cache_clear()
    problems = method_cards.check_card(card)
    method_cards.file_sha256.cache_clear()
    assert len(problems) == 1 and 'has changed' in problems[0]


def test_provenance_label():
    prov = method_cards.provenance('vp_simpol')
    assert prov['label'].startswith('vp_simpol@1 (')
    assert prov['parameters_match_card'] is True
    assert method_cards.provenance('not_a_method') == {'method': 'not_a_method', 'documented': False}


def test_api_results_carry_provenance():
    out = api.predict_properties(['CCO'], methods=['vp_nannoolal'])
    prov = out['provenance']
    assert set(prov['methods']) == set(out['methods_run'])
    assert prov['methods']['vp_nannoolal']['parameter_set_version'] == 1
    assert 'aerosysprop_version' in prov['software']
    aio = api.activity_coefficients(['O', 'CCO'], mole_fractions=[0.5, 0.5])
    assert aio['provenance']['methods']['aiomfac']['documented'] is True


def test_list_methods_reports_card_status():
    out = api.list_methods()
    simpol = next(m for m in out['single_component_methods'] if m['id'] == 'vp_simpol')
    assert simpol['card']['status'] in method_cards.STATUSES
    assert simpol['card']['references']
