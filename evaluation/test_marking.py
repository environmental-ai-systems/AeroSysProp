"""Regression tests for the evaluation's number reading and manifest paths.

    python -m pytest -q evaluation/test_marking.py
"""
import json

import pytest

from evaluation import mark


@pytest.mark.parametrize('text, want', [
    ('P ≈ 0.0649 atm ≈ 6,580 Pa', [0.0649, 6580.0]),        # thousands separator
    ('about 1,234,567.5 molecules', [1234567.5]),
    ('5,809 species and 4 MB', [5809.0, 4.0]),
    ('sets 1,2 and 3', [1.0, 2.0, 3.0]),                     # a list stays a list
    ('fractions 0.8,0.2', [0.8, 0.2]),
    ('12,34 and 1,2345', [12.0, 34.0, 1.0, 2345.0]),         # not three-digit groups
    ('log10 = −1.187, which is 0.065 atm', [-1.187, 0.065]),
    ('values -1.26, -1.22', [-1.26, -1.22]),
    ('8.6 × 10^2 Pa', [860.0]),
])
def test_numbers(text, want):
    assert [round(x, 6) for x in mark.numbers(text)] == want


def test_manifest_written_on_windows_is_readable(tmp_path):
    # A campaign run on Windows records the question set path with backslashes.
    camp = tmp_path / 'camp'
    (camp / 'transcripts').mkdir(parents=True)
    qfile = tmp_path / 'set' / 'questions.yaml'
    qfile.parent.mkdir()
    qfile.write_text('questions: []\n')
    (camp / 'manifest.json').write_text(json.dumps(
        {'question_set': {'path': str(qfile).replace('/', '\\')}}))
    assert mark.main([str(camp)]) == 0


@pytest.mark.parametrize('answer', [
    "The SMILES \"CC(C\" isn't valid. RDKit can't parse it.",
    'I can’t give a value because it is not valid SMILES.',
    'RDKit couldn’t parse it: the ring is never closed.',
    'That is not a valid SMILES string.',
    'RDKit rejects it because the parenthesis is never closed.',
])
def test_failure_wordings(answer):
    assert mark.FAILURE_RE.search(answer)


def _card_question(must):
    return {'id': 'prov-04', 'kind': 'open',
            'key': {'acceptable_methods': ['vp_evaporation'], 'must_report': must}}


def _transcript(answer, names):
    return {'status': 'complete', 'final_answer': answer,
            'tool_calls': [{'name': n, 'arguments': {}, 'is_error': False} for n in names]}


def test_card_question_answered_from_list_methods():
    # The card is the source: no prediction is needed, and years and counts need no unit.
    row = mark.mark(_transcript('Set 2 (Compernolle 2011) refits 805 groups.', ['list_methods']),
                    _card_question(['method_named', 'value_with_units']))
    assert row['method_ok'] == 1
    assert row.get('value_with_units') is None


def test_card_question_value_still_needs_units():
    row = mark.mark(_transcript('Set 2 gives -6.11 for pinic acid.', ['list_methods']),
                    _card_question(['value_with_units']))
    assert row['value_with_units'] == 0


def test_list_methods_is_not_a_prediction_elsewhere():
    q = {'id': 'lookup-01', 'kind': 'lookup', 'key': {'acceptable_methods': ['vp_simpol']}}
    assert mark.mark(_transcript('', ['list_methods']), q)['method_ok'] == 0


def test_judge_reference_fills_shared_facts():
    from evaluation import judge
    q = {'id': 'x', 'judge_reference': 'Check the card. {evaporation_sets}'}
    ref = judge.reference_for(q)
    assert ref.startswith('Check the card. EVAPORATION has two parameter sets')
    assert judge.reference_for({'id': 'y'}) == ''
