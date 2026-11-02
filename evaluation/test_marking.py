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
    # A campaign run on Windows records the question set as evaluation\questions.yaml.
    camp = tmp_path / 'camp'
    (camp / 'transcripts').mkdir(parents=True)
    (camp / 'manifest.json').write_text(json.dumps(
        {'question_set': {'path': 'evaluation\\questions.yaml'}}))
    assert mark.main([str(camp)]) == 0
