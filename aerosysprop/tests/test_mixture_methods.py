"""Tests for the mixture-method registry."""

import pytest

from aerosysprop.mixture_methods import (
    MIXTURE_METHODS, MixtureMethod, available_methods,
    method_by_id, method_by_name,
)


def test_registry_nonempty_and_has_aiomfac():
    assert MIXTURE_METHODS
    ids = {m.id for m in MIXTURE_METHODS}
    assert 'aiomfac' in ids


def test_all_entries_well_formed():
    for m in MIXTURE_METHODS:
        assert isinstance(m, MixtureMethod)
        assert m.id and m.name
        assert m.status in ('available', 'coming_soon')
        assert m.summary
        # available methods must document how to call them
        if m.status == 'available':
            assert m.python_example


def test_available_methods_filters_status():
    avail = available_methods()
    assert all(m.status == 'available' for m in avail)
    assert any(m.id == 'aiomfac' for m in avail)


def test_lookup_by_id_and_name():
    m = method_by_id('aiomfac')
    assert m.name == 'AIOMFAC activity coefficients'
    assert method_by_name(m.name).id == 'aiomfac'


def test_lookup_unknown_raises():
    with pytest.raises(KeyError):
        method_by_id('not_a_method')
    with pytest.raises(KeyError):
        method_by_name('Not A Method')


def test_aiomfac_has_diagnostics_notes():
    """The registry is where method-specific diagnostic caveats live so the
    docs page can render them generically."""
    m = method_by_id('aiomfac')
    assert m.diagnostics_notes
    joined = ' '.join(m.diagnostics_notes)
    assert 'BISULFATE_NEEDED' in joined
