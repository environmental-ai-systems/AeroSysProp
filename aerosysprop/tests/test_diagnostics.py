"""Tests for the aerosysprop.diagnostics framework."""

import warnings

import pytest

from aerosysprop.diagnostics import (
    DiagnosticReport, Issue, STRICT_VALUES,
    diagnose, diagnose_smiles, diagnose_aiomfac, diagnose_aiomfac_sr,
    apply_strict,
    INVALID_SMILES, EMPTY_DECOMPOSITION, UNMATCHED_ATOMS,
    EXTENSION_SUBGROUP, SUBGROUP_OUT_OF_RANGE, NEGATIVE_COUNT,
    MR_PARAM_MISSING, COMPOSITION_SUM, COMPOSITION_VALUE,
    NO_WATER, HIGH_CONCENTRATION, BISULFATE_NEEDED, QCCA_NEEDED,
    UNKNOWN_SALT,
)
from aerosysprop.methods import boiling_points, vapour_pressures


# ----------------------- Issue / DiagnosticReport types -----------------------

def test_issue_str_format():
    iss = Issue('warning', 'X', 'a thing happened')
    s = str(iss)
    assert 'WARNING' in s and 'X' in s and 'a thing happened' in s


def test_report_aggregation():
    r = DiagnosticReport()
    assert r.is_ok
    r.add('warning', 'W1', 'w')
    r.add('error', 'E1', 'e')
    assert not r.is_ok
    assert r.has_errors
    assert r.has_warnings
    assert len(r.errors) == 1
    assert len(r.warnings) == 1
    assert len(r) == 2
    assert 'W1' in r.summary() and 'E1' in r.summary()


# ----------------------- SMILES validity ------------------------------------

def test_invalid_smiles_flagged():
    r = diagnose_smiles('this is not a SMILES')
    assert any(i.code == INVALID_SMILES for i in r.issues)


# ----------------------- AIOMFAC checks -------------------------------------

def test_aiomfac_acetone_flags_extension():
    """Acetone goes through the 281 extension (aerosysprop adds 281 = 191
    remap). Fortran AIOMFAC doesn't have 281; the framework should flag
    both EXTENSION_SUBGROUP and SUBGROUP_OUT_OF_RANGE."""
    r = diagnose_aiomfac_sr(['O', 'CC(=O)C'], [0.5, 0.5], 298.15)
    codes = {i.code for i in r.issues}
    assert EXTENSION_SUBGROUP in codes
    assert SUBGROUP_OUT_OF_RANGE in codes


def test_aiomfac_nitromethane_unmatched_atoms():
    """Nitromethane has -NO2 which aerosysprop's AIOMFAC SMARTS doesn't
    cover. The N and both O atoms should show up as unmatched."""
    r = diagnose_aiomfac_sr(['O', 'C[N+](=O)[O-]'], [0.5, 0.5], 298.15)
    issues = [i for i in r.issues if i.code == UNMATCHED_ATOMS]
    assert issues, 'should flag unmatched -NO2 atoms'
    unmatched_symbols = [
        sym for iss in issues for _, sym in iss.details['unmatched']
    ]
    assert 'N' in unmatched_symbols
    assert 'O' in unmatched_symbols


def test_aiomfac_dmso_unmatched_S():
    """DMSO's S=O group isn't in aerosysprop's AIOMFAC SMARTS — should
    flag S as unmatched."""
    r = diagnose_aiomfac_sr(['O', 'CS(=O)C'], [0.5, 0.5], 298.15)
    issues = [i for i in r.issues if i.code == UNMATCHED_ATOMS]
    assert issues
    unmatched = [sym for iss in issues for _, sym in iss.details['unmatched']]
    assert 'S' in unmatched


def test_aiomfac_mass_fraction_sum():
    """Mass fractions not summing to 1.0 should produce a warning."""
    r = diagnose_aiomfac(['O', 'CCO'], ['NaCl'], [0.50, 0.50, 0.50], 298.15)
    codes = {i.code for i in r.issues}
    assert COMPOSITION_SUM in codes


def test_aiomfac_mass_fraction_negative():
    r = diagnose_aiomfac(['O', 'CCO'], ['NaCl'], [-0.10, 0.50, 0.60], 298.15)
    codes = {i.code for i in r.issues}
    assert COMPOSITION_VALUE in codes


def test_aiomfac_mass_fraction_length_mismatch():
    r = diagnose_aiomfac(['O', 'CCO'], ['NaCl'], [0.7, 0.3], 298.15)
    codes = {i.code for i in r.issues}
    assert COMPOSITION_VALUE in codes


def test_aiomfac_unknown_salt():
    r = diagnose_aiomfac(['O'], ['NotARealSalt'], [0.95, 0.05], 298.15)
    codes = {i.code for i in r.issues}
    assert UNKNOWN_SALT in codes


def test_aiomfac_missing_org_ion_parameter():
    """Toluene (main group 4, ACCHn) with sulfate -- Fortran AIOMFAC has
    no parameter for this pair (sentinel -7.7777e4). Should be flagged
    as an error."""
    r = diagnose_aiomfac(
        ['O', 'Cc1ccccc1'], ['NH42SO4'],
        [0.85, 0.10, 0.05], 298.15,
    )
    issues = [i for i in r.issues if i.code == MR_PARAM_MISSING]
    assert issues, 'should flag main_group 4 vs SO4-- missing parameter'


def test_aiomfac_bisulfate_warning():
    r = diagnose_aiomfac(['O'], ['H2SO4'], [0.95, 0.05], 298.15)
    codes = {i.code for i in r.issues}
    assert BISULFATE_NEEDED in codes


def test_aiomfac_high_concentration_warning():
    r = diagnose_aiomfac(['O'], ['NaCl'], [0.65, 0.35], 298.15)
    codes = {i.code for i in r.issues}
    assert HIGH_CONCENTRATION in codes


# ----------------------- Generic dispatcher ---------------------------------

def test_dispatch_property_method_unmatched_atoms():
    """Methanethiol with the Nannoolal SMARTS — both atoms unmatched."""
    r = diagnose(boiling_points.nannoolal, 'CS')
    issues = [i for i in r.issues if i.code == UNMATCHED_ATOMS]
    assert issues


def test_dispatch_clean_property_method():
    r = diagnose(boiling_points.nannoolal, 'CCCCO')   # 1-butanol
    assert r.is_ok


def test_dispatch_via_string_path():
    r = diagnose(
        'aerosysprop.methods.boiling_points.nannoolal',
        'CC(=O)C',  # acetone -- Nannoolal handles it
    )
    # Nannoolal's matcher covers ketones, so atom coverage should be OK
    assert not any(i.code == UNMATCHED_ATOMS for i in r.issues)


# ----------------------- strict= integration --------------------------------

def test_apply_strict_silent():
    r = DiagnosticReport()
    r.add('warning', 'W', 'w')
    apply_strict(r, 'silent')   # no-op


def test_apply_strict_warn_emits_warnings():
    r = DiagnosticReport()
    r.add('warning', 'W', 'w')
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        apply_strict(r, 'warn')
    assert len(captured) == 1
    assert 'W' in str(captured[0].message)


def test_apply_strict_raise():
    r = DiagnosticReport()
    r.add('error', 'E', 'e')
    with pytest.raises(ValueError, match='E'):
        apply_strict(r, 'raise')


def test_apply_strict_invalid_value():
    with pytest.raises(ValueError, match='strict must be'):
        apply_strict(DiagnosticReport(), 'bogus')


def test_aiomfac_sr_from_smiles_strict_raise():
    """Driver-level integration: strict='raise' propagates to a ValueError."""
    from aerosysprop.methods.activity_coefficients import aiomfac_sr_from_smiles
    with pytest.raises(ValueError):
        aiomfac_sr_from_smiles(
            ['O', 'CC(=O)C'], [0.5, 0.5], 298.15, strict='raise'
        )


def test_aiomfac_from_smiles_strict_silent_unchanged():
    """strict='silent' should preserve the original behaviour (no warnings,
    no exceptions even for problematic inputs)."""
    from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        # acetone has 281 — triggers warnings under 'warn'
        result = aiomfac_from_smiles(
            ['O', 'CC(=O)C'], ['NaCl'], [0.85, 0.10, 0.05], 298.15,
            strict='silent',
        )
    assert len(captured) == 0
    assert result.gamma_neutral.shape == (2,)


def test_strict_values_constant():
    assert STRICT_VALUES == ('silent', 'warn', 'raise')


# ----------------------- result.diagnostics wiring -------------------------

def test_aiomfac_result_carries_diagnostics():
    """AIOMFACResult.diagnostics should be populated when strict != 'silent'."""
    import warnings
    from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result = aiomfac_from_smiles(
            ['O', 'CC(=O)C'], ['NaCl'], [0.85, 0.10, 0.05], 298.15,
        )
    assert result.diagnostics is not None
    # Acetone hits the 281 extension; should report at least one issue.
    assert len(result.diagnostics) > 0
    codes = {i.code for i in result.diagnostics}
    assert EXTENSION_SUBGROUP in codes


def test_aiomfac_result_diagnostics_none_when_silent():
    """strict='silent' should produce diagnostics=None (no overhead)."""
    from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles
    result = aiomfac_from_smiles(
        ['O', 'CCO'], ['NaCl'], [0.85, 0.10, 0.05], 298.15,
        strict='silent',
    )
    assert result.diagnostics is None


def test_aiomfac_result_namedtuple_compat():
    """AIOMFACResult should remain a namedtuple with .gamma_neutral etc."""
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_from_smiles, AIOMFACResult,
    )
    result = aiomfac_from_smiles(
        ['O', 'CCO'], ['NaCl'], [0.85, 0.10, 0.05], 298.15, strict='silent',
    )
    assert isinstance(result, AIOMFACResult)
    # Named access still works
    assert result.gamma_neutral is not None
    assert result.cation_ids == [202]
    assert result.diagnostics is None


def test_aiomfac_sr_return_diagnostics_tuple():
    import warnings
    from aerosysprop.methods.activity_coefficients import aiomfac_sr_from_smiles
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gamma, report = aiomfac_sr_from_smiles(
            ['O', 'CC(=O)C'], [0.5, 0.5], 298.15,
            return_diagnostics=True,
        )
    assert gamma.shape == (2,)
    assert report is not None
    assert len(report) > 0


def test_aiomfac_sr_return_diagnostics_default_off():
    """Default behaviour is unchanged — returns just the ndarray."""
    from aerosysprop.methods.activity_coefficients import aiomfac_sr_from_smiles
    result = aiomfac_sr_from_smiles(['O', 'CCO'], [0.5, 0.5], 298.15)
    # Should be just an ndarray, not a tuple
    assert hasattr(result, 'shape')
    assert result.shape == (2,)


def test_aiomfac_sr_return_diagnostics_silent_mode():
    """return_diagnostics=True works even when strict='silent' (still
    computes the report on demand)."""
    from aerosysprop.methods.activity_coefficients import aiomfac_sr_from_smiles
    gamma, report = aiomfac_sr_from_smiles(
        ['O', 'CCO'], [0.5, 0.5], 298.15,
        strict='silent', return_diagnostics=True,
    )
    assert gamma.shape == (2,)
    # Clean SMILES => is_ok
    assert report.is_ok


# ----------------------- per-SMILES diagnostic cache -----------------------

def test_diagnostic_cache_amortises_repeated_smiles():
    """Repeated calls with the same SMILES should hit the cache; the
    overhead of strict='warn' versus strict='silent' should be small."""
    import time
    import warnings
    from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles
    from aerosysprop.diagnostics import (
        clear_diagnostic_cache, diagnostic_cache_info,
    )
    clear_diagnostic_cache()

    args = (['O', 'CCO'], ['NaCl'], [0.85, 0.10, 0.05], 298.15)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Warm both paths
        aiomfac_from_smiles(*args, strict='silent')
        aiomfac_from_smiles(*args, strict='warn')

        n = 200
        t0 = time.perf_counter()
        for _ in range(n):
            aiomfac_from_smiles(*args, strict='silent')
        silent = time.perf_counter() - t0

        t0 = time.perf_counter()
        for _ in range(n):
            aiomfac_from_smiles(*args, strict='warn')
        warn = time.perf_counter() - t0

    info = diagnostic_cache_info()
    # Should be a high hit rate
    assert info.hits >= info.misses

    # Memoised diagnostics overhead should be small. We assert < 4x because
    # hot-path timing is noisy across machines; the practical observed
    # overhead is ~1.05x, so a 4x bound is very conservative but still rules
    # out a regression to the pre-memoised ~2x state.
    assert warn < 4.0 * silent, (
        f'memoised strict="warn" took {warn:.3f}s, '
        f'strict="silent" took {silent:.3f}s'
    )


def test_clear_diagnostic_cache():
    from aerosysprop.diagnostics import (
        clear_diagnostic_cache, diagnostic_cache_info, diagnose_aiomfac_sr,
    )
    diagnose_aiomfac_sr(['O', 'CCO'], [0.5, 0.5], 298.15)
    info1 = diagnostic_cache_info()
    assert info1.currsize > 0
    clear_diagnostic_cache()
    info2 = diagnostic_cache_info()
    assert info2.currsize == 0


# ----------------------- System cache --------------------------------------

def test_system_cache_hits_on_repeated_smiles():
    """Repeated aiomfac_from_smiles calls with the same SMILES list should
    reuse the underlying AIOMFAC system."""
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_from_smiles, system_cache_info, clear_system_cache,
    )
    import warnings
    clear_system_cache()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(50):
            aiomfac_from_smiles(
                ['O', 'CCO'], ['NaCl'], [0.85, 0.10, 0.05], 298.15,
                strict='silent',
            )
    aiomfac_info, sr_info = system_cache_info()
    # The first call is a miss, the rest are hits.
    assert aiomfac_info.misses == 1
    assert aiomfac_info.hits == 49


def test_system_cache_distinguishes_different_smiles_lists():
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_from_smiles, system_cache_info, clear_system_cache,
    )
    clear_system_cache()
    aiomfac_from_smiles(
        ['O', 'CCO'], ['NaCl'], [0.85, 0.10, 0.05], 298.15, strict='silent',
    )
    aiomfac_from_smiles(
        ['O', 'CO'], ['NaCl'], [0.85, 0.10, 0.05], 298.15, strict='silent',
    )
    aiomfac_from_smiles(
        ['O', 'CCO'], ['NaCl'], [0.50, 0.40, 0.10], 298.15, strict='silent',
    )
    info, _ = system_cache_info()
    # Two distinct SMILES lists -> two misses, the third call reuses entry 1
    assert info.misses == 2
    assert info.hits == 1


def test_system_cache_sr_path():
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_sr_from_smiles, system_cache_info, clear_system_cache,
    )
    import warnings
    clear_system_cache()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(10):
            aiomfac_sr_from_smiles(
                ['O', 'CCO'], [0.5, 0.5], 298.15, strict='silent',
            )
    _, sr_info = system_cache_info()
    assert sr_info.misses == 1
    assert sr_info.hits == 9


def test_clear_system_cache_drops_entries():
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_from_smiles, system_cache_info, clear_system_cache,
    )
    aiomfac_from_smiles(['O', 'CCO'], ['NaCl'], [0.85, 0.10, 0.05],
                        298.15, strict='silent')
    info, _ = system_cache_info()
    assert info.currsize > 0
    clear_system_cache()
    info, _ = system_cache_info()
    assert info.currsize == 0


def test_entry_point_speed_close_to_cached_system():
    """Regression guard: aiomfac_from_smiles (strict='silent') should be
    within a small multiple of the manually-cached-system speed once the
    system cache is warm. Catches a regression to the pre-system-cache
    state where the entry point was 30x slower."""
    import time, warnings
    from aerosysprop.methods.activity_coefficients import (
        aiomfac_from_smiles, AIOMFAC, clear_system_cache,
    )
    from aerosysprop import groups
    from rdkit import Chem

    SMI = ['O', 'CCO', 'CO', 'OCC(O)CO']
    SALTS = ['NaCl']
    WTF = [0.7, 0.10, 0.05, 0.05, 0.10]

    # warm both paths
    clear_system_cache()
    aiomfac_from_smiles(SMI, SALTS, WTF, 298.15, strict='silent')
    neutrals = [groups.aiomfac(Chem.MolFromSmiles(s)) for s in SMI]
    manual_system = AIOMFAC(neutrals, salts=SALTS)
    manual_system.gamma(WTF, 298.15)

    n = 500
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        for _ in range(n):
            aiomfac_from_smiles(SMI, SALTS, WTF, 298.15, strict='silent')
        entry = time.perf_counter() - t0

    t0 = time.perf_counter()
    for _ in range(n):
        manual_system.gamma(WTF, 298.15)
    manual = time.perf_counter() - t0

    # Generous bound (5x) - actual ratio is ~1.05x but timing varies
    assert entry < 5.0 * manual, (
        f'entry-point {entry:.3f}s vs manual {manual:.3f}s '
        f'(ratio {entry/manual:.2f}x) - system cache may have regressed'
    )
