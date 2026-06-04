"""Documentation / user guide page.

Explains the cross-session database file, how diagnostics are reported, how
to drive AeroSysProp from your own Python code (rather than the app), and a
few other things worth knowing. Prediction-method details are intentionally
left to a separate page.
"""
import streamlit as st

from aerosysprop._ui import ensure_session_state, sidebar_session_loader
from aerosysprop.mixture_methods import MIXTURE_METHODS

ensure_session_state()

with st.sidebar:
    st.title('Documentation')
    st.divider()
    sidebar_session_loader()

st.header('Documentation')
st.markdown(
    'How AeroSysProp works as an application and as a Python library. '
    'This page covers the session database, the input diagnostics, and how '
    'to call the same calculations from your own code. '
    '*The prediction methods themselves are documented on a separate page.*'
)

(tab_overview, tab_db, tab_diag, tab_mixture, tab_python,
 tab_notes) = st.tabs([
    'Overview',
    'Session database',
    'Diagnostics',
    'Mixture methods',
    'Use from Python',
    'Other notes',
])


# ---------------------------------------------------------------------------
# Overview
# ---------------------------------------------------------------------------
with tab_overview:
    st.subheader('What the app does')
    st.markdown(
        """
AeroSysProp predicts thermophysical properties of organic molecules
relevant to atmospheric aerosol research, and activity coefficients of
mixtures. The sidebar navigates three working areas plus this guide:

- **Single-component properties** — predicts each compound independently
  (boiling points, vapour pressures, critical properties, densities) for
  a list of SMILES.
- **Mixture properties** — treats the loaded compounds as one coupled
  mixture. A family of methods (currently AIOMFAC activity coefficients;
  see the **Mixture methods** tab) where every species' result depends on
  the whole composition.
- **External services** — hands your compounds to third-party tools
  (e.g. OPERA) for comparison.

Everything the app does is also available as a plain Python library, so
you can embed the same calculations in scripts, notebooks, or numerical
simulations — see the **Use from Python** tab.
"""
    )


# ---------------------------------------------------------------------------
# Session database
# ---------------------------------------------------------------------------
with tab_db:
    st.subheader('The session database (.db file)')
    st.markdown(
        """
AeroSysProp can persist your work to a single **SQLite `.db` file** so it
survives between sessions. There is nothing running in the background —
the file *is* your saved state, and you carry it with you.

**Saving.** On the single-component and mixture pages, click *Save to
session database*. The app writes your results to a `.db` file and offers
it for download.

**Reloading.** Next time, upload that `.db` in the sidebar (*Session
database → Load existing session*). It is shared across every page for
the rest of the session: reload your compounds on the input tabs, and any
stored results are available again.

**What's stored.** Four tables:

| Table | Contents |
|---|---|
| `compounds` | each compound's name + SMILES (de-duplicated by SMILES) |
| `results` | single-component predictions, long-form: one row per compound × method × temperature × value |
| `mixtures` | one row per saved AIOMFAC mixture calculation (name, temperature, basis, salts) |
| `mixture_species` | per-species activity coefficients for each mixture (neutral / cation / anion / mean-ionic) |

Because results are stored *long-form* (one row per value, not one column
per method), you can keep adding methods to an existing `.db` over time
without changing its structure, and query any subset later.

**It's just SQLite.** The `.db` is a standard SQLite database — open it
with any SQLite browser, `sqlite3` on the command line, or pandas:
"""
    )
    st.code(
        """import sqlite3, pandas as pd

con = sqlite3.connect('aerosysprop_session.db')
compounds = pd.read_sql('SELECT * FROM compounds', con)
results   = pd.read_sql('SELECT * FROM results', con)
mixtures  = pd.read_sql('SELECT * FROM mixtures', con)
species   = pd.read_sql('SELECT * FROM mixture_species', con)
""",
        language='python',
    )
    st.markdown(
        'Or use the helper functions in `aerosysprop.db` '
        '(`list_compounds`, `list_mixtures`, `load_mixture_species`, '
        '`save_results`, `save_mixture_result`).'
    )


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
with tab_diag:
    st.subheader('Input diagnostics')
    st.markdown(
        """
Several failure modes in property/activity prediction are *silent*: the
SMARTS matcher may not recognise a functional group and quietly drop it,
a salt + organic combination may have no fitted interaction parameter, or
a composition may not sum to 1. AeroSysProp screens every input and
surfaces these as **diagnostics** so you are never given a
silently-wrong number without a flag.

In the app, diagnostics appear **inline with the results** (the mixture
page lists them in an expander). Each issue has a **severity** and a
**code**:

- **warning** — the calculation ran, but something may make it less
  reliable or differ from the reference implementation.
- **error** — the input is malformed or unsupported; treat the output
  with caution (it will not match the Fortran AIOMFAC reference).

**General codes** (apply to any SMILES / composition input)

| Code | Severity | Meaning |
|---|---|---|
| `INVALID_SMILES` | error | RDKit could not parse the SMILES. |
| `UNMATCHED_ATOMS` | warning | Some atoms have no SMARTS match — the method behaves as if those groups were absent (silent truncation). |
| `COMPOSITION_SUM` | warning | Mass / mole fractions don't sum to ≈ 1. |
| `COMPOSITION_VALUE` | error | A fraction is outside [0, 1] or the vector length is wrong. |
| `NO_WATER` | warning | An electrolyte mixture without water as a component. |

**Strict modes (in code).** When you call the methods from Python, a
`strict=` argument controls how issues surface:

- `'warn'` (default) — emit a Python warning and continue;
- `'raise'` — raise `ValueError` on any issue;
- `'silent'` — skip the pre-flight check entirely.

You can also pre-screen inputs without computing anything — see the
**Use from Python** tab.
"""
    )

    # Method-specific diagnostic caveats, pulled from the mixture-method
    # registry so each method documents its own codes.
    method_diags = [m for m in MIXTURE_METHODS if m.diagnostics_notes]
    if method_diags:
        st.markdown('**Method-specific notes**')
        for m in method_diags:
            st.markdown(f'*{m.name}*')
            for note in m.diagnostics_notes:
                st.markdown(f'- {note}')


# ---------------------------------------------------------------------------
# Mixture methods (registry-driven — grows as methods are added)
# ---------------------------------------------------------------------------
with tab_mixture:
    st.subheader('Mixture-property methods')
    st.markdown(
        'The **Mixture properties** page is a family of methods that act on a '
        'whole mixture rather than each compound independently. Each method '
        'below documents what it does, when to use it, and how to call it '
        'from Python. (Detailed equations and validation live in each '
        "method's own reference document.)"
    )
    for m in MIXTURE_METHODS:
        badge = '' if m.status == 'available' else '  *(coming soon)*'
        st.markdown(f'### {m.name}{badge}')
        st.markdown(m.summary)
        if m.when_to_use:
            st.markdown(f'**When to use.** {m.when_to_use}')
        if m.python_example:
            st.markdown('**From Python**')
            st.code(m.python_example, language='python')
        if m.reference_doc:
            st.caption(f'Detailed reference: `{m.reference_doc}`')
        st.divider()


# ---------------------------------------------------------------------------
# Use from Python
# ---------------------------------------------------------------------------
with tab_python:
    st.subheader('Driving AeroSysProp from your own code')
    st.markdown(
        'Everything in the app is a thin wrapper over the `aerosysprop` '
        'package. Install it once (`pip install -e .` from the repo) and '
        'import what you need.'
    )

    st.markdown('**Single-component properties**')
    st.code(
        """from rdkit import Chem
from aerosysprop.methods import boiling_points, vapour_pressures

mol = Chem.MolFromSmiles('CC(=O)O')          # acetic acid
bp  = boiling_points.nannoolal(mol)          # K
vp  = vapour_pressures.simpol(mol, 298.15)   # log10(atm)
""",
        language='python',
    )

    st.markdown('**Batch prediction** (the engine behind the single-component page)')
    st.code(
        """from rdkit import Chem
from aerosysprop.batch import run_all

mols = [('acetic_acid', Chem.MolFromSmiles('CC(=O)O')),
        ('ethanol',     Chem.MolFromSmiles('CCO'))]

results = run_all(
    mols,
    methods=['bp_nannoolal', 'vp_simpol', 'vp_nannoolal'],
    temperature=298.15,
    # one boiling-point source -> usual column names;
    # several -> '__bp<Source>'-suffixed columns for a sensitivity sweep
    bp_sources=['bp_nannoolal', 'bp_stein'],
)
""",
        language='python',
    )

    st.markdown(
        '**Mixture properties** — each mixture-property method documents its '
        'own Python entry points in the **Mixture methods** tab.'
    )

    st.markdown('**Pre-flight diagnostics without computing**')
    st.code(
        """from aerosysprop.diagnostics import diagnose_aiomfac

report = diagnose_aiomfac(['O', 'Cc1ccccc1'], ['NH42SO4'],
                          [0.85, 0.10, 0.05], 298.15)
report.is_ok            # False
report.has_errors       # True (toluene main group vs sulfate has no parameter)
for issue in report:
    print(issue.severity, issue.code, issue.message)
""",
        language='python',
    )

    st.markdown('**SMARTS-hit fingerprint matrix** (batch SMILES inspection)')
    st.code(
        """from aerosysprop.fingerprints import smarts_hits

df = smarts_hits(['O', 'CCO', 'OCC(O)CO'], method='aiomfac',
                 drop_zero_columns=True)   # pandas DataFrame: SMILES x subgroup
""",
        language='python',
    )

    st.markdown('**Saving / loading the session database from code**')
    st.code(
        """from aerosysprop import db

db.save_results('session.db', results, temperature=298.15)
db.list_compounds('session.db')
db.list_mixtures('session.db')
""",
        language='python',
    )

    st.info(
        'For numerical simulations that call AIOMFAC thousands–millions of '
        'times, the SMILES entry points cache the constructed system and the '
        'diagnostics, so repeated calls with the same compounds run at '
        '~0.1 ms each. See the performance notes in '
        '`aerosysprop/methods/activity_coefficients.md`.',
        icon='⚡',
    )


# ---------------------------------------------------------------------------
# Other notes
# ---------------------------------------------------------------------------
with tab_notes:
    st.subheader('Other things to note')
    st.markdown(
        """
**Loading compounds.** Each working page accepts SMILES three ways: typed
in (optionally `name SMILES` per line), uploaded as a file
(`.txt` `.csv` `.tsv` `.prop` one-or-two-column, `.sdf`, or a `.db`
session), or reloaded from a session database you've loaded in the
sidebar.

**Mixtures vs single components.** The single-component page predicts each
molecule independently. The mixture page is fundamentally different: it
treats the loaded compounds as *one mixture*, so every species' activity
coefficient depends on the whole composition. Water is the SMILES `O` and,
by AIOMFAC convention, is normally the first component.

**Salts.** On the mixture page, inorganic salts are chosen from the
AIOMFAC catalog (NaCl, (NH₄)₂SO₄, NH₄NO₃, …). They dissociate into ions;
activity coefficients are returned per ion plus a mean ionic value per
salt. Sulfuric acid / bisulfate systems are not yet supported (flagged by
the `BISULFATE_NEEDED` diagnostic).

**Reproducibility.** Activity-coefficient outputs are validated against
the public Fortran AIOMFAC reference to ~5 × 10⁻⁵ relative error across
thousands of mixtures. Where the Python port deviates (e.g. unsupported
subgroups, bisulfate), the diagnostics tell you.

**Methods documentation.** Details of the individual prediction methods
(Nannoolal, SIMPOL, EVAPORATION, AIOMFAC SR / MR / LR, …) — their
equations, parameter sources, and validation — are documented separately.
"""
    )
