"""'Add a method' page — a contributor guide for extending AeroSysProp with
your own single-component property methods or mixture-property methods.

Explains the repo layout, the registration points, whether the app updates
automatically, and a checklist for each kind of method.
"""
import streamlit as st

from aerosysprop._ui import ensure_session_state, sidebar_session_loader

ensure_session_state()

with st.sidebar:
    st.title('Add a method')
    st.divider()
    sidebar_session_loader()

st.header('Add a method')
st.markdown(
    'How to extend AeroSysProp with your own **single-component property** '
    'method or **mixture-property** method — from repo structure to '
    'registration to a checklist. This is contributor-level guidance; for '
    'using the existing methods see the **Documentation** page.'
)

tab_how, tab_repo, tab_prop, tab_mix, tab_update, tab_check = st.tabs([
    'How discovery works',
    'Repo structure',
    'Add a property method',
    'Add a mixture method',
    'Does the app auto-update?',
    'Checklist',
])


# ---------------------------------------------------------------------------
# How discovery works
# ---------------------------------------------------------------------------
with tab_how:
    st.subheader('How the app discovers methods')
    st.markdown(
        """
The app does **not** hardcode its method lists in the pages. Each kind of
method has a **registry** — a plain Python list, in code, where every
method is described by one entry (its id, display label, units, group,
output columns, …). The pages read that list when they run and build
themselves from it. "Registering" a method just means adding an entry to
the list.

There are two registries:

- **Single-component methods** → the `PROPERTY_METHODS` list in
  `aerosysprop/batch.py`. Each entry is a `PropertyMethod(...)`. The
  single-component page loops over this list and draws one checkbox per
  method, grouped by the entry's `group` field. (For backward
  compatibility the older lookup tables — `METHOD_META`, `METHOD_GROUPS`,
  `ALL_METHODS`, `_OUTPUT_KEYS` — are now *generated from* this one list,
  so you only edit the list.)
- **Mixture methods** → the `MIXTURE_METHODS` list in
  `aerosysprop/mixture_methods.py`. Each entry is a `MixtureMethod(...)`.
  The mixture page builds its method picker from it, and the
  Documentation page builds its "Mixture methods" tab from it.

So **adding a method is a code change**: you implement the calculation and
add one entry to the relevant registry. You do *not* wire it into the UI
by hand — once the entry is there, the page shows it automatically the
next time the script runs.
"""
    )


# ---------------------------------------------------------------------------
# Repo structure
# ---------------------------------------------------------------------------
with tab_repo:
    st.subheader('The parts of the repo you will touch')
    st.code(
        """aerosysprop/
├── groups.py                 # SMARTS group counting / aggregation helpers
├── data/                     # parameter + SMARTS files, loaded at import
│   └── __init__.py           #   exposes them as data.SOMETHING
├── methods/
│   ├── boiling_points.py     # single-component property methods live here
│   ├── vapour_pressures.py   #   (one module per property family)
│   ├── critical_properties.py
│   ├── liquid_densities.py
│   └── activity_coefficients.py   # AIOMFAC (a mixture method)
├── batch.py                  # registry + worker for single-component methods
├── mixture_methods.py        # registry for mixture-property methods
├── diagnostics.py            # optional input pre-flight checks
├── views/                    # the Streamlit pages (read the registries)
│   ├── single_component.py
│   └── mixture.py
└── tests/                    # pytest + fixture generators
""",
        language='text',
    )
    st.markdown(
        'New parameter tables or SMARTS go in `aerosysprop/data/` and are '
        'loaded in `data/__init__.py` so they are importable as '
        '`aerosysprop.data.YOUR_TABLE`.'
    )


# ---------------------------------------------------------------------------
# Add a property method
# ---------------------------------------------------------------------------
with tab_prop:
    st.subheader('Adding a single-component property method')
    st.markdown(
        'A property method takes one RDKit molecule (plus optional '
        'temperature / boiling point) and returns a scalar. Example: a new '
        'boiling-point method `bp_mymethod`.'
    )

    st.markdown('**1. Implement the function** in the relevant `methods/` module')
    st.code(
        """# aerosysprop/methods/boiling_points.py
from aerosysprop import data, groups
from aerosysprop.diagnostics import with_diagnostics

@with_diagnostics                      # optional: adds the strict= argument
def my_method(compound):
    counts = groups.nannoolal_primary(compound)
    return groups.aggregate_matches(counts, data.MY_TABLE) + 198.2
""",
        language='python',
    )

    st.markdown('**2. Register it** — add one `PropertyMethod` entry to '
                '`PROPERTY_METHODS` in `aerosysprop/batch.py`')
    st.code(
        """# aerosysprop/batch.py — append to the PROPERTY_METHODS registry.
# All the lookup tables (METHOD_META, ALL_METHODS, METHOD_GROUPS,
# _OUTPUT_KEYS, …) are derived from this list automatically.
PropertyMethod(
    id='bp_mymethod',
    label='My Method (2026)',
    unit='K',
    group='Boiling Point',              # appears under this sidebar group
    output_keys=('bp_mymethod_K',),     # result column name(s)
    # deps=frozenset({BP_ANY}),         # if it needs *a* boiling point
    # bp_dependent=True,                # if its output changes with the BP source
    # bp_short='MyBP',                  # only for boiling-point methods used as a source
)
""",
        language='python',
    )

    st.markdown('**3. Wire the calculation into the worker** `_run_one` in `batch.py`')
    st.code(
        """# inside _run_one(...), alongside the other boiling points:
if 'bp_mymethod' in methods:
    out['bp_mymethod_K'] = _call(boiling_points.my_method, mol)
""",
        language='python',
    )
    st.markdown(
        'If the method **depends on a boiling point** (like a vapour '
        'pressure), set `deps=frozenset({BP_ANY})`, add its keys to '
        '`_BP_DEPENDENT_KEYS`, and compute it inside the per-boiling-point '
        'loop of `_run_one` so it participates in the multi-source feature.'
    )

    st.markdown('**4. The UI updates itself.** The single-component page reads '
                '`METHOD_GROUPS` / `METHOD_META`, so a checkbox for your '
                'method now appears under its group automatically.')

    st.markdown('**5. Test it.** Add a case to `aerosysprop/tests/` and, if '
                'the method has reference values, a fixture.')


# ---------------------------------------------------------------------------
# Add a mixture method
# ---------------------------------------------------------------------------
with tab_mix:
    st.subheader('Adding a mixture-property method')
    st.markdown(
        'A mixture method acts on the whole mixture (its result for each '
        'species depends on the full composition). AIOMFAC is the first one; '
        'here is how to add another.'
    )

    st.markdown('**1. Implement the calculation** (e.g. a new `methods/` module)')
    st.code(
        """# aerosysprop/methods/my_mixture.py
def my_mixture_from_smiles(smiles, mass_fractions, T_K):
    ...
    return result   # however you choose to represent per-species output
""",
        language='python',
    )

    st.markdown('**2. Add a registry entry** in `aerosysprop/mixture_methods.py`')
    st.code(
        """MIXTURE_METHODS.append(MixtureMethod(
    id='mymodel',
    name='My activity model',
    status='available',                 # or 'coming_soon'
    summary='One or two sentences on what it does.',
    when_to_use='When this method applies.',
    python_example='from aerosysprop.methods.my_mixture import ...',
    reference_doc='aerosysprop/methods/my_mixture.md',
    diagnostics_notes=['`SOME_CODE` — what it means for this method.'],
))
""",
        language='python',
    )
    st.markdown(
        'This alone makes your method appear in the mixture-page **picker** '
        'and gain its own section in the Documentation **Mixture methods** '
        'tab.'
    )

    st.markdown('**3. Dispatch the computation** in `aerosysprop/views/mixture.py`')
    st.markdown(
        'The mixture page resolves the picked method to a registry entry '
        '(`method = method_by_name(...)`). The compute + results block is '
        'currently AIOMFAC-specific; branch on `method.id` to call your '
        'calculation and map its output into the per-species rows the page '
        'displays and saves:'
    )
    st.code(
        """if method.id == 'aiomfac':
    result = aiomfac_from_smiles(...)
elif method.id == 'mymodel':
    result = my_mixture_from_smiles(...)
# then build the species rows (neutral / cation / anion / mean) for the
# results table + session-database save.
""",
        language='python',
    )

    st.markdown('**4. (Optional) diagnostics.** Add input checks to '
                '`aerosysprop/diagnostics.py` if your method has its own '
                'failure modes; list the codes in the registry entry\'s '
                '`diagnostics_notes`.')

    st.markdown('**5. Test it** against a reference, ideally with a fixture '
                'generator under `aerosysprop/tests/`.')


# ---------------------------------------------------------------------------
# Does the app auto-update?
# ---------------------------------------------------------------------------
with tab_update:
    st.subheader('Does the Streamlit app update automatically?')
    st.markdown(
        """
**Short answer: the UI updates itself once the code is changed — but the
code change is on you.** There is no separate "register with the running
app" step.

- The pages are **built from the registries at runtime**, so as soon as
  your method is in `METHOD_META` / `METHOD_GROUPS` (property) or
  `MIXTURE_METHODS` (mixture), the checkbox / picker entry / docs section
  appear with no edits to the page code.
- Streamlit **re-runs the script on every file save**, so during
  development you usually see your change immediately after saving.
- **Restart the server** (stop and `streamlit run aerosysprop/app.py`
  again) when you have added a *new module*, new data loaded at import
  time, or changed imports — import-time changes are not always picked up
  by a hot re-run.
- If you installed the package (`pip install -e .`), an **editable
  install** points at your source tree, so your edits are seen without
  reinstalling. A non-editable install would need reinstalling.

So: edit + register + (save or restart) → the app reflects it. It will
not appear without the code change.
"""
    )


# ---------------------------------------------------------------------------
# Checklist
# ---------------------------------------------------------------------------
with tab_check:
    st.subheader('Checklists')

    st.markdown('**Single-component property method**')
    st.markdown(
        """
- [ ] Implement the function in the relevant `methods/` module (RDKit mol in, scalar out)
- [ ] Add any new SMARTS / parameter files to `data/` and load them in `data/__init__.py`
- [ ] Add one `PropertyMethod(...)` entry to `PROPERTY_METHODS` in `batch.py` (id, label, unit, group, output_keys, and `deps` / `bp_dependent` if it needs a boiling point)
- [ ] Wire the calculation into `_run_one` in `batch.py`
- [ ] For boiling-point-dependent methods, compute inside the per-BP loop of `_run_one` (the `bp_dependent=True` flag handles the column suffixing)
- [ ] Add a test (and a fixture if there are reference values)
- [ ] Save / restart the app and confirm the checkbox appears and runs
"""
    )

    st.markdown('**Mixture-property method**')
    st.markdown(
        """
- [ ] Implement the calculation (new `methods/` module or function)
- [ ] Add a `MixtureMethod(...)` entry to `MIXTURE_METHODS`
- [ ] Branch on `method.id` in `views/mixture.py` to dispatch the compute and build the per-species result rows
- [ ] (Optional) add input checks to `diagnostics.py`; list codes in `diagnostics_notes`
- [ ] (Optional) write a `<method>.md` reference and point `reference_doc` at it
- [ ] Add a test / fixture against a reference implementation
- [ ] Save / restart the app and confirm the method appears in the picker and the Mixture methods docs tab
"""
    )

    st.info(
        'Run the test suite with `pytest aerosysprop/tests/` before opening '
        'a pull request. New methods that reproduce a published / reference '
        'implementation should ship with a fixture validating them to a '
        'stated tolerance, matching how the existing methods are tested.',
        icon='✅',
    )
