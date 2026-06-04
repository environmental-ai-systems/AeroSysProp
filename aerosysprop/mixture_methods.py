"""
aerosysprop.mixture_methods
===========================
Registry of mixture-property methods.

"Mixture properties" is a family of methods that act on a whole mixture
rather than each component independently. AIOMFAC activity coefficients are
the first member; more (other activity-coefficient models, or other mixture
properties entirely) can be added by appending a ``MixtureMethod`` entry to
``MIXTURE_METHODS``.

Both the Streamlit mixture page and the in-app Documentation page read from
this registry, so adding a method here makes it appear in the method picker
and gain its own documentation section automatically — no edits to the UI or
docs pages required.

Detailed method-level documentation (equations, parameter provenance,
validation) lives separately, referenced via ``reference_doc``.
"""
from dataclasses import dataclass, field
from typing import List


@dataclass
class MixtureMethod:
    """Metadata + self-documentation for one mixture-property method."""
    id: str                       # short stable identifier, e.g. 'aiomfac'
    name: str                     # display name for the picker / docs
    status: str                   # 'available' | 'coming_soon'
    summary: str                  # one or two sentence markdown summary
    when_to_use: str              # markdown: when this method applies
    python_example: str           # code snippet (str) for the Documentation page
    reference_doc: str = ''       # path/link to detailed method documentation
    diagnostics_notes: List[str] = field(default_factory=list)
    # markdown bullet strings describing method-specific diagnostic caveats


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

MIXTURE_METHODS: List[MixtureMethod] = [
    MixtureMethod(
        id='aiomfac',
        name='AIOMFAC activity coefficients',
        status='available',
        summary=(
            'Activity coefficients of organic + water + electrolyte mixtures '
            'via AIOMFAC: short-range (modified UNIFAC) + middle-range '
            '(organic-ion / ion-ion) + long-range (Debye–Hückel). Validated '
            'against the public Fortran AIOMFAC reference.'
        ),
        when_to_use=(
            'Use for the non-ideal mixing of organics, water and common '
            'atmospheric inorganic salts (NaCl, (NH₄)₂SO₄, NH₄NO₃, …). '
            'Returns per-species activity coefficients and a mean ionic value '
            'per salt for a composition you specify.'
        ),
        python_example=(
            "from aerosysprop.methods.activity_coefficients import (\n"
            "    aiomfac_sr_from_smiles,   # non-electrolyte (organic + water)\n"
            "    aiomfac_from_smiles,      # full SR + MR + LR with salts\n"
            ")\n\n"
            "# Non-electrolyte: gamma of water + ethanol at 50:50 mole fraction\n"
            "g = aiomfac_sr_from_smiles(['O', 'CCO'], [0.5, 0.5], 298.15)\n\n"
            "# Electrolyte: water + ethanol + NaCl from mass fractions\n"
            "res = aiomfac_from_smiles(['O', 'CCO'], ['NaCl'],\n"
            "                          [0.75, 0.20, 0.05], 298.15)\n"
            "res.gamma_neutral   # mole-fraction basis (water, ethanol)\n"
            "res.gamma_cation    # molality basis (Na+)\n"
            "res.gamma_anion     # molality basis (Cl-)\n"
            "res.gamma_mean      # mean ionic gamma_pm per salt\n"
            "res.diagnostics     # the same report shown in the app\n"
        ),
        reference_doc='aerosysprop/methods/activity_coefficients.md',
        diagnostics_notes=[
            '`MR_PARAM_MISSING` — an organic main group has no fitted '
            'interaction parameter against one of the ions in the mixture.',
            '`EXTENSION_SUBGROUP` — an aerosysprop-only subgroup (281–285) '
            'that the reference Fortran AIOMFAC drops.',
            '`BISULFATE_NEEDED` — H₂SO₄ present; the bisulfate equilibrium '
            '(HSO₄⁻ ⇌ H⁺ + SO₄²⁻) is not yet modelled, so results will '
            'differ from the reference.',
            '`HIGH_CONCENTRATION` — a salt above ~30 % mass fraction, outside '
            'the well-fitted range.',
        ],
    ),
    # ---- add further mixture-property methods below ----
    # MixtureMethod(
    #     id='...', name='...', status='coming_soon', summary='...',
    #     when_to_use='...', python_example='...', reference_doc='...',
    # ),
]


def available_methods() -> List[MixtureMethod]:
    """Return the mixture methods that are ready to use."""
    return [m for m in MIXTURE_METHODS if m.status == 'available']


def method_by_id(method_id: str) -> MixtureMethod:
    """Look up a registry entry by id."""
    for m in MIXTURE_METHODS:
        if m.id == method_id:
            return m
    raise KeyError(f'unknown mixture method {method_id!r}')


def method_by_name(name: str) -> MixtureMethod:
    """Look up a registry entry by display name."""
    for m in MIXTURE_METHODS:
        if m.name == name:
            return m
    raise KeyError(f'unknown mixture method name {name!r}')
