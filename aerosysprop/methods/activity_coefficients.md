# Activity coefficients — AIOMFAC short-range (SR) implementation

This document describes the Python implementation of the AIOMFAC short-range
activity-coefficient model in AeroSysProp. The implementation is a faithful
port of the Fortran reference code released alongside Zuend *et al.* (2008,
2011) and maintained at
[github.com/andizuend/AIOMFAC](https://github.com/andizuend/AIOMFAC).

For atmospheric scientists who simply want activity coefficients, see
[Quick start](#quick-start). Non-electrolyte mixtures use the SR part only;
electrolyte (salt-containing) mixtures additionally use the middle-range
(MR) and long-range (LR, Debye–Hückel) parts.

## Table of contents

- [Scope](#scope)
- [Quick start](#quick-start)
- [API reference](#api-reference)
- [Theory](#theory)
- [Compound coverage](#compound-coverage)
- [Validation methodology](#validation-methodology)
- [Data files](#data-files)
- [Diagnostics](#diagnostics)
- [Known limitations](#known-limitations)
- [Data-fix log](#data-fix-log)
- [Citations](#citations)
- [Roadmap](#roadmap)

---

## Scope

AIOMFAC = **A**erosol **I**norganic–**O**rganic **M**ixtures **F**unctional
groups **A**ctivity **C**oefficients. The full model has three parts:

| Part | What it captures | Status |
|---|---|---|
| **SR** (short range, modified UNIFAC) | Entropic + energetic mixing of neutral functional groups | **✓ implemented** |
| **MR** (middle range) | Organic–ion and ion–ion specific interactions (binary cation-anion) | **✓ implemented** |
| **LR** (long range, Debye–Hückel) | Electrolyte screening at finite ionic strength | **✓ implemented** |
| **MR ext.** (triple-ion Qcca, cation-cation Rcc) | Bisulfate / carbonate equilibria specific interactions | not yet |
| **Dissociation equilibria** | Partial dissociation of bisulfate (HSO4⁻ ⇌ H⁺ + SO4²⁻) and carbonate | not yet |

For a non-electrolyte mixture (only neutral organics + water), the MR and LR
contributions are identically zero and **the SR part alone gives the full
AIOMFAC result**. For salt-containing mixtures the full SR + MR + LR is
used. The currently unsupported pieces — triple-ion / cation-cation terms
and the bisulfate equilibrium — only matter for systems containing NH4⁺
together with H⁺ and/or HSO4⁻; everything else is fully supported.

---

## Quick start

### Non-electrolyte mixtures (SR-only)

```python
from aerosysprop.methods.activity_coefficients import aiomfac_sr_from_smiles

# γ_water, γ_ethanol for a 50:50 water-ethanol mixture at 25 °C
gamma = aiomfac_sr_from_smiles(['O', 'CCO'], [0.5, 0.5], 298.15)
print(gamma)            # [1.26196 1.18692]
```

`aiomfac_sr_from_smiles` returns γ on a mole-fraction basis with
pure-component reference state (γ_i → 1 as x_i → 1).

### Electrolyte mixtures (SR + MR + LR)

```python
from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles

# water + ethanol + NaCl from mass fractions
result = aiomfac_from_smiles(
    smiles=['O', 'CCO'],
    salts=['NaCl'],
    mass_fractions=[0.75, 0.20, 0.05],
    T_K=298.15,
)
print(result.gamma_neutral)   # mole-fraction basis, [γ_water, γ_ethanol]
print(result.gamma_cation)    # molality basis, [γ_Na+]
print(result.gamma_anion)     # molality basis, [γ_Cl-]
print(result.gamma_mean)      # mean ionic γ_± per salt, [γ_±(NaCl)]
```

For neutrals, γ is on a mole-fraction basis with pure-component reference;
for ions, γ is on a molality basis with reference state of infinite dilution
in pure water (per the standard atmospheric / electrolyte convention). Salt
names follow the AIOMFAC catalog (`'NaCl'`, `'KCl'`, `'NH4Cl'`, `'Na2SO4'`,
`'NH42SO4'` for ammonium sulfate, `'NaNO3'`, `'NH4NO3'`, `'KNO3'`, `'MgCl2'`,
`'Mg(NO3)2'`, `'MgSO4'`, `'CaSO4'`, `'Ca(NO3)2'`, `'Ca(Cl)2'`, `'K2SO4'`,
`'HCl'`, `'HNO3'`).

### Cached system (many compositions)

For repeated evaluation at many compositions, build the system once and
reuse — the per-system setup (subgroup decomposition, parameter lookup,
pure-component reference, MR matrix assembly) is then cached:

```python
from aerosysprop.methods.activity_coefficients import AIOMFAC
from aerosysprop import groups
from rdkit import Chem

neutrals = [groups.aiomfac(Chem.MolFromSmiles(smi)) for smi in ['O', 'CCO']]
system = AIOMFAC(neutrals, salts=['NaCl'])
for w_NaCl in [0.01, 0.05, 0.10]:
    mass_fracs = [0.90 - w_NaCl, 0.10, w_NaCl]   # water, ethanol, NaCl
    result = system.gamma(mass_fracs, 298.15)
    print(w_NaCl, result.gamma_mean[0])
```

---

## API reference

### `AIOMFAC_SR(components)`

Cached system for SR activity-coefficient evaluation. Build once for a fixed
component list; evaluate at any composition / temperature.

| Argument | Type | Meaning |
|---|---|---|
| `components` | sequence of `{int: number}` | One dict per component: subgroup-ID → count (int or float; SMARTS matching can yield fractional counts for ring-shared atoms) |

Subgroup IDs follow the AIOMFAC numbering convention (1–200 for neutral
organic subgroups). Use `aerosysprop.groups.aiomfac(rdkit_mol)` to produce
these dicts from SMILES.

**Methods:**

| Method | Returns |
|---|---|
| `ln_gamma(mole_fractions, T_K)` | `ndarray` of ln γᵢ in component order |
| `gamma(mole_fractions, T_K)` | `ndarray` of γᵢ |

### `aiomfac_sr_from_smiles(smiles, mole_fractions, T_K)`

Convenience wrapper. Decomposes each SMILES via `groups.aiomfac()`, builds an
`AIOMFAC_SR` system, and returns γᵢ.

| Argument | Type | Meaning |
|---|---|---|
| `smiles` | sequence of `str` | One SMILES per component |
| `mole_fractions` | sequence of `float` | Same length as `smiles` |
| `T_K` | `float` | Temperature in Kelvin |

Returns: `numpy.ndarray` of γᵢ, shape `(len(smiles),)`.

### `AIOMFAC(neutral_components, salts=None)`

Full AIOMFAC system: SR + MR + LR. Handles non-electrolyte and electrolyte
mixtures uniformly; for the non-electrolyte case it reduces to the SR-only
result automatically.

| Argument | Type | Meaning |
|---|---|---|
| `neutral_components` | sequence of `{int: number}` | Subgroup-count dict per neutral component, water-first by convention |
| `salts` | sequence of `str`, optional | Salt names from the AIOMFAC catalog (`'NaCl'`, etc.) |

**Method:**

`.gamma(mass_fractions, T_K)` returns an `AIOMFACResult` named tuple with
fields:

| Field | Type | Meaning |
|---|---|---|
| `gamma_neutral` | `ndarray` | mole-fraction-basis γᵢ for each neutral component |
| `gamma_cation` | `ndarray` | molality-basis γ for each distinct cation, indexed by `cation_ids` |
| `gamma_anion` | `ndarray` | molality-basis γ for each distinct anion, indexed by `anion_ids` |
| `gamma_mean` | `ndarray` | mean ionic γ_± per salt, parallel to `salt_names`. Defined as (γ_c^ν_c × γ_a^ν_a)^(1/(ν_c+ν_a)) |
| `cation_ids`, `anion_ids`, `salt_names` | lists | Parallel lookups for the indices above |

### `aiomfac_from_smiles(smiles, salts, mass_fractions, T_K)`

SMILES-friendly entry point for the full SR + MR + LR calculation.

| Argument | Type | Meaning |
|---|---|---|
| `smiles` | sequence of `str` | SMILES for each neutral component (water first) |
| `salts` | sequence of `str` | Salt names from the AIOMFAC catalog |
| `mass_fractions` | sequence of `float` | All components: [neutral_1, ..., salt_1, ...], sum to 1 |
| `T_K` | `float` | Temperature in Kelvin |

Returns: `AIOMFACResult` (see above).

---

## Theory

The SR part of AIOMFAC is a modified UNIFAC model. For each component *i*:

```
ln γᵢ^SR  =  ln γᵢ^C  +  ln γᵢ^R  −  ln γᵢ^R,ref
```

### Combinatorial (entropy) contribution

We use the Voutsas–Tassios (1997) variant of the UNIFAC combinatorial term
that the Fortran reference uses in `SRgcomb`:

```
ln γᵢ^C  =  ln(rᵢ/r̄)  +  1  −  rᵢ/r̄  −  (z/2) qᵢ [ ln(VQR) + 1 − VQR ]
```

with

```
rᵢ  = Σₖ νₖᵢ Rₖ              # van der Waals volume of component i
qᵢ  = Σₖ νₖᵢ Qₖ              # van der Waals surface area
r̄   = Σⱼ xⱼ rⱼ
q̄   = Σⱼ xⱼ qⱼ
VQR = (rᵢ/r̄) · (q̄/qᵢ)
z   = 10  (coordination number, hardcoded in AIOMFAC)
```

Subgroup Rₖ and Qₖ values come from `data.AIOMFAC_RI[k]` and
`data.AIOMFAC_QI[k]`; the matcher provides νₖᵢ via `groups.aiomfac()`.

### Residual (energetic) contribution

Standard UNIFAC residual form, summed over subgroups present in component *i*:

```
ln γᵢ^R  =  Σₖ νₖᵢ Γₖ
```

with the group-level residual

```
Γₖ  =  Qₖ [ 1  −  ln(Sₖ)  −  Σₘ θₘ Ψₖₘ / Sₘ ]

Sₖ  =  Σₘ θₘ Ψₘₖ
θₖ  =  Qₖ Xₖ / Σⱼ Qⱼ Xⱼ        # subgroup surface-area fraction
Xₖ  =  S1ₖ / Σⱼ S1ⱼ            # subgroup mole fraction in mixture
S1ₖ =  Σᵢ νₖᵢ xᵢ
```

The interaction parameter Ψₘₖ links **main groups** (groups of related
subgroups) rather than subgroups directly:

```
Ψₘₖ  =  exp(−aₘ(m), m(k) / T)
```

where `aₘ(m), m(k)` is the asymmetric main-group interaction parameter from
`data.AIOMFAC_SR_INTERACTIONS` (rebuilt from `Fortran ModSRparam.f90:ARR`).
Self-interactions are zero by convention.

### Reference state

AIOMFAC uses the pure-component reference state for non-electrolytes
(γᵢ → 1 as xᵢ → 1). Achieved by computing `ln Γₖ` again with **xⱼ = δᵢⱼ**
(pure i) and subtracting:

```
ln γᵢ^R,ref  =  Σₖ νₖᵢ Γₖ⁽ⁱ⁾,    Γₖ⁽ⁱ⁾ ≡ Γₖ |_(x = eᵢ)
```

The combinatorial part is zero at the pure-component limit by construction —
no reference subtraction needed there.

### Temperature dependence

The default AIOMFAC parameterisation uses

```
Ψₘₖ(T) = exp(−aₘₖ / T)
```

A 3-parameter form (Ganbavale *et al.* 2015) is also defined in the Fortran
source but is only activated for specific system IDs. The Python
implementation currently uses the 1-parameter form, matching every dataset
released in `AIOMFAC-web` v3.04 and prior — which is the active set in this
repo.

### MR (middle-range) — organic–ion + ion–ion

The MR part captures specific (chemically distinct) interactions between
neutral organic groups and ions, and between cations and anions, as a
function of ionic strength I. From Zuend *et al.* (2008, 2011), implemented
following Fortran `ModMRpart.f90:LR_MR_activity`. All ionic-strength
dependences use the Pitzer-like form `exp(−ω·√I)`.

**Org-ion interactions** (2-parameter):

```
B_kc(I) = b1_kc + b2_kc · exp(−1.2 · √I)        # neutral main group k ↔ cation c
B_ka(I) = b1_ka + b2_ka · exp(−1.2 · √I)        # neutral main group k ↔ anion a
```

ω = 1.2 is fixed for all org-ion pairs. The pre-factor b1 captures the
infinite-dilution part, b2 modulates the short-I rise.

**Ion-ion binary interactions:**

```
B_ca(I)  = b1_ca + b2_ca · exp(−ω_ca · √I)
Cn_ca(I) = b3_ca + c2_ca · exp(−ω2_ca · √I)     # 3-body part
```

ω_ca defaults to 0.8 with per-pair overrides; ω2_ca defaults to 0.6 with
per-pair overrides — all stored in `aiomfac_mr_ion_ion_interactions_c1.matrix`
and `_c3.matrix`. Each (cation, anion) pair has six parameters total.

The MR Gibbs-energy contribution and its derivative w.r.t. ionic strength
are summed via the formulas in Fortran `LR_MR_activity()` (lines 239-336).
Outputs:

- **gnmrln(i)** for each neutral component: built from main-group-level
  contributions weighted by ITAB (subgroup composition).
- **gcmrln(c)** for each cation: org-ion + ion-ion contributions weighted
  by ion charge.
- **gamrln(a)** for each anion: same structure.

### LR (long-range Debye-Hückel)

Pitzer–Debye–Hückel form, mole-fraction-based for water solvent. Following
Fortran `LR_MR_activity` lines 156-195:

```
A_DH = 1.327757e5 · √ρ_w / (ε_w · T)^1.5
B_DH = 6.359696 · √ρ_w / (ε_w · T)^0.5
where ρ_w = 997 kg/m³, ε_w = 78.54 (T = 298.15 K reference)

bb  = 1 + B_DH · √I
TermI = bb − 1/bb − 2 log(bb)

gnlrln_i = M_i · (2 A_DH / B_DH³) · TermI       # neutrals scaled by molar mass
gclrln_c = −Z_c² · A_DH · √I / bb                # cations
galrln_a = −Z_a² · A_DH · √I / bb                # anions
```

### Reference-state conversion for ions

The LR and MR contributions for ions come out on a mole-fraction basis but
ion activity coefficients in the electrolyte literature are conventionally
reported on a molality basis with reference state of infinite dilution in
water. The conversion factor (Fortran line 348) is subtracted from the
ln γ of each ion:

```
T_molal = log(M_water / M̄_solvent + M_water · Σ molalities_ions)
ln γ_c,molality = ln γ_c,SR + ln γ_c,MR + ln γ_c,LR − T_molal
ln γ_a,molality = ln γ_a,SR + ln γ_a,MR + ln γ_a,LR − T_molal
```

where M̄_solvent is the mean molar mass of the salt-free solvent mixture.

### Mean ionic activity coefficient

For a salt with stoichiometry ν_c cations + ν_a anions (e.g., NaCl has
ν_c = ν_a = 1, Na₂SO₄ has ν_c = 2 ν_a = 1):

```
γ_± = (γ_c^ν_c · γ_a^ν_a)^(1 / (ν_c + ν_a))
```

Returned in `AIOMFACResult.gamma_mean`.

---

## Compound coverage

### Non-electrolyte (Phase 1 SR)

The SR regression suite (`aerosysprop/tests/test_activity_coefficients.py`)
validates the Python output against Fortran AIOMFAC for **1520 fixtures**,
spanning **30 distinct main groups**:

| Fixture group | Count | Notes |
|---|---|---|
| Curated water-binaries | 79 | Diverse hand-picked compounds across functional-group classes |
| MCM331 unique-signature water-binaries | 1432 | Every closed-shell species in `tests/SMILES/MCM331.xml` whose subgroup decomposition is unique and Fortran-compatible |
| Multi-component mixes | 9 | 3-to-4-component mixtures from different functional-group classes |

**Functional-group classes exercised** (curated set):
alcohols (primary / secondary / tertiary), diols, polyols, mono- and
dicarboxylic acids, hydroxy- and keto-acids, esters, ethers, PEG, aliphatic
and aromatic aldehydes, ketones, simple aromatics, sugars, alkenes, amines,
amides, nitriles, halides, hydrocarbons, plus 12 α-pinene-oxidation
products. Cross-class multi-component mixes include water + ethanol +
glycerol, water + acid + ester, dicarboxylic-acid + alcohol, and a
water + pinic-acid + pinonic-acid + levoglucosan SOA proxy.

**MCM coverage:** of 5809 species in MCM331.xml, 2193 are radicals (open
shell, AIOMFAC not applicable), 1641 hit the aerosysprop 281+ extension
subgroups (no Fortran counterpart — see [Known limitations](#known-limitations)),
and 1968 are Fortran-compatible. After de-duplicating by subgroup signature,
1432 unique decompositions are tested.

**Tolerance:** every fixture matches to **5 × 10⁻⁵ relative error** (the
precision of the Fortran printed output, which is 5–6 significant figures).

### Electrolyte (Phase 2 SR + MR + LR)

The electrolyte regression suite
(`aerosysprop/tests/test_activity_coefficients_electrolytes.py`) validates
the full SR + MR + LR pipeline against Fortran AIOMFAC for **43 fixtures**
spanning every common atmospheric salt, all combinations of multiple
cations and anions, multi-organic + multi-salt mixtures, divalent
cations, and edge concentrations. Breakdown:

| Category | # fixtures | What it exercises |
|---|---|---|
| Water + single salt (alkali / ammonium chlorides, nitrates, sulfates) | 9 | Individual cation/anion parameter sets at multiple concentrations |
| Water + divalent cation single salt (MgCl₂, CaCl₂, MgSO₄) | 3 | Z² scaling in LR and MR for doubly-charged cations |
| Water + acid (HCl, HNO₃) | 2 | H⁺ as a cation |
| Water + 2 salts | 6 | Cross-coupling: same cation with two different anions, two cations with same anion, full 2×2 ion-ion pair activation |
| Water + 3-4 salts | 3 | All-NH₄⁺ with three anions, mixed alkali / ammonium / nitrate / sulfate / chloride |
| Water + mixed divalent + monovalent | 2 | Mg²⁺ ↔ Ca²⁺ ↔ Cl⁻; NH4⁺ + Mg²⁺ + SO4²⁻ |
| Water + single organic + single salt | 6 | Org-ion B_kc and B_ka coupling for various organic functional groups |
| Water + multiple organics + single salt | 4 | Multi-organic mixing with shared ion environment (incl. a 4-organic case) |
| Water + multiple organics + multiple salts | 5 | True atmospheric SOA + inorganic mixtures (urban, marine, biomass-burning, α-pinene-SOA proxies) |
| Edge concentrations | 3 | Very dilute (10⁻⁴), high concentration up to 40%, near saturation |
| High-organic-count stress tests | 6 | 5, 10, and 20 organics with one or two salts — realistic atmospheric organic-rich mixtures with all functional groups co-present in a single phase |

Notable atmospheric aerosol proxies validated:

- **urban_aerosol_proxy**: water + ethanol + acetic-acid + (NH₄)₂SO₄ + NH₄NO₃
- **marine_SOA_seasalt_sulfate**: water + pinonic-acid + glucose + (NH₄)₂SO₄ + NaCl
- **biomass_burning_proxy**: water + levoglucosan + acetic-acid + NaCl + (NH₄)₂SO₄
- **alpha_pinene_SOA_mixed_inorg**: water + pinic-acid + pinonic-acid + (NH₄)₂SO₄ + NH₄NO₃
- **20_organic_atm_SOA_sulfate_nitrate**: water + 20 organics (terpene SOA + dicarboxylic acids + sugars + alcohols + aromatics) + (NH₄)₂SO₄ + NH₄NO₃ — the most demanding single mixture in the suite, exercising 9 main groups simultaneously with multiple cations and anions

Every composition point matches Fortran to **5 × 10⁻⁵ relative error** for
every species (neutrals, cations, anions, and mean ionic γ_±).

Total Phase 1 + Phase 2 pytest count: **3152 cases** (1524 SR fixtures × 2
input variants + 2 sanity + 50 electrolyte fixtures × 2 input variants +
2 sanity), running in ~10 s.

**Ion ordering convention:** The Python API returns `gamma_cation` and
`gamma_anion` in *first-encountered order* across the user-supplied salt
list — matching Fortran AIOMFAC's output convention. If you pass salts as
`['NH42SO4', 'NH4NO3']`, the anions come back in order `[SO4²⁻, NO3⁻]`.
The `cation_ids` and `anion_ids` fields of `AIOMFACResult` give the
explicit ID lookup table.

---

## Validation methodology

Reference outputs are generated by running the Fortran AIOMFAC executable on
inputs containing the same SMILES-derived subgroup decompositions the Python
code uses. Both sides therefore see identical inputs; only the activity-
coefficient maths is independent. The two test variants:

- **`test_sr_matches_fortran`** — feeds subgroup dicts pulled from the JSON
  fixture directly into `AIOMFAC_SR`. Isolates the SR maths.
- **`test_sr_smiles_pipeline_matches_fortran`** — feeds the SMILES from the
  JSON fixture through `aiomfac_sr_from_smiles()`. Catches drift in
  `groups.aiomfac()` or in any pre-processing.

### Regenerating fixtures

The full pipeline (Python → Fortran → JSON) is in
`aerosysprop/tests/generate_aiomfac_sr_fixtures.py`. To refresh after
parameter changes or matcher updates:

```bash
# macOS only: gfortran (Homebrew gcc) needs the Xcode SDK to find libSystem
export SDKROOT=$(xcrun --show-sdk-path)

# Build the Fortran reference if not already built
# (see https://github.com/andizuend/AIOMFAC)

python aerosysprop/tests/generate_aiomfac_sr_fixtures.py
```

Build expects the Fortran executable at
`/Users/davidtopping/Documents/GitHub/AIOMFAC/FortranCode/AIOMFAC-web.out`
— edit `FORTRAN_EXE` in the generator if yours is elsewhere.

Skipped compounds (decomposition incompatible with Fortran) are logged to
`aerosysprop/tests/aiomfac_sr_skipped.txt` with reasons.

The generator takes ~60 s for the full 1520-fixture build. For fast
iteration on the curated set only, set `INCLUDE_MCM = False` at the top of
the script.

### Running the tests

```bash
pytest aerosysprop/tests/test_activity_coefficients.py
```

Expected: **3042 passed in ~10 s** (1520 fixtures × 2 test variants + 2
sanity tests).

---

## Data files

All AIOMFAC parameter data lives in `aerosysprop/data/`:

| File | Contents | Format |
|---|---|---|
| `aiomfac_main.data` | Per-subgroup main group, molar mass, Rₖ, Qₖ | `subgroup_id  main_group  mass  R  Q  # name` |
| `aiomfac_sr_interactions.matrix` | Main-group interaction parameters aₘₖ (Ψ = exp(−aₘₖ/T)) | 76×76 matrix, row & col labels are main-group IDs; asymmetric |
| `aiomfac_unifac.smarts` | SMARTS patterns for organic subgroups | `subgroup_id  SMARTS` |
| `aiomfac_hctail.smarts` | SMARTS for alcohol hydrocarbon-tail logic | as above |
| `aiomfac_ions.smarts` | SMARTS for ion species (Phase 2) | as above |
| `aiomfac_ions.data` | Per-ion charges (Phase 2) | `ion_id  charge` |
| `aiomfac_salts.data` | Salt → cation/anion decompositions (Phase 2) | tabular |
| `aiomfac_mr_*.matrix` | Middle-range interaction parameters (Phase 2) | as above |

All these are loaded at package import time by `aerosysprop/data/__init__.py`
and exposed as module-level `AIOMFAC_*` dicts and nested dicts.

---

## Diagnostics

Several failure modes in AeroSysProp are silent — the matcher may drop
functional groups it has no SMARTS for, the AIOMFAC extension subgroups
281-285 silently disappear at the Fortran reference level, an organic
main group may lack a parameter against a salt's ion, mass fractions may
not sum to 1, etc. The :mod:`aerosysprop.diagnostics` module flags every
known issue class.

### Surface via warnings (default) or raise

The AIOMFAC entry points accept a ``strict=`` argument
(default ``'warn'``):

```python
from aerosysprop.methods.activity_coefficients import aiomfac_from_smiles

# default: emits warnings.warn for each detected issue, then proceeds
result = aiomfac_from_smiles(['O', 'CC(=O)C'], ['NaCl'], [0.85, 0.10, 0.05], 298.15)
# UserWarning: [WARNING] EXTENSION_SUBGROUP 'CC(=O)C' decomposes to extension subgroup(s) [281] ...
# UserWarning: [ERROR  ] SUBGROUP_OUT_OF_RANGE 'CC(=O)C' produces subgroup IDs outside Fortran ...

# raise: ValueError if anything is flagged
aiomfac_from_smiles(..., strict='raise')   # raises ValueError

# silent: original behaviour, useful inside hot loops where you've
# pre-validated inputs
aiomfac_from_smiles(..., strict='silent')
```

### Pre-flight check without computing

To screen inputs for ANY aerosysprop method without actually running the
calculation, use ``diagnose()``:

```python
from aerosysprop.diagnostics import diagnose, diagnose_aiomfac
from aerosysprop.methods import boiling_points, vapour_pressures

# Will the Nannoolal boiling-point method see the whole molecule?
report = diagnose(boiling_points.nannoolal, 'CS(=O)C')   # DMSO
report.is_ok                  # False
for issue in report:
    print(issue.severity, issue.code, issue.message)
# warning UNMATCHED_ATOMS  2 heavy atom(s) in 'CS(=O)C' are not matched
#   by any SMARTS pattern; the method will silently behave as if these
#   atoms (and their attached functional groups) were absent.

# Full AIOMFAC pre-flight in one call:
report = diagnose_aiomfac(
    smiles=['O', 'Cc1ccccc1'],          # water + toluene
    salts=['NH42SO4'],
    mass_fractions=[0.85, 0.10, 0.05],
    T_K=298.15,
)
# Flags: main_group 4 has no MR parameter against anion 261 (SO4--)
```

### Issue codes

| Severity | Code | Detected when... |
|---|---|---|
| error | `INVALID_SMILES` | RDKit can't parse the SMILES |
| error | `EMPTY_DECOMPOSITION` | matcher returned no subgroups |
| warning | `UNMATCHED_ATOMS` | one or more heavy atoms have no SMARTS match (silent truncation) |
| warning | `EXTENSION_SUBGROUP` | subgroup 281–285 (aerosysprop's 2012 extension; Fortran drops these) |
| error | `SUBGROUP_OUT_OF_RANGE` | subgroup ID outside 1..265 |
| warning | `FRACTIONAL_COUNT` | non-integer subgroup count (ring-shared atoms) |
| error | `NEGATIVE_COUNT` | negative count (matcher bug) |
| error | `MR_PARAM_MISSING` | required (main group, ion) MR parameter has Fortran's `-7.7777e4` "not estimated" sentinel |
| warning | `COMPOSITION_SUM` | mass / mole fractions don't sum to ≈1 |
| error | `COMPOSITION_VALUE` | values outside [0, 1] or wrong vector length |
| warning | `NO_WATER` | electrolyte mixture without water as a neutral |
| warning | `HIGH_CONCENTRATION` | a salt is above 30% mass fraction (parameter set may not have been fitted here) |
| warning | `BISULFATE_NEEDED` | H₂SO₄ in mix (partial dissociation not modelled) |
| warning | `QCCA_NEEDED` | NH₄⁺ + H⁺ both present (triple-ion terms not modelled) |
| error | `UNKNOWN_SALT` | salt name not in the AIOMFAC catalog |

### SMARTS hit matrix for inspection

For bulk inspection of how the matchers see a list of SMILES — for
example, sanity-checking matcher decompositions across a batch, or
fingerprint-style comparison — use the companion
:mod:`aerosysprop.fingerprints` module:

```python
from aerosysprop.fingerprints import smarts_hits, available_methods

available_methods()
# ['aiomfac', 'aiomfac_hctail', 'aiomfac_ions', 'simpol', 'nannoolal_primary',
#  'nannoolal_secondary', 'evaporation', 'evaporation2', 'stein_and_brown',
#  'girolami', 'le_bas', 'schroeder', 'composition', 'moller', 'fp4', 'maccs']

df = smarts_hits(
    smiles=['O', 'CCO', 'OCC(O)CO', 'Cc1ccccc1'],
    method='aiomfac',
    drop_zero_columns=True,
)
#             1   2   3   9   11  14  16
# O           0   0   0   0   0   0   1
# CCO         1   1   0   0   0   1   0
# OCC(O)CO    0   2   1   0   0   3   0
# Cc1ccccc1   1   0   0   5   1   0   0
```

For side-by-side comparison across multiple SMARTS sets, use
``smarts_hits_all(smiles, methods=[...])`` which returns a single
DataFrame with a ``MultiIndex(method, id)`` on the columns.

### Diagnostics on the result object

When you call `aiomfac_from_smiles` (full SR+MR+LR) with the default
`strict='warn'` or `strict='raise'`, the same diagnostic report that
drives the warning/raise behaviour is also attached to the result as
`result.diagnostics`. This lets you introspect every issue after the
fact without re-running the pre-flight check:

```python
result = aiomfac_from_smiles(
    smiles=['O', 'CC(=O)C'], salts=['NaCl'],
    mass_fractions=[0.85, 0.10, 0.05], T_K=298.15,
)
result.gamma_neutral       # numeric output as before
result.diagnostics         # DiagnosticReport object (or None if strict='silent')
result.diagnostics.errors  # list of Issue objects with severity='error'
for issue in result.diagnostics:
    print(issue.code, issue.message)
```

The SR entry point follows the same pattern but via an opt-in
`return_diagnostics=True` argument (its default return type stays a
bare numpy array for backward compatibility):

```python
gamma, report = aiomfac_sr_from_smiles(
    ['O', 'CCO'], [0.5, 0.5], 298.15,
    return_diagnostics=True,
)
report.is_ok
```

### Performance in hot loops

AIOMFAC is built to be called from numerical simulation frameworks
(partitioning solvers, Monte Carlo, droplet growth models, equilibrium
iterators) that call the method **tens of thousands to millions of
times** with a fixed SMILES list but varying composition. Two layers
of caching make the SMILES entry points as fast as the manually cached
system for that pattern:

1. **System cache** — `aiomfac_from_smiles` and `aiomfac_sr_from_smiles`
   memoise the constructed `AIOMFAC` / `AIOMFAC_SR` system by
   `(smiles tuple, salts tuple)`. After the first call, subsequent
   calls with the same SMILES list skip SMILES parsing, AIOMFAC
   subgroup decomposition, and MR parameter table assembly. 128-entry
   LRU.

2. **Diagnostic cache** — the per-SMILES portion of the pre-flight
   check (matcher + atom-coverage SMARTS) is also memoised by SMILES
   so `strict='warn'` adds essentially no overhead on a warm cache.
   4096-entry LRU.

Concrete numbers (M1, single-thread, 10-organic + 1-salt mixture):

| Call style | Time per call | Throughput |
|---|---|---|
| ``aiomfac_from_smiles(..., strict='silent')`` | 0.092 ms | ~11k calls/s |
| ``aiomfac_from_smiles(..., strict='warn')`` (default) | 0.102 ms | ~10k calls/s |
| ``AIOMFAC(...).gamma(...)`` (manual cached system, reference) | 0.090 ms | ~11k calls/s |

So the default behaviour is **within 10% of the theoretical floor** for
this access pattern. The first call pays a one-time ~3 ms cost (system
construction); every call after that is free.

For 1 million calls in a numerical simulation:

- Pre-system-cache implementation: ~50 minutes (3 ms × 1M)
- Current implementation, default settings: ~100 seconds (0.10 ms × 1M)

#### Cache controls

```python
from aerosysprop.methods.activity_coefficients import (
    system_cache_info,    # ((aiomfac_info, aiomfac_sr_info))
    clear_system_cache,   # drop all cached AIOMFAC / AIOMFAC_SR systems
)
from aerosysprop.diagnostics import (
    diagnostic_cache_info,   # cache stats for per-SMILES diagnostics
    clear_diagnostic_cache,  # drop all diagnostic-analysis entries
)

ai, sr = system_cache_info()
print(ai)   # CacheInfo(hits=4099, misses=1, maxsize=128, currsize=1)
```

#### Workloads with unbounded distinct SMILES

The 128-entry LRU is sized for typical simulation patterns (a handful
of distinct mixtures). For screening / fingerprinting workloads with
millions of *distinct* SMILES lists, the cache itself grows to its
LRU bound and old entries are evicted automatically — no leak, but
the cache effectively reverts to "one-shot" behaviour and each call
pays full construction cost. Either call ``clear_system_cache()``
periodically to keep memory tight, set ``strict='silent'`` to skip
the diagnostic memoisation entirely, or — if you have a specific
fast set of compounds you want to keep hot — construct the
``AIOMFAC`` system once explicitly and reuse it directly:

```python
from aerosysprop.methods.activity_coefficients import AIOMFAC
from aerosysprop import groups
from rdkit import Chem

neutrals = [groups.aiomfac(Chem.MolFromSmiles(s)) for s in smiles]
system   = AIOMFAC(neutrals, salts=['NaCl'])
for x in composition_grid:
    result = system.gamma(x, 298.15)
```

This is also the recommended pattern when integrating AIOMFAC into a
larger framework where you want explicit control over which mixture
state is held in memory.

### Inspecting reports

`DiagnosticReport` is iterable and exposes structured access:

```python
report = diagnose_aiomfac([...], [...], [...], 298.15)
report.is_ok          # True iff zero issues
report.has_errors
report.has_warnings
report.errors         # list[Issue]
report.warnings       # list[Issue]
report.summary()      # multi-line human-readable string
for issue in report:
    issue.code        # str — e.g., 'EXTENSION_SUBGROUP'
    issue.severity    # 'error' | 'warning' | 'info'
    issue.message     # str
    issue.details     # dict — structured context (e.g., {'subgroups': [281]})
```

## Known limitations

### Compounds the matcher can't decompose cleanly

`groups.aiomfac()` silently drops some functional groups it lacks SMARTS for:

| Compound | What's dropped | Workaround |
|---|---|---|
| Nitromethane, nitrobenzene | -NO₂ | None — uses CH₃ / aromatic only |
| Methanethiol, DMSO | -SH, S=O | None |
| Methyl nitrate, formyl peroxynitrate | -ONO₂, peroxynitrate | None |
| Formamide | Everything (returns `{}`) | None |
| Many MCM hydroperoxides | -OOH | Filtered by 281+ extension check; ~1641 MCM species affected |

For the *dropped* compounds, Python and Fortran *agree* on the truncated
representation (both see the same dict), so the regression suite still
passes, but the resulting γ is for the partial molecule. Use case dictates
whether this is acceptable.

### `281`–`285` extension subgroups

`aerosysprop/data/__init__.py` adds "Mark Barley 2012" subgroups 281–285
(renamed from 191/231/261/531/561) to handle MCM compounds, but these don't
exist in Fortran AIOMFAC (its `topsubno = 265`). For affected compounds
(acetone, methyl hydroperoxide, peracetic acid, etc.) the matcher emits
subgroup IDs the Fortran reference silently drops. The fixture generator
filters these out before running Fortran.

### Glycerol triol-parser bug (FIXED 2026-06-02)

`groups.aiomfac()`'s triol branch (alcohols_alone == 3) used to compute a
broken `c_atoms_between_oh` for compactly-packed triols, producing a
negative count for subgroup 142 (CH₂[alc]). Glycerol was the canonical
example. The formula has been replaced with `max(distances) − (alcohols_alone − 2)`,
which gives 0 for glycerol, 1 for 1,2,4-butanetriol, 2 for 1,3,5-pentanetriol.
Glycerol now decomposes correctly to `{150: 2, 151: 1, 153: 3}` and is
exercised in the regression suite (binary water + glycerol, plus
water + glycerol + NaCl in Phase 2). Three additional MCM compounds that
previously triggered the same code path are now Fortran-compatible.

### Single-temperature validation

All fixtures are at T = 298.15 K. The temperature law
Ψ(T) = exp(−a/T) is well-tested by physics, but multi-temperature
regression cases would extend coverage. Adding them is a small generator
edit (see the `BINARY_X_SWEEP` constant for the composition sweep — an
analogous `TEMPERATURE_SWEEP` would extend each fixture to multiple Ts).

### Unparameterised organic-ion pairs

Not every (organic main group, ion) pair has fitted parameters in the
AIOMFAC parameter set. The Fortran source flags those cells with the
sentinel value `-7.7777e4`. When you build an `AIOMFAC` system whose
organics contribute a main group lacking a parameter for one of the ions
in the mixture, Python now raises `ValueError` with a clear diagnostic —
matching Fortran's error code 1/2 behaviour. The most common case is
**main group 4 (ACCHn — aromatic-attached alkyls)** paired with SO4²⁻ or
NO3⁻: toluene, xylene, m-cresol, and similar compounds with `c-CH3` /
`c-CH2-` substituents will fail when combined with sulfate or nitrate
salts. Substitute another aromatic (vanillin, phenol, catechol,
benzaldehyde, benzoic acid — all using main groups 3/8/9/13/20) if you
need that combination.

### Bisulfate equilibrium (H2SO4 / HSO4⁻)

Sulfuric acid partially dissociates in solution:
`H2SO4 → H⁺ + HSO4⁻ ⇌ 2 H⁺ + SO4²⁻`. The Fortran code solves this
equilibrium iteratively (`SubModDissociationEquil.f90`), producing
*both* HSO4⁻ and SO4²⁻ ion populations at every composition. The current
Python port treats H2SO4 as fully dissociated into 2 H⁺ + SO4²⁻ and so
will disagree with Fortran for any system involving H2SO4 (or generally
any case with H⁺ + HSO4⁻). Same caveat applies to the carbonate
equilibrium (CO3²⁻ ⇌ HCO3⁻ ⇌ H2CO3). Mark as Phase 3.

### MR triple-ion (Qcca) and cation-cation (Rcc) terms

The Fortran MR also includes Qcca (cation-cation-anion 3-body) and Rcc
(cation-cation binary) interaction terms. These are only activated for
mixtures containing both NH4⁺ and H⁺ (plus optionally HSO4⁻); they're
otherwise zero. The Python port skips them — fine for any mixture without
this specific cation combination. Combined with the bisulfate-equilibrium
limitation, this means **systems with NH4⁺/H⁺/HSO4⁻ simultaneously
present are out of scope** for the current Python AIOMFAC.

### Ions in the MR catalog

The Python port handles every cation/anion the Fortran MRinteractcoeff
subroutine assigns parameters for: cations Li⁺ (201), Na⁺ (202), K⁺ (203),
NH4⁺ (204), H⁺ (205), Ca²⁺ (221), Mg²⁺ (223); and anions Cl⁻ (242),
Br⁻ (243), I⁻ (244), NO3⁻ (245), IO3⁻ (246), OH⁻ (247), HSO4⁻ (248),
ClO4⁻ (249), HCO3⁻ (250), Hmalo⁻ (251), Hglut⁻ (252), Hsucc⁻ (253),
MeOS⁻ (254), EtOS⁻ (255), IsopOS⁻ (256), SO4²⁻ (261), CO3²⁻ (262),
malo²⁻ (263), glut²⁻ (264), succ²⁻ (265). Ions outside this set are
unparameterised in Fortran too (no Phase 2 capability is being lost).

---

## Data-fix log

When the Python port was first validated against Fortran on a curated set,
several latent issues in the bundled aerosysprop data were uncovered. All
have been fixed at source; the regression suite locks them in.

| Date | File / function | Issue | Fix |
|---|---|---|---|
| 2026-05-22 | `groups.py:1558` | `aiomfac()` returned subgroup 14 (old OH, main group 5) for alcohols; Fortran rejects with error 13 when paired with main group 68 | Append a 14→153 conversion at the end of the function when any of subgroups 149–152 are present |
| 2026-05-22 | `aiomfac_main.data` | Subgroup 151 main group: 68 → **52** (Fortran NKTAB has separate `CH[OH]` for sec/tert alcohols) | Patched lines |
| 2026-05-22 | `aiomfac_main.data` | Subgroup 152 main group: 68 → **52** | Patched lines |
| 2026-05-22 | `aiomfac_main.data` | Subgroup 137 main group: 65 → **20** | Patched line |
| 2026-05-22 | `aiomfac_sr_interactions.matrix` | Rows 52, 53, 54, 56, 58, 60, 61, 63, 64, 65, 66, 67, 68, 69 were near-duplicates of standard rows (1/5/9/13/20). Inherited "pre-fit" placeholder values from an older AIOMFAC version, never refreshed when the Fortran code refitted those parameters | Regenerated the whole 76×76 matrix from `AIOMFAC/FortranCode/ModSRparam.f90:set_SR_group_interaction_param` |
| 2026-05-24 | `aiomfac_main.data` | Subgroup 18 (CH₃CO) Q: 1.448 → **1.488** (one-digit typo affected all ketones — pyruvic, pinonic, MEK) | Patched line |
| 2026-05-24 | `aiomfac_main.data` | Subgroup 154 (CH₂OCH₂ PEG) R: 1.5927 → **1.38143**, Q: 1.32 → **3.0** (Fortran applies a PEG-specific override in `SRsystm`) | Patched line |
| 2026-05-22 | `aerosysprop/data/__init__.py:91` | `_read_matrix` exception handler used `e.args += 'string'` (str + tuple TypeError), masking the real parsing error | Changed to `e.args += ('string',)` |
| 2026-06-02 | `aiomfac_mr_*.matrix` (org-ion b1, b2 + ion-ion b1, b2, b3, c2, c3) | All MR matrices were stale (Na⁺/Cl⁻ b1 was 0.034 in Python vs 0.054 in Fortran, etc.). | Regenerated all 6 MR matrices from `AIOMFAC/FortranCode/ModMRpart.f90:MRdata + MRinteractcoeff` (76 main groups × 27 ions for org-ion; 7 cations × 20 anions for ion-ion). Also added the previously-missing `aiomfac_mr_ion_ion_interactions_c1.matrix` for per-pair `omegaTAB` overrides. |
| 2026-06-02 | `aiomfac_salts.data` | NH4Cl mass was 62.491 g/mol; correct is 53.491 (NH4⁺ 18.04 + Cl⁻ 35.45). Old typo. | Patched. |
| 2026-06-02 | `groups.py:1487-1530` (triol branch) | The `c_atoms_between_oh` formula returned a negative value for compactly-packed triols like glycerol (OHs at C1, C2, C3 — SMARTS 100 matches multiple times, fall-through gate fails). | Replaced with the chemically correct `max(distances) − (alcohols_alone − 2)`. Glycerol now decomposes to `{150:2, 151:1, 153:3}`; 3 additional MCM triols become Fortran-compatible. |
| 2026-06-02 | `aiomfac_mr_org_ion_interactions_b1.matrix`, `b2.matrix` + `AIOMFAC` class | Extractor was silently converting Fortran's "-7.7777e4" not-estimated sentinels to 0, so Python computed wrong-but-finite gammas instead of erroring when mixtures used unparameterised organic-ion pairs. | Sentinel now preserved in the data files; `AIOMFAC.__init__` raises `ValueError` with a clear diagnostic if any required organic-ion slot is unparameterised — matching Fortran's error code 1/2 behaviour. |

Subgroups **244, 246, 247** also have R/Q mismatches vs Fortran but they are
ions whose R/Q values are not used by the MR/LR equations (MR uses per-ion
charges and interaction parameters, not R/Q). No Phase 2 impact.

---

## Citations

Please cite the original AIOMFAC papers when using this implementation:

- Zuend, A., Marcolli, C., Luo, B. P., and Peter, T. (2008). A thermodynamic
  model of mixed organic–inorganic aerosols to predict activity coefficients.
  *Atmospheric Chemistry and Physics*, **8**, 4559–4593.
  [doi:10.5194/acp-8-4559-2008](https://doi.org/10.5194/acp-8-4559-2008).

- Zuend, A., Marcolli, C., Booth, A. M., Lienhard, D. M., Soonsin, V.,
  Krieger, U. K., Topping, D. O., McFiggans, G., Peter, T., and Seinfeld,
  J. H. (2011). New and extended parameterization of the thermodynamic model
  AIOMFAC: calculation of activity coefficients for organic–inorganic
  mixtures containing carboxyl, hydroxyl, carbonyl, ether, ester, alkenyl,
  alkyl, and aromatic functional groups.
  *Atmospheric Chemistry and Physics*, **11**, 9155–9206.
  [doi:10.5194/acp-11-9155-2011](https://doi.org/10.5194/acp-11-9155-2011).

For the UNIFAC functional-group methodology this builds on:

- Fredenslund, A., Jones, R. L., and Prausnitz, J. M. (1975).
  Group-contribution estimation of activity coefficients in nonideal liquid
  mixtures. *AIChE Journal*, **21**, 1086–1099.

- Voutsas, E. C., and Tassios, D. P. (1997). Analysis of the UNIFAC-type
  group-contribution models at the highly dilute region.
  *Industrial & Engineering Chemistry Research*, **36**, 4965–4972.

For atmospheric-functional-group extensions (peroxy/organonitrate, not yet
matched here):

- Compernolle, S., Ceulemans, K., and Müller, J.-F. (2009). Influence of
  non-ideality on condensation to aerosol.
  *Atmospheric Chemistry and Physics*, **9**, 1325–1337.

---

## Roadmap

| Phase | Scope | Status |
|---|---|---|
| **1** | SR (modified UNIFAC) for neutral organic + water mixtures | **complete** |
| **2** | MR (organic–ion, ion–ion binary) + LR (Debye–Hückel) for electrolyte mixtures | **complete** |
| **3a** | Bisulfate / carbonate dissociation equilibria (H2SO4 ⇌ H⁺ + HSO4⁻ ⇌ 2 H⁺ + SO4²⁻); MR triple-ion Qcca and cation-cation Rcc terms (for NH4⁺/H⁺/HSO4⁻ systems) | not started |
| **3b** | Extend `groups.aiomfac()` to emit AIOMFAC subgroups 155–172 (Compernolle 2009 peroxy / organonitrate / PAN) so the ~1641 MCM species currently filtered can be modelled properly | not started |
| **4** | Temperature-sweep validation (multi-T fixtures) | not started |
