# AeroSysProp: guide for AI agents

AeroSysProp estimates properties of organic molecules relevant to atmospheric aerosol from their structure (SMILES), using published group-contribution methods, and computes AIOMFAC activity coefficients for mixtures of water, organics and inorganic salts. Every number it returns is an estimate. Treat it as a model output with known limits, not as a measured value.

## Choosing a tool

| Task | Tool |
|---|---|
| Find method IDs, units and dependencies | `list_methods` |
| Check a SMILES or see which methods cover a molecule | `diagnose_smiles` |
| Boiling point, vapour pressure, critical properties, liquid density | `predict_properties` |
| Non-ideal mixing of water, organics and salts | `activity_coefficients` (with `list_salts`) |
| Molecules held in a file | `list_data_files`, then `load_smiles_file` |
| Results saved earlier from the AeroSysProp app | `read_session_database` |
| Independent QSAR cross-check (if the OPERA container runs) | `opera_predict` |

## Units

Boiling and critical temperatures are in K, critical pressure in bar, liquid density in g/cm3. Vapour pressure is **log10 of the saturation vapour pressure in atmospheres** at the requested temperature, so -6 means 1e-6 atm. To convert to Pa, compute 10^value × 101325. Temperatures are always kelvin. OPERA uses different units (°C, log10 mmHg); do not mix them with AeroSysProp values without converting.

## Reading diagnostics

Results carry diagnostics per method. The important codes:

* `UNMATCHED_ATOMS`: the method's functional groups do not cover some heavy atoms. The method behaves as though those atoms (and the groups attached to them) were absent. Treat the value as unreliable and say so.
* `INVALID_SMILES`: RDKit could not parse the input. Check the SMILES rather than guessing.
* AIOMFAC codes such as `MR_PARAM_MISSING`, `EXTENSION_SUBGROUP`, `HIGH_CONCENTRATION` and `BISULFATE_NEEDED` are described in `list_methods` under the mixture method, and in the resource `aerosysprop://docs/aiomfac`.

`not_computed` lists outputs a method could not produce for that molecule (the value is `null`).

## Good practice

* **Compare methods.** Vapour pressure estimates for multifunctional compounds can differ by orders of magnitude between methods. When asked for "the" vapour pressure, run several (`methods: ["vapour_pressure"]`) and report the range and which methods agree, rather than silently choosing one.
* **Boiling point source matters.** Nannoolal, Myrdal & Yalkowsky vapour pressures, Nannoolal critical properties and the Schroeder, Le Bas and Tyn & Calus densities all take a boiling point as input. Pass several `bp_sources` to see the sensitivity.
* **EVAPORATION and SIMPOL** do not need a boiling point and were developed with atmospheric oxygenated organics in mind, which makes them useful comparisons for aerosol-relevant compounds.
* **Batch your calls.** One call with 50 SMILES is far better than 50 calls.
* **Mixtures.** Water is `"O"`. With salts, give `mass_fractions` for the neutral components first, then the salts, in the same order as the `smiles` and `salts` lists. Neutral-component activity coefficients are on a mole-fraction basis; ion coefficients are on a molality basis.

## Provenance

Every result includes a `provenance` block. Each method has a label such as `vp_simpol@1 (3f9a1c2b7d40)`: the method ID, its parameter set version and a fingerprint of the parameter files. The method's card, summarised under `card` in `list_methods`, gives its publication status. When reporting values in anything a person may cite, include the provenance label, and say when a method is `published-with-modifications` or `in-house` rather than presenting it as the literature method. A `PARAMETERS_CHANGED` warning means the installed parameter files differ from the approved version.

## What this server will not do

It cannot write files, save results or change a session database, and it can only read files inside the directory the user configured. If a user wants results saved, return them and let the user save them through the AeroSysProp app or their own code.
