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
| More than about 50 molecules, or a whole file | `start_batch_run`, then `batch_run_status` and `read_batch_results` |
| Earlier batch runs | `list_batch_runs` |

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

Some methods have more than one parameter set: the same equations with different group definitions. `list_methods` lists them under `parameter_sets`. EVAPORATION has two; set 2 is the default. Request a set with `"<id>@<n>"` (for example `"vp_evaporation@1"`, or both sets to compare them). Output columns always carry the set, as in `vp_evaporation_log10atm__set2`, and provenance labels read `vp_evaporation@2 (...)`. Always say which set a value came from. The old ID `vp_evaporation2` still works but returns a `DEPRECATED_METHOD_ID` warning.

## Batch runs

`predict_properties` returns every value, which is only practical for a few tens of molecules. For larger sets use `start_batch_run` (up to 100,000 molecules; the full MCM of about 5,800 species takes around 40 s on 4 cores). It writes a new run folder containing `results.csv`, `diagnostics.csv`, `run.json` (the full provenance manifest, with a fingerprint of each output file) and `ro-crate-metadata.json` (the same provenance in the RO-Crate community standard). Then:

1. Poll `batch_run_status` until `status` is `complete`. The summary gives a `coverage` statement, `caveats`, `reporting_requirements`, the range of every column for **clean** molecules (closed-shell, no method flagged) and for **all** molecules, and `method_spread`: how far methods disagree, computed on clean molecules.
2. Use `read_batch_results` for detail: `view: "flagged"` (optionally one `method`), `view: "largest_spread"` with a `group` exactly as named in `method_spread`, or `names` / `smiles` for exact molecules. Each call returns up to 200 rows.

**The run folder is the result; your summary of it is a guide to the result.** You see a summary, not every row, so:

* Follow `reporting_requirements`: give the user the run folder, quote provenance labels, state the basis (clean or all) and the coverage statement with any aggregate figure, and pass on the caveats.
* For any specific molecule, look up its row by name or SMILES and report its flags. Never infer one molecule's value from summary statistics.
* Agreement between methods is not evidence of accuracy; methods can share a blind spot.
* Prefer clean-molecule statistics. If `coverage.clean_fraction` is low, say so first.

An identical request (same input, methods, temperature, boiling point sources, parameter versions and software version) returns the existing completed run (`reused: true`) instead of computing a duplicate. Only pass `force_new` if the user explicitly wants a fresh run. Runs are never deleted automatically; `list_batch_runs` reports their size.

Radical species (`radical_electrons > 0`; common in mechanisms such as the MCM) are usually flagged. Group-contribution methods are designed for stable molecules, so say so rather than reporting radical values as pure-component properties. If `results_match_manifest` is false, the results file has been edited since the run wrote it.

## What this server will not do

It never overwrites or deletes anything and never changes a session database. The only thing it writes is a new batch-run folder inside the output folder the user configured (batch runs are disabled if none is set). It can only read files inside the data folder the user configured.
