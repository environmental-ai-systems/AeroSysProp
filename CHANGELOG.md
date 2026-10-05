# Changelog

All notable changes to AeroSysProp are recorded here. Version numbers
follow [semantic versioning](https://semver.org/). Every prediction is also
labelled with its method ID and parameter set version, recorded in the
method cards (see `METHODS_GOVERNANCE.md`), so a change to a method's
parameters is visible even when the package version is unchanged.

## [1.0.0] (unreleased)

First release of AeroSysProp as a successor to UManSysProp v1.0 (Topping
et al., 2016, https://doi.org/10.5194/gmd-9-899-2016). The changes below
are relative to UManSysProp.

### Core library

- Replaced OpenBabel and Pybel with RDKit for SMILES parsing and SMARTS
  matching, and removed Python 2 code. Predictions were checked against
  stored reference values throughout the migration.
- Loaders for plain text, comma or tab separated files, SQLite databases
  and SDF files, with automatic detection of the SMILES column
  (`aerosysprop.io`).
- Parallel batch processing of large compound sets with a progress bar
  and CSV output (`aerosysprop.batch`).
- Configurable boiling point source for the vapour pressure, critical
  property and density methods, including several sources at once to
  test sensitivity to that choice.
- Registries for single-component methods (`PropertyMethod`) and mixture
  methods (`aerosysprop.mixture_methods`), so a new method appears in the
  app, the documentation page and the API without changes elsewhere.

### Activity coefficients

- AIOMFAC activity coefficients from SMILES: the short-range modified
  UNIFAC part, the middle-range organic-ion and ion-ion parts and the
  long-range Debye-Hückel part, for mixtures of organics, water and common
  atmospheric salts. Validated against the Fortran AIOMFAC reference for
  more than 1500 atmospheric and laboratory compounds.
- Optional bisulfate dissociation equilibrium (HSO4- to H+ and SO4--),
  with a solver protected against non-finite activity coefficients.
- Cached system objects and cache controls for repeated calls inside
  numerical models, with guidance on throughput and thread safety.

### Diagnostics and inspection

- Diagnostics that flag problems the methods would otherwise pass over
  silently (`aerosysprop.diagnostics`): atoms or functional groups that no
  SMARTS pattern assigns, AIOMFAC subgroups missing from or outside the
  Fortran reference, fractional or negative group counts, missing
  interaction parameters for a mixture, and compositions that do not sum
  to one. Each SMILES or mixture can be checked before a calculation, or
  the issues reported as warnings or errors during it.
- SMARTS hit matrices for a batch of molecules across any pattern set,
  for checking how molecules are broken into groups
  (`aerosysprop.fingerprints`).
- An issue template for reporting incorrect group assignments.

### Method governance and provenance

- A versioned method card for every method (`aerosysprop/methods/cards/`)
  recording its references, publication status, documented departures
  from the source papers, known limitations, validation and review state.
- A SHA-256 fingerprint of every parameter file, checked against the cards
  in continuous integration, so any change to SMARTS patterns or
  coefficients fails the build until the card is updated.
- A prediction is identified by its method ID, parameter set version and
  software version with git commit.
- Support for more than one parameter set per method. EVAPORATION now has
  two: set 1, the original port, and set 2 (the default), which revises
  the alcohol terms to match the published method. `vp_evaporation2` is
  kept as a deprecated alias for set 2.
- A documented process for proposing, reviewing and recording method
  changes (`METHODS_GOVERNANCE.md`), with code owners and a pull request
  template.

### Web application

- A Streamlit multipage app with pages for single-component properties,
  mixture activity coefficients, external services, in-app documentation
  and a guide to adding a method.
- A session database (SQLite) that keeps compounds, single-component
  results and mixture calculations between pages and between sessions,
  and can be reloaded or exported to CSV (`aerosysprop.db`).
- File upload with downloadable example files, and diagnostics shown
  alongside results.

### Containers and external tools

- A Docker image and Compose file for running the web app with one
  command.
- Optional integration with OPERA as a separately installed or
  containerised service, with an asynchronous job API, an endpoint
  selector and documentation. AeroSysProp does not bundle OPERA and works
  fully without it.

### Interfaces for AI agents

- An agent-facing API that takes and returns plain values with explicit
  units, returns diagnostics with every value and reports errors with a
  stable code and a hint for fixing the call (`aerosysprop.api`).
- A read-only Model Context Protocol server (`aerosysprop-mcp`) over
  stdio or HTTP. File access is confined to a configured data folder.
- Batch runs for large molecule sets that write a run folder with
  results, diagnostics, a provenance manifest with file fingerprints and
  RO-Crate metadata, so the folder can be archived as it stands
  (`aerosysprop.runs`).
- A guidance level switch (`AEROSYSPROP_GUIDANCE`) that changes only the
  explanatory text served to agents, never the values.

### Testing

- Regression tests for every property method, the IO loaders, AIOMFAC
  (non-electrolyte and electrolyte), the EVAPORATION parameter sets, the
  API, the MCP server, batch runs and the method cards.
- Continuous integration on Python 3.10 and 3.12.

### Known issues

- EVAPORATION parameter set 2 also changes 13 alpha-hydroxy acids in the
  MCM by about 0.25 log units, which appears to be unintended. It is
  under review by the core team (see `METHODS_GOVERNANCE.md`).
- Two legacy test files use exact floating point comparisons that can
  fail in the last digit on some platforms. They run in CI but do not
  block it.
