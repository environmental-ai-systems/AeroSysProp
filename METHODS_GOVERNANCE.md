# Governance of prediction methods

AeroSysProp's value depends on users knowing exactly which technique, with which parameters, produced a number. Several of our methods have no single published version, carry additions for atmospheric compounds, or have had their parameters revised while the model structure stayed the same. This document sets out how we record that and who decides.

## Core team

Changes to prediction methods and their parameters are the responsibility of the core management team:

* @simonom

Membership is listed in `.github/CODEOWNERS`, which makes core team approval a requirement for merging any change to `aerosysprop/methods/`, `aerosysprop/data/` or `aerosysprop/groups.py`.

## What is versioned

Every prediction is identified by three things.

1. **Method ID** (for example `vp_evaporation`). This stays the same as long as the model equations are unchanged. A change to the equations is a new method with a new ID.
2. **Parameter set version** (an integer). This covers the SMARTS patterns that assign functional groups and the coefficient files. Any change that can alter a result, including a single SMARTS pattern, requires a new parameter set version. Earlier parameter sets are kept so that past results can be reproduced.
3. **Software version** (package version plus git commit). This covers changes to code, including group-counting logic in `groups.py` and coefficients held in code rather than data files.

Every result from the agent API carries a provenance label such as `vp_simpol@1 (3f9a1c2b7d40)`, together with the software version.

## Method cards

Each method has a card in `aerosysprop/methods/cards/<method_id>.yaml` recording:

* **status**: `published` (faithful implementation of a published method), `published-with-modifications` (published method with documented additions or interpretations) or `in-house` (no published version);
* **references** with DOIs, each marked `verified` once a team member has checked it;
* **parameter set**: its version, the files it comprises and their SHA-256 fingerprints, and where it came from;
* **history** of earlier parameter sets;
* **implementation notes**: every known deviation from the source paper;
* **group interpretation**: places where the SMARTS encode a judgement about a paper's group definitions rather than a literal transcription;
* **known limitations** and **validation** evidence;
* **review** state (`draft` or `approved`), approver and date.

## Fingerprints

A SHA-256 fingerprint is a 64-character code computed from a file's exact contents. Any change to the file, however small, gives a completely different code. The fingerprint of each parameter file is recorded in its card when the parameter set is approved, and the test suite recomputes it on every run. If a file no longer matches its card, the tests fail and name the file. This stops parameters changing silently.

## Changing a method's parameters

1. Open a pull request describing the change and the evidence for it (for example, comparison with measurements or with the source paper's worked examples).
2. In the method card, move the current `parameter_set` into `history`, increase `version` by one, describe the change under `origin`, and set the fingerprints of changed files to `null`.
3. Run `python -m aerosysprop.method_cards --record` to fill in the new fingerprints, then `python -m aerosysprop.method_cards` to check every card.
4. Regenerate the affected reference outputs in the test fixtures and explain the differences in the pull request.
5. A core team member reviews and, on approval, sets `review.state: approved`, `approved_by` and `approved_on`.

Fixes that do not change any intended value (comments, formatting of code) need no new parameter version, but still need review if they touch the protected paths.

## Group assignment and community reports

Most published methods define functional groups in words and tables; our SMARTS patterns are an interpretation. Two kinds of error are possible. Atoms that no group matches are detected automatically and reported as `UNMATCHED_ATOMS`. Atoms assigned to the wrong group, or counted twice, are silent, so we rely on tests and on the community.

Anyone who believes a molecule is assigned incorrectly should open a "Group assignment discrepancy" issue with the molecule, the expected groups and the page or table of the paper that supports them. The core team decides whether the report leads to a new parameter set, a documented interpretation in the card, or no change, and records the outcome in the card.

Where a source paper gives worked examples or group assignments for specific molecules, these should become group-count tests, so that disputed interpretations are settled by tests rather than by discussion.

## Open items for the core team

The draft cards record the following for decision:

* Whether `vp_evaporation2` should become parameter set 2 of `vp_evaporation` (a breaking change to method IDs and output columns).
* The Schroeder and Le Bas density temperature correction always uses the Nannoolal boiling point, while the critical properties follow the chosen boiling point source.
* The origin of the atmospheric additions to the Nannoolal and Stein & Brown group sets.
* The critical property citation (labelled Nannoolal et al. 2004; the method appears to be Nannoolal, Rarey and Ramjugernath 2007).
* Verification of every reference and DOI, and the default status of each method.
