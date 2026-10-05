# Contributing to AeroSysProp

AeroSysProp is a community tool. It is not tied to any single institution,
and its usefulness depends on people outside the original team reporting
problems, checking methods against the literature and adding new ones.
Contributions of every size are welcome, from a typo fix to a new method.

By taking part you agree to follow our [code of conduct](CODE_OF_CONDUCT.md).

## Ways to contribute

**Report a problem.** Open an issue with the SMILES, the method, what you
expected and what you got. The output of the diagnostics
(`aerosysprop.diagnostics.diagnose`) is often the quickest way to show what
went wrong.

**Report a group assignment.** If a molecule seems to be split into the
wrong functional groups, use the *Group assignment discrepancy* issue
template. Include the expected groups and the page or table of the source
paper that supports them, so the report can be turned into a test.

**Check a method against its source.** Every method card in
`aerosysprop/methods/cards/` lists its references, each marked `verified`
once someone has checked it, along with documented departures from the
paper. Checking a reference, or adding a worked example from a paper as a
test, is one of the most useful things a new contributor can do.

**Improve the code or documentation.** Bug fixes, tests, documentation,
the web app and the agent interface all follow the ordinary process below.

**Add or change a method.** See [below](#changes-to-methods-and-parameters).

## The ordinary process

1. For anything beyond a small fix, open an issue first so the change can
   be discussed before you spend time on it.
2. Fork the repository and create a branch from `main`.
3. Set up a development environment:

   ```bash
   pip install -e ".[dev]"
   ```

4. Make your change, with tests where it changes behaviour.
5. Run the checks that continuous integration runs:

   ```bash
   python -m aerosysprop.method_cards
   pytest aerosysprop/tests
   ```

6. Open a pull request and complete the template.

Please use UK spelling in user-facing text. `AGENTS.md` sets out the
conventions for the codebase, and applies to human contributors as much as
to AI coding agents.

## Changes to methods and parameters

A user must always be able to tell which technique, with which parameters,
produced a number. Changes to the files under `aerosysprop/methods/`,
`aerosysprop/data/` and `aerosysprop/groups.py` therefore follow the
process in [METHODS_GOVERNANCE.md](METHODS_GOVERNANCE.md): a new parameter
set version, an updated method card with new fingerprints, regenerated
reference values, the evidence for the change, and approval from the core
team.

Earlier parameter sets are kept, so past results stay reproducible. A
proposal that changes results is not a problem in itself; a change that
cannot be traced is.

## Adding a method

A new single-component method needs an implementation in
`aerosysprop/methods/`, its SMARTS and coefficient files in
`aerosysprop/data/`, one entry in the `PROPERTY_METHODS` registry in
`aerosysprop/batch.py`, a method card and tests. A new mixture method is an
entry in `aerosysprop/mixture_methods.py`. Once registered, a method appears
in the web app, the in-app documentation and the agent API without further
changes. The *Add a method* page in the web app walks through the steps.

## Joining the core team

The core team reviews changes to methods and parameters and decides open
questions recorded in the method cards. Membership is listed in
`METHODS_GOVERNANCE.md` and `docs/CODEOWNERS`. We aim for a team drawn from
several institutions. Regular contributors who have reviewed or improved
methods are invited to join, and anyone interested is welcome to say so in
an issue.

## Credit

Contributors are acknowledged in the release notes. Those who make
substantial contributions to a release are invited to be listed as authors
in `CITATION.cff`, which sets how the software is cited.
