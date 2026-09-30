# This file is part of AeroSysProp.
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
"""
aerosysprop.method_cards
========================
Method cards: the versioned record of every prediction technique.

Each method has a YAML card in ``aerosysprop/methods/cards/`` recording its
publication status, references, the parameter set in use (with a SHA-256
fingerprint of every parameter file), documented deviations from the source
papers, and review state. See ``METHODS_GOVERNANCE.md`` for the process.

A prediction is identified by three things:

* the method ID (stable while the model equations are unchanged);
* the parameter set version (increases whenever SMARTS patterns or
  coefficients change); and
* the software version (package version and git commit), which covers
  changes to the code itself.

Command line, for the core team::

    python -m aerosysprop.method_cards            # check every card
    python -m aerosysprop.method_cards --record   # fill in fingerprints left as null
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

import yaml

from aerosysprop import __version__

PACKAGE_DIR = Path(__file__).parent
CARDS_DIR = PACKAGE_DIR / 'methods' / 'cards'

STATUSES = ('published', 'published-with-modifications', 'in-house')
REVIEW_STATES = ('draft', 'approved')
REQUIRED_FIELDS = ('id', 'name', 'property', 'status', 'model', 'review')


def parameter_sets(card: dict) -> list:
    """All parameter sets a card offers. A card has either ``parameter_set``
    (one set) or ``parameter_sets`` (several, exactly one marked default)."""
    if 'parameter_sets' in card:
        return list(card['parameter_sets'])
    return [card['parameter_set']] if 'parameter_set' in card else []


def default_set(card: dict) -> dict:
    sets = parameter_sets(card)
    for ps in sets:
        if ps.get('default'):
            return ps
    return sets[0] if sets else {}


def select_set(card: dict, version=None) -> Optional[dict]:
    if version is None:
        return default_set(card)
    for ps in parameter_sets(card):
        if ps.get('version') == version:
            return ps
    return None


def _split_key(key: str):
    """'vp_evaporation@2' -> ('vp_evaporation', 2); 'vp_simpol' -> ('vp_simpol', None)."""
    if '@' in key:
        mid, _, v = key.partition('@')
        return mid, int(v) if v.isdigit() else v
    return key, None


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def load_cards() -> dict:
    """Return all method cards keyed by method ID."""
    cards = {}
    for path in sorted(CARDS_DIR.glob('*.yaml')):
        card = yaml.safe_load(path.read_text(encoding='utf-8'))
        card['_path'] = path
        cards[card['id']] = card
    return cards


def card_for(method_id: str):
    return load_cards().get(method_id)


# ---------------------------------------------------------------------------
# Fingerprints
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def file_sha256(relpath: str) -> str:
    """SHA-256 of a parameter file, given relative to the package directory.

    Line endings are normalised so a checkout on Windows gives the same
    fingerprint as one on Linux or macOS.
    """
    raw = (PACKAGE_DIR / relpath).read_bytes().replace(b'\r\n', b'\n')
    return hashlib.sha256(raw).hexdigest()


def parameter_set_fingerprint(files: dict) -> str:
    """One fingerprint for a whole parameter set, combining its files."""
    h = hashlib.sha256()
    for rel in sorted(files):
        h.update(f'{rel}:{files[rel]}\n'.encode())
    return h.hexdigest()


def check_card(card: dict) -> list:
    """Return a list of problems with one card (empty if it is sound)."""
    problems = []
    cid = card.get('id', card.get('_path', '?'))
    for f in REQUIRED_FIELDS:
        if f not in card:
            problems.append(f'{cid}: missing field {f!r}')
    if card.get('status') not in STATUSES:
        problems.append(f'{cid}: status must be one of {STATUSES}')
    if (card.get('review') or {}).get('state') not in REVIEW_STATES:
        problems.append(f'{cid}: review.state must be one of {REVIEW_STATES}')
    sets = parameter_sets(card)
    if not sets:
        problems.append(f'{cid}: missing parameter_set or parameter_sets')
    if 'parameter_sets' in card and sum(1 for ps in sets if ps.get('default')) != 1:
        problems.append(f'{cid}: exactly one of parameter_sets must be marked default: true')
    versions = [ps.get('version') for ps in sets]
    if len(set(map(str, versions))) != len(versions):
        problems.append(f'{cid}: parameter set versions must be unique')
    for pset in sets:
        if not isinstance(pset.get('version'), int):
            problems.append(f'{cid}: parameter set version must be an integer')
        for rel, recorded in (pset.get('files') or {}).items():
            path = PACKAGE_DIR / rel
            if not path.is_file():
                problems.append(f'{cid}: parameter file {rel} does not exist')
            elif recorded is None:
                problems.append(f'{cid}: no fingerprint recorded for {rel}')
            elif file_sha256(rel) != recorded:
                problems.append(
                    f'{cid}: {rel} has changed since parameter set '
                    f'{pset.get("version")} was recorded. If the change is '
                    'intended, declare a new parameter set version in the card.')
    return problems


def check_cards() -> list:
    problems = []
    for card in load_cards().values():
        problems.extend(check_card(card))
    return problems


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def _git_commit():
    """Short commit hash when running from a git checkout, else None."""
    repo = PACKAGE_DIR.parent
    if not (repo / '.git').exists():
        return None
    try:
        out = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, timeout=2)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def software_provenance() -> dict:
    return {'aerosysprop_version': __version__, 'git_commit': _git_commit()}


def provenance(method_id: str) -> dict:
    """Provenance record for one method, suitable for attaching to results.

    ``method_id`` may name a parameter set, as in ``'vp_evaporation@1'``;
    otherwise the card's default set is used.
    """
    mid, version = _split_key(method_id)
    card = card_for(mid)
    if card is None:
        return {'method': method_id, 'documented': False}
    pset = select_set(card, version)
    if pset is None:
        return {'method': method_id, 'documented': False,
                'warning': f'{mid} has no parameter set {version!r} in its card.'}
    method_id = mid
    files = pset.get('files') or {}
    current = {rel: file_sha256(rel) for rel in files if (PACKAGE_DIR / rel).is_file()}
    matches = bool(files) and all(files[r] == current.get(r) for r in files)
    fp = parameter_set_fingerprint(files) if files else None
    out = {
        'method': method_id,
        'documented': True,
        'status': card['status'],
        'parameter_set_version': pset['version'],
        'parameter_fingerprint': fp[:12] if fp else None,
        'parameters_match_card': matches,
        'review_state': card['review']['state'],
        'label': f'{method_id}@{pset["version"]}' + (f' ({fp[:12]})' if fp else ''),
    }
    if not matches and files:
        out['warning'] = ('Parameter files differ from the version recorded in '
                          'the method card; results may not be reproducible.')
    return out


def summary(method_id: str) -> dict:
    """Card information for method catalogues (e.g. the agent API)."""
    card = card_for(method_id)
    if card is None:
        return {'documented': False}
    refs = [r.get('citation') for r in card['model'].get('references', [])]
    sets = parameter_sets(card)
    return {
        'documented': True,
        'status': card['status'],
        'references': refs,
        'parameter_set_version': default_set(card).get('version'),
        **({'parameter_sets': [{'version': ps['version'], 'default': bool(ps.get('default')),
                                'description': ps.get('description') or ps.get('origin')}
                               for ps in sets]} if len(sets) > 1 else {}),
        'known_limitations': card.get('known_limitations', []),
        'deviations_from_source': card.get('implementation_notes', []),
        'review_state': card['review']['state'],
    }


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def _record_missing() -> list:
    """Fill fingerprints written as ``null`` in the cards, preserving the rest
    of each file (comments included). Existing fingerprints are never changed:
    to change one, declare a new parameter set version and set it to null."""
    changed = []
    for card in load_cards().values():
        path = card['_path']
        text = path.read_text(encoding='utf-8')
        new = text
        files = {rel: rec for ps in parameter_sets(card) for rel, rec in (ps.get('files') or {}).items()}
        for rel, recorded in files.items():
            if recorded is None and (PACKAGE_DIR / rel).is_file():
                pattern = re.compile(r'^(\s+' + re.escape(rel) + r':\s*)null\s*$', re.M)
                new = pattern.sub(lambda m: m.group(1) + file_sha256(rel), new, count=1)
        if new != text:
            path.write_text(new, encoding='utf-8')
            changed.append(path.name)
    load_cards.cache_clear()
    return changed


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog='python -m aerosysprop.method_cards',
                                     description='Check or record method card fingerprints.')
    parser.add_argument('--record', action='store_true',
                        help='Fill in fingerprints that are null in the cards.')
    args = parser.parse_args(argv)
    if args.record:
        for name in _record_missing():
            print(f'recorded fingerprints in {name}')
    problems = check_cards()
    for p in problems:
        print('PROBLEM:', p)
    cards = load_cards()
    print(f'{len(cards)} cards checked, {len(problems)} problem(s).')
    for cid, card in cards.items():
        for ps in parameter_sets(card):
            key = f'{cid}@{ps["version"]}' if len(parameter_sets(card)) > 1 else cid
            prov = provenance(key)
            mark = ' (default)' if ps.get('default') else ''
            print(f'  {prov["label"] + mark:45s} {card["status"]:30s} review: {card["review"]["state"]}')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
