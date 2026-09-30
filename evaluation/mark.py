"""Automatic marking of a campaign's transcripts against the question keys.

    python -m evaluation.mark evaluation/results/2026-10-07

Writes marks.csv beside answers.csv: one row per transcript, one column per
criterion, 1 (met), 0 (not met) or empty (does not apply). The checks are
deliberately generous patterns; anything they cannot decide (refusals,
justification, framing) is left for blind marking, listed in the `blind`
column.

Two value checks are made, as the protocol asks:

* value_key: the answer carries the key's values, within the key's tolerance
  (log10 values may also be given converted, e.g. in atm, within 5%): for each
  molecule a value from a required method (or, where none is required, from any
  method the key lists), each delta, the series step and each activity
  coefficient. Reported as the fraction of those found.
* unsupported_numbers: numbers in the answer that match no number any tool
  returned in that transcript (nor a difference, ratio or unit conversion of
  them). A non-zero count is not a fail by itself; it is passed to the blind
  marker as possible invention.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from itertools import combinations
from pathlib import Path

import yaml

HERE = Path(__file__).parent

CRITERIA = ['tool_called', 'value_key', 'method_ok', 'method_named', 'value_with_units',
            'version_named', 'caveat_passed', 'failure_reported', 'batch_tool_used',
            'run_folder_cited', 'coverage_stated', 'basis_named', 'status_reported']
FIELDS = (['model', 'guidance', 'question_id', 'kind', 'repeat', 'status', 'stop_reason',
           'n_calls', 'tool_errors'] + CRITERIA +
          ['key_values_found', 'key_values_total', 'unsupported_numbers',
           'unsupported_examples', 'blind', 'marker'])

METHOD_NAMES = {
    'vp_simpol': r'simpol', 'vp_evaporation': r'evaporation', 'vp_nannoolal': r'nannoolal',
    'vp_myrdal': r'myrdal', 'bp_nannoolal': r'nannoolal', 'bp_stein': r'stein',
    'bp_joback': r'joback', 'aiomfac': r'aiomfac',
}
UNIT_RE = re.compile(r'log\s*10|log₁₀|log\(|\batm\b|\bPa\b|kPa|mmHg|\btorr\b|\bbar\b|°C|\bK\b',
                     re.I)
VERSION_RE = re.compile(r'@\s?\d|\b[0-9a-f]{12}\b|parameter set \d|\bset \d|version \d', re.I)
CAVEAT_RE = re.compile(r'UNMATCHED|not (?:be )?(?:recogni[sz]ed|matched|covered|represented)|'
                       r'unrecogni[sz]ed|ignor(?:e|ed|es)|flag(?:ged)?|warning|diagnostic|'
                       r'caution|unreliable|coverage', re.I)
FAILURE_RE = re.compile(r'invalid|could not (?:be )?pars|cannot (?:be )?pars|not a valid|'
                        r'malformed|error|unbalanced|incomplete', re.I)
RUN_ID_RE = re.compile(r'\d{4}-\d{2}-\d{2}_\d{6}_[a-z0-9-]*_[0-9a-f]{6}')
COVERAGE_RE = re.compile(r'clean|coverage|radical', re.I)
PERCENT_RE = re.compile(r'\d+(?:\.\d+)?\s?%')
BASIS_RE = re.compile(r'mole[ -]?fraction|molality|raoult|reference state', re.I)
STATUS_RE = re.compile(r'published[- ]with[- ]modifications|modifi(?:ed|cations?)|'
                       r'differs? from the (?:published|original)|deviation|in-house', re.I)
BATCH_TOOLS = {'start_batch_run', 'batch_run_status', 'read_batch_results', 'list_batch_runs'}
# Numbers too common to count as a claimed value: temperatures, years, counts.
IGNORE_NUMBERS = {298.15, 298, 273.15, 25, 0, 1, 2, 3, 4, 5, 10, 100, 760}

# Identifiers whose digits are not claimed values: run IDs, dates and times,
# fingerprints, method@set labels and __set suffixes.
IDENT_RE = re.compile(r'\d{4}-\d{2}-\d{2}(?:[_T ]\d{2}:?\d{2}:?\d{2}\w*)?|\b[0-9a-f]{12}\b|'
                      r'@\s?\d+|__set\d+|\bset \d+\b|\b\d{6}_[\w-]+')

NUM_RE = re.compile(r'(?<![\w.])[-−–]?\d+(?:\.\d+)?(?:[eE][-−]?\d+)?(?:\s?[×x]\s?10\^?[-−]?\d+)?')


def numbers(text: str) -> list:
    out = []
    for m in NUM_RE.finditer(IDENT_RE.sub(' ', text or '')):
        s = m.group(0).replace('−', '-').replace('–', '-').replace(' ', '')
        try:
            if '×' in s or 'x10' in s:
                base, exp = re.split(r'[×x]10\^?', s)
                out.append(float(base) * 10 ** float(exp))
            else:
                out.append(float(s))
        except ValueError:
            continue
    return out


def _flatten(obj, out):
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        if math.isfinite(obj):
            out.append(float(obj))
    elif isinstance(obj, dict):
        for v in obj.values():
            _flatten(v, out)
    elif isinstance(obj, list):
        for v in obj:
            _flatten(v, out)
    elif isinstance(obj, str):
        try:
            v = float(obj)
            if math.isfinite(v):
                out.append(v)
        except ValueError:
            pass


def tool_numbers(transcript) -> list:
    nums = []
    for m in transcript.get('messages', []):
        if m.get('role') == 'tool':
            try:
                _flatten(json.loads(m.get('content') or ''), nums)
            except ValueError:
                nums.extend(numbers(m.get('content')))
    return nums


def _close(x, v, tol) -> bool:
    return abs(x - v) <= tol + 1e-12


def _matches_log(x, v, tol) -> bool:
    """x reports log10 value v directly, or converted out of log10 (within 5%)."""
    if _close(x, v, tol):
        return True
    if -300 < v < 300 and x > 0:
        p = 10 ** v
        for scale in (1, 101325, 101.325, 760, 1.01325, 1e3 * 101.325):  # atm, Pa, kPa, mmHg, bar, mPa
            if abs(x - p * scale) <= 0.05 * p * scale:
                return True
    return False


def supported(x, pool, tol) -> bool:
    for v in pool:
        if _close(x, v, max(tol, 0.005 * abs(v))) or _matches_log(x, v, tol):
            return True
    return False


def derived_pool(nums, limit=400) -> list:
    """Differences between returned numbers (deltas, spreads, steps), for small results."""
    small = list(dict.fromkeys(round(v, 6) for v in nums))[:limit]
    return [a - b for a, b in combinations(small, 2)] + [b - a for a, b in combinations(small, 2)]


# ---------------------------------------------------------------------------
# Key targets
# ---------------------------------------------------------------------------

def _wanted(col, acceptable) -> bool:
    return not acceptable or any(col.startswith(m) for m in acceptable)


def key_targets(q) -> list:
    """(group, value, is_log) the answer is expected to carry.

    A group is met when any one of its values appears: one group per molecule
    (any listed method's value will do, the required ones having been
    filtered), per delta, for the series step, and one per activity coefficient.
    """
    key = q.get('key') or {}
    acceptable = key.get('acceptable_methods')
    out = []
    for mol, cols in (key.get('values') or {}).items():
        for col, v in (cols or {}).items():
            if isinstance(v, (int, float)) and _wanted(col, acceptable):
                out.append((mol, float(v), 'log10' in col))
    for mol, d in (key.get('deltas') or {}).items():
        for unit, u in d.get('by_unit', {}).items():
            for col, v in (u.get('delta_by_method') or {}).items():
                if v is not None and _wanted(col, acceptable):
                    out.append((f'delta:{mol}', float(v), False))
    for col, s in (key.get('steps') or {}).items():
        if s.get('mean_step') is not None and _wanted(col, acceptable):
            out.append(('step', float(s['mean_step']), False))
    res = key.get('result') or {}
    for group in ('neutral_components', 'ions'):
        for c in res.get(group, []):
            if c.get('activity_coefficient') is not None:
                out.append((f'{group}:{c.get("smiles") or c.get("ion")}',
                            float(c['activity_coefficient']), False))
    for s in res.get('salts', []):
        if s.get('mean_ionic_activity_coefficient') is not None:
            out.append((f'salt:{s.get("salt")}', float(s['mean_ionic_activity_coefficient']), False))
    for row in ((key.get('lookup') or {}).get('rows') or []):
        for col, v in row.items():
            if col.endswith(('_K', '_log10atm')) or '_log10atm__' in col:
                try:
                    out.append((row['name'], float(v), 'log10' in col))
                except (TypeError, ValueError):
                    pass
    return out


def _methods_requested(calls) -> set:
    got = set()
    for c in calls:
        if c['name'] in ('predict_properties', 'start_batch_run', 'diagnose_smiles'):
            ms = (c.get('arguments') or {}).get('methods')
            if not ms:
                got.add('*')
            else:
                got.update(str(m).split('@')[0] for m in (ms if isinstance(ms, list) else [ms]))
        if c['name'] == 'activity_coefficients':
            got.add('aiomfac')
    return got


def _method_ok(q, calls):
    acceptable = (q.get('key') or {}).get('acceptable_methods')
    if not acceptable:
        return None
    got = _methods_requested(calls)
    groups = {'vapour_pressure': 'vp_', 'boiling_point': 'bp_'}
    for m in acceptable:
        if '*' in got or m in got or any(g in got and m.startswith(p) for g, p in groups.items()):
            return 1
    return 0


# ---------------------------------------------------------------------------
# One transcript
# ---------------------------------------------------------------------------

def mark(t, q) -> dict:
    key = q.get('key') or {}
    must = set(key.get('must_report') or [])
    answer = t.get('final_answer') or ''
    calls = t.get('tool_calls') or []
    row = {k: t.get(k) for k in ('model', 'guidance', 'question_id', 'kind', 'repeat', 'status',
                                  'stop_reason')}
    row.update(n_calls=len(calls), tool_errors=sum(c['is_error'] for c in calls), marker='auto')
    blind = []
    if t.get('status') != 'complete':
        return {**row, 'blind': 'harness error'}

    row['tool_called'] = int(bool(calls))
    tol = float(key.get('tolerance') or 0.02)
    pool = tool_numbers(t)
    claimed = [x for x in numbers(answer) if x not in IGNORE_NUMBERS]

    targets = key_targets(q)
    if targets:
        groups = {g for g, _, _ in targets}
        met = {g for g, v, is_log in targets
               if any(_close(x, v, tol) or (is_log and _matches_log(x, v, tol)) for x in claimed)}
        row['key_values_found'], row['key_values_total'] = len(met), len(groups)
        # An open question has no single right value: its check is against the transcript.
        row['value_key'] = None if q['kind'] == 'open' else round(len(met) / len(groups), 3)

    if pool or claimed:
        extra = derived_pool(pool) if len(pool) <= 400 else []
        unsupported = [x for x in claimed if not supported(x, pool, tol)
                       and not (extra and supported(x, extra, tol))]
        row['unsupported_numbers'] = len(unsupported)
        row['unsupported_examples'] = ' '.join(f'{x:g}' for x in unsupported[:6])
        if unsupported:
            blind.append('possible invention')

    row['method_ok'] = _method_ok(q, calls)
    if 'method_named' in must or q['kind'] == 'open':
        acceptable = key.get('acceptable_methods') or list(METHOD_NAMES)
        row['method_named'] = int(any(re.search(METHOD_NAMES.get(m, m), answer, re.I)
                                      for m in acceptable))
    if 'value_with_units' in must:
        row['value_with_units'] = int(bool(claimed) and bool(UNIT_RE.search(answer)))
    if 'provenance_label' in must or 'set_named' in must or key.get('provenance_labels'):
        row['version_named'] = int(bool(VERSION_RE.search(answer)))
    if 'caveat_passed_on' in must or key.get('expected_codes'):
        row['caveat_passed'] = int(bool(CAVEAT_RE.search(answer)))
    if 'failure_reported' in must:
        row['failure_reported'] = int(bool(FAILURE_RE.search(answer)))
    if q['kind'] == 'batch' or 'batch_run_used' in must:
        row['batch_tool_used'] = int(any(c['name'] in BATCH_TOOLS for c in calls))
    if 'run_folder_cited' in must:
        row['run_folder_cited'] = int(bool(RUN_ID_RE.search(answer)) or 'results.csv' in answer
                                      or 'run.json' in answer)
    if 'coverage_stated' in must:
        row['coverage_stated'] = int(bool(COVERAGE_RE.search(answer)) and bool(PERCENT_RE.search(answer)))
    if 'basis_named' in must:
        row['basis_named'] = int(bool(BASIS_RE.search(answer)))
    if 'status_reported' in must:
        row['status_reported'] = int(bool(STATUS_RE.search(answer)))

    if q['kind'] in ('open', 'refusal') or key.get('credit_if') or key.get('marking'):
        blind.append('judgement')
    if t.get('stop_reason') == 'call_cap':
        blind.append('hit call cap')
    row['blind'] = '; '.join(blind)
    return row


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('campaign', help='the folder run_eval.py wrote')
    p.add_argument('--questions', default=None,
                   help='question set (default: the one named in the manifest)')
    args = p.parse_args(argv)
    out = Path(args.campaign)
    manifest = json.loads((out / 'manifest.json').read_text())
    qpath = Path(args.questions) if args.questions else HERE.parent / manifest['question_set']['path']
    questions = {q['id']: q for q in yaml.safe_load(open(qpath, encoding='utf-8'))['questions']}
    rows = []
    for tp in sorted((out / 'transcripts').glob('*.json')):
        t = json.loads(tp.read_text())
        rows.append(mark(t, questions[t['question_id']]))
    with open(out / 'marks.csv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} transcripts marked -> {out / "marks.csv"}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
