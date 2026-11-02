"""Automatic marking of a campaign's transcripts against the question keys.

    python -m evaluation.mark evaluation/results/2026-10-07

Writes marks.csv beside answers.csv: one row per transcript, one column per
criterion, 1 (met), 0 (not met) or empty (does not apply). Where the key holds
the answer (coverage figures, version labels, the reference basis, the run the
tools reported) a check compares against it rather than looking for a word.
Anything they cannot decide (refusals, justification, framing) is left for
blind marking, listed in the `blind` column.

Value checks, as the protocol asks:

* value_key: the answer carries the key's values, within the key's tolerance
  (log10 values may also be given converted, e.g. in atm, within 5%): for each
  molecule a value from a required method (or, where none is required, from any
  method the key lists), each delta, the series step and each activity
  coefficient. Reported as the fraction of those found.
* value_key_all: the same, but every expected value counts separately, so an
  answer giving one method of the three asked for scores a third.
* unsupported_numbers: numbers in the answer that match no number any tool
  returned in that transcript (nor a difference, ratio or unit conversion of
  them). A non-zero count is not a fail by itself; it is passed to the blind
  marker as possible invention. Whole numbers of four or more digits that
  appear anywhere in the text the model was given (the question, or a tool's
  output such as the method cards' references) also count as supported, so a
  correctly cited "Pankow & Asher 2008" is not flagged but an invented year is.
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

CRITERIA = ['tool_called', 'value_key', 'value_key_all', 'method_ok', 'method_named',
            'value_with_units', 'version_named', 'caveat_passed', 'no_false_warning',
            'failure_reported', 'batch_tool_used', 'run_folder_cited', 'coverage_stated',
            'basis_named', 'status_reported']
FIELDS = (['model', 'guidance', 'question_id', 'kind', 'repeat', 'status', 'stop_reason',
           'n_calls', 'tool_errors'] + CRITERIA +
          ['key_values_found', 'key_values_total', 'unsupported_numbers',
           'unsupported_examples', 'blind', 'marker'])

METHOD_NAMES = {
    'vp_simpol': r'simpol', 'vp_evaporation': r'evaporation', 'vp_nannoolal': r'nannoolal',
    'vp_myrdal': r'myrdal', 'bp_nannoolal': r'nannoolal', 'bp_stein': r'stein',
    'bp_joback': r'joback', 'aiomfac': r'aiomfac',
}
UNIT_RE = re.compile(r'log\s*10|log₁₀|log\(|\batm\b|\bPa\b|kPa|mmHg|\btorr\b|\bbar\b|°C|\bK\b|'
                     r'\bk?g\s*(?:/|per|\s)\s*(?:c?m|mL)', re.I)
VERSION_RE = re.compile(r'@\s?\d|\b[0-9a-f]{12}\b|parameter set \d|\bset \d|version \d', re.I)
# Phrases that pass a coverage warning on. Bare words such as "coverage" or
# "diagnostic" are not enough: "coverage is complete" is not a warning.
WARN_RE = re.compile(r'UNMATCHED|not (?:be )?(?:recogni[sz]ed|matched|covered|represented|'
                     r'accounted for)|unrecogni[sz]ed|ignor(?:e|ed|es|ing)\b|flagged|warning|'
                     r'caution|unreliable|less reliable|not reliable|'
                     r'(?:incomplete|partial|limited|poor) coverage|'
                     r'coverage (?:problem|issue|gap|warning|limitation)|missing groups?|'
                     r'no (?:group|parameter)s? for|'
                     r'not (?:designed|intended|parameteri[sz]ed|suitable|valid|meant) for|'
                     r'(?:does|do) not (?:typically )?(?:support|handle|cover)', re.I)
# A warning phrase in a sentence with one of these says nothing is wrong
# ("no warnings were raised"), so that sentence is not a warning.
# "any ... problems" is the question being restated or a hypothetical, not a finding.
NEGATED_RE = re.compile(r'\b(?:no|zero|any)\W+(?:[\w‑-]+\W+){0,2}(?:warnings?|flags?|diagnostics?|'
                        r'issues?|problems?|limitations?|gaps?|concerns?|unmatched)|'
                        r'\bwithout\W+(?:any\W+)?(?:[\w‑-]+\W+)?(?:warnings?|flags?|issues?|problems?)|'
                        r'\bnot flagged|\b(?:user|you) (?:asked|is asking|are asking)|'
                        r'\basking (?:if|whether|about)|'
                        r'\bnone (?:of \w+ )?(?:were|was|are|is) flagged|'
                        r'\bnone of (?:the|these) (?:\w+\W+){0,3}(?:raised|reported|flagged|returned)|'
                        r'\bfull coverage|'
                        r'\bno atoms? (?:were |was )?(?:left )?unmatched|'
                        r'\ball atoms (?:were |are )?matched|\bfully (?:covered|matched|recogni[sz]ed)|'
                        r'\bclean\b', re.I)
FAILURE_RE = re.compile(r'invalid|could not (?:be )?(?:pars|read|interpret)|'
                        r'cannot (?:be )?(?:pars|read|interpret)|not a valid|malformed|unbalanced|'
                        r'unclosed|(?:parse|parsing|syntax) error|'
                        r'incomplete (?:smiles|ring|branch|bracket|parenthes)|'
                        # a tool's service being down (opera-01)
                        r'not reachable|unreachable|unavailable|not running|'
                        r'could not (?:be )?(?:reach|connect)|cannot (?:be )?(?:reach|connect)', re.I)
RUN_ID_RE = re.compile(r'\d{4}-\d{2}-\d{2}_\d{6}_[a-z0-9-]*_[0-9a-f]{6}')
COVERAGE_RE = re.compile(r'clean|coverage|radical', re.I)
PERCENT_RE = re.compile(r'\d+(?:\.\d+)?\s?%')
BASIS_RE = re.compile(r'mole[ -]?fraction|molality|raoult|reference state', re.I)
NEUTRAL_BASIS_RE = re.compile(r'mole[ -]?fraction|raoult', re.I)
ION_BASIS_RE = re.compile(r'molality|molal\b|infinite dilution', re.I)
TEMP_RE = re.compile(r'298(?:\.15)?\s*K\b|25\s*°\s*C|273\.15\s*K\b', re.I)
LIST_ITEM_RE = re.compile(r'^\s*(?:[-*•+]|\d+[.)])\s')
FOLDER_RE = re.compile(r'runs/([0-9a-f]{12})\b')
STATUS_RE = re.compile(r'published[- ]with[- ]modifications|modifi(?:ed|cations?)|'
                       r'differs? from the (?:published|original)|deviation|in-house', re.I)
# Listing earlier runs is not using one.
BATCH_TOOLS = {'start_batch_run', 'batch_run_status', 'read_batch_results'}
# Numbers too common to count as a claimed value: temperatures, years, counts.
IGNORE_NUMBERS = {298.15, 298, 273.15, 25, 0, 1, 2, 3, 4, 5, 10, 100, 760}

# Identifiers whose digits are not claimed values: run IDs, dates and times,
# fingerprints, method@set labels and __set suffixes.
IDENT_RE = re.compile(r'\d{4}-\d{2}-\d{2}(?:[_T ]\d{2}:?\d{2}:?\d{2}\w*)?|\b[0-9a-f]{12}\b|'
                      r'@\s?\d+|__set\d+|\bset \d+\b|\b\d{6}_[\w-]+')

# A comma between digit groups of exactly three ("6,580 Pa") is a thousands
# separator, read as one number; "1,2" or "0.8,0.2" stay separate numbers.
NUM_RE = re.compile(r'(?<![\w.])[-−–]?(?:\d{1,3}(?:,\d{3})+(?![\d,])|\d+)(?:\.\d+)?'
                    r'(?:[eE][-−]?\d+)?(?:\s?[×x]\s?10\^?[-−]?\d+)?')


# Models write minus signs as any of these (gpt-oss uses the non-breaking hyphen).
MINUS_RE = re.compile('[\u2010\u2011\u2012\u2212\ufe63\uff0d]')


def numbers(text: str) -> list:
    out = []
    for m in NUM_RE.finditer(IDENT_RE.sub(' ', MINUS_RE.sub('-', text or ''))):
        s = m.group(0).replace('−', '-').replace('–', '-').replace(' ', '').replace(',', '')
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


# Citation years, page ranges and pattern numbers. Integers only, so an example
# value in a tool description cannot make an invented property value look supported.
TEXT_INT_RE = re.compile(r'(?<![\w.])\d{4,}(?![\w]|\.\d)')


def text_integers(transcript) -> set:
    texts = [transcript.get('question') or '']
    texts += [m.get('content') or '' for m in transcript.get('messages', [])
              if m.get('role') == 'tool']
    return {float(s) for t in texts for s in TEXT_INT_RE.findall(str(t))}


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
# Reporting checks
# ---------------------------------------------------------------------------

def _sentences(text) -> list:
    return [s for s in re.split(r'(?<=[.!?])\s+|\n+', text or '') if s.strip()]


def warns(text) -> bool:
    """A sentence carries a warning phrase and does not say that nothing is wrong."""
    return any(WARN_RE.search(s) and not NEGATED_RE.search(s) for s in _sentences(text))


def units_with_values(answer) -> bool:
    """A reported number has its unit on the same line, in its table's header, or
    in the line introducing its list. A restated temperature does not count."""
    lines = TEMP_RE.sub(' ', answer or '').splitlines()
    for i, line in enumerate(lines):
        if not any(x not in IGNORE_NUMBERS for x in numbers(line)):
            continue
        context, j = [line], i
        if line.lstrip().startswith('|'):
            while j > 0 and lines[j - 1].lstrip().startswith('|'):
                j -= 1
            context.append(lines[j])
        elif LIST_ITEM_RE.match(line):
            while j > 0 and (LIST_ITEM_RE.match(lines[j - 1]) or not lines[j - 1].strip()):
                j -= 1
            if j > 0:
                context.append(lines[j - 1])
        if any(UNIT_RE.search(c) for c in context):
            return True
    return False


def version_quoted(answer, key) -> bool:
    """The answer quotes one of the key's provenance labels: its fingerprint, the
    method with the right @version, or EVAPORATION's set number."""
    labels = key.get('provenance_labels') or {}
    if not labels:  # set questions whose key holds no label
        return bool(VERSION_RE.search(answer))
    low = answer.lower()
    for label in labels.values():
        m = re.match(r'(\w+)@(\d+)\s*\(([0-9a-f]{12})\)', label)
        if not m:
            continue
        mid, n, fp = m.groups()
        if fp in low:
            return True
        for name in (mid, METHOD_NAMES.get(mid, mid)):
            if re.search(rf'{name}[\w.]*`?\s*@\s*{n}\b', answer, re.I):
                return True
        if 'evaporation' in mid and re.search(rf'(?<![a-z])set\s*{n}\b', answer, re.I):
            return True
    return False


def coverage_correct(answer, key) -> bool:
    """The answer gives one of the run's coverage figures: a clean, flagged or
    radical share (overall or per method) within one percentage point, or the
    count out of the total."""
    cov = key.get('coverage') or {}
    if not cov:
        return bool(COVERAGE_RE.search(answer)) and bool(PERCENT_RE.search(answer))
    shares = [cov.get(k) for k in ('clean_fraction', 'flagged_fraction', 'radical_fraction')]
    for f in (cov.get('flagged_fraction_by_method') or {}).values():
        shares += [f, 1 - f]
    shares = [100 * f for f in shares if f is not None]
    for m in re.finditer(r'(\d+(?:\.\d+)?)\s?%', answer):
        if any(abs(float(m.group(1)) - f) <= 1.0 for f in shares):
            return True
    counts = {cov.get(k) for k in ('clean', 'flagged', 'radical_species')} - {None}
    for m in re.finditer(r'(\d+)\s*(?:of|/|out of)\s*(?:the\s*)?(\d+)', answer):
        if int(m.group(2)) == cov.get('molecules') and int(m.group(1)) in counts:
            return True
    return False


def basis_correct(answer, key) -> bool:
    """Mole fraction (Raoult) for neutral components, molality for ions."""
    res = key.get('result') or {}
    need = []
    if res.get('neutral_components'):
        need.append(NEUTRAL_BASIS_RE)
    if res.get('ions') or res.get('salts'):
        need.append(ION_BASIS_RE)
    if not need:
        return bool(BASIS_RE.search(answer))
    return all(r.search(answer) for r in need)


def runs_reported(transcript) -> set:
    """Run IDs and run folder names the tools reported in this conversation."""
    text = ' '.join(str(m.get('content') or '') for m in transcript.get('messages', [])
                    if m.get('role') == 'tool')
    return set(RUN_ID_RE.findall(text)) | set(FOLDER_RE.findall(text))


# ---------------------------------------------------------------------------
# Key targets
# ---------------------------------------------------------------------------

def _method_ids(names) -> set:
    """Method IDs for a list of method IDs and property groups ("vapour_pressure"),
    resolved the way the tools resolve them."""
    from aerosysprop import api, batch
    ids = set()
    for n in names:
        group = api._GROUP_ALIASES.get(api._slug(str(n)))
        ids.update(batch.METHOD_GROUPS[group] if group else [str(n).split('@')[0]])
    return ids


def _column_prefixes(names) -> tuple:
    """Output columns a list of methods or groups produces, e.g. critical_properties
    gives tc_nannoolal_K and pc_nannoolal_bar."""
    from aerosysprop import batch
    by_id = {m.id: m for m in batch.PROPERTY_METHODS}
    return tuple(k for i in _method_ids(names) for k in (by_id[i].output_keys if i in by_id else [i]))


def _wanted(col, acceptable) -> bool:
    return not acceptable or col.startswith(_column_prefixes(acceptable))


def key_targets(q) -> list:
    """(group, value, is_log, column) the answer is expected to carry.

    A group is met when any one of its values appears: one group per molecule
    (any listed method's value will do, the required ones having been
    filtered), per delta, for the series step, and one per activity coefficient.
    The column separates a group's values for value_key_all.
    """
    key = q.get('key') or {}
    acceptable = key.get('acceptable_methods')
    out = []
    for mol, cols in (key.get('values') or {}).items():
        for col, v in (cols or {}).items():
            if isinstance(v, (int, float)) and _wanted(col, acceptable):
                out.append((mol, float(v), 'log10' in col, col))
    for mol, d in (key.get('deltas') or {}).items():
        for unit, u in d.get('by_unit', {}).items():
            for col, v in (u.get('delta_by_method') or {}).items():
                if v is not None and _wanted(col, acceptable):
                    out.append((f'delta:{mol}', float(v), False, col))
    for col, s in (key.get('steps') or {}).items():
        if s.get('mean_step') is not None and _wanted(col, acceptable):
            out.append(('step', float(s['mean_step']), False, col))
    res = key.get('result') or {}
    for group in ('neutral_components', 'ions'):
        for c in res.get(group, []):
            if c.get('activity_coefficient') is not None:
                out.append((f'{group}:{c.get("smiles") or c.get("ion")}',
                            float(c['activity_coefficient']), False, ''))
    for s in res.get('salts', []):
        if s.get('mean_ionic_activity_coefficient') is not None:
            out.append((f'salt:{s.get("salt")}', float(s['mean_ionic_activity_coefficient']), False, ''))
    for row in ((key.get('lookup') or {}).get('rows') or []):
        for col, v in row.items():
            if col.endswith(('_K', '_log10atm')) or '_log10atm__' in col:
                try:
                    out.append((row['name'], float(v), 'log10' in col, col))
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
                got.update(_method_ids(ms if isinstance(ms, list) else [ms]))
        if c['name'] == 'activity_coefficients':
            got.add('aiomfac')
    return got


def _method_ok(q, calls):
    acceptable = (q.get('key') or {}).get('acceptable_methods')
    if not acceptable:
        return None
    got = _methods_requested(calls)
    return int('*' in got or bool(got & _method_ids(acceptable)))


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

    # A call the server rejected is not a tool result.
    row['tool_called'] = int(any(not c['is_error'] for c in calls))
    tol = float(key.get('tolerance') or 0.02)
    pool = tool_numbers(t)
    claimed = [x for x in numbers(answer) if x not in IGNORE_NUMBERS]

    targets = key_targets(q)
    if targets:
        # Values in K, bar or g/cm3 may be rounded: allow 0.5%, as supported() does.
        found = {(g, col) for g, v, is_log, col in targets
                 if any(_close(x, v, tol if is_log else max(tol, 0.005 * abs(v)))
                        or (is_log and _matches_log(x, v, tol)) for x in claimed)}
        groups = {g for g, _, _, _ in targets}
        met = {g for g, _ in found}
        every = {(g, col) for g, _, _, col in targets}
        row['key_values_found'], row['key_values_total'] = len(met), len(groups)
        # An open question has no single right value: its check is against the transcript.
        if q['kind'] != 'open':
            row['value_key'] = round(len(met) / len(groups), 3)
            row['value_key_all'] = round(len(found) / len(every), 3)

    if pool or claimed:
        extra = derived_pool(pool) if len(pool) <= 400 else []
        quoted = text_integers(t)
        unsupported = [x for x in claimed if x not in quoted and not supported(x, pool, tol)
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
        row['value_with_units'] = int(units_with_values(answer))
    if ('provenance_label' in must or 'set_named' in must or key.get('provenance_labels')) \
            and not key.get('invalid_smiles'):
        row['version_named'] = int(version_quoted(answer, key))
    if 'caveat_passed_on' in must or key.get('expected_codes'):
        row['caveat_passed'] = int(warns(answer))
    if 'no_false_warning' in must:
        row['no_false_warning'] = int(not warns(answer))
    if 'failure_reported' in must:
        row['failure_reported'] = int(bool(FAILURE_RE.search(answer)))
    if q['kind'] == 'batch' or 'batch_run_used' in must:
        row['batch_tool_used'] = int(any(c['name'] in BATCH_TOOLS and not c['is_error']
                                         for c in calls))
    if 'run_folder_cited' in must:
        # Models print IDs with non-breaking hyphens; the run must be one a tool reported.
        cited = answer.replace('\u2011', '-').replace('\u2010', '-')
        row['run_folder_cited'] = int(any(r in cited for r in runs_reported(t)))
    if 'coverage_stated' in must:
        row['coverage_stated'] = int(coverage_correct(answer, key))
    if 'basis_named' in must:
        row['basis_named'] = int(basis_correct(answer, key))
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
    qpath = Path(args.questions) if args.questions else HERE.parent / manifest['question_set']['path'].replace('\\', '/')
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
