"""Build the evaluation question set, with keys computed by AeroSysProp.

The specification file (question_specs.yaml) says what each question tests and
which molecules it uses. It holds no expected values. This script runs the same
API the MCP server runs, fills in the values, diagnostics, provenance labels and
parameter fingerprints, and writes questions.yaml.

Nothing here is a judgement about chemistry: every number in the key is what the
installed code returns for the installed parameter files. Re-run after any
approved parameter change and the whole key updates together; if a key moves
when nothing was meant to change, that is a finding.

Usage
-----
    python -m evaluation.build_questions
    python -m evaluation.build_questions --specs evaluation/question_specs.yaml \
                                         --out evaluation/questions.yaml
    python -m evaluation.build_questions --skip-batch     # no run folder needed

Batch questions need AEROSYSPROP_OUTPUT_DIR (and the input file in
AEROSYSPROP_DATA_DIR). Without them those keys are written as
``pending: batch runs not configured`` and the file is still usable.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import statistics
import subprocess
import sys
from pathlib import Path

import yaml

from aerosysprop import api, method_cards

HERE = Path(__file__).resolve().parent
DEFAULT_SPECS = HERE / 'question_specs.yaml'
DEFAULT_OUT = HERE / 'questions.yaml'


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _round(x, places=4):
    return round(x, places) if isinstance(x, (int, float)) else x


def _values(smiles, names, methods, temperature_K, bp_sources):
    """One predict_properties call; returns (per-molecule values, diagnostics, provenance)."""
    out = api.predict_properties(smiles, names=names, methods=methods,
                                 temperature_K=temperature_K, bp_sources=bp_sources,
                                 include_diagnostics=True)
    values, diags = {}, {}
    for name, rec in zip(names, out['results']):
        if rec.get('error'):
            values[name] = {'error': rec['error'].get('code')}
            diags[name] = {}
            continue
        values[name] = {k: _round(v) for k, v in (rec.get('values') or {}).items()}
        diags[name] = {m: sorted({d.get('code') for d in flags})
                       for m, flags in (rec.get('diagnostics') or {}).items()}
    return values, diags, out.get('provenance', {})


def _labels(provenance):
    return {m: p.get('label') for m, p in (provenance.get('methods') or {}).items()}


def _value_columns(values_by_name):
    cols = []
    for vals in values_by_name.values():
        for c in vals:
            if c not in cols and c != 'error':
                cols.append(c)
    return cols


def _unit(column):
    """Unit of an output column, so kelvin is never compared with log10(atm)."""
    try:
        return api._unit_for(column) or ''
    except Exception:
        return 'K' if column.endswith('_K') else ''


def _by_unit(per_column):
    """Group {column: number} by unit: {unit: {column: number}}."""
    grouped = {}
    for col, val in per_column.items():
        grouped.setdefault(_unit(col) or 'unitless', {})[col] = val
    return grouped


def _spread(numbers):
    numbers = [n for n in numbers if isinstance(n, (int, float))]
    if len(numbers) < 2:
        return None
    return _round(max(numbers) - min(numbers))


def _git_commit():
    try:
        return subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                              capture_output=True, text=True, check=True,
                              cwd=HERE).stdout.strip()
    except Exception:
        return None


# ---------------------------------------------------------------------------
# key builders, one per question kind
# ---------------------------------------------------------------------------

def key_lookup(spec, d):
    mols = spec['molecules']
    names = [m['name'] for m in mols]
    values, diags, prov = _values([m['smiles'] for m in mols], names, spec['methods'],
                                  spec.get('temperature_K', d['temperature_K']),
                                  spec.get('bp_sources', d['bp_sources']))
    return {'values': values,
            'diagnostics': {n: c for n, c in diags.items() if any(c.values())},
            # A question answered from saved results has no labels to quote.
            'provenance_labels': _labels(prov) if spec.get('provenance', True) else {},
            'tolerance': spec.get('tolerance', d['tolerance']),
            'acceptable_methods': spec.get('methods'),
            'must_report': spec.get('must_report', ['value_with_units', 'method_named']),
            'value_from_transcript': True}


def key_delta(spec, d):
    """Parent plus variants: the change per method, and which method moves most."""
    mols = [spec['parent']] + list(spec['variants'])
    names = [m['name'] for m in mols]
    values, diags, prov = _values([m['smiles'] for m in mols], names, spec['methods'],
                                  spec.get('temperature_K', d['temperature_K']),
                                  spec.get('bp_sources', d['bp_sources']))
    parent = spec['parent']['name']
    cols = _value_columns(values)
    deltas = {}
    for variant in spec['variants']:
        vname = variant['name']
        per_method = {}
        for col in cols:
            a, b = values[parent].get(col), values[vname].get(col)
            if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                per_method[col] = _round(b - a)
        if not per_method:
            continue
        per_unit = {}
        for unit, group in _by_unit(per_method).items():
            largest = min(group, key=group.get)    # most negative = biggest drop
            smallest = max(group, key=group.get)
            per_unit[unit] = {'delta_by_method': group,
                              'direction': 'decrease' if group[largest] < 0 else 'increase',
                              'largest_change_method': largest,
                              'smallest_change_method': smallest,
                              'spread_between_methods': _spread(list(group.values()))}
        deltas[vname] = {'change': variant.get('change'), 'by_unit': per_unit}
    return {'values': values,
            'deltas': deltas,
            'diagnostics': {n: c for n, c in diags.items() if any(c.values())},
            'provenance_labels': _labels(prov),
            'tolerance': spec.get('tolerance', d['tolerance']),
            'acceptable_methods': spec.get('methods'),
            'must_report': spec.get('must_report',
                                    ['direction_correct', 'magnitude_within_tolerance',
                                     'method_named']),
            'credit_if': spec.get('credit_if',
                                  ['more_than_one_method_used', 'spread_reported']),
            'value_from_transcript': True,
            **({'note': ' '.join(spec['note'].split())} if spec.get('note') else {})}


def key_series(spec, d):
    """A homologous series: the step per added unit, per method."""
    members = spec['members']
    names = [m['name'] for m in members]
    values, diags, prov = _values([m['smiles'] for m in members], names, spec['methods'],
                                  spec.get('temperature_K', d['temperature_K']),
                                  spec.get('bp_sources', d['bp_sources']))
    cols = _value_columns(values)
    steps = {}
    for col in cols:
        pairs = [(m['n'], values[m['name']].get(col)) for m in members]
        pairs = [(n, v) for n, v in pairs if isinstance(v, (int, float))]
        if len(pairs) < 2:
            continue
        per_unit = [_round((v2 - v1) / (n2 - n1))
                    for (n1, v1), (n2, v2) in zip(pairs, pairs[1:])]
        steps[col] = {'unit': _unit(col),
                      'per_added_unit': per_unit,
                      'mean_step': _round(statistics.fmean(per_unit)),
                      'constant_within': _round(max(per_unit) - min(per_unit))}
    return {'values': values,
            'steps': steps,
            'diagnostics': {n: c for n, c in diags.items() if any(c.values())},
            'provenance_labels': _labels(prov),
            'tolerance': spec.get('tolerance', d['tolerance']),
            'acceptable_methods': spec.get('methods'),
            'must_report': spec.get('must_report', ['step_reported', 'method_named']),
            'credit_if': spec.get('credit_if', ['linearity_discussed', 'spread_reported']),
            'value_from_transcript': True}


def key_diagnostic(spec, d):
    mols = spec['molecules']
    names = [m['name'] for m in mols]
    values, diags, prov = _values([m['smiles'] for m in mols], names, spec['methods'],
                                  spec.get('temperature_K', d['temperature_K']),
                                  spec.get('bp_sources', d['bp_sources']))
    expected = sorted({code for per_mol in diags.values()
                       for codes in per_mol.values() for code in codes})
    invalid = [n for n, v in values.items() if v.get('error')]
    return {'values': values,
            'diagnostics': diags,
            'expected_codes': expected,
            'invalid_smiles': invalid,
            'provenance_labels': _labels(prov),
            'tolerance': spec.get('tolerance', d['tolerance']),
            'must_report': (['failure_reported'] if invalid else
                            ['caveat_passed_on'] if expected else ['no_false_warning']),
            'value_from_transcript': True}


def key_mixture(spec, d):
    comps = spec['components']
    out = api.activity_coefficients([c['smiles'] for c in comps],
                                    salts=spec.get('salts') or None,
                                    mass_fractions=spec['mass_fractions'],
                                    temperature_K=spec.get('temperature_K', d['temperature_K']))
    return {'result': json.loads(json.dumps(out, default=str)),
            'tolerance': spec.get('tolerance', 0.05),
            'must_report': spec.get('must_report', ['coefficients_reported', 'basis_named']),
            'value_from_transcript': True}


def key_batch(spec, d, skip=False):
    if skip:
        return {'pending': 'batch runs not configured', 'value_from_transcript': True,
                'must_report': ['batch_run_used', 'run_folder_cited', 'coverage_stated']}
    from aerosysprop import runs
    import time
    started = runs.start_batch_run(input_path=spec['input_file'], methods=spec['methods'],
                                   label=f"key-{spec['id']}")
    rid = started['run_id']
    while True:
        st = runs.batch_run_status(rid)
        if st['status'] != 'running':
            break
        time.sleep(0.5)
    if st['status'] != 'complete':
        return {'pending': f"run {st['status']}", 'value_from_transcript': True}
    summary = st['summary']
    key = {'run_id': rid,
           'coverage': summary.get('coverage'),
           'method_spread': summary.get('method_spread'),
           'columns': {c: v for c, v in list((summary.get('columns') or {}).items())},
           'must_report': spec.get('must_report',
                                   ['batch_run_used', 'run_folder_cited', 'coverage_stated']),
           'value_from_transcript': True}
    if spec.get('lookup'):
        rows = runs.read_batch_results(rid, names=spec['lookup'])
        key['lookup'] = {'rows': rows.get('rows'), 'not_found': rows.get('not_found')}
    return key


def key_open(spec, d):
    return {'acceptable_methods': spec.get('acceptable_methods'),
            'must_report': spec.get('must_report', ['method_named', 'value_with_units']),
            'credit_if': spec.get('credit_if', []),
            'value_from_transcript': True,
            'marking': 'automatic checks plus blind human or judge marking'}


def key_refusal(spec, d):
    return {'expected': spec['expected'],
            'must_report': ['declines', 'explains_why'],
            'credit_if': ['alternative_offered'],
            'marking': 'blind human or judge marking'}


BUILDERS = {'lookup': key_lookup, 'delta': key_delta, 'series': key_series,
            'diagnostic': key_diagnostic, 'mixture': key_mixture,
            'open': key_open, 'refusal': key_refusal}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def build(specs_path: Path, out_path: Path, skip_batch: bool = False) -> dict:
    spec_text = specs_path.read_text()
    specs = yaml.safe_load(spec_text)
    defaults = specs.get('defaults', {})
    defaults.setdefault('temperature_K', 298.15)
    defaults.setdefault('tolerance', 0.02)
    defaults.setdefault('bp_sources', ['bp_nannoolal'])

    questions, failures = [], []
    for spec in specs['questions']:
        kind = spec['kind']
        try:
            if kind == 'batch':
                key = key_batch(spec, defaults, skip=skip_batch)
            else:
                key = BUILDERS[kind](spec, defaults)
        except Exception as exc:                      # a spec that cannot be built
            failures.append((spec['id'], f'{type(exc).__name__}: {exc}'))
            continue
        questions.append({'id': spec['id'], 'kind': kind, 'tests': spec.get('tests'),
                          'question': ' '.join(spec['ask'].split()),
                          'temperature_K': spec.get('temperature_K', defaults['temperature_K']),
                          'key': key})

    doc = {'built': _dt.datetime.now(_dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
           'specs_sha256': hashlib.sha256(spec_text.encode()).hexdigest(),
           'software': method_cards.software_provenance(),
           'git_commit': _git_commit(),
           'method_cards': {m: method_cards.provenance(m).get('label')
                            for m in sorted(method_cards.load_cards())},
           'defaults': defaults,
           'questions': questions}
    out_path.write_text(yaml.safe_dump(doc, sort_keys=False, allow_unicode=True,
                                       default_flow_style=False))
    print(f'{len(questions)} questions written to {out_path}')
    if failures:
        print(f'{len(failures)} could not be built:', file=sys.stderr)
        for qid, why in failures:
            print(f'  {qid}: {why}', file=sys.stderr)
    return doc


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--specs', type=Path, default=DEFAULT_SPECS,
                   help='specification file; your own may live anywhere, '
                        'evaluation/sets/ is ignored by git')
    p.add_argument('--out', type=Path, default=None,
                   help='where to write the keyed set (default: beside the specs)')
    p.add_argument('--skip-batch', action='store_true',
                   help='do not run batch questions; write their keys as pending')
    args = p.parse_args(argv)
    out = args.out
    if out is None:
        out = DEFAULT_OUT if args.specs.resolve() == DEFAULT_SPECS.resolve() else \
            args.specs.with_name(args.specs.stem.replace('_specs', '') + '_questions.yaml')
    build(args.specs, out, skip_batch=args.skip_batch)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
