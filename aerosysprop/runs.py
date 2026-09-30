# This file is part of AeroSysProp.
#
# AeroSysProp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
"""
aerosysprop.runs
================
Batch runs for large molecule sets, designed for agents.

A batch run computes properties for up to ``AEROSYSPROP_MAX_BATCH`` molecules
and writes a *run folder* instead of returning every row::

    <AEROSYSPROP_OUTPUT_DIR>/2026-09-21_143005_mcm-vp_a1b2c3/
        results.csv       one row per molecule: identifiers, values, flags
        diagnostics.csv   one row per diagnostic issue
        run.json          provenance manifest: request, method labels,
                          software version, timings, summary and the
                          SHA-256 fingerprint of each output file
        ro-crate-metadata.json
                          the same provenance in the community-standard
                          RO-Crate format (Process Run Crate profile), so
                          the folder can be archived or deposited as is

The caller gets back a compact summary and can read slices of the results on
demand. The run folder is the result; a summary of it, including an agent's
account, is a guide to the result. Summaries therefore carry a coverage
statement and caveats that must accompany any aggregate figure.

Write rules
-----------
* Runs are written only inside ``AEROSYSPROP_OUTPUT_DIR``; if it is unset,
  batch runs are disabled.
* Every run creates a new folder. Nothing is ever overwritten or deleted, and
  the session database is never touched.
* Only one run executes at a time per server process.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import secrets
import statistics
import sys
import threading
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from aerosysprop import api, batch, guidance, method_cards
from aerosysprop.api import ApiError, json_safe

OUTPUT_DIR_ENV = 'AEROSYSPROP_OUTPUT_DIR'
MAX_BATCH = int(os.environ.get('AEROSYSPROP_MAX_BATCH', '100000'))
MAX_ROWS_PER_READ = 200
SCHEMA_VERSION = 1
RESULTS_FILE, DIAGNOSTICS_FILE, MANIFEST_FILE = 'results.csv', 'diagnostics.csv', 'run.json'
CRATE_FILE = 'ro-crate-metadata.json'
PROCESS_RUN_CRATE = 'https://w3id.org/ro/wfrun/process/0.5'
_RUN_ID = re.compile(r'^\d{4}-\d{2}-\d{2}_\d{6}_[a-z0-9-]{0,40}_[0-9a-f]{6}$')

_jobs: dict = {}
_lock = threading.Lock()

__all__ = ['start_batch_run', 'batch_run_status', 'read_batch_results', 'list_batch_runs',
           'OUTPUT_DIR_ENV', 'MAX_BATCH']


# ---------------------------------------------------------------------------
# Paths and small helpers
# ---------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def _output_root() -> Path:
    root = os.environ.get(OUTPUT_DIR_ENV)
    if not root:
        raise ApiError(
            'BATCH_RUNS_DISABLED',
            'Batch runs are disabled because no output folder is configured.',
            hint=f'Start the server with {OUTPUT_DIR_ENV} set to the folder where '
                 'run folders should be written.')
    path = Path(root).expanduser().resolve()
    if not path.is_dir():
        raise ApiError('BATCH_RUNS_DISABLED', f'{OUTPUT_DIR_ENV} does not point to a folder.')
    return path


def _run_dir(run_id: str) -> Path:
    if not isinstance(run_id, str) or not _RUN_ID.match(run_id):
        raise ApiError('UNKNOWN_RUN', f'{run_id!r} is not a valid run ID.',
                       hint='Use a run_id returned by start_batch_run or list_batch_runs.')
    d = _output_root() / run_id
    if not (d / MANIFEST_FILE).is_file():
        raise ApiError('UNKNOWN_RUN', f'No run {run_id!r} in the output folder.',
                       hint='Call list_batch_runs to see available runs.')
    return d


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _write_manifest(d: Path, manifest: dict):
    # The manifest is the only file ever rewritten, and only by the run that
    # created it: written to a temporary name, then renamed into place.
    tmp = d / (MANIFEST_FILE + '.tmp')
    tmp.write_text(json.dumps(json_safe(manifest), indent=2, ensure_ascii=False), encoding='utf-8')
    os.replace(tmp, d / MANIFEST_FILE)


def _read_manifest(d: Path) -> dict:
    return json.loads((d / MANIFEST_FILE).read_text(encoding='utf-8'))


def _slug(label) -> str:
    return re.sub(r'[^a-z0-9]+', '-', (label or '').lower()).strip('-')[:40]


def _pid_alive(pid) -> Optional[bool]:
    if not isinstance(pid, int) or sys.platform == 'win32':
        return None  # cannot check safely on Windows
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


# ---------------------------------------------------------------------------
# Starting a run
# ---------------------------------------------------------------------------

def start_batch_run(input_path: Optional[str] = None, smiles=None, names=None, methods=None,
                    temperature_K=298.15, bp_sources=None,
                    label: Optional[str] = None, requested_by: Optional[str] = None,
                    smiles_column: Optional[int] = None, name_column: Optional[int] = None,
                    table: str = 'compounds', force_new: bool = False) -> dict:
    """Start a batch run in the background and return its run ID at once.

    If a completed run with an identical request already exists (same input,
    methods, temperature, boiling point sources, parameter versions and
    software version), it is returned instead of computing a duplicate; pass
    ``force_new=True`` to run again regardless.

    Give either ``input_path`` (a file in the data directory) or ``smiles``
    (with optional ``names``). Poll :func:`batch_run_status` for progress.

    Diagnostics are always computed for batch runs: the summary's split
    between clean and flagged molecules depends on them, so switching them
    off would make every closed-shell molecule look clean.
    """
    include_diagnostics = True
    root = _output_root()
    if (input_path is None) == (smiles is None):
        raise ApiError('INVALID_ARGUMENT', 'Give exactly one of input_path or smiles.')
    if input_path is not None:
        p = api._resolve_in_root(input_path)
        pairs = api._read_molecules(p, smiles_column, name_column, table)
        input_info = {'kind': 'file', 'path': input_path, 'sha256': _sha256(p),
                      'molecules_read': len(pairs)}
    else:
        smi = api._as_str_list(smiles, 'smiles')
        nm = api._as_str_list(names, 'names') if names is not None else list(smi)
        if len(nm) != len(smi):
            raise ApiError('INVALID_ARGUMENT', 'names must be the same length as smiles.')
        pairs = list(zip(nm, smi))
        input_info = {'kind': 'inline', 'molecules': len(pairs)}
    if not pairs:
        raise ApiError('EMPTY_INPUT', 'No molecules were found in the input.')
    if len(pairs) > MAX_BATCH:
        raise ApiError('TOO_MANY_MOLECULES',
                       f'{len(pairs)} molecules requested; the batch limit is {MAX_BATCH}.',
                       hint='Split the input, or raise AEROSYSPROP_MAX_BATCH.')
    method_ids, requested_sets, request_notes = api._normalise_request(methods)
    bp = api._normalise_bp_sources(bp_sources)
    T, notes = api._check_temperature(temperature_K)
    resolved = [m for m in batch.ALL_METHODS
                if m in batch.resolve_methods(method_ids or batch.ALL_METHODS, bp_sources=bp)]
    _, sets = batch.normalise_parameter_sets(resolved, requested_sets)
    provenance, prov_notes = api._provenance(resolved, sets)
    request_key = _request_key(pairs, input_info, resolved, T, bp, provenance, sets)
    if not force_new:
        existing = _find_identical(root, request_key)
        if existing is not None:
            return {'run_id': existing.name, 'status': 'complete', 'reused': True,
                    'molecules': len(pairs), 'folder': str(existing), 'methods_run': resolved,
                    'next': ('An identical completed run already exists, so no new run was made. '
                             'Call batch_run_status with this run_id for its summary.')}

    with _lock:
        busy = [rid for rid, j in _jobs.items() if j['status'] == 'running']
        if busy:
            raise ApiError('BUSY', f'Run {busy[0]} is still running.',
                           hint='Wait for it to finish (batch_run_status), then start the next.')
        stamp = datetime.now().strftime('%Y-%m-%d_%H%M%S')
        run_id = f'{stamp}_{_slug(label)}_{secrets.token_hex(3)}'
        d = root / run_id
        d.mkdir(exist_ok=False)  # a new folder every time; never reuse one
        manifest = {
            'schema_version': SCHEMA_VERSION,
            'run_id': run_id,
            'label': label,
            'status': 'running',
            'created': _now(),
            'finished': None,
            'requested_by': requested_by,
            'pid': os.getpid(),
            'request': {
                'input': input_info,
                'methods_requested': method_ids or 'all',
                'methods_run': resolved,
                'parameter_sets': {k: list(v) for k, v in sets.items()},
                'temperature_K': T,
                'boiling_point_sources': bp,
                'include_diagnostics': bool(include_diagnostics),
                'request_key': request_key,
            },
            'provenance': provenance,
            'warnings': request_notes + notes + prov_notes,
            'outputs': {},
            'summary': None,
        }
        _write_manifest(d, manifest)
        _jobs[run_id] = {'status': 'running', 'done': 0, 'total': len(pairs),
                         'started': time.time()}
    threading.Thread(target=_execute, name=f'aerosysprop-run-{run_id}', daemon=True,
                     args=(d, pairs, method_ids, T, bp, include_diagnostics, manifest, sets)).start()
    return {'run_id': run_id, 'status': 'running', 'reused': False, 'molecules': len(pairs),
            'folder': str(d), 'methods_run': resolved,
            'next': 'Call batch_run_status with this run_id until status is complete.'}


def _request_key(pairs, input_info, resolved, T, bp, provenance, sets=None) -> str:
    """Fingerprint of everything that determines a run's results."""
    h = hashlib.sha256()
    for n, smi in pairs:
        h.update(f'{n}\t{smi}\n'.encode())
    key = {'molecules': h.hexdigest(), 'methods': resolved, 'temperature_K': T,
           'parameter_sets': {k: list(v) for k, v in (sets or {}).items()},
           'bp_sources': bp, 'software': provenance['software'],
           'method_labels': {m: p.get('label') for m, p in provenance['methods'].items()},
           'parameters_match_cards': all(p.get('parameters_match_card', False)
                                         for p in provenance['methods'].values()
                                         if p.get('documented'))}
    return hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()


def _find_identical(root: Path, request_key: str) -> Optional[Path]:
    """A completed, unaltered run with the same request key, if one exists."""
    for d in sorted((p for p in root.iterdir() if p.is_dir() and _RUN_ID.match(p.name)),
                    key=lambda p: p.name, reverse=True):
        try:
            m = _read_manifest(d)
        except (OSError, ValueError):
            continue
        if m.get('status') != 'complete' or m['request'].get('request_key') != request_key:
            continue
        outputs = m.get('outputs', {})
        if all((d / f).is_file() and _sha256(d / f) == o.get('sha256')
               for f, o in outputs.items()):
            return d
    return None


def _folder_size(d: Path) -> int:
    return sum(f.stat().st_size for f in d.rglob('*') if f.is_file())


def _execute(d: Path, pairs, method_ids, T, bp, include_diagnostics, manifest, sets=None):
    run_id = manifest['run_id']
    job = _jobs[run_id]

    def progress(done, total):
        job['done'], job['total'] = done, total

    try:
        names = [n for n, _ in pairs]
        smiles = [s for _, s in pairs]
        resolved, records, columns = api._compute(smiles, names, method_ids, T, bp,
                                                  include_diagnostics, progress,
                                                  parameter_sets=sets)
        value_cols = [c for c in batch.resolved_output_columns(set(resolved), bp, sets)
                      if c in columns]
        value_cols += [c for c in columns if c not in value_cols]
        n_diag = _write_outputs(d, records, value_cols)
        manifest['outputs'] = {
            RESULTS_FILE: {'rows': len(records), 'sha256': _sha256(d / RESULTS_FILE),
                           'columns': {c: columns.get(c, '') for c in value_cols}},
            DIAGNOSTICS_FILE: {'rows': n_diag, 'sha256': _sha256(d / DIAGNOSTICS_FILE)},
        }
        manifest['summary'] = _summarise(records, value_cols, columns)
        manifest['diagnostic_codes'] = DIAGNOSTIC_CODES
        manifest['status'] = 'complete'
    except Exception as e:  # recorded in the manifest; nothing is left half-written
        manifest['status'] = 'failed'
        manifest['error'] = f'{type(e).__name__}: {e}'[:2000]
    finally:
        manifest['finished'] = _now()
        manifest['elapsed_s'] = round(time.time() - job['started'], 1)
        try:
            _write_crate(d, manifest)
        except Exception as e:  # the crate is a convenience; never lose the run over it
            manifest.setdefault('warnings', []).append(
                {'severity': 'warning', 'code': 'CRATE_NOT_WRITTEN', 'message': str(e)[:500]})
        _write_manifest(d, manifest)
        job['status'] = manifest['status']


def _write_outputs(d: Path, records, value_cols) -> int:
    head = ['name', 'input_smiles', 'canonical_smiles', 'formula', 'radical_electrons'] + value_cols + \
           ['flagged_methods', 'not_computed', 'error']
    tmp = d / (RESULTS_FILE + '.partial')
    with open(tmp, 'x', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(head)
        for r in records:
            vals = r.get('values', {})
            w.writerow([r['name'], r['input_smiles'], r.get('canonical_smiles', ''),
                        r.get('formula', ''), r.get('radical_electrons', '')] +
                       ['' if vals.get(c) is None else repr(vals[c]) for c in value_cols] +
                       [';'.join(sorted(r.get('diagnostics', {}))),
                        ';'.join(r.get('not_computed', [])),
                        r['error']['code'] if 'error' in r else ''])
    os.rename(tmp, d / RESULTS_FILE)
    n = 0
    with open(d / DIAGNOSTICS_FILE, 'x', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        # Compact form: SMILES are in results.csv, and the full wording of each
        # code is given once in run.json (diagnostic_codes).
        w.writerow(['name', 'method', 'severity', 'code', 'detail'])
        for r in records:
            for method, issues in sorted(r.get('diagnostics', {}).items()):
                for i in issues:
                    w.writerow([r['name'], method, i['severity'], i['code'], _detail(i)])
                    n += 1
    return n


def _detail(issue) -> str:
    unmatched = (issue.get('details') or {}).get('unmatched')
    if unmatched:
        return f'{len(unmatched)} unmatched heavy atom(s): ' + ','.join(str(a[1]) for a in unmatched)
    return issue['message'][:160]


DIAGNOSTIC_CODES = {
    'UNMATCHED_ATOMS': 'The method\'s functional groups do not cover these heavy atoms; the '
                       'method behaves as though they (and their groups) were absent.',
    'INVALID_SMILES': 'RDKit could not parse the SMILES.',
}


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def _column_groups(value_cols):
    """Group columns that estimate the same quantity in the same unit, so the
    spread between methods (or boiling point sources) can be reported."""
    key_to_method = {k: m for m in batch.PROPERTY_METHODS for k in m.output_keys}
    groups = {}
    for c in value_cols:
        base = c.split('__', 1)[0]  # strip __bp<source> and __set<n> suffixes
        m = key_to_method.get(base)
        if m is None:
            continue
        quantity = m.group if len(m.output_keys) == 1 else base.rsplit('_', 1)[0].split('_')[0]
        groups.setdefault(f'{quantity} ({api._unit_for(c)})', []).append(c)
    return {g: cols for g, cols in groups.items() if len(cols) > 1}


def _spreads(rows, cols):
    out = []
    for r in rows:
        vals = [r[c] for c in cols if r.get(c) is not None]
        if len(vals) >= 2:
            out.append((max(vals) - min(vals), r))
    return out


def _stats(xs):
    return {'n': len(xs), 'min': min(xs) if xs else None,
            'median': statistics.median(xs) if xs else None, 'max': max(xs) if xs else None}


def _is_clean(r) -> bool:
    """Closed-shell molecule with no diagnostic warning from any method."""
    return not r.get('radical_electrons') and not r.get('diagnostics')


CAVEATS = [
    'Values are group-contribution estimates, not measurements.',
    'Agreement between methods does not establish accuracy: methods can share '
    'the same blind spot (for example, a functional group none of them represents).',
    'Aggregate figures hide molecule-level caveats. Statistics are given for '
    '"clean" molecules (closed-shell, no method flagged) and for all molecules; '
    'quote the clean figures, and the coverage, whenever reporting aggregates.',
    'Radical species are outside the design range of these methods; their values '
    'should not be reported as pure-component properties.',
    'The run folder (results.csv, diagnostics.csv, run.json) is the result. A '
    'summary is a guide to it; check specific molecules against results.csv.',
]

REPORTING_REQUIREMENTS = [
    'Give the user the run folder path.',
    'Quote provenance labels (method@version) for any value reported.',
    'When reporting an aggregate, state its basis (clean or all molecules) and '
    'the coverage figures.',
    'For a specific molecule, read its row (read_batch_results with names or '
    'smiles) rather than inferring from the summary, and report any flags on it.',
]


def _summarise(records, value_cols, units) -> dict:
    valid = [r for r in records if 'error' not in r]
    clean = [r for r in valid if _is_clean(r)]

    def rows_of(rs):
        return [dict(name=r['name'], canonical_smiles=r['canonical_smiles'],
                     flagged=sorted(r.get('diagnostics', {})),
                     radical=bool(r.get('radical_electrons')), **r['values']) for r in rs]
    rows_all, rows_clean = rows_of(valid), rows_of(clean)

    stats = {}
    for c in value_cols:
        stats[c] = {'unit': units.get(c, ''),
                    'clean': _stats([r[c] for r in rows_clean if r.get(c) is not None]),
                    'all': _stats([r[c] for r in rows_all if r.get(c) is not None])}

    flags = Counter(m for r in valid for m in r.get('diagnostics', {}))
    n = len(valid) or 1
    radicals = sum(1 for r in valid if r.get('radical_electrons'))
    flagged = sum(1 for r in valid if r.get('diagnostics'))
    coverage = {
        'molecules': len(records),
        'computed': len(valid),
        'invalid_smiles': len(records) - len(valid),
        'clean': len(clean),
        'clean_fraction': round(len(clean) / n, 3),
        'flagged': flagged,
        'flagged_fraction': round(flagged / n, 3),
        'radical_species': radicals,
        'radical_fraction': round(radicals / n, 3),
        'flagged_fraction_by_method': {m: round(k / n, 3) for m, k in flags.most_common()},
        'statement': (f'{len(clean)} of {len(valid)} molecules ({100 * len(clean) / n:.0f}%) are '
                      f'clean. {flagged} ({100 * flagged / n:.0f}%) have at least one method '
                      f'flagged, and {radicals} ({100 * radicals / n:.0f}%) are radicals.'),
    }

    spreads = {}
    for group, cols in _column_groups(value_cols).items():
        sp_clean, sp_all = _spreads(rows_clean, cols), _spreads(rows_all, cols)
        if not sp_all:
            continue
        vals = sorted(x for x, _ in sp_clean)
        top = sorted(sp_clean, key=lambda t: -t[0])[:10]
        spreads[group] = {
            'columns': cols,
            'basis': 'clean molecules',
            'molecules_compared': len(sp_clean),
            'median_spread': statistics.median(vals) if vals else None,
            'p90_spread': vals[int(0.9 * (len(vals) - 1))] if vals else None,
            'median_spread_all_molecules': statistics.median(sorted(x for x, _ in sp_all)),
            'largest': [{'name': r['name'], 'smiles': r['canonical_smiles'], 'spread': x,
                         'values': {c: r.get(c) for c in cols}} for x, r in top],
        }
    caveats = list(CAVEATS)
    if coverage['clean_fraction'] < 0.5:
        caveats.insert(0, f'Fewer than half the molecules are clean '
                          f'({coverage["clean_fraction"]:.0%}); aggregates over all molecules '
                          'are dominated by values the methods are not designed for.')
    return {
        'coverage': coverage,
        'caveats': caveats,
        'reporting_requirements': REPORTING_REQUIREMENTS,
        'columns': stats,
        'method_spread': spreads,
    }


def _write_crate(d: Path, m: dict):
    """Write ro-crate-metadata.json (RO-Crate 1.1, Process Run Crate 0.5)."""
    sw = m['provenance']['software']
    req = m['request']
    action_id = '#run-' + m['run_id']
    params = [
        {'@id': '#param-methods', '@type': 'PropertyValue', 'name': 'methods',
         'value': ', '.join(req['methods_run'])},
        {'@id': '#param-temperature', '@type': 'PropertyValue', 'name': 'temperature',
         'value': req['temperature_K'], 'unitText': 'K'},
        {'@id': '#param-bp-sources', '@type': 'PropertyValue', 'name': 'boiling point sources',
         'value': ', '.join(req['boiling_point_sources'])},
    ]
    if req.get('parameter_sets'):
        params.append({'@id': '#param-parameter-sets', '@type': 'PropertyValue',
                       'name': 'parameter sets',
                       'value': ', '.join(f'{k}@{v}' for k, vs in req['parameter_sets'].items()
                                          for v in vs)})
    inp = req['input']
    params.append({'@id': '#param-input', '@type': 'PropertyValue', 'name': 'input',
                   'value': (f'{inp["path"]} (sha256 {inp["sha256"]})' if inp['kind'] == 'file'
                             else f'{inp["molecules"]} SMILES given inline')})
    methods = [{'@id': f'#method-{mid}', '@type': 'PropertyValue', 'name': mid,
                'value': pv.get('label'), 'description': f'status: {pv.get("status")}; '
                f'review: {pv.get("review_state")}'}
               for mid, pv in m['provenance']['methods'].items()]
    files = [f for f in (RESULTS_FILE, DIAGNOSTICS_FILE, MANIFEST_FILE) if (d / f).exists()
             or f == MANIFEST_FILE]
    file_entities = [{'@id': RESULTS_FILE, '@type': 'File', 'name': 'Predicted properties',
                      'encodingFormat': 'text/csv'},
                     {'@id': DIAGNOSTICS_FILE, '@type': 'File', 'name': 'Diagnostic warnings',
                      'encodingFormat': 'text/csv'},
                     {'@id': MANIFEST_FILE, '@type': 'File', 'name': 'AeroSysProp run manifest',
                      'encodingFormat': 'application/json',
                      'description': 'Full provenance, summary and SHA-256 fingerprints of the outputs.'}]
    file_entities = [f for f in file_entities if f['@id'] in files]
    graph = [
        {'@id': CRATE_FILE, '@type': 'CreativeWork',
         'conformsTo': {'@id': 'https://w3id.org/ro/crate/1.1'}, 'about': {'@id': './'}},
        {'@id': './', '@type': 'Dataset', 'conformsTo': {'@id': PROCESS_RUN_CRATE},
         'name': f'AeroSysProp batch run {m["run_id"]}',
         'description': m.get('label') or 'AeroSysProp batch run',
         'datePublished': m['finished'] or m['created'],
         'hasPart': [{'@id': f['@id']} for f in file_entities],
         'mentions': {'@id': action_id}},
        {'@id': PROCESS_RUN_CRATE, '@type': 'CreativeWork', 'name': 'Process Run Crate',
         'version': '0.5'},
        {'@id': action_id, '@type': 'CreateAction', 'name': f'AeroSysProp batch run {m["run_id"]}',
         'startTime': m['created'], 'endTime': m['finished'],
         'actionStatus': {'@id': 'http://schema.org/CompletedActionStatus'
                          if m['status'] == 'complete' else 'http://schema.org/FailedActionStatus'},
         'instrument': {'@id': '#aerosysprop'},
         'object': [{'@id': p['@id']} for p in params + methods],
         'result': [{'@id': f['@id']} for f in file_entities if f['@id'] != MANIFEST_FILE]},
        {'@id': '#aerosysprop', '@type': 'SoftwareApplication', 'name': 'AeroSysProp',
         'softwareVersion': sw.get('aerosysprop_version'),
         'description': f'git commit {sw.get("git_commit")}' if sw.get('git_commit') else None},
    ] + params + methods + file_entities
    if m['status'] != 'complete' and m.get('error'):
        graph[3]['error'] = m['error']
    crate = {'@context': 'https://w3id.org/ro/crate/1.1/context',
             '@graph': [{k: v for k, v in e.items() if v is not None} for e in graph]}
    (d / CRATE_FILE).write_text(json.dumps(json_safe(crate), indent=2), encoding='utf-8')


# ---------------------------------------------------------------------------
# Status, reading and listing
# ---------------------------------------------------------------------------

def _served_summary(summary: dict) -> dict:
    """The summary as served at the current guidance level. run.json always
    keeps the full summary; at the bare level the coverage statement, caveats
    and reporting requirements are left out and the statistics are unchanged."""
    if not guidance.is_bare():
        return summary
    served = {k: v for k, v in summary.items() if k not in ('caveats', 'reporting_requirements')}
    served['coverage'] = {k: v for k, v in summary['coverage'].items() if k != 'statement'}
    return served


def batch_run_status(run_id: str) -> dict:
    """Status of a run; once complete, includes its summary and provenance."""
    d = _run_dir(run_id)
    m = _read_manifest(d)
    status = m['status']
    out = {'run_id': run_id, 'status': status, 'folder': str(d), 'created': m['created']}
    if status == 'running':
        job = _jobs.get(run_id)
        if job is not None:
            out['progress'] = {'done': job['done'], 'total': job['total'],
                               'elapsed_s': round(time.time() - job['started'], 1)}
        elif _pid_alive(m.get('pid')) is False:
            out['status'] = 'interrupted'
            out['note'] = ('The server process running this job has stopped; the run did not '
                           'finish. Start it again with start_batch_run.')
        else:
            out['note'] = 'Running in another server process; progress is not visible here.'
        return out
    out.update({'finished': m.get('finished'), 'elapsed_s': m.get('elapsed_s')})
    if status == 'failed':
        out['error'] = m.get('error')
        return out
    out['summary'] = _served_summary(m['summary'])
    out['provenance'] = {'software': m['provenance']['software'],
                         'methods': {k: v.get('label') for k, v in m['provenance']['methods'].items()}}
    out['outputs'] = {k: {'rows': v['rows']} for k, v in m['outputs'].items()}
    out['warnings'] = m.get('warnings', [])
    return out


def _load_rows(d: Path, m: dict):
    cols = list(m['outputs'][RESULTS_FILE]['columns'])
    rows = []
    with open(d / RESULTS_FILE, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            for c in cols:
                r[c] = float(r[c]) if r[c] != '' else None
            rows.append(r)
    return cols, rows


def read_batch_results(run_id: str, view: str = 'rows', columns=None,
                       method: Optional[str] = None, name_contains: Optional[str] = None,
                       group: Optional[str] = None, limit: int = 50, offset: int = 0,
                       names=None, smiles=None) -> dict:
    """Read a slice of a completed run. Read-only.

    ``view``: ``'rows'`` (all rows, paged), ``'flagged'`` (rows with a
    diagnostic warning, optionally for one ``method``) or ``'largest_spread'``
    (rows where methods disagree most within ``group``, as named in the
    summary's method_spread). ``names`` and ``smiles`` select exact
    molecules (SMILES are compared after canonicalisation); any not found are
    listed in ``not_found``.
    """
    d = _run_dir(run_id)
    m = _read_manifest(d)
    if m['status'] != 'complete':
        raise ApiError('RUN_NOT_COMPLETE', f'Run {run_id} is {m["status"]}.',
                       hint='Call batch_run_status and wait until status is complete.')
    matches = _sha256(d / RESULTS_FILE) == m['outputs'][RESULTS_FILE]['sha256']
    cols, rows = _load_rows(d, m)
    not_found = []
    if names is not None or smiles is not None:
        wanted_names = api._as_str_list(names, 'names') if names is not None else []
        wanted_smiles = {}
        for smi in (api._as_str_list(smiles, 'smiles') if smiles is not None else []):
            mol = api._parse_smiles(smi)
            wanted_smiles[api.Chem.MolToSmiles(mol) if mol is not None else smi] = smi
        by_name = {r['name'] for r in rows}
        by_smiles = {r['canonical_smiles'] for r in rows}
        not_found = [n for n in wanted_names if n not in by_name] + \
                    [orig for can, orig in wanted_smiles.items() if can not in by_smiles]
        rows = [r for r in rows if r['name'] in wanted_names
                or r['canonical_smiles'] in wanted_smiles]
    if name_contains:
        rows = [r for r in rows if name_contains.lower() in r['name'].lower()]
    if view == 'rows':
        pass
    elif view == 'flagged':
        rows = [r for r in rows if r['flagged_methods'] and
                (method is None or method in r['flagged_methods'].split(';'))]
    elif view == 'largest_spread':
        groups = _column_groups(cols)
        if group not in groups:
            raise ApiError('INVALID_ARGUMENT', f'Unknown group {group!r}.',
                           details={'groups': list(groups)})
        rows = [r for s, r in sorted(_spreads(rows, groups[group]), key=lambda t: -t[0])]
        for r in rows:
            vals = [r[c] for c in groups[group] if r.get(c) is not None]
            r['spread'] = max(vals) - min(vals)
    else:
        raise ApiError('INVALID_ARGUMENT', f'Unknown view {view!r}.',
                       details={'views': ['rows', 'flagged', 'largest_spread']})
    keep = None
    if columns:
        wanted = api._as_str_list(columns, 'columns')
        unknown = [c for c in wanted if c not in cols]
        if unknown:
            raise ApiError('INVALID_ARGUMENT', f'Unknown columns {unknown}.',
                           details={'columns': cols})
        keep = ['name', 'canonical_smiles', 'radical_electrons'] + wanted + \
               ['flagged_methods', 'not_computed', 'spread']
    limit = max(1, min(int(limit), MAX_ROWS_PER_READ))
    offset = max(0, int(offset))
    page = rows[offset:offset + limit]
    if keep:
        page = [{k: r[k] for k in keep if k in r} for r in page]
    return {'run_id': run_id, 'view': view, 'total_matching': len(rows), 'offset': offset,
            'not_found': not_found,
            'returned': len(page), 'rows': json_safe(page),
            'units': m['outputs'][RESULTS_FILE]['columns'],
            'results_match_manifest': matches,
            'note': None if matches else
            'results.csv has been changed since the run wrote it; values may not be original.'}


def list_batch_runs(limit: int = 20) -> dict:
    """List runs in the output folder, newest first."""
    root = _output_root()
    runs = []
    for d in sorted((p for p in root.iterdir() if p.is_dir() and _RUN_ID.match(p.name)),
                    key=lambda p: p.name, reverse=True)[:max(1, min(int(limit), 200))]:
        try:
            m = _read_manifest(d)
        except (OSError, ValueError):
            continue
        runs.append({'run_id': m['run_id'], 'label': m.get('label'), 'status': m['status'],
                     'created': m['created'],
                     'molecules': ((m.get('summary') or {}).get('coverage') or {}).get('molecules')
                     or m['request']['input'].get('molecules')
                     or m['request']['input'].get('molecules_read'),
                     'methods_run': m['request']['methods_run'],
                     'temperature_K': m['request']['temperature_K'],
                     'size_mb': round(_folder_size(d) / 1e6, 2)})
    all_runs = [p for p in root.iterdir() if p.is_dir() and _RUN_ID.match(p.name)]
    total = sum(_folder_size(p) for p in all_runs)
    return {'output_folder': str(root), 'runs': runs, 'total_runs': len(all_runs),
            'total_size_mb': round(total / 1e6, 1),
            'note': 'Runs are never deleted automatically. To remove one, delete its folder by hand.'}
