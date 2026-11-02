"""Ask each model each question under each guidance level, and record everything.

For every model, guidance level, question and repeat the harness starts a fresh
AeroSysProp MCP server as a subprocess, gives the model the question and the
server's tool list, passes the model's tool calls to the server until it stops
asking (or hits the call cap), and writes a transcript. Models are served by
Ollama on this machine, or with --backend anthropic by the Anthropic API (needs
the anthropic package and ANTHROPIC_API_KEY). Both backends get the same system
prompt, the same tools, the same call cap and a fresh server per conversation,
and write transcripts in the same shape, so mark.py and judge.py read either.

    python -m evaluation.run_eval --models qwen3:8b,llama3.1:8b --out evaluation/results/pilot
    python -m evaluation.run_eval --models qwen3:4b --questions-only lookup-01,diag-01 \\
        --repeats 1 --out evaluation/results/smoke
    python -m evaluation.run_eval --backend anthropic --models claude-sonnet-5-5 \\
        --questions-only lookup-01,diag-01 --repeats 1 --out evaluation/results/claude-smoke

Every conversation starts from the same state: an output folder holding one
completed run of the batch questions' request over mcm_subset.txt, so "that run"
means the same thing whatever order the questions come in. Each conversation
gets its own copy under an opaque name (tools report the folder path, so it must
not name the model or condition); copies the model did not add to are removed.

The run can be stopped and started again: finished transcripts are skipped.
Marking is separate (evaluation/mark.py).
"""
from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import json
import logging
import os
import platform
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import yaml

HERE = Path(__file__).parent
REPO = HERE.parent

SYSTEM_PROMPT = ('You answer questions about molecular properties using the AeroSysProp tools. '
                 'Do not answer from memory. If a tool reports a problem with a molecule, '
                 'say so in your answer.')
# 'model' (the default) keeps each model's own recommended sampling (its Ollama
# defaults, recorded in the manifest) with a fixed seed: Qwen3, for one, warns
# that greedy decoding makes its thinking models repeat themselves, and in the
# Spark smoke test qwen3:4b looped for 12,000 tokens at temperature 0.
# 'greedy' pins temperature 0 for every model.
SAMPLING_MODES = {'greedy': {'temperature': 0, 'top_p': 1, 'seed': 7},
                  'model': {'seed': 7}}
NUM_CTX = 32768
NUM_PREDICT = 8192  # output cap per turn, so a repetition loop cannot run for minutes
MAX_CALLS = 12
REQUEST_TIMEOUT_S = 1800
SEED_INPUT = 'mcm_subset.txt'
ANSWER_FIELDS = ['model', 'guidance', 'question_id', 'kind', 'repeat', 'status', 'stop_reason',
                 'tool_calls', 'tools_used', 'tool_errors', 'wall_seconds', 'prompt_tokens',
                 'completion_tokens', 'max_prompt_tokens', 'context_near_limit',
                 'turns_hit_output_cap', 'final_answer']


# ---------------------------------------------------------------------------
# Ollama
# ---------------------------------------------------------------------------

def _http(host, path, body=None, timeout=60):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(f'{host.rstrip("/")}{path}', data=data,
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


def model_info(host, model) -> dict:
    show = _http(host, '/api/show', {'model': model})
    tags = {m['name']: m for m in _http(host, '/api/tags')['models']}
    details = show.get('details', {})
    return {'name': model, 'digest': tags.get(model, {}).get('digest'),
            'default_parameters': show.get('parameters'),
            'family': details.get('family'), 'parameter_size': details.get('parameter_size'),
            'quantisation': details.get('quantization_level'), 'format': details.get('format'),
            'capabilities': show.get('capabilities'),
            'context_length': next((v for k, v in show.get('model_info', {}).items()
                                    if k.endswith('.context_length')), None)}


def unload(host, model):
    try:
        _http(host, '/api/generate', {'model': model, 'keep_alive': 0})
    except (urllib.error.URLError, OSError):
        pass


# ---------------------------------------------------------------------------
# One conversation
# ---------------------------------------------------------------------------

def _ollama_tools(tools) -> list:
    return [{'type': 'function',
             'function': {'name': t.name, 'description': t.description or '',
                          'parameters': t.input_schema}} for t in tools]


def _result_text(result) -> str:
    return '\n'.join(c.text for c in result.content if getattr(c, 'text', None) is not None)


async def converse(host, model, question, env, sampling, max_calls=MAX_CALLS,
                   num_ctx=NUM_CTX, num_predict=NUM_PREDICT) -> dict:
    """Run one conversation. Returns the transcript fields it produces."""
    from mcp import Client, StdioServerParameters

    params = StdioServerParameters(command=sys.executable,
                                   args=['-m', 'aerosysprop.mcp_server'], env=env)
    calls, turns = [], []
    stop, final = None, ''
    t0 = time.time()
    async with Client(params) as client:
        # The server's own instructions follow the fixed prompt, as MCP clients
        # such as Claude Desktop pass them; they differ between guidance levels.
        instructions = client.instructions or ''
        system = SYSTEM_PROMPT + (f'\n\n{instructions}' if instructions else '')
        messages = [{'role': 'system', 'content': system},
                    {'role': 'user', 'content': question}]
        tools = (await client.list_tools()).tools
        tool_defs = _ollama_tools(tools)
        while True:
            body = {'model': model, 'messages': messages, 'tools': tool_defs, 'stream': False,
                    'options': {**sampling, 'num_ctx': num_ctx, 'num_predict': num_predict}}
            t = time.time()
            reply = await asyncio.to_thread(_http, host, '/api/chat', body, REQUEST_TIMEOUT_S)
            msg = reply.get('message', {})
            turns.append({'seconds': round(time.time() - t, 2),
                          'prompt_tokens': reply.get('prompt_eval_count'),
                          'completion_tokens': reply.get('eval_count'),
                          'done_reason': reply.get('done_reason')})
            messages.append(msg)
            requested = msg.get('tool_calls') or []
            if not requested:
                cut = reply.get('done_reason') == 'length'
                stop, final = ('output_cap' if cut else 'answered'), (msg.get('content') or '').strip()
                break
            if len(calls) + len(requested) > max_calls:
                stop, final = 'call_cap', (msg.get('content') or '').strip()
                break
            for tc in requested:
                fn = tc.get('function', {})
                name, args = fn.get('name', ''), fn.get('arguments') or {}
                if isinstance(args, str):
                    try:
                        args = json.loads(args)
                    except ValueError:
                        args = {'_unparsed': args}
                t = time.time()
                try:
                    res = await client.call_tool(name, args)
                    text, is_error = _result_text(res), bool(res.is_error)
                except Exception as e:  # unknown tool, protocol error: tell the model
                    text, is_error = f'Error: {type(e).__name__}: {e}', True
                calls.append({'name': name, 'arguments': args, 'is_error': is_error,
                              'seconds': round(time.time() - t, 2), 'result_chars': len(text)})
                tool_msg = {'role': 'tool', 'tool_name': name, 'content': text}
                if tc.get('id'):
                    tool_msg['tool_call_id'] = tc['id']
                messages.append(tool_msg)
    prompt = [x['prompt_tokens'] or 0 for x in turns]
    return {
        'server_instructions': instructions,
        'messages': messages, 'tools_offered': tool_defs, 'tool_calls': calls, 'turns': turns,
        'stop_reason': stop, 'final_answer': final,
        'wall_seconds': round(time.time() - t0, 2),
        'prompt_tokens': sum(prompt), 'max_prompt_tokens': max(prompt, default=0),
        'completion_tokens': sum(x['completion_tokens'] or 0 for x in turns),
        # Ollama truncates an over-long prompt silently; flag any turn that came close.
        'context_near_limit': max(prompt, default=0) >= 0.95 * num_ctx,
        'turns_hit_output_cap': sum(1 for x in turns if x['done_reason'] == 'length'),
    }


# ---------------------------------------------------------------------------
# Anthropic API
# ---------------------------------------------------------------------------
# Current Claude models think adaptively by default and reject non-default
# temperature, top_p and top_k, so the model's own defaults are the only
# setting, as with --sampling model for local models. There is no seed.
# Thinking is requested with display "summarized" so the transcript keeps a
# summary of it; that changes what is returned, not how the model reasons.
# Thinking blocks are passed back unmodified, as the API requires.
ANTHROPIC_THINKING = {'type': 'adaptive', 'display': 'summarized'}


def anthropic_client():
    try:
        import anthropic
    except ImportError:
        raise SystemExit('The anthropic backend needs the anthropic package: pip install anthropic')
    if not os.environ.get('ANTHROPIC_API_KEY'):
        raise SystemExit('Set ANTHROPIC_API_KEY to use the anthropic backend.')
    return anthropic, anthropic.Anthropic(max_retries=8)


def _anthropic_tools(tools) -> list:
    return [{'name': t.name, 'description': t.description or '',
             'input_schema': t.input_schema} for t in tools]


def _block_dict(b) -> dict:
    return b.model_dump(mode='json', exclude_none=True) if hasattr(b, 'model_dump') else dict(b)


async def converse_anthropic(api, model, question, env, effort=None, max_calls=MAX_CALLS,
                             max_tokens=NUM_PREDICT) -> dict:
    """One conversation through the Anthropic API. Same loop, prompt and caps as
    converse(); messages are also stored in the Ollama-like shape the markers read."""
    from mcp import Client, StdioServerParameters

    params = StdioServerParameters(command=sys.executable,
                                   args=['-m', 'aerosysprop.mcp_server'], env=env)
    calls, turns = [], []
    stop, final = None, ''
    served_model = None
    t0 = time.time()
    async with Client(params) as client:
        instructions = client.instructions or ''
        system = SYSTEM_PROMPT + (f'\n\n{instructions}' if instructions else '')
        tools = (await client.list_tools()).tools
        tool_defs = _anthropic_tools(tools)
        api_messages = [{'role': 'user', 'content': question}]   # what the API is sent
        messages = [{'role': 'system', 'content': system},       # what the markers read
                    {'role': 'user', 'content': question}]
        raw = []                                                  # every response, verbatim
        while True:
            body = {'model': model, 'max_tokens': max_tokens, 'system': system,
                    'messages': api_messages, 'tools': tool_defs,
                    'thinking': ANTHROPIC_THINKING}
            if effort:
                body['output_config'] = {'effort': effort}
            t = time.time()
            reply = await asyncio.to_thread(api.messages.create, **body)
            served_model = reply.model
            u = reply.usage
            cached = (getattr(u, 'cache_read_input_tokens', 0) or 0) + \
                     (getattr(u, 'cache_creation_input_tokens', 0) or 0)
            turns.append({'seconds': round(time.time() - t, 2),
                          'prompt_tokens': (u.input_tokens or 0) + cached,
                          'completion_tokens': u.output_tokens,
                          # 'length' as Ollama reports it, so the same checks apply
                          'done_reason': 'length' if reply.stop_reason == 'max_tokens'
                                         else reply.stop_reason,
                          'stop_reason': reply.stop_reason,
                          'request_id': getattr(reply, '_request_id', None)})
            blocks = [_block_dict(b) for b in reply.content]
            raw.append(blocks)
            api_messages.append({'role': 'assistant', 'content': reply.content})
            text = '\n'.join(b.text for b in reply.content if b.type == 'text').strip()
            thinking = '\n'.join(b.thinking for b in reply.content
                                 if b.type == 'thinking' and b.thinking).strip()
            requested = [b for b in reply.content if b.type == 'tool_use']
            msg = {'role': 'assistant', 'content': text}
            if thinking:
                msg['thinking'] = thinking
            if requested:
                msg['tool_calls'] = [{'id': b.id, 'function': {'name': b.name,
                                                               'arguments': b.input}}
                                     for b in requested]
            messages.append(msg)
            if not requested:
                stop = 'output_cap' if reply.stop_reason == 'max_tokens' else (
                    'refusal' if reply.stop_reason == 'refusal' else 'answered')
                final = text
                break
            if len(calls) + len(requested) > max_calls:
                stop, final = 'call_cap', text
                break
            results = []
            for b in requested:
                name, args = b.name, b.input or {}
                t = time.time()
                try:
                    res = await client.call_tool(name, args)
                    out, is_error = _result_text(res), bool(res.is_error)
                except Exception as e:  # unknown tool, protocol error: tell the model
                    out, is_error = f'Error: {type(e).__name__}: {e}', True
                calls.append({'name': name, 'arguments': args, 'is_error': is_error,
                              'seconds': round(time.time() - t, 2), 'result_chars': len(out)})
                results.append({'type': 'tool_result', 'tool_use_id': b.id,
                                'content': out, 'is_error': is_error})
                messages.append({'role': 'tool', 'tool_name': name, 'content': out,
                                 'tool_call_id': b.id})
            api_messages.append({'role': 'user', 'content': results})
    prompt = [x['prompt_tokens'] or 0 for x in turns]
    return {
        'server_instructions': instructions,
        'messages': messages, 'tools_offered': tool_defs, 'tool_calls': calls, 'turns': turns,
        'raw_responses': raw, 'served_model': served_model,
        'stop_reason': stop, 'final_answer': final,
        'wall_seconds': round(time.time() - t0, 2),
        'prompt_tokens': sum(prompt), 'max_prompt_tokens': max(prompt, default=0),
        'completion_tokens': sum(x['completion_tokens'] or 0 for x in turns),
        'context_near_limit': False,  # the context is far larger than any conversation here
        'turns_hit_output_cap': sum(1 for x in turns if x['done_reason'] == 'length'),
    }


# ---------------------------------------------------------------------------
# Campaign
# ---------------------------------------------------------------------------

def _slug(model) -> str:
    return re.sub(r'[^A-Za-z0-9.]+', '-', model).strip('-')


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git(*args) -> str:
    try:
        return subprocess.run(['git', *args], cwd=REPO, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ''


def _gpu() -> str:
    try:
        return subprocess.run(['nvidia-smi', '--query-gpu=name,driver_version',
                               '--format=csv,noheader'], capture_output=True, text=True,
                              timeout=20).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ''


def build_seed(out: Path, data_dir: Path, questions) -> Path:
    """The state every conversation starts from: one completed batch run."""
    seed = out / 'seed'
    if seed.is_dir() and any(seed.iterdir()):
        return seed
    batch_q = [q for q in questions if q['kind'] == 'batch']
    code = _seed_script(SEED_INPUT)
    seed.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, 'AEROSYSPROP_DATA_DIR': str(data_dir),
           'AEROSYSPROP_OUTPUT_DIR': str(seed)}
    p = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f'Could not build the seed run:\n{p.stderr[-2000:]}')
    if batch_q:
        print(f'seed run written to {seed}', flush=True)
    return seed


def _seed_script(input_file) -> str:
    # The same request build_questions makes for the batch keys: every vapour
    # pressure method over the subset, at the default temperature.
    return (
        'import time, warnings\n'
        'warnings.simplefilter("ignore")\n'
        'from aerosysprop import runs\n'
        f'r = runs.start_batch_run(input_path={input_file!r}, methods=["vapour_pressure"],\n'
        '                          label="earlier run")\n'
        'while runs.batch_run_status(r["run_id"])["status"] == "running":\n'
        '    time.sleep(0.2)\n'
        'assert runs.batch_run_status(r["run_id"])["status"] == "complete"\n')


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--models', required=True,
                   help='comma-separated Ollama tags, or Claude model IDs with --backend anthropic')
    p.add_argument('--backend', choices=['ollama', 'anthropic'], default='ollama')
    p.add_argument('--effort', choices=['low', 'medium', 'high', 'xhigh', 'max'], default=None,
                   help="anthropic only: thinking effort; omit for the model's default")
    p.add_argument('--conditions', default='full,bare')
    p.add_argument('--questions', default=str(HERE / 'questions.yaml'))
    p.add_argument('--questions-only', default='', help='comma-separated question ids')
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--out', required=True, help='campaign folder, e.g. evaluation/results/2026-10-07')
    p.add_argument('--data-dir', default=str(HERE / 'data'))
    p.add_argument('--host', default='http://localhost:11434')
    p.add_argument('--num-ctx', type=int, default=NUM_CTX)
    p.add_argument('--num-predict', type=int, default=NUM_PREDICT)
    p.add_argument('--sampling', choices=sorted(SAMPLING_MODES), default='model')
    p.add_argument('--max-calls', type=int, default=MAX_CALLS)
    p.add_argument('--order-seed', type=int, default=20261001)
    args = p.parse_args(argv)
    if args.backend == 'anthropic' and args.sampling != 'model':
        p.error('the Anthropic API accepts only the model\'s default sampling: use --sampling model')
    if args.effort and args.backend != 'anthropic':
        p.error('--effort applies to --backend anthropic only')
    # The MCP client logs every tool error; transcripts record them already.
    logging.getLogger('mcp').setLevel(logging.ERROR)

    models = [m.strip() for m in args.models.split(',') if m.strip()]
    conditions = [c.strip() for c in args.conditions.split(',') if c.strip()]
    qdoc = yaml.safe_load(open(args.questions, encoding='utf-8'))
    questions = qdoc['questions']
    if args.questions_only:
        wanted = [x.strip() for x in args.questions_only.split(',')]
        unknown = set(wanted) - {q['id'] for q in questions}
        if unknown:
            raise SystemExit(f'Unknown question ids: {", ".join(sorted(unknown))}')
        questions = [q for q in questions if q['id'] in wanted]

    out = Path(args.out).resolve()
    (out / 'transcripts').mkdir(parents=True, exist_ok=True)
    (out / 'runs').mkdir(exist_ok=True)
    data_dir = Path(args.data_dir).resolve()
    seed = build_seed(out, data_dir, questions)

    api = None
    if args.backend == 'anthropic':
        anthropic, api = anthropic_client()
        serving = {'software': 'anthropic-api', 'sdk_version': anthropic.__version__}
        server_label = f'anthropic api (sdk {anthropic.__version__})'
    else:
        try:
            ollama_version = _http(args.host, '/api/version')['version']
        except (urllib.error.URLError, OSError) as e:
            raise SystemExit(f'Cannot reach Ollama at {args.host}: {e}')
        serving = {'software': 'ollama', 'version': ollama_version, 'host': args.host}
        server_label = f'ollama {ollama_version}'
    from aerosysprop import __version__
    manifest_path = out / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    fixed = {'sampling_mode': args.sampling, 'num_ctx': args.num_ctx,
             'num_predict': args.num_predict, 'max_tool_calls': args.max_calls,
             'order_seed': args.order_seed, 'backend': args.backend, 'effort': args.effort}
    changed = {k: (manifest[k], v) for k, v in fixed.items() if k in manifest and manifest[k] != v}
    if manifest.get('question_set', {}).get('sha256') not in (None, _sha256(args.questions)):
        changed['question_set'] = ('a different file', args.questions)
    if changed:  # resuming must not mix settings within one campaign
        raise SystemExit(f'{out} was started with other settings: {changed}. '
                         'Use a new --out folder.')
    manifest.update({
        'started': manifest.get('started') or datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'aerosysprop': {'version': __version__, 'commit': _git('rev-parse', '--short', 'HEAD'),
                        'dirty': bool(_git('status', '--porcelain', '--untracked-files=no'))},
        'question_set': {'path': Path(os.path.relpath(args.questions, REPO)).as_posix(),
                         'sha256': _sha256(args.questions),
                         'specs_sha256': qdoc.get('specs_sha256'),
                         'keys_built_at_commit': qdoc.get('git_commit'),
                         'ids': [q['id'] for q in questions]},
        'conditions': conditions, 'repeats': args.repeats,
        'system_prompt': SYSTEM_PROMPT,
        'system_prompt_note': "followed by the server's instructions for the condition",
        'sampling_mode': args.sampling,
        'sampling': SAMPLING_MODES[args.sampling], 'num_ctx': args.num_ctx,
        'num_predict': args.num_predict,
        'max_tool_calls': args.max_calls, 'order_seed': args.order_seed,
        'backend': args.backend, 'effort': args.effort,
        'anthropic_thinking': ANTHROPIC_THINKING if args.backend == 'anthropic' else None,
        'serving': serving,
        'machine': {'platform': platform.platform(), 'machine': platform.machine(),
                    'gpu': _gpu(), 'python': platform.python_version()},
        'seed_state': sorted(p.name for p in seed.iterdir()),
    })
    manifest.setdefault('models', {})
    for m in models:
        manifest['models'][m] = (model_info(args.host, m) if args.backend == 'ollama' else
                                 {'name': m, 'provider': 'anthropic',
                                  'effort': args.effort or 'model default',
                                  'thinking': ANTHROPIC_THINKING,
                                  'sampling': 'model defaults (the API rejects others)'})
    manifest_path.write_text(json.dumps(manifest, indent=1))

    total = len(models) * len(conditions) * len(questions) * args.repeats
    done = 0
    for model in models:
        jobs = [(c, q, r) for c in conditions for q in questions for r in range(1, args.repeats + 1)]
        random.Random(f'{args.order_seed}:{model}').shuffle(jobs)
        info = manifest['models'][model]
        for cond, q, rep in jobs:
            done += 1
            stem = f'{_slug(model)}__{cond}__{q["id"]}__r{rep}'
            tpath = out / 'transcripts' / f'{stem}.json'
            if tpath.exists() and json.loads(tpath.read_text()).get('status') == 'complete':
                continue
            # Tools report the output folder's path, so its name must not give
            # away the model or the condition: an opaque name, mapped in the transcript.
            run_dir = out / 'runs' / hashlib.sha256(f'{args.order_seed}:{stem}'.encode()).hexdigest()[:12]
            if run_dir.exists():
                shutil.rmtree(run_dir)
            shutil.copytree(seed, run_dir)
            env = {'AEROSYSPROP_GUIDANCE': cond, 'AEROSYSPROP_DATA_DIR': str(data_dir),
                   'AEROSYSPROP_OUTPUT_DIR': str(run_dir), 'PYTHONWARNINGS': 'ignore'}
            record = {'model': model, 'quantisation': info.get('quantisation'),
                      'model_digest': info.get('digest'),
                      'server': server_label, 'backend': args.backend, 'guidance': cond,
                      'question_id': q['id'], 'kind': q['kind'], 'repeat': rep,
                      'question': q['question'].strip(),
                      'sampling_mode': args.sampling,
                      'sampling': {**SAMPLING_MODES[args.sampling], 'context': args.num_ctx,
                                   'num_predict': args.num_predict,
                                   'model_defaults': info.get('default_parameters')},
                      'max_tool_calls': args.max_calls,
                      'aerosysprop': manifest['aerosysprop'],
                      'started': datetime.now(timezone.utc).isoformat(timespec='seconds')}
            if args.backend == 'anthropic':
                record['effort'] = args.effort or 'model default'
                record['thinking'] = ANTHROPIC_THINKING
            try:
                if args.backend == 'anthropic':
                    record.update(asyncio.run(converse_anthropic(
                        api, model, q['question'].strip(), env, args.effort,
                        args.max_calls, args.num_predict)))
                else:
                    record.update(asyncio.run(converse(
                        args.host, model, q['question'].strip(), env,
                        SAMPLING_MODES[args.sampling], args.max_calls, args.num_ctx,
                        args.num_predict)))
                record['status'] = 'complete'
            except Exception as e:  # recorded, not fatal: the campaign carries on
                record.update({'status': 'error', 'error': f'{type(e).__name__}: {e}'[:2000]})
            new_runs = sorted({p.name for p in run_dir.iterdir()} - set(manifest['seed_state']))
            record['output_folder'] = Path(os.path.relpath(run_dir, out)).as_posix()
            record['new_run_folders'] = new_runs
            if not new_runs:
                shutil.rmtree(run_dir)
            tpath.write_text(json.dumps(record, indent=1, default=str))
            print(f'[{done}/{total}] {stem}: {record["status"]} '
                  f'{record.get("stop_reason") or record.get("error", "")[:80]} '
                  f'calls={len(record.get("tool_calls", []))} {record.get("wall_seconds", "")}s',
                  flush=True)
        if args.backend == 'ollama':
            unload(args.host, model)

    write_answers(out)
    manifest['finished'] = datetime.now(timezone.utc).isoformat(timespec='seconds')
    manifest_path.write_text(json.dumps(manifest, indent=1))
    return 0


def write_answers(out: Path):
    rows = []
    for t in sorted((out / 'transcripts').glob('*.json')):
        r = json.loads(t.read_text())
        calls = r.get('tool_calls', [])
        rows.append({
            'model': r['model'], 'guidance': r['guidance'], 'question_id': r['question_id'],
            'kind': r['kind'], 'repeat': r['repeat'], 'status': r['status'],
            'stop_reason': r.get('stop_reason'), 'tool_calls': len(calls),
            'tools_used': ';'.join(c['name'] for c in calls),
            'tool_errors': sum(c['is_error'] for c in calls),
            'wall_seconds': r.get('wall_seconds'), 'prompt_tokens': r.get('prompt_tokens'),
            'completion_tokens': r.get('completion_tokens'),
            'max_prompt_tokens': r.get('max_prompt_tokens'),
            'context_near_limit': r.get('context_near_limit'),
            'turns_hit_output_cap': r.get('turns_hit_output_cap'),
            'final_answer': r.get('final_answer', '')})
    with open(out / 'answers.csv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=ANSWER_FIELDS)
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    sys.exit(main())
