"""Check that a locally served model can call tools at all.

A model that describes a tool call in prose, rather than returning it in the
structured slot its chat template provides, cannot take part in the evaluation.
That is worth finding out in one second rather than after a long run.

Uses only the standard library, and talks to Ollama's chat endpoint. Nothing
here touches AeroSysProp: the tool below is a stand-in with an obvious answer,
so the only thing being tested is the model's tool-calling machinery.

    python evaluation/check_tool_calling.py --model qwen3:8b
    python evaluation/check_tool_calling.py --model llama3.1:8b --host http://localhost:11434
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request

TOOLS = [{
    'type': 'function',
    'function': {
        'name': 'get_vapour_pressure',
        'description': ('Return the vapour pressure of a molecule in log10(atm). '
                        'Takes a SMILES string. Always use this tool rather than '
                        'answering from memory.'),
        'parameters': {
            'type': 'object',
            'properties': {
                'smiles': {'type': 'string', 'description': 'SMILES, e.g. CCO for ethanol'},
                'temperature_K': {'type': 'number', 'description': 'temperature in kelvin'},
            },
            'required': ['smiles'],
        },
    },
}]

QUESTION = 'What is the vapour pressure of ethanol (SMILES CCO) at 298.15 K?'
ANSWER = {'smiles': 'CCO', 'temperature_K': 298.15, 'vapour_pressure_log10atm': -1.0189}


def chat(host: str, model: str, messages: list, tools=None, timeout=180) -> dict:
    body = {'model': model, 'messages': messages, 'stream': False,
            'options': {'temperature': 0}}
    if tools:
        body['tools'] = tools
    req = urllib.request.Request(f'{host.rstrip("/")}/api/chat',
                                 data=json.dumps(body).encode(),
                                 headers={'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--model', required=True, help='an Ollama model tag, e.g. qwen3:8b')
    p.add_argument('--host', default='http://localhost:11434')
    p.add_argument('--verbose', action='store_true', help='print the raw replies')
    args = p.parse_args(argv)

    messages = [{'role': 'user', 'content': QUESTION}]
    try:
        first = chat(args.host, args.model, messages, TOOLS)
    except urllib.error.HTTPError as exc:  # before URLError, which it subclasses
        detail = exc.read().decode(errors='replace')[:200]
        print(f'FAIL  Ollama refused the request ({exc.code}): {detail}')
        print(f'      has the model been pulled? try: ollama pull {args.model}')
        return 2
    except urllib.error.URLError as exc:
        print(f'FAIL  cannot reach Ollama at {args.host}: {exc}')
        print('      is the Ollama service running? try: ollama list')
        return 2

    msg = first.get('message', {})
    if args.verbose:
        print(json.dumps(first, indent=1)[:2000])
    calls = msg.get('tool_calls') or []

    if not calls:
        text = (msg.get('content') or '').strip()
        print(f'FAIL  {args.model} returned no tool call.')
        if 'get_vapour_pressure' in text:
            print('      It described the call in prose instead of the structured slot,')
            print('      which the harness cannot parse. This model cannot take part.')
        else:
            print('      It answered from memory. Model answered:')
        print(f'      {text[:300]}')
        return 1

    call = calls[0]['function']
    raw_args = call.get('arguments')
    parsed = json.loads(raw_args) if isinstance(raw_args, str) else (raw_args or {})
    print(f'ok    tool call: {call.get("name")}({json.dumps(parsed)})')
    if call.get('name') != 'get_vapour_pressure':
        print('WARN  called a tool that does not exist')
        return 1
    if 'CCO' not in str(parsed.get('smiles', '')):
        print('WARN  the SMILES argument does not look right')

    # Second turn: hand back a result and see whether the model uses it.
    messages += [msg, {'role': 'tool', 'name': 'get_vapour_pressure',
                       'content': json.dumps(ANSWER)}]
    second = chat(args.host, args.model, messages, TOOLS)
    answer = (second.get('message', {}) or {}).get('content', '').strip()
    if args.verbose:
        print(json.dumps(second, indent=1)[:2000])
    print(f'ok    final answer: {answer[:300]}')
    if '-1.0' in answer or '1.0189' in answer or '1.02' in answer:
        print(f'PASS  {args.model} calls tools and uses the result. Ready for the evaluation.')
        return 0
    print(f'WARN  {args.model} called the tool but did not report the value it was given.')
    print('      Usable, but expect it to score badly on reporting what the tool returned.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
