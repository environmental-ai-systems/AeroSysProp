# Evaluation tools: testing AI agents that use AeroSysProp

These scripts measure what an AI model does when it answers questions through the AeroSysProp MCP server, and whether the documentation layer served to agents changes that. The design of the study they were written for is described in [docs/agent-evaluation.md](../docs/agent-evaluation.md). The study's own question set and results are not here; they will be released with the paper describing AeroSysProp. Everything needed to write, key, run and mark a question set of your own is.

| File | What it does |
| --- | --- |
| `example_specs.yaml` | A worked example question, and the template for your own set |
| `build_questions.py` | Runs AeroSysProp on a specification and writes the keyed question set |
| `check_tool_calling.py` | Checks that a model can make structured tool calls at all |
| `run_eval.py` | Asks every question of every model at both guidance levels and saves the transcripts |
| `mark.py` | Automatic marking of the transcripts against the keys |
| `judge.py` | Sends the answers the automatic checks cannot decide to a judge model, blind |
| `analyse.py` | Paired comparison of the two guidance levels |
| `make_session_db.py` | Writes a small web-app session database for questions about saved sessions |
| `data/mcm_subset.txt` | 583 MCM species (every tenth), the input for batch questions |

## A worked example

`example_specs.yaml` holds one question: the vapour pressure of MBTCA (3-methylbutane-1,2,3-tricarboxylic acid) at 298.15 K from SIMPOL.1 and EVAPORATION, and how far apart the two estimates are. The specification says what the question asks, which molecule and methods it uses, and what the answer must report. It contains no expected value. Build the key:

```bash
pip install -e ".[all]"
python -m evaluation.build_questions --specs evaluation/example_specs.yaml --skip-batch
```

This writes `evaluation/example_questions.yaml`. The key holds the values AeroSysProp returned, the diagnostics (none for this molecule), the method versions with their parameter fingerprints, the tolerance and what must be reported:

```yaml
key:
  values:
    MBTCA:
      vp_evaporation_log10atm__set2: -10.7215
      vp_simpol_log10atm: -12.0801
  diagnostics: {}
  provenance_labels:
    vp_evaporation@2: vp_evaporation@2 (d3935c1d1b7a)
    vp_simpol: vp_simpol@1 (44ca9771f3fb)
  tolerance: 0.02
  must_report: [value_with_units, method_named]
  value_from_transcript: true
```

The two methods differ by about 1.4 orders of magnitude for this molecule, so a good answer gives both values with units, names the methods and states the gap. An answer that quotes one number, or rounds both until the gap disappears, is marked down. `value_from_transcript` means every number in the answer is also checked against what the tools actually returned in that conversation, so a value no tool produced counts as invented.

Then run it against a model and mark it. One repeat is enough for a first look:

```bash
python evaluation/check_tool_calling.py --model qwen3:8b
python -m evaluation.run_eval --models qwen3:8b --questions evaluation/example_questions.yaml \
    --repeats 1 --out evaluation/results/example
python -m evaluation.mark evaluation/results/example
```

`marks.csv` in the campaign folder has one row per conversation, with each automatic check and a `blind` column naming anything that still needs a person or a judge model.

## Writing your own set

Copy the example into `evaluation/sets/`, which git ignores, so your questions and keys stay yours:

```bash
mkdir -p evaluation/sets
cp evaluation/example_specs.yaml evaluation/sets/my_specs.yaml
# edit the molecules, methods and wording, add questions, then:
python -m evaluation.build_questions --specs evaluation/sets/my_specs.yaml --skip-batch
```

The builder writes `evaluation/sets/my_questions.yaml`. You never type an expected value: re-run the builder after any change to the methods and every key updates together. If a key moves when nothing was meant to change, that is worth knowing.

Eight kinds of question are supported:

| Kind | What the key holds | How it is marked |
| --- | --- | --- |
| `lookup` | The values returned for named methods | Value within tolerance, method named, units |
| `delta` | Parent and variants, the change per method, which method moves most | Direction, size, whether the spread between methods was reported |
| `series` | The step per added carbon, per method | Step within tolerance, whether linearity was discussed |
| `diagnostic` | The diagnostic codes the tool reports | Whether the warning reached the answer, and no false warnings |
| `mixture` | The AIOMFAC result | Coefficients within tolerance, basis named |
| `batch` | The run summary, coverage and any exact lookups | Batch run used, folder cited, coverage stated |
| `open` | Acceptable methods and what must be reported | Automatic checks plus blind marking |
| `refusal` | What the model should decline to do | Blind marking |

A `delta` question takes a `parent` molecule and a list of `variants`, each with a `change` describing what was added. A `series` question takes a list of `members`, a `mixture` question takes `components`, `mass_fractions` and optionally `salts`, and a `batch` question takes an `input_file` in the data folder. Batch questions need `AEROSYSPROP_DATA_DIR` and `AEROSYSPROP_OUTPUT_DIR` set; leave out `--skip-batch` to build them. Questions about a saved web-app session need `python -m evaluation.make_session_db` first. For questions that need judgement (`open`, `refusal`, or any with `credit_if`), a `judge_reference` field gives the judge model the facts a good answer should rest on; `{card_facts}` and `{evaporation_sets}` in it insert the standard notes on method cards and on the EVAPORATION parameter sets.

Questions whose ID starts with `prov-` ask about method provenance, and `mark.py` credits an answer read from `list_methods` as well as from a prediction.

A comparison of two models on ten questions about your own compounds, or of one model before and after you change a tool description, is a quick way to find out whether an agent can be trusted with them.

## Running a full campaign

`run_eval.py` asks every question of every model at both guidance levels (`--conditions full,bare`), three times each by default (`--repeats`), in an order shuffled per model. Each conversation gets a fresh MCP server with `AEROSYSPROP_GUIDANCE` set, the same system prompt followed by the server's own instructions for that level, and at most 12 tool calls. Every conversation starts from an output folder already holding one completed batch run over `mcm_subset.txt`, so questions that refer to "the last run" mean the same thing whatever order they come in.

Local models are served by [Ollama](https://ollama.com). Sampling defaults to each model's recommended settings with a fixed seed (`--sampling model`), because greedy decoding makes some thinking models loop; `--sampling greedy` uses temperature 0 instead. A model that fails `check_tool_calling.py` cannot take part, and that is recorded as a result: the harness never repairs a model's arguments.

Claude models run through the Anthropic API in the same loop, with the same tools and no web search:

```bash
pip install anthropic
export ANTHROPIC_API_KEY=...        # PowerShell: $env:ANTHROPIC_API_KEY = "..."
python -m evaluation.run_eval --backend anthropic --models claude-sonnet-5-5 \
    --questions evaluation/sets/my_questions.yaml --out evaluation/results/claude
```

Run each backend into its own campaign folder. The harness writes `manifest.json` (software version and commit, question set checksum, models with their digests and quantisation, sampling), one transcript per conversation, and `answers.csv`. Stopping and restarting the same command carries on where it left off; restarting with different settings is refused.

After `mark.py`, `judge.py` sends the undecided answers to a judge model with the model and guidance level hidden, and `analyse.py` gives the paired comparison of the two levels:

```bash
python -m evaluation.judge export  evaluation/results/claude   # blind, shuffled items
python -m evaluation.judge preview evaluation/results/claude --n 2
python -m evaluation.judge submit  evaluation/results/claude
python -m evaluation.judge collect evaluation/results/claude
python -m evaluation.judge merge   evaluation/results/claude
python -m evaluation.analyse evaluation/results/claude
```

A judge model should not be the only check on answers from its own family: have a person, or a second judge, look at those.

## Tests

```bash
python -m pytest evaluation/test_marking.py
```
