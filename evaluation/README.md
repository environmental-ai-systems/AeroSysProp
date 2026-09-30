# Evaluation: does the guidance layer change what a model does?

Two files start the question set:

- `question_specs.yaml` says what each question tests and which molecules it uses. It holds **no expected values**.
- `build_questions.py` runs AeroSysProp and writes `questions.yaml`, with values, diagnostics, provenance labels and fingerprints filled in.

Keys are therefore generated, never typed. Re-run the builder after any approved
parameter change and every key updates together. If a key moves when nothing was
meant to change, that is a finding rather than a nuisance.

## Building the set

```bash
pip install -e ".[all]"
python -m evaluation.build_questions --skip-batch      # no run folder needed
```

For the batch questions, set the two folders first and drop the species file in
the input folder:

```bash
export AEROSYSPROP_DATA_DIR=~/aerosysprop-inputs
export AEROSYSPROP_OUTPUT_DIR=~/aerosysprop-runs
python -m evaluation.build_questions
```

`mcm_subset.txt` is a tab separated `name<TAB>SMILES` file. A few hundred species
is enough; taking every tenth species from the MCM gives a mixture of clean
molecules, radicals and flagged cases, which is what the coverage questions need.

## Kinds of question

| Kind | What the key holds | How it is marked |
| --- | --- | --- |
| `lookup` | The values returned for named methods | Value within tolerance, method used, units |
| `delta` | Parent and variants, the change per method, which method moves most | Direction, magnitude, whether the spread was reported |
| `series` | The step per added carbon, per method | Step within tolerance, whether linearity was discussed |
| `diagnostic` | The diagnostic codes the tool reports | Whether the caveat reached the answer, and no false warnings |
| `mixture` | The AIOMFAC result | Coefficients within tolerance, basis named |
| `batch` | The run summary, coverage and any exact lookups | Batch run used, folder cited, coverage stated |
| `open` | Acceptable methods and what must be reported | Automatic checks plus blind marking |
| `refusal` | What the model should decline to do | Blind marking |

Every key carries `value_from_transcript: true`. The marker checks the answer
against the tool result in that transcript as well as against the key, so a
number no tool returned is caught as invention whatever path the model took.

## Running the evaluation

Models are served by [Ollama](https://ollama.com) on the same machine. First
check each model can call tools with list arguments, as AeroSysProp's tools need:

```bash
python evaluation/check_tool_calling.py --model qwen3:8b
```

A model that fails this check cannot take part, and that is reported as a
result rather than worked around: the harness never repairs a model's
arguments. On the Spark (Ollama 0.35.0), `llama3.1:8b` fails it by sending
lists as strings such as `"['CCO']"`.

Then run the campaign and mark it:

```bash
python -m evaluation.run_eval --models qwen3:8b,granite4.1:8b --out evaluation/results/2026-10-07
python -m evaluation.mark evaluation/results/2026-10-07
```

`run_eval.py` asks every question of every model under both guidance levels
(`--conditions full,bare`), three times each (`--repeats`), in an order shuffled
per model. Each conversation gets a fresh MCP server with `AEROSYSPROP_GUIDANCE`
set, the fixed system prompt followed by the server's own instructions for that
level (as MCP clients such as Claude Desktop pass them), and up to 12 tool
calls. Sampling defaults to each model's own recommended settings with a fixed
seed (`--sampling model`), because greedy decoding makes some thinking models
loop; `--sampling greedy` pins temperature 0 instead. It writes
`manifest.json` (software, commit, question set checksum, models with their
digests and quantisation, sampling), one transcript per conversation, and
`answers.csv`. Stopping and restarting the same command carries on where it
left off; restarting with different settings is refused.

Every conversation starts from the same state: an output folder already holding
one completed run of the batch request over `mcm_subset.txt`, so questions that
refer to "that run" or "the last batch" mean the same thing whatever order they
come in. Folders are given opaque names because tools report their paths.

`mark.py` writes `marks.csv`: the automatic checks (a tool was called, values
against the key, method, version named, caveat passed on, run folder cited,
coverage stated, numbers no tool returned), with a `blind` column listing what
still needs a human or judge.

For a quick test, `--questions-only lookup-01,diag-01 --repeats 1`.

## Your own question set

The shipped set is a reference, not a limit. To test the compounds you care
about, copy the template and point the builder at it:

```bash
mkdir -p evaluation/sets
cp evaluation/example_specs.yaml evaluation/sets/my_specs.yaml
# edit the molecules and questions, then:
python -m evaluation.build_questions --specs evaluation/sets/my_specs.yaml --skip-batch
# writes evaluation/sets/my_questions.yaml
```

You never type an expected value: the builder computes the keys the same way it
does for the shipped set, and stamps them with the commit and the parameter
fingerprints they came from. Run the harness against your set exactly as against
the shipped one.

`evaluation/sets/` and `evaluation/results/` are ignored by git, so your
questions and your runs stay yours. The shipped specification, its keyed set and
the builder are listed in `docs/CODEOWNERS`, because they define what "the tool
answers well" means for everyone; your own sets need nobody's approval.

A local comparison is worth running for its own sake, separately from this study.
Two models, or one model before and after you change a tool description, on ten
questions about your own compounds, tells you quickly whether an agent can be
trusted with them.

## Before the set is frozen

1. The core team reviews and edits every question. Claude drafted them, and
   Claude is one of the models under test.
2. The team adds a held-back set of its own, which Claude never sees. If the
   results are close, those are the questions to report most prominently.
3. Freeze: record the `specs_sha256` from `questions.yaml` in the run manifest.

The `delta` and `series` questions are the scientifically interesting ones. They
ask what happens to a vapour pressure when one functional group is added or a
chain is lengthened, and whether the answer checks that change against a method
built to capture it. EVAPORATION carries interaction terms for groups on
adjacent carbons, SIMPOL is additive, and the boiling-point-dependent methods
inherit the error of whichever boiling point fed them. The two hexanediols in
`delta-diol-position-03` differ only in where the second hydroxyl sits, so a
purely additive method cannot separate them at all: that question exists to see
whether a model notices.
