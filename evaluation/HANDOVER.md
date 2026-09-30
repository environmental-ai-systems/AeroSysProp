# Handover: running the agent evaluation on the Blackwell desktop (1 Oct 2026)

Read this first if you are a fresh session picking up the evaluation. It
replaces the handover written for the DGX Spark on 30 Sep 2026. The protocol is
the Claude doc "Agent evaluation protocol"
(https://claude.ai/code/artifact/9d095b52-c622-4b9f-bfbb-f8bd89828998); where it
and this file disagree, this file records the later decision (see "Decisions").

## What we are testing

Whether the documentation layer in AeroSysProp changes what an AI model does
with it. One install serves two conditions through `AEROSYSPROP_GUIDANCE=full`
or `bare`. Values, diagnostics, provenance labels, error messages and run
folders are identical in both, so a difference in the answers can only have
entered through the documentation. The second question is who benefits: do
small local models gain most from good self-description?

## Where things stand

Everything below is built, tested on the Spark and pushed.

Repository `environmental-ai-systems/AeroSysProp`:

- **The guidance switch** is on `main` (pull request #5, merged 1 Oct 2026).
  At `bare`: one-line tool descriptions, no titles or argument descriptions
  (types, defaults, enums and bounds kept, so validation is identical), one-line
  server instructions, no guide resources, `list_methods` gives IDs, output
  columns and units only, `predict_properties` drops its `interpretation`
  paragraph, and `batch_run_status` drops the coverage statement, caveats and
  reporting requirements (`run.json` always keeps them). An unknown level stops
  the server. `aerosysprop/tests/test_guidance.py` compares the levels byte for
  byte, including through a real server process per level; 4,001 tests pass.
- **`agent-eval`**: the question set, keys and builder, with the switch
  merged in (the same commit as on `main`, so no conflicts later), plus:
  - `evaluation/run_eval.py`, the harness (see `evaluation/README.md`).
  - `evaluation/mark.py`, the automatic marking.
  - `evaluation/check_tool_calling.py`, now with a list argument like the real
    tools, which is what catches `llama3.1:8b`.

Work on `agent-eval`; do not commit to `main`.

## What the Spark established

- AeroSysProp runs on ARM (aarch64, GB10): 3,983 of the existing tests pass,
  RDKit 2026.03.6 from conda-forge.
- **The keys travel between architectures.** A full rebuild on ARM matched the
  x86-built `questions.yaml` in all 899 numbers to within 7e-16 (six AIOMFAC
  coefficients differ in the last digit; everything else is bit-identical). The
  0.02 log10 atm tolerance stands. Any exact-equality comparison of mixture keys
  across machines would fail; tolerance comparisons pass.
- Smoke tests of the harness found and fixed: the condition leaking through the
  output folder path (folders now have opaque names), a 12,000-token repetition
  loop from `qwen3:4b` at temperature 0 (now an 8,192-token cap per turn, and
  model-default sampling), and two marking false alarms.
- Early signs the questions discriminate (single runs, not results): `qwen3:4b`
  told the user it had overwritten a run folder, which the tool cannot do
  (`refuse-01`); and gave up after an invalid-SMILES error while stating the
  correct SMILES itself (`diag-01`, bare).

## Models

Sixteen take part; all pass `check_tool_calling.py` on the Spark (Ollama
0.35.0). IDs are the Spark's `ollama list`: compare them on the desktop, and
record any that differ (an upstream update).

| Model | ID | Size | Check |
| --- | --- | --- | --- |
| qwen3:4b | 359d7dd4bcda | 2.5 GB | pass |
| qwen3:8b | 500a1f067a9f | 5.2 GB | pass (did not repeat the value once) |
| qwen3:14b | bdbd181c33f2 | 9.3 GB | pass |
| qwen3:32b | 030ee887880f | 20 GB | pass |
| qwen3:30b-a3b | ad815644918f | 18 GB | pass |
| qwen3.8:27b | aaee06c39dcf | 17 GB | pass |
| gpt-oss:20b | 17052f91a42e | 13 GB | pass |
| llama3.1:70b | 711a9e8463af | 42 GB | pass |
| granite4.1:3b | 6fd349357287 | 2.1 GB | pass |
| granite4.1:8b | 444af1c4b2fe | 5.3 GB | pass (did not repeat the value once) |
| granite4.1:30b | 3f3e5df8a021 | 17 GB | pass |
| gemma4:12b | 6114515d63c1 | 8.0 GB | pass |
| gemma4:26b | 001e5dafc3c7 | 18 GB | pass |
| gemma4:31b | 17ba34c06c80 | 20 GB | pass |
| mistral-small3.2:24b | 5a408ab55df5 | 15 GB | pass |
| lfm2.5:8b | 9cf756159fc2 | 5.2 GB | pass |

Excluded: `llama3.1:8b` (46e0c10c039e) sends list arguments as strings such as
`"['CCO']"`; the server rejects them and the model never corrects itself.
Report it as unable to take part. Claude is also to be tested, but the harness
only talks to Ollama so far (see "Still to do").

## Setting up the desktop (x86 Linux, RTX Pro 6000 Blackwell, 96 GB)

The whole full run happens on this one machine, so hardware does not confound
the size ladder. Nothing needs copying from the Spark.

1. `nvidia-smi` shows the card with a CUDA 12.8+ driver. Install Ollama
   (`curl -fsSL https://ollama.com/install.sh | sh`, needs sudo, so the user
   runs it) and check `journalctl -u ollama` shows an "inference compute …
   library=CUDA" line naming the card. Prefer Ollama 0.35.0 to match the Spark.
2. Miniforge to `~/miniforge3`, then
   `~/miniforge3/bin/conda create -n aerosysprop -c conda-forge python=3.11 rdkit`.
   Nothing in the system Python; call `~/miniforge3/envs/aerosysprop/bin/python`
   directly.
3. GitHub access for that machine (its own SSH key, or `gh auth login`), clone,
   `git switch agent-eval`, `python -m pip install -e ".[all]" pytest`. Set
   `git config user.name/user.email` before committing.
4. Tests: `pytest -q aerosysprop/tests --ignore=aerosysprop/tests/test_main.py
   --ignore=aerosysprop/tests/test_vapour_pressure.py` (4,001 pass). Keys:
   `AEROSYSPROP_DATA_DIR=$PWD/evaluation/data AEROSYSPROP_OUTPUT_DIR=~/aerosysprop-runs
   python -m evaluation.build_questions`; only build stamps, batch run IDs and
   last-digit AIOMFAC values should differ. Then `git restore evaluation/questions.yaml`.
5. Pull the 16 models, compare IDs with the table, and run
   `python evaluation/check_tool_calling.py --model <tag>` on each: the Spark's
   results do not carry over to different hardware.
6. Pilot, one repeat, then look at transcripts and marks before going further:
   `python -m evaluation.run_eval --models <tags> --repeats 1 --out evaluation/results/pilot`
   then `python -m evaluation.mark evaluation/results/pilot`.
7. Full run with `--repeats 3` into a dated folder under `evaluation/results/`,
   nothing else on the GPU. It resumes if interrupted and refuses to mix
   settings. Back the folder up afterwards: git ignores it.

Timing: on the Spark the full run was estimated at about four days (the dense
thinking models dominate; `qwen3:32b` generates at about 10 tokens/s and
`llama3.1:70b` at about 5). The desktop's memory bandwidth is about six times
the Spark's, so expect roughly a day; measure it in the pilot. The smoke test
took 30 to 90 s per conversation for `qwen3:4b`, most of it thinking.

## Decisions

Earlier ones still stand: keys are generated, never typed; prescribed and open
questions are marked differently; Claude drafted the questions and is under
test, so the core team reviews them and a held-back set by David and Simon is
reported prominently; applicability sections for method cards are deferred and
would become a third guidance level.

Made on the Spark, 30 Sep to 1 Oct 2026:

- **Sampling:** each model's own recommended settings with a fixed seed
  (`--sampling model`, the default), not temperature 0. Qwen warns greedy
  decoding makes its thinking models loop, and it did. The protocol doc still
  says temperature 0 and needs updating.
- **System prompt:** the fixed prompt followed by the server's instructions for
  the level, as MCP clients pass them. The protocol doc says the prompt is
  identical in both conditions and needs updating. The guide resources are not
  reachable through the harness in either condition.
- **Seeded state:** every conversation starts with one completed batch run in
  its output folder, so "that run" (`batch-03` to `batch-06`) and "the last
  batch" (`refuse-01`) mean the same thing whatever the order.
- **Tool-calling failures are results:** the harness never repairs a model's
  arguments.
- **Where to run:** the full run on the Blackwell desktop, with the machine to
  itself; the Spark was for setup and debugging.

## Still to do

1. Update the protocol doc for the decisions above.
2. A Claude backend in `run_eval.py` (Anthropic API in the same loop; needs an
   API key), or run Claude through Claude Desktop, and record which in the
   manifest.
3. Blind marking: `marks.csv` lists what needs it in the `blind` column
   (refusals, justification, framing, possible invention). Not built yet: a
   shuffled, model- and condition-hidden export for markers or a judge model.
4. Carried over: core team review of the 46 questions (pinic and levulinic acid
   SMILES especially); the held-back question set; whether batch runs default
   to a folder when `AEROSYSPROP_OUTPUT_DIR` is unset (currently off); branch
   protection on `main` (waiting on GitHub Education); EVAPORATION set 2's diol
   correction on 13 MCM alpha-hydroxy acids (a parameter set 3 is a core team
   decision).
