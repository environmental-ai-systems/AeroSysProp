# Evaluating AI agents that use AeroSysProp

AeroSysProp can be called by AI agents through its MCP server (see "Use with AI agents (MCP)" in the [README](../README.md)). Every tool comes with a documentation layer written for agents. This covers tool descriptions, an agent guide, method cards, and the caveats and reporting requirements in batch run summaries. This page describes a study of whether that layer changes what a model does with the tool, and which models benefit from it most.

## The question

An agent that calls a property method can still give a poor answer. It might quote a value for a molecule the method could not fully describe, choose a method that does not suit the compound, round a vapour pressure until the difference between two methods disappears, or report a number that no tool returned. The documentation layer exists to prevent these errors. The study asks how much it actually helps, and whether small models that run on local hardware gain more from it than large hosted models.

## Two guidance levels, identical values

The server reads one environment variable at start-up. `AEROSYSPROP_GUIDANCE=full`, the default, serves the whole documentation layer. `AEROSYSPROP_GUIDANCE=bare` serves tool names, argument names and units only, with no guide resource, no method card content in `list_methods`, and no coverage statement, caveats or reporting requirements in batch summaries. Values, diagnostics, provenance labels, error messages and run folders are identical at both levels, and the test suite checks this byte for byte. Every question is asked of every model at both levels, so any difference between the answers can only have come from the documentation. For normal use, leave the variable unset.

## Questions and keys

The question set covers the ways an agent can go wrong with a property tool. There are eight kinds of question:

- lookups of a named method's value;
- the change in a property when one functional group is added to a molecule;
- the step per carbon along a homologous series;
- molecules whose functional groups a method cannot fully assign, where the diagnostic has to reach the answer;
- AIOMFAC activity coefficients for mixtures with salts;
- batch runs, where the run folder and the coverage of the molecule set must be reported;
- open questions, where the model chooses the method itself;
- requests the model should decline.

No expected value is typed by hand. A builder script runs AeroSysProp on the specification of each question and writes the key, stamped with the software commit and the SHA-256 fingerprint of every parameter set used. If an approved change to a method's parameters alters a key, every key is rebuilt together. Keys for the questions that do not need a batch run, built on x86 and ARM machines, agree to within 1e-15.

## Running and marking

The harness starts a fresh MCP server for each conversation and gives every model the same system prompt and the same tools. It allows at most 12 tool calls, and every question is repeated three times in a shuffled order. Local models are served through Ollama on a single machine, so hardware does not confound comparisons between model sizes. Claude models are run through the Anthropic API in the same loop, with no web search. Each transcript records the model, its quantisation and digest, the sampling settings, the guidance level, the software version and commit, timings and token counts.

Marking happens in two stages. Automatic checks come first. They confirm that a tool was called and that values match the key within a stated tolerance (0.02 log10 units for vapour pressures). They also check that the method and its parameter set were named, that diagnostics and caveats were passed on, and that the run folder and coverage were reported. Any number that no tool returned in that conversation is treated as invented. A value rounded more coarsely than the tolerance is not credited. Answers that need judgement, such as the reasoning behind a choice of method or a refusal, are then marked by a judge model or a person, with the model and the guidance level hidden and the answers shuffled.

## Safeguards

Claude drafted the first version of the questions and is also one of the models tested. The core team therefore reviews every question before the set is frozen, and writes a held-back set that Claude never sees. That set stays private so that it remains unseen by any model. A Claude model is also one of the judges, so Claude's own answers are checked by a person or a second judge. A model that cannot make a structured tool call is recorded as unable to take part, and its calls are never repaired.

## Your own questions

The evaluation tools also work with a question set you write yourself about the compounds you care about. The builder computes the keys for your set exactly as for the shipped set, so you never type an expected value. Two models compared on ten questions about your own molecules, or one model before and after a change to a tool description, quickly shows whether an agent can be trusted with them.

## Availability

The harness, the marking code and the question set will be released together with the paper describing AeroSysProp. Until then they are kept private, so that the questions do not reach models in training before the results are reported. The guidance switch is part of the released package from v1.0.0.
