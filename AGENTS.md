# AGENTS.md

Guidance for AI coding agents working on the AeroSysProp codebase. (Agents that want to *use* AeroSysProp should connect to the MCP server instead; see `aerosysprop/agent_guide.md`.)

## Layout

* `aerosysprop/methods/`: the published property methods. Each must remain a faithful implementation of its source paper; any deviation is documented in the docstring.
* `aerosysprop/batch.py`: `PROPERTY_METHODS` is the single registry of single-component methods. Adding a method means one `PropertyMethod` entry plus its wiring in `_run_one`.
* `aerosysprop/diagnostics.py`: pre-flight checks. `_METHOD_REGISTRY` maps each method to its SMARTS set.
* `aerosysprop/api.py`: the agent-facing layer (plain JSON in and out). Add new capabilities here first.
* `aerosysprop/mcp_server.py`: MCP tools, each a thin wrapper over one `api` function.
* `aerosysprop/app.py`, `views/`: Streamlit app.

## Rules

* `api.py` never writes. The only writing path is `runs.py`, which creates a new run folder inside `AEROSYSPROP_OUTPUT_DIR` per run and never overwrites or deletes; keep it that way. Never pass `output_csv` to `run_all` from `api.py`. Database reads use `read_only=True`; nothing writes to a session database.
* Worker processes started from the API or the MCP server must use the `spawn` context (`api._SPAWN`); forking inside a multi-threaded server is unsafe.
* Batch-run summaries must keep the coverage statement, caveats and clean/all split; anything an agent sees in place of the full results must say what it leaves out.
* File access from the API goes through `_resolve_in_root`, which confines it to `AEROSYSPROP_DATA_DIR`.
* Nothing in `mcp_server.py` may print to stdout (it carries the stdio protocol).
* API outputs must be JSON-safe: pass values through `json_safe`, never return RDKit objects, numpy types, tuples or NaN.
* When adding a single-component method, also add it to `_DIAG_PATH` in `api.py` if it has a diagnostic.
* Use UK spelling in user-facing text.

## Tests

```bash
pip install -e ".[dev,mcp]"
pytest aerosysprop/tests          # full suite
pytest aerosysprop/tests/test_api.py aerosysprop/tests/test_mcp_server.py
```

`test_main.py` and `test_vapour_pressure.py` use relative paths and must be run from `aerosysprop/tests/`. Three of their assertions compare floats exactly and can fail in the last digit on some platforms; this predates the agent API.
