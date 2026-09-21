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

* The API and MCP server are **read-only**. Do not add tools that write files or databases, and never pass `output_csv` to `run_all` from `api.py`. Database reads use `read_only=True`.
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
