# OPERA prediction service (phase 2)

A small HTTP service that wraps the [OPERA](https://github.com/NIEHS/OPERA) QSAR
command-line tool, so the AeroSysProp app can request predictions over the
network and never needs OPERA — or MATLAB Runtime — installed itself.

> **Optional & externally managed.** OPERA is a third-party tool with its own
> licence and release cycle; AeroSysProp never requires it. If you already have
> OPERA installed you don't need this container at all — and on Apple Silicon you
> don't configure Rosetta yourself. For the full "do I already have OPERA? / do I
> need emulation?" decision tree, see
> [`EXTERNAL_TOOLS.md`](EXTERNAL_TOOLS.md).

```
 app (native or container)  ──POST /predict {smiles:[...]}──▶  opera container
                            ◀──────── JSON records ──────────  (OPERA CL + MCR)
```

- The app finds the service via the **`OPERA_URL`** environment variable
  (default `http://localhost:8601`). Native app → `localhost:8601`; app in the
  compose network → `http://opera:8601`.
- OPERA itself emits CSV/TSV; **the wrapper parses it and returns JSON** (or CSV
  passthrough with `format: "csv"`). See [`server.py`](server.py).

## HTTP API

| Endpoint | Purpose |
|---|---|
| `GET /health` | Liveness + readiness (`opera_available`). |
| `POST /predict` | **Synchronous** run — convenient for small batches / curl. Body `{smiles, names?, endpoints?, format?}`; returns `{n, columns, records}` (or CSV). |
| `POST /jobs` | **Asynchronous** run for large inputs — returns a `{id, status:"queued", ...}` immediately and runs OPERA in a background thread. Same body as `/predict`. |
| `GET /jobs/{id}` | Poll status/progress: `status` (queued/running/done/error/cancelled), `n_done`/`n_total` (descriptors processed), `elapsed`, `message` (latest OPERA log line); `records`+`columns` when done. Use `?results=false` to omit records while polling. |
| `DELETE /jobs/{id}` | Request cancellation of a running job (kills the OPERA process group). |

The app uses `/jobs` for *Run OPERA* so large files (thousands of compounds) get
a live progress bar and Cancel instead of one long request that would time out.
Endpoint selection (`endpoints: ["logP","BP",...]`) maps to OPERA's `-e`; omit
for all endpoints.

## Files

| File | Purpose |
|---|---|
| `server.py` | FastAPI wrapper: `/health`, `/predict`, and the async `/jobs` API. Runs OPERA via subprocess, parses output, returns JSON. |
| `requirements.txt` | Python deps for the wrapper (fastapi, uvicorn, pandas, pydantic). |
| `Dockerfile` | Image: Ubuntu + OPERA CL (MATLAB Runtime bundled) + the wrapper. x86_64. Verified against OPERA v2.9.2. |

## Status: verified end-to-end

This has been **built and run successfully**, including under x86 emulation
(Rosetta) on an Apple Silicon Mac:

- OPERA **v2.9.2** Linux CL installs and runs; the `_mcr` release **bundles
  MATLAB Runtime v912 (R2022a)** — there is no separate MathWorks download.
- Install lands at OPERA's default `/usr/local/bin/OPERA`; launcher
  `application/run_OPERA.sh <mcr_dir>` (the Dockerfile installs to the default
  path to avoid OPERA's custom-folder `installdir` check).
- The wrapper drives OPERA and returns JSON: a 2-molecule `/predict` returned
  150 columns (all endpoints), e.g. ethanol BP 78 °C, acetic acid BP 118 °C.
- Rough timings under Rosetta: MCR install ≈ 1.5 min (build time); a cold
  `/predict` ≈ 60–75 s (MCR start + PaDEL descriptors dominate; OPERA CL spins up
  a fresh process per request).

> The single thing you supply is the OPERA release tarball (it is ~2 GB and
> NIEHS-distributed, so it is not committed here) — see the steps below.

## Build and run

1. **Download the OPERA Linux CL release** and drop it into `opera/vendor/`
   (gitignored):
   ```
   opera/vendor/OPERA2.9_CL_mcr.tar.xz
   ```
   from <https://github.com/NIEHS/OPERA/releases> (~1.96 GB; bundles MCR v912).
2. **Build** (x86_64; emulated via Rosetta on Apple Silicon — automatic):
   ```bash
   docker compose --profile opera build opera
   ```
3. **Run and check readiness:**
   ```bash
   docker compose --profile opera up -d opera
   curl -s localhost:8601/health      # -> {"opera_available": true, ...}
   ```
4. **Smoke-test a prediction** (first call is a cold start, ~60–75 s):
   ```bash
   curl -s -X POST localhost:8601/predict \
        -H 'Content-Type: application/json' \
        -d '{"smiles":["CCO","CC(=O)O"],"names":["ethanol","acetic_acid"]}' | head -c 400
   ```
5. **Use it from the app.** Start the app (native or container) and open
   **External Services → OPERA**. The tab auto-detects the service; when it is
   up you get a ✅ and the *Run OPERA* button is enabled.

If emulated OPERA is too slow for your needs, the same setup works with the
`opera` service running on a Linux host/CI — only `OPERA_URL` changes, not the
app or the wrapper.

### Adapting to other OPERA versions / endpoint selection
The wrapper assumes `run_OPERA.sh <mcr> -s <input.smi> -o <output.csv>`, which
computes **all** endpoints by default. To pin a different version, change the
tarball name + `OPERA_CMD`/`MCR_ROOT` in the Dockerfile. To narrow endpoints,
pass `endpoints` in the `/predict` body (mapped to OPERA's `-e`), or set
`OPERA_DEFAULT_ARGS`. OPERA's endpoint codes (from `help.txt`): `logP, BP, MP,
logVP, WS, HL, AOH, BCF, BioDeg, RB, KM, KOA, logKoc, RT, pKa, logD, ER, AR,
CATMoS, FuB, Clint, Caco2` and groups `PC, EF, Tox, ADME, StrP`.

## Local wrapper test (no OPERA needed)

You can exercise the wrapper's HTTP layer without OPERA installed — `/health`
will simply report `opera_available: false` and `/predict` returns 503:

```bash
pip install -r opera/requirements.txt
uvicorn server:app --app-dir opera --port 8601
curl -s localhost:8601/health
```
