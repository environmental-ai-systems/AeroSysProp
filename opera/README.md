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

## Understanding the output & methodology

OPERA predictions are **QSAR-model estimates, not measurements**. For each
endpoint the output carries reliability indicators you should check before using
a value:

- `*_pred` — the predicted value (units per endpoint; see the in-app endpoint
  descriptions or OPERA's bundled `help.txt`).
- `*_predRange` — a prediction interval.
- `AD_*` — applicability-domain flag (`1` = inside the model's domain, `0` =
  outside → treat with caution).
- `AD_index_*` — continuous similarity-to-training-set measure (0–1).
- `Conf_index_*` — accuracy/confidence estimate (0–1).

For the model methodologies, training data, validation and intended use, consult
the **OPERA documentation** and the primary references (also in `help.txt`):

- Mansouri et al., *SAR & QSAR Env. Res.* (2016) — <https://doi.org/10.1080/1062936X.2016.1253611>
- Mansouri et al., *J. Cheminform.* (2018) — <https://doi.org/10.1186/s13321-018-0263-1>
- Mansouri et al., *J. Cheminform.* (2019) — <https://doi.org/10.1186/s13321-019-0384-1>
- CATMoS acute toxicity, Mansouri et al., *EHP* (2021) — <https://doi.org/10.1289/EHP8495>
- CompTox Chemistry Dashboard — <https://comptox.epa.gov/dashboard>
- Project home — <https://github.com/NIEHS/OPERA>

The app shows a short version of this note on the OPERA tab.

## Calling OPERA from Python

AeroSysProp's own property/AIOMFAC methods are a Python library (see the main
README's *Quick start*). The **external services are reached over HTTP**, so any
language can call them. Point at the running service — native app:
`http://localhost:8601`; inside the compose network: `http://opera:8601`; or a
remote host you manage.

Synchronous (small batches):

```python
import requests
BASE = "http://localhost:8601"
r = requests.post(f"{BASE}/predict", timeout=600, json={
    "smiles": ["CCO", "CC(=O)O"],
    "names":  ["ethanol", "acetic_acid"],
    # "endpoints": ["BP", "VP", "logP", "WS"],   # optional subset (-> OPERA -e)
})
r.raise_for_status()
data = r.json()                        # {"n", "columns", "records"}
print(data["records"][0]["BP_pred"])   # -> 78
```

Asynchronous (recommended for large inputs — won't time out):

```python
import time, requests
BASE = "http://localhost:8601"

jid = requests.post(f"{BASE}/jobs", timeout=30,
                    json={"smiles": smiles_list, "names": names}).json()["id"]

while True:
    s = requests.get(f"{BASE}/jobs/{jid}",
                     params={"results": "false"}, timeout=30).json()
    print(s["status"], s["n_done"], "/", s["n_total"], "-", s.get("message"))
    if s["status"] in ("done", "error", "cancelled"):
        break
    time.sleep(5)

if s["status"] == "done":
    import pandas as pd
    full = requests.get(f"{BASE}/jobs/{jid}", timeout=120).json()
    df = pd.DataFrame(full["records"])
    print(df[["MoleculeID", "BP_pred", "LogVP_pred"]])

# requests.delete(f"{BASE}/jobs/{jid}")   # cancel a running job
```

Notes: add `"format": "csv"` to a `/predict` body to get CSV text (`r.text`)
instead of JSON. A `503` means OPERA isn't installed in the service — check
`/health` (`opera_available`). The same pattern works for any future HTTP
service added to `aerosysprop/services.py`.

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
