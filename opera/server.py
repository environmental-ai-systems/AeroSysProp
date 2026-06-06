"""
opera.server
============
A thin HTTP wrapper around the OPERA command-line tool, so the AeroSysProp app
can run OPERA predictions over HTTP — without the app needing OPERA (or MATLAB
Runtime) installed.

The app stays decoupled: it only needs this service's URL (``OPERA_URL``). The
service shells out to OPERA internally, parses OPERA's CSV output, and returns
JSON (or CSV passthrough on request).

Endpoints
---------
``GET    /health``        liveness + readiness (is the OPERA command present?)
``GET    /``              basic service info
``POST   /predict``       run OPERA synchronously; returns JSON (small batches)
``POST   /jobs``          submit an asynchronous job; returns a job id immediately
``GET    /jobs/{id}``     poll a job's status/progress (``?results=false`` to omit
                          the records while polling); returns records when done
``DELETE /jobs/{id}``     request cancellation of a running job

Why async jobs
--------------
OPERA is a batch tool with a heavy MATLAB-Runtime cold start; large inputs
(thousands of molecules) can run for many minutes. A synchronous request would
hold a connection open the whole time and time out. The job API runs OPERA in a
background thread and reports progress (descriptors calculated / total, plus the
latest OPERA log line), so the UI can show a live status bar and offer cancel.

Configuration (environment variables)
--------------------------------------
``OPERA_CMD``           command that launches OPERA, e.g.
                        ``"/usr/local/bin/OPERA/application/run_OPERA.sh /usr/local/bin/OPERA/v912"``.
``OPERA_DEFAULT_ARGS``  extra OPERA flags always appended (shlex-split).
``OPERA_TIMEOUT``       per-job timeout in seconds (default 7200 = 2 h).
``OPERA_JOB_RETENTION`` seconds to keep finished jobs in memory (default 3600).

NOTE: single-worker deployment assumed — the job registry is in-process memory.
(The bundled Dockerfile runs one uvicorn worker.)
"""
from __future__ import annotations

import glob
import json
import os
import shlex
import shutil
import signal
import subprocess
import tempfile
import threading
import time
import uuid
from typing import Dict, List, Optional

import pandas as pd
from fastapi import FastAPI, HTTPException
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

OPERA_CMD = os.environ.get("OPERA_CMD", "OPERA_CL")
OPERA_DEFAULT_ARGS = os.environ.get("OPERA_DEFAULT_ARGS", "")
OPERA_TIMEOUT = int(os.environ.get("OPERA_TIMEOUT", "7200"))
JOB_RETENTION_S = int(os.environ.get("OPERA_JOB_RETENTION", "3600"))

app = FastAPI(
    title="OPERA prediction service",
    summary="HTTP wrapper around the OPERA QSAR command-line tool.",
    version="0.2.0",
)


# ---------------------------------------------------------------------------
# Request model
# ---------------------------------------------------------------------------

class PredictRequest(BaseModel):
    smiles: List[str] = Field(..., description="SMILES strings to predict.")
    names: Optional[List[str]] = Field(
        None,
        description="Optional identifiers aligned with `smiles`. Defaults to "
                    "mol_1, mol_2, ... when omitted.",
    )
    endpoints: Optional[List[str]] = Field(
        None,
        description="Optional OPERA endpoint selection (passed via `-e`); "
                    "omit for all endpoints.",
    )
    format: str = Field(
        "json",
        description="Response format for /predict: 'json' or 'csv'.",
    )


# ---------------------------------------------------------------------------
# OPERA invocation helpers
# ---------------------------------------------------------------------------

def _opera_executable() -> str:
    parts = shlex.split(OPERA_CMD)
    return parts[0] if parts else OPERA_CMD


def _opera_available() -> bool:
    exe = _opera_executable()
    if os.path.sep in exe:
        return os.path.isfile(exe) and os.access(exe, os.X_OK)
    from shutil import which
    return which(exe) is not None


def _write_smi(smiles: List[str], names: Optional[List[str]], path: str) -> None:
    """Write a tab-separated `SMILES<TAB>ID` per line input file for OPERA."""
    if names is None:
        names = [f"mol_{i+1}" for i in range(len(smiles))]
    if len(names) != len(smiles):
        raise ValueError("`names` must align with `smiles`.")
    with open(path, "w") as fh:
        for smi, name in zip(smiles, names):
            fh.write(f"{smi}\t{name}\n")


def _build_argv(smi_path: str, out_path: str, endpoints: Optional[List[str]]) -> List[str]:
    argv = shlex.split(OPERA_CMD)
    argv += ["-s", smi_path, "-o", out_path]
    if OPERA_DEFAULT_ARGS:
        argv += shlex.split(OPERA_DEFAULT_ARGS)
    if endpoints:
        argv += ["-e", *endpoints]
    return argv


def _read_opera_output(out_path: str, work_dir: str) -> pd.DataFrame:
    """Locate and parse OPERA's output table (explicit -o path, else any
    tabular file produced in the working directory)."""
    candidates = []
    if os.path.isfile(out_path):
        candidates.append(out_path)
    else:
        for ext in ("*.csv", "*.tsv", "*.txt"):
            candidates += glob.glob(os.path.join(work_dir, ext))
    if not candidates:
        raise FileNotFoundError(
            "OPERA produced no .csv/.tsv/.txt output in %s" % work_dir
        )
    target = max(candidates, key=os.path.getmtime)
    sep = "," if target.lower().endswith(".csv") else "\t"
    return pd.read_csv(target, sep=sep)


def _records_from_df(df: pd.DataFrame):
    """JSON-safe records: pandas' to_json emits null for NaN and converts numpy
    scalars to native Python (Starlette's JSON serializer rejects NaN/Inf)."""
    return json.loads(df.to_json(orient="records"))


def _count_descriptor_lines(work_dir: str) -> int:
    """Best-effort count of molecules whose descriptors OPERA has written so far
    (the descriptor CSV grows as PaDEL/CDK process molecules). Header excluded."""
    best = 0
    for pat in ("*_PadelDesc.csv", "*_CDKDesc.csv"):
        for f in glob.glob(os.path.join(work_dir, pat)):
            try:
                with open(f, "rb") as fh:
                    n = sum(1 for _ in fh)
                best = max(best, max(0, n - 1))
            except OSError:
                pass
    return best


def _tail_line(path: str) -> str:
    try:
        with open(path, "r", errors="replace") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        return lines[-1] if lines else ""
    except OSError:
        return ""


def _tail(path: str, n: int = 4000) -> str:
    try:
        with open(path, "r", errors="replace") as fh:
            return fh.read()[-n:]
    except OSError:
        return ""


def _kill_process_group(proc: "subprocess.Popen") -> None:
    """Terminate OPERA and its MATLAB-Runtime children (started in a new
    session, so we can signal the whole process group)."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(os.getpgid(proc.pid), sig)
        except Exception:
            try:
                proc.terminate() if sig == signal.SIGTERM else proc.kill()
            except Exception:
                pass
        try:
            proc.wait(timeout=10)
            return
        except Exception:
            continue


# ---------------------------------------------------------------------------
# Synchronous run (used by /predict for small batches)
# ---------------------------------------------------------------------------

def _run_opera(req: PredictRequest):
    if not req.smiles:
        raise HTTPException(status_code=422, detail="`smiles` is empty.")
    if not _opera_available():
        raise HTTPException(
            status_code=503,
            detail="OPERA command not found: %r. Set OPERA_CMD and ensure OPERA "
                   "+ MATLAB Runtime are installed." % _opera_executable(),
        )
    with tempfile.TemporaryDirectory(prefix="opera_") as work:
        smi_path = os.path.join(work, "input.smi")
        out_path = os.path.join(work, "opera_results.csv")
        try:
            _write_smi(req.smiles, req.names, smi_path)
        except ValueError as e:
            raise HTTPException(status_code=422, detail=str(e))
        argv = _build_argv(smi_path, out_path, req.endpoints)
        try:
            proc = subprocess.run(argv, cwd=work, capture_output=True,
                                  text=True, timeout=OPERA_TIMEOUT)
        except subprocess.TimeoutExpired:
            raise HTTPException(status_code=504,
                                detail="OPERA timed out after %ds." % OPERA_TIMEOUT)
        if proc.returncode != 0:
            raise HTTPException(status_code=500,
                                detail="OPERA exited %d. stderr:\n%s"
                                       % (proc.returncode, proc.stderr[-4000:]))
        try:
            df = _read_opera_output(out_path, work)
        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            raise HTTPException(status_code=500,
                                detail="Could not read OPERA output: %s" % e)
        return df


# ---------------------------------------------------------------------------
# Asynchronous job registry
# ---------------------------------------------------------------------------

_JOBS: Dict[str, dict] = {}
_LOCK = threading.Lock()


def _update(job_id: str, **kw) -> None:
    with _LOCK:
        j = _JOBS.get(job_id)
        if j is not None:
            j.update(kw)


def _prune_jobs() -> None:
    now = time.time()
    with _LOCK:
        stale = [jid for jid, j in _JOBS.items()
                 if j["status"] in ("done", "error", "cancelled")
                 and now - j.get("finished", now) > JOB_RETENTION_S]
        for jid in stale:
            _JOBS.pop(jid, None)


def _public_job(job_id: str, include_results: bool = True) -> dict:
    with _LOCK:
        j = _JOBS.get(job_id)
        if j is None:
            raise HTTPException(status_code=404, detail="No such job: %s" % job_id)
        out = {k: v for k, v in j.items() if not k.startswith("_")}
    if not include_results:
        out = {k: v for k, v in out.items() if k not in ("records", "columns")}
    return out


def _run_job(job_id: str, req: PredictRequest) -> None:
    work = tempfile.mkdtemp(prefix="opera_job_")
    log_path = os.path.join(work, "opera_run.log")
    out_path = os.path.join(work, "opera_results.csv")
    smi_path = os.path.join(work, "input.smi")
    n_total = len(req.smiles)
    try:
        _write_smi(req.smiles, req.names, smi_path)
        argv = _build_argv(smi_path, out_path, req.endpoints)
        _update(job_id, status="running", started=time.time(),
                message="Starting OPERA…")
        with open(log_path, "w") as logf:
            proc = subprocess.Popen(argv, cwd=work, stdout=logf,
                                    stderr=subprocess.STDOUT,
                                    start_new_session=True)
        with _LOCK:
            if job_id in _JOBS:
                _JOBS[job_id]["_proc"] = proc
        start = time.time()
        while True:
            rc = proc.poll()
            n_done = min(_count_descriptor_lines(work), n_total)
            _update(job_id, n_done=n_done,
                    message=_tail_line(log_path) or "Running…",
                    elapsed=round(time.time() - start, 1))
            with _LOCK:
                cancel = bool(_JOBS.get(job_id, {}).get("_cancel"))
            if cancel:
                _kill_process_group(proc)
                _update(job_id, status="cancelled", message="Cancelled by user.",
                        finished=time.time())
                return
            if rc is not None:
                break
            if time.time() - start > OPERA_TIMEOUT:
                _kill_process_group(proc)
                _update(job_id, status="error", message="Timed out.",
                        error="OPERA timed out after %ds." % OPERA_TIMEOUT,
                        finished=time.time())
                return
            time.sleep(2)
        if proc.returncode != 0:
            _update(job_id, status="error", message="OPERA failed.",
                    error="OPERA exited %d.\n%s" % (proc.returncode,
                                                    _tail(log_path)),
                    finished=time.time())
            return
        try:
            df = _read_opera_output(out_path, work)
        except Exception as e:
            _update(job_id, status="error", message="Output unreadable.",
                    error="Could not read OPERA output: %s\n%s"
                          % (e, _tail(log_path)),
                    finished=time.time())
            return
        _update(job_id, status="done", n_done=int(len(df)), n=int(len(df)),
                columns=list(df.columns), records=_records_from_df(df),
                message="Done — %d molecules predicted." % len(df),
                finished=time.time())
    except Exception as e:  # pragma: no cover - defensive
        _update(job_id, status="error", error=str(e), message="Failed.",
                finished=time.time())
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.get("/")
def root():
    return {
        "service": "opera",
        "opera_cmd": OPERA_CMD,
        "opera_available": _opera_available(),
        "endpoints": ["/health", "/predict", "/jobs"],
    }


@app.get("/health")
def health():
    return {
        "status": "ok",
        "opera_available": _opera_available(),
        "opera_cmd": OPERA_CMD,
    }


@app.post("/predict")
def predict(req: PredictRequest):
    """Synchronous prediction — convenient for small batches and curl tests.
    For large inputs use the /jobs API so the request doesn't time out."""
    df = _run_opera(req)
    if req.format.lower() == "csv":
        return PlainTextResponse(df.to_csv(index=False), media_type="text/csv")
    return {"n": int(len(df)), "columns": list(df.columns),
            "records": _records_from_df(df)}


@app.post("/jobs")
def submit_job(req: PredictRequest):
    """Submit an asynchronous OPERA job; returns a job id immediately."""
    if not req.smiles:
        raise HTTPException(status_code=422, detail="`smiles` is empty.")
    if not _opera_available():
        raise HTTPException(
            status_code=503,
            detail="OPERA command not found: %r." % _opera_executable(),
        )
    _prune_jobs()
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        _JOBS[job_id] = {
            "id": job_id, "status": "queued", "message": "Queued.",
            "n_total": len(req.smiles), "n_done": 0, "elapsed": 0.0,
            "created": time.time(), "records": None, "columns": None,
            "error": None,
        }
    threading.Thread(target=_run_job, args=(job_id, req), daemon=True).start()
    return _public_job(job_id, include_results=False)


@app.get("/jobs/{job_id}")
def get_job(job_id: str, results: bool = True):
    """Poll a job. Pass ``?results=false`` to omit the (potentially large)
    records/columns while polling; fetch once more with results=true when done."""
    return _public_job(job_id, include_results=results)


@app.delete("/jobs/{job_id}")
def cancel_job(job_id: str):
    """Request cancellation of a queued/running job."""
    with _LOCK:
        j = _JOBS.get(job_id)
        if j is None:
            raise HTTPException(status_code=404, detail="No such job: %s" % job_id)
        if j["status"] in ("queued", "running"):
            j["_cancel"] = True
            j["message"] = "Cancelling…"
    return _public_job(job_id, include_results=False)
