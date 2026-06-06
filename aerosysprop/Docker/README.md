# Running AeroSysProp in Docker

This folder holds the container image for the **AeroSysProp web app** (the
Streamlit UI: single-component properties + AIOMFAC mixture activity
coefficients). The app image is intentionally small and self-contained — all
dependencies install from PyPI wheels, so there is no conda environment and no
GitHub clone at build time.

The OPERA predictive tool is **not** part of this image. The plan is to run
OPERA as a *separate* container and have the app talk to it over the
docker-compose network (see [Phase 2](#phase-2-adding-opera) below).

---

## Prerequisites

- **Docker Desktop** (or Docker Engine + the Compose plugin). On macOS install
  Docker Desktop; nothing else is required on the host.
- Run all commands **from the repository root** (where `docker-compose.yml`
  lives), because the build context is the repo root.

## Phase 1 — run the app (now)

```bash
# from the repository root
docker compose up --build app        # foreground; Ctrl-C to stop
# or in the background:
docker compose up -d --build app
```

Then open <http://localhost:8501>.

**Stop / close the app:**

```bash
docker compose down                  # stop AND remove the container + network
# or stop without removing (faster restart later):
docker compose stop app
docker compose start app             # bring it back up
```

`docker compose down` is the clean shutdown — the built image stays cached for
next time. To remove that too: `docker image rm aerosysprop-app:latest`.

Without compose (equivalent):

```bash
docker build -f aerosysprop/Docker/Dockerfile -t aerosysprop-app .
docker run --rm -p 8501:8501 aerosysprop-app
```

### Notes

- **Architecture.** The image builds natively for your host — `arm64` on Apple
  Silicon, `amd64` on Intel. That is the fast path for the app. OPERA (phase 2)
  is x86_64-only and stays in its own container, so the two never need to share
  an architecture.
- **Health.** The container reports healthy once Streamlit's
  `/_stcore/health` endpoint returns `ok` (≈20–30 s after start).
- **Live code editing.** The image bakes the source in. For a fast dev loop you
  can instead mount the working tree over it:
  ```bash
  docker run --rm -p 8501:8501 -v "$PWD/aerosysprop:/app/aerosysprop" aerosysprop-app
  ```

---

## Phase 2 — adding OPERA

Goal: two containers, networked by `docker-compose.yml` — the **app** in its own
image and **OPERA** in its own image. Nothing about the app image changes.

The single biggest unknown to resolve *first*: **does an OPERA Linux build run
on this Mac at all?** OPERA is a compiled MATLAB-Runtime command-line tool
distributed as x86_64. On Apple Silicon it runs under x86 emulation, which is
usually fine but can be slow or flaky for MATLAB-Runtime apps. So the sequence
is:

1. **Prove it runs.** Build/obtain an `./opera/Dockerfile` (MATLAB Runtime +
   the OPERA command-line tool), then confirm it predicts on a sample `.smi`:
   ```bash
   OPERA2.x_CL -s /exchange/input.smi -o /exchange/opera_results/
   ```
2. **Wire it to the app.** Two clean options:
   - **Shared file drop-box** (simplest): both services mount
     `./.docker-exchange` at `/exchange`. The app writes `input.smi` there and
     reads the results CSV back. (The current "Run locally" button in the
     External Services page already exports the `.smi`; it would point here.)
   - **Thin HTTP wrapper** (cleanest decoupling): wrap the OPERA CLI in a small
     Flask/FastAPI endpoint inside the OPERA container; the app POSTs SMILES and
     gets the CSV back. The app only needs the service URL.

The `opera` service is sketched (commented out) at the bottom of
`docker-compose.yml`, pinned to `platform: linux/amd64`. Uncomment it — and the
`./.docker-exchange` volume on the `app` service — once step 1 succeeds.

> If x86-emulated OPERA proves too painful on this Mac, the same compose setup
> works with OPERA running on a Linux host/CI instead — only the service's
> address changes, not the app.
