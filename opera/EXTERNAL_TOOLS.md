# External tools (e.g. OPERA) — optional & externally managed

**AeroSysProp is fully functional on its own.** Third-party predictive tools such
as [OPERA](https://github.com/NIEHS/OPERA) are *optional* integrations: they are
**externally-managed software** with their own licences, installers, and release
cycles. AeroSysProp can hand them SMILES and display their results alongside its
own predictions for comparison — but it does **not** bundle, own, vendor, or
require them. You can ignore this entire area and lose no AeroSysProp capability.

This separation is deliberate. It means:

- The app, the property methods, and AIOMFAC never depend on any external tool
  being present, installed, or reachable.
- Each external tool is yours to install, license, update, and run wherever you
  like — on this machine, another machine, a server, or CI.
- Adding/removing an integration is a registry entry (`aerosysprop/services.py`),
  not a change to the core.

---

## Three ways to use OPERA — pick one

### 0. Don't use it (the default)
Nothing to install or configure. Everything in AeroSysProp works without OPERA.

### 1. You already have OPERA (you manage the install)
AeroSysProp just feeds it SMILES — it doesn't care where or how OPERA is installed.

- **File hand-off (zero setup):** *External Services → OPERA → download the
  `.smi`*, run it through your own OPERA, read the results back:
  ```bash
  ./run_OPERA.sh <mcr_dir> -s aerosysprop_for_opera.smi -o predictions.csv
  ```
- **Live hand-off (one-click in the app):** put OPERA behind the thin HTTP wrapper
  ([`opera/server.py`](server.py)) on a machine you control and set the
  **`OPERA_URL`** environment variable to it. The app's *Run OPERA* button then
  calls it directly. The wrapper is only a shim around *your* existing
  `run_OPERA.sh`; it duplicates nothing.

### 2. You don't have OPERA and want it containerised
Our [`opera/Dockerfile`](Dockerfile) builds a self-contained OPERA service from
the official NIEHS release — which **bundles MATLAB Runtime v912**, so there is no
MathWorks account or separate download. See [`README.md`](README.md):
```bash
docker compose --profile opera up --build opera   # app then reaches it at :8601
```

> In all cases the app stays decoupled: it only needs a URL (`OPERA_URL`, default
> `http://localhost:8601`) or the downloadable `.smi`. Run the app natively or in
> its own container — either works.

---

## The x86 / Rosetta question (Apple Silicon only)

OPERA is compiled for **x86_64** and ships with the **MATLAB Runtime**. Apple
Silicon Macs are **arm64**, so OPERA cannot run *natively* on an M-series Mac.
This only matters if you intend to run OPERA itself on such a Mac.

| Your situation | What you do | Rosetta / emulation? |
|---|---|---|
| **No OPERA** (fine without it) | nothing | n/a |
| **OPERA already running somewhere x86** (Intel Mac, Linux box/server, VM, CI) | file hand-off, or set `OPERA_URL` to your HTTP wrapper | **No** — your arm64 app just talks to it over the network; architectures need not match |
| **Want OPERA in a container on an Apple-Silicon Mac** | `docker compose --profile opera up --build opera` | **Yes, automatic** — see below |
| **Want OPERA in a container on Intel mac / Linux** | same command | **No** — image is native amd64 |

### If you already have OPERA installed
You do **not** need Rosetta, our Dockerfile, or any emulation on your side. OPERA
runs wherever you installed it (it must be an x86 OS — Intel Mac, Linux, Windows,
a VM/server). Expose it via the file hand-off or the HTTP wrapper and point
`OPERA_URL` at it. The arm64 app communicates with the x86 OPERA over HTTP; the
mismatch is irrelevant.

### If you want our container to handle it (Apple Silicon)
You don't configure anything. The OPERA image is pinned to
`platform: linux/amd64` in `docker-compose.yml`, so Docker Desktop runs it under
x86 emulation. Recent Docker Desktop uses **Rosetta** for this
(*Settings → General → “Use Rosetta for x86/amd64 emulation”*, on by default),
which is far faster than QEMU and gives the MATLAB Runtime a realistic chance of
running well. The `platform:` pin plus Docker Desktop's Rosetta do all the work.

Check emulation is available:
```bash
docker run --rm --platform linux/amd64 ubuntu uname -m   # expect: x86_64
```
If that prints `x86_64`, you're ready — `docker compose --profile opera up --build
opera` will build and run OPERA emulated, and the app talks to it normally.

### Intel Mac or Linux host
The OPERA image is native `amd64` — **no emulation, no Rosetta** — and the exact
same commands apply.

---

## Container lifecycle (starting & stopping)

The app and OPERA are **separate containers with independent lifecycles** — this
is deliberate, so a natively-run app can still use OPERA. The practical
consequence to be aware of:

- **Closing the app does *not* stop OPERA** when the app runs natively
  (`streamlit run …`). The OPERA container keeps running until you stop it:
  ```bash
  docker compose --profile opera down     # stop & remove the OPERA service
  # …or stop it without removing (faster restart later):
  docker compose stop opera
  ```
- **Running both via compose** ties their lifecycles together — a single
  `docker compose down` stops both:
  ```bash
  docker compose --profile opera up        # app + OPERA together
  docker compose down                       # stops both
  ```
- The OPERA container is idle-cheap (it only uses CPU during a prediction) but
  holds a ~4.5 GB image and its memory, so stop it when you're done.

The app reinforces this: when it detects the OPERA service, the
**External Services → OPERA** tab shows the stop command inline.

## Large inputs

OPERA is a batch tool with a heavy MATLAB-Runtime cold start, so big files
(thousands of compounds) can run for many minutes. The service exposes an
**asynchronous job API** (`POST /jobs`, `GET /jobs/{id}`, `DELETE /jobs/{id}`)
and the app uses it: clicking *Run OPERA* submits a background job and shows a
**live progress bar** (descriptors calculated / total, elapsed, latest OPERA log
line) with a **Cancel** button — so a 10,000-compound run won't time out or
freeze the UI. See `opera/README.md` for the API.

## Bottom line
- OPERA (and any external tool) is **opt-in**; the core never needs it.
- **Already have it?** Point the app at it (`OPERA_URL`) or use the `.smi`
  export — no Rosetta, no container required.
- **Don't have it?** Our `--profile opera` container builds it (MCR bundled) and,
  on Apple Silicon, Docker + Rosetta handle the x86 emulation automatically.
