# syntax=docker/dockerfile:1
#
# AeroSysProp — Streamlit application image
# =========================================
# A small, self-contained image that runs the AeroSysProp web UI (single
# component properties + AIOMFAC mixture activity coefficients).
#
# All runtime dependencies (rdkit, numpy, pandas, streamlit, tqdm) install
# cleanly from PyPI wheels, so there is no need for conda or a GitHub clone —
# the build copies the local source tree straight in.
#
# Build context is the REPOSITORY ROOT (see docker-compose.yml), so the COPY
# paths below are relative to the repo root, not to this Docker/ folder.
#
#   # from the repo root:
#   docker compose up --build app          # -> http://localhost:8501
#   # or without compose:
#   docker build -f aerosysprop/Docker/Dockerfile -t aerosysprop-app .
#   docker run --rm -p 8501:8501 aerosysprop-app
#
# Note on architecture: this image builds natively for your host (arm64 on
# Apple Silicon, amd64 on Intel). The future OPERA service is x86_64-only and
# lives in its own container — the two talk over the docker-compose network, so
# they do NOT need to share a CPU architecture.

FROM python:3.12-slim

LABEL maintainer="david.topping@manchester.ac.uk"
LABEL org.opencontainers.image.title="AeroSysProp app"
LABEL org.opencontainers.image.description="Streamlit UI for AeroSysProp single-component properties and AIOMFAC mixture activity coefficients"
LABEL org.opencontainers.image.source="https://github.com/loftytopping/AeroSysProp"

# RDKit's manylinux wheels dlopen a couple of X11 shared libraries at import
# time even for headless use; install just those.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Streamlit is a long-lived server: keep Python output unbuffered and run it
# headless on all interfaces inside the container.
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

WORKDIR /app

# Copy the packaging metadata and the package itself, then install. Runtime
# dependencies are declared in pyproject.toml's [project.dependencies]; the
# Streamlit interface is the optional `ui` extra, so install with `.[ui]`.
COPY pyproject.toml ./
COPY aerosysprop ./aerosysprop
RUN pip install ".[ui]"

EXPOSE 8501

# Streamlit exposes a readiness endpoint; use it for the container healthcheck.
HEALTHCHECK --interval=30s --timeout=5s --start-period=25s --retries=3 \
    CMD python -c "import sys,urllib.request; sys.exit(0 if urllib.request.urlopen('http://localhost:8501/_stcore/health', timeout=3).read()==b'ok' else 1)"

# The entry script bootstraps sys.path, so running it directly works whether or
# not the package is installed.
CMD ["streamlit", "run", "aerosysprop/app.py"]
