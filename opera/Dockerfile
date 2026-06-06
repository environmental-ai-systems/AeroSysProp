# syntax=docker/dockerfile:1
#
# OPERA prediction service (x86_64)
# =================================
# Wraps the OPERA command-line QSAR tool in a small HTTP service (opera/server.py)
# so the AeroSysProp app can request predictions over the network and never needs
# OPERA or MATLAB Runtime itself.
#
# The build steps below are VERIFIED against OPERA v2.9.2 (Linux CL, `_mcr` build,
# which bundles MATLAB Runtime v912 / R2022a — no MathWorks download needed).
#
# You supply the release tarball (it is large + NIEHS-distributed, so it is not in
# this repo). Download it and drop it in opera/vendor/ before building:
#
#   opera/vendor/OPERA2.9_CL_mcr.tar.xz
#   from https://github.com/NIEHS/OPERA/releases  (~1.96 GB)
#
# Then, from the repository root:
#   docker compose --profile opera build opera
#   docker compose --profile opera up opera
#   curl localhost:8601/health           # -> {"opera_available": true, ...}
#
# OPERA is x86_64 only; on Apple Silicon this image runs under emulation (Rosetta)
# via the `platform: linux/amd64` pin in docker-compose.yml. See
# opera/EXTERNAL_TOOLS.md for the full emulation story.

FROM --platform=linux/amd64 ubuntu:22.04

LABEL org.opencontainers.image.title="AeroSysProp OPERA service"
LABEL org.opencontainers.image.description="HTTP wrapper around the OPERA QSAR CLI (v2.9.2, MCR v912 bundled)"

ENV DEBIAN_FRONTEND=noninteractive

# Build/runtime dependencies:
#  - unzip + xz-utils: the OPERA installer unpacks its bundled MCR with unzip
#    (without it the installer wrongly falls back to a network download and fails)
#  - default-jre: OPERA computes PaDEL/CDK molecular descriptors in Java
#  - libX*/libgtk/etc: MATLAB Runtime (R2022a) shared-library dependencies
#  - python3/pip: the FastAPI wrapper
RUN apt-get update && apt-get install -y --no-install-recommends \
        unzip xz-utils ca-certificates default-jre \
        libxt6 libx11-6 libxext6 libxmu6 libxpm4 libxtst6 libxcomposite1 \
        libxrandr2 libxcursor1 libxdamage1 libxfixes3 libxft2 libxinerama1 \
        libxrender1 libgl1 libglu1-mesa libasound2 libcairo2 libpango-1.0-0 \
        libgtk2.0-0 libnss3 libxss1 \
        python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install OPERA + its bundled MATLAB Runtime. We install to OPERA's *default*
# location (/usr/local/bin/OPERA): OPERA records its expected install dir and
# refuses to run from a non-default path, so installing to the default avoids any
# runtime patching. The tarball is bind-mounted for this step only (BuildKit), so
# the ~2 GB archive never becomes an image layer; only the installed tree is kept.
RUN --mount=type=bind,source=vendor/OPERA2.9_CL_mcr.tar.xz,target=/tmp/opera.tar.xz \
    set -eux; \
    mkdir -p /tmp/ox; \
    tar xJf /tmp/opera.tar.xz -C /tmp/ox; \
    INST="$(find /tmp/ox -name '*Installer.install' | head -1)"; \
    cd "$(dirname "$INST")"; \
    chmod +x ./*.install; \
    ./*.install -mode silent -agreeToLicense yes -destinationFolder /usr/local/bin/OPERA; \
    rm -rf /tmp/ox

# Verified install layout: launcher + MATLAB Runtime root.
ENV MCR_ROOT=/usr/local/bin/OPERA/v912
ENV OPERA_CMD="/usr/local/bin/OPERA/application/run_OPERA.sh /usr/local/bin/OPERA/v912"
# Extra flags always passed to OPERA. Empty => all endpoints (OPERA's default).
ENV OPERA_DEFAULT_ARGS=""

# The HTTP wrapper.
WORKDIR /srv
COPY requirements.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt
COPY server.py ./

ENV PYTHONUNBUFFERED=1
EXPOSE 8601

HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD python3 -c "import sys,urllib.request; sys.exit(0 if urllib.request.urlopen('http://localhost:8601/health', timeout=3).status==200 else 1)"

CMD ["uvicorn", "server:app", "--host", "0.0.0.0", "--port", "8601"]
