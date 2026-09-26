# Image for running WattEG on a cloud executor (the `gcb` profile).
#
# - The environment is baked in: every task runs `pixi run --frozen ...`, and solving per task would
#   repeat minutes of work thousands of times.
# - /opt/pipeline is load-bearing: conf/gcb.config symlinks each task's ${projectDir} to it. Change
#   it here and change it there.
# - The pipeline is cloned at WATTEG_REF. ALWAYS pass a commit SHA; `main` makes the image silently
#   non-reproducible. Rebuild whenever pixi.lock changes, including the pysceptre pin.
# - pysceptre's scripts/ (sceptre_io.py, which watteg-prepare-sim-input and make_test_data import)
#   is not shipped in its wheel, so it is fetched at the commit pinned in pixi.lock and put on
#   PYTHONPATH. Without it every run fails at PREPARE_SIM_INPUT with "No module named sceptre_io".
# - Keep this file under ~7 KB: Wave rejects a containerfile over 10,000 bytes base64-encoded.
#
#   wave --containerfile containers/Dockerfile --build-arg WATTEG_REF=<sha> --await

FROM ubuntu:22.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl bzip2 git \
    && rm -rf /var/lib/apt/lists/*
ARG PIXI_VERSION=0.76.2
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION="v${PIXI_VERSION}" PIXI_NO_PATH_UPDATE=1 bash
ENV PATH="/root/.pixi/bin:${PATH}"
ARG WATTEG_REF=main
RUN git clone https://github.com/EngreitzLab/WattEG.git /opt/pipeline \
    && cd /opt/pipeline \
    && git checkout "${WATTEG_REF}" \
    && git rev-parse HEAD > /opt/pipeline/.image-built-from
WORKDIR /opt/pipeline
RUN pixi install --frozen
RUN rev=$(grep -o 'pysceptre.git?tag=[^#]*#[0-9a-f]*' pixi.lock | head -1 | sed 's/.*#//') \
    && test -n "${rev}" \
    && git clone https://github.com/broadinstitute/pysceptre.git /opt/pysceptre \
    && git -C /opt/pysceptre checkout "${rev}"
RUN rm -rf /root/.cache/rattler /root/.cache/pixi /root/.cache/uv /root/.pixi/cache

FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /root/.pixi/bin/pixi /root/.pixi/bin/pixi
COPY --from=builder /opt/pipeline /opt/pipeline
COPY --from=builder /opt/pysceptre/scripts /opt/pysceptre-scripts
ENV PATH="/root/.pixi/bin:${PATH}"
ENV PYTHONPATH=/opt/pysceptre-scripts
WORKDIR /opt/pipeline
RUN pixi run --frozen python -c "import watteg, pysceptre, sceptre_io; from pysceptre.pipeline.discovery import run_discovery_ntcells_complement; print('watteg', watteg.__version__)" \
    && pixi run --frozen watteg-prepare-sim-input --help > /dev/null \
    && echo "built from $(cat /opt/pipeline/.image-built-from)"
ENV WATTEG_ROOT=/opt/pipeline
CMD ["/bin/bash"]
