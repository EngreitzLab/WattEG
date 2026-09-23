# Image for running WattEG on a cloud executor (the `gcb` profile).
#
# /opt/pipeline is load-bearing. conf/gcb.config's beforeScript symlinks each task's ${projectDir}
# here; without it `pixi run` finds a manifest with no environment beside it and every task dies
# with "pixi: command not found". Change it here and change it there.
#
# The environment is baked in because every task shells out to `pixi run --frozen`, and because it
# cannot be rebuilt from a conda spec alone: sceptre and ondisc install from pinned commits and
# compile C++. `pixi run setup` is what does that.
#
# BUILD:
#   wave --containerfile containers/Dockerfile --await
#   docker build --build-arg WATTEG_REF=<sha> -t <registry>/watteg:<tag> -f containers/Dockerfile .
#
# ALWAYS PASS A COMMIT SHA. The default is `main`, which makes the image non-reproducible.
# REBUILD WHENEVER pixi.lock CHANGES -- nothing enforces it, and a stale image fails at the first
# task needing the new package with an R error about a missing library.
#
# Keep this file under ~7 KB: Wave rejects a containerfile over 10,000 bytes base64-encoded.
#
# Two stages, not one: a `RUN rm -rf` in a later layer records the deletion but still ships the
# bytes. Measured -- single stage 2.14 GiB, two stages 0.68 GiB.

FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# ca-certificates: stock ubuntu:22.04 ships none, so every HTTPS call fails.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl bzip2 git \
    && rm -rf /var/lib/apt/lists/*

# Pinned: an older pixi cannot read a newer lockfile format.
ARG PIXI_VERSION=0.76.2
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION="v${PIXI_VERSION}" PIXI_NO_PATH_UPDATE=1 bash
ENV PATH="/root/.pixi/bin:${PATH}"

ARG WATTEG_REF=main
RUN git clone https://github.com/EngreitzLab/WattEG.git /opt/pipeline \
    && cd /opt/pipeline \
    && git checkout "${WATTEG_REF}" \
    && git rev-parse HEAD > /opt/pipeline/.image-built-from

WORKDIR /opt/pipeline

# Mirrors workflow/slurm_executor/00_setup_env.sbatch: the two must not drift, or a cloud run and a
# Sherlock run stop being comparable.
RUN pixi install --frozen \
    && pixi install --frozen --environment build \
    && pixi run --frozen setup \
    && pixi run --frozen check-api

RUN rm -rf /opt/pipeline/.pixi/envs/build /root/.cache/rattler /root/.cache/pixi /root/.pixi/cache

# ---- runtime -------------------------------------------------------------------------------
# Same base and the SAME PATHS: conda environments embed absolute paths.

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# pixi comes too: process scripts invoke `pixi run` themselves.
COPY --from=builder /root/.pixi/bin/pixi /root/.pixi/bin/pixi
COPY --from=builder /opt/pipeline /opt/pipeline

ENV PATH="/root/.pixi/bin:${PATH}"
WORKDIR /opt/pipeline

# Fail the BUILD rather than the first task of a 5,000-task sweep.
RUN pixi run --frozen check-api \
    && echo "built from $(cat /opt/pipeline/.image-built-from)"

ENV WATTEG_ROOT=/opt/pipeline

CMD ["/bin/bash"]
