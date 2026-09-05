# Image for running WattEG on a cloud executor (the `gcb` profile).
#
# WHY THE ENVIRONMENT IS BAKED IN RATHER THAN SOLVED AT RUN TIME
#
# Every task runs `pixi run --frozen --manifest-path ${projectDir}/pixi.toml Rscript ...`, so the
# solved environment has to exist on disk before the first task starts. Solving it per task would
# be minutes of work repeated thousands of times, and would need outbound network from every worker.
#
# It also CANNOT be reproduced from a conda specification alone, which rules out Wave's
# conda-recipe mode: sceptre is installed from a pinned commit with a local patch applied, ondisc
# from another pinned commit, and both compile C++. `pixi run setup` is what does that, so the image
# has to run it at build time. See src/install_sceptre.R and src/install_ondisc.R.
#
# WHY /opt/pipeline
#
# conf/gcb.config's beforeScript symlinks the task's ${projectDir} to /opt/pipeline. Google Batch
# stages pipeline source into a per-task clone path that has nothing to do with where this image's
# environment lives, so without that symlink `pixi run` finds a manifest with no environment beside
# it and every task dies with "pixi: command not found". The path is therefore load-bearing: change
# it here and change it there.
#
# BUILDING
#
# The pipeline is CLONED at a pinned ref rather than COPYed from a build context. Two reasons:
#
#   * Wave can then build it from the Dockerfile alone, with no context upload. The repository is
#     public, so no credentials are involved.
#   * The image records exactly which commit it contains. A context copy records nothing, and this
#     image has to be rebuilt whenever pixi.lock changes -- being able to ask an image what it is
#     matters more than saving a clone.
#
#   wave --containerfile containers/Dockerfile --await
#   docker build --build-arg WATTEG_REF=<sha> -t <registry>/watteg:<tag> -f containers/Dockerfile .
#
# ALWAYS PASS A COMMIT SHA for anything but a throwaway build. The default is `main`, which makes
# the image silently non-reproducible -- two builds a week apart would contain different pipelines.
#
# Wave's default (ephemeral) build expires roughly 24 hours after it completes, which is fine for a
# one-off validation run and useless for anything citable. For a durable image add `--freeze` and a
# `--build-repo` that Seqera has push credentials for; conf/gcb.config has the longer note.
#
# THE IMAGE MUST BE REBUILT WHENEVER pixi.lock CHANGES. There is no check that enforces this: a
# stale image simply fails at the first task that needs the new package, with an R error naming a
# missing library rather than anything about the image. nanoparquet, added 2026-09-04, is the most
# recent example -- CONSOLIDATE_REPLICATES fails without it.

# ---- stage 1: build the environment --------------------------------------------------------
#
# MULTI-STAGE, because layer arithmetic is otherwise easy to get wrong. A `RUN rm -rf` in a later
# layer records the deletion but the bytes still ship in the earlier one -- measured here, not
# assumed: a single-stage version that pruned 1.7 GB in its own RUN produced an image byte-identical
# in size to the unpruned one (2.14 GiB both ways, same 2,104 MiB layer, two extra no-op layers).
#
# With two stages the final image contains only what is COPYed, so anything the builder creates and
# then removes cannot ship regardless of how the RUNs are arranged. That is a structural guarantee
# rather than a discipline someone has to remember, and it is what actually worked -- compressed
# sizes read off the registry manifest:
#
#   single stage, no prune      7 layers   2.14 GiB
#   single stage, prune in RUN  9 layers   2.14 GiB   <- no effect whatsoever
#   two stages                  6 layers   0.68 GiB   <- 3.1x smaller
#
# The remaining 0.68 GiB is 629 MiB of /opt/pipeline, 33 MiB of pixi binary, and the 28 MiB ubuntu
# base. The base is NOT worth optimising: debian:bookworm-slim measures 26.9 MiB against ubuntu's
# 28.4 MiB, 1.5 MiB of 700, and swapping a validated base for 0.2% is a bad trade. Nothing here
# needs a heavy R base image either -- conda-forge's r-base carries R's shared libraries inside the
# environment, which is why stage 2 can be a bare ubuntu with only ca-certificates added.

FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# ca-certificates is not optional: a stock ubuntu:22.04 ships none at all, so every HTTPS call --
# pixi fetching packages, install_sceptre.R fetching a tarball from GitHub -- dies with
# "No CA certificates were loaded from the system".
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl bzip2 git \
    && rm -rf /var/lib/apt/lists/*

# Pinned rather than latest. pixi reads the lockfile format, and an older pixi cannot read a newer
# one -- the module system's 0.53.0 fails outright on this repo's format-v7 lock. Pinning keeps the
# image reproducible and matches what Sherlock runs.
ARG PIXI_VERSION=0.76.2
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION="v${PIXI_VERSION}" PIXI_NO_PATH_UPDATE=1 bash
ENV PATH="/root/.pixi/bin:${PATH}"

# Pin this to a commit SHA for any build whose output matters. `main` is a moving target.
ARG WATTEG_REF=main
RUN git clone https://github.com/EngreitzLab/WattEG.git /opt/pipeline \
    && cd /opt/pipeline \
    && git checkout "${WATTEG_REF}" \
    && git rev-parse HEAD > /opt/pipeline/.image-built-from

WORKDIR /opt/pipeline

# Mirrors workflow/slurm_executor/00_setup_env.sbatch, deliberately: the two ways of building this
# environment must not drift, or a cloud run and a Sherlock run stop being comparable.
#
#   --frozen  use pixi.lock as committed, never re-solve. A worker silently re-solving would defeat
#             the lockfile and would need network it may not have.
#   build env carries the compilers, remotes, and the packages R CMD INSTALL demands from sceptre's
#             DESCRIPTION even though the runtime path never calls them.
#   setup     installs the patched sceptre and ondisc into .pixi/rlibs, compiling both.
RUN pixi install --frozen \
    && pixi install --frozen --environment build \
    && pixi run --frozen setup \
    && pixi run --frozen check-api

# Dropped before the copy. .pixi/envs/build is 1.7 GB against 1.6 GB for the runtime environment
# and 7.5 MB for the compiled sceptre and ondisc in .pixi/rlibs -- measured inside the image. No
# task calls a compiler, and the compiled packages live in rlibs rather than inside the build
# environment, so nothing the runtime needs goes with it.
RUN rm -rf /opt/pipeline/.pixi/envs/build /root/.cache/rattler /root/.cache/pixi /root/.pixi/cache

# ---- stage 2: the runtime image ------------------------------------------------------------
#
# Same base as the builder, and the SAME PATHS. Conda environments embed absolute paths, so
# /opt/pipeline has to be /opt/pipeline in both stages or the copied environment will not activate.
#
# conf/gcb.config's beforeScript symlinks the task's ${projectDir} here. Google Batch stages
# pipeline source into a per-task clone path unrelated to where this image's environment lives, so
# without that symlink `pixi run` finds a manifest with no environment beside it and every task dies
# with "pixi: command not found". The path is load-bearing: change it here and change it there.

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# The pixi binary comes too: every process script invokes `pixi run --frozen --manifest-path ...`
# itself, so the runtime needs pixi even though the environment is already solved. That is why the
# shell-hook/entrypoint form of this pattern does not apply -- Nextflow calls pixi, not the image.
COPY --from=builder /root/.pixi/bin/pixi /root/.pixi/bin/pixi
COPY --from=builder /opt/pipeline /opt/pipeline

ENV PATH="/root/.pixi/bin:${PATH}"
WORKDIR /opt/pipeline

# Prove the runtime image works on its own terms: check-api runs in the default environment and
# reaches sceptre's unexported slots, so a copy that missed something fails the BUILD rather than
# the first task of a 6,000-task sweep.
RUN pixi run --frozen check-api \
    && echo "built from $(cat /opt/pipeline/.image-built-from)"

ENV WATTEG_ROOT=/opt/pipeline

CMD ["/bin/bash"]
