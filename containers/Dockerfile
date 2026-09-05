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

FROM ubuntu:22.04

# glibc 2.35 here, against the 2.17 on Sherlock's CentOS 7 nodes. That difference is the whole
# reason the Sherlock runner also uses a container -- conda-forge's r-base 4.4 needs >= 2.28 and
# glibc is a virtual package, so no userspace manager can supply it. See conf/sherlock.config.

ENV DEBIAN_FRONTEND=noninteractive

# ca-certificates is not optional: a stock ubuntu:22.04 ships none at all, and every HTTPS call --
# pixi fetching packages, install_sceptre.R fetching a tarball from GitHub -- dies with
# "No CA certificates were loaded from the system".
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        bzip2 \
        git \
    && rm -rf /var/lib/apt/lists/*

# Pinned rather than latest. pixi reads the lockfile format, and an older pixi cannot read a newer
# one -- the module system's 0.53.0 fails outright on this repo's format-v7 lock. Pinning here keeps
# the image reproducible and matches what Sherlock runs.
ARG PIXI_VERSION=0.76.2
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION="v${PIXI_VERSION}" PIXI_NO_PATH_UPDATE=1 bash
ENV PATH="/root/.pixi/bin:${PATH}"

# Pin this to a commit SHA for any build whose output matters. `main` is a moving target.
ARG WATTEG_REF=main
RUN git clone https://github.com/EngreitzLab/WattEG.git /opt/pipeline \
    && cd /opt/pipeline \
    && git checkout "${WATTEG_REF}" \
    && git rev-parse HEAD > /opt/pipeline/.image-built-from \
    && echo "built from $(cat /opt/pipeline/.image-built-from)"

WORKDIR /opt/pipeline

# Mirrors workflow/slurm_executor/00_setup_env.sbatch exactly, deliberately: the two ways of
# building this environment must not drift, or a cloud run and a Sherlock run stop being comparable.
#
#   --frozen  use pixi.lock as committed, never re-solve. A worker silently re-solving would defeat
#             the lockfile and would need network it may not have.
#   build env carries the compilers, remotes, and the packages R CMD INSTALL demands from sceptre's
#             DESCRIPTION even though the runtime path never calls them.
#   setup     installs the patched sceptre and ondisc into .pixi/rlibs, compiling both.
# ONE RUN, deliberately. Docker layers are additive: a `RUN rm -rf` in a LATER layer records the
# deletion but the bytes still ship in the earlier layer, so the image does not shrink at all. That
# was measured, not assumed -- an earlier version of this file pruned in separate RUN steps and
# produced an image byte-identical in size to the unpruned one (2.14 GiB both ways, same 2,104 MiB
# layer, just two extra no-op layers). Everything that creates large intermediate state and
# everything that removes it has to live in the same layer.
#
# Mirrors workflow/slurm_executor/00_setup_env.sbatch, deliberately: the two ways of building this
# environment must not drift, or a cloud run and a Sherlock run stop being comparable.
#
#   --frozen  use pixi.lock as committed, never re-solve. A worker silently re-solving would defeat
#             the lockfile and would need network it may not have.
#   build env carries the compilers, remotes, and the packages R CMD INSTALL demands from sceptre's
#             DESCRIPTION even though the runtime path never calls them.
#   setup     installs the patched sceptre and ondisc into .pixi/rlibs, compiling both.
#
# Then, in the same layer:
#
#   check-api      once before the prune, proving the build is sound
#   rm envs/build  1.7 GB against 1.6 GB for the runtime environment and 7.5 MB for the compiled
#                  sceptre and ondisc in .pixi/rlibs -- measured inside the image. No task calls a
#                  compiler. Safe because the compiled packages live in rlibs, not in the build
#                  environment, and both environments resolve the same r-base from one lockfile.
#   rm caches      package caches are only needed while installing
#   check-api      AGAIN, after the prune, in the default environment. It reaches sceptre's
#                  unexported slots, so it fails the BUILD loudly if the prune took something the
#                  runtime needed, rather than shipping an image that dies on the first task.
#
# This matters more than a one-off download: every Batch worker pulls this image, and spot
# preemption makes that far more frequent than the task count suggests -- a single 68-CPU-hour
# validation run took 21 interruptions.
RUN pixi install --frozen \
    && pixi install --frozen --environment build \
    && pixi run --frozen setup \
    && pixi run --frozen check-api \
    && rm -rf /opt/pipeline/.pixi/envs/build \
    && rm -rf /root/.cache/rattler /root/.cache/pixi /root/.pixi/cache \
    && pixi run --frozen check-api

ENV WATTEG_ROOT=/opt/pipeline

CMD ["/bin/bash"]
