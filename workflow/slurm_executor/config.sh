# Shared settings for the SLURM executor scripts. Sourced by every sbatch script in this
# directory, so that parameters live in exactly one place.
#
# This is not itself submittable -- it is sourced, not run.
#
# The values below mirror config/config.yml. They are duplicated here rather than parsed out of
# the YAML because parsing YAML in bash is not worth the fragility, but they must be kept in
# step: if you change one, change the other.

set -euo pipefail

## PATHS ==========================================================================================

# Resolved from this file's own location, not $PWD, so the scripts work regardless of where
# sbatch was invoked from. BASH_SOURCE[0] is config.sh even when sourced.
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${CONFIG_DIR}/../.." && pwd)"

# The sample being validated: the old pipeline was run on this one with the size-factor shuffle
# already removed, so its outputs are directly comparable to the refactored code.
SAMPLE="day0_grna20_no_shuffle"
OLD_DIR="${REPO_ROOT}/results/${SAMPLE}"
SCEPTRE_OBJECT="${OLD_DIR}/sceptre_object.rds"

# All refactored-pipeline output lands here, keeping it separate from the old results.
OUTDIR="${REPO_ROOT}/results/refactor"
PREPARED_DIR="${OUTDIR}/prepared"
SPLITS_DIR="${OUTDIR}/splits"
SMOKE_DIR="${OUTDIR}/smoke"
COMPARISON_DIR="${OUTDIR}/comparison"

## WHICH NULL MODEL THE SIMULATION TESTS AGAINST =================================================
#
# This is part of the output path, and that is load-bearing rather than tidy-mindedness.
#
# The null-model configuration does not change the *shape* of a task's output -- the same
# reps x pairs rows come out under every configuration -- so nothing about a finished file
# distinguishes one from another. Step 04's resume check compares row counts, and on 2026-08-13 it
# therefore skipped all 1,000 tasks of the null_fit rerun (job 38882207) because the as_is output
# from job 38849611 was sitting at the same paths with exactly the right number of rows. The array
# reported 1,000 successes and recomputed nothing.
#
# Putting the configuration in the directory name makes that collision impossible: two
# configurations can no longer address the same file. The `provenance` sidecar written next to each
# output in step 04 is the second layer, and catches the case a directory name cannot -- the same
# configuration re-run against a *different* null-model bundle.
#
#   as_is     inherited real-data cache. Biased low by ~0.03; only for reproducing job 38849611.
#   null_fit  per-simulation null models from step 02b. The correct configuration; production default.
#   cleared   refits inside every call. The faithful reference, 4.3x the cost. Never for production.
#
# See docs/status.md, "Settled -- which per-gene null model does the simulation test against?".
NULL_MODEL_CONFIG="null_fit"

SIM_DIR="${OUTDIR}/sim_${NULL_MODEL_CONFIG}"
POWER_DIR="${OUTDIR}/power_${NULL_MODEL_CONFIG}"

# The as_is baseline from job 38849611, kept deliberately: it is the only direct measurement of how
# much the inherited cache moved real numbers at 100 replicates. Step 09 compares against it, so it
# is named explicitly rather than reached through SIM_DIR, which now follows NULL_MODEL_CONFIG.
AS_IS_SIM_DIR="${OUTDIR}/sim_as_is"
AS_IS_POWER_DIR="${OUTDIR}/power_as_is"

# Per-replicate null models (step 02b), and the null_fit validation subset (steps 08-09). The
# subset writes to its own tree so it can never be mistaken for, or concatenated with, either the
# as_is baseline or the full null_fit run.
NULL_MODELS_DIR="${OUTDIR}/null_models"
NULL_MODELS="${PREPARED_DIR}/null_precomputations.rds"
NULL_FIT_SUBSET_DIR="${OUTDIR}/sim_null_fit_subset"

LOG_DIR="${REPO_ROOT}/logs/refactor"

## PIPELINE PARAMETERS ============================================================================

# Effect sizes as a *fractional decrease* in expression (0.15 = 15% knockdown).
#
# Only 0.15 for the validation run: it is the only effect size the old pipeline covers, so it is
# the only one that can be compared against. A typical production sweep is 5% to 50%:
#
#   EFFECT_SIZES=(0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50)
#
# Adding effect sizes costs CPU linearly (~635 CPU-hours each) but need not cost wall clock --
# see the N_SPLITS budget below, which trades splits against effect sizes to keep everything in
# one wave.
EFFECT_SIZES=(0.15)

NUM_REPLICATES=100     # Monte-Carlo simulations per pair
SEED=20250812          # base RNG seed
GUIDE_SPREAD_C=0.65    # guide knockdown ~ Beta(mean es, sd c*es*(1-es)); docs/methods.md
# fixed: each gene's realised mean over the perturbed cells pinned to es; random: es on average
# (docs/methods.md). Both write to SIM_DIR, and step 04's provenance check redoes rather than reuses
# an output made under the other one, so give a second estimand its own OUTDIR.
ESTIMAND=fixed

## NULL MODELS (step 02b) =========================================================================
#
# Simulations per array task when fitting null models. One fit covers every gene (~237) on one
# simulation and takes ~20 minutes, so 1 per task keeps wall clock at ~20 minutes for the whole set
# instead of ~33 hours serial. Raise it only if the array-task budget is tight.
REPS_PER_NULL_CHUNK=1

## NULL_FIT VALIDATION SUBSET (steps 08-09) =======================================================
#
# How many splits to re-run under null_fit to measure what the inherited cache cost us, without
# paying for a full second sweep.
#
# The comparison is exactly paired: seeds derive from (seed, target, rep, effect_size), which does
# not involve the null model, so the simulated count matrices are bit-identical between the as_is
# run and this one. Both configurations skip the GLM fit -- one from the inherited cache, one from
# the bundle -- so RNG consumption matches too. Every difference in power is therefore attributable
# to the null model, with zero Monte-Carlo noise between the two runs. That is why 20 splits are
# decisive where 53 pairs at 5 simulations were only suggestive.
NULL_FIT_SUBSET_SPLITS=20

# Parallel tasks per effect size. Does not affect results -- seeds derive from
# (seed, target, rep, effect_size), so the pipeline is invariant to the split layout. It decides
# wall clock only.
#
# THE BUDGET THAT ACTUALLY BINDS
#
# The simulation is submitted as one 2D array over (split x effect size), so the task count is
# N_SPLITS * ${#EFFECT_SIZES[@]}. That product is capped by the QOS, not by anything about this
# pipeline:
#
#   sacctmgr show qos ->  owners  MaxSubmitJobsPU=3000  MaxTRESPU=cpu=16384
#                         normal  MaxSubmitJobsPU=2000  MaxTRESPU=cpu=512
#
# Array tasks count individually against MaxSubmitJobsPU, so a 10-effect-size sweep at 1000
# splits (10,000 tasks) is rejected at submission. MAX_ARRAY_TASKS below leaves headroom under
# the 3000 cap for the odd interactive job.
#
# Because ~all tasks run concurrently, wall clock is roughly one task's duration regardless of
# how the budget is divided -- so spend it.
#
# MEASURED cost, fitted on all 999 real tasks of job 38887744 (null_fit + gRNA precomputation
# reuse): 1.140 s per target + 0.5561 s per pair, per replicate. Over 3,026 targets and 34,886
# pairs at 100 simulations that is 635 CPU-hours per effect size predicted against 634 measured,
# down from 999 before the gRNA reuse and against the old pipeline's 1,308 -- a 2.1x speedup.
#
# Do not re-derive this from the step 3 smoke test. That gave 9.68 s per target per simulation and
# implied ~814 CPU-hours; it was measured before the gRNA reuse and on 36 targets.
#
#   effect sizes   N_SPLITS   tasks   per task    wall clock (all effect sizes)
#   1              1000       1000    ~36 min     ~40 min
#   1              2000       2000    ~18 min     ~20 min
#   2              1400       2800    ~26 min     ~30 min
#   10             290        2900    ~131 min    ~2.2 h      (6,350 CPU-hours!)
#
# Rule of thumb: N_SPLITS ~= MAX_ARRAY_TASKS / number of effect sizes. And note the last row:
# a full 5-50% sweep is 8,000+ CPU-hours simulated naively. See docs/status.md on two-stage
# simulation allocation, and the cascading-effect-size idea, for how to avoid paying it.
#
# Two other ceilings: whole targets stay together, so N_SPLITS can never usefully exceed the
# 2,798 targets; and below ~5 minutes per task the fixed costs (container start, R startup, the
# 18 MB + 50 MB reads) start to dominate. See docs/status.md on preferring more splits over more
# rep chunks.
N_SPLITS=1000

# Headroom below the QOS MaxSubmitJobsPU of 3000. check_array_budget() enforces it.
MAX_ARRAY_TASKS=2900

# Total tasks in the step 4 array, and the (split, effect size) decomposition it implies.
# Ordering is split-major within each effect size, so tasks 1..N_SPLITS are the first effect
# size, N_SPLITS+1..2*N_SPLITS the second, and so on. That makes it possible to submit or
# resubmit a single effect size by array range.
N_EFFECT_SIZES=${#EFFECT_SIZES[@]}
TOTAL_SIM_TASKS=$(( N_SPLITS * N_EFFECT_SIZES ))

# Refuse to submit an array the scheduler will reject, and say what to do about it.
check_array_budget() {
  if [[ "${TOTAL_SIM_TASKS}" -gt "${MAX_ARRAY_TASKS}" ]]; then
    echo "ERROR: ${N_SPLITS} splits x ${N_EFFECT_SIZES} effect sizes = ${TOTAL_SIM_TASKS} array tasks," >&2
    echo "       which exceeds the ${MAX_ARRAY_TASKS}-task budget (QOS MaxSubmitJobsPU is 3000)." >&2
    echo "       Lower N_SPLITS in config.sh to about $(( MAX_ARRAY_TASKS / N_EFFECT_SIZES )) and re-run step 02." >&2
    echo "       Wall clock is unaffected: fewer, longer tasks all still run concurrently." >&2
    exit 1
  fi
}

## SLURM ==========================================================================================

# engreitz alone has 288 cores and was queueing for a long time; owners adds a few thousand, at
# the cost of being preemptible and capped at 2 days. Every request here is well under that, and
# step 4 skips already-complete tasks on resubmission, so a preempted task costs only its own
# work.
PARTITION="engreitz,owners"

# The simulation array is deliberately NOT throttled -- see the --array line in
# 04_run_power_simulation.sbatch. Throttling to %50 turned a ~20 minute run into a ~6 hour one.
# Slurm's MaxArraySize here is 1,000,000, so the array size itself is not a constraint.

## ENVIRONMENT ====================================================================================
#
# WHY EVERYTHING RUNS IN A CONTAINER
#
# Sherlock is CentOS 7 (glibc 2.17) on every node, login and compute alike. conda-forge dropped
# glibc 2.17 support after CentOS 7 reached end of life, so the r-base 4.4 builds this project
# pins require glibc >= 2.28 and pixi refuses to install them:
#
#   × Cannot install environment 'default'
#     Virtual package '__glibc >=2.28' does not match any of the available virtual packages on
#     your machine: [__glibc=2.17=0, ...]
#
# This is not something pixi can fix. glibc is a *virtual* package: conda packages link against
# the host's C library through the system dynamic loader, so pixi can only detect the host glibc,
# never supply a newer one. Refusing to install binaries that would fail at load time is the
# correct behaviour. The same is true of any userspace version manager -- mise, asdf, conda --
# for the same reason.
#
# So the environment runs inside an Ubuntu 22.04 container (glibc 2.35), which does supply its
# own C library. The alternatives were to re-solve the lockfile against packages old enough for
# CentOS 7 (discards the pinned solve, and r-base 4.4 may have no such build), or to abandon pixi
# for Sherlock's own R module (drops the lockfile as a description of what ran). The container
# keeps pixi.lock authoritative, and Nextflow supports Apptainer natively, so this carries
# straight into the Nextflow step rather than being throwaway scaffolding.
#
# The image is a stock ubuntu:22.04 with nothing added -- the pixi binary is statically linked,
# so it is simply bind-mounted in from $GROUP_HOME. To rebuild it:
#
#   APPTAINER_CACHEDIR=$SCRATCH/apptainer_cache \
#     apptainer pull $GROUP_HOME/software/containers/ubuntu2204.sif docker://ubuntu:22.04

CONTAINER_IMAGE="${GROUP_HOME}/software/containers/ubuntu2204.sif"

# Bind mounts. Apptainer binds $HOME and the current directory by default but not $GROUP_HOME or
# $SCRATCH, and every path this pipeline touches -- repo, inputs, outputs, pixi cache -- lives
# under one of those two.
#
# The third bind is CA certificates. A stock ubuntu:22.04 image ships no ca-certificates package
# at all -- /etc/ssl/certs does not exist -- so pixi cannot open an HTTPS connection and dies
# with "No CA certificates were loaded from the system". Rather than build a derived image just
# to add them, the host's bundle is mounted in and pointed at with SSL_CERT_FILE. This matters
# only while downloading; once .pixi/ is populated the runtime steps need no network.
HOST_CA_BUNDLE="/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem"
CONTAINER_CA_BUNDLE="/etc/ssl/certs/ca-certificates.crt"

CONTAINER_BINDS=(
  --bind "${GROUP_HOME}:${GROUP_HOME}"
  --bind "${SCRATCH}:${SCRATCH}"
  --bind "${HOST_CA_BUNDLE}:${CONTAINER_CA_BUNDLE}:ro"
)

# Bring pixi onto PATH.
#
# The module system only carries pixi/0.53.0, which reads lockfile format v6 and earlier. The
# committed pixi.lock is format v7, so 0.53.0 fails outright:
#
#   × Failed to load lock file ... The lock file version is 7, but only up to including
#     version 6 is supported by the current version.
#
# Rather than re-solve the lockfile with the older pixi -- which would discard the pinned solve
# that makes this pipeline reproducible -- a current pixi is installed as a standalone binary
# under $GROUP_HOME. It is a single static executable and needs no root.
#
# To reinstall or upgrade it:
#   PIXI_HOME=$GROUP_HOME/software/pixi curl -fsSL https://pixi.sh/install.sh | bash
#
# Note that installer appends a PATH line to ~/.bashrc. Remove it: $GROUP_HOME is NFS, and PATH
# entries on slow filesystems are a documented cause of hung logins and held jobs on Sherlock.
# PATH is set explicitly below instead, which is also what makes this work identically in a
# batch job, where ~/.bashrc may not be sourced at all.
PIXI_BIN_DIR="${GROUP_HOME}/software/pixi/bin"

load_pixi() {
  if [[ ! -x "${PIXI_BIN_DIR}/pixi" ]]; then
    echo "ERROR: no pixi at ${PIXI_BIN_DIR}/pixi." >&2
    echo "       Install it with:" >&2
    echo "         PIXI_HOME=${GROUP_HOME}/software/pixi curl -fsSL https://pixi.sh/install.sh | bash" >&2
    echo "       then delete the PATH line it appends to ~/.bashrc." >&2
    echo "       The pixi/0.53.0 module is not a substitute: it cannot read a v7 pixi.lock." >&2
    exit 1
  fi
  if [[ ! -f "${CONTAINER_IMAGE}" ]]; then
    echo "ERROR: container image not found at ${CONTAINER_IMAGE}." >&2
    echo "       Build it with:" >&2
    echo "         APPTAINER_CACHEDIR=\$SCRATCH/apptainer_cache \\" >&2
    echo "           apptainer pull ${CONTAINER_IMAGE} docker://ubuntu:22.04" >&2
    exit 1
  fi

  # Pixi caches downloaded packages under $HOME/.cache by default. $HOME is 15 GB and
  # NFS-backed, and a conda-forge R stack will fill a meaningful fraction of it, so the cache
  # is redirected to scratch. Do this before any pixi invocation.
  export PIXI_CACHE_DIR="${SCRATCH}/pixi_cache"
  export XDG_CACHE_HOME="${SCRATCH}/xdg_cache"
  export RATTLER_CACHE_DIR="${PIXI_CACHE_DIR}"

  # Apptainer's own cache and build scratch, for the same reason.
  export APPTAINER_CACHEDIR="${SCRATCH}/apptainer_cache"
  export APPTAINER_TMPDIR="${SCRATCH}/apptainer_tmp"

  # HOME inside the container. Anything that insists on writing to a home directory -- pixi's own
  # state, R's session files -- lands here rather than in the real $HOME, which is 15 GB.
  export CONTAINER_HOME="${SCRATCH}/container_home"

  mkdir -p "${PIXI_CACHE_DIR}" "${XDG_CACHE_HOME}" "${CONTAINER_HOME}" \
           "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

  echo "pixi:       ${PIXI_BIN_DIR}/pixi ($("${PIXI_BIN_DIR}/pixi" --version 2>/dev/null))"
  echo "container:  ${CONTAINER_IMAGE}"
}

# Run a command inside the container, inside the pixi environment, from the repo root.
#
# --frozen means "use pixi.lock exactly as committed, never re-solve". That is the whole point of
# the lockfile: a compute node silently re-solving the environment would defeat reproducibility,
# and would also need network access it may not have.
#
# --cleanenv starts from an empty environment rather than inheriting the host's. Without it the
# host's Lmod variables, PATH and R settings leak into the container and can shadow the pixi
# environment -- the failure mode is an R that half-works, which is worse than one that does not
# start. Everything the run genuinely needs is passed explicitly with --env.
container_run() {
  # HOME is set with --home, not --env: Apptainer refuses the latter outright
  # ("Overriding HOME environment variable with APPTAINERENV_HOME is not permitted").
  apptainer exec --cleanenv \
    "${CONTAINER_BINDS[@]}" \
    --home "${CONTAINER_HOME}" \
    --pwd "${REPO_ROOT}" \
    --env "PIXI_CACHE_DIR=${PIXI_CACHE_DIR}" \
    --env "RATTLER_CACHE_DIR=${RATTLER_CACHE_DIR}" \
    --env "XDG_CACHE_HOME=${XDG_CACHE_HOME}" \
    --env "SSL_CERT_FILE=${CONTAINER_CA_BUNDLE}" \
    --env "SSL_CERT_DIR=/etc/ssl/certs" \
    --env "REQUESTS_CA_BUNDLE=${CONTAINER_CA_BUNDLE}" \
    --env "CURL_CA_BUNDLE=${CONTAINER_CA_BUNDLE}" \
    --env "PATH=${PIXI_BIN_DIR}:/usr/local/bin:/usr/bin:/bin" \
    "${CONTAINER_IMAGE}" "$@"
}

pixi_run() {
  cd "${REPO_ROOT}"
  container_run pixi run --frozen "$@"
}

# Same, but for the environments that are not the default one (currently only `build`, which
# carries the compiler toolchain used to compile sceptre).
pixi_run_env() {
  local env_name="$1"; shift
  cd "${REPO_ROOT}"
  container_run pixi run --frozen --environment "${env_name}" "$@"
}

## HELPERS ========================================================================================

log_header() {
  echo "=================================================================="
  echo "  $*"
  echo "  job ${SLURM_JOB_ID:-<none>}${SLURM_ARRAY_TASK_ID:+ task ${SLURM_ARRAY_TASK_ID}}"
  echo "  host $(hostname)  started $(date '+%Y-%m-%d %H:%M:%S')"
  echo "=================================================================="
}

# Fail early and legibly rather than letting an R script error on a missing input.
require_file() {
  local path="$1" what="${2:-input}"
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: ${what} is missing or empty: ${path}" >&2
    echo "       Has the preceding step in workflow/slurm_executor/ been run?" >&2
    exit 1
  fi
}

require_dir_nonempty() {
  local path="$1" what="${2:-directory}"
  if [[ ! -d "${path}" ]] || [[ -z "$(ls -A "${path}" 2>/dev/null)" ]]; then
    echo "ERROR: ${what} is missing or empty: ${path}" >&2
    exit 1
  fi
}

# Report what the job actually used, so the resource requests can be tightened from measurement
# rather than guesswork. See docs/status.md, Comparison 3.
report_usage() {
  echo
  echo "------------------------------------------------------------------"
  echo "  finished $(date '+%Y-%m-%d %H:%M:%S')"
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "  resource usage (sacct may lag a few seconds):"
    sacct -j "${SLURM_JOB_ID}" --format=JobID%20,Elapsed,MaxRSS,State 2>/dev/null || true
  fi
  echo "------------------------------------------------------------------"
}
