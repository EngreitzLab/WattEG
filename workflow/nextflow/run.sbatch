#!/bin/bash
#SBATCH --job-name=power_nextflow
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/nextflow/%x-%A.out
#SBATCH --error=logs/nextflow/%x-%A.err
#
# Drive the Nextflow pipeline from inside a job.
#
# The driver is a JVM that lives for the whole run, submitting and reaping Slurm jobs as it goes.
# That is exactly the kind of long-lived process that must not sit on a login node -- Sherlock's
# login nodes are shared and reap this sort of thing. It needs little CPU but a walltime at least
# as long as the pipeline, hence 2 days against 8 GB.
#
# Usage:
#   sbatch workflow/nextflow/run.sbatch                                  # full run
#   PARAMS_FILE=config/test.yml sbatch workflow/nextflow/run.sbatch      # the fast path
#   sbatch workflow/nextflow/run.sbatch -resume                          # pick up where it stopped
#
# Anything passed on the command line is appended to the `nextflow run` invocation, so -resume,
# -entry, --n_splits 100 and friends all work.
#
# The test settings live in a params file, not a profile: Nextflow gives -params-file precedence
# over profile params, and this script always passes one, so a `-profile test` would have been
# silently overridden by config/config.yml and produced a full-size run. Profiles choose where a run
# executes; params files choose what it runs.

set -euo pipefail

# Slurm copies the submitted script into its spool directory, so $BASH_SOURCE points there rather
# than into the repository -- `cd` on it lands somewhere unwritable. Ask the controller where the
# script actually came from instead. Same idiom as workflow/slurm_executor/*.sbatch.
_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
REPO_ROOT="$(cd "$(dirname "${_cmd:-${BASH_SOURCE[0]}}")/../.." && pwd)"
cd "${REPO_ROOT}"
mkdir -p logs/nextflow

ml biology nextflow/25.04.7

# Everything Nextflow writes goes to Lustre, not to $HOME.
#
# NXF_HOME defaults to ~/.nextflow, where it caches plugins and assets; NXF_WORK defaults to
# ./work, which would put every task's staged inputs and outputs on the same filesystem as the
# repo. $HOME is 15 GB and NFS-backed, and a 1,000-task run with a 448 MB input object staged per
# task would fill it. $SCRATCH is 100 TB.
export NXF_HOME="${GROUP_HOME}/software/nextflow_home"
export NXF_WORK="${SCRATCH}/nextflow_work/watteg"
export NXF_TEMP="${SCRATCH}/nextflow_tmp"
export NXF_OPTS='-Xms1g -Xmx4g'

# Apptainer's own caches, for the same reason.
export APPTAINER_CACHEDIR="${SCRATCH}/apptainer_cache"
export APPTAINER_TMPDIR="${SCRATCH}/apptainer_tmp"

mkdir -p "${NXF_HOME}" "${NXF_WORK}" "${NXF_TEMP}" "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

echo "=================================================================="
echo "  Nextflow driver"
echo "  job ${SLURM_JOB_ID:-<none>}  host $(hostname)  $(date '+%Y-%m-%d %H:%M:%S')"
echo "  nextflow: $(nextflow -v 2>&1)"
echo "  work:     ${NXF_WORK}"
echo "=================================================================="

# -profile sherlock unless the caller supplies its own -profile.
PROFILE_ARGS=()
if [[ ! " $* " =~ " -profile " ]]; then
  PROFILE_ARGS=(-profile sherlock)
fi

PARAMS_FILE="${PARAMS_FILE:-config/config.yml}"
if [[ ! -s "${PARAMS_FILE}" ]]; then
  echo "ERROR: params file not found: ${PARAMS_FILE}" >&2
  exit 1
fi
echo "  params:   ${PARAMS_FILE}"

nextflow run "${REPO_ROOT}" \
    "${PROFILE_ARGS[@]}" \
    -params-file "${PARAMS_FILE}" \
    -ansi-log false \
    "$@"

echo
echo "finished $(date '+%Y-%m-%d %H:%M:%S')"
