#!/bin/bash
#
# Submit the step 2b null-model array, and step 2c to merge it.
#
# Same reason as 04_submit.sh: the array range depends on NUM_REPLICATES / REPS_PER_NULL_CHUNK, and
# #SBATCH directives are read before config.sh is sourced, so the range cannot be a directive
# without going stale. It is computed here and passed on the command line, which overrides one.
#
# The merge is submitted with --dependency=afterok on the array, so the two can be launched together
# and the merge simply will not run if any chunk fails.
#
# Usage:
#   workflow/slurm_executor/02b_submit.sh              # array + merge
#   workflow/slurm_executor/02b_submit.sh --dry-run    # show what would be submitted
#   workflow/slurm_executor/02b_submit.sh --no-merge   # array only
#
# Any other arguments are passed through to sbatch.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

DRY_RUN=0
MERGE=1
PASSTHROUGH=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)  DRY_RUN=1; shift ;;
    --no-merge) MERGE=0; shift ;;
    *)          PASSTHROUGH+=("$1"); shift ;;
  esac
done

# Ceiling division: with 100 simulations at 1 per chunk that is 100 tasks, and a non-multiple leaves a
# short final chunk, which 02b handles.
n_tasks=$(( (NUM_REPLICATES + REPS_PER_NULL_CHUNK - 1) / REPS_PER_NULL_CHUNK ))

echo "Null models"
echo "  simulations:        ${NUM_REPLICATES}"
echo "  reps per chunk:    ${REPS_PER_NULL_CHUNK}"
echo "  array tasks:       ${n_tasks}"
echo "  chunk directory:   ${NULL_MODELS_DIR}"
echo "  merged output:     ${NULL_MODELS}"
echo "  seed:              ${SEED}"
echo

mkdir -p "${LOG_DIR}"

# ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}, not "${PASSTHROUGH[@]}": config.sh sets `set -u`, and bash
# 4.2 (what CentOS 7 ships) treats an empty array as unset, so the plain form fails outright with
# "PASSTHROUGH[@]: unbound variable" whenever this script is invoked with no extra arguments -- which
# is the normal case. Same reasoning as 04_submit.sh, where this was already fixed; see the longer
# note there on why the `+` form rather than `-`.
fit_cmd=(sbatch --array="1-${n_tasks}" ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}
         "${CONFIG_DIR}/02b_fit_null_models.sbatch")
merge_cmd_base=(sbatch ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"} "${CONFIG_DIR}/02c_merge_null_models.sbatch")

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "would run: ${fit_cmd[*]}"
  [[ "${MERGE}" -eq 1 ]] && echo "would run: ${merge_cmd_base[*]} --dependency=afterok:<arrayid>"
  exit 0
fi

out="$("${fit_cmd[@]}")"
echo "${out}"
array_id="$(echo "${out}" | awk '{print $NF}')"

if [[ "${MERGE}" -eq 1 ]]; then
  # afterok on the whole array: the merge refuses an incomplete set anyway, but not launching it at
  # all gives a clearer failure than a merge that errors on a gap.
  sbatch --dependency="afterok:${array_id}" ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"} \
    "${CONFIG_DIR}/02c_merge_null_models.sbatch"
fi
