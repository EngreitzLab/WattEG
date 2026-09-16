#!/bin/bash
#
# Submit the step 4 simulation array.
#
# This wrapper exists because the array range depends on N_SPLITS x the number of effect sizes,
# and #SBATCH directives are read before the script runs, so they cannot reference config.sh.
# Rather than keep a literal --array in sync by hand -- a footgun that silently skips work when
# it drifts -- the range is computed here and passed on the command line, which overrides any
# directive.
#
# Usage:
#   workflow/slurm_executor/04_submit.sh                # everything
#   workflow/slurm_executor/04_submit.sh --dry-run      # show what would be submitted
#   workflow/slurm_executor/04_submit.sh --effect-size 0.15
#                                                       # just that one effect size's range
#
# Any other arguments are passed through to sbatch, e.g. --dependency=afterok:12345.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

DRY_RUN=0
ONLY_ES=""
PASSTHROUGH=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)     DRY_RUN=1; shift ;;
    --effect-size) ONLY_ES="$2"; shift 2 ;;
    *)             PASSTHROUGH+=("$1"); shift ;;
  esac
done

check_array_budget

## WORK OUT THE ARRAY RANGE =======================================================================

if [[ -n "${ONLY_ES}" ]]; then
  # Split-major ordering means one effect size is a contiguous range.
  es_idx=-1
  for i in "${!EFFECT_SIZES[@]}"; do
    [[ "${EFFECT_SIZES[$i]}" == "${ONLY_ES}" ]] && es_idx=$i
  done
  if [[ "${es_idx}" -lt 0 ]]; then
    echo "ERROR: effect size ${ONLY_ES} is not in EFFECT_SIZES (${EFFECT_SIZES[*]})." >&2
    exit 1
  fi
  array_start=$(( es_idx * N_SPLITS + 1 ))
  array_end=$(( (es_idx + 1) * N_SPLITS ))
else
  array_start=1
  array_end=${TOTAL_SIM_TASKS}
fi

n_tasks=$(( array_end - array_start + 1 ))

## REPORT =========================================================================================

echo "step 4 -- power simulation"
echo "  effect sizes:   ${EFFECT_SIZES[*]}${ONLY_ES:+  (submitting only ${ONLY_ES})}"
echo "  splits:         ${N_SPLITS}"
echo "  simulations:     ${NUM_REPLICATES}"
echo "  array:          ${array_start}-${array_end}  (${n_tasks} tasks)"
echo "  partition:      ${PARTITION}"
echo

# The array is deliberately unthrottled: with ~all tasks running concurrently the wall clock is
# one task rather than the sum. Throttling to %50 turned a ~20 minute run into a ~6 hour one.
if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "(dry run) would submit:"
  # The `+` form for the same reason as at the real sbatch call below -- here it also keeps the
  # dry-run line free of the stray double space an empty expansion would leave.
  echo "  sbatch --array=${array_start}-${array_end} ${PASSTHROUGH[@]+${PASSTHROUGH[*]} }${CONFIG_DIR}/04_run_power_simulation.sbatch"
  exit 0
fi

# Splits must exist and match N_SPLITS, or the array indexes into the wrong thing.
require_dir_nonempty "${SPLITS_DIR}" "splits directory from step 02"
n_available=$(find "${SPLITS_DIR}" -name 'split_*.tsv' | wc -l)
if [[ "${n_available}" -ne "${N_SPLITS}" ]]; then
  echo "ERROR: ${n_available} split files present but N_SPLITS=${N_SPLITS}." >&2
  echo "       Re-run step 02 first." >&2
  exit 1
fi

# ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}, not "${PASSTHROUGH[@]-}": config.sh sets `set -u`, and
# bash 4.2 (what CentOS 7 ships) treats an empty array as unset, so both forms avoid the
# unbound-variable error -- but the `-` form substitutes a single *empty-string* argument, which
# sbatch then tries to open as the batch script ("sbatch: error: Unable to open file"). The `+`
# form expands to nothing at all when the array is empty, and preserves word boundaries when it
# is not.
jobid=$(sbatch --parsable \
  --array="${array_start}-${array_end}" \
  ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"} \
  "${CONFIG_DIR}/04_run_power_simulation.sbatch")

echo "submitted ${jobid}"
echo
echo "  watch:    squeue --me -j ${jobid}"
echo "  failures: sacct -j ${jobid} --format=JobID%16,State,Elapsed,MaxRSS | grep -v COMPLETED"
echo
echo "Completed tasks are skipped on resubmission, so re-running this script after a partial"
echo "failure or a preemption picks up only what is missing."
