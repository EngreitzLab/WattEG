#!/bin/bash
#SBATCH --job-name=power_07_compare
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 7 -- validate the refactored pipeline against the old one.
#
# This is the whole point of the exercise: the old pipeline was run on this sample with the
# size-factor shuffle already removed, so old and new implement the same statistical model and
# their power estimates should agree up to Monte-Carlo noise.
#
# The comparison cannot be exact. Seeding moved inside the simulation loop during the refactor, so
# the RNG streams differ, and the old pipeline has no set.seed() at all -- two runs of it do not
# agree with each other either. See workflow/compare_old_new.R for what is checked exactly and
# what is checked distributionally.
#
# Exits non-zero if any check fails, so the job state reflects the verdict.
#
# Usage:  sbatch workflow/slurm_executor/07_compare_old_new.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 7: compare refactored output against the old pipeline"

load_pixi
mkdir -p "${COMPARISON_DIR}"

# Only 0.15 exists in the old run, so that is the only effect size that can be validated.
ES="0.15"

OLD_POWER="${OLD_DIR}/power_analysis/power_analysis_results_es_${ES}.tsv"

# WHICH SIDE TO COMPARE AGAINST, AND WHY IT IS THE as_is ONE
#
# The old pipeline inherits @response_precomputations from the real discovery analysis, exactly as
# the refactored code did before null_fit. So the like-for-like comparison -- the one that isolates
# the refactor rather than the refactor plus a deliberate change of null model -- is old against
# as_is. That is the run that gates this job's exit status.
#
# The comparison against null_fit is also run, because null_fit is what the pipeline produces now
# and a reader will ask. It carries a known, measured offset of +0.0063 mean power (see
# docs/status.md, Consequences), which is well inside compare_old_new.R's 0.02 median-offset
# tolerance, so it should also pass -- but it is reported separately and does not gate.
AS_IS_POWER="${AS_IS_POWER_DIR}/power_es${ES}.tsv"
NULL_FIT_POWER="${POWER_DIR}/power_es${ES}.tsv"

require_file "${OLD_POWER}"      "old pipeline power results"
require_file "${AS_IS_POWER}"    "as_is power results (job 38849611 via step 05)"
require_file "${NULL_FIT_POWER}" "null_fit power results from step 05"

echo "old:       ${OLD_POWER}"
echo "as_is:     ${AS_IS_POWER}    <- gates the verdict"
echo "null_fit:  ${NULL_FIT_POWER}  <- reported, does not gate"
echo "out:       ${COMPARISON_DIR}"
echo

# config.sh sets `set -e`, which would abort the job before the verdict could be reported. The
# non-zero exit is the expected signal for a failed comparison, not an error to be swallowed, so
# it is captured and re-raised at the end.
set +e

echo "=================================================================="
echo "  OLD vs as_is  -- like for like, both inherit the real-data cache"
echo "=================================================================="
pixi_run Rscript workflow/compare_old_new.R \
  --old "${OLD_POWER}" \
  --new "${AS_IS_POWER}" \
  --outdir "${COMPARISON_DIR}/old_vs_as_is" \
  --min-correlation 0.95
status=$?

echo
echo "=================================================================="
echo "  OLD vs null_fit -- what the pipeline produces now (informational)"
echo "=================================================================="
pixi_run Rscript workflow/compare_old_new.R \
  --old "${OLD_POWER}" \
  --new "${NULL_FIT_POWER}" \
  --outdir "${COMPARISON_DIR}/old_vs_null_fit" \
  --min-correlation 0.95
null_fit_status=$?
echo
echo "null_fit comparison exit status: ${null_fit_status} (informational, does not gate)"

set -e

echo
if [[ ${status} -eq 0 ]]; then
  echo "VALIDATION PASSED -- the refactored pipeline agrees with the old one."
  echo "Nextflow development (Step 7 in docs/status.md) can begin."
else
  echo "VALIDATION FAILED -- see ${COMPARISON_DIR}/old_vs_as_is/comparison_report.txt"
  echo "Do not start Nextflow work until the discrepancy is understood."
fi

report_usage
exit ${status}
