#!/bin/bash
#SBATCH --job-name=power_16_covariate
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Section 6 of paper/outline.md: can a pair's power curve be predicted from covariates, for pairs
# that were never simulated? And figures.md's Figure 5c: are the errors spread evenly, or
# concentrated?
#
# WHY IT NEEDS ITS OWN CURVE FIT
#
# The design-space scan (step 14) only fits three- and four-point grids, because reduced designs are
# what it is testing. Here `k` is the thing being PREDICTED, so it should be the best estimate
# available -- all six effect sizes at 100 replicates. That fit does not otherwise exist, so it is
# made here, once, from the scan's own 100-replicate power tables.
#
# WHAT THIS DOES AND DOES NOT ESTABLISH
#
# Held-out pairs are the cheap, in-regime test: fit on most pairs, predict the rest, and see how
# wrong it is. It costs nothing because day0's sweep already measured every pair.
#
# It does NOT establish the `trans` claim. `trans` genes are not drawn from the element's
# neighbourhood, so their expression distribution differs by construction -- and expression is the
# covariate carrying most of the signal. Held-out `cis` performance therefore says nothing about
# `trans`, which needs its own measured subset (paper/experiments.md step 5). Nor does it establish
# cross-dataset transfer, which needs the coefficients refitted on a second dataset and compared.
#
# Usage:
#   sbatch workflow/slurm_executor/16_covariate_model.sbatch
#   sbatch --dependency=afterok:<scan job> workflow/slurm_executor/16_covariate_model.sbatch

set -euo pipefail

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

# SCAN_SAMPLE, not SAMPLE -- config.sh exports SAMPLE for the sbatch runner's own dataset.
SCAN_SAMPLE="${SCAN_SAMPLE:-day0}"
SWEEP="${SCAN_SWEEP:-/home/users/emattei/oak/projects/element-gene-power-analysis/power_sweep/${SCAN_SAMPLE}/${SCAN_SAMPLE}}"
OUT="${SCAN_OUT:-${SWEEP}/reduced_design}"

log_header "Step 16: covariate model for ${SCAN_SAMPLE}"
load_pixi

THRESHOLD="${SWEEP}/prepared/discovery_threshold.txt"
REFERENCE="${SWEEP}/power_summary.tsv"
FULL_POWER="${OUT}/power/100_1"
SIX_POINT="${OUT}/curves_six_point/six_point_100.tsv"

for f in "${THRESHOLD}" "${REFERENCE}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR: missing ${f}" >&2
    exit 1
  fi
done

n_full=$(ls "${FULL_POWER}"/power_es*.tsv 2>/dev/null | wc -l)
if [[ "${n_full}" -ne 6 ]]; then
  echo "ERROR: need all six 100-replicate power tables in ${FULL_POWER}, found ${n_full}." >&2
  echo "       Run 14_design_space_scan.sbatch first -- it produces them." >&2
  exit 1
fi

mkdir -p "$(dirname "${SIX_POINT}")"

## ---- the reference k: all six effect sizes, 100 replicates -------------------------------------
#
# Kept out of ${OUT}/curves so the step-15 report, which globs that directory and counts what it
# finds, does not mistake a six-point fit for one of the reduced designs it is scoring.
if [[ ! -s "${SIX_POINT}" ]]; then
  files=""
  for es in 0.05 0.1 0.15 0.2 0.25 0.5; do
    files="${files}${files:+,}${FULL_POWER}/power_es${es}.tsv"
  done
  echo ">>> fitting the six-point reference curve"
  pixi_run Rscript src/fit_power_curve.R \
      --power "${files}" \
      --threshold-file "${THRESHOLD}" \
      --out "${SIX_POINT}"
else
  echo "  six-point fit already present: ${SIX_POINT}"
fi

## ---- the covariate model, held out, with residual structure -----------------------------------
echo
echo ">>> fitting log(k) ~ log(pert_cells) + log(expression)"
pixi_run Rscript src/fit_covariate_model.R \
    --curves "${SIX_POINT}" \
    --reference "${REFERENCE}" \
    --threshold-file "${THRESHOLD}" \
    --out-summary "${OUT}/covariate_model.txt" \
    --out-residuals "${OUT}/covariate_residuals.tsv"

echo
echo "Wrote ${OUT}/covariate_model.txt and covariate_residuals.tsv"
echo
echo "Read the decile table first. If the interval width each decile needs is roughly constant,"
echo "one global calibration is defensible and section 6 can carry per-pair intervals. If it moves"
echo "a lot, the residuals are structured -- which is the case docs/status.md warns is worse for a"
echo "downstream classifier than uniform noise, and it means calibration has to be conditional."
