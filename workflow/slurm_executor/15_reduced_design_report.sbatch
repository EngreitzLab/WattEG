#!/bin/bash
#SBATCH --job-name=power_15_report
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 3 of paper/experiments.md: score every design the scan produced against the full sweep.
#
# Runs after 14_design_space_scan.sbatch and is meant to be submitted with
# `--dependency=afterok:<scan job>`, so it never scores a partial scan. That matters more than it
# looks: the scan writes one curve file per design and this script globs the directory, so a scan
# that died halfway would silently produce a report over whichever designs happened to finish --
# which would look like a complete answer with some designs quietly missing.
#
# The guard below is belt and braces for the case where it is run by hand without the dependency.
#
# Usage:
#   sbatch --dependency=afterok:<scan job id> workflow/slurm_executor/15_reduced_design_report.sbatch
#   SCAN_SAMPLE=day2 sbatch ... (same variable as the scan; see the note there about SAMPLE)

set -euo pipefail

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

# SCAN_SAMPLE, not SAMPLE -- config.sh above exports SAMPLE for the sbatch runner's own dataset.
SCAN_SAMPLE="${SCAN_SAMPLE:-day0}"
SWEEP="${SCAN_SWEEP:-/home/users/emattei/oak/projects/element-gene-power-analysis/power_sweep/${SCAN_SAMPLE}/${SCAN_SAMPLE}}"
OUT="${SCAN_OUT:-${SWEEP}/reduced_design}"

log_header "Step 15: scoring reduced designs for ${SCAN_SAMPLE}"
load_pixi

CURVES="${OUT}/curves"
REFERENCE="${SWEEP}/power_summary.tsv"
THRESHOLD="${SWEEP}/prepared/discovery_threshold.txt"

for f in "${REFERENCE}" "${THRESHOLD}"; do
  if [[ ! -s "${f}" ]]; then
    echo "ERROR: missing required input ${f}" >&2
    exit 1
  fi
done

# How many designs the scan was supposed to produce: one curve per (grid, reps, seed), with the
# 100-replicate row run at a single seed because subsampling to 100 of 100 is the identity.
# Checking the count rather than just "is the directory non-empty" is the difference between
# noticing a truncated scan and reporting on one.
n_curves=$(find "${CURVES}" -name '*.tsv' 2>/dev/null | wc -l)
n_grids=4
n_reps_multi=4          # 10, 20, 30, 50
n_seeds=3
expected=$(( n_grids * (n_reps_multi * n_seeds + 1) ))
echo "  curves found:    ${n_curves}"
echo "  curves expected: ${expected}"
if [[ "${n_curves}" -lt "${expected}" ]]; then
  echo "ERROR: the scan looks incomplete (${n_curves} of ${expected} designs)." >&2
  echo "       Re-run 14_design_space_scan.sbatch -- it skips finished cells, so it resumes." >&2
  exit 1
fi

echo
echo ">>> scoring"
pixi_run Rscript src/compare_reduced_design.R \
    --curves "${CURVES}" \
    --reference "${REFERENCE}" \
    --threshold-file "${THRESHOLD}" \
    --out "${OUT}/reduced_design_summary.tsv"

echo
echo "Wrote ${OUT}/reduced_design_summary.tsv"
echo
echo "Read it as: MAE is held-out prediction error against the measured value, and noise_floor is"
echo "the binomial standard error of that measurement. A design whose MAE is close to the floor is"
echo "as good as measuring; mdes_ceiling_* is what the reference's own noise allows, so a design"
echo "near the ceiling is measurement-limited rather than model-limited."
