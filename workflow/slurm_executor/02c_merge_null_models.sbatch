#!/bin/bash
#SBATCH --job-name=power_02c_merge_null
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 2c -- merge the per-simulation null-model chunks from step 2b into one bundle.
#
# Cheap: a bundle entry is 11 named coefficients plus a theta scalar, so 100 simulations x ~237 genes
# is a few hundred KB. The work is in the checking, not the merging -- merge_null_models.R refuses a
# set with a gap, a duplicate simulation, a mismatched seed or gene list, or coefficients that are
# identical between replicates.
#
# Usage:  sbatch workflow/slurm_executor/02c_merge_null_models.sbatch
#         (02b_submit.sh submits it with --dependency=afterok on the array)

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 2c: merge null-model chunks"

load_pixi
require_dir_nonempty "${NULL_MODELS_DIR}" "null-model chunks from step 02b"

n_chunks=$(find "${NULL_MODELS_DIR}" -name 'chunk_*.rds' | wc -l)
n_expected=$(( (NUM_REPLICATES + REPS_PER_NULL_CHUNK - 1) / REPS_PER_NULL_CHUNK ))
echo "  chunks found:    ${n_chunks} (expected ${n_expected})"
echo "  simulations:      ${NUM_REPLICATES}"
echo "  output:          ${NULL_MODELS}"
echo

if [[ "${n_chunks}" -ne "${n_expected}" ]]; then
  echo "ERROR: ${n_chunks} chunks present but ${n_expected} expected." >&2
  echo "       Re-run the missing step 02b tasks; completed chunks are skipped." >&2
  exit 1
fi

pixi_run Rscript src/merge_null_models.R \
  --inputs "${NULL_MODELS_DIR}/chunk_*.rds" \
  --reps "${NUM_REPLICATES}" \
  --out "${NULL_MODELS}"

## CHECKS =========================================================================================

echo
require_file "${NULL_MODELS}" "merged null models"
ls -la "${NULL_MODELS}" | sed 's/^/  /'

report_usage
