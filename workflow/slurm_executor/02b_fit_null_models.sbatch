#!/bin/bash
#SBATCH --job-name=power_02b_null_models
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A_%a.out
#SBATCH --error=logs/refactor/%x-%A_%a.err
#
# Step 2b -- fit the per-gene null model on a null simulation of each replicate.
#
# WHY IT IS A STEP OF ITS OWN
#
# sceptre skips fitting a gene's null model when @response_precomputations already holds an entry for
# it. The template used to inherit those entries from the *real* discovery analysis, so simulated
# counts were tested against real-data coefficients -- which measurably understates power (7 of 265
# discovery calls flipped against a faithful refit, all in the same direction; see docs/status.md).
#
# A null simulation applies no knockdown, so its fit depends on neither the target nor the effect
# size. One fit per (gene, simulation) therefore serves every target the gene pairs with (~147 on
# average), every split, and the whole effect-size sweep. Fitting inside the simulation array would
# instead pay for it once per task -- 1,000 times over.
#
# WHY AN ARRAY
#
# Measured: one simulation fits all 237 genes over 586,309 cells in 1,025s. Serial, 100 simulations
# would be ~28.5 CPU-hours in one 28-hour job; as an array of one simulation per task it is ~17
# minutes of wall clock for the same CPU. The chunks are merged by step 02c.
#
# RESOURCES
#
# Measured peak R heap is 6.6 GB -- more than the ~2.2 GB the count matrix (237 x 586,309 doubles)
# and the all-ones effect-size matrix account for, the rest being sceptre's own working copies. 10G
# leaves ~1.5x headroom. Scale it with the gene count if a dataset has many more than 237 tested genes,
# since the count matrix is genes x cells. 2 hours against a measured 17 minutes covers a slow node.
#
# Usage:  sbatch --array=1-100 workflow/slurm_executor/02b_fit_null_models.sbatch
#         (or let 02b_submit.sh compute the range from config.sh)

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: no SLURM_ARRAY_TASK_ID -- this script must be submitted as an array." >&2
  echo "       Use: workflow/slurm_executor/02b_submit.sh" >&2
  exit 1
fi

## THIS TASK'S REPLICATE RANGE ====================================================================

rep_offset=$(( (SLURM_ARRAY_TASK_ID - 1) * REPS_PER_NULL_CHUNK ))
reps_here=${REPS_PER_NULL_CHUNK}

# The last chunk is short whenever NUM_REPLICATES is not a multiple of REPS_PER_NULL_CHUNK. Without
# this the merge would find simulations beyond NUM_REPLICATES and reject the whole set.
if (( rep_offset + reps_here > NUM_REPLICATES )); then
  reps_here=$(( NUM_REPLICATES - rep_offset ))
fi
if (( reps_here <= 0 )); then
  echo "Task ${SLURM_ARRAY_TASK_ID} has no simulations to fit (offset ${rep_offset} of ${NUM_REPLICATES}); nothing to do."
  exit 0
fi

log_header "Step 2b: null models for simulations $((rep_offset + 1))-$((rep_offset + reps_here))"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"        "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"

mkdir -p "${NULL_MODELS_DIR}"
chunk_out="${NULL_MODELS_DIR}/$(printf 'chunk_%04d.rds' "${rep_offset}")"

echo "  simulations:  $((rep_offset + 1))-$((rep_offset + reps_here)) of ${NUM_REPLICATES}"
echo "  seed:        ${SEED}  (must match step 04)"
echo "  output:      ${chunk_out}"
echo

# Resume: a completed chunk is left alone, so re-running the array after a partial failure picks up
# only what is missing. Same contract as step 04.
if [[ -s "${chunk_out}" ]]; then
  echo ">>> chunk already present, skipping"
  ls -la "${chunk_out}" | sed 's/^/  /'
  report_usage
  exit 0
fi

pixi_run Rscript src/fit_null_models.R \
  --sim-input "${PREPARED_DIR}/sim_input.rds" \
  --sceptre-template "${PREPARED_DIR}/sceptre_template.rds" \
  --grna-targets "${PREPARED_DIR}/grna_targets.tsv" \
  --reps "${reps_here}" \
  --rep-offset "${rep_offset}" \
  --seed "${SEED}" \
  --out "${chunk_out}"

require_file "${chunk_out}" "null-model chunk"
ls -la "${chunk_out}" | sed 's/^/  /'

report_usage
