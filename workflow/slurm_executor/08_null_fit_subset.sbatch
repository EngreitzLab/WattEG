#!/bin/bash
#SBATCH --job-name=power_08_null_fit_subset
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=logs/refactor/%x-%A_%a.out
#SBATCH --error=logs/refactor/%x-%A_%a.err
#
# Step 8 -- re-run a subset of splits under null_fit, to measure what the inherited cache cost.
#
# WHY A SUBSET IS ENOUGH
#
# The comparison against step 4's output is *exactly paired*. Seeds derive from
# (seed, target, rep, effect_size), which does not involve the null model, so the simulated count
# matrices here are bit-identical to the ones step 4 drew. Both configurations also skip the GLM fit
# -- step 4 from the inherited cache, this one from the bundle -- so RNG consumption matches. Every
# difference in power is therefore attributable to the null model, with no Monte-Carlo noise between
# the two runs.
#
# That is why NULL_FIT_SUBSET_SPLITS splits at 100 simulations settle a question that 53 pairs at 5
# simulations could only raise, at ~1/50th the cost of a second full sweep.
#
# Only the first effect size is run: the bias is a property of the null model, not of the knockdown.
#
# Usage:  sbatch --array=1-20 workflow/slurm_executor/08_null_fit_subset.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: no SLURM_ARRAY_TASK_ID -- submit as an array, e.g." >&2
  echo "       sbatch --array=1-${NULL_FIT_SUBSET_SPLITS} workflow/slurm_executor/08_null_fit_subset.sbatch" >&2
  exit 1
fi

if (( SLURM_ARRAY_TASK_ID > NULL_FIT_SUBSET_SPLITS )); then
  echo "ERROR: task ${SLURM_ARRAY_TASK_ID} exceeds NULL_FIT_SUBSET_SPLITS=${NULL_FIT_SUBSET_SPLITS}." >&2
  exit 1
fi

ES="${EFFECT_SIZES[0]}"

# Splits 1..NULL_FIT_SUBSET_SPLITS, selected by position exactly as step 04 does -- split_pairs.R
# pads the filename to the width of --n-splits, so the index cannot be formatted directly.
SPLIT="$(ls "${SPLITS_DIR}"/split_*.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")"

log_header "Step 8: null_fit subset, split ${SLURM_ARRAY_TASK_ID}/${NULL_FIT_SUBSET_SPLITS}, effect size ${ES}"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"        "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"
require_file "${NULL_MODELS}"                       "null_precomputations.rds from step 02c"
require_file "${SPLIT}"                             "split file"

out_dir="${NULL_FIT_SUBSET_DIR}/es${ES}"
mkdir -p "${out_dir}"
out="${out_dir}/$(basename "${SPLIT}")"

echo "  split file:  ${SPLIT}"
echo "  pairs:       $(($(wc -l < "${SPLIT}") - 1))"
echo "  null models: ${NULL_MODELS}"
echo "  output:      ${out}"
echo

# Same resume contract as step 04.
if [[ -s "${out}" ]]; then
  echo ">>> output already present, skipping"
  report_usage
  exit 0
fi

pixi_run Rscript src/run_power_simulation.R \
  --sim-input "${PREPARED_DIR}/sim_input.rds" \
  --sceptre-template "${PREPARED_DIR}/sceptre_template.rds" \
  --pairs "${SPLIT}" \
  --grna-targets "${PREPARED_DIR}/grna_targets.tsv" \
  --null-precomputations "${NULL_MODELS}" \
  --effect-size "${ES}" \
  --reps "${NUM_REPLICATES}" \
  --rep-offset 0 \
  --guide-sd "${GUIDE_SD}" \
  --seed "${SEED}" \
  --out "${out}"

require_file "${out}" "null_fit simulation output"

report_usage
