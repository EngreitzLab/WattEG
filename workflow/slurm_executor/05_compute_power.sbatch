#!/bin/bash
#SBATCH --job-name=power_05_compute
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# 8 GB, not the 32 GB this asked for until 2026-08-14: measured 2.1 GB peak and 42 seconds over
# 1,000 simulation files at 34,886 pairs x 100 simulations (job 38978285). The old figure was never
# measured and queued the job behind almost everything on a shared partition.
#
# Step 5 -- collapse the per-simulation results into one power estimate per pair.
#
#   power = mean(p_value < threshold & log_2_fold_change < 0)
#
# reported with a Wilson score interval, one row per element-gene pair.
#
# Memory: compute_power.R rbinds every simulation row into one data frame. At ~35,000 pairs x 100
# simulations that is ~3.5M rows, which fits comfortably in 32 GB but not in the 8 GB the
# simulation tasks use.
#
# Usage:  sbatch workflow/slurm_executor/05_compute_power.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 5: compute power per pair"

load_pixi
require_file "${PREPARED_DIR}/discovery_threshold.txt" "discovery_threshold.txt from step 01"
mkdir -p "${POWER_DIR}"

for es in "${EFFECT_SIZES[@]}"; do
  sim_dir="${SIM_DIR}/es${es}"
  require_dir_nonempty "${sim_dir}" "simulation outputs for effect size ${es}"

  # Every split must have produced an output, or the power estimates would silently cover only
  # the pairs that happened to finish.
  n_sim=$(find "${sim_dir}" -name 'split_*.tsv' | wc -l)
  echo
  echo ">>> effect size ${es}: ${n_sim} of ${N_SPLITS} split outputs present"
  if [[ "${n_sim}" -ne "${N_SPLITS}" ]]; then
    echo "ERROR: ${n_sim} simulation outputs but ${N_SPLITS} splits were submitted." >&2
    echo "       Missing splits:" >&2
    for f in "${SPLITS_DIR}"/split_*.tsv; do
      [[ -s "${sim_dir}/$(basename "${f}")" ]] || echo "         $(basename "${f}")" >&2
    done
    echo "       Resubmit those array tasks (step 04) before computing power." >&2
    exit 1
  fi

  # --simulations takes a comma-separated list. Built with paste rather than by interpolating a
  # glob, so the ordering is stable and the argument is a single token.
  sim_list=$(ls "${sim_dir}"/split_*.tsv | paste -sd, -)
  out="${POWER_DIR}/power_es${es}.tsv"

  pixi_run Rscript src/compute_power.R \
    --simulations "${sim_list}" \
    --threshold-file "${PREPARED_DIR}/discovery_threshold.txt" \
    --conf-level 0.95 \
    --out "${out}"

  require_file "${out}" "power output for es ${es}"

  echo
  echo "  wrote $(($(wc -l < "${out}") - 1)) pairs to ${out}"

  # n_reps below NUM_REPLICATES means simulation rows were dropped for a missing p_value or fold
  # change. A handful is unremarkable; a systematic shortfall is not.
  awk -F'\t' -v expected="${NUM_REPLICATES}" '
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    { n++; if ($col["n_reps"] + 0 < expected) short++ }
    END { printf "  pairs with fewer than %d usable simulations: %d of %d\n", expected, short + 0, n }
  ' "${out}"
done

report_usage
