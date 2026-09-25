#!/bin/bash
#SBATCH --job-name=power_04_sim
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=logs/refactor/%x-%A_%a.out
#SBATCH --error=logs/refactor/%x-%A_%a.err
#
# Step 4 -- the real simulation. One array task per (split, effect size).
#
# SUBMIT THIS WITH 04_submit.sh, NOT WITH sbatch DIRECTLY.
#
# The array range depends on N_SPLITS x the number of effect sizes, and #SBATCH directives are
# read by sbatch before the script runs, so they cannot reference a variable sourced at runtime.
# 04_submit.sh computes the range and passes it as --array, which overrides any directive here.
# There is deliberately no --array directive, so a bare `sbatch` fails loudly rather than
# silently running a stale range.
#
# WHY (split x effect size) AND NOT A LOOP OVER EFFECT SIZES
#
# Looping over effect sizes inside each task would share the per-task setup -- the 18 MB
# sim_input.rds and 50 MB sceptre_template.rds reads, plus container and R startup -- across all
# of them. But that setup is ~1 second against a task measured in tens of minutes, while the loop
# multiplies each task's duration by the number of effect sizes and divides parallelism by the
# same factor. With thousands of idle cores on `owners`, trading a little repeated I/O for full
# parallelism is much the better deal: a 10-effect-size sweep goes from ~11 hours to ~70 minutes.
#
# --cpus-per-task=1 is deliberate: sceptre is called with parallel = FALSE and all the
# parallelism comes from the array fan-out. Asking for more cores would idle them.
#
# --time=4h against a measured ~49 min at 1000 splits. The headroom is not superstition: the
# per-target cost varies 4x with pair count (10.3s to 41.4s for 2 reps, measured in the step 3
# smoke test), so a split that happens to collect several large targets runs well over the
# median. Raise it if you lower N_SPLITS -- at 290 splits a task is ~2.8 hours.
#
# Resubmitting: completed tasks are skipped (see the resume check below), so re-running
# 04_submit.sh after a partial failure is safe and picks up only what is missing. To redo
# specific tasks:
#   sbatch --array=17,42,103 workflow/slurm_executor/04_run_power_simulation.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: no SLURM_ARRAY_TASK_ID -- this script must be submitted as an array." >&2
  echo "       Use: workflow/slurm_executor/04_submit.sh" >&2
  exit 1
fi

## DECOMPOSE THE TASK ID INTO (split, effect size) ================================================
#
# Split-major within each effect size: tasks 1..N_SPLITS are effect size 1, N_SPLITS+1..2N are
# effect size 2, and so on. So a single effect size can be resubmitted by array range.

if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${TOTAL_SIM_TASKS}" ]]; then
  echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} exceeds ${TOTAL_SIM_TASKS} = ${N_SPLITS} splits x ${N_EFFECT_SIZES} effect sizes." >&2
  echo "       The array range and config.sh have diverged; resubmit with 04_submit.sh." >&2
  exit 1
fi

idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
split_idx=$(( idx % N_SPLITS + 1 ))
es_idx=$(( idx / N_SPLITS ))
ES="${EFFECT_SIZES[${es_idx}]}"

log_header "Step 4: split ${split_idx}/${N_SPLITS}, effect size ${ES}"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"        "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"
require_dir_nonempty "${SPLITS_DIR}"                "splits directory from step 02"

# The per-simulation null models (NULL_MODELS comes from config.sh -- do not redefine it here, or it
# will drift from the path step 02c writes). Without these, sceptre either refits inside every call
# (4.3x) or -- if the template still carries the inherited real-data cache -- silently tests
# simulated counts against real-data coefficients, which understates power. run_power_simulation.R
# refuses the latter.
require_file "${NULL_MODELS}" "null_precomputations.rds from step 02b"

n_available=$(find "${SPLITS_DIR}" -name 'split_*.tsv' | wc -l)
if [[ "${n_available}" -ne "${N_SPLITS}" ]]; then
  echo "ERROR: ${n_available} split files present but N_SPLITS=${N_SPLITS} in config.sh." >&2
  echo "       Re-run step 02 after changing N_SPLITS." >&2
  exit 1
fi

# By position, not by formatting the index: split_pairs.R pads filenames to the width of
# --n-splits, so 1000 splits give split_0001.tsv but 280 would give split_001.tsv.
SPLIT=$(ls "${SPLITS_DIR}"/split_*.tsv | sed -n "${split_idx}p")
SPLIT_NAME=$(basename "${SPLIT}" .tsv)
n_pairs=$(($(wc -l < "${SPLIT}") - 1))
echo "split file:  ${SPLIT}"
echo "pairs:       ${n_pairs}"

## SIMULATE =======================================================================================

out_dir="${SIM_DIR}/es${ES}"
mkdir -p "${out_dir}"
out="${out_dir}/${SPLIT_NAME}.tsv"
provenance="${out}.provenance"

# What this output was produced from. Row count cannot distinguish two runs that differ only in
# their null model -- see the long note in config.sh -- so the resume check below requires this
# string to match as well, not just the line count.
NULL_MODELS_MD5=$(md5sum "${NULL_MODELS}" | cut -d' ' -f1)
want_provenance="config=${NULL_MODEL_CONFIG} null_models_md5=${NULL_MODELS_MD5} reps=${NUM_REPLICATES} seed=${SEED} guide_spread_c=${GUIDE_SPREAD_C} estimand=${ESTIMAND} es=${ES}"

# Resume support. A task that already produced a complete output is skipped, so the array can be
# resubmitted wholesale after a partial failure without redoing finished work.
#
# "Complete" means both:
#   * the right row count -- reps x pairs data rows plus a header. A task killed mid-write, or
#     preempted on `owners`, leaves a short file, which is caught here rather than silently
#     under-counting simulations in compute_power.R.
#   * matching provenance. An output produced under a different null-model configuration, a
#     different null-model bundle, or different seed/reps is NOT resumable, and re-running it is
#     the whole point of the resubmission. A missing sidecar means the file predates this check,
#     so it cannot be vouched for and is redone.
if [[ -s "${out}" ]]; then
  expected=$(( NUM_REPLICATES * n_pairs + 1 ))
  actual=$(wc -l < "${out}")
  have_provenance=""
  [[ -s "${provenance}" ]] && have_provenance=$(cat "${provenance}")

  if [[ "${actual}" -eq "${expected}" && "${have_provenance}" == "${want_provenance}" ]]; then
    echo
    echo ">>> already complete (${actual} lines, provenance matches), skipping"
    report_usage
    exit 0
  fi

  echo
  if [[ "${actual}" -ne "${expected}" ]]; then
    echo ">>> existing output is incomplete (${actual} of ${expected} lines), redoing"
  elif [[ -z "${have_provenance}" ]]; then
    echo ">>> existing output has no provenance sidecar, redoing"
    echo "      want: ${want_provenance}"
  else
    echo ">>> existing output was produced differently, redoing"
    echo "      have: ${have_provenance}"
    echo "      want: ${want_provenance}"
  fi
  rm -f "${out}" "${provenance}"
fi

echo
echo ">>> ${NUM_REPLICATES} reps at effect size ${ES} -> ${out}"
pixi_run Rscript src/run_power_simulation.R \
  --sim-input "${PREPARED_DIR}/sim_input.rds" \
  --sceptre-template "${PREPARED_DIR}/sceptre_template.rds" \
  --pairs "${SPLIT}" \
  --grna-targets "${PREPARED_DIR}/grna_targets.tsv" \
  --null-precomputations "${NULL_MODELS}" \
  --effect-size "${ES}" \
  --reps "${NUM_REPLICATES}" \
  --rep-offset 0 \
  --guide-spread-c "${GUIDE_SPREAD_C}" \
  --estimand "${ESTIMAND}" \
  --seed "${SEED}" \
  --out "${out}"

require_file "${out}" "simulation output"
echo "    wrote $(($(wc -l < "${out}") - 1)) rows"

# Written only after the output is confirmed present, so an interrupted task never leaves a
# sidecar vouching for a file that does not exist or is short.
printf '%s\n' "${want_provenance}" > "${provenance}"
echo "    provenance: ${want_provenance}"

report_usage
