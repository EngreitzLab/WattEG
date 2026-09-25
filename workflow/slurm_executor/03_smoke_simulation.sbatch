#!/bin/bash
#SBATCH --job-name=power_03_smoke
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 3 -- smoke test of the simulation, before committing ~300 CPU-hours to step 4.
#
# Three splits at 2 replicates. This exercises every part the full array does -- sceptre loads,
# the emptied template deserialises, the perturbation matrices subset, counts are drawn, sceptre
# is called, a TSV lands -- but finishes in minutes.
#
# It also gives the first real per-target timing on this sample, which is what step 4's --time
# request should be based on.
#
# Usage:  sbatch workflow/slurm_executor/03_smoke_simulation.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

SMOKE_SPLITS=3
SMOKE_REPS=2
EFFECT_SIZE="${EFFECT_SIZES[0]}"

log_header "Step 3: smoke test -- ${SMOKE_SPLITS} splits x ${SMOKE_REPS} reps at effect size ${EFFECT_SIZE}"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"       "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"
require_dir_nonempty "${SPLITS_DIR}"                "splits directory from step 02"

mkdir -p "${SMOKE_DIR}"

# Collect the outputs this run writes, and check only those. Globbing ${SMOKE_DIR} instead would pick
# up leftovers from an earlier run, and because split_pairs.R pads the split index to the width of
# --n-splits, a previous run at a different N_SPLITS leaves *differently named* files that a glob
# cannot tell apart: 280 splits give smoke_split_001_*.tsv, 1000 give smoke_split_0001_*.tsv. That is
# not hypothetical -- it made the simulation-count check below fail against 250-row files from a
# 280-split run while every file this run wrote was correct. Step 02 avoids the same trap by
# `rm -rf`-ing its output directory; this keeps prior smoke output instead, since it is cheap to keep
# and occasionally useful to compare against.
smoke_outputs=()

for i in $(seq 1 "${SMOKE_SPLITS}"); do
  # Selected by position rather than by formatting the index: split_pairs.R pads the filename to
  # the width of --n-splits, so 280 splits give split_001.tsv but 1000 would give split_0001.tsv.
  split=$(ls "${SPLITS_DIR}"/split_*.tsv | sed -n "${i}p")
  out="${SMOKE_DIR}/smoke_$(basename "${split}" .tsv)_es${EFFECT_SIZE}.tsv"

  echo
  echo ">>> smoke split ${i}/${SMOKE_SPLITS}: $(basename "${split}")"
  pixi_run Rscript src/run_power_simulation.R \
    --sim-input "${PREPARED_DIR}/sim_input.rds" \
    --sceptre-template "${PREPARED_DIR}/sceptre_template.rds" \
    --pairs "${split}" \
    --grna-targets "${PREPARED_DIR}/grna_targets.tsv" \
    --effect-size "${EFFECT_SIZE}" \
    --reps "${SMOKE_REPS}" \
    --rep-offset 0 \
    --guide-spread-c "${GUIDE_SPREAD_C}" \
    --seed "${SEED}" \
    --out "${out}"
  smoke_outputs+=("${out}")
done

## CHECKS =========================================================================================

echo
echo ">>> smoke outputs"
for f in "${smoke_outputs[@]}"; do
  require_file "${f}" "smoke output"
  printf '  %-44s %6d rows\n' "$(basename "${f}")" "$(($(wc -l < "${f}") - 1))"
done

first="${smoke_outputs[0]}"

echo
echo ">>> columns"
head -1 "${first}" | tr '\t' '\n' | sed 's/^/    /'

echo
echo ">>> first rows"
head -3 "${first}" | cut -c1-200 | sed 's/^/    /'

echo
echo ">>> sanity"
# Each pair must appear exactly SMOKE_REPS times, and p-values must be usable -- an all-NA
# p_value column is the signature of the simulated matrix reaching sceptre in the wrong class,
# which fails silently rather than erroring (see the note in lib/simulate.R).
awk -F'\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    next
  }
  {
    n++
    key = $col["grna_target"] SUBSEP $col["response_id"]
    reps[key]++
    if ($col["p_value"] == "NA" || $col["p_value"] == "") na_p++
    if ($col["log_2_fold_change"] == "NA" || $col["log_2_fold_change"] == "") na_lfc++
  }
  END {
    for (k in reps) { npairs++; if (reps[k] != REPS) bad++ }
    printf "    rows                    %d\n", n
    printf "    distinct pairs          %d\n", npairs
    printf "    pairs with != %d reps    %d\n", REPS, bad + 0
    printf "    NA p_value              %d\n", na_p + 0
    printf "    NA log_2_fold_change    %d\n", na_lfc + 0
    if (bad > 0)          { print  "    ERROR: simulation count per pair is wrong"; exit 1 }
    if (na_p >= n)        { print  "    ERROR: every p_value is NA"; exit 1 }
  }
' REPS="${SMOKE_REPS}" "${smoke_outputs[@]}"

echo
echo "Smoke test passed."
echo
echo "Sizing step 4 from the per-target timings above is only valid if step 4 runs the same null-model"
echo "configuration as this run did. With --null-precomputations (null_fit) a call is far cheaper than"
echo "without it (the pipeline refits per call, ~4.3x); check which one this run used before using its"
echo "timings, or step 4's --time will be over-requested by roughly that factor."

report_usage
