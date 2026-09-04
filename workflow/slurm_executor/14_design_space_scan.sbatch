#!/bin/bash
#SBATCH --job-name=power_14_scan
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Steps 2 and 3 of paper/experiments.md: does a reduced sweep reproduce a full one?
#
# WHY THIS COSTS NO SIMULATION
#
# The question is whether three effect sizes at fewer replicates give the same per-pair curves as
# six at 100. That is answered by SUBSAMPLING a sweep that already exists, not by running smaller
# ones. Replicates are i.i.d. draws, so a random subset of 30 is distributed exactly like a fresh
# 30-replicate run -- and subsampling makes the comparison paired, on the same pairs and the same
# simulated data, which a rerun could not be. So this is minutes of arithmetic over the stored
# per-replicate output, not CPU-days of simulation.
#
# WHAT IT PRODUCES
#
#   power/<reps>_<seed>/power_es<es>.tsv   power recomputed from a subsample
#   curves/<grid>_<reps>_<seed>.tsv        the fitted one-parameter curve per pair
#
# and leaves the comparison against the full 100-replicate reference to the analysis that reads
# them. The reference itself is the sweep's own power_summary.tsv -- step 1 needs nothing built.
#
# THE DESIGN SPACE, and why it is a scan rather than one comparison
#
# Reporting that "3 effect sizes at 30 replicates works" invites the question "why 30?". A scan
# answers it with a curve: where the fit breaks is more informative than one point that happens to
# hold. Replicate counts and grids are both varied, and each (grid, reps) cell is repeated over
# several subsample seeds so the result carries a spread rather than a single draw.
#
# Grids are dataset-specific and matter more than the number of points -- day0 is poorly powered
# and 84% of its pairs transition between 5% and 25%, so a grid that brackets that transition is
# expected to beat wtc11's best grid here. That expectation is exactly what this measures.
#
# Usage:  sbatch workflow/slurm_executor/14_design_space_scan.sbatch
#         SCAN_SAMPLE=day2 sbatch workflow/slurm_executor/14_design_space_scan.sbatch

set -euo pipefail

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

# SCAN_SAMPLE, not SAMPLE. config.sh above already exports SAMPLE (=day0_grna20_no_shuffle, the
# sbatch runner's own dataset), so `${SAMPLE:-day0}` silently inherits that instead of defaulting --
# the first run of this script scanned a sweep that does not exist and died on a missing threshold
# file. Same trap as the `scratch` shadowing in conf/sherlock.config: a local name colliding with an
# inherited one, reported as something else entirely. Do not rename these back.
SCAN_SAMPLE="${SCAN_SAMPLE:-day0}"
SWEEP="${SCAN_SWEEP:-/home/users/emattei/oak/projects/element-gene-power-analysis/power_sweep/${SCAN_SAMPLE}/${SCAN_SAMPLE}}"
OUT="${SCAN_OUT:-${SWEEP}/reduced_design}"

# All six, so any three or four of them can be selected as a fit grid and the rest held out.
ALL_ES=(0.05 0.1 0.15 0.2 0.25 0.5)

# Replicate counts. 100 is included as the no-subsampling control: at 100 the subsample IS the full
# set, so that row isolates the effect of the grid from the effect of the replicate count.
REPS=(${REPS_LIST:-10 20 30 50 100})

# Repeats per (grid, reps) cell. One draw would give a number with no error bar, and the spread is
# the point -- it is what says whether a design is reliably good or got lucky.
SEEDS=(${SEEDS_LIST:-1 2 3})

log_header "Step 14: reduced-design scan for ${SCAN_SAMPLE}"
load_pixi

mkdir -p "${OUT}/power" "${OUT}/curves"
echo "  sweep:  ${SWEEP}"
echo "  out:    ${OUT}"

THRESHOLD="${SWEEP}/prepared/discovery_threshold.txt"
if [[ ! -s "${THRESHOLD}" ]]; then
  echo "ERROR: no discovery threshold at ${THRESHOLD}" >&2
  exit 1
fi

## ---- recompute power at each (effect size, replicate count, seed) ------------------------------
#
# 100 replicates needs no subsampling and no repeats: every seed would select the same full set, so
# it is computed once and reused, which also keeps it byte-identical to the sweep's own output.
for reps in "${REPS[@]}"; do
  for seed in "${SEEDS[@]}"; do
    if [[ "${reps}" -eq 100 && "${seed}" -ne "${SEEDS[0]}" ]]; then continue; fi
    dir="${OUT}/power/${reps}_${seed}"
    mkdir -p "${dir}"
    for es in "${ALL_ES[@]}"; do
      target="${dir}/power_es${es}.tsv"
      [[ -s "${target}" ]] && continue          # resumable: re-running skips finished cells
      sub_args=(--subsample-reps "${reps}" --subsample-seed "${seed}")
      [[ "${reps}" -eq 100 ]] && sub_args=()
      echo ">>> power: es=${es} reps=${reps} seed=${seed}"
      pixi_run Rscript src/compute_power.R \
          --simulations "${SWEEP}/sim/es${es}" \
          --threshold-file "${THRESHOLD}" \
          --out "${target}" \
          "${sub_args[@]}"
    done
  done
done

## ---- fit the curve on each grid, so the held-out effect sizes can be predicted -----------------
#
# Four grids: wtc11's best (0.05/0.1/0.25) as the cross-dataset comparison, two that bracket day0's
# own transition, and a four-point grid to separate "more points" from "better placed points".
GRIDS=(
  "0.05,0.1,0.25"
  "0.05,0.15,0.25"
  "0.1,0.2,0.25"
  "0.05,0.1,0.2,0.25"
)

for grid in "${GRIDS[@]}"; do
  for reps in "${REPS[@]}"; do
    for seed in "${SEEDS[@]}"; do
      if [[ "${reps}" -eq 100 && "${seed}" -ne "${SEEDS[0]}" ]]; then continue; fi
      dir="${OUT}/power/${reps}_${seed}"
      name="$(echo "${grid}" | tr ',' '-')_${reps}_${seed}"
      target="${OUT}/curves/${name}.tsv"
      [[ -s "${target}" ]] && continue
      files=""
      for es in ${grid//,/ }; do
        files="${files}${files:+,}${dir}/power_es${es}.tsv"
      done
      echo ">>> curve: grid=${grid} reps=${reps} seed=${seed}"
      pixi_run Rscript src/fit_power_curve.R \
          --power "${files}" \
          --threshold-file "${THRESHOLD}" \
          --out "${target}"
    done
  done
done

echo
echo "Done. Power tables in ${OUT}/power, fitted curves in ${OUT}/curves."
echo "The reference to compare against is ${SWEEP}/power_summary.tsv (100 replicates, six points)."
