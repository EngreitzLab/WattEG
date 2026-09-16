#!/bin/bash
#SBATCH --job-name=power_10_grna_precomp
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 10 -- validate the gRNA precomputation reuse: exact, and worth having.
#
# patches/0001-reuse-grna-precomputation.patch lets the pipeline fit the gRNA null model once per
# target instead of once per (target, simulation). sceptre regresses perturbation status on the
# covariates and draws its synthetic assignments from the fitted probabilities; both inputs are
# fixed for the whole simulation loop and neither depends on the simulated counts, so every
# simulation after the first was refitting an identical model.
#
# The claim is that reuse is *exact*, not approximate. This job tests it the only way that settles
# it: run the same split twice under the same seed, once with reuse and once with
# --no-grna-precomp-reuse, and compare the outputs byte for byte. The fit consumes no RNG, so
# skipping it cannot shift the random stream and the p-values must match exactly. Anything less
# than an exact match means the two configurations are not analysing the same thing, and the reuse
# must not be used until that is understood.
#
# The wall-clock difference between the two runs is the first measurement of what the patch buys.
# It is a single-split measurement on one node, so treat it as an order of magnitude rather than a
# number to quote -- and note it is exactly the kind of figure this project has already had to
# retract once (see docs/status.md on the 59 % claim).
#
# Usage:  sbatch workflow/slurm_executor/10_grna_precomp_equivalence.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

# One split at a handful of replicates. Simulations are what the reuse acts on -- it removes
# (reps - 1) of every target's `reps` fits -- so the saving only becomes visible above 1, and 5 is
# enough to see it without paying for 100.
CHECK_REPS=5
EFFECT_SIZE="${EFFECT_SIZES[0]}"
CHECK_DIR="${OUTDIR}/grna_precomp_check"

log_header "Step 10: gRNA precomputation reuse -- exactness and cost, ${CHECK_REPS} reps at effect size ${EFFECT_SIZE}"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"        "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"
require_dir_nonempty "${SPLITS_DIR}"                "splits directory from step 02"

mkdir -p "${CHECK_DIR}"

split=$(ls "${SPLITS_DIR}"/split_*.tsv | sed -n 1p)
echo "split:          $(basename "${split}") ($(($(wc -l < "${split}") - 1)) pairs)"

# Match the configuration the production array runs. Without the null bundle sceptre refits the
# per-gene null model inside every call, which is ~4.3x the work and would swamp the gRNA fit we
# are trying to measure -- the comparison would still be valid, but the timing would be
# uninformative about the configuration that actually runs.
# Expanded below as ${null_args[@]+"${null_args[@]}"}, not "${null_args[@]}". config.sh sets `set -u`
# and this cluster runs bash 4.2, where expanding an empty array under `set -u` is an "unbound
# variable" error; bash only fixed that in 4.4.
null_args=()
if [[ -s "${NULL_MODELS}" ]]; then
  null_args=(--null-precomputations "${NULL_MODELS}")
  echo "null models:    ${NULL_MODELS}"
else
  echo "null models:    none on disk -- sceptre will refit the per-gene null model in every call,"
  echo "                so the timings below overstate the denominator. Exactness is unaffected."
fi

run_config() {
  local label="$1"; shift
  local out="${CHECK_DIR}/${label}_es${EFFECT_SIZE}.tsv"
  echo
  echo ">>> ${label}"
  local started=${SECONDS}
  pixi_run Rscript src/run_power_simulation.R \
    --sim-input "${PREPARED_DIR}/sim_input.rds" \
    --sceptre-template "${PREPARED_DIR}/sceptre_template.rds" \
    --pairs "${split}" \
    --grna-targets "${PREPARED_DIR}/grna_targets.tsv" \
    ${null_args[@]+"${null_args[@]}"} \
    --effect-size "${EFFECT_SIZE}" \
    --reps "${CHECK_REPS}" \
    --rep-offset 0 \
    --guide-sd "${GUIDE_SD}" \
    --seed "${SEED}" \
    --out "${out}" \
    "$@"
  echo "    elapsed: $((SECONDS - started))s -> ${out}"
  echo $((SECONDS - started)) > "${CHECK_DIR}/${label}.seconds"
}

## RUN BOTH CONFIGURATIONS =========================================================================

run_config refit --no-grna-precomp-reuse
run_config reuse

## EXACTNESS =======================================================================================

echo
echo ">>> exactness"
refit_out="${CHECK_DIR}/refit_es${EFFECT_SIZE}.tsv"
reuse_out="${CHECK_DIR}/reuse_es${EFFECT_SIZE}.tsv"
require_file "${refit_out}" "refit output"
require_file "${reuse_out}" "reuse output"

if cmp -s "${refit_out}" "${reuse_out}"; then
  echo "    IDENTICAL: the two outputs match byte for byte."
  exact=1
else
  echo "    DIFFERENT -- the reuse is not exact. First differing lines:"
  diff "${refit_out}" "${reuse_out}" | head -20 | sed 's/^/      /'
  exact=0
fi

## COST ============================================================================================

echo
echo ">>> cost"
refit_s=$(cat "${CHECK_DIR}/refit.seconds")
reuse_s=$(cat "${CHECK_DIR}/reuse.seconds")
printf '    refit every simulation   %5ds\n' "${refit_s}"
printf '    reuse per target        %5ds\n' "${reuse_s}"
if [[ "${reuse_s}" -gt 0 ]]; then
  awk -v a="${refit_s}" -v b="${reuse_s}" 'BEGIN {
    printf "    speedup                 %5.2fx  (%+.1f%% wall clock)\n", a / b, 100 * (b - a) / a
  }'
fi
echo
echo "    At ${CHECK_REPS} simulations the reuse removes $((CHECK_REPS - 1)) of every target'\''s ${CHECK_REPS} gRNA fits."
echo "    NUM_REPLICATES is ${NUM_REPLICATES} in production, so the fraction removed there is larger"
echo "    ($((NUM_REPLICATES - 1))/${NUM_REPLICATES}); the saving per call is what this measures, and it does not grow."

if [[ "${exact}" -ne 1 ]]; then
  echo
  echo "FAILED: reuse changed the results. Do not run the production array with it until this is" >&2
  echo "        explained. --no-grna-precomp-reuse restores the previous behaviour." >&2
  report_usage
  exit 1
fi

echo
echo "Passed: reuse is exact, and the timings above are the first measurement of its cost."
report_usage
