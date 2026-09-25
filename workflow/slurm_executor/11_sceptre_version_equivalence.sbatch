#!/bin/bash
#SBATCH --job-name=power_11_sceptre_version
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 11 -- does moving the sceptre pin from v0.10.3 to 0.99.0 change the numbers?
#
# The pin moved from v0.10.3 (9faa4373) to 0.99.0 (3ba046b8). 0.99.0's NEWS describes it as a
# Bioconductor-prep release with no statistical change, and reading the diff agrees: the one RNG
# difference on the association path is in `partition_response_ids()`, where `set.seed(4)` became
# `withr::with_seed(4, ...)` -- and that whole block sits inside `if (parallel)`, which this
# pipeline never enters because it calls `run_discovery_analysis(parallel = FALSE)`.
#
# That is a reading of the source, not a measurement. This job measures it, because the old-vs-new
# comparison (step 07) compares against results produced under a different sceptre, and if the
# version alone moves p-values then that comparison is confounded and has to be reinterpreted.
#
# The comparison is exactly paired. Seeds derive from (seed, target, rep, effect_size), which does
# not involve the sceptre version, so the simulated count matrices are bit-identical between the two
# runs. Both runs use --no-grna-precomp-reuse, so the patch is out of the picture and the only
# difference is the sceptre build.
#
# For that to be possible, stock v0.10.3 is installed into its own library: our patch does not
# apply to it, and unpatched sceptre has no `grna_precomputations` argument. Passing an empty
# SCEPTRE_PATCH_DIR is what asks install_sceptre.R for a stock build; run_power_simulation.R only
# passes the argument when it has a cache, so with reuse off it runs against stock sceptre.
#
# Usage:  sbatch workflow/slurm_executor/11_sceptre_version_equivalence.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

OLD_SHA="9faa4373d1b72164b630fa03c597e52e668a4e03"   # v0.10.3, the previous pin
OLD_LIB="${SCRATCH}/sceptre_stock_v0.10.3"
EMPTY_PATCH_DIR="${SCRATCH}/sceptre_no_patches"
CHECK_REPS=3
EFFECT_SIZE="${EFFECT_SIZES[0]}"
CHECK_DIR="${OUTDIR}/sceptre_version_check"

log_header "Step 11: sceptre v0.10.3 vs 0.99.0, ${CHECK_REPS} reps at effect size ${EFFECT_SIZE}"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"        "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"
require_dir_nonempty "${SPLITS_DIR}"                "splits directory from step 02"

mkdir -p "${CHECK_DIR}" "${OLD_LIB}" "${EMPTY_PATCH_DIR}"
split=$(ls "${SPLITS_DIR}"/split_*.tsv | sed -n 1p)
echo "split:      $(basename "${split}") ($(($(wc -l < "${split}") - 1)) pairs)"
echo "old lib:    ${OLD_LIB}"

# Expanded below as ${null_args[*]-}, not ${null_args[*]}. config.sh sets `set -u` and this cluster
# runs bash 4.2, where expanding an empty array under `set -u` is an "unbound variable" error.
null_args=()
if [[ -s "${NULL_MODELS}" ]]; then
  null_args=(--null-precomputations "${NULL_MODELS}")
fi

## INSTALL STOCK v0.10.3 ALONGSIDE =================================================================

# Inline env assignments rather than --env forwarding: container_run runs with --cleanenv and a
# fixed --env list, and pixi sets R_LIBS_USER itself from [activation.env], so overriding it has to
# happen inside the command.
echo
echo ">>> installing stock sceptre v0.10.3 into ${OLD_LIB}"
container_run pixi run --frozen --environment build bash -c "
  R_LIBS_USER='${OLD_LIB}' \
  SCEPTRE_SHA='${OLD_SHA}' \
  SCEPTRE_REF='v0.10.3' \
  SCEPTRE_PATCH_DIR='${EMPTY_PATCH_DIR}' \
  Rscript src/install_sceptre.R
"

echo
echo ">>> versions in play"
container_run pixi run --frozen bash -c "
  echo -n '  main library:  '
  Rscript -e 'cat(as.character(utils::packageVersion(\"sceptre\")), \"/\", utils::packageDescription(\"sceptre\")\$LocalPatches, \"\n\")'
  echo -n '  side library:  '
  R_LIBS_USER='${OLD_LIB}' Rscript -e 'cat(as.character(utils::packageVersion(\"sceptre\")), \"/\", utils::packageDescription(\"sceptre\")\$LocalPatches, \"\n\")'
"

## RUN THE SAME SPLIT UNDER EACH VERSION ===========================================================

run_under() {
  local label="$1" lib="$2"
  local out="${CHECK_DIR}/${label}_es${EFFECT_SIZE}.tsv"
  echo
  echo ">>> ${label}"
  local started=${SECONDS}
  container_run pixi run --frozen bash -c "
    ${lib:+R_LIBS_USER='${lib}'} Rscript src/run_power_simulation.R \
      --sim-input '${PREPARED_DIR}/sim_input.rds' \
      --sceptre-template '${PREPARED_DIR}/sceptre_template.rds' \
      --pairs '${split}' \
      --grna-targets '${PREPARED_DIR}/grna_targets.tsv' \
      ${null_args[*]-} \
      --effect-size '${EFFECT_SIZE}' \
      --reps '${CHECK_REPS}' \
      --rep-offset 0 \
      --guide-spread-c '${GUIDE_SPREAD_C}' \
      --seed '${SEED}' \
      --out '${out}' \
      --no-grna-precomp-reuse
  "
  echo "    elapsed: $((SECONDS - started))s -> ${out}"
}

run_under "v0.10.3_stock" "${OLD_LIB}"
run_under "v0.99.0_patched" ""

## COMPARE =========================================================================================

old_out="${CHECK_DIR}/v0.10.3_stock_es${EFFECT_SIZE}.tsv"
new_out="${CHECK_DIR}/v0.99.0_patched_es${EFFECT_SIZE}.tsv"
require_file "${old_out}" "v0.10.3 output"
require_file "${new_out}" "0.99.0 output"

echo
echo ">>> comparison"
if cmp -s "${old_out}" "${new_out}"; then
  echo "    IDENTICAL: the sceptre version bump does not change any number."
  echo "    Step 07's old-vs-new comparison is therefore not confounded by the bump."
  exit_code=0
else
  echo "    DIFFERENT. The version bump moves results; step 07 must account for it."
  echo
  echo "    row counts: $(($(wc -l < "${old_out}") - 1)) vs $(($(wc -l < "${new_out}") - 1))"
  echo "    p_value summary of the paired difference:"
  container_run pixi run --frozen Rscript -e "
    a <- read.delim('${old_out}'); b <- read.delim('${new_out}')
    k <- c('grna_target', 'response_id', 'rep')
    k <- k[k %in% names(a)]
    m <- merge(a, b, by = k, suffixes = c('_old', '_new'))
    cat('      paired rows:', nrow(m), '\n')
    d <- abs(m\$p_value_old - m\$p_value_new)
    cat('      max |dp|:   ', formatC(max(d, na.rm = TRUE), format = 'e', digits = 3), '\n')
    cat('      spearman:   ', formatC(suppressWarnings(cor(m\$p_value_old, m\$p_value_new,
        method = 'spearman', use = 'complete.obs')), digits = 4), '\n')
  " | sed 's/^/  /'
  exit_code=1
fi

echo
echo "Note: this isolates the sceptre version only -- both runs had the gRNA precomputation reuse"
echo "off. Step 10 isolates the patch. Step 07 covers the whole refactor."
report_usage
exit "${exit_code}"
