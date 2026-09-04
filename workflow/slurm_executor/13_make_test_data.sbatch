#!/bin/bash
#SBATCH --job-name=power_13_test_data
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Build the synthetic sceptre object, then run the whole pipeline on it.
#
# The generator validates its own output against read_sceptre_object()'s requirements, but that only
# proves the object is well formed. What matters is that the five pipeline steps run on it, because
# the point of the object is to let CI exercise the sceptre path without lab data. So this builds it
# and then runs prepare -> split -> null models -> simulate -> power end to end, at a scale that
# finishes in minutes.
#
# MOI and RESAMPLING select which kind of screen to synthesise. They exist because these are the
# two settings the pipeline INHERITS from the object rather than setting, so the only way to show it
# handles them is to build an object that carries them:
#
#   MOI=high RESAMPLING=default        crt          the default, and what every real dataset is
#   MOI=low  RESAMPLING=default        permutations low MOI resolves to permutations on its own
#   MOI=high RESAMPLING=permutations   permutations separates the mechanism from the MOI
#
# The third matters because the first two confound them: without it, a permutations failure could
# not be told apart from a low-MOI failure.
#
# Usage:  sbatch workflow/slurm_executor/13_make_test_data.sbatch
#         MOI=low sbatch workflow/slurm_executor/13_make_test_data.sbatch
#         MOI=high RESAMPLING=permutations sbatch workflow/slurm_executor/13_make_test_data.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

MOI="${MOI:-high}"
RESAMPLING="${RESAMPLING:-default}"

log_header "Step 13: synthetic test data (moi=${MOI}, resampling=${RESAMPLING}), and the pipeline on it"

load_pixi

# Variant-suffixed, so running one does not clobber another's object or outputs and the three can
# run concurrently.
VARIANT="${MOI}_${RESAMPLING}"
OBJ="${REPO_ROOT}/tests/data/sceptre_object_${VARIANT}.rds"
WORK="${SCRATCH}/synthetic_pipeline_check_${VARIANT}"
rm -rf "${WORK}"
mkdir -p "${WORK}" "$(dirname "${OBJ}")"

echo
echo ">>> generating"
pixi_run Rscript src/make_test_data.R --out "${OBJ}" \
    --moi "${MOI}" --resampling-mechanism "${RESAMPLING}"

echo
echo ">>> step 1: prepare"
pixi_run Rscript src/prepare_sim_input.R --sceptre-object "${OBJ}" --outdir "${WORK}/prepared"

# The point of the variants: check the pipeline REPORTS the screen it was actually handed. A run
# that silently analysed a low-MOI permutations screen as though it were high-MOI CRT would produce
# plausible numbers, so this asserts rather than prints.
echo
echo ">>> checking the reported analysis mode"
cat "${WORK}/prepared/analysis_mode.tsv"
expected_moi="${MOI}"
expected_mech="${RESAMPLING}"
if [[ "${expected_mech}" == "default" ]]; then
  expected_mech=$([[ "${MOI}" == "low" ]] && echo permutations || echo crt)
fi
got_moi=$(awk -F'\t' '$1=="moi"{print $2}' "${WORK}/prepared/analysis_mode.tsv")
got_mech=$(awk -F'\t' '$1=="resampling_mechanism"{print $2}' "${WORK}/prepared/analysis_mode.tsv")
if [[ "${got_moi}" != "${expected_moi}" || "${got_mech}" != "${expected_mech}" ]]; then
  echo "FAIL: built moi=${expected_moi}/${expected_mech} but the pipeline reports ${got_moi}/${got_mech}" >&2
  exit 1
fi
echo "  ok: reports ${got_moi} MOI, ${got_mech}"

# n_splits must not exceed the number of targets; the synthetic object has far fewer than a real
# one, which is exactly the kind of mismatch a CI run should catch rather than inherit.
n_targets=$(( $(cut -f1 "${WORK}/prepared/pairs.tsv" | tail -n +2 | sort -u | wc -l) ))
n_splits=$(( n_targets < 4 ? n_targets : 4 ))
echo "  targets: ${n_targets}, using ${n_splits} splits"

echo
echo ">>> step 2: split"
pixi_run Rscript src/split_pairs.R --pairs "${WORK}/prepared/pairs.tsv" \
  --n-splits "${n_splits}" --outdir "${WORK}/splits"

echo
echo ">>> step 2b/2c: null models (2 replicates)"
pixi_run Rscript src/fit_null_models.R \
  --sim-input "${WORK}/prepared/sim_input.rds" \
  --sceptre-template "${WORK}/prepared/sceptre_template.rds" \
  --grna-targets "${WORK}/prepared/grna_targets.tsv" \
  --reps 2 --rep-offset 0 --seed 1 --out "${WORK}/null_chunk.rds"
pixi_run Rscript src/merge_null_models.R \
  --inputs "${WORK}/null_chunk.rds" --reps 2 --out "${WORK}/null_precomputations.rds"

echo
echo ">>> step 4: simulate (2 replicates, first split)"
first_split=$(ls "${WORK}/splits"/split_*.tsv | head -1)
pixi_run Rscript src/run_power_simulation.R \
  --sim-input "${WORK}/prepared/sim_input.rds" \
  --sceptre-template "${WORK}/prepared/sceptre_template.rds" \
  --pairs "${first_split}" \
  --grna-targets "${WORK}/prepared/grna_targets.tsv" \
  --null-precomputations "${WORK}/null_precomputations.rds" \
  --effect-size 0.5 --reps 2 --rep-offset 0 --guide-sd 0.13 --seed 1 \
  --out "${WORK}/sim.tsv"

echo
echo ">>> step 5: power"
pixi_run Rscript src/compute_power.R \
  --simulations "${WORK}/sim.tsv" \
  --threshold-file "${WORK}/prepared/discovery_threshold.txt" \
  --conf-level 0.95 --out "${WORK}/power.tsv"

echo
echo ">>> step 6: summarize"
pixi_run Rscript src/summarize_power.R \
  --power "${WORK}/power.tsv" --power-threshold 0.8 --out "${WORK}/summary.tsv"

echo
echo ">>> results"
printf '  %-16s %s rows\n' "power.tsv" "$(( $(wc -l < "${WORK}/power.tsv") - 1 ))"
printf '  %-16s %s rows\n' "summary.tsv" "$(( $(wc -l < "${WORK}/summary.tsv") - 1 ))"
echo "  object: $(du -h "${OBJ}" | cut -f1)"

echo
echo "SYNTHETIC PIPELINE CHECK PASSED -- every step runs without lab data."

report_usage
