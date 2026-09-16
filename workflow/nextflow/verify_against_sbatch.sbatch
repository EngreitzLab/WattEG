#!/bin/bash
#SBATCH --job-name=power_nf_verify
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/nextflow/%x-%A.out
#SBATCH --error=logs/nextflow/%x-%A.err
#
# Check the Nextflow runner against the sbatch runner, output by output.
#
# The two runners exist side by side on purpose: the sbatch scripts produced every measurement in
# docs/status.md, so they are the reference until Nextflow reproduces them. This is the check that
# decides when that has happened. Run it after a Nextflow run and it says, per artefact, whether the
# two agree.
#
# Byte-identity is the bar for everything except the null-model bundle, which the two runners may
# cover to different depths (a test run fits 2 simulations, a production run 100). For that one the
# comparison is over the shared simulations, and it should still be exact -- null models are seeded
# from (seed, rep) alone, so simulation 1 is the same object either way.
#
# Usage:
#   sbatch workflow/nextflow/verify_against_sbatch.sbatch [nextflow_outdir]
#
# Defaults to the test tree; pass results/nextflow for a production run.

set -uo pipefail

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
REPO_ROOT="$(cd "$(dirname "${_cmd:-${BASH_SOURCE[0]}}")/../.." && pwd)"
cd "${REPO_ROOT}"
mkdir -p logs/nextflow

source workflow/slurm_executor/config.sh

NF_OUT="${1:-${REPO_ROOT}/results/nextflow_test}"
NF_DIR="${NF_OUT}/${SAMPLE}"

log_header "Verify: Nextflow (${NF_OUT}) against sbatch (${OUTDIR})"
load_pixi

status=0

echo
echo ">>> prepared/ -- byte comparison"
for f in pairs.tsv grna_targets.tsv discovery_threshold.txt sim_input.rds sceptre_template.rds; do
  a="${PREPARED_DIR}/${f}"
  b="${NF_DIR}/prepared/${f}"
  if [[ ! -s "${b}" ]]; then
    printf '  %-26s SKIP (not produced by this run)\n' "${f}"
  elif cmp -s "${a}" "${b}"; then
    printf '  %-26s IDENTICAL\n' "${f}"
  else
    printf '  %-26s DIFFERS\n' "${f}"
    status=1
  fi
done

echo
echo ">>> splits/ -- byte comparison of every split"
if [[ -d "${NF_DIR}/splits" ]]; then
  n_diff=$(diff -rq "${SPLITS_DIR}" "${NF_DIR}/splits" 2>&1 | wc -l)
  n_files=$(find "${NF_DIR}/splits" -name 'split_*.tsv' | wc -l)
  if [[ "${n_diff}" -eq 0 ]]; then
    echo "  ${n_files} splits, all IDENTICAL"
  else
    echo "  ${n_diff} difference(s) across ${n_files} splits"
    diff -rq "${SPLITS_DIR}" "${NF_DIR}/splits" 2>&1 | head -5 | sed 's/^/    /'
    status=1
  fi
else
  echo "  SKIP (not produced by this run)"
fi

echo
echo ">>> null_precomputations.rds -- shared simulations"
nf_null="${NF_DIR}/prepared/null_precomputations.rds"
if [[ -s "${nf_null}" && -s "${NULL_MODELS}" ]]; then
  pixi_run Rscript workflow/compare_null_models.R \
    --a "${NULL_MODELS}" \
    --b "${nf_null}" | sed 's/^/  /'
  [[ "${PIPESTATUS[0]}" -ne 0 ]] && status=1
else
  echo "  SKIP (one side missing)"
fi

echo
echo ">>> sim/ -- per-simulation rows, over the simulations both runners covered"
#
# NOT a byte comparison, and the reason matters.
#
# The simulation is seeded from (seed, target, rep, effect_size), so both runners draw the same
# random numbers. What they do not share is the arithmetic: tasks land on whatever node the
# scheduler gives them, and this cluster spans six CPU generations whose vectorised BLAS kernels
# differ in the last bit. Those last bits are then amplified by sceptre's skew-normal tail
# approximation, so a p-value of 1e-81 on one node can be 1e-80 on another while the fold change
# agrees to five significant figures past the fifteenth.
#
# Demanding byte-identity here would therefore fail on hardware, not on the pipeline. What must
# hold is the thing power is computed from: whether each simulation clears the discovery threshold.
# That is checked exactly, and any disagreement is a real defect.
ES="${EFFECT_SIZES[0]}"
nf_sim_dir="${NF_DIR}/sim/es${ES}"
sb_sim_dir="${SIM_DIR}/es${ES}"
THRESHOLD=$(cat "${PREPARED_DIR}/discovery_threshold.txt" 2>/dev/null || echo 0)

if [[ -d "${nf_sim_dir}" && -d "${sb_sim_dir}" ]]; then
  # The pair list is built first rather than counted inside the pipeline below: a `while` on the
  # right of a pipe runs in a subshell, so anything it counts is lost when the subshell exits.
  pairs_file="$(mktemp)"
  for f in "${nf_sim_dir}"/split_*.tsv; do
    [[ -s "${f}" ]] || continue
    base="$(basename "${f}")"
    counterpart="${sb_sim_dir}/${base%%_es*}.tsv"
    [[ -s "${counterpart}" ]] || continue
    printf '%s\t%s\n' "${counterpart}" "${f}" >> "${pairs_file}"
  done
  n_checked=$(wc -l < "${pairs_file}")

  while IFS=$'\t' read -r sb nf; do
    awk -F'\t' -v thr="${THRESHOLD}" -v label="$(basename "${nf}")" '
      FNR==NR { if (FNR>1) { k=$1"|"$12; p[k]=$6+0; fc[k]=$7+0 }; next }
      FNR>1 {
        k=$1"|"$12
        if (!(k in p)) next
        n++
        rp = (p[k]==0) ? 0 : ($6+0 - p[k]) / p[k];  if (rp<0) rp=-rp
        rf = (fc[k]==0) ? 0 : ($7+0 - fc[k]) / fc[k]; if (rf<0) rf=-rf
        if (rp > maxp) maxp = rp
        if (rf > maxf) maxf = rf
        if (((p[k] < thr) ? 1 : 0) != ((($6+0) < thr) ? 1 : 0)) flips++
      }
      END {
        printf "  %-34s %4d rows | max rel diff: p %.2g, fold_change %.2g | significance flips: %d\n",
               label, n, maxp, maxf, flips+0
        if (flips+0 > 0) exit 1
      }' "${sb}" "${nf}" || echo "__FLIP__"
  done < "${pairs_file}" | tee "${pairs_file}.out"

  if grep -q "__FLIP__" "${pairs_file}.out" 2>/dev/null; then
    echo "  FAIL -- a simulation changed significance between runners."
    status=1
  elif [[ "${n_checked}" -eq 0 ]]; then
    echo "  SKIP (no overlapping split outputs)"
  else
    echo "  ${n_checked} split(s) compared, all agreeing on every significance call"
  fi
  rm -f "${pairs_file}" "${pairs_file}.out"
else
  echo "  SKIP (one side missing)"
fi

echo
if [[ "${status}" -eq 0 ]]; then
  echo "VERIFY PASSED -- the two runners agree on everything compared."
else
  echo "VERIFY FAILED -- see above."
fi

report_usage
exit ${status}
