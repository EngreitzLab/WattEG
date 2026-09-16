#!/bin/bash
#SBATCH --job-name=power_01_prepare
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 1 -- turn the sceptre object into the small inputs the simulation needs.
#
# Runs once per sample. This is the step that makes everything downstream cheap: it emits a
# ~59 MB sim_input.rds in place of the 427 MB of objects every parallel task used to read.
#
# Memory: measured at 8.5 GB MaxRSS on day0_grna20_no_shuffle (job 38846800, 70s wall), against
# the 7.7 GB docs/status.md recorded for the same data. 16 GB is that plus headroom.
#
# Usage:  sbatch workflow/slurm_executor/01_prepare_sim_input.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 1: prepare_sim_input for sample ${SAMPLE}"

load_pixi
require_file "${SCEPTRE_OBJECT}" "sceptre object"
mkdir -p "${PREPARED_DIR}"

echo "input:   ${SCEPTRE_OBJECT}  ($(du -h "${SCEPTRE_OBJECT}" | cut -f1))"
echo "outdir:  ${PREPARED_DIR}"
echo

pixi_run Rscript src/prepare_sim_input.R \
  --sceptre-object "${SCEPTRE_OBJECT}" \
  --outdir "${PREPARED_DIR}"

## CHECKS =========================================================================================

echo
echo ">>> outputs"
for f in sim_input.rds sceptre_template.rds pairs.tsv grna_targets.tsv discovery_threshold.txt; do
  require_file "${PREPARED_DIR}/${f}" "expected output ${f}"
  printf '  %-24s %s\n' "${f}" "$(du -h "${PREPARED_DIR}/${f}" | cut -f1)"
done

echo
echo ">>> pair count"
# Minus one for the header row. The old pipeline produced 34,886 pairs for this sample, so a
# materially different number here means the two runs are not testing the same thing.
new_pairs=$(($(wc -l < "${PREPARED_DIR}/pairs.tsv") - 1))
echo "  QC-passing pairs (new):  ${new_pairs}"
old_power="${OLD_DIR}/power_analysis/power_analysis_results_es_0.15.tsv"
if [[ -s "${old_power}" ]]; then
  echo "  pairs in old output:     $(($(wc -l < "${old_power}") - 1))"
fi

echo
echo ">>> discovery p-value threshold"
# The threshold decides which simulated simulations count as detections, so a mismatch against the
# old pipeline would shift every power estimate. The old code took max(p_value) over the
# significant discoveries; the new code derives it from @discovery_result inside the object.
# Cross-check the two here, in awk, so it costs nothing and needs no R.
new_threshold=$(cat "${PREPARED_DIR}/discovery_threshold.txt")
echo "  derived from @discovery_result (new):  ${new_threshold}"

old_csv="${OLD_DIR}/sceptre_discovery_results.csv"
if [[ -s "${old_csv}" ]]; then
  old_threshold=$(awk -F',' '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        gsub(/^"|"$/, "", $i)
        if ($i == "p_value")     p = i
        if ($i == "significant") s = i
      }
      next
    }
    p && s {
      sig = $s; gsub(/^"|"$/, "", sig)
      if (sig == "TRUE" || sig == "True" || sig == "true") {
        if ($p + 0 > max) max = $p + 0
      }
    }
    END { if (max > 0) printf "%.10g", max; else print "none" }
  ' "${old_csv}")
  echo "  max p among significant (old CSV):     ${old_threshold}"
  echo
  echo "  These should agree. If they do not, do not proceed -- every downstream power"
  echo "  number depends on this threshold."
else
  echo "  (old discovery CSV not found at ${old_csv}; skipping cross-check)"
fi

report_usage
