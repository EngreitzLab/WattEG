#!/usr/bin/env Rscript
# What did the inherited null-model cache cost us? Compare as_is against null_fit, pair by pair.
#
# WHY THIS COMPARISON IS EXACT
#
# Seeds derive from (seed, target, rep, effect_size) and not from anything about the null model, so
# the simulated count matrices behind the two runs are bit-identical. Both configurations skip the
# per-gene GLM fit -- as_is from the cache inherited from the real discovery analysis, null_fit from
# fit_null_models.R -- so RNG consumption matches as well. Every difference reported here is
# therefore attributable to the null model alone, with no Monte-Carlo noise between the two runs.
#
# That is what makes a subset of splits informative: this is a paired comparison of the same draws,
# not two independent estimates that have to be separated from sampling noise.
#
# WHAT IT REPORTS
#
# Power is `mean(p_value < threshold & log_2_fold_change < 0)` per pair -- the same definition
# compute_power.R uses, including the fold-change condition, because a significant p-value with a
# positive fold change is not a detected knockdown.
#
# Usage:
#   compare_null_fit.R --as-is 'sim/es0.15/split_*.tsv' --null-fit 'sim_null_fit/es0.15/split_*.tsv' \
#     --threshold-file prepared/discovery_threshold.txt --out comparison/null_fit_vs_as_is.tsv

suppressPackageStartupMessages({
  library(optparse)
})

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  lib <- if (length(file_arg) == 1) {
    file.path(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))), "bin", "lib")
  } else {
    "lib"
  }
  source(file.path(lib, "cli.R"))
})

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--as-is", type = "character", default = NULL, dest = "as_is",
              help = "Glob or comma-separated list of step 04 outputs (inherited cache)."),
  make_option("--null-fit", type = "character", default = NULL, dest = "null_fit",
              help = "Glob or comma-separated list of step 08 outputs (fitted null models)."),
  make_option("--threshold-file", type = "character", default = NULL, dest = "threshold_file",
              help = "discovery_threshold.txt from prepare_sim_input.R."),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per pair present in both runs.")
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Compare per-pair power between the inherited cache and fitted null models."
))
require_options(opts, c("as_is", "null_fit", "threshold_file", "out"))

expand_inputs <- function(spec, label) {
  paths <- if (grepl("[*?[]", spec)) sort(Sys.glob(spec)) else trimws(strsplit(spec, ",")[[1]])
  paths <- paths[nzchar(paths)]
  if (length(paths) == 0) stop("--", label, " matched no files: ", spec, call. = FALSE)
  paths
}

threshold <- suppressWarnings(as.numeric(readLines(opts$threshold_file, warn = FALSE)[1]))
if (is.na(threshold) || threshold <= 0 || threshold >= 1) {
  stop("Could not read a usable threshold from ", opts$threshold_file, ".", call. = FALSE)
}
log_step("Threshold: ", format(threshold, digits = 6))

## LOAD ============================================================================================

read_run <- function(spec, label) {
  paths <- expand_inputs(spec, label)
  log_step(label, ": ", length(paths), " file(s)")
  rows <- lapply(paths, function(p) {
    d <- read_tsv_file(p, required_columns = c("response_id", "grna_target", "p_value",
                                               "log_2_fold_change", "rep"))
    d[, c("response_id", "grna_target", "p_value", "log_2_fold_change", "rep")]
  })
  do.call(rbind, rows)
}

as_is <- read_run(opts$as_is, "as_is")
null_fit <- read_run(opts$null_fit, "null_fit")

# Per-pair power, exactly as compute_power.R defines it.
power_by_pair <- function(d) {
  called <- d$p_value < threshold & d$log_2_fold_change < 0
  key <- paste(d$grna_target, d$response_id, sep = "|")
  data.frame(pair = names(tapply(called, key, mean)),
             power = as.numeric(tapply(called, key, mean)),
             reps = as.integer(tapply(called, key, length)),
             row.names = NULL, stringsAsFactors = FALSE)
}

pa <- power_by_pair(as_is)
pn <- power_by_pair(null_fit)

merged <- merge(pa, pn, by = "pair", suffixes = c("_as_is", "_null_fit"))
if (nrow(merged) == 0) {
  stop("No pairs are present in both runs. Do the two globs cover the same splits?", call. = FALSE)
}
log_step("Pairs in both runs: ", nrow(merged), " (as_is ", nrow(pa), ", null_fit ", nrow(pn), ")")

# The comparison is only paired if both runs used the same replicates for a pair. Unequal counts mean
# one side is truncated, and the power values would differ for that reason rather than the null model.
unequal <- merged$reps_as_is != merged$reps_null_fit
if (any(unequal)) {
  stop(sum(unequal), " pair(s) have different replicate counts between the runs, e.g. ",
       merged$pair[which(unequal)[1]], " (", merged$reps_as_is[which(unequal)[1]], " vs ",
       merged$reps_null_fit[which(unequal)[1]], "). The comparison would not be paired.",
       call. = FALSE)
}

## PER-REPLICATE CALL FLIPS ========================================================================
#
# The aggregate p-value statistics from cache_experiment.R looked reassuring while the call-level
# ones did not, so this reports calls, not correlations.

both <- merge(
  transform(as_is, key = paste(grna_target, response_id, rep, sep = "|")),
  transform(null_fit, key = paste(grna_target, response_id, rep, sep = "|")),
  by = "key", suffixes = c("_as_is", "_null_fit")
)
call_as_is <- both$p_value_as_is < threshold & both$log_2_fold_change_as_is < 0
call_null_fit <- both$p_value_null_fit < threshold & both$log_2_fold_change_null_fit < 0
flips <- call_as_is != call_null_fit

cat("\n=========== per-replicate discovery calls ===========\n")
cat(sprintf("  observations compared:      %d\n", nrow(both)))
cat(sprintf("  called by as_is:            %d (%.2f %%)\n",
            sum(call_as_is), 100 * mean(call_as_is)))
cat(sprintf("  called by null_fit:         %d (%.2f %%)\n",
            sum(call_null_fit), 100 * mean(call_null_fit)))
cat(sprintf("  flips:                      %d (%.3f %%)\n", sum(flips), 100 * mean(flips)))
cat(sprintf("    null_fit only:            %d\n", sum(call_null_fit & !call_as_is)))
cat(sprintf("    as_is only:               %d\n", sum(call_as_is & !call_null_fit)))

# One-directional flips mean bias; a split means noise. At 5 replicates on one target this was 7-0.
# A binomial test on the direction says whether that survives at scale.
n_flip <- sum(flips)
if (n_flip > 0) {
  n_null_only <- sum(call_null_fit & !call_as_is)
  bt <- stats::binom.test(n_null_only, n_flip, p = 0.5)
  cat(sprintf("  direction: %d of %d flips favour null_fit (binomial p = %.3g)\n",
              n_null_only, n_flip, bt$p.value))
  cat("  A one-directional split is bias; an even one is noise.\n")
}

cat("\n=========== p-value agreement, for reference ===========\n")
d_p <- abs(both$p_value_as_is - both$p_value_null_fit)
cat(sprintf("  max|diff| = %.4g   median|diff| = %.4g   spearman = %.4f\n",
            max(d_p, na.rm = TRUE), stats::median(d_p, na.rm = TRUE),
            stats::cor(both$p_value_as_is, both$p_value_null_fit, method = "spearman",
                       use = "complete.obs")))
cat(sprintf("  within a factor of 10 of the threshold: %d of %d\n",
            sum(both$p_value_null_fit > threshold / 10 & both$p_value_null_fit < threshold * 10,
                na.rm = TRUE), nrow(both)))
cat("  Reported only to show that these aggregates hide the call-level disagreement above.\n")

## PER-PAIR POWER ==================================================================================

merged$power_difference <- merged$power_null_fit - merged$power_as_is
d <- merged$power_difference

cat("\n=========== per-pair power ===========\n")
cat(sprintf("  pairs:                      %d at %d replicates\n",
            nrow(merged), stats::median(merged$reps_as_is)))
cat(sprintf("  mean power, as_is:          %.4f\n", mean(merged$power_as_is)))
cat(sprintf("  mean power, null_fit:       %.4f\n", mean(merged$power_null_fit)))
cat(sprintf("  mean difference:            %+.4f  (null_fit - as_is)\n", mean(d)))
cat(sprintf("  mean |difference|:          %.4f\n", mean(abs(d))))
cat(sprintf("  max |difference|:           %.4f\n", max(abs(d))))
cat(sprintf("  pairs differing at all:     %d (%.1f %%)\n",
            sum(d != 0), 100 * mean(d != 0)))
cat(sprintf("  pairs where null_fit higher: %d\n", sum(d > 0)))
cat(sprintf("  pairs where as_is higher:    %d\n", sum(d < 0)))
if (any(d != 0)) {
  st <- stats::wilcox.test(merged$power_null_fit, merged$power_as_is, paired = TRUE)
  cat(sprintf("  paired signed-rank p:       %.3g\n", st$p.value))
}

# The deliverable is a per-pair claim of measured power >= 0.8, so what matters is not the mean shift
# but how many pairs change side of that line. See docs/output.md.
for (cut in c(0.8, 0.5)) {
  a <- sum(merged$power_as_is >= cut)
  b <- sum(merged$power_null_fit >= cut)
  moved <- sum((merged$power_as_is >= cut) != (merged$power_null_fit >= cut))
  cat(sprintf("  pairs at power >= %.1f:      as_is %d, null_fit %d (%d pairs change side)\n",
              cut, a, b, moved))
}

## WRITE ===========================================================================================

merged <- merged[order(-abs(merged$power_difference)), ]
write_tsv_file(merged[, c("pair", "reps_as_is", "power_as_is", "power_null_fit",
                          "power_difference")], opts$out)
log_step("Wrote ", opts$out)

cat("\nInterpretation: because the two runs share their simulated counts, a non-zero mean\n")
cat("difference is bias rather than noise, and its sign says which way the inherited cache\n")
cat("moved the reported power.\n")
