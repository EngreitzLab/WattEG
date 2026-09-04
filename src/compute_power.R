#!/usr/bin/env Rscript
#
# Turn per-replicate simulation results into a power estimate per element-gene pair.
#
# Power is the fraction of replicates in which sceptre would have called the association:
#
#   power = mean(p_value < threshold & log_2_fold_change < 0)
#
# The log_2_fold_change condition makes this one-sided -- a replicate only counts if the simulated
# perturbation reduced expression.
#
# Because that is a binomial proportion over a finite number of replicates, the point estimate is
# reported with a Wilson score interval. See docs/choosing-num-replicates.md for how to pick the
# replicate count and how to read the intervals.
#
# Three fixes relative to compute_power_from_simulations.R:
#
#  * It no longer needs a separate discovery-results file. That file was an undeclared fourth
#    manual input which no rule produced; the threshold comes from @discovery_result inside the
#    sceptre object, written out by prepare_sim_input.R (verified identical on sample1: both give
#    7.71e-4).
#  * `group_by(rep) %>% group_by(grna_target, response_id)` -- the first grouping was silently
#    discarded by the second, so it was dead code.
#  * `max()` over an empty set of significant discoveries returned -Inf and silently reported zero
#    power for every pair. The threshold is now validated.
#
# Usage:
#   compute_power.R --simulations combined_es0.15.tsv \
#     --threshold-file discovery_threshold.txt --out power_es0.15.tsv

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    normalizePath(".")
  }
  # lib/ sits one level up from src/, but a workflow engine that stages every input into one flat
  # task directory collapses that, so look beside the script and in it as well. Kept in step with
  # lib_dirs() in cli.R, which resolves the rest of the library the same way.
  for (dir in unique(c(file.path(dirname(here), "lib"), file.path(here, "lib"), here))) {
    if (file.exists(file.path(dir, "cli.R"))) {
      source(file.path(dir, "cli.R"))
      return(invisible(NULL))
    }
  }
  stop("Cannot find lib/cli.R relative to ", here, call. = FALSE)
})

source_lib("stats.R")

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--simulations", type = "character", default = NULL, dest = "simulations",
              help = paste("Per-replicate simulation results, one or more TSVs (comma-separated).",
                           "Concatenated per-split/per-chunk outputs of run_power_simulation.R.")),
  make_option("--threshold-file", type = "character", default = NULL, dest = "threshold_file",
              help = paste("File containing the nominal p-value threshold, as written by",
                           "prepare_sim_input.R. Mutually exclusive with --alpha.")),
  make_option("--alpha", type = "double", default = NULL, dest = "alpha",
              help = paste("Nominal p-value threshold, set explicitly. Use this only when the",
                           "sceptre object has no discovery result to derive one from -- the",
                           "derived threshold reflects the multiple-testing correction actually",
                           "applied to the real data, whereas a bare alpha does not.")),
  make_option("--conf-level", type = "double", default = 0.95, dest = "conf_level",
              help = "Confidence level for the Wilson interval [default %default]."),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per pair."),
  make_option("--subsample-reps", type = "integer", default = NULL, dest = "subsample_reps",
              help = paste("Compute power from this many randomly chosen replicates instead of",
                           "all of them. For the reduced-design study in broadinstitute/WattEG-paper;",
                           "has no place in a production run, where every replicate that was paid",
                           "for should be used.")),
  make_option("--subsample-seed", type = "integer", default = 1L, dest = "subsample_seed",
              help = paste("Seed for --subsample-reps [default %default]. Vary it to repeat the",
                           "same reduced design over different subsets and report the spread,",
                           "which is what turns a single comparison into an error bar."))
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Compute power per pair from per-replicate simulation results."
))
require_options(opts, c("simulations", "out"))

if (is.null(opts$threshold_file) == is.null(opts$alpha)) {
  stop("Provide exactly one of --threshold-file or --alpha.", call. = FALSE)
}
if (opts$conf_level <= 0 || opts$conf_level >= 1) {
  stop("--conf-level must be strictly between 0 and 1.", call. = FALSE)
}

## THRESHOLD =======================================================================================

threshold <- if (!is.null(opts$alpha)) {
  log_step("Using explicit threshold --alpha ", opts$alpha)
  opts$alpha
} else {
  value <- suppressWarnings(as.numeric(readLines(opts$threshold_file, warn = FALSE)[1]))
  log_step(sprintf("Using threshold %.6g from %s", value, opts$threshold_file))
  value
}
if (!is.finite(threshold) || threshold <= 0 || threshold > 1) {
  stop("The p-value threshold (", threshold, ") is not a usable probability. ",
       "A non-finite value usually means the source discovery results contained no significant ",
       "pair; pass --alpha to set the threshold explicitly.", call. = FALSE)
}

## LOAD ============================================================================================

paths <- trimws(strsplit(opts$simulations, ",", fixed = TRUE)[[1]])

# A directory expands to the TSVs inside it. The Nextflow runner passes an explicit file list, which
# is what makes the replicate-count check below meaningful -- it validates the files it was handed.
# But a sweep publishes 1,000 files per effect size, and naming them all on a command line to
# re-aggregate an existing sweep (the reduced-design study does this hundreds of times) is a
# 100 KB argument list for no benefit. Sorted, so the row order does not depend on readdir order.
expanded <- unlist(lapply(paths, function(p) {
  if (dir.exists(p)) sort(list.files(p, pattern = "[.]tsv$", full.names = TRUE)) else p
}), use.names = FALSE)
if (length(expanded) == 0) {
  stop("--simulations matched no files: ", opts$simulations, call. = FALSE)
}
paths <- expanded
required <- c("grna_target", "response_id", "p_value", "log_2_fold_change", "rep")
sims <- do.call(rbind, lapply(paths, read_tsv_file, required_columns = required))
log_step("Read ", nrow(sims), " replicate rows from ", length(paths), " file(s)")

# Replicates with no fold change estimate carry no information about whether the association would
# have been called, so they are dropped -- which means n_reps can fall below --reps and has to be
# reported per pair rather than assumed.
n_before <- nrow(sims)
sims <- sims[!is.na(sims$log_2_fold_change) & !is.na(sims$p_value), , drop = FALSE]
if (nrow(sims) < n_before) {
  log_step("Dropped ", n_before - nrow(sims), " replicate row(s) with a missing p_value or ",
           "log_2_fold_change")
}
if (nrow(sims) == 0) {
  stop("No usable replicate rows remain.", call. = FALSE)
}

if (anyDuplicated(sims[, c("grna_target", "response_id", "rep")])) {
  n_dup <- sum(duplicated(sims[, c("grna_target", "response_id", "rep")]))
  stop(n_dup, " duplicated (grna_target, response_id, rep) row(s). Overlapping --rep-offset ",
       "values across chunks would double-count replicates.", call. = FALSE)
}

# --- optional: keep only a subsample of the replicates ------------------------------------------
#
# For the reduced-design study in broadinstitute/WattEG-paper, which asks whether fewer replicates
# and fewer effect sizes reproduce the full design. That question is answered by SUBSAMPLING an existing
# sweep rather than by running smaller ones: replicates are i.i.d. draws, so a random subset of 30
# is distributed exactly like a fresh 30-replicate run, and subsampling makes the comparison paired
# -- the same pairs, the same simulated data -- which a rerun would not be.
#
# The SAME replicate ids are kept for every pair, deliberately. A replicate index is arbitrary, so
# one shared subset is the faithful analogue of "we ran 30 replicates"; drawing an independent
# subset per pair would average over more of the replicate noise than a real 30-replicate run does
# and would flatter the reduced design.
#
# Sampling from the observed ids rather than 1:reps because rows with a missing fold change have
# already been dropped above, so a pair can have fewer than --reps replicates present.
if (!is.null(opts$subsample_reps)) {
  available <- sort(unique(sims$rep))
  if (opts$subsample_reps > length(available)) {
    stop("--subsample-reps ", opts$subsample_reps, " exceeds the ", length(available),
         " replicate(s) present.", call. = FALSE)
  }
  set.seed(opts$subsample_seed)
  keep <- sort(sample(available, opts$subsample_reps))
  sims <- sims[sims$rep %in% keep, , drop = FALSE]
  log_step("Subsampled to ", length(keep), " of ", length(available), " replicates ",
           "(seed ", opts$subsample_seed, "): ", nrow(sims), " rows remain")
  if (nrow(sims) == 0) {
    stop("No replicate rows remain after subsampling.", call. = FALSE)
  }
}

## AGGREGATE =======================================================================================

pair_key <- paste(sims$grna_target, sims$response_id, sep = "\r")
success <- sims$p_value < threshold & sims$log_2_fold_change < 0

agg <- function(values, fun) vapply(split(values, pair_key), fun, numeric(1))

n_reps <- agg(rep(1, nrow(sims)), sum)
successes <- agg(as.numeric(success), sum)
interval <- wilson_interval(successes, n_reps, opts$conf_level)

keys <- do.call(rbind, strsplit(names(n_reps), "\r", fixed = TRUE))
power <- data.frame(
  grna_target = keys[, 1],
  response_id = keys[, 2],
  power = successes / n_reps,
  power_ci_low = interval$low,
  power_ci_high = interval$high,
  n_reps = n_reps,
  mean_log_2_fold_change = agg(sims$log_2_fold_change, mean),
  stringsAsFactors = FALSE
)

# Carried through when the simulation reported them, so the output is self-describing.
for (column in c("num_pert_cells", "average_expression_all_cells")) {
  if (column %in% colnames(sims)) {
    power[[sub("^num_", "mean_", column)]] <- agg(sims[[column]], mean)
  }
}
if ("effect_size" %in% colnames(sims)) {
  effect_sizes <- unique(sims$effect_size)
  if (length(effect_sizes) > 1) {
    stop("The input mixes effect sizes (", paste(effect_sizes, collapse = ", "),
         "). Run compute_power.R once per effect size.", call. = FALSE)
  }
  power$effect_size <- effect_sizes
}

power <- power[order(-power$power, power$grna_target, power$response_id), ]
rownames(power) <- NULL

## WRITE ===========================================================================================

write_tsv_file(power, opts$out)
log_step("Wrote ", nrow(power), " pairs to ", opts$out)

reps_range <- range(power$n_reps)
message(sprintf("  replicates per pair: %d-%d", reps_range[1], reps_range[2]))
message(sprintf("  mean power: %.3f | pairs at 0: %d | at 1: %d | in (0.1,0.9): %d of %d",
                mean(power$power), sum(power$power == 0), sum(power$power == 1),
                sum(power$power > 0.1 & power$power < 0.9), nrow(power)))
message(sprintf("  median 95%% CI width: %.3f  (widest %.3f)",
                stats::median(power$power_ci_high - power$power_ci_low),
                max(power$power_ci_high - power$power_ci_low)))
if (min(power$n_reps) < 100 && is.null(opts$subsample_reps)) {
  message("  note: with fewer than ~100 replicates a per-pair estimate is coarse; see ",
          "docs/choosing-num-replicates.md")
}
