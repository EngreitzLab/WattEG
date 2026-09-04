#!/usr/bin/env Rscript
#
# Step 3 of paper/experiments.md: does a reduced sweep reproduce the full one?
#
# Reads the fitted curves produced by the design-space scan (14_design_space_scan.sbatch), predicts
# the effect sizes each design HELD OUT, and scores those predictions against the full
# six-point / 100-replicate sweep. One row out per (grid, replicates, seed).
#
# TWO QUESTIONS, KEPT SEPARATE BECAUSE THEY CAN FAIL INDEPENDENTLY
#
#   1. Interpolation. Does the fitted curve predict power at an effect size it never saw?
#      Scored as mean absolute error against the measured value, with the binomial noise floor of
#      that measurement reported alongside -- an MAE below the floor is not "better than perfect",
#      it means the comparison is noise-limited and the design is as good as measuring.
#
#   2. Minimum detectable effect size. The actual deliverable: one sentence per pair. Scored as
#      exact agreement with the reference MDES and agreement within one grid step.
#
# THE REFERENCE IS NOT TRUTH, AND THAT IS QUANTIFIED RATHER THAN ASSERTED
#
# The six-point MDES is itself an estimate from 100 binomial draws per point, so some disagreement
# is the REFERENCE being wrong, not the fit. Reporting raw agreement without that ceiling overstates
# the fit's error. So the reference is bootstrapped: for each pair and effect size, resample
# successes ~ Binomial(n_reps, measured power), recompute the MDES, and repeat. The resulting
# self-agreement -- how often a resampled reference agrees with the original -- is the ceiling any
# design could reach. A design at 89 % against a reference whose own ceiling is 91 % is close to
# measurement-limited; the same 89 % against a ceiling of 99 % is not.
#
# Parametric rather than resampling the stored replicates: the per-replicate files are 3 GB per
# effect size and the binomial is exactly the sampling distribution of the statistic, so drawing
# from it is equivalent and cheap.
#
# Usage:
#   compare_reduced_design.R --curves <dir> --reference power_summary.tsv \
#     --threshold-file discovery_threshold.txt --out reduced_design_summary.tsv

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    normalizePath(".")
  }
  for (dir in unique(c(file.path(dirname(here), "lib"), file.path(here, "lib"), here))) {
    if (file.exists(file.path(dir, "cli.R"))) {
      source(file.path(dir, "cli.R"))
      return(invisible(NULL))
    }
  }
  stop("Cannot find lib/cli.R relative to ", here, call. = FALSE)
})

suppressPackageStartupMessages(library(optparse))
source_lib("stats.R")

option_list <- list(
  make_option("--curves", type = "character", default = NULL, dest = "curves",
              help = "Directory of fitted-curve TSVs from 14_design_space_scan.sbatch."),
  make_option("--reference", type = "character", default = NULL, dest = "reference",
              help = "power_summary.tsv from the full six-point, 100-replicate sweep."),
  make_option("--threshold-file", type = "character", default = NULL, dest = "threshold_file",
              help = "The sweep's discovery_threshold.txt, for deriving z."),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per (grid, replicates, seed)."),
  make_option("--power-threshold", type = "double", default = 0.8, dest = "power_threshold",
              help = "Power level the MDES is defined at [default %default]."),
  make_option("--bootstrap", type = "integer", default = 200L, dest = "bootstrap",
              help = paste("Bootstrap draws for the reference's own noise ceiling [default",
                           "%default]. 0 skips it."))
)
opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Score reduced sweep designs against the full one."
))
require_options(opts, c("curves", "reference", "threshold_file", "out"))

ALL_ES <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)

# effect_label() from lib/stats.R, the same function summarize_power.R and fit_power_curve.R use to
# build these column names. Reimplementing it here would be one more place for the 0.2 -> "2" bug to
# come back.
es_label <- effect_label

## LOAD ============================================================================================

threshold <- as.numeric(readLines(opts$threshold_file, warn = FALSE)[1])
if (!is.finite(threshold) || threshold <= 0 || threshold >= 1) {
  stop("Unusable threshold in ", opts$threshold_file, call. = FALSE)
}
# Same derivation as fit_power_curve.R:176 for a run that passes neither --z nor --fit-z, which is
# what the scan does. Duplicated rather than read from the curve files because z is not among their
# columns; the self-check below is what stops the two drifting apart silently.
z <- stats::qnorm(1 - threshold)
log_step("z = ", format(z, digits = 6), " (from threshold ", format(threshold, digits = 6), ")")

reference <- read_tsv_file(opts$reference,
                           required_columns = c("grna_target", "response_id",
                                                "min_detectable_effect_size"))
log_step("Reference: ", nrow(reference), " pairs")

ref_key <- paste(reference$grna_target, reference$response_id, sep = "\r")

# Measured power and replicate count per effect size, as matrices in ALL_ES order.
ref_power <- vapply(ALL_ES, function(es) reference[[paste0("power_at_effect_size_", es_label(es))]],
                    numeric(nrow(reference)))
ref_nreps <- vapply(ALL_ES, function(es) reference[[paste0("power_at_effect_size_", es_label(es),
                                                           "_n_reps")]],
                    numeric(nrow(reference)))

## THE MDES RULE, IN ONE PLACE =====================================================================
#
# Qualifying requires clearing the power threshold at that effect size AND every larger one tested.
# Taking the first effect size that clears would let noise at a single point bias every pair toward
# looking more detectable than it is -- the rule summarize_power.R already applies, restated here so
# the bootstrap uses exactly the same definition as the reference it is bootstrapping.
mdes_from_power <- function(power_mat, es = ALL_ES) {
  n <- nrow(power_mat)
  out <- rep(NA_real_, n)
  # Vectorised transcription of summarize_power.R:174-182. Walk down from the largest effect size;
  # the first KNOWN failure disqualifies everything weaker, and the answer is the lowest clearing
  # effect size reached before it. An NA is untested, not a failure, so it neither clears nor
  # breaks -- getting that wrong would silently truncate the suffix at the first missing point.
  alive <- rep(TRUE, n)
  for (j in rev(seq_along(es))) {
    p <- power_mat[, j]
    known <- !is.na(p)
    out[alive & known & p >= opts$power_threshold] <- es[j]
    alive <- alive & !(known & p < opts$power_threshold)
  }
  out
}

ref_mdes <- reference$min_detectable_effect_size

# The reference's MDES recomputed from its own power columns must reproduce the column
# summarize_power.R wrote. If it does not, this script's rule has drifted from the pipeline's and
# every agreement number below would be measuring that drift instead of the design.
check_mdes <- mdes_from_power(ref_power)
disagree <- sum(!is.na(check_mdes) != !is.na(ref_mdes) |
                  (!is.na(check_mdes) & !is.na(ref_mdes) & check_mdes != ref_mdes))
if (disagree > 0) {
  stop("Recomputing the reference MDES from its own power columns disagrees with ",
       "min_detectable_effect_size for ", disagree, " pair(s). This script's rule and ",
       "summarize_power.R's have diverged; fix before trusting any number below.", call. = FALSE)
}
log_step("MDES rule reproduces the reference exactly (", nrow(reference), " pairs)")

## THE REFERENCE'S OWN CEILING =====================================================================

ceiling_exact <- NA_real_
ceiling_within1 <- NA_real_
if (opts$bootstrap > 0) {
  log_step("Bootstrapping the reference's own MDES, ", opts$bootstrap, " draws")
  set.seed(1L)
  agree_exact <- numeric(opts$bootstrap)
  agree_within1 <- numeric(opts$bootstrap)
  for (b in seq_len(opts$bootstrap)) {
    resampled <- ref_power
    for (j in seq_along(ALL_ES)) {
      resampled[, j] <- stats::rbinom(nrow(ref_power), ref_nreps[, j], ref_power[, j]) /
        ref_nreps[, j]
    }
    boot_mdes <- mdes_from_power(resampled)
    cmp <- compare_mdes(boot_mdes, ref_mdes)
    agree_exact[b] <- cmp$exact
    agree_within1[b] <- cmp$within1
  }
  ceiling_exact <- mean(agree_exact)
  ceiling_within1 <- mean(agree_within1)
  log_step(sprintf("  ceiling: %.1f%% exact, %.1f%% within one grid step",
                   100 * ceiling_exact, 100 * ceiling_within1))
}

## SCORE EACH DESIGN ===============================================================================

curve_files <- sort(list.files(opts$curves, pattern = "[.]tsv$", full.names = TRUE))
if (length(curve_files) == 0) {
  stop("No curve files in ", opts$curves, call. = FALSE)
}
log_step("Scoring ", length(curve_files), " design(s)")

rows <- lapply(curve_files, function(path) {
  # <grid with - separators>_<reps>_<seed>.tsv, e.g. 0.05-0.1-0.25_30_2.tsv
  stem <- sub("[.]tsv$", "", basename(path))
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) {
    log_step("  skipping unparseable name: ", basename(path))
    return(NULL)
  }
  grid <- as.numeric(strsplit(parts[1], "-", fixed = TRUE)[[1]])
  reps <- as.integer(parts[2])
  seed <- as.integer(parts[3])
  held_out <- setdiff(ALL_ES, grid)

  curves <- read_tsv_file(path, required_columns = c("grna_target", "response_id", "k"))
  # Align to the reference by pair, keeping only pairs the fit succeeded on.
  idx <- match(paste(curves$grna_target, curves$response_id, sep = "\r"), ref_key)
  keep <- !is.na(idx) & is.finite(curves$k)
  curves <- curves[keep, , drop = FALSE]
  idx <- idx[keep]

  # Fitted power from the model, not from the file's own columns, so held-out effect sizes -- which
  # the scan did not ask fit_power_curve.R to predict -- are available too.
  fitted_at <- function(es) stats::pnorm(curves$k * (-log(1 - es)) - z)

  # SELF-CHECK: on a FITTED effect size the file already carries the answer, so recomputing it here
  # must reproduce that column. This is what catches z or the formula drifting from
  # fit_power_curve.R, which is the one thing duplicating them risks.
  probe <- paste0("fitted_power_at_effect_size_", es_label(grid[1]))
  if (probe %in% colnames(curves)) {
    delta <- max(abs(fitted_at(grid[1]) - curves[[probe]]), na.rm = TRUE)
    if (!is.finite(delta) || delta > 1e-8) {
      stop("Recomputed fitted power differs from ", probe, " by ", format(delta),
           " in ", basename(path), ". z or the curve formula has drifted from ",
           "fit_power_curve.R.", call. = FALSE)
    }
  }

  # --- question 1: interpolation onto the held-out effect sizes --------------------------------
  errs <- unlist(lapply(held_out, function(es) {
    abs(fitted_at(es) - ref_power[idx, match(es, ALL_ES)])
  }), use.names = FALSE)
  floors <- unlist(lapply(held_out, function(es) {
    p <- ref_power[idx, match(es, ALL_ES)]
    sqrt(p * (1 - p) / ref_nreps[idx, match(es, ALL_ES)])
  }), use.names = FALSE)

  # --- question 2: does the reduced design reproduce the reference MDES? -----------------------
  fitted_power_all <- vapply(ALL_ES, fitted_at, numeric(nrow(curves)))
  fit_mdes <- mdes_from_power(fitted_power_all)
  cmp <- compare_mdes(fit_mdes, ref_mdes[idx])

  data.frame(
    grid = parts[1], n_points = length(grid), reps = reps, seed = seed,
    held_out = paste(held_out, collapse = ","),
    pairs_fitted = nrow(curves),
    mae = mean(errs, na.rm = TRUE),
    p90_error = as.numeric(stats::quantile(errs, 0.9, na.rm = TRUE)),
    noise_floor = mean(floors, na.rm = TRUE),
    mdes_exact = cmp$exact,
    mdes_within1 = cmp$within1,
    mdes_ceiling_exact = ceiling_exact,
    mdes_ceiling_within1 = ceiling_within1,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
summary_df <- summary_df[order(summary_df$mae), , drop = FALSE]
rownames(summary_df) <- NULL

write_tsv_file(summary_df, opts$out)
log_step("Wrote ", nrow(summary_df), " design(s) to ", opts$out)

message("\nBest designs by held-out MAE:")
top <- utils::head(summary_df, 10)
for (i in seq_len(nrow(top))) {
  message(sprintf("  grid %-18s reps %3d seed %d | MAE %.4f (floor %.4f, %.2fx) | MDES %.1f%% exact, %.1f%% +-1",
                  top$grid[i], top$reps[i], top$seed[i], top$mae[i], top$noise_floor[i],
                  top$mae[i] / top$noise_floor[i], 100 * top$mdes_exact[i],
                  100 * top$mdes_within1[i]))
}
if (is.finite(ceiling_exact)) {
  message(sprintf("\nReference's own ceiling: %.1f%% exact, %.1f%% within one step. ",
                  100 * ceiling_exact, 100 * ceiling_within1),
          "A design cannot beat this, and one close to it is measurement-limited rather than ",
          "model-limited.")
}
