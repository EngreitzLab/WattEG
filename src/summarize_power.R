#!/usr/bin/env Rscript
#
# Merge the per-effect-size power tables into one wide table, one row per element-gene pair.
#
# Without this the pipeline leaves you with a separate file per effect size and no combined view.
# This recovers the useful half of the old format_sceptre_output.R -- a single table with one
# power column per effect size -- without that script's dependencies on `distances`,
# `guide_targets` and `features` files, none of which the pipeline produces.
#
# It also derives the quantity the analysis is usually actually after: the smallest tested effect
# size at which a pair reaches a target power. See --power-threshold.
#
# Usage:
#   summarize_power.R --power power_es0.15.tsv,power_es0.2.tsv --out power_summary.tsv

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
  make_option("--power", type = "character", default = NULL, dest = "power",
              help = "Comma-separated compute_power.R outputs, one per effect size."),
  make_option("--power-threshold", type = "double", default = 0.8, dest = "power_threshold",
              help = paste("Power level used to derive min_detectable_effect_size: the smallest",
                           "tested effect size at which a pair reaches it [default %default].")),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per pair."),
  make_option("--sim-input", type = "character", default = NULL, dest = "sim_input",
              help = paste("sim_input.rds from prepare_sim_input.R. Optional. Supplied, the gene's",
                           "negative-binomial dispersion and normalised mean are joined onto each",
                           "pair, which makes this table self-sufficient for covariate modelling.",
                           "Omitted, those columns are absent and everything else is unchanged."))
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Merge per-effect-size power tables into one wide table."
))
require_options(opts, c("power", "out"))

if (opts$power_threshold <= 0 || opts$power_threshold > 1) {
  stop("--power-threshold must be in (0, 1].", call. = FALSE)
}

## LOAD ============================================================================================

paths <- trimws(strsplit(opts$power, ",", fixed = TRUE)[[1]])
required <- c("grna_target", "response_id", "power", "effect_size")

tables <- lapply(paths, function(path) {
  df <- read_tsv_file(path, required_columns = required)
  effect_sizes <- unique(df$effect_size)
  if (length(effect_sizes) != 1) {
    stop(path, " contains ", length(effect_sizes), " effect sizes; expected exactly one.",
         call. = FALSE)
  }
  df
})
effect_sizes <- vapply(tables, function(df) df$effect_size[1], numeric(1))

if (anyDuplicated(effect_sizes)) {
  stop("Two or more inputs share an effect size (",
       paste(effect_sizes[duplicated(effect_sizes)], collapse = ", "), ").", call. = FALSE)
}

# Ascending, so min_detectable_effect_size below can take the first that clears the threshold and
# the columns read left to right from weakest to strongest perturbation.
ordering <- order(effect_sizes)
tables <- tables[ordering]
effect_sizes <- effect_sizes[ordering]
log_step("Merging ", length(tables), " effect size(s): ", paste(effect_sizes, collapse = ", "))


## MERGE ===========================================================================================

# Outer join on the pair, so a pair tested at only some effect sizes still appears (with NA
# elsewhere) rather than being silently dropped.
pairs <- unique(do.call(rbind, lapply(tables, function(df) df[, c("grna_target", "response_id")])))
pairs <- pairs[order(pairs$grna_target, pairs$response_id), , drop = FALSE]
rownames(pairs) <- NULL
key_of <- function(df) paste(df$grna_target, df$response_id, sep = "\r")
pair_keys <- key_of(pairs)

out <- pairs

# Columns shared across effect sizes: taken from the first table that has the pair. Flagged if they
# disagree between effect sizes, which would mean the inputs do not describe the same experiment.
shared <- c("mean_pert_cells", "average_expression_all_cells")
for (column in shared) {
  if (!any(vapply(tables, function(df) column %in% colnames(df), logical(1)))) next
  values <- matrix(NA_real_, nrow = nrow(pairs), ncol = length(tables))
  for (i in seq_along(tables)) {
    if (!column %in% colnames(tables[[i]])) next
    values[, i] <- tables[[i]][[column]][match(pair_keys, key_of(tables[[i]]))]
  }
  spread <- apply(values, 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < 2) 0 else diff(range(v)) / max(abs(v), 1e-12)
  })
  if (any(spread > 1e-6, na.rm = TRUE)) {
    message("  note: '", column, "' differs between effect sizes for ",
            sum(spread > 1e-6, na.rm = TRUE), " pair(s); reporting the first non-missing value.")
  }
  out[[column]] <- apply(values, 1, function(v) if (all(is.na(v))) NA_real_ else v[!is.na(v)][1])
}

# Per-gene values from sim_input, joined here rather than carried through the simulation.
#
# WHY HERE. The theory behind the power curve is SE^2 ~ (1/n_pert_cells) * (1/mu + 1/theta), so any
# covariate model of power needs the gene's dispersion as well as its expression. Only
# `average_expression_all_cells` reaches this table through the per-replicate output, so every such
# analysis has had to load a 16 MB sim_input.rds to find the other half. Joining it here costs one
# read and makes the summary self-sufficient.
#
# It is done at this step, and not by having run_power_simulation.R emit a `dispersion` column,
# because that would only help sweeps run AFTER the change -- an existing sweep's per-replicate
# output is already written, and re-running it to add a per-gene constant would cost thousands of
# CPU-hours. Summarising again from stored power tables costs seconds.
#
# `gene_mean` comes along because it is the mu in that formula and is NOT the same quantity as
# `average_expression_all_cells`: they correlate at r = 0.9999 but differ by a scale factor (median
# ratio 0.84 on day0), since one is size-factor normalised and the other is not. Adding only the
# dispersion would leave the formula still needing the RDS, which is the problem this solves.
#
# Optional throughout: with no --sim-input the columns are simply absent, so nothing downstream
# that predates them can break.
if (!is.null(opts$sim_input)) {
  sim <- readRDS(opts$sim_input)
  rd <- sim$row_data
  missing_cols <- setdiff(c("mean", "dispersion"), colnames(rd))
  if (length(missing_cols) > 0) {
    stop("--sim-input's row_data has no ", paste(missing_cols, collapse = ", "), " column(s). ",
         "See lib/sim_input.R for the contract.", call. = FALSE)
  }
  gi <- match(out$response_id, rownames(rd))
  n_missing <- sum(is.na(gi))
  if (n_missing > 0) {
    log_step("note: ", n_missing, " of ", nrow(out), " pair(s) have a gene absent from ",
             "--sim-input; their gene_mean and dispersion are NA")
  }
  # `dispersion` is 1/theta, not theta -- built that way in lib/simulate.R and consumed as
  # rnbinom(size = 1/dispersion). The name is kept identical to sim_input's so the two cannot be
  # read as different quantities; the theory term is therefore `1/gene_mean + dispersion`.
  out$gene_mean <- rd$mean[gi]
  out$dispersion <- rd$dispersion[gi]

  # average_expression_all_cells comes from here too, and only from here, now that
  # run_power_simulation.R no longer repeats it once per replicate. If an older power table still
  # carries it, the value already set from `shared` above is kept -- they are the same number, and
  # not overwriting it keeps a re-summarised old sweep byte-identical.
  if (!"average_expression_all_cells" %in% colnames(out) ||
        all(is.na(out$average_expression_all_cells))) {
    out$average_expression_all_cells <- rd$average_expression_all_cells[gi]
  }
  log_step("Joined per-gene mean, dispersion and expression from ", opts$sim_input)
}

# Per-effect-size columns.
power_columns <- character(0)
for (i in seq_along(tables)) {
  df <- tables[[i]]
  idx <- match(pair_keys, key_of(df))
  label <- effect_label(effect_sizes[i])
  base <- paste0("power_at_effect_size_", label)
  power_columns <- c(power_columns, base)

  out[[base]] <- df$power[idx]
  for (from in c("power_ci_low", "power_ci_high", "n_reps")) {
    if (from %in% colnames(df)) {
      suffix <- switch(from, power_ci_low = "_ci_low", power_ci_high = "_ci_high",
                       n_reps = "_n_reps")
      out[[paste0(base, suffix)]] <- df[[from]][idx]
    }
  }
}

## DERIVED =========================================================================================

# Smallest tested effect size from which a pair is detectable. NA means no tested effect size
# qualified, which is a statement about the effect sizes you ran, not proof the pair is
# undetectable -- so the largest effect size tested is reported alongside it for context.
#
# TWO CHOICES ARE BAKED IN HERE AND BOTH CHANGE HOW THE COLUMN CAN BE USED.
#
# 1. THE SUFFIX RULE. A pair qualifies at effect size e only if it clears the threshold at e *and*
#    at every larger effect size tested. Taking the first effect size that clears, in isolation,
#    lets Monte-Carlo noise win: at 100 replicates a pair whose true power is 0.75 clears 0.8
#    around a third of the time, so across six effect sizes a spurious early clear is likely, and
#    the error is one-directional -- it always reports the pair as more detectable than it is.
#    "Detectable from e upwards" is also the claim the column gets used to make.
#
#    Effect sizes a pair was not tested at are unknown, not failures, so they cannot block a
#    suffix and are skipped.
#
# 2. THREE BASES, because the interval on power carries through to an interval on this. Power is
#    monotone in effect size, so thresholding a *lower* bound on power yields a *larger* effect
#    size -- the direction inverts:
#
#      min_detectable_effect_size           from `power`          point estimate
#      min_detectable_effect_size_ci_low    from `power_ci_high`  optimistic edge
#      min_detectable_effect_size_ci_high   from `power_ci_low`   conservative edge
#
#    For "this pair was powered well enough that a non-significant result means something", the
#    conservative edge is the one to use: the smallest knockdown the data can *certify* the assay
#    would have caught. See docs/output.md.
min_detectable <- function(columns) {
  values <- as.matrix(out[, columns, drop = FALSE])
  apply(values, 1, function(row) {
    best <- NA_real_
    # Walk down from the largest effect size. The first known failure disqualifies everything
    # weaker than it, so the answer is the lowest clearing effect size reached before that.
    for (i in rev(seq_along(row))) {
      if (is.na(row[i])) next
      if (row[i] < opts$power_threshold) break
      best <- effect_sizes[i]
    }
    best
  })
}

out$min_detectable_effect_size <- min_detectable(power_columns)

ci_low_columns <- paste0(power_columns, "_ci_low")
ci_high_columns <- paste0(power_columns, "_ci_high")
have_intervals <- all(c(ci_low_columns, ci_high_columns) %in% colnames(out))
if (have_intervals) {
  out$min_detectable_effect_size_ci_low <- min_detectable(ci_high_columns)
  out$min_detectable_effect_size_ci_high <- min_detectable(ci_low_columns)
} else {
  message("  note: no power_ci_low / power_ci_high columns in the inputs; ",
          "reporting the point estimate only. A negative result cannot be certified from it.")
}

out$max_effect_size_tested <- max(effect_sizes)

out <- out[order(out$min_detectable_effect_size, na.last = TRUE,
                 decreasing = FALSE), , drop = FALSE]
rownames(out) <- NULL

## WRITE ===========================================================================================

write_tsv_file(out, opts$out)
log_step("Wrote ", nrow(out), " pairs x ", length(effect_sizes), " effect size(s) to ", opts$out)

message(sprintf("  reaching power >= %.2f at some tested effect size: %d of %d (%.0f%%)",
                opts$power_threshold, sum(!is.na(out$min_detectable_effect_size)), nrow(out),
                100 * mean(!is.na(out$min_detectable_effect_size))))
if (have_intervals) {
  # The gap between these two is the cost of insisting on a certified negative rather than a
  # point-estimate one, and it is the number to quote when reporting how many pairs the analysis
  # can actually say anything about.
  certified <- sum(!is.na(out$min_detectable_effect_size_ci_high))
  message(sprintf("  certifiable (power_ci_low >= %.2f): %d of %d pairs (%.0f%%)",
                  opts$power_threshold, certified, nrow(out), 100 * certified / nrow(out)))
  stricter <- sum(out$min_detectable_effect_size_ci_high > out$min_detectable_effect_size,
                  na.rm = TRUE)
  message(sprintf("  of those, %d land at a larger effect size than the point estimate suggests",
                  stricter))
}
for (i in seq_along(effect_sizes)) {
  column <- power_columns[i]
  values <- out[[column]]
  message(sprintf("  %-24s mean power %.3f | >= %.2f in %d of %d pairs",
                  column, mean(values, na.rm = TRUE), opts$power_threshold,
                  sum(values >= opts$power_threshold, na.rm = TRUE), sum(!is.na(values))))
}
