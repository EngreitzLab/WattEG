#!/usr/bin/env Rscript
#
# Fit each pair's power curve across effect sizes, so a sweep needs three effect sizes rather than
# six and the minimum detectable effect size stops being snapped to whichever values were run.
#
# PROTOTYPE. Validated on one dataset only -- see Step 11 in docs/status.md for the held-out test
# and its caveats. Do not use it to replace a measured sweep without repeating that test on the
# dataset at hand.
#
# THE MODEL, AND WHY IT IS NOT AN ARBITRARY CURVE FIT
#
# For a Wald-type test at a fixed rejection threshold, power for one pair is
#
#   power(effect_size) = Phi(beta / SE - z),   beta = -log(1 - effect_size),  z = qnorm(1 - alpha)
#
# `beta` is the effect on the scale the test works on, `SE` collects everything pair-specific
# (perturbed cells, expression, dispersion) and `z` is fixed by the discovery threshold and shared
# by every pair. So on the probit scale the curve is a straight line through -z with slope 1/SE:
# one free parameter per pair. Three effect sizes therefore leave two degrees of freedom to check
# the fit rather than just enough to force it.
#
# That is fitted here as a binomial GLM with a probit link, no intercept, and -z as an offset. The
# GLM is used rather than least squares on qnorm(power) because it handles 0/100 and 100/100
# correctly: a saturated point still carries information about the slope, where qnorm() would be
# infinite and the point would have to be discarded.
#
# Usage:
#   fit_power_curve.R --power power_es0.05.tsv,power_es0.1.tsv,power_es0.25.tsv \
#     --threshold-file prepared/discovery_threshold.txt --out power_curves.tsv

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

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--power", type = "character", default = NULL, dest = "power",
              help = "Comma-separated compute_power.R outputs, one per effect size (>= 2)."),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per pair."),
  make_option("--power-threshold", type = "double", default = 0.8, dest = "power_threshold",
              help = "Power level defining the minimum detectable effect size [default %default]."),
  make_option("--threshold-file", type = "character", default = NULL, dest = "threshold_file",
              help = paste("File holding the discovery p-value threshold, as written by",
                           "prepare_sim_input.R. z is taken as qnorm(1 - threshold).")),
  make_option("--z", type = "double", default = NULL, dest = "z",
              help = "Set z directly instead of deriving it from --threshold-file."),
  make_option("--fit-z", action = "store_true", default = FALSE, dest = "fit_z",
              help = paste("Estimate one z shared across all pairs by profiling total deviance,",
                           "instead of taking it from the threshold. Use when the threshold and",
                           "the observed curves disagree -- see docs/status.md, Step 11.")),
  make_option("--predict-at", type = "character", default = NULL, dest = "predict_at",
              help = paste("Comma-separated effect sizes to report fitted power at, in addition to",
                           "the ones supplied. These need not have been simulated.")),
  make_option("--conf-level", type = "double", default = 0.95, dest = "conf_level",
              help = "Confidence level for the interval on the effect size [default %default].")
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Fit a one-parameter power curve per pair across effect sizes."
))
require_options(opts, c("power", "out"))

if (opts$power_threshold <= 0 || opts$power_threshold >= 1) {
  stop("--power-threshold must be strictly between 0 and 1.", call. = FALSE)
}
if (is.null(opts$z) && is.null(opts$threshold_file) && !opts$fit_z) {
  stop("Supply one of --z, --threshold-file or --fit-z: the curve is anchored at -z and cannot ",
       "be fitted without it.", call. = FALSE)
}

## LOAD ============================================================================================

paths <- trimws(strsplit(opts$power, ",", fixed = TRUE)[[1]])
if (length(paths) < 2) {
  stop("At least two effect sizes are needed to fit a curve; got ", length(paths), ".",
       call. = FALSE)
}
required <- c("grna_target", "response_id", "power", "n_reps", "effect_size")

tables <- lapply(paths, function(path) {
  df <- read_tsv_file(path, required_columns = required)
  sizes <- unique(df$effect_size)
  if (length(sizes) != 1) {
    stop(path, " contains ", length(sizes), " effect sizes; expected exactly one.", call. = FALSE)
  }
  df
})
effect_sizes <- vapply(tables, function(df) df$effect_size[1], numeric(1))
if (anyDuplicated(effect_sizes)) {
  stop("Two or more inputs share an effect size.", call. = FALSE)
}
ordering <- order(effect_sizes)
tables <- tables[ordering]
effect_sizes <- effect_sizes[ordering]
log_step("Fitting across ", length(effect_sizes), " effect sizes: ",
         paste(effect_sizes, collapse = ", "))

# One row per (pair, effect size), which is the shape the GLM wants.
long <- do.call(rbind, lapply(tables, function(df) {
  data.frame(key = paste(df$grna_target, df$response_id, sep = "\r"),
             grna_target = df$grna_target, response_id = df$response_id,
             effect_size = df$effect_size,
             # Successes rather than the proportion: compute_power.R reports power = successes /
             # n_reps, and n_reps can differ between effect sizes for the same pair when a
             # simulation produced no fold-change estimate.
             successes = round(df$power * df$n_reps), n_reps = df$n_reps,
             stringsAsFactors = FALSE)
}))
long <- long[!is.na(long$successes) & long$n_reps > 0, , drop = FALSE]
long$beta <- -log(1 - long$effect_size)

## z ===============================================================================================

#' Total residual deviance over all pairs at a candidate z, used to profile z out of the model.
total_deviance <- function(z, split_data) {
  sum(vapply(split_data, function(d) {
    fit <- try(stats::glm(cbind(successes, n_reps - successes) ~ beta - 1, data = d,
                          family = stats::binomial(link = "probit"),
                          offset = rep(-z, nrow(d))), silent = TRUE)
    if (inherits(fit, "try-error")) NA_real_ else stats::deviance(fit)
  }, numeric(1)), na.rm = TRUE)
}

split_data <- split(long[, c("beta", "successes", "n_reps")], long$key)
log_step(length(split_data), " pairs, ",
         round(mean(vapply(split_data, nrow, integer(1))), 2), " effect sizes each on average")

if (opts$fit_z) {
  # Coarse grid then a local refinement: the deviance is smooth in z and this avoids depending on
  # an optimiser for a one-dimensional problem. Profiled on a subsample -- z is shared by every
  # pair, so a few thousand of them determine it as well as all of them.
  sample_keys <- names(split_data)
  if (length(sample_keys) > 2000) {
    sample_keys <- sample_keys[seq(1, length(sample_keys), length.out = 2000)]
  }
  subset_data <- split_data[sample_keys]
  coarse <- seq(1, 5, by = 0.25)
  dev_coarse <- vapply(coarse, total_deviance, numeric(1), split_data = subset_data)
  best <- coarse[which.min(dev_coarse)]
  fine <- seq(max(0.25, best - 0.25), best + 0.25, by = 0.05)
  dev_fine <- vapply(fine, total_deviance, numeric(1), split_data = subset_data)
  z <- fine[which.min(dev_fine)]
  log_step("Profiled z = ", z, " over ", length(subset_data), " pairs")
  if (!is.null(opts$threshold_file) && file.exists(opts$threshold_file)) {
    thr <- as.numeric(readLines(opts$threshold_file, warn = FALSE)[1])
    message(sprintf("  note: the threshold %.3g implies z = %.3f; the fitted value is %.3f.",
                    thr, stats::qnorm(1 - thr), z))
  }
} else if (!is.null(opts$z)) {
  z <- opts$z
  log_step("Using z = ", z, " as supplied")
} else {
  thr <- as.numeric(readLines(opts$threshold_file, warn = FALSE)[1])
  if (is.na(thr) || thr <= 0 || thr >= 1) {
    stop("Could not read a p-value threshold from ", opts$threshold_file, ".", call. = FALSE)
  }
  z <- stats::qnorm(1 - thr)
  log_step("Using z = ", round(z, 4), " from threshold ", thr)
}

## FIT =============================================================================================

# Slope k = 1 / SE. Larger k means the pair is easier to detect, so the effect size solving
# Phi(k * beta - z) = power_threshold is
#
#   beta* = (qnorm(power_threshold) + z) / k,   effect_size* = 1 - exp(-beta*)
#
# and the direction inverts: a *lower* bound on k gives a *larger* minimum detectable effect size,
# which is the conservative end and the one a false-negative argument needs. Same inversion as
# min_detectable_effect_size_ci_high in summarize_power.R -- see docs/output.md.
target_q <- stats::qnorm(opts$power_threshold)
zq <- stats::qnorm(1 - (1 - opts$conf_level) / 2)

effect_size_at <- function(k) {
  # NA rather than a number when the curve never reaches the target power at any effect size below
  # a complete knockdown: k <= 0 means power does not increase with effect size in this fit.
  out <- rep(NA_real_, length(k))
  ok <- !is.na(k) & k > 0
  out[ok] <- 1 - exp(-(target_q + z) / k[ok])
  out[!is.na(out) & (out <= 0 | out >= 1)] <- NA_real_
  out
}

fit_one <- function(d) {
  n_points <- nrow(d)
  saturated <- all(d$successes == d$n_reps)
  empty <- all(d$successes == 0)
  fit <- try(stats::glm(cbind(successes, n_reps - successes) ~ beta - 1, data = d,
                        family = stats::binomial(link = "probit"),
                        offset = rep(-z, n_points)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(c(k = NA_real_, se_k = NA_real_, deviance = NA_real_, df = NA_real_,
             n_points = n_points, status = 3))
  }
  co <- summary(fit)$coefficients
  c(k = unname(co[1, 1]), se_k = unname(co[1, 2]),
    deviance = stats::deviance(fit), df = stats::df.residual(fit),
    n_points = n_points,
    # 0 fitted, 1 every point saturated, 2 every point empty: the last two are reported rather
    # than dropped, because "detectable everywhere" and "detectable nowhere" are real answers.
    status = if (saturated) 1 else if (empty) 2 else 0)
}

fits <- do.call(rbind, lapply(split_data, fit_one))
out <- data.frame(key = rownames(fits), fits, row.names = NULL, stringsAsFactors = FALSE)

keys <- do.call(rbind, strsplit(out$key, "\r", fixed = TRUE))
out$grna_target <- keys[, 1]
out$response_id <- keys[, 2]
out$key <- NULL

out$min_detectable_effect_size <- effect_size_at(out$k)
out$min_detectable_effect_size_ci_low <- effect_size_at(out$k + zq * out$se_k)
out$min_detectable_effect_size_ci_high <- effect_size_at(out$k - zq * out$se_k)

# Fitted power at the effect sizes supplied, so the fit can be compared against what was measured,
# plus any extra ones asked for.
predict_at <- effect_sizes
if (!is.null(opts$predict_at)) {
  extra <- as.numeric(trimws(strsplit(opts$predict_at, ",", fixed = TRUE)[[1]]))
  if (anyNA(extra) || any(extra <= 0 | extra >= 1)) {
    stop("--predict-at must be effect sizes strictly between 0 and 1.", call. = FALSE)
  }
  predict_at <- sort(unique(c(predict_at, extra)))
}
for (es in predict_at) {
  label <- sub("\\.", "_", sub("^0\\.", "", format(es * 100, trim = TRUE)))
  out[[paste0("fitted_power_at_effect_size_", label)]] <-
    stats::pnorm(out$k * (-log(1 - es)) - z)
}

column_order <- c("grna_target", "response_id", "n_points", "k", "se_k", "deviance", "df",
                 "status", "min_detectable_effect_size", "min_detectable_effect_size_ci_low",
                 "min_detectable_effect_size_ci_high")
out <- out[, c(column_order, setdiff(colnames(out), column_order)), drop = FALSE]
out <- out[order(out$min_detectable_effect_size, na.last = TRUE), , drop = FALSE]
rownames(out) <- NULL

## WRITE ===========================================================================================

write_tsv_file(out, opts$out)
log_step("Wrote ", nrow(out), " pairs to ", opts$out)

message(sprintf("  z used                        %.4f", z))
message(sprintf("  fitted                        %d of %d pairs",
                sum(out$status == 0), nrow(out)))
message(sprintf("  saturated at every point      %d", sum(out$status == 1)))
message(sprintf("  zero successes everywhere     %d", sum(out$status == 2)))
message(sprintf("  fit failed                    %d", sum(out$status == 3)))
message(sprintf("  reaching power >= %.2f below a complete knockdown: %d (%.0f%%)",
                opts$power_threshold, sum(!is.na(out$min_detectable_effect_size)),
                100 * mean(!is.na(out$min_detectable_effect_size))))

# Residual deviance against its degrees of freedom is the model check the third effect size buys.
# Much above 1 per degree of freedom means the straight line on the probit scale does not describe
# these curves, and the fit should not be trusted in place of a measured sweep.
ok <- out$status == 0 & !is.na(out$deviance) & out$df > 0
if (any(ok)) {
  ratio <- out$deviance[ok] / out$df[ok]
  message(sprintf("  deviance per df: median %.2f, 90th %.2f  (>> 1 means the model does not fit)",
                  stats::median(ratio), stats::quantile(ratio, 0.9)))
}
