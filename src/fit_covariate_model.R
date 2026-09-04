#!/usr/bin/env Rscript
#
# Section 6 of paper/outline.md: predict a pair's power curve from covariates instead of measuring
# it, which is the only route to pairs that were never simulated.
#
# The theory says SE^2 ~ (1 / n_pert_cells) * (1 / mu + 1 / theta), so k = 1/SE should be
# predictable from the perturbed-cell count and the gene's expression, both of which the pipeline
# already reports per pair. This fits that, and then does the two things docs/status.md says have to
# happen before it can carry a paper.
#
# WHY A HELD-OUT SPLIT AND NOT JUST R^2
#
# R^2 on the pairs the model was fitted on says nothing about a pair the model has never seen, which
# is the entire use case. So the pairs are split, the model is fitted on one part, and the error is
# reported on the other. Held-out pairs are the cheap, in-regime version of the question -- `trans`
# pairs are the expensive, out-of-regime version, and good held-out performance here does NOT
# license the `trans` claim, because `trans` genes are not drawn from the element's neighbourhood
# and expression is the covariate carrying most of the signal.
#
# WHY THE RESIDUALS ARE BROKEN DOWN BY DECILE
#
# This is figures.md's Figure 5c, and it is not decoration. A x1.28 error in k is about +-0.20 in
# power mid-transition, and it is MODEL error, not sampling noise: it does not shrink with more
# simulation. docs/status.md's concern is that it is "likely concentrated in particular genes rather
# spread evenly", and if so:
#
#   * one global interval width is the wrong shape -- too narrow for some pairs, too wide for others
#   * for scE2G / ENCODE-rE2G training labels, mislabelled negatives concentrated in a non-random
#     subset are worse for a classifier than uniform noise
#
# So the deliverable is not the residual sd. It is whether the residual sd and bias move across
# deciles of expression and of perturbed-cell count. Reporting the global number alone would hide
# exactly the thing that matters.
#
# WHY THE INTERVALS ARE CALIBRATED, NOT TAKEN FROM THE REGRESSION
#
# A prediction interval from the fit assumes the residuals are sampling noise. They are not. So the
# honest version calibrates on measured pairs: take the empirical residual quantiles and widen each
# prediction by them. That gives intervals that are merely wider than measured ones rather than
# falsely narrow -- and quantifies what a measured sweep is buying, which is the point.
#
# Usage:
#   fit_covariate_model.R --curves six_point_100.tsv --reference power_summary.tsv \
#     --out-summary covariate_model.txt --out-residuals covariate_residuals.tsv

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

option_list <- list(
  make_option("--curves", type = "character", default = NULL, dest = "curves",
              help = paste("Fitted curves with a `k` column. Use the six-point, 100-replicate fit:",
                           "k is the thing being predicted, so it should be the best available",
                           "estimate, not one from a reduced design.")),
  make_option("--reference", type = "character", default = NULL, dest = "reference",
              help = "power_summary.tsv, for mean_pert_cells and average_expression_all_cells."),
  make_option("--threshold-file", type = "character", default = NULL, dest = "threshold_file",
              help = paste("The sweep's discovery_threshold.txt. Needed to turn an error in k into",
                           "an error in power: z = qnorm(1 - threshold) sets how steep every",
                           "pair's curve is at its transition.")),
  make_option("--out-summary", type = "character", default = NULL, dest = "out_summary",
              help = "Text report: coefficients, held-out error, per-decile residuals."),
  make_option("--out-residuals", type = "character", default = NULL, dest = "out_residuals",
              help = "Per-pair TSV: covariates, fitted and predicted k, residual, decile."),
  make_option("--holdout-frac", type = "double", default = 0.2, dest = "holdout_frac",
              help = "Fraction of pairs held out [default %default]."),
  make_option("--seed", type = "integer", default = 20250812L, dest = "seed",
              help = "Seed for the held-out split [default %default].")
)
opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Fit and validate the covariate model for the per-pair power curve."
))
require_options(opts, c("curves", "reference", "threshold_file", "out_summary",
                        "out_residuals"))

# Same derivation as fit_power_curve.R:176, which is where the curves' own k came from.
threshold <- as.numeric(readLines(opts$threshold_file, warn = FALSE)[1])
if (!is.finite(threshold) || threshold <= 0 || threshold >= 1) {
  stop("Unusable threshold in ", opts$threshold_file, call. = FALSE)
}
z <- stats::qnorm(1 - threshold)
log_step("z = ", format(z, digits = 6))

## LOAD AND JOIN ===================================================================================

curves <- read_tsv_file(opts$curves, required_columns = c("grna_target", "response_id", "k"))
reference <- read_tsv_file(opts$reference,
                           required_columns = c("grna_target", "response_id", "mean_pert_cells",
                                                "average_expression_all_cells"))

idx <- match(paste(curves$grna_target, curves$response_id, sep = "\r"),
             paste(reference$grna_target, reference$response_id, sep = "\r"))
df <- data.frame(
  grna_target = curves$grna_target,
  response_id = curves$response_id,
  k = curves$k,
  pert_cells = reference$mean_pert_cells[idx],
  expression = reference$average_expression_all_cells[idx],
  stringsAsFactors = FALSE
)
if ("status" %in% colnames(curves)) df$status <- curves$status

log_step("Joined ", nrow(df), " pairs")

# Both covariates are logged, so both must be strictly positive; k must be a usable finite slope.
# Dropping these is not a modelling choice -- log(0) is undefined -- but the count is reported
# because a large drop would mean the model is being fitted on a biased subset.
usable <- is.finite(df$k) & df$k > 0 &
  is.finite(df$pert_cells) & df$pert_cells > 0 &
  is.finite(df$expression) & df$expression > 0
log_step("Usable for the model: ", sum(usable), " of ", nrow(df),
         " (dropped ", sum(!usable), ")")
df <- df[usable, , drop = FALSE]
if (nrow(df) < 100) {
  stop("Only ", nrow(df), " usable pairs; not enough to fit or validate.", call. = FALSE)
}

df$log_k <- log(df$k)
df$log_pert <- log(df$pert_cells)
df$log_expr <- log(df$expression)

## FIT, IN SAMPLE ==================================================================================

full_fit <- stats::lm(log_k ~ log_pert + log_expr, data = df)
co <- stats::coef(full_fit)
r2 <- summary(full_fit)$r.squared
resid_sd <- stats::sd(stats::residuals(full_fit))

# Each covariate alone, because the joint R^2 hides which one is doing the work. On wtc11 expression
# alone explained 83.7 % and perturbed cells alone 2.3 % -- not because cell count is unimportant
# (its exponent matches theory) but because it barely varies across pairs. Whether that holds here
# is a fact about this dataset, and a reviewer will ask.
r2_pert <- summary(stats::lm(log_k ~ log_pert, data = df))$r.squared
r2_expr <- summary(stats::lm(log_k ~ log_expr, data = df))$r.squared

## HELD-OUT ========================================================================================

set.seed(opts$seed)
test <- sample(nrow(df), max(1L, round(opts$holdout_frac * nrow(df))))
train_fit <- stats::lm(log_k ~ log_pert + log_expr, data = df[-test, , drop = FALSE])
pred_log_k <- stats::predict(train_fit, newdata = df[test, , drop = FALSE])
held_resid <- df$log_k[test] - pred_log_k

# A residual on the log scale is a multiplicative error in k, which is the interpretable quantity:
# exp(|resid|) is the factor by which k is wrong.
held_factor <- exp(stats::quantile(abs(held_resid), c(0.5, 0.9), names = FALSE))

# What a multiplicative error in k costs in POWER, which is what anyone actually reads. Derived
# rather than approximated, and worth writing out because the algebra is what makes the number
# interpretable:
#
#   power(es) = Phi(k * beta - z)
#
# Evaluate at each pair's OWN transition midpoint, where the curve is steepest and the error is
# therefore worst: power = 0.5 means k_true * beta = z, i.e. beta = z / k_true. Substituting a
# predicted k_pred = k_true * exp(residual),
#
#   power_pred = Phi(k_true * exp(r) * z / k_true - z) = Phi(z * (exp(r) - 1))
#
# so the error is |Phi(z * (exp(r) - 1)) - 0.5| and depends only on the residual and z. Note this is
# a worst case per pair, not an average over effect sizes: away from the midpoint the curve is
# flatter and the same error in k costs less power.
power_err <- abs(stats::pnorm(z * (exp(held_resid) - 1)) - 0.5)

## RESIDUAL STRUCTURE -- FIGURE 5c =================================================================

df$fitted_log_k <- stats::fitted(full_fit)
df$residual <- stats::residuals(full_fit)

decile_of <- function(x) {
  breaks <- stats::quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE)
  cut(x, breaks = unique(breaks), include.lowest = TRUE, labels = FALSE)
}
df$expr_decile <- decile_of(df$log_expr)
df$pert_decile <- decile_of(df$log_pert)

by_decile <- function(decile_col, label) {
  parts <- split(df$residual, df[[decile_col]])
  data.frame(
    covariate = label,
    decile = as.integer(names(parts)),
    n = vapply(parts, length, integer(1)),
    mean_residual = vapply(parts, mean, numeric(1)),
    sd_residual = vapply(parts, stats::sd, numeric(1)),
    # The 5th-95th residual range is the width a calibrated interval would need in this decile. If
    # it moves across deciles, one global width is the wrong shape.
    q05 = vapply(parts, function(v) stats::quantile(v, 0.05, names = FALSE), numeric(1)),
    q95 = vapply(parts, function(v) stats::quantile(v, 0.95, names = FALSE), numeric(1)),
    stringsAsFactors = FALSE
  )
}
deciles <- rbind(by_decile("expr_decile", "log_expression"),
                 by_decile("pert_decile", "log_pert_cells"))
rownames(deciles) <- NULL

# The headline for Figure 5c: does the interval width a decile needs vary across deciles? A ratio
# near 1 means one global calibration is defensible; well above 1 means it is not.
width <- deciles$q95 - deciles$q05
width_ratio <- tapply(width, deciles$covariate, function(w) max(w) / min(w))

## WRITE ===========================================================================================

write_tsv_file(df[, c("grna_target", "response_id", "k", "pert_cells", "expression",
                      "log_k", "fitted_log_k", "residual", "expr_decile", "pert_decile")],
               opts$out_residuals)

lines <- c(
  "COVARIATE MODEL FOR THE PER-PAIR POWER CURVE",
  "",
  sprintf("pairs fitted: %d", nrow(df)),
  "",
  "log(k) = a + b * log(perturbed cells) + c * log(expression)",
  sprintf("  a = %+.4f", co[["(Intercept)"]]),
  sprintf("  b = %+.4f   (theory: 0.500 -- the sqrt(n) law)", co[["log_pert"]]),
  sprintf("  c = %+.4f   (theory: 0 to 0.5)", co[["log_expr"]]),
  sprintf("  R^2 = %.4f, residual sd of log k = %.4f  ->  k within x%.3f",
          r2, resid_sd, exp(resid_sd)),
  "",
  "Each covariate alone (which one is doing the work):",
  sprintf("  log(perturbed cells) only: R^2 = %.4f", r2_pert),
  sprintf("  log(expression) only:      R^2 = %.4f", r2_expr),
  "",
  sprintf("HELD OUT (%.0f%% of pairs, n = %d, seed %d)",
          100 * opts$holdout_frac, length(test), opts$seed),
  sprintf("  median |error| in k: x%.3f", held_factor[1]),
  sprintf("  p90    |error| in k: x%.3f", held_factor[2]),
  sprintf("  implied |error| in power near the transition: median %.3f, p90 %.3f",
          stats::median(power_err), stats::quantile(power_err, 0.9, names = FALSE)),
  "",
  "RESIDUAL STRUCTURE (figures.md Figure 5c)",
  "  A single global interval width is only defensible if the width each decile needs is",
  "  roughly constant. Ratio of widest to narrowest 5-95% residual range:",
  sprintf("    by %s: %.2fx", names(width_ratio), width_ratio),
  "",
  "  Model error does not shrink with more replicates, so these widths are a floor on any",
  "  interval a prediction can honestly carry.",
  ""
)
writeLines(lines, opts$out_summary)
message(paste(lines, collapse = "\n"))
log_step("Wrote ", opts$out_summary, " and ", opts$out_residuals)
