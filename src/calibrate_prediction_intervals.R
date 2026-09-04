#!/usr/bin/env Rscript
#
# Section 6, the part that decides whether prediction is publishable: can a predicted power carry an
# honest interval, and is that interval narrow enough to be worth anything?
#
# THE QUESTION, AND WHY COVERAGE IS THE TEST
#
# The covariate model predicts k to within about x1.2, which is +-0.16 in power at the median and
# +-0.30 at p90. A point prediction with that error is not a per-pair claim. A prediction plus an
# honest interval is -- but only if the interval actually contains the truth as often as it says.
#
# So this calibrates intervals on measured pairs and then checks COVERAGE on pairs the calibration
# never saw: does the nominal 90 % interval contain the measured power for ~90 % of held-out pairs?
# Coverage alone is not enough, because [0, 1] covers everything and says nothing, so width is
# reported alongside -- and then the only question that matters operationally:
#
#   HOW MANY PAIRS COULD BE CERTIFIED FROM A PREDICTION, VERSUS FROM MEASUREMENT?
#
# Certifying a negative means asserting power was at least 0.8, which from an interval means its
# LOWER bound clears 0.8. If prediction certifies a useful fraction of what measurement does, section
# 6 stands on its own. If it certifies almost nothing, prediction is for `trans` pairs and design
# counterfactuals, and measurement remains the per-pair method. Either answer is publishable; only
# the first one supports the stronger claim.
#
# WHY THE CALIBRATION IS CONDITIONAL ON EXPRESSION
#
# Measured on day0: the residual width the model needs varies 8.6x across expression deciles, and
# that survives both the log-linear surrogate and the theoretically-derived form. So a single global
# width is the wrong shape -- it would over-cover mid-expression pairs and under-cover the tails,
# and under-coverage is the dangerous direction: it certifies negatives that were never powered.
# Both are computed here so the difference is measured rather than argued.
#
# WHAT WOULD INVALIDATE THIS
#
# Calibrating and evaluating on the same pairs guarantees coverage and proves nothing, so the split
# is strict: the model, the strata boundaries and the residual quantiles all come from the training
# pairs only. Nothing about the held-out pairs informs the interval they are scored against.
#
# Usage:
#   calibrate_prediction_intervals.R --curves six_point_100.tsv --reference power_summary.tsv \
#     --threshold-file discovery_threshold.txt --sim-input sim_input.rds \
#     --out-summary intervals.txt --out-pairs interval_pairs.tsv

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
              help = "Six-point 100-replicate fitted curves, providing the target k."),
  make_option("--reference", type = "character", default = NULL, dest = "reference",
              help = "power_summary.tsv: covariates and the measured power to score against."),
  make_option("--threshold-file", type = "character", default = NULL, dest = "threshold_file",
              help = "discovery_threshold.txt, for z."),
  make_option("--sim-input", type = "character", default = NULL, dest = "sim_input",
              help = paste("sim_input.rds. Supplied, the theoretically-derived model is used",
                           "(log(1/mu + 1/theta)), which fits better than the log-linear",
                           "surrogate; omitted, the surrogate is used.")),
  make_option("--out-summary", type = "character", default = NULL, dest = "out_summary",
              help = "Text report: coverage, width, and certification rates."),
  make_option("--out-pairs", type = "character", default = NULL, dest = "out_pairs",
              help = "Held-out pairs with their predicted intervals."),
  make_option("--conf-level", type = "double", default = 0.9, dest = "conf_level",
              help = "Nominal interval coverage [default %default]."),
  make_option("--n-strata", type = "integer", default = 10L, dest = "n_strata",
              help = "Expression strata for conditional calibration [default %default]."),
  make_option("--holdout-frac", type = "double", default = 0.3, dest = "holdout_frac",
              help = "Fraction held out [default %default]."),
  make_option("--power-threshold", type = "double", default = 0.8, dest = "power_threshold",
              help = "Power level a certified negative must clear [default %default]."),
  make_option("--seed", type = "integer", default = 20250812L, dest = "seed",
              help = "Seed for the split [default %default].")
)
opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Calibrate and validate prediction intervals for per-pair power."
))
require_options(opts, c("curves", "reference", "threshold_file", "out_summary", "out_pairs"))

ALL_ES <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)

threshold <- as.numeric(readLines(opts$threshold_file, warn = FALSE)[1])
z <- stats::qnorm(1 - threshold)
log_step("z = ", format(z, digits = 6))

## LOAD ============================================================================================

curves <- read_tsv_file(opts$curves, required_columns = c("grna_target", "response_id", "k"))
reference <- read_tsv_file(opts$reference,
                           required_columns = c("grna_target", "response_id", "mean_pert_cells",
                                                "average_expression_all_cells"))
idx <- match(paste(curves$grna_target, curves$response_id, sep = "\r"),
             paste(reference$grna_target, reference$response_id, sep = "\r"))

df <- data.frame(
  grna_target = curves$grna_target, response_id = curves$response_id,
  k = curves$k,
  pert_cells = reference$mean_pert_cells[idx],
  expression = reference$average_expression_all_cells[idx],
  stringsAsFactors = FALSE
)
for (es in ALL_ES) {
  df[[paste0("measured_", effect_label(es))]] <-
    reference[[paste0("power_at_effect_size_", effect_label(es))]][idx]
}

keep <- is.finite(df$k) & df$k > 0 & is.finite(df$pert_cells) & df$pert_cells > 0 &
  is.finite(df$expression) & df$expression > 0
df <- df[keep, , drop = FALSE]
df$log_k <- log(df$k)
df$log_pert <- log(df$pert_cells)
df$log_expr <- log(df$expression)

# The theoretically-derived predictor when dispersion is available, since it fits better and is not
# an arbitrary functional form. The residual STRUCTURE is the same either way, which is why the
# conditional calibration below is needed regardless of which is used.
if (!is.null(opts$sim_input)) {
  rd <- readRDS(opts$sim_input)$row_data
  gi <- match(df$response_id, rownames(rd))
  ok <- is.finite(rd$mean[gi]) & rd$mean[gi] > 0 & is.finite(rd$dispersion[gi])
  df <- df[ok, , drop = FALSE]
  gi <- gi[ok]
  df$predictor <- log(1 / rd$mean[gi] + rd$dispersion[gi])
  model_name <- "theoretical: log(k) ~ log(pert) + log(1/mu + 1/theta)"
} else {
  df$predictor <- df$log_expr
  model_name <- "surrogate: log(k) ~ log(pert) + log(expression)"
}
log_step(nrow(df), " pairs; model = ", model_name)

## SPLIT, THEN CALIBRATE ON TRAIN ONLY =============================================================

set.seed(opts$seed)
test_i <- sample(nrow(df), round(opts$holdout_frac * nrow(df)))
train <- df[-test_i, , drop = FALSE]
test <- df[test_i, , drop = FALSE]

fit <- stats::lm(log_k ~ log_pert + predictor, data = train)
train$resid <- stats::residuals(fit)
test$pred_log_k <- stats::predict(fit, newdata = test)

alpha <- (1 - opts$conf_level) / 2
qs <- c(alpha, 1 - alpha)

# Global calibration: one pair of residual quantiles for everything.
global_q <- stats::quantile(train$resid, qs, names = FALSE)

# Conditional calibration: strata on expression, with boundaries taken from TRAIN only. Test pairs
# outside the training range are clamped into the end strata rather than dropped -- a real
# prediction has to answer for them, and dropping them would flatter the coverage.
breaks <- stats::quantile(train$log_expr, probs = seq(0, 1, length.out = opts$n_strata + 1),
                          names = FALSE)
breaks[1] <- -Inf
breaks[length(breaks)] <- Inf
train$stratum <- cut(train$log_expr, breaks = breaks, labels = FALSE)
test$stratum <- cut(test$log_expr, breaks = breaks, labels = FALSE)

strat_q <- t(vapply(split(train$resid, train$stratum),
                    function(v) stats::quantile(v, qs, names = FALSE), numeric(2)))
rownames(strat_q) <- names(split(train$resid, train$stratum))

## SCORE ON HELD-OUT PAIRS =========================================================================

power_at <- function(k, es) stats::pnorm(k * (-log(1 - es)) - z)

# CALIBRATE IN POWER SPACE, NOT IN k SPACE.
#
# A first version of this calibrated the residuals of log k and then mapped the resulting k interval
# through the curve. That answers the wrong question. An interval on k claims to cover the true k;
# what a per-pair claim needs is an interval that covers the MEASURED power, and three distinct
# errors sit between those: the model's error in k, the curve's own misfit (deviance/df reaches 9.5
# at high expression), and the binomial noise in the measurement itself.
#
# It also fails structurally at saturated effect sizes. At 0.5 nearly every pair has power near 1,
# so a x1.17 range in k maps to a power interval of width 0.000 -- which cannot cover a measured
# value that scatters, no matter how well k was calibrated. Measured coverage there was 85 % against
# a nominal 90 % for reasons that had nothing to do with the model.
#
# So the residual that gets calibrated is `measured power - predicted power`, per effect size, which
# absorbs all three error sources by construction. Quantiles still come from the training pairs
# only, so held-out coverage remains a real test rather than a tautology.
power_residuals <- function(rows, es) {
  k_hat <- exp(rows$pred_log_k)
  rows[[paste0("measured_", effect_label(es))]] - power_at(k_hat, es)
}
train$pred_log_k <- stats::fitted(fit)

score <- function(mode, label) {
  rows <- lapply(ALL_ES, function(es) {
    meas <- test[[paste0("measured_", effect_label(es))]]
    p_hat <- power_at(exp(test$pred_log_k), es)
    tr <- power_residuals(train, es)
    if (identical(mode, "global")) {
      q <- stats::quantile(tr, qs, na.rm = TRUE, names = FALSE)
      lo_off <- rep(q[1], nrow(test)); hi_off <- rep(q[2], nrow(test))
    } else {
      per <- t(vapply(split(tr, train$stratum),
                      function(v) stats::quantile(v, qs, na.rm = TRUE, names = FALSE), numeric(2)))
      key <- as.character(test$stratum)
      lo_off <- per[key, 1]; hi_off <- per[key, 2]
    }
    # Clamped to [0, 1]: a power interval reaching outside that is not wrong, but reporting it
    # would overstate the width and understate how often the bound is uninformative.
    p_lo <- pmax(0, pmin(1, p_hat + lo_off))
    p_hi <- pmax(0, pmin(1, p_hat + hi_off))
    inside <- meas >= p_lo & meas <= p_hi
    data.frame(calibration = label, effect_size = es,
               coverage = mean(inside, na.rm = TRUE),
               median_width = stats::median(p_hi - p_lo, na.rm = TRUE),
               # Certification from the prediction: the LOWER bound must clear the power
               # threshold. Compared against what measurement itself certifies on the same pairs.
               certified_pred = mean(p_lo >= opts$power_threshold, na.rm = TRUE),
               certified_meas = mean(meas >= opts$power_threshold, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

global <- score("global", "global")
cond <- score("conditional", "conditional")

# Per-stratum coverage under global calibration is the direct evidence for going conditional: if a
# single width were adequate, coverage would be flat across strata.
per_stratum <- do.call(rbind, lapply(sort(unique(test$stratum)), function(s) {
  sel <- test$stratum == s
  tr <- power_residuals(train, 0.15)
  gq <- stats::quantile(tr, qs, na.rm = TRUE, names = FALSE)
  cq <- stats::quantile(tr[train$stratum == s], qs, na.rm = TRUE, names = FALSE)
  p_hat <- power_at(exp(test$pred_log_k[sel]), 0.15)
  meas <- test[[paste0("measured_", effect_label(0.15))]][sel]
  covg <- function(q) mean(meas >= pmax(0, pmin(1, p_hat + q[1])) &
                             meas <= pmax(0, pmin(1, p_hat + q[2])), na.rm = TRUE)
  data.frame(stratum = s, n = sum(sel),
             coverage_global = covg(gq), coverage_conditional = covg(cq),
             stringsAsFactors = FALSE)
}))

## WRITE ===========================================================================================

out_pairs <- test[, c("grna_target", "response_id", "k", "pert_cells", "expression",
                      "pred_log_k", "stratum")]
for (es in ALL_ES) {
  tr <- power_residuals(train, es)
  per <- t(vapply(split(tr, train$stratum),
                  function(v) stats::quantile(v, qs, na.rm = TRUE, names = FALSE), numeric(2)))
  key <- as.character(test$stratum)
  p_hat <- power_at(exp(test$pred_log_k), es)
  out_pairs[[paste0("pred_power_", effect_label(es))]] <- p_hat
  out_pairs[[paste0("pred_power_lo_", effect_label(es))]] <- pmax(0, pmin(1, p_hat + per[key, 1]))
  out_pairs[[paste0("pred_power_hi_", effect_label(es))]] <- pmax(0, pmin(1, p_hat + per[key, 2]))
  out_pairs[[paste0("measured_", effect_label(es))]] <- test[[paste0("measured_", effect_label(es))]]
}
write_tsv_file(out_pairs, opts$out_pairs)

fmt <- function(d) {
  vapply(seq_len(nrow(d)), function(i) sprintf(
    "    es %-5s  coverage %5.1f%%  median width %.3f  certified: pred %5.1f%% vs measured %5.1f%%",
    d$effect_size[i], 100 * d$coverage[i], d$median_width[i],
    100 * d$certified_pred[i], 100 * d$certified_meas[i]), character(1))
}

lines <- c(
  "CALIBRATED PREDICTION INTERVALS FOR PER-PAIR POWER",
  "",
  sprintf("model: %s", model_name),
  sprintf("pairs: %d train / %d held out (seed %d)", nrow(train), nrow(test), opts$seed),
  sprintf("nominal coverage: %.0f%%   expression strata: %d",
          100 * opts$conf_level, opts$n_strata),
  sprintf("log-k residual quantiles (informational): %+.4f to %+.4f",
          global_q[1], global_q[2]),
  "",
  "GLOBAL calibration (one width for every pair):",
  fmt(global),
  "",
  "CONDITIONAL calibration (width per expression stratum):",
  fmt(cond),
  "",
  sprintf("Per-stratum coverage at effect size 0.15 (nominal %.0f%%):", 100 * opts$conf_level),
  "    stratum      n   global  conditional",
  sprintf("    %7d %6d  %5.1f%%      %5.1f%%", per_stratum$stratum, per_stratum$n,
          100 * per_stratum$coverage_global, 100 * per_stratum$coverage_conditional),
  "",
  "HOW TO READ THIS",
  "  Coverage far below nominal means the interval lies: it would certify negatives that were",
  "  never powered. Coverage far above nominal with a wide interval means it is honest but",
  "  useless. The certification columns are the operational answer -- what fraction of pairs a",
  "  PREDICTION can certify, against what MEASUREMENT certifies on the same pairs. If those are",
  "  close, prediction stands alone; if prediction certifies far fewer, it is for trans pairs and",
  "  design counterfactuals and measurement stays the per-pair method.",
  ""
)
writeLines(lines, opts$out_summary)
message(paste(lines, collapse = "\n"))
log_step("Wrote ", opts$out_summary, " and ", opts$out_pairs)
