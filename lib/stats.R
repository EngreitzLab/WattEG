## Pure statistical and formatting helpers.
##
## These live here rather than inside the scripts that use them so they can be tested without
## invoking a command-line interface. Each is used by exactly one script today; that is not a reason
## to inline them, because both have been a source of silent, wrong output before -- see the
## comments below and docs/status.md, "Correctness fixes".

#' Wilson score interval for a binomial proportion.
#'
#' Preferred over the normal approximation `p +/- z*sqrt(p(1-p)/n)` because these estimates live at
#' the boundaries, where that approximation degenerates: 0 successes out of 100 gives the interval
#' [0, 0], asserting the power is certainly zero. Wilson gives [0, 0.037], which is what 100
#' replicates without a success actually supports.
#'
#' This matters beyond tidiness: `power_ci_low` is the column thresholded to certify a negative, so
#' an interval that collapses at the boundary would certify claims the data cannot support.
wilson_interval <- function(successes, n, conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  phat <- successes / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  halfwidth <- (z / denom) * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))
  list(
    low = pmax(0, center - halfwidth),
    high = pmin(1, center + halfwidth)
  )
}

#' Column-name fragment for an effect size: 0.15 -> "15", 0.125 -> "12_5".
#'
#' Trailing zeros are trimmed only when the value actually has a fractional part, so 0.125 becomes
#' "12_5" rather than "12.50". Trimming unconditionally is wrong and quietly so: it would turn 0.2
#' into "2" and 0.5 into "5", i.e. a column called power_at_effect_size_2 for a 20 % knockdown.
#' That bug shipped once; see docs/status.md.
#'
#' The decimal point then becomes an underscore (12.5 -> "12_5") to keep every column name
#' snake_case and free of dots, which several downstream readers treat as name separators.
effect_label <- function(effect_size) {
  formatted <- format(effect_size * 100, trim = TRUE, scientific = FALSE)
  if (grepl(".", formatted, fixed = TRUE)) {
    formatted <- sub("0+$", "", formatted)
    formatted <- sub("\\.$", "", formatted)
  }
  gsub(".", "_", formatted, fixed = TRUE)
}

#' Agreement between two minimum-detectable-effect-size vectors, on a shared effect-size grid.
#'
#' Used to ask whether a reduced sweep design reproduces the full design's MDES per pair, and --
#' with the reference bootstrapped against itself -- how much agreement the reference's own noise
#' allows in the first place. Both callers need the same rule, which is why it lives here rather
#' than inside the comparison script.
#'
#' `NA` means "never reached the power threshold at any tested effect size", which is a real answer
#' rather than missing data: two NAs agree, and an NA against a value disagrees. Treating NA as
#' missing and dropping those pairs would silently score the designs on an easier subset -- the
#' pairs that reached the threshold at all.
#'
#' `within1` allows a disagreement of one step ON THE GRID, not one unit of effect size, because the
#' grid is uneven (0.25 to 0.5 is a bigger jump than 0.05 to 0.1). A pair whose true MDES sits near
#' a grid boundary flips between neighbouring steps on noise alone, so exact agreement understates
#' how well a design does; reporting both is what separates "wrong" from "one step out".
#'
#' @param fitted,reference MDES vectors of equal length, with values drawn from `grid` or NA
#' @param grid the effect sizes, ascending
compare_mdes <- function(fitted, reference, grid = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)) {
  if (length(fitted) != length(reference)) {
    stop("fitted and reference must have the same length.", call. = FALSE)
  }
  n <- length(fitted)
  if (n == 0) {
    return(list(exact = NA_real_, within1 = NA_real_, n = 0L))
  }

  both_na <- is.na(fitted) & is.na(reference)
  both_val <- !is.na(fitted) & !is.na(reference)

  exact <- both_na
  exact[both_val] <- fitted[both_val] == reference[both_val]

  # Positions on the grid, so "one step" means one grid position rather than one unit.
  fi <- match(fitted, grid)
  ri <- match(reference, grid)
  if (any(!is.na(fitted) & is.na(fi)) || any(!is.na(reference) & is.na(ri))) {
    stop("MDES values must come from `grid` or be NA.", call. = FALSE)
  }
  within1 <- both_na
  within1[both_val] <- abs(fi[both_val] - ri[both_val]) <= 1

  list(exact = mean(exact), within1 = mean(within1), n = n)
}
