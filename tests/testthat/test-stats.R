test_that("wilson_interval matches published values at the boundaries", {
  # The boundaries are the whole reason this is not the normal approximation, so they are what the
  # test pins. Reference values are the standard Wilson score interval at 95 %.
  zero <- wilson_interval(0, 100)
  expect_equal(zero$low, 0)
  expect_equal(zero$high, 0.0370, tolerance = 1e-3)

  full <- wilson_interval(100, 100)
  expect_equal(full$low, 0.9630, tolerance = 1e-3)
  expect_equal(full$high, 1)

  half <- wilson_interval(50, 100)
  expect_equal(half$low, 0.4038, tolerance = 1e-3)
  expect_equal(half$high, 0.5962, tolerance = 1e-3)
})

test_that("wilson_interval does not collapse where the normal approximation does", {
  # p +/- z*sqrt(p(1-p)/n) gives [0, 0] for 0 successes, which would assert that power is certainly
  # zero and let a pair be called underpowered with no uncertainty at all.
  expect_gt(wilson_interval(0, 100)$high, 0)
  expect_lt(wilson_interval(100, 100)$low, 1)
})

test_that("wilson_interval is vectorised and stays within [0, 1]", {
  res <- wilson_interval(c(0, 1, 50, 99, 100), 100)
  expect_length(res$low, 5)
  expect_true(all(res$low >= 0 & res$low <= 1))
  expect_true(all(res$high >= 0 & res$high <= 1))
  expect_true(all(res$low <= res$high))
})

test_that("wilson_interval widens as the confidence level rises", {
  narrow <- wilson_interval(50, 100, conf_level = 0.80)
  wide   <- wilson_interval(50, 100, conf_level = 0.99)
  expect_lt(wide$low, narrow$low)
  expect_gt(wide$high, narrow$high)
})

test_that("wilson_interval narrows as replicates increase", {
  expect_lt(
    diff(unlist(wilson_interval(500, 1000)[c("low", "high")])),
    diff(unlist(wilson_interval(50, 100)[c("low", "high")]))
  )
})

test_that("effect_label keeps whole percentages whole", {
  # The bug this pins: trimming trailing zeros unconditionally turned 0.2 into "2" and 0.5 into
  # "5", so a 20 % knockdown was reported in a column named power_at_effect_size_2.
  expect_equal(effect_label(0.2), "20")
  expect_equal(effect_label(0.5), "50")
  expect_equal(effect_label(0.15), "15")
  expect_equal(effect_label(0.05), "5")
  expect_equal(effect_label(0.1), "10")
})

test_that("effect_label keeps a real fractional part", {
  expect_equal(effect_label(0.125), "12_5")
  expect_equal(effect_label(0.025), "2_5")
})

test_that("effect_label emits no dots, which downstream readers treat as separators", {
  for (es in c(0.05, 0.1, 0.125, 0.15, 0.2, 0.25, 0.5)) {
    expect_false(grepl(".", effect_label(es), fixed = TRUE), info = paste("effect size", es))
  }
})

test_that("effect_label is distinct across a realistic sweep", {
  sweep <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)
  expect_equal(length(unique(vapply(sweep, effect_label, character(1)))), length(sweep))
})

# compare_mdes() scores whether a reduced sweep design reproduces the full design's minimum
# detectable effect size. The NA handling is the part worth pinning: NA means "never reached the
# power threshold at any tested effect size", which is an answer, not missing data.

test_that("compare_mdes counts exact and within-one-grid-step agreement", {
  grid <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)

  # identical -> everything agrees
  x <- c(0.05, 0.15, 0.5)
  expect_equal(compare_mdes(x, x, grid)$exact, 1)
  expect_equal(compare_mdes(x, x, grid)$within1, 1)

  # one grid step out: 0.1 vs 0.15 are adjacent, so not exact but within one
  res <- compare_mdes(c(0.1), c(0.15), grid)
  expect_equal(res$exact, 0)
  expect_equal(res$within1, 1)

  # two steps out: 0.05 vs 0.15
  res <- compare_mdes(c(0.05), c(0.15), grid)
  expect_equal(res$within1, 0)

  # "one step" is a GRID step, not one unit of effect size: 0.25 and 0.5 are adjacent on this grid
  # despite being 0.25 apart, while 0.15 and 0.25 are two steps despite being 0.1 apart.
  expect_equal(compare_mdes(c(0.25), c(0.5), grid)$within1, 1)
  expect_equal(compare_mdes(c(0.15), c(0.25), grid)$within1, 0)
})

test_that("compare_mdes treats NA as an answer, not as missing data", {
  grid <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)

  # both never reach the threshold -> they agree
  expect_equal(compare_mdes(NA_real_, NA_real_, grid)$exact, 1)
  expect_equal(compare_mdes(NA_real_, NA_real_, grid)$within1, 1)

  # one reaches it and the other does not -> disagreement, and not rescued by within1
  expect_equal(compare_mdes(NA_real_, 0.5, grid)$exact, 0)
  expect_equal(compare_mdes(NA_real_, 0.5, grid)$within1, 0)
  expect_equal(compare_mdes(0.5, NA_real_, grid)$within1, 0)

  # NA pairs are counted in the denominator: dropping them would score designs on the easier subset
  # of pairs that reached the threshold at all.
  res <- compare_mdes(c(0.1, NA, 0.2), c(0.1, 0.5, 0.2), grid)
  expect_equal(res$n, 3L)
  expect_equal(res$exact, 2 / 3)
})

test_that("compare_mdes rejects values that are not on the grid", {
  expect_error(compare_mdes(0.07, 0.1, c(0.05, 0.1)), "must come from")
  expect_error(compare_mdes(c(0.05, 0.1), 0.1, c(0.05, 0.1)), "same length")
})
