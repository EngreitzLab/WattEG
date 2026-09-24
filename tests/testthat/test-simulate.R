test_that("center_effect_size_matrix puts each gene's perturbed mean on its target effect size", {
  # This pin is what defines the estimand: power at a FIXED element effect (docs/methods.md).
  # Control cells come in at exactly 1 -- other guides have no effect on the tested genes -- and
  # must leave at exactly 1.
  set.seed(1)
  n_genes <- 3
  n_cells <- 200
  gene_effect_sizes <- c(0.85, 0.90, 0.95)
  pert_status <- rep(c(1, 0), each = n_cells / 2)

  mat <- matrix(1, nrow = n_genes, ncol = n_cells)
  mat[, pert_status == 1] <- runif(n_genes * n_cells / 2, 0.5, 1.2)
  centred <- center_effect_size_matrix(mat, pert_status, gene_effect_sizes)

  expect_equal(rowMeans(centred[, pert_status == 1, drop = FALSE]), gene_effect_sizes)
  expect_true(all(centred[, pert_status == 0] == 1))
})

test_that("center_effect_size_matrix preserves shape and never goes negative", {
  set.seed(2)
  pert_status <- rep(c(1, 0), each = 25)
  mat <- matrix(1, nrow = 4, ncol = 50)
  mat[, pert_status == 1] <- runif(4 * 25, 0, 0.2)
  centred <- center_effect_size_matrix(mat, pert_status, rep(0.1, 4))

  expect_equal(dim(centred), dim(mat))
  expect_true(all(centred >= 0))
  expect_equal(rowMeans(centred[, pert_status == 1, drop = FALSE]), rep(0.1, 4))
})

test_that("center_effect_size_matrix refuses control cells that are not exactly 1", {
  # A control cell off 1 means either the matrix is still in perturbed-then-control order (the
  # cell-order mask then calls perturbed cells "control") or control noise has crept back in.
  mat <- matrix(c(0.8, 0.9, 1, 1.02), nrow = 1)
  expect_error(center_effect_size_matrix(mat, c(1, 1, 0, 0), 0.85), "exactly 1")
})

test_that("create_effect_size_matrix returns genes x cells", {
  # Orientation is genes in rows, cells in columns. run_power_simulation.R relies on it directly --
  # `es_mat[, restore_cell_order]` reorders columns by cell, and center_effect_size_matrix subsets
  # columns by perturbation status. A transposed matrix would pair each cell with another gene's
  # effect, which is the class of bug that produced the size-factor shuffle.
  set.seed(3)
  n_cells <- 40
  gene_effect_sizes <- c(0.8, 0.9)
  # 0 = no guide, 1..2 = targeting guides, 3..4 = control guides
  grna_pert_status <- sample(0:4, n_cells, replace = TRUE)
  pert_guides <- c("g1", "g2")

  mat <- create_effect_size_matrix(grna_pert_status, pert_guides, gene_effect_sizes, guide_sd = 0.13)

  expect_equal(nrow(mat), length(gene_effect_sizes))
  expect_equal(ncol(mat), n_cells)
})

test_that("create_effect_size_matrix leaves guide-free cells unperturbed", {
  set.seed(4)
  # Cells 1 and 2 carry no guide; 3-4 carry targeting guides, 5-6 control guides.
  grna_pert_status <- c(0, 0, 1, 2, 3, 4)
  mat <- create_effect_size_matrix(grna_pert_status, c("g1", "g2"), c(0.8, 0.9), guide_sd = 0.13)

  # Row 1 of the guide table is the no-effect row, so a cell with no guide gets exactly 1 for
  # every gene. Cells are columns.
  expect_equal(unname(mat[, 1]), c(1, 1))
  expect_equal(unname(mat[, 2]), c(1, 1))
  expect_false(all(mat[, 3] == 1))
  # Cells 5-6 carry guides that belong to something else: no effect on these genes, exactly 1.
  expect_equal(unname(mat[, 5]), c(1, 1))
  expect_equal(unname(mat[, 6]), c(1, 1))
})

test_that("create_effect_size_matrix never returns a negative effect", {
  # guide_sd large enough that the normal draws go negative and must be clamped.
  set.seed(5)
  mat <- create_effect_size_matrix(sample(0:4, 100, replace = TRUE), c("g1", "g2"),
                                   c(0.5, 0.5), guide_sd = 2)
  expect_true(all(mat >= 0))
})

test_that("build_dispersion_vector inverts theta, in the order asked for", {
  precomp <- list(
    GENE_A = list(theta = 4),
    GENE_B = list(theta = 2),
    GENE_C = list(theta = 10)
  )
  disp <- build_dispersion_vector(precomp, c("GENE_C", "GENE_A"))

  expect_equal(names(disp), c("GENE_C", "GENE_A"))
  expect_equal(unname(disp), c(1 / 10, 1 / 4))
})

test_that("build_dispersion_vector errors rather than recycling when a gene is missing", {
  # The original bug: dispersions came back in a list column with NULL holes, and unlist() silently
  # shortened the vector so rnbinom recycled -- every gene after the first gap was simulated with
  # another gene's dispersion. A hard error is the point of this function.
  precomp <- list(GENE_A = list(theta = 4))
  expect_error(
    build_dispersion_vector(precomp, c("GENE_A", "GENE_MISSING")),
    "no entry in @response_precomputations"
  )
})

test_that("build_dispersion_vector errors on a non-finite dispersion", {
  expect_error(
    build_dispersion_vector(list(GENE_A = list(theta = 0)), "GENE_A"),
    "non-finite dispersion"
  )
  expect_error(
    build_dispersion_vector(list(GENE_A = list(theta = NA_real_)), "GENE_A"),
    "non-finite dispersion"
  )
})

test_that("baseline_expression's fitted model is exp(X %*% coefs), gene by gene", {
  # The quantity the whole simulation rests on: the expected count sceptre's own null model gives
  # each cell. Checked against a per-gene loop rather than against a second matrix expression, so
  # a transpose error in tcrossprod() cannot agree with itself.
  set.seed(11)
  n_genes <- 4
  n_cells <- 25
  p <- 3
  coefs <- matrix(rnorm(n_genes * p, sd = 0.3), nrow = n_genes,
                  dimnames = list(paste0("g", 1:n_genes), NULL))
  covariates <- cbind(1, matrix(rnorm(n_cells * (p - 1)), nrow = n_cells))

  x <- list(genes = rownames(coefs), cells = paste0("c", seq_len(n_cells)),
            fitted_coefs = coefs,
            row_data = data.frame(mean = runif(n_genes), dispersion = runif(n_genes),
                                  row.names = rownames(coefs)),
            col_data = data.frame(size_factors = runif(n_cells, 0.5, 2)))

  got <- baseline_expression(x, covariates, model = "fitted")
  expect_equal(dim(got), c(n_genes, n_cells))
  for (i in seq_len(n_genes)) {
    expect_equal(got[i, ], as.vector(exp(covariates %*% coefs[i, ])))
  }
})

test_that("baseline_expression's size_factor model is the outer product it always was", {
  set.seed(12)
  x <- list(genes = c("a", "b"), cells = paste0("c", 1:6),
            fitted_coefs = NULL,
            row_data = data.frame(mean = c(2, 5), dispersion = c(0.1, 0.2),
                                  row.names = c("a", "b")),
            col_data = data.frame(size_factors = c(0.5, 1, 1.5, 2, 2.5, 3)))

  expect_equal(baseline_expression(x, model = "size_factor"),
               outer(c(2, 5), c(0.5, 1, 1.5, 2, 2.5, 3)))
})

test_that("baseline_expression refuses a mismatched or missing design rather than guessing", {
  x <- list(genes = "a", cells = c("c1", "c2"),
            fitted_coefs = matrix(0.1, nrow = 1, ncol = 3, dimnames = list("a", NULL)),
            row_data = data.frame(mean = 1, dispersion = 0.1, row.names = "a"),
            col_data = data.frame(size_factors = c(1, 1)))

  expect_error(baseline_expression(x, model = "fitted"), "needs the covariate matrix")
  # Wrong number of cells, and wrong number of covariates: both would otherwise recycle or
  # broadcast into a plausible-looking matrix.
  expect_error(baseline_expression(x, matrix(0, nrow = 5, ncol = 3), model = "fitted"),
               "5 rows but there are 2 cells")
  expect_error(baseline_expression(x, matrix(0, nrow = 2, ncol = 2), model = "fitted"),
               "different designs")

  no_coefs <- x
  no_coefs$fitted_coefs <- NULL
  expect_error(baseline_expression(no_coefs, matrix(0, nrow = 2, ncol = 3), model = "fitted"),
               "carries no fitted_coefs")
})

test_that("draw_counts centres on the baseline it is handed, not on a normalised mean", {
  # With the effect size fixed at 1 and a large cell count, the simulated mean has to track the
  # baseline. This is what makes the two expression models distinguishable at all.
  set.seed(13)
  n_cells <- 4000
  baseline <- rbind(rep(5, n_cells), rep(0.5, n_cells))
  x <- list(genes = c("a", "b"), cells = paste0("c", seq_len(n_cells)),
            row_data = data.frame(mean = c(999, 999), dispersion = c(0.05, 0.05),
                                  row.names = c("a", "b")),
            col_data = data.frame(size_factors = rep(1, n_cells)))

  counts <- draw_counts(x, matrix(1, nrow = 2, ncol = n_cells), baseline)
  expect_equal(dim(counts), c(2L, n_cells))
  # row_data$mean is deliberately absurd: if it leaked into the draw this would fail loudly.
  expect_equal(unname(rowMeans(counts)), c(5, 0.5), tolerance = 0.05)
})

test_that("draw_counts rejects a baseline of the wrong shape", {
  x <- list(genes = c("a", "b"), cells = c("c1", "c2", "c3"),
            row_data = data.frame(mean = c(1, 1), dispersion = c(0.1, 0.1),
                                  row.names = c("a", "b")),
            col_data = data.frame(size_factors = rep(1, 3)))
  es <- matrix(1, nrow = 2, ncol = 3)
  expect_error(draw_counts(x, es, matrix(1, nrow = 2, ncol = 2)), "baseline is 2 x 2")
})

test_that("build_fitted_coefs_matrix reads the other half of the same fit", {
  precomps <- list(
    g1 = list(fitted_coefs = c(0.1, 0.2), theta = 5),
    g2 = list(fitted_coefs = c(0.3, 0.4), theta = 8),
    g3 = list(fitted_coefs = c(0.5, 0.6), theta = 2)
  )
  got <- build_fitted_coefs_matrix(precomps, c("g3", "g1"))
  expect_equal(rownames(got), c("g3", "g1"))
  expect_equal(got["g1", ], c(0.1, 0.2))
  expect_equal(got["g3", ], c(0.5, 0.6))

  # The same guard build_dispersion_vector has, for the same reason: a gene with no entry used to
  # become a NULL hole that unlist() dropped, shifting every later gene.
  expect_error(build_fitted_coefs_matrix(precomps, c("g1", "absent")),
               "no entry in @response_precomputations")
  precomps$g2$fitted_coefs <- c(0.3, NA)
  expect_error(build_fitted_coefs_matrix(precomps, c("g1", "g2")), "non-finite fitted coefficient")
})


## PROOF TESTS: the effect-size path the simulation actually runs ==================================
##
## Everything below drives the production functions (create_guide_pert_status ->
## simulate_effect_sizes) on the shared fixture in tests/fixtures/, which the Python suite reads
## byte for byte. The previous version of this block rebuilt the create -> reorder -> centre order by
## hand, so a regression in run_power_simulation.R would still have passed; it now calls the same
## simulate_effect_sizes() the simulation calls.
##
## What is being proven is the estimand decided 2026-09-24 (docs/methods.md): power at a FIXED
## element effect. In every replicate the realised mean effect across the perturbed cells equals the
## requested one; the target's guides still differ from each other; every other cell is exactly 1.

test_that("simulate_effect_sizes pins each gene's realised mean to the requested effect", {
  fx <- fixture_target()
  is_pert <- fx$pert_status == 1
  for (es in c(0.05, 0.15, 0.5)) {
    target <- c(g1 = 1 - es, g2 = 1 - es, g3 = 1 - es)
    set.seed(100)
    realised <- replicate(50, {
      m <- simulate_effect_sizes(fx$status, fx$pert_status, fx$restore, fx$pert_guides,
                                 target, guide_sd = 0.13)
      expect_true(all(m[, !is_pert] == 1))
      rowMeans(m[, is_pert, drop = FALSE])
    })
    expect_equal(unname(c(realised)), rep(1 - es, length(realised)), tolerance = 1e-12)
  }
})

test_that("centring in the wrong order is refused, and skipping it leaves the mean free to vary", {
  fx <- fixture_target()
  is_pert <- fx$pert_status == 1
  target <- c(g1 = 0.85)

  # The order every version ran until 2026-09-21, going back to the original DC-TAP code: centre
  # while the columns are still perturbed-then-control. The cell-order mask then labels some
  # perturbed cells "control", and those carry a targeting effect rather than 1. That used to shift
  # the wrong columns in silence; now it cannot run at all.
  set.seed(7)
  block_order <- create_effect_size_matrix(fx$status, fx$pert_guides, target, guide_sd = 0.13)
  expect_error(center_effect_size_matrix(block_order, fx$pert_status, target), "exactly 1")

  # Without centring the realised mean wobbles from replicate to replicate by about
  # guide_sd * sqrt(sum n_g^2) / sum n_g -- that is the other estimand, random guide effects.
  # Centring is what removes it.
  set.seed(8)
  uncentred <- replicate(200, {
    m <- create_effect_size_matrix(fx$status, fx$pert_guides, target, guide_sd = 0.13)
    mean(m[, fx$restore, drop = FALSE][1, is_pert])
  })
  expect_gt(sd(uncentred), 0.01)
})

test_that("the guide status is each cell's own guide, in cell order", {
  fx <- fixture_target()
  n_t <- length(fx$pert_guides)
  expect_length(fx$status, length(fx$cells))
  st <- fx$status[fx$restore]

  ok <- vapply(seq_along(fx$cells), function(j) {
    carried <- fx$carried[[j]]
    if (fx$pert_status[j] == 1) {
      st[j] %in% match(intersect(carried, fx$pert_guides), fx$pert_guides)
    } else if (length(carried) == 0) {
      st[j] == 0
    } else {
      st[j] > n_t && fx$other_guides[st[j] - n_t] %in% carried
    }
  }, logical(1))
  expect_true(all(ok), info = paste("wrong status for cells:",
                                    paste(utils::head(fx$cells[!ok]), collapse = ", ")))

  # Cells carrying the same guide get the same effect within a replicate.
  set.seed(13)
  m <- simulate_effect_sizes(fx$status, fx$pert_status, fx$restore, fx$pert_guides,
                             c(g1 = 0.85), guide_sd = 0.13)
  is_pert <- fx$pert_status == 1
  spread <- tapply(m[1, is_pert], st[is_pert], function(v) diff(range(v)))
  expect_true(all(spread == 0))
})

test_that("pinning the mean keeps the guide-to-guide spread", {
  # Centring adds one constant per gene; it must not flatten the guides into one value. With n_g
  # perturbed cells on guide g and N in total, the expected within-replicate variance across the
  # perturbed cells is guide_sd^2 * (1 - sum(n_g^2) / N^2).
  fx <- fixture_target()
  is_pert <- fx$pert_status == 1
  st <- fx$status[fx$restore][is_pert]
  n_g <- tabulate(st, nbins = length(fx$pert_guides))
  expected <- 0.13^2 * (1 - sum(n_g^2) / sum(n_g)^2)

  set.seed(9)
  within <- replicate(400, {
    v <- simulate_effect_sizes(fx$status, fx$pert_status, fx$restore, fx$pert_guides,
                               c(g1 = 0.85), guide_sd = 0.13)[1, is_pert]
    mean((v - mean(v))^2)
  })
  expect_equal(mean(within), expected, tolerance = 0.1)
})

test_that("no control cell can index a targeting row, even when the last target guide is unused", {
  fx <- fixture_target()
  # t8 is listed last and carried by no cell (tests/fixtures/README.md). Offsetting control statuses
  # by the highest target index any cell carried (7) put guide o001 on t8's targeting row.
  expect_false("t8" %in% unlist(fx$carried))
  st_ctrl <- fx$status[fx$restore][fx$pert_status == 0]
  expect_true(all(st_ctrl == 0 | st_ctrl > length(fx$pert_guides)))

  # With no guide-to-guide spread every effect is exact: perturbed cells at the target, all other
  # cells at 1. On the old offset some control cells came out at 0.5.
  m <- create_effect_size_matrix(fx$status, fx$pert_guides, c(g1 = 0.5), guide_sd = 0)
  m <- m[, fx$restore, drop = FALSE]
  expect_true(all(m[1, fx$pert_status == 1] == 0.5))
  expect_true(all(m[1, fx$pert_status == 0] == 1))
})

test_that("a lone perturbed cell carrying several target guides still gets one status per cell", {
  # One perturbed cell carrying both target guides: the subset used to drop to a vector, skip the
  # conversion, and hand back 7 statuses for 6 cells without an error.
  grna_perts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 3), j = c(2, 2, 1, 3, 5), x = 1, dims = c(4, 6),
    dimnames = list(c("t1", "t2", "o1", "o2"), paste0("c", 1:6))
  )
  pert_status <- c(0L, 1L, 0L, 0L, 0L, 0L)
  set.seed(10)
  status <- create_guide_pert_status(pert_status, grna_perts, c("t1", "t2"))

  expect_length(status, 6)
  in_cell_order <- status[order(cell_order(pert_status))]
  expect_true(in_cell_order[2] %in% c(1, 2))
  expect_equal(in_cell_order[c(1, 3, 5)], c(3, 4, 3))   # o1 -> 2 + 1, o2 -> 2 + 2
  expect_equal(in_cell_order[c(4, 6)], c(0, 0))
})

test_that("strong knockdowns are pinned too, and a pin that cannot be reached is an error", {
  fx <- fixture_target()
  is_pert <- fx$pert_status == 1
  for (es in c(0.7, 0.9)) {
    set.seed(11)
    realised <- replicate(30, {
      m <- simulate_effect_sizes(fx$status, fx$pert_status, fx$restore, fx$pert_guides,
                                 c(g1 = 1 - es, g2 = 1 - es), guide_sd = 0.13)
      expect_true(all(m >= 0))
      rowMeans(m[, is_pert, drop = FALSE])
    })
    expect_equal(unname(c(realised)), rep(1 - es, length(realised)), tolerance = 1e-12)
  }

  # One shift of -0.2 sends two cells below zero; clamping them leaves the mean at 0.233, not 0.1.
  # Allowed a single pass it must say so rather than return a matrix that misses the target.
  mat <- matrix(c(0, 0, 0.9, 1, 1), nrow = 1)
  expect_error(center_effect_size_matrix(mat, c(1, 1, 1, 0, 0), 0.1, max_iter = 1L),
               "Could not pin")
  expect_equal(mean(center_effect_size_matrix(mat, c(1, 1, 1, 0, 0), 0.1)[1, 1:3]), 0.1,
               tolerance = 1e-12)
})

test_that("baseline_expression refuses covariates in a different order than the coefficients", {
  coefs <- matrix(c(0.1, 0.5), nrow = 1, dimnames = list("a", c("(Intercept)", "log(n_umis)")))
  covariates <- cbind("(Intercept)" = 1, "log(n_umis)" = c(1, 2, 3))
  x <- list(genes = "a", cells = paste0("c", 1:3), fitted_coefs = coefs,
            row_data = data.frame(mean = 1, dispersion = 0.1, row.names = "a"),
            col_data = data.frame(size_factors = c(1, 1, 1)))

  expect_no_error(baseline_expression(x, covariates, model = "fitted"))
  expect_error(baseline_expression(x, covariates[, 2:1], model = "fitted"),
               "name their columns differently")

  # Each cell keeps its own covariates: permuting the cells permutes the baseline identically.
  perm <- c(3, 1, 2)
  expect_equal(baseline_expression(x, covariates[perm, , drop = FALSE], model = "fitted"),
               baseline_expression(x, covariates, model = "fitted")[, perm, drop = FALSE])
})

test_that("subset_genes refuses duplicate gene ids", {
  x <- list(genes = c("a", "b"), cells = "c1", fitted_coefs = NULL,
            row_data = data.frame(mean = c(1, 2), dispersion = c(0.1, 0.2),
                                  row.names = c("a", "b")))
  expect_error(subset_genes(x, c("a", "a")), "duplicate gene ids")
  expect_equal(subset_genes(x, c("b", "a"))$genes, c("b", "a"))
})
