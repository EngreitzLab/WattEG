test_that("center_effect_size_matrix puts each gene's perturbed mean on its target effect size", {
  # Drawing per-guide effect sizes and clamping negatives to 0 biases the mean upwards, so the
  # centring step is what makes "effect size 0.15" mean a 15 % knockdown on average rather than
  # something slightly weaker.
  set.seed(1)
  n_genes <- 3
  n_cells <- 200
  gene_effect_sizes <- c(0.85, 0.90, 0.95)
  pert_status <- rep(c(1, 0), each = n_cells / 2)

  mat <- matrix(runif(n_genes * n_cells, 0.5, 1.2), nrow = n_genes)
  centred <- center_effect_size_matrix(mat, pert_status, gene_effect_sizes)

  expect_equal(rowMeans(centred[, pert_status == 1, drop = FALSE]), gene_effect_sizes)
  expect_equal(rowMeans(centred[, pert_status == 0, drop = FALSE]), rep(1, n_genes))
})

test_that("center_effect_size_matrix preserves shape and never goes negative", {
  set.seed(2)
  mat <- matrix(runif(4 * 50, 0, 0.2), nrow = 4)
  pert_status <- rep(c(1, 0), each = 25)
  centred <- center_effect_size_matrix(mat, pert_status, rep(0.1, 4))

  expect_equal(dim(centred), dim(mat))
  expect_true(all(centred >= 0))
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
