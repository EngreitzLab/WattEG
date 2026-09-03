# materialize_odm_response_matrix() and attach_response_odm() have no real sceptre object to test
# against in CI, so they are exercised directly here against a small synthetic odm.

make_test_odm <- function(mat) {
  skip_if_not_installed("ondisc")
  h5_path <- withr::local_tempfile(fileext = ".h5", .local_envir = parent.frame())
  ondisc::create_odm_from_r_matrix(methods::as(mat, "CsparseMatrix"), h5_path,
                                    chunk_size = min(5L, ncol(mat)))
}

# odm only preserves row (gene) names -- cell barcodes live in @covariate_data_frame elsewhere in
# the sceptre object, not in the response matrix's dimnames -- so every matrix built for these
# tests needs rownames (create_odm_from_r_matrix requires them) but never colnames.
random_count_matrix <- function(n_gene, n_cell, density = 0.3, seed = 1) {
  withr::with_seed(seed, {
    m <- matrix(0L, n_gene, n_cell)
    n_nonzero <- round(n_gene * n_cell * density)
    idx <- sample.int(n_gene * n_cell, n_nonzero)
    m[idx] <- sample.int(20, n_nonzero, replace = TRUE)
    rownames(m) <- paste0("gene", seq_len(n_gene))
    m
  })
}

# A minimal sceptre_object that passes read_sceptre_object()'s earlier checks (functs_called,
# grna_integration_strategy) so attach_response_odm() is reached. @integer_id defaults to the
# odm's own, i.e. "matches" -- tests that want a mismatch override it explicitly.
make_test_sceptre_object <- function(response_odm) {
  skip_if_not_installed("sceptre")
  so <- methods::new("sceptre_object")
  so@functs_called <- c(assign_grnas = TRUE, run_qc = TRUE)
  so@grna_integration_strategy <- "union"
  so@integer_id <- response_odm@integer_id
  so@response_matrix <- list(response_odm)
  so
}

test_that("materialize_odm_response_matrix reproduces the source matrix exactly", {
  mat <- random_count_matrix(12, 8)
  odm <- make_test_odm(mat)
  materialized <- materialize_odm_response_matrix(odm)

  expect_equal(unname(as.matrix(materialized)), unname(mat))
  expect_equal(rownames(materialized), rownames(mat))
})

test_that("materialize_odm_response_matrix handles all-zero rows and columns", {
  # compute_expression_stats() already handles rows/columns with no nonzero entries (that is what
  # usable_gene / log_geomeans == -Inf is for); this just checks the odm read path does not drop
  # or misplace them on the way in.
  mat <- random_count_matrix(10, 6, density = 0.5)
  mat[3, ] <- 0L  # all-zero gene
  mat[, 5] <- 0L  # all-zero cell

  odm <- make_test_odm(mat)
  materialized <- materialize_odm_response_matrix(odm)

  expect_equal(unname(as.matrix(materialized)), unname(mat))
  expect_equal(unname(Matrix::rowSums(materialized)[3]), 0)
  expect_equal(unname(Matrix::colSums(materialized)[5]), 0)
})

test_that("project_and_check_odm_size passes for a small matrix and stops above the limit", {
  odm <- make_test_odm(random_count_matrix(20, 10, density = 0.5))

  expect_silent(project_and_check_odm_size(odm, limit_gb = ODM_MATERIALIZATION_LIMIT_GB))
  expect_error(
    project_and_check_odm_size(odm, limit_gb = 0),
    "projected to need"
  )
})

test_that("read_sceptre_object requires --response-odm for an odm-backed object", {
  odm <- make_test_odm(random_count_matrix(6, 4))
  so <- make_test_sceptre_object(odm)
  rds_path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(so, rds_path)

  expect_error(read_sceptre_object(rds_path), "--response-odm")
})

test_that("read_sceptre_object rejects a --response-odm with a mismatched integer_id", {
  odm <- make_test_odm(random_count_matrix(6, 4))
  other_odm <- make_test_odm(random_count_matrix(5, 3))  # a different backing file/integer_id
  so <- make_test_sceptre_object(odm)
  so@integer_id <- other_odm@integer_id + 1L  # guaranteed not to match either odm below

  rds_path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(so, rds_path)
  odm_path <- odm@h5_file

  expect_error(
    read_sceptre_object(rds_path, response_odm_fp = odm_path),
    "distinct IDs"
  )
})

test_that("read_sceptre_object reconnects a matching --response-odm end to end", {
  # This is the failure mode the whole fix exists for: odm's C++ external pointer does not
  # survive a saveRDS()/readRDS() round trip, even within the same R session, so reading a
  # just-deserialized odm without reconnecting first errors with "external pointer is not valid".
  mat <- random_count_matrix(8, 5)
  odm <- make_test_odm(mat)
  so <- make_test_sceptre_object(odm)

  rds_path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(so, rds_path)
  odm_path <- odm@h5_file

  reloaded <- read_sceptre_object(rds_path, response_odm_fp = odm_path)
  result <- get_response_matrix(reloaded)

  expect_false(methods::is(result, "odm"))
  expect_equal(unname(as.matrix(result)), unname(mat))
})

# sceptre_analysis_mode() reads which test a screen was configured for. There is no
# resampling_mechanism slot -- set_analysis_parameters() collapses that string to the boolean
# @run_permutations -- so these pin the mapping in both directions, plus the MOI reporting that
# explains how a screen reaches permutations without anyone asking for it.

test_that("sceptre_analysis_mode reports crt when run_permutations is FALSE", {
  skip_if_not_installed("sceptre")
  so <- methods::new("sceptre_object")
  so@run_permutations <- FALSE
  so@low_moi <- FALSE

  mode <- sceptre_analysis_mode(so)
  expect_false(mode$run_permutations)
  expect_equal(mode$resampling_mechanism, "crt")
  expect_equal(mode$moi, "high")
})

test_that("sceptre_analysis_mode reports permutations when run_permutations is TRUE", {
  skip_if_not_installed("sceptre")
  so <- methods::new("sceptre_object")
  so@run_permutations <- TRUE
  so@low_moi <- TRUE

  mode <- sceptre_analysis_mode(so)
  expect_true(mode$run_permutations)
  expect_equal(mode$resampling_mechanism, "permutations")
  expect_equal(mode$moi, "low")
})

test_that("format_analysis_mode emits one key-value line per setting", {
  skip_if_not_installed("sceptre")
  so <- methods::new("sceptre_object")
  so@run_permutations <- TRUE
  so@low_moi <- TRUE

  lines <- format_analysis_mode(sceptre_analysis_mode(so))
  expect_equal(lines, c("resampling_mechanism\tpermutations",
                        "run_permutations\ttrue",
                        "moi\tlow"))
})
