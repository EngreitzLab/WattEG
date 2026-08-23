## Read the pipeline's single input -- a sceptre object -- and pull out everything it needs.
##
## This is the one file that touches sceptre's internal S4 slots. The pipeline is not using
## sceptre's public API here: it swaps the response matrix for simulated counts, narrows the
## discovery pairs to one perturbation at a time, and reads the cached negative-binomial
## precomputations. None of that is a documented interface, so bin/check_sceptre_api.R asserts
## every slot exists and the sceptre version is pinned. Keeping the access in one place means a
## sceptre upgrade breaks in one file rather than five.
##
## Everything the old pipeline took as a separate input file is derived here instead:
##   gene_grna_group_pairs.rds     -> qc_passing_pairs()      (@discovery_pairs_with_info)
##   grna_groups_table.rds         -> grna_target_map()       (@grna_target_data_frame)
##   discovery results (rds/csv)   -> discovery_threshold()   (@discovery_result)

suppressPackageStartupMessages(library(Matrix))

#' Read and validate a sceptre object.
#'
#' @param path path to the .rds
#' @param require_discovery whether @discovery_result must be populated (needed only to derive
#'   the p-value threshold)
#' @param response_odm_fp path to the response matrix's backing .odm file. Only needed when the
#'   object's response matrix is odm-backed (out-of-core); ignored otherwise. Required rather than
#'   read from the odm's own stored path for the same reason sceptre's own
#'   `read_ondisc_backed_sceptre_object()` requires it: the object may have been created on a
#'   different machine or moved since, and the stored path is not guaranteed to still resolve.
read_sceptre_object <- function(path, require_discovery = FALSE, response_odm_fp = NULL) {
  if (!file.exists(path)) {
    stop("Sceptre object not found: ", path, call. = FALSE)
  }
  so <- readRDS(path)

  if (!methods::is(so, "sceptre_object")) {
    stop(path, " is a ", paste(class(so), collapse = "/"), ", not a sceptre_object.", call. = FALSE)
  }

  # assign_grnas() populates the assignment lists the perturbation matrices are built from, and
  # run_qc() adds the pass_qc column the pairs are filtered on.
  for (funct in c("assign_grnas", "run_qc")) {
    if (!isTRUE(so@functs_called[[funct]])) {
      stop("`", funct, "()` has not been called on this sceptre object (functs_called is FALSE). ",
           "The power analysis needs it.", call. = FALSE)
    }
  }

  if (require_discovery && !isTRUE(so@functs_called[["run_discovery_analysis"]])) {
    stop("`run_discovery_analysis()` has not been called on this sceptre object, so there is no ",
         "discovery result to take the significance threshold from. Either supply an object that ",
         "has it, or pass --alpha explicitly.", call. = FALSE)
  }

  # `grna_group` in @discovery_pairs_with_info holds the perturbation *target* only under the
  # "union" integration strategy; under "singleton" it holds individual gRNA ids, and every pair
  # in this pipeline would then be silently mis-keyed.
  if (!identical(so@grna_integration_strategy, "union")) {
    stop("grna_integration_strategy is '", so@grna_integration_strategy,
         "' but this pipeline requires 'union': the target-level perturbation matrix and the ",
         "discovery pairs are both keyed by target.", call. = FALSE)
  }

  if (methods::is(so@response_matrix[[1]], "odm")) {
    so@response_matrix[[1]] <- attach_response_odm(so, response_odm_fp)
  }

  so
}

#' Reconnect an odm-backed response matrix to its backing file after a readRDS() round trip.
#'
#' odm objects hold a C++ external pointer to their backing HDF5 file, and external pointers are
#' never preserved across serialization -- readRDS() always hands back a dead one, even within the
#' same R process. Mirrors sceptre's own `read_ondisc_backed_sceptre_object()`: reconnect from a
#' freshly supplied path, then check `@integer_id` against the sceptre object's own, which is the
#' same defense sceptre uses against silently pairing an object with the wrong backing file.
attach_response_odm <- function(so, response_odm_fp) {
  if (is.null(response_odm_fp)) {
    stop("This sceptre object's response matrix is odm-backed (out-of-core), which needs the ",
         "backing .odm file's path to reconnect -- pass --response-odm.", call. = FALSE)
  }
  if (!requireNamespace("ondisc", quietly = TRUE)) {
    stop("This sceptre object's response matrix is odm-backed, which requires the `ondisc` ",
         "package. It is not installed.", call. = FALSE)
  }
  if (!file.exists(response_odm_fp)) {
    stop("--response-odm file not found: ", response_odm_fp, call. = FALSE)
  }

  response_odm <- ondisc::initialize_odm_from_backing_file(response_odm_fp)
  if (!identical(so@integer_id, response_odm@integer_id)) {
    stop("The sceptre object and --response-odm (", response_odm_fp, ") have distinct IDs; ",
         "the sceptre object was likely built from a different backing .odm file.", call. = FALSE)
  }
  response_odm
}

#' Cell barcodes, in the column order of the response matrix.
sceptre_cells <- function(so) {
  cells <- rownames(so@covariate_data_frame)
  if (is.null(cells)) {
    stop("@covariate_data_frame has no rownames, so cell barcodes cannot be recovered. ",
         "The response matrix has no column names either, so there is nothing to fall back on.",
         call. = FALSE)
  }
  cells
}

#' The discovery pairs that passed QC, as grna_target / response_id.
#'
#' Normalises sceptre's `grna_group` column name to `grna_target`, which is what the rest of the
#' pipeline and the output tables use. Filtering to pass_qc here rather than inside the
#' simulation means the splits are balanced over work that will actually run: on sample1 that is
#' 32,386 of 33,809 pairs.
qc_passing_pairs <- function(so) {
  pairs <- so@discovery_pairs_with_info

  key <- intersect(c("grna_group", "grna_target"), colnames(pairs))[1]
  if (is.na(key)) {
    stop("@discovery_pairs_with_info has neither a 'grna_group' nor a 'grna_target' column; ",
         "found: ", paste(colnames(pairs), collapse = ", "), ".", call. = FALSE)
  }
  if (!"pass_qc" %in% colnames(pairs)) {
    stop("@discovery_pairs_with_info has no 'pass_qc' column; has run_qc() been called?",
         call. = FALSE)
  }

  keep <- which(pairs$pass_qc)   # which() drops NA, which logical indexing would keep as NA rows
  out <- data.frame(
    grna_target = as.character(pairs[[key]][keep]),
    response_id = as.character(pairs$response_id[keep]),
    stringsAsFactors = FALSE
  )
  if (nrow(out) == 0) {
    stop("No discovery pairs passed QC, so there is nothing to run a power analysis on.",
         call. = FALSE)
  }
  out
}

#' gRNA id -> target mapping. Replaces the old grna_groups_table.rds input.
grna_target_map <- function(so) {
  gt <- so@grna_target_data_frame
  required <- c("grna_id", "grna_target")
  missing_cols <- setdiff(required, colnames(gt))
  if (length(missing_cols) > 0) {
    stop("@grna_target_data_frame is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". Found: ", paste(colnames(gt), collapse = ", "), ".", call. = FALSE)
  }
  data.frame(
    grna_id = as.character(gt$grna_id),
    grna_target = as.character(gt$grna_target),
    stringsAsFactors = FALSE
  )
}

#' The nominal p-value threshold that counts as a discovery.
#'
#' Taken from the real discovery results: the largest nominal p-value that survived multiple
#' testing correction. This is what the original compute_power_from_simulations.R did, except it
#' read it from a separate discovery-results file; @discovery_result already holds the same table
#' (verified against sample1's sceptre_discovery_results.csv: both give 7.71e-4).
#'
#' Note `which(significant)` rather than `significant == TRUE`: pairs that failed QC have
#' p_value = NA and significant = NA (1,423 of them on sample1), and plain logical indexing would
#' pull those NAs into the max().
discovery_threshold <- function(so) {
  dr <- so@discovery_result
  if (nrow(dr) == 0) {
    stop("@discovery_result is empty, so no significance threshold can be derived. ",
         "Pass --alpha to set one explicitly.", call. = FALSE)
  }
  for (column in c("p_value", "significant")) {
    if (!column %in% colnames(dr)) {
      stop("@discovery_result has no '", column, "' column.", call. = FALSE)
    }
  }

  significant_p <- dr$p_value[which(dr$significant)]
  if (length(significant_p) == 0) {
    stop("No pair in @discovery_result is significant, so the largest significant p-value is ",
         "undefined. The original code returned -Inf here and silently reported zero power for ",
         "every pair. Pass --alpha to set the threshold explicitly instead.", call. = FALSE)
  }

  threshold <- max(significant_p)
  if (!is.finite(threshold) || threshold <= 0 || threshold > 1) {
    stop("Derived p-value threshold ", threshold, " is not a usable probability.", call. = FALSE)
  }
  threshold
}

#' An empty dgRMatrix with the same shape and dimnames as `original`.
empty_like <- function(original) {
  new("dgRMatrix",
      j = integer(0),
      p = integer(nrow(original) + 1L),
      x = numeric(0),
      Dim = dim(original),
      Dimnames = list(rownames(original), colnames(original)))
}

#' Strip the two count matrices the power simulation never needs.
#'
#' Every parallel task loads this object, so its size is paid n_splits x n_effect_sizes times
#' over. Measured on sample1 (total object 1,774 MB):
#'
#'   @response_matrix  1,095 MB  overwritten with simulated counts on every rep
#'   @grna_matrix        417 MB  raw gRNA counts, only ever read by assign_grnas()
#'
#' Dropping @grna_matrix was verified empirically, not just by reading sceptre's source: running
#' run_discovery_analysis() on the same target under the same seed with and without it produces
#' an identical result table (all.equal TRUE on p_value, log_2_fold_change and every other
#' column). It is safe because assign_grnas() and run_qc() have already been called, so
#' skip_assign_grnas_and_run_qc() is a no-op and nothing else in the discovery path touches it.
#'
#' The shape and dimnames are preserved because sceptre computes n_cells from
#' ncol(get_grna_matrix()) in some code paths.
slim_sceptre_object <- function(so) {
  so@response_matrix <- list(empty_like(so@response_matrix[[1]]))
  if (length(so@grna_matrix) >= 1 && !methods::is(so@grna_matrix[[1]], "odm")) {
    so@grna_matrix <- list(empty_like(so@grna_matrix[[1]]))
  }

  # Drop the per-gene null models fitted during the *real* discovery analysis.
  #
  # sceptre skips fitting a gene's null model whenever @response_precomputations already holds an
  # entry for it -- exact memoization in a real analysis, where the gene's counts are identical in
  # every one of the ~147 pairs it appears in. In a simulation the counts are redrawn every replicate,
  # so carrying these entries over means testing simulated counts against coefficients fitted to real
  # counts. Measured, that understates power: 7 of 265 discovery calls flipped against a faithful
  # refit, all 7 in the same direction (docs/status.md).
  #
  # Cleared here rather than in prepare_sim_input.R because it must happen *after*
  # build_dispersion_vector() has read the slot -- the dispersions must keep coming from the real
  # data, since they set the noise the simulation exists to reproduce. Only the null model moves.
  # fit_null_models.R supplies the replacement, per replicate.
  so@response_precomputations <- list()

  so
}

## RESPONSE MATRIX ACCESS ==========================================================================
##
## The only place the pipeline reads actual expression values. Isolated here because it is also
## the only place that needs to change to support out-of-core (ondisc / `odm`) objects:
## everything downstream works from the per-gene means and per-cell size factors this produces,
## and the simulated matrix handed back to sceptre is always a small in-memory dgRMatrix.
##
## odm response matrices are materialized into an ordinary in-memory sparse matrix here, once, and
## nothing downstream (compute_expression_stats(), the CRT simulation) needs to change: the result
## is bit-identical to what they would compute had the sceptre object carried an in-memory matrix
## from the start. This works because the pipeline only ever needs *all* genes at once for the
## poscounts size factors (a per-cell reduction over the whole gene set), and that turns out to be
## small even at whole-transcriptome scale: ~1.8 GB sparse for 38,606 genes x 147,856 cells,
## measured on a real DC-TAP-seq object (2.6% density). ODM_MATERIALIZATION_LIMIT_GB exists as a
## backstop against the genuinely-500-GB case this comment used to warn about, not because the
## common case is expected to be close to it.
##
## odm's `[` only accepts a single row index, never a range -- there is no bulk row-chunk read --
## so this reads one gene at a time. It is still a one-time, ~2 ms/gene cost (measured), paid once
## per sample rather than once per simulation replicate.
##
## The odm is reconnected to its backing file in read_sceptre_object() (see attach_response_odm()
## there), not here: by the time get_response_matrix() runs, `so@response_matrix[[1]]` is either
## an already-connected odm or an ordinary in-memory matrix.

ODM_MATERIALIZATION_LIMIT_GB <- 16

#' Project the in-memory size of fully materializing an odm, from a random sample of rows.
#'
#' Sampling rather than reading every row up front, so a dataset that *would* blow the limit is
#' rejected without first paying the cost of the full read it is about to refuse to do.
#'
#' @param limit_gb stop() if the projection exceeds this. Sized far above what any whole-gene,
#'   whole-cell matrix we have measured needs (~1.8 GB); it exists to catch the case
#'   materialization was never meant for, not to police ordinary datasets.
project_and_check_odm_size <- function(m, limit_gb = ODM_MATERIALIZATION_LIMIT_GB, sample_n = 200) {
  n_gene <- nrow(m)
  n_cell <- ncol(m)
  sample_idx <- if (n_gene <= sample_n) seq_len(n_gene) else sort(sample.int(n_gene, sample_n))
  sampled_nnz <- vapply(sample_idx, function(i) sum(m[i, ] != 0), integer(1))

  # 12 bytes/nonzero: a 4-byte column index plus an 8-byte double, the marginal cost of one
  # triplet entry in the sparse matrix this is about to become (row/column pointers are
  # O(n_gene + n_cell) and negligible against O(nnz) for a real single-cell matrix).
  projected_gb <- mean(sampled_nnz) * n_gene * 12 / 1024^3
  if (projected_gb > limit_gb) {
    stop(sprintf(
      paste("Materializing this odm response matrix is projected to need ~%.1f GiB (%d genes x",
            "%d cells, sampled density %.4f), over the %.0f GiB limit. This is the case",
            "out-of-core storage exists for -- the fix here always reads every gene, since the",
            "poscounts size factors are a per-cell reduction over the whole gene set -- so there",
            "is no partial-read option. Raise ODM_MATERIALIZATION_LIMIT_GB deliberately if the",
            "host actually has the memory for it."),
      projected_gb, n_gene, n_cell, mean(sampled_nnz) / n_cell, limit_gb), call. = FALSE)
  }
  invisible(projected_gb)
}

#' Read every row of a (reconnected) odm into one in-memory sparse matrix.
#'
#' One row at a time -- see the "RESPONSE MATRIX ACCESS" note above for why -- assembled via
#' triplets rather than rbind(), which would copy the growing matrix on every row.
materialize_odm_response_matrix <- function(m) {
  n_gene <- nrow(m)
  n_cell <- ncol(m)
  gene_ids <- rownames(m)
  cell_ids <- colnames(m)

  i_parts <- vector("list", n_gene)
  j_parts <- vector("list", n_gene)
  x_parts <- vector("list", n_gene)
  for (g in seq_len(n_gene)) {
    row <- m[g, ]
    nz <- which(row != 0)
    if (length(nz) == 0L) next
    i_parts[[g]] <- rep.int(g, length(nz))
    j_parts[[g]] <- nz
    x_parts[[g]] <- row[nz]
  }

  Matrix::sparseMatrix(
    i = unlist(i_parts, use.names = FALSE),
    j = unlist(j_parts, use.names = FALSE),
    x = as.double(unlist(x_parts, use.names = FALSE)),
    dims = c(n_gene, n_cell),
    dimnames = list(gene_ids, cell_ids),
    repr = "R"  # dgRMatrix, matching what the non-odm path already hands compute_expression_stats()
  )
}

#' Fetch the in-memory response matrix, or fail with a specific message.
#'
#' For an odm-backed object (already reconnected by read_sceptre_object()) this checks the
#' materialized size is sane and returns an ordinary in-memory sparse matrix built from it -- so
#' callers never need to know whether the object started out-of-core.
get_response_matrix <- function(so) {
  m <- so@response_matrix[[1]]

  if (methods::is(m, "odm")) {
    project_and_check_odm_size(m)
    return(materialize_odm_response_matrix(m))
  }
  if (!methods::is(m, "sparseMatrix") && !is.matrix(m)) {
    stop("@response_matrix[[1]] is a ", paste(class(m), collapse = "/"),
         "; expected a sparse Matrix, a base matrix, or an odm.", call. = FALSE)
  }
  m
}
