## sim_input: the container the power simulation reads.
##
## This replaces SingleCellExperiment. The old pipeline stored simulation inputs in an SCE,
## but only ever used it as a container: a per-gene table, a per-cell table, two sparse
## perturbation matrices, row/column subsetting, and dimnames. Nothing on the live code path
## touched rowRanges or any genomic range, so none of the Bioconductor stack was earning its
## place. See docs/development.md.
##
## Structure:
##   genes     character, gene ids (was rownames(sce))
##   cells     character, cell barcodes (was colnames(sce))
##   row_data  data.frame, one row per gene, rownames == genes
##               mean                          size-factor-normalised mean; drives the simulation
##                                             only under the legacy `size_factor` baseline
##               dispersion                    1/theta from sceptre's cached precomputations
##               average_expression_all_cells  raw mean, reported in the output only
##   fitted_coefs  matrix, genes x covariates, rownames == genes
##               sceptre's own Poisson-GLM coefficients for the gene, from
##               @response_precomputations. exp(coefs %*% t(covariate_matrix)) is the expected
##               count the discovery test's null model gives each cell, and is what the
##               simulation draws from. See baseline_expression() in simulate.R.
##   col_data  data.frame, one row per cell, aligned to `cells` by position (no rownames)
##               size_factors                  poscounts size factors
##               <categorical covariates>      factors, e.g. batch_factor, replicate_factor;
##                                             whichever one --cell-batches names is used to
##                                             stratify control sampling
##               pert                          0/1, added by add_pert_status()
##   perts     list of sparse matrices, columns aligned to `cells` by position (no colnames)
##               grna_perts  one row per gRNA
##               cre_perts   one row per perturbation target
##
## Cell barcodes are stored exactly once, in `cells`. Putting them on col_data's rownames and on
## both perturbation matrices' colnames as well meant serialising 586,309 barcodes four times;
## everything indexes by position instead, which is also faster than matching names.
##
## Note the count matrix is deliberately absent: the simulation draws counts from
## row_data$mean and row_data$dispersion, so carrying the real counts through every parallel
## task only cost memory and deserialisation time.

suppressPackageStartupMessages(library(Matrix))

PERT_LEVELS <- c("grna_perts", "cre_perts")

new_sim_input <- function(genes, cells, row_data, col_data, perts, fitted_coefs = NULL) {
  x <- list(
    genes = as.character(genes),
    cells = as.character(cells),
    row_data = row_data,
    col_data = col_data,
    perts = perts,
    fitted_coefs = fitted_coefs
  )
  class(x) <- "sim_input"
  validate_sim_input(x)
  x
}

validate_sim_input <- function(x) {
  if (!inherits(x, "sim_input")) {
    stop("Expected a sim_input object, got ", paste(class(x), collapse = "/"), ".", call. = FALSE)
  }

  if (anyNA(x$genes) || anyDuplicated(x$genes)) {
    stop("sim_input$genes must be unique and non-missing.", call. = FALSE)
  }
  if (anyNA(x$cells) || anyDuplicated(x$cells)) {
    stop("sim_input$cells must be unique and non-missing.", call. = FALSE)
  }

  if (!identical(rownames(x$row_data), x$genes)) {
    stop("sim_input$row_data rownames must equal sim_input$genes, in the same order.", call. = FALSE)
  }
  # col_data is aligned to x$cells positionally and carries no rownames, for the same reason
  # the perturbation matrices carry no colnames -- see the note below.
  if (nrow(x$col_data) != length(x$cells)) {
    stop("sim_input$col_data has ", nrow(x$col_data), " rows but there are ", length(x$cells),
         " cells.", call. = FALSE)
  }

  # Optional only so that a sim_input written before the fitted baseline existed still loads;
  # run_power_simulation.R refuses to use the fitted model without it rather than falling back
  # silently to a different one.
  if (!is.null(x$fitted_coefs)) {
    if (!is.matrix(x$fitted_coefs)) {
      stop("sim_input$fitted_coefs must be a matrix.", call. = FALSE)
    }
    if (!identical(rownames(x$fitted_coefs), x$genes)) {
      stop("sim_input$fitted_coefs rownames must equal sim_input$genes, in the same order.",
           call. = FALSE)
    }
    if (any(!is.finite(x$fitted_coefs))) {
      stop("sim_input$fitted_coefs has non-finite entries.", call. = FALSE)
    }
  }

  required_row <- c("mean", "dispersion", "average_expression_all_cells")
  missing_row <- setdiff(required_row, colnames(x$row_data))
  if (length(missing_row) > 0) {
    stop("sim_input$row_data is missing column(s): ", paste(missing_row, collapse = ", "), ".",
         call. = FALSE)
  }
  if (!"size_factors" %in% colnames(x$col_data)) {
    stop("sim_input$col_data is missing the 'size_factors' column.", call. = FALSE)
  }

  # The simulation cannot draw from a non-finite mean or dispersion, and rnbinom() would
  # silently return NA rather than failing. Catch it here, where the message can be useful.
  for (column in c("mean", "dispersion")) {
    bad <- !is.finite(x$row_data[[column]])
    if (any(bad)) {
      stop(sum(bad), " gene(s) have a non-finite '", column, "', including ",
           paste(utils::head(x$genes[bad], 3), collapse = ", "), ".", call. = FALSE)
    }
  }
  if (any(!is.finite(x$col_data$size_factors)) || any(x$col_data$size_factors <= 0)) {
    stop("sim_input$col_data$size_factors must all be finite and positive.", call. = FALSE)
  }

  missing_levels <- setdiff(PERT_LEVELS, names(x$perts))
  if (length(missing_levels) > 0) {
    stop("sim_input$perts is missing: ", paste(missing_levels, collapse = ", "), ".", call. = FALSE)
  }
  # Columns are aligned to x$cells *positionally* and carry no colnames: storing 586,309 cell
  # barcodes on each perturbation matrix as well as on cells and col_data meant serialising them
  # four times over. Only the dimension is checked here; every accessor below indexes by position.
  for (level in PERT_LEVELS) {
    if (ncol(x$perts[[level]]) != length(x$cells)) {
      stop("sim_input$perts$", level, " has ", ncol(x$perts[[level]]), " columns but there are ",
           length(x$cells), " cells.", call. = FALSE)
    }
    if (is.null(rownames(x$perts[[level]]))) {
      stop("sim_input$perts$", level, " must have rownames (gRNA ids or target names).",
           call. = FALSE)
    }
  }

  invisible(x)
}

#' Restrict to a set of genes, preserving the order given.
#'
#' The perturbation matrices are per-gRNA and per-target, so they are untouched.
subset_genes <- function(x, genes) {
  unknown <- setdiff(genes, x$genes)
  if (length(unknown) > 0) {
    stop(length(unknown), " gene(s) not present in sim_input, including ",
         paste(utils::head(unknown, 3), collapse = ", "), ".", call. = FALSE)
  }
  x$genes <- as.character(genes)
  x$row_data <- x$row_data[x$genes, , drop = FALSE]
  if (!is.null(x$fitted_coefs)) {
    x$fitted_coefs <- x$fitted_coefs[x$genes, , drop = FALSE]
  }
  x
}

#' Restrict to a set of cells, preserving the order given.
#'
#' Takes either cell names or integer positions. Positions are preferred in hot paths: matching
#' names against 586,309 barcodes is not free, and the perturbation matrices are aligned
#' positionally anyway.
subset_cells <- function(x, cells) {
  if (is.numeric(cells)) {
    idx <- as.integer(cells)
    if (any(idx < 1L) || any(idx > length(x$cells))) {
      stop("Cell index out of range.", call. = FALSE)
    }
  } else {
    idx <- match(as.character(cells), x$cells)
    if (anyNA(idx)) {
      stop(sum(is.na(idx)), " cell(s) not present in sim_input, including ",
           paste(utils::head(cells[is.na(idx)], 3), collapse = ", "), ".", call. = FALSE)
    }
  }
  x$cells <- x$cells[idx]
  x$col_data <- x$col_data[idx, , drop = FALSE]
  x$perts <- lapply(x$perts, function(m) m[, idx, drop = FALSE])
  x
}

#' Add the perturbation status of one target as col_data$pert.
#'
#' Stored as integer 0/1. The old code stored a factor and then compared it with
#' `pert_status == 1`, which worked only via factor-to-character coercion; an integer removes
#' that coupling. Values are thresholded at > 0 because the perturbation matrices are
#' indicator matrices, so any nonzero entry means "perturbed".
add_pert_status <- function(x, target, level = PERT_LEVELS) {
  level <- match.arg(level)
  pert_matrix <- x$perts[[level]]
  row <- match(target, rownames(pert_matrix))
  if (is.na(row)) {
    stop("Target '", target, "' is not a row of perts$", level, ".", call. = FALSE)
  }
  # Already in x$cells order (columns are positionally aligned), so no reordering by name.
  x$col_data$pert <- as.integer(pert_matrix[row, ] > 0)
  x
}

n_genes <- function(x) length(x$genes)
n_cells <- function(x) length(x$cells)

print.sim_input <- function(x, ...) {
  cat("sim_input:", n_genes(x), "genes x", n_cells(x), "cells\n")
  cat("  row_data:", paste(colnames(x$row_data), collapse = ", "), "\n")
  cat("  col_data:", paste(colnames(x$col_data), collapse = ", "), "\n")
  cat("  fitted_coefs:",
      if (is.null(x$fitted_coefs)) "absent (legacy sim_input)" else
        paste(ncol(x$fitted_coefs), "covariates"), "\n")
  for (level in names(x$perts)) {
    cat("  perts$", level, ": ", nrow(x$perts[[level]]), " rows\n", sep = "")
  }
  invisible(x)
}
