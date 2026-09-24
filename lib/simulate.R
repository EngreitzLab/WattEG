## Simulate Perturb-seq counts with a specified effect size, including guide-to-guide
## variability in that effect size.
##
## Ported from R/power_simulations_fun.R. fit_negbinom_deseq2() is not carried over: nothing
## called it, and it was the only reason the pipeline depended on DESeq2.
##
## It departs from the original deliberately, and docs/methods.md ("What simulated power means")
## is the specification. In short:
##   - each cell keeps its own size factor and covariates: the original shuffled size factors
##     across cells, pairing one cell's perturbation status with another's library size;
##   - the expected counts are sceptre's own fitted model, exp(X . beta), not mean * size factor;
##   - the realised mean effect over the perturbed cells is pinned to the requested one in every
##     replicate (power at a fixed element effect). The original meant to, but centred before
##     reordering and so shifted the wrong columns;
##   - guides that belong to anything else have an effect of exactly 1; the original gave them
##     N(1, guide_sd), which double-counted noise already in the fitted dispersion.
## Each of these changes the RNG stream, so no seed reproduces the original draw for draw.

suppressPackageStartupMessages(library(Matrix))

#' A gene's expected counts per cell, before any perturbation.
#'
#' This decides what "unperturbed" means, and therefore what power is measured against. There
#' are two ways to produce it and they are not equivalent.
#'
#' `fitted` (the default) is `exp(X %*% coefs)`: the expected count sceptre's own null model
#' gives that cell for that gene, so the simulation and the test that judges it are on one scale
#' by construction.
#'
#' `size_factor` is `mean_i * sf_j`, which is what this pipeline did until 2026-09-21. It is kept
#' so the sweeps already run can be reproduced, and for no other reason. Three things are wrong
#' with it, in increasing order of weight:
#'
#'   * It mixes two models -- the dispersion comes from sceptre's negative binomial and the level
#'     from a DESeq2-style normalisation, for the same gene.
#'   * It gets the level wrong. `row_data$mean` sits 16 % below the mean sceptre's model implies;
#'     multiplying by the size factor recovers most of that and leaves the simulated genes about
#'     4 % low on day0. The residual is a dropped covariance term: `mean_i` is a mean of ratios,
#'     and `E[x*sf] = E[x]E[sf] + Cov(x, sf)`.
#'   * It gets the SHAPE wrong, which the level hides. sceptre's mean varies with every covariate;
#'     `mean_i * sf_j` varies with one scalar per cell. Measured against the real counts on day0
#'     over 60 genes and 567,690 cells, the fitted baseline predicts the zero fraction to 0.0007
#'     against 0.0053, and 99.5 % of the observed count variance against 86.5 %.
#'
#' Note the variance error pulls the OPPOSITE way from the level error -- less variance inflates
#' power where less expression deflates it -- so which way the change moves power is not obvious
#' and has not been measured.
#'
#' @param x sim_input, already subset to the genes and cells being simulated
#' @param covariate_matrix cells x covariates, rows aligned to x$cells. Required by `fitted`.
#' @return a dense genes x cells matrix of expected counts
baseline_expression <- function(x, covariate_matrix = NULL,
                                model = c("fitted", "size_factor")) {
  model <- match.arg(model)
  if (model == "size_factor") {
    # outer() replaces the original matrix(rep(...)) + sweep() pair: same values, one allocation
    # instead of three.
    return(outer(x$row_data$mean, x$col_data$size_factors))
  }

  if (is.null(x$fitted_coefs)) {
    stop("This sim_input carries no fitted_coefs, so the 'fitted' baseline cannot be formed. ",
         "Re-run prepare_sim_input.R, or pass --expression-model size_factor to reproduce the ",
         "pre-2026-09-21 behaviour.", call. = FALSE)
  }
  if (is.null(covariate_matrix)) {
    stop("The 'fitted' baseline needs the covariate matrix.", call. = FALSE)
  }
  if (nrow(covariate_matrix) != length(x$cells)) {
    stop("covariate_matrix has ", nrow(covariate_matrix), " rows but there are ",
         length(x$cells), " cells; they must be aligned.", call. = FALSE)
  }
  if (ncol(covariate_matrix) != ncol(x$fitted_coefs)) {
    stop("covariate_matrix has ", ncol(covariate_matrix), " columns but fitted_coefs has ",
         ncol(x$fitted_coefs), "; they come from different designs.", call. = FALSE)
  }
  # The product below pairs coefficients with covariates by POSITION, so a matching count is not
  # enough: the same covariates in a different order would give every cell a plausible, wrong
  # expected count. Both carry the design's column names (e.g. "log(response_n_umis)").
  if (!identical(colnames(covariate_matrix), colnames(x$fitted_coefs))) {
    stop("covariate_matrix and fitted_coefs name their columns differently (",
         paste(utils::head(colnames(covariate_matrix), 3), collapse = ", "), " vs ",
         paste(utils::head(colnames(x$fitted_coefs), 3), collapse = ", "),
         "); the product would pair the wrong coefficient with each covariate.", call. = FALSE)
  }
  # tcrossprod(A, B) is A %*% t(B): (genes x p) against (cells x p) gives genes x cells.
  exp(tcrossprod(x$fitted_coefs, covariate_matrix))
}

#' Simulate counts for one perturbation.
#'
#' @param x sim_input, already subset to the genes being tested
#' @param effect_size_mat genes x cells multiplier, columns in the same order as x$cells
#' @param baseline genes x cells expected counts from baseline_expression(). Computed once per
#'   target rather than once per simulation: it does not depend on the effect size or the draw.
#' @return a dense genes x cells matrix of counts
draw_counts <- function(x, effect_size_mat, baseline) {
  gene_dispersions <- x$row_data$dispersion

  n_gene <- length(gene_dispersions)
  n_cell <- length(x$cells)

  if (!identical(dim(effect_size_mat), c(n_gene, n_cell))) {
    stop("effect_size_mat is ", paste(dim(effect_size_mat), collapse = " x "),
         " but ", n_gene, " x ", n_cell, " was expected.", call. = FALSE)
  }
  if (!identical(dim(baseline), c(n_gene, n_cell))) {
    stop("baseline is ", paste(dim(baseline), collapse = " x "),
         " but ", n_gene, " x ", n_cell, " was expected.", call. = FALSE)
  }

  # mu[i, j] = baseline[i, j] * effect_size[i, j]. Each cell keeps its own baseline: the
  # effect-size matrix is indexed by cell, so shuffling would pair one cell's perturbation status
  # with another cell's expected expression.
  mu <- baseline * effect_size_mat

  # size = theta = 1 / dispersion. mu is consumed column-major, and `size` recycles over the
  # gene index, which is why gene_dispersions must be exactly n_gene long -- see
  # build_dispersion_vector() for the guard that guarantees it.
  counts <- rnbinom(n_cell * n_gene, mu = mu, size = 1 / gene_dispersions)

  matrix(counts, nrow = n_gene, ncol = n_cell, dimnames = list(x$genes, x$cells))
}

#' Convert simulated counts into the matrix class sceptre requires.
#'
#' sceptre needs CSR (`dgRMatrix`) here and there is no way around it: this pipeline assigns
#' `@response_matrix` directly, bypassing the normalisation `import_data()` would have applied,
#' and the consumers are unforgiving. `load_row()` (sceptre's s4_helpers.R) dispatches on `odm`
#' and `dgRMatrix` with *no else branch*, so handing it a dense matrix returns NULL silently
#' rather than raising an error; and with `run_permutations = FALSE` the CRT path reads
#' `response_matrix@j/@p/@x` directly, which only a dgRMatrix has.
#'
#' So density does not decide the representation here -- but it does decide the cost, and it is
#' worth logging. The original built a dgCMatrix via `Matrix(..., sparse = TRUE)` and then
#' re-built it as dgRMatrix with `as(..., "RsparseMatrix")`; this constructs the dgRMatrix in a
#' single pass instead. Rows are iterated because the gene count is small (single or double
#' digits) and iterating rows yields CSR's required within-row column ordering for free.
#'
#' @param counts dense genes x cells matrix from draw_counts()
#' @param report_density log the realised nonzero fraction
as_sceptre_response_matrix <- function(counts, report_density = TRUE) {
  n_gene <- nrow(counts)
  n_cell <- ncol(counts)

  nonzero_cols <- lapply(seq_len(n_gene), function(g) which(counts[g, ] != 0))
  per_row <- lengths(nonzero_cols)

  if (isTRUE(report_density)) {
    density <- sum(per_row) / (as.numeric(n_gene) * n_cell)
    message(sprintf(
      "Simulated matrix %d x %d, density %.1f%% (%s storage would be smaller; sceptre requires CSR)",
      n_gene, n_cell, 100 * density, if (density > 0.5) "dense" else "sparse"
    ))
  }

  values <- unlist(
    lapply(seq_len(n_gene), function(g) counts[g, nonzero_cols[[g]]]),
    use.names = FALSE
  )

  new("dgRMatrix",
      j = as.integer(unlist(nonzero_cols, use.names = FALSE)) - 1L,
      p = c(0L, as.integer(cumsum(per_row))),
      x = as.numeric(values),
      Dim = c(n_gene, n_cell),
      Dimnames = list(rownames(counts), colnames(counts)))
}

#' Build the per-gene dispersion vector from sceptre's cached precomputations.
#'
#' The original wrote `rowData(sce)$dispersion <- dispersion_values[rownames(sce)]`, which
#' produced a *list* column with a NULL entry for any gene sceptre had no precomputation for.
#' `unlist()` then dropped those NULLs, shortening the vector, and rnbinom() recycled it -- so
#' every gene after the first gap would have been simulated with another gene's dispersion, with
#' no warning. On sample1 this never fired (all 239 tested genes have precomputations, though 19
#' of the 292 genes in the response matrix do not), but the failure is silent, so it is worth a
#' hard guard rather than a comment.
#'
#' @param precomputations sceptre_object@response_precomputations
#' @param genes genes to build the vector for
#' @return named numeric vector, one finite value per gene, in the order given
build_dispersion_vector <- function(precomputations, genes) {
  missing_genes <- setdiff(genes, names(precomputations))
  if (length(missing_genes) > 0) {
    stop(length(missing_genes), " gene(s) have no entry in @response_precomputations and so no ",
         "dispersion estimate, including: ",
         paste(utils::head(missing_genes, 5), collapse = ", "),
         ". Re-run sceptre's precomputation, or drop these genes from the discovery pairs.",
         call. = FALSE)
  }

  theta <- vapply(precomputations[genes], function(p) p$theta, numeric(1))
  dispersion <- 1 / theta

  if (any(!is.finite(dispersion))) {
    bad <- genes[!is.finite(dispersion)]
    stop(length(bad), " gene(s) have a non-finite dispersion (theta of 0 or NA), including: ",
         paste(utils::head(bad, 5), collapse = ", "), ".", call. = FALSE)
  }

  stats::setNames(dispersion, genes)
}

#' Build the per-gene coefficient matrix from sceptre's cached precomputations.
#'
#' The sibling of build_dispersion_vector(), reading the other half of the same fit. Taking both
#' from one fit is the point: until 2026-09-21 the dispersion came from sceptre's negative
#' binomial and the expression level from a DESeq2-style normalisation, so a simulated gene's
#' noise and its level came from different models of the same data.
#'
#' @param precomputations sceptre_object@response_precomputations
#' @param genes genes to build the matrix for
#' @return genes x covariates matrix, rownames == genes
build_fitted_coefs_matrix <- function(precomputations, genes) {
  missing_genes <- setdiff(genes, names(precomputations))
  if (length(missing_genes) > 0) {
    stop(length(missing_genes), " gene(s) have no entry in @response_precomputations and so no ",
         "fitted coefficients, including: ",
         paste(utils::head(missing_genes, 5), collapse = ", "),
         ". Re-run sceptre's precomputation, or drop these genes from the discovery pairs.",
         call. = FALSE)
  }

  coefs <- do.call(rbind, lapply(precomputations[genes], function(p) p$fitted_coefs))
  rownames(coefs) <- genes

  if (any(!is.finite(coefs))) {
    bad <- genes[apply(!is.finite(coefs), 1, any)]
    stop(length(bad), " gene(s) have a non-finite fitted coefficient, including: ",
         paste(utils::head(bad, 5), collapse = ", "), ".", call. = FALSE)
  }
  coefs
}

## GUIDE-LEVEL VARIABILITY =========================================================================

#' Pick one expressed guide at random per cell, from a CSC perturbation matrix.
#'
#' Returns 0-based row indices, or -1 for a cell with no guide (the caller adds 1, mapping it to
#' the no-effect row).
sample_guide <- function(pert_status) {
  num_cols <- dim(pert_status)[[2]]
  return_vector <- integer(num_cols)

  for (col_idx in seq_len(num_cols)) {
    start_idx <- pert_status@p[col_idx] + 1
    end_idx <- pert_status@p[col_idx + 1]
    column_non_zeros <- end_idx - (start_idx - 1)

    if (column_non_zeros > 0) {
      selected_i <- pert_status@i[start_idx + sample(0:(column_non_zeros - 1), 1)]
      return_vector[col_idx] <- selected_i
    } else {
      return_vector[col_idx] <- -1L
    }
  }

  return_vector
}

#' Collapse a perturbation matrix to one status per cell, choosing randomly among multiples.
convert_pert_mat_to_vector <- function(pert_mat) {
  sample_guide(pert_mat) + 1
}

#' The cell ordering create_guide_pert_status() produces: perturbed cells then control cells.
#'
#' `es_mat[, order(cell_order(pert_status))]` restores x$cells order. Computed once per target
#' rather than once per rep -- the original re-derived it by name inside the rep loop.
cell_order <- function(pert_status) {
  c(which(pert_status == 1), which(pert_status == 0))
}

#' Per-cell gRNA perturbation status: which guide (if any) each cell carries.
#'
#' 0 is no guide, 1..length(pert_guides) is one of this target's guides (by its position in
#' `pert_guides`), and anything above that is a guide belonging to something else. Control-cell
#' statuses are offset by `length(pert_guides)`, the number of target guides -- not by the highest
#' index any perturbed cell happened to carry. Offsetting by that maximum, as this did until
#' 2026-09-24, put a control guide on a *targeting* row of the effect table whenever the
#' last-listed target guide was carried by no perturbed cell.
#'
#' The subsets keep their dimensions (`drop = FALSE`) and are always converted. A subset that
#' collapsed to a vector -- one perturbed cell carrying several target guides, say -- used to skip
#' the conversion and return a status vector of the wrong length, which the positional reorder
#' downstream then misread without an error.
#'
#' IMPORTANT: the result is ordered perturbed-cells-first-then-control-cells (each block in
#' ascending cell position), which is *not* the order of x$cells. Do not reorder it by hand: go
#' through simulate_effect_sizes(), which applies the permutation in the one order that is
#' correct. The original relied on cell barcodes as names for this; positions are cheaper and
#' cannot silently mismatch.
create_guide_pert_status <- function(pert_status, grna_perts, pert_guides) {
  grnas_pert_cells <- convert_pert_mat_to_vector(
    grna_perts[pert_guides, pert_status == 1, drop = FALSE]
  )
  grnas_ctrl_cells <- convert_pert_mat_to_vector(
    grna_perts[!rownames(grna_perts) %in% pert_guides, pert_status == 0, drop = FALSE]
  )

  ctrl_perts <- grnas_ctrl_cells > 0
  grnas_ctrl_cells[ctrl_perts] <- grnas_ctrl_cells[ctrl_perts] + length(pert_guides)

  status <- c(grnas_pert_cells, grnas_ctrl_cells)
  if (length(status) != length(pert_status)) {
    stop("create_guide_pert_status() built ", length(status), " statuses for ",
         length(pert_status), " cells.", call. = FALSE)
  }
  status
}

#' Effect-size matrix with guide-to-guide variability among the target's guides.
#'
#' Each of the target's guides draws its own effect size around the requested relative expression,
#' with standard deviation guide_sd, independently for every gene; negative draws clamp to 0.
#'
#' Every other guide has an effect of exactly 1. Those guides belong to other elements, or are
#' non-targeting, and have no business moving this gene. The dispersion the counts are drawn with
#' was fitted to real cells that already carry their real guides, so it already contains whatever
#' those guides do. Drawing an extra N(1, guide_sd) multiplier for them -- as every version did
#' until 2026-09-24, going back to the original DC-TAP code -- counted that noise twice: a gene
#' with theta 146 came back from its own simulated null data with theta 42, and power was
#' understated for highly expressed, low-dispersion genes.
create_effect_size_matrix <- function(grna_pert_status, pert_guides, gene_effect_sizes, guide_sd) {
  n_pert_guides <- length(pert_guides)
  n_ctrl_guides <- max(0L, max(grna_pert_status) - n_pert_guides)
  n_genes <- length(gene_effect_sizes)

  # matrix() because vapply() returns a plain vector when there is a single guide.
  guide_effect_sizes_pert <- matrix(
    vapply(gene_effect_sizes, FUN = rnorm, n = n_pert_guides, sd = guide_sd,
           FUN.VALUE = numeric(n_pert_guides)),
    nrow = n_pert_guides, ncol = n_genes
  )
  guide_effect_sizes_pert[guide_effect_sizes_pert < 0] <- 0
  guide_effect_sizes_ctrl <- matrix(1, nrow = n_ctrl_guides, ncol = n_genes)

  # Row 1 is the no-effect row, used by cells carrying no guide.
  guide_effect_sizes <- rbind(1, guide_effect_sizes_pert, guide_effect_sizes_ctrl)

  t(guide_effect_sizes[grna_pert_status + 1, , drop = FALSE])
}

#' Pin each gene's realised mean effect over the perturbed cells to the requested effect size.
#'
#' This is what makes the simulated power the power at a FIXED element effect (estimand A, decided
#' 2026-09-24; see docs/methods.md): in every replicate, the cell-weighted mean effect across the
#' perturbed cells equals the requested one. The guides still differ from each other within a
#' replicate; only the replicate-to-replicate wobble of their mean is removed. It is not a
#' correction for clamping, which is what this docstring used to claim: at es 0.15 a guide clamps
#' with probability ~3e-11.
#'
#' For each gene the result is pmax(v + c, 0), where v are the perturbed cells' (already clamped)
#' effects and c is the one constant that puts the mean exactly on the target. Where nothing would
#' clamp -- every realistic case up to es 0.5 -- that is the plain shift c = target - mean(v), to a
#' few ulps. At
#' strong knockdowns a shift down pushes some cells below zero, so c is solved exactly rather than
#' approached: see pin_to_mean(). An earlier version shifted, clamped and repeated; that converges
#' only linearly once most cells clamp, and at es >= ~0.99 ran out of iterations and stopped a task
#' on a pin that exists.
#'
#' Control cells must already be exactly 1 (create_effect_size_matrix() guarantees it), so there is
#' nothing to centre for them; this checks it rather than assuming it.
center_effect_size_matrix <- function(effect_size_mat, pert_status, gene_effect_sizes,
                                      tol = 1e-12) {
  # The check is for real failures (NaN, a logic error), not rounding, so its tolerance grows with
  # the number of perturbed cells the way summation error does. The mean is taken with mean(), which
  # is compensated, not rowMeans(), which on platforms without long double accumulates n * eps.
  if (length(pert_status) != ncol(effect_size_mat)) {
    stop("pert_status has ", length(pert_status), " entries but the effect-size matrix has ",
         ncol(effect_size_mat), " cells.", call. = FALSE)
  }
  is_pert <- pert_status == 1

  if (any(effect_size_mat[, !is_pert] != 1)) {
    stop("Control cells must carry an effect of exactly 1; found other values. ",
         "See create_effect_size_matrix().", call. = FALSE)
  }

  if (any(is_pert)) {
    pert <- effect_size_mat[, is_pert, drop = FALSE]
    for (g in seq_len(nrow(pert))) {
      pert[g, ] <- pin_to_mean(pert[g, ], gene_effect_sizes[[g]])
    }
    gap <- gene_effect_sizes - apply(pert, 1, mean)
    allowed <- max(tol, 4 * ncol(pert) * .Machine$double.eps)
    if (any(!is.finite(gap)) || any(abs(gap) >= allowed)) {
      stop("Could not pin the realised mean effect to the requested one (largest miss ",
           signif(max(abs(gap)), 3), ").", call. = FALSE)
    }
    effect_size_mat[, is_pert] <- pert
  }
  effect_size_mat
}

#' The values pmax(v + c, 0) whose mean is exactly `target`, for the one constant c that does it.
#'
#' f(c) = mean(pmax(v + c, 0)) is continuous, non-decreasing and piecewise linear, with a kink where
#' each value hits zero; for any target > 0 it has a root. Sorting v in decreasing order, if exactly
#' the top j values stay positive then f(c) = (sum of those j + j * c) / n, so c = (n * target -
#' sum) / j. The right j is the largest one whose kink, f(-s_j), is still at or below the target.
#'
#' The closed form carries the rounding of a long cumulative sum, and perturbed-cell effects are
#' heavily tied (each cell takes one of a few guide values), so that error adds up rather than
#' cancelling: at tens of thousands of cells it reached ~1e-12. One correction step on the linear
#' piece the solution lies on, using the compensated mean(), removes it.
pin_to_mean <- function(v, target) {
  n <- length(v)
  s <- sort(v, decreasing = TRUE)
  cs <- cumsum(s)
  f_at_kink <- (cs - seq_len(n) * s) / n    # f(-s_j): the mean when the j-th largest just hits 0
  j <- max(which(f_at_kink <= target))      # f(-s_1) = 0, so j >= 1 whenever target >= 0
  out <- pmax(v + (n * target - cs[j]) / j, 0)
  pos <- which(out > 0)   # which(), not a logical mask: a NaN input must reach the caller's check
  if (length(pos) > 0) {
    out[pos] <- pmax(out[pos] + (target - mean(out)) * n / length(pos), 0)
  }
  out
}

#' One replicate's effect-size matrix, in cell order, pinned to the requested effect.
#'
#' The only path from a guide status to the matrix draw_counts() consumes. It exists so that the
#' simulation and its tests run the SAME sequence: the order of these steps is exactly what was
#' wrong from the original DC-TAP code until 2026-09-21 (centring before the reorder shifted the
#' wrong columns), and a test that re-implements the order by hand cannot catch a regression in it.
#'
#' @param grna_pert_status from create_guide_pert_status(), perturbed-then-control order
#' @param pert_status 1/0 per cell, in cell order
#' @param restore_cell_order order(cell_order(pert_status)), computed once per target
#' @return genes x cells, in cell order
simulate_effect_sizes <- function(grna_pert_status, pert_status, restore_cell_order,
                                  pert_guides, gene_effect_sizes, guide_sd) {
  es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                      gene_effect_sizes = gene_effect_sizes, guide_sd = guide_sd)
  es_mat <- es_mat[, restore_cell_order, drop = FALSE]
  center_effect_size_matrix(es_mat, pert_status = pert_status,
                            gene_effect_sizes = gene_effect_sizes)
}
