#!/usr/bin/env Rscript
#
# Turn a sceptre object into the small inputs the power simulation needs.
#
# This is the step that makes the rest of the pipeline cheap. Previously every one of the
# n_splits x n_effect_sizes parallel tasks read *two* objects that each carried the full
# gene-by-cell count matrix -- the sceptre object (1,774 MB in memory on sample1, of which
# 1,095 MB is the response matrix) and perturb_sce.rds -- and the simulation read neither
# matrix. It draws counts from the per-gene mean and dispersion, and sceptre's response matrix
# is overwritten on every rep anyway.
#
# So this script front-loads the work once per sample and emits:
#   sim_input.rds         per-gene and per-cell statistics + the two perturbation matrices
#   sceptre_template.rds  the sceptre object with the response matrix emptied out
#   pairs.tsv             the QC-passing discovery pairs
#   grna_targets.tsv      the gRNA -> target mapping
#
# The last two replace the gene_grna_group_pairs.rds and grna_groups_table.rds inputs, which
# duplicated data already inside the sceptre object.
#
# Usage:
#   prepare_sim_input.R --sceptre-object results/sample1/sceptre_object.rds --outdir prepared/

suppressPackageStartupMessages(library(Matrix))

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    normalizePath(".")
  }
  # lib/ sits one level up from src/, but a workflow engine that stages every input into one flat
  # task directory collapses that, so look beside the script and in it as well. Kept in step with
  # lib_dirs() in cli.R, which resolves the rest of the library the same way.
  for (dir in unique(c(file.path(dirname(here), "lib"), file.path(here, "lib"), here))) {
    if (file.exists(file.path(dir, "cli.R"))) {
      source(file.path(dir, "cli.R"))
      return(invisible(NULL))
    }
  }
  stop("Cannot find lib/cli.R relative to ", here, call. = FALSE)
})
source_lib("sim_input.R", "simulate.R", "sceptre_io.R")

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--sceptre-object", type = "character", default = NULL, dest = "sceptre_object",
              help = "Input sceptre object (.rds). The only required input."),
  make_option("--response-odm", type = "character", default = NULL, dest = "response_odm",
              help = paste("Path to the response matrix's backing .odm file. Only needed if",
                           "--sceptre-object's response matrix is odm-backed (out-of-core);",
                           "ignored otherwise. See lib/sceptre_io.R:attach_response_odm().")),
  make_option("--outdir", type = "character", default = NULL, dest = "outdir",
              help = "Directory for all outputs. Overridden by the individual --out-* options."),
  make_option("--out-sim-input", type = "character", default = NULL, dest = "out_sim_input",
              help = "Output path for sim_input.rds."),
  make_option("--out-sceptre-template", type = "character", default = NULL,
              dest = "out_sceptre_template",
              help = "Output path for sceptre_template.rds."),
  make_option("--out-pairs", type = "character", default = NULL, dest = "out_pairs",
              help = "Output path for the QC-passing pairs TSV."),
  make_option("--out-grna-targets", type = "character", default = NULL, dest = "out_grna_targets",
              help = "Output path for the gRNA-to-target TSV."),
  make_option("--out-threshold", type = "character", default = NULL, dest = "out_threshold",
              help = paste("Output path for the nominal p-value threshold derived from",
                           "@discovery_result. Written so compute_power.R does not have to load",
                           "the whole sceptre object just to read one number.")),
  make_option("--out-analysis-mode", type = "character", default = NULL, dest = "out_analysis_mode",
              help = paste("Output path for the resampling mechanism and MOI read off the sceptre",
                           "object. Same reason as --out-threshold: which test produced a set of",
                           "power numbers should be readable without deserialising the object.")),
  make_option("--all-genes", action = "store_true", default = FALSE, dest = "all_genes",
              help = paste("Keep every gene in the response matrix, not just those in QC-passing",
                           "pairs. Only useful for inspection; the simulation never tests the rest.")),
  make_option("--no-compress", action = "store_false", default = TRUE, dest = "compress",
              help = paste("Write the .rds outputs uncompressed. Compression is on by default",
                           "because it measured both smaller *and* faster to read on sample1:",
                           "59 MB / 1.09s compressed vs 322 MB / 1.20s uncompressed. It costs",
                           "~6s extra at write time, once per sample, against a read paid by",
                           "every simulation task."))
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Derive the power-simulation inputs from a sceptre object."
))
require_options(opts, "sceptre_object")

resolve_output <- function(explicit, default_name) {
  if (!is.null(explicit)) return(explicit)
  if (is.null(opts$outdir)) {
    stop("Provide either --outdir or an explicit path for every --out-* option.", call. = FALSE)
  }
  file.path(opts$outdir, default_name)
}

out_sim_input <- resolve_output(opts$out_sim_input, "sim_input.rds")
out_template <- resolve_output(opts$out_sceptre_template, "sceptre_template.rds")
out_pairs <- resolve_output(opts$out_pairs, "pairs.tsv")
out_grna_targets <- resolve_output(opts$out_grna_targets, "grna_targets.tsv")
out_threshold <- resolve_output(opts$out_threshold, "discovery_threshold.txt")
out_analysis_mode <- resolve_output(opts$out_analysis_mode, "analysis_mode.tsv")

## HELPERS =========================================================================================

#' Build a cell-by-perturbation indicator matrix from one of sceptre's assignment lists.
#'
#' Carried over from create_sce_object.R, which built it in a single sparse allocation rather
#' than by growing a matrix (commit b3cf7ab). Values are forced to 1: these are indicators, and
#' a cell assigned the same gRNA twice must not count double.
#'
#' @param assignments list of integer cell indices, one element per gRNA or target
#' @param cell_index if the indices are positions within `cells_in_use` rather than absolute
#'   cell positions, pass `cells_in_use` to map them
build_assignment_matrix <- function(assignments, n_cols, row_names, col_names, cell_index = NULL) {
  n_rows <- length(assignments)
  assignment_lengths <- lengths(assignments)

  if (sum(assignment_lengths) == 0L) {
    return(sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                        dims = c(n_rows, n_cols), dimnames = list(row_names, col_names)))
  }

  row_idx <- rep.int(seq_len(n_rows), assignment_lengths)
  col_idx <- unlist(assignments, use.names = FALSE)

  if (!is.null(cell_index)) {
    mapped <- rep(NA_integer_, length(col_idx))
    valid <- !is.na(col_idx) & col_idx >= 1L & col_idx <= length(cell_index)
    mapped[valid] <- cell_index[col_idx[valid]]
    col_idx <- mapped
  }

  keep <- !is.na(col_idx) & col_idx >= 1L & col_idx <= n_cols
  if (any(!keep)) {
    message("Dropping ", sum(!keep), " assignment index/indices outside the cell range.")
  }

  out <- sparseMatrix(i = row_idx[keep], j = as.integer(col_idx[keep]),
                      x = rep.int(1, sum(keep)),
                      dims = c(n_rows, n_cols), dimnames = list(row_names, col_names))
  if (length(out@x) > 0L) out@x[] <- 1
  out
}

#' Per-gene and per-cell expression statistics.
#'
#' The only function in the pipeline that reads actual expression values, and therefore the only
#' one an ondisc/odm backend would need to reimplement (see sceptre_io.R).
#'
#' Computes DESeq2-style "poscounts" size factors without ever densifying the matrix. DESeq2
#' itself is not used: DESeqDataSetFromMatrix() coerces to dense, which at 292 x 586,309 would be
#' 1.4 GB and pointless. That is why bioconductor-deseq2 is not a dependency.
compute_expression_stats <- function(counts) {
  n_cell <- ncol(counts)

  # Column-compressed, because every statistic below is a per-cell reduction. The input from
  # sceptre is dgRMatrix (row-compressed), so this conversion is unavoidable.
  started <- proc.time()[["elapsed"]]
  counts <- drop0(as(counts, "CsparseMatrix"))
  if (any(counts@x < 0)) {
    stop("Count matrix contains negative values.", call. = FALSE)
  }
  log_resources("convert to CSC", started)

  # Geometric mean per gene over all cells, computed on the log scale. Genes that are zero
  # everywhere get -Inf and are excluded from the per-cell medians below.
  log_counts_x <- log(counts@x)
  gene_of_nonzero <- counts@i + 1L
  log_geomeans <- as.numeric(
    tapply(log_counts_x, factor(gene_of_nonzero, levels = seq_len(nrow(counts))), sum,
           default = 0)
  ) / n_cell
  row_totals <- Matrix::rowSums(counts)
  log_geomeans[row_totals == 0] <- -Inf
  usable_gene <- is.finite(log_geomeans)

  # Ratio of each observed count to its gene's geometric mean, on the log scale. Fully
  # vectorised over all nonzeros; the original recomputed log() and the validity mask separately
  # inside a per-cell closure, 586,309 times.
  started <- proc.time()[["elapsed"]]
  ratio <- log_counts_x - log_geomeans[gene_of_nonzero]
  ratio[!usable_gene[gene_of_nonzero]] <- NA_real_

  # The nonzeros are already grouped by cell (that is what CSC means), so the per-cell median
  # needs no sort by cell -- just a median within each contiguous slice.
  slice_start <- counts@p[-length(counts@p)] + 1L
  slice_end <- counts@p[-1L]
  size_factors <- vapply(seq_len(n_cell), function(k) {
    if (slice_start[k] > slice_end[k]) return(NA_real_)
    v <- ratio[slice_start[k]:slice_end[k]]
    v <- v[!is.na(v)]
    if (length(v) == 0L) return(NA_real_)
    stats::median.default(v)
  }, numeric(1))
  size_factors <- exp(size_factors)
  log_resources("size factors", started)

  if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
    n_bad <- sum(!is.finite(size_factors) | size_factors <= 0)
    stop(n_bad, " of ", n_cell, " cells have no usable size factor (no nonzero counts in any ",
         "gene with a finite geometric mean). Filter these cells out before running the power ",
         "analysis.", call. = FALSE)
  }

  # Raw mean, reported in the output. Kept distinct from `mean` below, which is normalised and is
  # what the simulation draws from.
  average_expression_all_cells <- Matrix::rowSums(counts) / n_cell

  # Size-factor-normalised mean: divide each column by its size factor, then average.
  normalized <- counts
  normalized@x <- normalized@x / rep.int(size_factors, diff(counts@p))
  normalized_mean <- Matrix::rowSums(normalized) / n_cell

  list(
    size_factors = size_factors,
    average_expression_all_cells = average_expression_all_cells,
    normalized_mean = normalized_mean,
    density = length(counts@x) / (as.numeric(nrow(counts)) * n_cell)
  )
}

## MAIN ============================================================================================

total_started <- proc.time()[["elapsed"]]

log_step("Reading sceptre object: ", opts$sceptre_object)
started <- proc.time()[["elapsed"]]
so <- read_sceptre_object(opts$sceptre_object, response_odm_fp = opts$response_odm)
log_resources("readRDS", started)

# Which test this screen was run under, read off the object rather than assumed. The power
# simulation re-runs this same test, so it decides what the power numbers mean -- and it is set
# upstream and inherited silently. Logged here, and written to analysis_mode.tsv below, because
# three different upstream choices land on permutations without ever saying so: an explicit
# request, low MOI with the default, or an odm-backed response matrix. See
# sceptre_analysis_mode() in lib/sceptre_io.R.
analysis_mode <- sceptre_analysis_mode(so)
log_step("Resampling mechanism: ", analysis_mode$resampling_mechanism,
         " (MOI: ", analysis_mode$moi, ")")
if (analysis_mode$run_permutations) {
  log_step("NOTE: this screen used permutations, not the CRT path every measurement in ",
           "docs/status.md was made under. Power is still correct for THIS screen -- the ",
           "simulation re-runs the test the screen actually ran -- but the cost model does not ",
           "apply and cross-dataset comparisons are comparing across tests.")
}

cells <- sceptre_cells(so)
n_cell <- length(cells)
log_step("Cells: ", n_cell, " (cells_in_use: ", length(so@cells_in_use), ")")

pairs <- qc_passing_pairs(so)
log_step("QC-passing pairs: ", nrow(pairs), " across ",
         length(unique(pairs$grna_target)), " targets")

grna_targets <- grna_target_map(so)
log_step("gRNA-to-target rows: ", nrow(grna_targets))

## Expression statistics -------------------------------------------------------------------------
response_matrix <- get_response_matrix(so)
if (ncol(response_matrix) != n_cell) {
  stop("Response matrix has ", ncol(response_matrix), " columns but @covariate_data_frame has ",
       n_cell, " rows; cell barcodes cannot be aligned.", call. = FALSE)
}

log_step("Computing expression statistics over ", nrow(response_matrix), " x ", n_cell)
stats_out <- compute_expression_stats(response_matrix)
log_step(sprintf("Response matrix density: %.1f%%", 100 * stats_out$density))

## Per-gene table ---------------------------------------------------------------------------------
all_genes <- rownames(response_matrix)
if (is.null(all_genes)) {
  stop("The response matrix has no rownames, so genes cannot be identified.", call. = FALSE)
}

tested_genes <- unique(pairs$response_id)
unknown_genes <- setdiff(tested_genes, all_genes)
if (length(unknown_genes) > 0) {
  stop(length(unknown_genes), " gene(s) appear in the discovery pairs but not in the response ",
       "matrix, including: ", paste(utils::head(unknown_genes, 5), collapse = ", "), ".",
       call. = FALSE)
}

genes <- if (opts$all_genes) all_genes else intersect(all_genes, tested_genes)
log_step("Genes kept: ", length(genes), " of ", length(all_genes),
         if (opts$all_genes) " (--all-genes)" else " (those in QC-passing pairs)")

# Fails loudly if any kept gene has no cached precomputation. The old code produced a list column
# with NULL holes here, which unlist() silently dropped, shifting every later gene's dispersion.
dispersion <- build_dispersion_vector(so@response_precomputations, genes)

gene_idx <- match(genes, all_genes)
row_data <- data.frame(
  mean = stats_out$normalized_mean[gene_idx],
  dispersion = as.numeric(dispersion),
  average_expression_all_cells = stats_out$average_expression_all_cells[gene_idx],
  row.names = genes
)

## Per-cell table ---------------------------------------------------------------------------------
# The batch column is *not* assumed to be called "batch": on sample1 the covariates are
# batch_factor and replicate_factor, so the old `"batch" %in% colnames(...)` test never fired and
# the batch-stratified control sampling was unreachable. Every non-numeric covariate is carried
# through instead, and --cell-batches names the one to use at simulation time.
# No rownames: with 586,309 cells, storing the barcodes here as well as in $cells and on both
# perturbation matrices serialised them four times over.
col_data <- data.frame(size_factors = stats_out$size_factors)
covariates <- so@covariate_data_frame
candidate_batches <- names(which(vapply(
  covariates, function(column) is.factor(column) || is.character(column), logical(1)
)))
for (column in candidate_batches) {
  # Factor, not character: the levels are stored once and the per-cell values become a 4-byte
  # integer each instead of a pointer to a string.
  col_data[[column]] <- as.factor(covariates[[column]])
}
if (length(candidate_batches) > 0) {
  log_step("Categorical covariates available for --cell-batches: ",
           paste(candidate_batches, collapse = ", "))
} else {
  log_step("No categorical covariates found; --cell-batches will not be usable.")
}

## Perturbation matrices --------------------------------------------------------------------------
log_step("Building perturbation matrices")
started <- proc.time()[["elapsed"]]

individual_grna_assignments <- so@initial_grna_assignment_list
max_grna_idx <- suppressWarnings(max(unlist(individual_grna_assignments, use.names = FALSE)))
if (is.finite(max_grna_idx) && max_grna_idx > n_cell) {
  stop("@initial_grna_assignment_list contains cell index ", max_grna_idx, " but there are only ",
       n_cell, " cells; these indices were assumed to be absolute cell positions.", call. = FALSE)
}
grna_perts <- build_assignment_matrix(
  assignments = individual_grna_assignments,
  n_cols = n_cell,
  row_names = names(individual_grna_assignments),
  col_names = NULL
)

# grna_group_idxs indexes positions *within cells_in_use*, not absolute cell positions, so it
# needs the cells_in_use mapping.
target_assignments <- so@grna_assignments$grna_group_idxs
cre_perts <- build_assignment_matrix(
  assignments = target_assignments,
  n_cols = n_cell,
  row_names = names(target_assignments),
  col_names = NULL,
  cell_index = so@cells_in_use
)
log_resources("perturbation matrices", started)

missing_targets <- setdiff(unique(pairs$grna_target), rownames(cre_perts))
if (length(missing_targets) > 0) {
  stop(length(missing_targets), " target(s) in the QC-passing pairs have no entry in ",
       "@grna_assignments$grna_group_idxs, including: ",
       paste(utils::head(missing_targets, 5), collapse = ", "), ".", call. = FALSE)
}

## Assemble and write ------------------------------------------------------------------------------
sim <- new_sim_input(
  genes = genes,
  cells = cells,
  row_data = row_data,
  col_data = col_data,
  perts = list(grna_perts = grna_perts, cre_perts = cre_perts)
)
print(sim)

log_step("Writing outputs")
started <- proc.time()[["elapsed"]]
for (path in unique(dirname(c(out_sim_input, out_template, out_pairs, out_grna_targets,
                              out_threshold, out_analysis_mode)))) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

saveRDS(sim, out_sim_input, compress = opts$compress)
saveRDS(slim_sceptre_object(so), out_template, compress = opts$compress)
write_tsv_file(pairs, out_pairs)
write_tsv_file(grna_targets, out_grna_targets)

# Unconditional, unlike the threshold below: an object with no discovery result still has a
# resampling mechanism, and a missing file would be indistinguishable from an old run that
# predates this output.
writeLines(format_analysis_mode(analysis_mode), out_analysis_mode)

# The significance threshold that a simulated replicate has to beat, taken from the real discovery
# results already inside the object. Written as a plain number so the downstream step does not have
# to deserialise the sceptre object for it. If the object has no discovery result, the threshold has
# to be supplied explicitly downstream via --alpha.
threshold <- tryCatch(discovery_threshold(so), error = function(e) {
  log_step("No usable discovery threshold in the object (", conditionMessage(e),
           ") -- compute_power.R will need --alpha.")
  NA_real_
})
if (!is.na(threshold)) {
  log_step(sprintf("Discovery p-value threshold: %.6g (largest significant nominal p-value)",
                   threshold))
  writeLines(format(threshold, digits = 17), out_threshold)
}
log_resources("write outputs", started)

report_size <- function(path) {
  sprintf("  %-24s %7.1f MB", basename(path), file.size(path) / 1024^2)
}
message("Output sizes:")
message(report_size(out_sim_input))
message(report_size(out_template))
message(report_size(out_pairs))
message(report_size(out_grna_targets))
if (file.exists(out_threshold)) message(report_size(out_threshold))
message(sprintf("  %-24s %7.1f MB  (input, for comparison)",
                basename(opts$sceptre_object), file.size(opts$sceptre_object) / 1024^2))

log_resources("TOTAL", total_started)
