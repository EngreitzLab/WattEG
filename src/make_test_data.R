#!/usr/bin/env Rscript
#
# Generate a small synthetic sceptre object the whole pipeline can run on.
#
# WHY THIS EXISTS
#
# Every other test in this repository is a pure-function test. Nothing exercises the sceptre path --
# the precomputations, the discovery threshold, the QC-passing pair table -- without a 448 MB object
# of unpublished lab data, which cannot go in CI and cannot be shared. This produces an object with
# the same shape and the same slots, small enough to build in under a minute.
#
# WHAT IT HAS TO SATISFY
#
# The pipeline's entry point, read_sceptre_object() in lib/sceptre_io.R, is strict, and rightly so.
# The object must have assign_grnas() and run_qc() recorded in @functs_called, a
# grna_integration_strategy of "union", pairs carrying pass_qc, and -- because the significance
# threshold is derived from the real discovery results rather than from a nominal alpha --
# @discovery_result with at least one significant pair. That last one is the reason this script
# simulates actual knockdowns rather than pure noise: an object with no discoveries is not merely
# uninteresting, it makes discovery_threshold() error, by design.
#
# Usage:
#   make_test_data.R --out tests/data/sceptre_object.rds
#   make_test_data.R --out obj.rds --n-genes 40 --n-cells 4000 --seed 1

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

source_lib("sceptre_io.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(sceptre)
})

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output path for the synthetic sceptre object (.rds)."),
  make_option("--n-genes", type = "integer", default = 30L, dest = "n_genes",
              help = "Number of genes [default %default]."),
  make_option("--n-cells", type = "integer", default = 3000L, dest = "n_cells",
              help = "Number of cells [default %default]."),
  make_option("--n-targets", type = "integer", default = 12L, dest = "n_targets",
              help = "Number of perturbation targets [default %default]."),
  make_option("--grnas-per-target", type = "integer", default = 3L, dest = "grnas_per_target",
              help = "gRNAs per target [default %default]."),
  make_option("--n-nt-grnas", type = "integer", default = 6L, dest = "n_nt_grnas",
              help = "Non-targeting control gRNAs [default %default]."),
  make_option("--true-knockdown", type = "double", default = 0.6, dest = "true_knockdown",
              help = paste("Fractional decrease applied to the genes each target really",
                           "regulates. Large on purpose: the object is useless unless some pair",
                           "comes out significant [default %default].")),
  make_option("--genes-per-target", type = "integer", default = 2L, dest = "genes_per_target",
              help = "Genes each target really regulates [default %default]."),
  make_option("--seed", type = "integer", default = 20250812L, dest = "seed",
              help = "RNG seed [default %default]."),
  # The two axes the pipeline inherits from the object rather than setting, so they are the two it
  # can only be shown to handle by building an object that carries them. See
  # sceptre_analysis_mode() in lib/sceptre_io.R for why they are worth varying.
  make_option("--moi", type = "character", default = "high", dest = "moi",
              help = paste("'high' or 'low' [default %default]. High MOI makes the control group",
                           "the complement; low MOI makes it the NT cells, a path no real dataset",
                           "has yet taken through this pipeline.")),
  make_option("--resampling-mechanism", type = "character", default = "default",
              dest = "resampling_mechanism",
              help = paste("'crt', 'permutations' or 'default' [default %default]. 'default'",
                           "resolves to crt for high MOI and permutations for low MOI, which is",
                           "the implicit coupling worth being able to override here."))
)
opts <- parse_args(OptionParser(option_list = option_list))
require_options(opts, "out")

if (!opts$moi %in% c("high", "low")) {
  stop("--moi must be 'high' or 'low', not '", opts$moi, "'.", call. = FALSE)
}
if (!opts$resampling_mechanism %in% c("crt", "permutations", "default")) {
  stop("--resampling-mechanism must be 'crt', 'permutations' or 'default', not '",
       opts$resampling_mechanism, "'.", call. = FALSE)
}

init_seed(opts$seed)

## SYNTHETIC COUNTS ================================================================================

n_genes <- opts$n_genes
n_cells <- opts$n_cells
n_targets <- opts$n_targets

gene_ids <- sprintf("GENE%03d", seq_len(n_genes))
cell_ids <- sprintf("CELL%05d", seq_len(n_cells))
targets <- sprintf("chr1:%d-%d", seq_len(n_targets) * 10000, seq_len(n_targets) * 10000 + 300)

# Built target-major and paired positionally with the target vector, so the two cannot drift.
grna_target_data_frame <- rbind(
  expand.grid(i = seq_len(opts$grnas_per_target), grna_target = targets,
              stringsAsFactors = FALSE) |>
    within({
      grna_id <- sprintf("%s_g%d", grna_target, i)
      rm(i)
    }),
  data.frame(grna_target = "non-targeting",
             grna_id = sprintf("NT_g%d", seq_len(opts$n_nt_grnas)),
             stringsAsFactors = FALSE)
)[, c("grna_id", "grna_target")]
rownames(grna_target_data_frame) <- NULL
grna_ids <- grna_target_data_frame$grna_id
stopifnot(!anyDuplicated(grna_ids))

log_step("Synthesising ", n_genes, " genes x ", n_cells, " cells, ",
         n_targets, " targets, ", length(grna_ids), " gRNAs")

# A VARIABLE number of gRNAs per cell, which is not decoration.
#
# The obvious simplification -- exactly one gRNA per cell -- produces an object sceptre refuses:
# grna_n_nonzero is then constant at 1, so it is collinear with the intercept and
# convert_covariate_df_to_design_matrix() rejects the formula as containing "redundant
# information". Real high-MOI data has a spread, and the covariate only carries information
# because of it.
n_grnas_per_cell <- 1L + rpois(n_cells, lambda = 1.5)
cell_grnas <- lapply(n_grnas_per_cell, function(k) {
  sample(seq_along(grna_ids), min(k, length(grna_ids)))
})

grna_matrix <- sparseMatrix(
  i = unlist(cell_grnas),
  j = rep(seq_len(n_cells), lengths(cell_grnas)),
  x = rpois(sum(lengths(cell_grnas)), lambda = 40) + 5,
  dims = c(length(grna_ids), n_cells),
  dimnames = list(grna_ids, cell_ids)
)

# Cell x target membership: a cell is perturbed for a target if it carries any of that target's
# gRNAs. With several gRNAs per cell a cell can belong to more than one target, exactly as in the
# real data -- which is also why the control group is the complement rather than the NT cells.
cell_targets <- lapply(cell_grnas, function(idx) {
  unique(grna_target_data_frame$grna_target[idx])
})
is_perturbed_for <- function(target) vapply(cell_targets, function(t) target %in% t, logical(1))

# Gene expression: a spread of means so some genes are well powered and some are not, which is the
# situation the analysis exists to distinguish.
gene_mean <- exp(runif(n_genes, log(5), log(300)))
gene_size <- runif(n_genes, 2, 20)   # NB size parameter; 1/size is the dispersion

# The pairs each target genuinely regulates. Assigned round-robin so every target has some, and no
# gene is regulated by so many targets that its control group disappears.
true_pairs <- do.call(rbind, lapply(seq_len(n_targets), function(k) {
  genes <- gene_ids[((k - 1) * opts$genes_per_target + seq_len(opts$genes_per_target) - 1) %%
                      n_genes + 1]
  data.frame(grna_target = targets[k], response_id = genes, stringsAsFactors = FALSE)
}))

# Draw counts. Each cell's mean is the gene mean, scaled down where that cell is perturbed for a
# target that really regulates the gene.
effect <- matrix(1, nrow = n_genes, ncol = n_cells, dimnames = list(gene_ids, cell_ids))
target_membership <- vapply(targets, is_perturbed_for, logical(n_cells))
for (i in seq_len(nrow(true_pairs))) {
  g <- match(true_pairs$response_id[i], gene_ids)
  hit <- which(target_membership[, true_pairs$grna_target[i]])
  effect[g, hit] <- 1 - opts$true_knockdown
}
log_step("gRNAs per cell: median ", stats::median(n_grnas_per_cell),
         ", range ", min(n_grnas_per_cell), "-", max(n_grnas_per_cell),
         "; perturbed cells per target: median ", stats::median(colSums(target_membership)))

size_factor <- exp(rnorm(n_cells, 0, 0.25))
counts <- matrix(
  rnbinom(n_genes * n_cells,
          mu = as.vector(effect * outer(gene_mean, size_factor)),
          size = rep(gene_size, times = n_cells)),
  nrow = n_genes, dimnames = list(gene_ids, cell_ids)
)
response_matrix <- as(as(counts, "dgCMatrix"), "CsparseMatrix")

## THROUGH SCEPTRE =================================================================================
#
# Built with sceptre's own API rather than by fabricating slots. The point of this object is that it
# behaves like a real one, and only sceptre knows what that means -- a hand-assembled object would
# drift from the library the moment the pin moves.

extra_covariates <- data.frame(batch = factor(sample(c("b1", "b2"), n_cells, replace = TRUE)))

# Every (target, gene) combination, so the object carries both true and null pairs -- the null ones
# are what make a discovery threshold meaningful.
discovery_pairs <- expand.grid(
  grna_target = targets,
  response_id = gene_ids,
  stringsAsFactors = FALSE
)[, c("grna_target", "response_id")]

log_step("Importing into sceptre (", nrow(discovery_pairs), " candidate pairs)")

log_step("MOI: ", opts$moi, ", resampling mechanism: ", opts$resampling_mechanism)

so <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = opts$moi,
  extra_covariates = extra_covariates
)

so <- set_analysis_parameters(
  sceptre_object = so,
  discovery_pairs = discovery_pairs,
  side = "left",
  grna_integration_strategy = "union",
  resampling_mechanism = opts$resampling_mechanism
)

log_step("Assigning gRNAs")
so <- assign_grnas(so, method = "thresholding", threshold = 1)

log_step("Running QC")
so <- run_qc(so)

log_step("Running discovery analysis")
so <- run_discovery_analysis(so, parallel = FALSE)

## VALIDATE AGAINST THE PIPELINE'S OWN CONTRACT ====================================================
#
# Checked here rather than left for the first pipeline step to discover: a generator that emits an
# object the pipeline rejects is worse than no generator, because the failure surfaces somewhere
# else entirely.

out <- opts$out
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(so, out)
log_step("Wrote ", out, " (", format(structure(file.size(out), class = "object_size"),
                                     units = "auto"), ")")

cat("\n>>> validating against the pipeline's requirements\n")
check <- read_sceptre_object(out, require_discovery = TRUE)
pairs <- qc_passing_pairs(check)
threshold <- discovery_threshold(check)
n_sig <- sum(check@discovery_result$significant, na.rm = TRUE)
n_precomp <- length(check@response_precomputations)

cat("  QC-passing pairs:        ", nrow(pairs), " of ", nrow(discovery_pairs), "\n", sep = "")
cat("  significant discoveries: ", n_sig, "\n", sep = "")
cat("  discovery threshold:     ", format(threshold, digits = 6), "\n", sep = "")
cat("  response precomputations:", n_precomp, "\n")
cat("  true pairs planted:      ", nrow(true_pairs), "\n", sep = "")

if (nrow(pairs) < n_targets) {
  stop("Only ", nrow(pairs), " pairs passed QC, fewer than the ", n_targets, " targets. Raise ",
       "--n-cells so each target has enough perturbed cells.", call. = FALSE)
}
if (length(check@response_precomputations) == 0) {
  stop("@response_precomputations is empty, so no dispersions can be derived.", call. = FALSE)
}

cat("\nOK -- this object satisfies every requirement read_sceptre_object() enforces.\n")
