#!/usr/bin/env Rscript
# Fit the per-gene null model that the simulation tests against, once per replicate.
#
# WHY THIS STEP EXISTS
#
# Before testing a pair, sceptre fits a Poisson GLM of the gene's counts on the cell covariates
# (perform_response_precomputation) -- the null model the perturbed cells are compared against. It
# skips that fit whenever @response_precomputations already holds an entry for the response_id.
#
# sceptre_template.rds used to inherit 272 such entries from the *real* discovery analysis, covering
# all 237 genes with QC-passing pairs, so every simulated count was tested against coefficients
# fitted to real counts. Measured, that understates power: 7 of 265 discovery calls flipped against a
# faithful refit and all 7 flipped the same way. See docs/status.md.
#
# The faithful alternative -- refitting inside every call -- costs 4.3x, because each pair in a call
# is a distinct gene and so the refit scales with pairs. This script buys the fidelity back. A gene's
# null model is fitted on counts with *no knockdown anywhere*, which makes it independent of both the
# target and the effect size, so one fit per (gene, simulation) serves every target the gene pairs
# with (~147 on average) and the whole effect-size sweep. That is why this is a separate step rather
# than something run_power_simulation.R does inline: targets are spread over ~1,000 array tasks, and
# fitting per task would pay for it ~1,000 times instead of once.
#
# Output is small. One entry is 11 named coefficients plus a theta scalar, so 100 simulations x ~272
# genes is a few hundred KB.
#
# Usage:
#   fit_null_models.R --sim-input sim_input.rds --sceptre-template sceptre_template.rds \
#     --grna-targets grna_targets.tsv --reps 100 --seed 20250812 --out null_precomputations.rds

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(sceptre)
})

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
source_lib("sim_input.R", "pert_input.R", "simulate.R", "sceptre_io.R")

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--sim-input", type = "character", default = NULL, dest = "sim_input",
              help = "sim_input.rds from prepare_sim_input.R."),
  make_option("--sceptre-template", type = "character", default = NULL, dest = "sceptre_template",
              help = "sceptre_template.rds from prepare_sim_input.R."),
  make_option("--grna-targets", type = "character", default = NULL, dest = "grna_targets",
              help = "grna_targets.tsv from prepare_sim_input.R."),
  make_option("--reps", type = "integer", default = NULL, dest = "reps",
              help = "Number of simulations to fit. Must cover every rep the simulation will run."),
  make_option("--rep-offset", type = "integer", default = 0L, dest = "rep_offset",
              help = paste("Number of simulations already covered by an earlier chunk, matching",
                           "run_power_simulation.R's --rep-offset [default %default].")),
  make_option("--seed", type = "integer", default = NULL, dest = "seed",
              help = paste("RNG seed. Must be the same --seed the simulation uses, or the null",
                           "models will not correspond to the simulations they are used for.")),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output RDS: a list of @response_precomputations, named by replicate.")
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Fit the per-gene null model on a null simulation of each replicate."
))
require_options(opts, c("sim_input", "sceptre_template", "grna_targets", "reps", "seed", "out"))
if (opts$reps < 1) stop("--reps must be at least 1.", call. = FALSE)

init_seed(opts$seed)

## LOAD ============================================================================================

total_started <- proc.time()[["elapsed"]]

log_step("Loading inputs")
started <- proc.time()[["elapsed"]]
sim <- readRDS(opts$sim_input)
validate_sim_input(sim)
template <- readRDS(opts$sceptre_template)
log_resources("load inputs", started)

grna_map <- read_tsv_file(opts$grna_targets, required_columns = c("grna_id", "grna_target"))

all_pairs <- template@discovery_pairs_with_info
pair_key <- intersect(c("grna_group", "grna_target"), colnames(all_pairs))[1]
all_pairs <- all_pairs[which(all_pairs$pass_qc), , drop = FALSE]

# One pair per gene is enough. A precomputation depends only on the gene's counts and the cell
# covariates, never on which target the pair is with, so testing each gene once produces the same
# coefficients as testing it against all ~147 of its targets -- at 237 tests per call instead of
# 34,886. The pairs are chosen deterministically (first occurrence in the object's own row order) so
# the output does not depend on how the table happens to be sorted.
representative <- all_pairs[!duplicated(all_pairs$response_id), , drop = FALSE]
genes <- representative$response_id
log_step("Genes to fit: ", length(genes), " (one representative pair each, from ",
         nrow(all_pairs), " QC-passing pairs)")

if (length(genes) == 0) {
  stop("No QC-passing pairs in the sceptre object, so there are no null models to fit.",
       call. = FALSE)
}

# The inherited cache is what this step exists to replace. Warn rather than error: prepare_sim_input.R
# clears it, but an older prepared/ directory will still have it, and silently fitting on top of it
# would reproduce exactly the bug being fixed.
if (length(template@response_precomputations) > 0) {
  log_step("NOTE: template carries ", length(template@response_precomputations),
           " inherited precomputations; clearing them before fitting")
  template@response_precomputations <- list()
}

## FIT =============================================================================================

# Every gene, so the result serves any split. Cells come from an arbitrary target's pert_input
# because with complement controls every cell is in use regardless of target -- the null simulation
# applies no knockdown, so which target supplied the cell set does not affect the counts. Asserted
# below rather than assumed.
anchor_target <- representative[[pair_key]][1]
anchor_guides <- grna_map$grna_id[grna_map$grna_target == anchor_target]
if (length(anchor_guides) == 0) {
  stop("No gRNAs map to target '", anchor_target, "' in ", opts$grna_targets, ".", call. = FALSE)
}

set.seed(derive_seed(opts$seed, anchor_target, 0L, 0))
pert_object <- pert_input(sim, anchor_target, level = "cre_perts")
gene_object <- subset_genes(pert_object, genes)
n_cells <- length(pert_object$col_data$pert)

if (!setequal(gene_object$genes, genes)) {
  stop("subset_genes() returned ", length(gene_object$genes), " genes where ", length(genes),
       " were requested.", call. = FALSE)
}

null_template <- template
null_template@discovery_pairs_with_info <- representative

log_step("Fitting ", length(genes), " null models x ", opts$reps, " simulations over ",
         n_cells, " cells")

# No knockdown: relative expression is 1 everywhere. Built once -- it does not vary by simulation, and
# at 237 x 586,309 doubles it is ~1.1 GB, so rebuilding it per simulation would be pure waste.
null_effect <- matrix(1, nrow = length(gene_object$genes), ncol = n_cells)

precomputations <- vector("list", opts$reps)
rep_ids <- integer(opts$reps)

for (rep_local in seq_len(opts$reps)) {
  rep_id <- opts$rep_offset + rep_local
  rep_started <- proc.time()[["elapsed"]]

  # Seeded from the simulation alone -- deliberately not from (target, rep, effect_size) like the
  # simulation, because the whole point is that a null fit belongs to a simulation and not to any
  # target or effect size. NULL_FIT_TARGET_KEY keeps this stream disjoint from every target's stream.
  set.seed(derive_seed(opts$seed, NULL_FIT_TARGET_KEY, rep_id, 0))

  counts <- draw_counts(gene_object, null_effect)

  obj <- null_template
  obj@response_matrix <- list(as_sceptre_response_matrix(counts, report_density = FALSE))
  obj@response_precomputations <- list()

  # run_discovery_analysis() is the only public route to a precomputation: it fits the null model as
  # a side effect and leaves it in the slot. The test results are discarded -- we want the fit.
  #
  # Both wrappers are needed: sceptre's "consider parallel = TRUE" note was cat()ed to stdout up to
  # v0.10.3 but is message()d to stderr from 0.99.0, and capture.output() only catches the former.
  suppressMessages(invisible(utils::capture.output(
    obj <- run_discovery_analysis(sceptre_object = obj, parallel = FALSE, print_progress = FALSE)
  )))

  fitted <- obj@response_precomputations
  missing_genes <- setdiff(genes, names(fitted))
  if (length(missing_genes) > 0) {
    stop("Simulation ", rep_id, ": no null model was produced for ", length(missing_genes),
         " gene(s), including: ", paste(utils::head(missing_genes, 5), collapse = ", "),
         ". A gene without one would silently fall back to being refitted per call.",
         call. = FALSE)
  }

  precomputations[[rep_local]] <- fitted
  rep_ids[rep_local] <- rep_id
  log_resources(paste0("simulation ", rep_id, " (", length(fitted), " genes)"), rep_started)

  rm(counts, obj, fitted)
}

names(precomputations) <- as.character(rep_ids)

## WRITE ===========================================================================================

# Recorded alongside the fits so run_power_simulation.R can refuse a mismatched file rather than
# testing against null models fitted for different simulations or a different dataset.
out <- list(
  precomputations = precomputations,
  seed = opts$seed,
  reps = rep_ids,
  genes = genes,
  n_cells = n_cells
)
saveRDS(out, opts$out)
log_step("Wrote ", opts$out, " (", length(precomputations), " simulations x ", length(genes),
         " genes)")
log_resources("total", total_started)
