#!/usr/bin/env Rscript
#
# Run the power simulation for one split, one effect size, one chunk of reps.
#
# For each perturbation target in the split, and for each rep, this simulates a count matrix
# under a specified effect size and asks sceptre whether it would have called the association.
# The fraction of reps in which it would is the power, computed downstream by compute_power.R.
#
# The unit of work is (split, effect size, rep chunk). Chunking reps is what makes this
# schedulable: the cost is one run_discovery_analysis() call per (target, rep), which on sample1
# is 2,798 targets x 100 reps = 279,800 calls per effect size. Trapping all 100 reps inside one
# task -- as the original did -- meant no amount of --cores could help.
#
# Usage:
#   run_power_simulation.R --sim-input sim_input.rds --sceptre-template sceptre_template.rds \
#     --pairs split_01.tsv --grna-targets grna_targets.tsv --effect-size 0.15 \
#     --reps 20 --rep-offset 0 --seed 1 --out sim_01_es0.15_chunk1.tsv

suppressPackageStartupMessages({
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
  make_option("--pairs", type = "character", default = NULL, dest = "pairs",
              help = "One split from split_pairs.R (columns grna_target, response_id)."),
  make_option("--grna-targets", type = "character", default = NULL, dest = "grna_targets",
              help = "grna_targets.tsv from prepare_sim_input.R."),
  make_option("--null-precomputations", type = "character", default = NULL,
              dest = "null_precomputations",
              help = paste("null_precomputations.rds from fit_null_models.R: the per-gene null",
                           "model fitted on a null simulation of each replicate. Without it,",
                           "sceptre refits the null model inside every call, which is faithful but",
                           "costs 4.3x. See docs/status.md.")),
  make_option("--effect-size", type = "double", default = NULL, dest = "effect_size",
              help = paste("Effect size as a *fractional decrease* in expression, e.g. 0.15 for",
                           "a 15%% knockdown. Converted internally to a relative expression",
                           "level of 1 - effect_size.")),
  make_option("--reps", type = "integer", default = NULL, dest = "reps",
              help = "Number of simulation reps to run in this chunk."),
  make_option("--rep-offset", type = "integer", default = 0L, dest = "rep_offset",
              help = paste("Number of reps already covered by earlier chunks. The reported `rep`",
                           "column is rep_offset + 1..reps, keeping it unique across chunks",
                           "[default %default].")),
  make_option("--guide-sd", type = "double", default = 0.13, dest = "guide_sd",
              help = paste("Standard deviation of the per-gRNA effect size around the target",
                           "effect size, i.e. guide-to-guide variability [default %default].",
                           "Was hardcoded at 0.13 in the original.")),
  make_option("--n-control-cells", type = "integer", default = NULL, dest = "n_control_cells",
              help = paste("Sample this many control cells per target instead of using every",
                           "non-perturbed cell. Unset means use all of them, which is the",
                           "original behaviour. This is the single largest runtime lever: on",
                           "sample1 the median target has 396 perturbed cells against 430,684",
                           "controls, so every rep simulates 586,309 cells to test ~396.")),
  make_option("--cell-batches", type = "character", default = NULL, dest = "cell_batches",
              help = paste("Column of col_data to stratify control sampling by, so the controls",
                           "match the perturbed cells' composition (e.g. batch_factor). Only",
                           "meaningful with --n-control-cells.")),
  make_option("--seed", type = "integer", default = NULL, dest = "seed",
              help = "RNG seed. Required: results are stochastic and must be reproducible."),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per (pair, rep)."),
  make_option("--no-grna-precomp-reuse", action = "store_false", default = TRUE,
              dest = "grna_precomp_reuse",
              help = paste("Refit the gRNA null model inside every call instead of once per",
                           "target. sceptre regresses perturbation status on the covariates to get",
                           "the probabilities it draws synthetic assignments from; that fit does",
                           "not depend on the counts, so it is identical for every rep of a",
                           "target. Reuse is on by default and is exact -- this flag exists to",
                           "demonstrate that, by showing the p-values are unchanged.")),
  make_option("--gc-every", type = "integer", default = 0L, dest = "gc_every",
              help = paste("Call gc() every N reps. 0 disables it [default %default]. The",
                           "original called gc() on every rep to contain allocation churn that",
                           "the reduced copying in draw_counts() largely removes."))
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Simulate power for one split x effect size x rep chunk."
))
require_options(opts, c("sim_input", "sceptre_template", "pairs", "grna_targets",
                        "effect_size", "reps", "seed", "out"))

if (opts$effect_size <= 0 || opts$effect_size >= 1) {
  stop("--effect-size must be a fractional decrease strictly between 0 and 1 (got ",
       opts$effect_size, ").", call. = FALSE)
}
if (opts$reps < 1) stop("--reps must be at least 1.", call. = FALSE)
if (!is.null(opts$cell_batches) && is.null(opts$n_control_cells)) {
  stop("--cell-batches only applies when --n-control-cells is set.", call. = FALSE)
}

init_seed(opts$seed)

# The config gives a percentage decrease; the simulation multiplies expression by a relative
# level. This conversion was inline at sceptre_power_analysis.R:43.
relative_expression <- 1 - opts$effect_size

## LOAD ============================================================================================

total_started <- proc.time()[["elapsed"]]

log_step("Loading inputs")
started <- proc.time()[["elapsed"]]
sim <- readRDS(opts$sim_input)
validate_sim_input(sim)
template <- readRDS(opts$sceptre_template)
log_resources("load inputs", started)

# Which test this screen runs, read off the template rather than assumed. Logged on every task so
# a simulation output can always be traced to the test that produced it, and used just below to
# decide whether the gRNA precomputation cache is worth building at all.
analysis_mode <- sceptre_analysis_mode(template)
log_step("Resampling mechanism: ", analysis_mode$resampling_mechanism,
         " (MOI: ", analysis_mode$moi, ")")

use_grna_precomp <- isTRUE(opts$grna_precomp_reuse) && !analysis_mode$run_permutations
if (isTRUE(opts$grna_precomp_reuse) && analysis_mode$run_permutations) {
  log_step("Skipping the gRNA precomputation: it is a CRT-path optimisation and this screen ",
           "uses permutations, which would build the cache once per target and discard it.")
}

split_pairs <- read_tsv_file(opts$pairs, required_columns = c("grna_target", "response_id"))
grna_map <- read_tsv_file(opts$grna_targets, required_columns = c("grna_id", "grna_target"))

# --- the per-gene null model ------------------------------------------------------------------
# sceptre skips fitting a gene's null model whenever @response_precomputations already holds an entry
# for it. Three states are possible and only two are wanted, so they are resolved here rather than
# left to whatever the template happens to carry:
#
#   --null-precomputations given   fits from fit_null_models.R, one set per replicate. Correct and
#                                  cheap; this is the intended path.
#   no flag, empty template slot   sceptre refits inside every call. Faithful, 4.3x dearer.
#   no flag, populated slot        the bug this guard exists for. The template used to inherit 272
#                                  fits from the *real* discovery analysis, so simulated counts were
#                                  tested against real-data coefficients -- which measurably
#                                  understates power. It is silent, so it is refused.
null_precomp <- NULL
if (!is.null(opts$null_precomputations)) {
  null_bundle <- readRDS(opts$null_precomputations)
  if (!is.list(null_bundle) || is.null(null_bundle$precomputations)) {
    stop("--null-precomputations ", opts$null_precomputations, " is not a bundle from ",
         "fit_null_models.R.", call. = FALSE)
  }
  if (!identical(as.integer(null_bundle$seed), as.integer(opts$seed))) {
    stop("--null-precomputations was fitted with --seed ", null_bundle$seed, " but this run uses ",
         opts$seed, ". The null models would not correspond to these replicates.", call. = FALSE)
  }
  wanted <- as.character(opts$rep_offset + seq_len(opts$reps))
  absent <- setdiff(wanted, names(null_bundle$precomputations))
  if (length(absent) > 0) {
    stop("--null-precomputations covers replicates ",
         paste(range(null_bundle$reps), collapse = "-"), " but this chunk needs ",
         paste(range(as.integer(wanted)), collapse = "-"), "; missing ", length(absent), ".",
         call. = FALSE)
  }
  null_precomp <- null_bundle$precomputations
  log_step("Null models: ", length(null_precomp), " replicates from ",
           opts$null_precomputations)
  template@response_precomputations <- list()
} else if (length(template@response_precomputations) > 0) {
  stop("The sceptre template carries ", length(template@response_precomputations),
       " inherited @response_precomputations from the real discovery analysis. Testing simulated ",
       "counts against real-data coefficients understates power (see docs/status.md). Either pass ",
       "--null-precomputations from fit_null_models.R, or re-run prepare_sim_input.R, which now ",
       "clears the slot.", call. = FALSE)
} else {
  log_step("Null models: refitted inside every call (no --null-precomputations; ~4.3x cost)")
}

# The split file says *which* targets and genes this task owns; the sceptre-format pair rows
# (with n_nonzero_trt / n_nonzero_cntrl / pass_qc, which sceptre needs) come from the object
# itself rather than being reconstructed.
all_pairs <- template@discovery_pairs_with_info
pair_key <- intersect(c("grna_group", "grna_target"), colnames(all_pairs))[1]
all_pairs <- all_pairs[which(all_pairs$pass_qc), , drop = FALSE]

targets <- unique(split_pairs$grna_target)
log_step("Split covers ", length(targets), " targets / ", nrow(split_pairs), " pairs")
log_step("Effect size ", opts$effect_size, " (relative expression ", relative_expression,
         "), reps ", opts$rep_offset + 1L, "-", opts$rep_offset + opts$reps)
if (!is.null(opts$n_control_cells)) {
  log_step("Sampling ", opts$n_control_cells, " control cells per target",
           if (!is.null(opts$cell_batches)) paste0(" stratified by ", opts$cell_batches) else "")
}

## SIMULATE ========================================================================================

results <- vector("list", length(targets) * opts$reps)
result_idx <- 0L
target_timings <- numeric(0)

for (target in targets) {
  target_started <- proc.time()[["elapsed"]]

  target_pairs <- all_pairs[all_pairs[[pair_key]] == target, , drop = FALSE]
  if (nrow(target_pairs) == 0) {
    message("Skipping target '", target, "': no QC-passing pairs in the sceptre object.")
    next
  }
  pert_genes <- intersect(split_pairs$response_id[split_pairs$grna_target == target],
                          target_pairs$response_id)
  target_pairs <- target_pairs[target_pairs$response_id %in% pert_genes, , drop = FALSE]

  pert_guides <- grna_map$grna_id[grna_map$grna_target == target]
  if (length(pert_guides) == 0) {
    stop("No gRNAs map to target '", target, "' in ", opts$grna_targets, ".", call. = FALSE)
  }

  # --- cells for this target -------------------------------------------------------------------
  # Seed the per-target random work (control-cell sampling, and picking one gRNA per cell in
  # create_guide_pert_status) from a target-level key. rep = 0 is reserved for this setup, with
  # reps numbered from 1. Without this the setup would inherit whatever RNG state the previous
  # target's rep loop happened to leave behind, which differs depending on how reps were chunked.
  set.seed(derive_seed(opts$seed, target, 0L, opts$effect_size))

  pert_object <- if (is.null(opts$n_control_cells)) {
    pert_input(sim, target, level = "cre_perts")
  } else {
    pert_input_sampled(sim, target, level = "cre_perts",
                       n_ctrl = opts$n_control_cells, cell_batches = opts$cell_batches)
  }
  pert_status <- pert_object$col_data$pert
  n_pert_cells <- sum(pert_status == 1)
  if (n_pert_cells == 0) {
    message("Skipping target '", target, "': no perturbed cells.")
    next
  }

  # --- guide-level perturbation status (constant across reps) ----------------------------------
  grna_perts <- as(pert_object$perts$grna_perts, "CsparseMatrix")
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts, pert_guides)

  # create_guide_pert_status() returns perturbed-then-control order; this permutation puts the
  # effect-size matrix back into cell order. Hoisted out of the rep loop: the original rebuilt it
  # by matching cell barcodes on every single rep.
  restore_cell_order <- order(cell_order(pert_status))

  gene_object <- subset_genes(pert_object, target_pairs$response_id)
  effect_sizes <- stats::setNames(
    rep(relative_expression, length(gene_object$genes)), gene_object$genes
  )

  # --- sceptre object, everything except the response matrix set once --------------------------
  target_template <- template
  target_template@discovery_pairs_with_info <- target_pairs

  if (!is.null(opts$n_control_cells)) {
    # With sampled controls the object's cell indexing no longer matches the simulated matrix, so
    # it has to be rebuilt. Note the covariates are reordered to the *simulated column order*
    # (perturbed cells first): the original subset them with `rownames %in% colnames(...)`, which
    # preserves the original order instead and would have misaligned covariates with counts. That
    # never surfaced because n_ctrl was hardcoded to FALSE.
    kept <- match(gene_object$cells, rownames(target_template@covariate_data_frame))
    if (anyNA(kept)) {
      stop("Could not align sampled cells to @covariate_data_frame rownames.", call. = FALSE)
    }
    target_template@covariate_data_frame <- target_template@covariate_data_frame[kept, , drop = FALSE]
    target_template@covariate_matrix <- target_template@covariate_matrix[kept, , drop = FALSE]
    target_template@cells_in_use <- seq_along(kept)
    target_template@grna_assignments$grna_group_idxs[[target]] <- seq_len(n_pert_cells)
  }

  # --- the gRNA null model, fitted once for this target ----------------------------------------
  # sceptre regresses perturbation status on the covariates and draws the synthetic assignments
  # from the fitted probabilities. Both inputs are fixed for the whole rep loop -- the covariate
  # matrix and cells_in_use come from the template, and the assignments from grna_group_idxs -- and
  # neither depends on the simulated counts, so the fit is identical on every rep. Unpatched sceptre
  # refits it inside all `opts$reps` calls; `grna_precomputations` hands it the same answer instead.
  # See patches/0001-reuse-grna-precomputation.patch.
  #
  # Note this must come *after* the sampled-control block above, which rewrites both the covariate
  # matrix and grna_group_idxs.
  #
  # Only the fit is reused. crt_index_sampler_fast() still runs per call, so the synthetic
  # assignments are redrawn every rep and the resampling distribution is untouched.
  #
  # The cache holds regression coefficients, not per-cell fitted probabilities: one target's
  # probabilities would be one double per cell, and the whole point is that this stays negligible
  # next to the count matrix. sceptre reconstructs the probabilities from them bit-identically.
  #
  # Gated on the CRT path as well as on the flag. The patch injects `grna_precomputations` into
  # args_to_pass only on the CRT branch -- the permutation branch shares the same do.call and was
  # deliberately given no ignored formal -- so on a permutations screen this cache is built once
  # per target and then discarded. Computing it anyway is strictly wasted work, and the only way
  # to avoid it used to be passing --no-grna-precomp-reuse, which nobody has a reason to suspect
  # they need. `use_grna_precomp` is resolved once, outside the target loop.
  grna_precomp <- if (use_grna_precomp) {
    compute_grna_precomputations(target_template, analysis = "discovery_analysis",
                                 grna_targets = target)
  } else {
    NULL
  }

  for (rep_local in seq_len(opts$reps)) {
    rep_id <- opts$rep_offset + rep_local

    # Seed per (target, rep, effect size) rather than once per task, so the estimate for a given
    # pair and rep does not depend on how the work was partitioned across splits or rep chunks.
    # See derive_seed() in lib/cli.R.
    set.seed(derive_seed(opts$seed, target, rep_id, opts$effect_size))

    es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                        gene_effect_sizes = effect_sizes, guide_sd = opts$guide_sd)
    es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status,
                                       gene_effect_sizes = effect_sizes)
    es_mat <- es_mat[, restore_cell_order, drop = FALSE]

    counts <- draw_counts(gene_object, es_mat)

    sceptre_use <- target_template
    sceptre_use@response_matrix <- list(as_sceptre_response_matrix(counts, report_density = FALSE))

    # This replicate's null models, fitted on a null simulation of the same replicate. Passing the
    # full set rather than this target's genes is deliberate: sceptre looks entries up by
    # response_id, so the extra ones are inert, and subsetting would cost a match() per rep for no
    # benefit. Left empty, sceptre refits each gene here instead.
    if (!is.null(null_precomp)) {
      sceptre_use@response_precomputations <- null_precomp[[as.character(rep_id)]]
    }

    # run_discovery_analysis() emits a "consider parallel = TRUE" note on every call, and
    # print_progress = FALSE does not silence it -- it is gated on `parallel` alone. Left unchecked
    # that is one line per (target, rep): 279,800 per effect size on sample1.
    #
    # Both wrappers are needed, and which one does the work depends on the sceptre version. Up to
    # v0.10.3 the note was cat()ed to stdout, which capture.output() catches; sceptre 0.99.0 turned
    # most cat() calls into message(), which goes to stderr and slips straight past
    # capture.output(). Keeping both means the pin can move either way without silently
    # reintroducing 279,800 lines of log.
    # `grna_precomputations` is only passed when there is a cache to pass. Unpatched sceptre has no
    # such argument and would reject it as unused, so building the argument list conditionally is
    # what lets --no-grna-precomp-reuse run against a stock sceptre -- which is how the patch's
    # effect can be isolated from the sceptre version's. check_sceptre_api.R still asserts the
    # patch is present, so this is graceful degradation, not a silent fallback.
    discovery_args <- list(sceptre_object = sceptre_use, parallel = FALSE, print_progress = FALSE)
    if (!is.null(grna_precomp)) {
      discovery_args$grna_precomputations <- grna_precomp
    }
    suppressMessages(invisible(utils::capture.output(
      sceptre_use <- do.call(run_discovery_analysis, discovery_args)
    )))
    discovery_result <- get_result(sceptre_object = sceptre_use,
                                   analysis = "run_discovery_analysis")

    discovery_result$num_pert_cells <- n_pert_cells
    discovery_result$rep <- rep_id
    discovery_result$effect_size <- opts$effect_size

    result_idx <- result_idx + 1L
    results[[result_idx]] <- discovery_result

    rm(es_mat, counts, sceptre_use)
    if (opts$gc_every > 0L && rep_local %% opts$gc_every == 0L) gc(verbose = FALSE)
  }

  elapsed <- proc.time()[["elapsed"]] - target_started
  target_timings <- c(target_timings, elapsed)
  log_step(sprintf("target %s: %d pairs, %d cells, %d reps in %.1fs (%.2fs/rep)",
                   target, nrow(target_pairs), length(gene_object$cells), opts$reps,
                   elapsed, elapsed / opts$reps))
}

## COMBINE AND WRITE ===============================================================================

results <- results[seq_len(result_idx)]
if (result_idx == 0L) {
  stop("No targets produced results for this split.", call. = FALSE)
}
combined <- do.call(rbind, results)

# Keep only what a reader needs. sceptre's get_result() returns more than this, and at 100
# replicates x 34,886 pairs x 6 effect sizes the extra columns were 39 % of a 3 GB output:
#
#   fold_change                   2^log_2_fold_change -- the same number twice
#   se_fold_change                read by nothing downstream
#   significant                   sceptre's own call at ITS threshold, not the discovery threshold
#                                 this pipeline tests against. compute_power.R recomputes it from
#                                 p_value and log_2_fold_change, so keeping it invites the wrong
#                                 column being believed
#   average_expression_all_cells  a per-GENE constant repeated once per replicate; summarize_power.R
#                                 now joins it from sim_input.rds, where it lives once
#
# Kept deliberately, despite also being constant per pair or per file:
#
#   effect_size                   compute_power.R refuses an input that mixes effect sizes, and that
#                                 check reads this column. 3 % of the output for a guard against
#                                 silently averaging power across knockdown levels
#   pass_qc, n_nonzero_trt/cntrl  diagnostics. The first thing to look at when a pair's power is
#                                 surprising, and there is nowhere else to recover them from
keep <- c("grna_target", "response_id", "p_value", "log_2_fold_change", "rep", "effect_size",
          "num_pert_cells", "pass_qc", "n_nonzero_trt", "n_nonzero_cntrl")
present <- intersect(keep, colnames(combined))
dropped <- setdiff(colnames(combined), present)
if (length(dropped) > 0) {
  log_step("Dropping ", length(dropped), " column(s) nothing downstream reads: ",
           paste(dropped, collapse = ", "))
}

# as.data.frame() FIRST. sceptre's get_result() returns a data.table, and `dt[, present]` does not
# mean "select the columns named in `present`" there -- data.table evaluates the expression in the
# frame's own scope, so it looks for a column literally called "present" and errors. Coercing makes
# the base-R idiom mean what it says, and nothing downstream needs a data.table: write.table takes
# either.
combined <- as.data.frame(combined)[, present, drop = FALSE]

write_tsv_file(combined, opts$out)
log_step("Wrote ", nrow(combined), " rows to ", opts$out)

if (length(target_timings) > 0) {
  message(sprintf("[timing] per-target: min %.1fs median %.1fs max %.1fs | per-rep median %.2fs",
                  min(target_timings), stats::median(target_timings), max(target_timings),
                  stats::median(target_timings) / opts$reps))
}
log_resources("TOTAL", total_started)
