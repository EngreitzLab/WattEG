#!/usr/bin/env Rscript
# Merge per-chunk null-model bundles from fit_null_models.R into one file.
#
# Fitting all simulations in one process takes ~20 minutes per simulation, so ~33 hours for 100. The
# fits are independent, so the work is run as an array over simulation chunks and merged here.
#
# The merge is where the invariants get checked, because a silently wrong bundle is worse than a
# missing one: run_power_simulation.R looks each simulation up by name, so a gap, a duplicate or a
# gene set that differs between chunks would surface as some tasks refitting and others not --
# exactly the inconsistency this whole mechanism exists to remove.
#
# Usage:
#   merge_null_models.R --inputs 'dir/chunk_*.rds' --reps 100 --out null_precomputations.rds
#   merge_null_models.R --inputs a.rds,b.rds --reps 100 --out null_precomputations.rds

suppressPackageStartupMessages({
  library(optparse)
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

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--inputs", type = "character", default = NULL, dest = "inputs",
              help = paste("Chunk files from fit_null_models.R: either a comma-separated list or a",
                           "single glob (quote it so the shell does not expand it).")),
  make_option("--reps", type = "integer", default = NULL, dest = "reps",
              help = paste("Total simulations expected after merging. Required: it is the only way",
                           "to notice that a chunk failed and its simulations are simply absent.")),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output RDS, in the same shape a single fit_null_models.R run produces.")
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Merge per-chunk null-model bundles into one."
))
require_options(opts, c("inputs", "reps", "out"))

paths <- if (grepl("[*?[]", opts$inputs)) {
  sort(Sys.glob(opts$inputs))
} else {
  trimws(strsplit(opts$inputs, ",", fixed = TRUE)[[1]])
}
paths <- paths[nzchar(paths)]
if (length(paths) == 0) {
  stop("--inputs matched no files: ", opts$inputs, call. = FALSE)
}
missing_paths <- paths[!file.exists(paths)]
if (length(missing_paths) > 0) {
  stop(length(missing_paths), " input file(s) do not exist, including: ",
       paste(utils::head(missing_paths, 3), collapse = ", "), call. = FALSE)
}
log_step("Merging ", length(paths), " chunk(s)")

## MERGE ===========================================================================================

precomputations <- list()
seed <- NULL
genes <- NULL
n_cells <- NULL

for (path in paths) {
  bundle <- readRDS(path)
  if (!is.list(bundle) || is.null(bundle$precomputations)) {
    stop(path, " is not a bundle from fit_null_models.R.", call. = FALSE)
  }

  # Every chunk must describe the same dataset and the same base seed, or the simulations are not
  # comparable and the merged file would be a mixture of two experiments.
  if (is.null(seed)) {
    seed <- bundle$seed
    genes <- bundle$genes
    n_cells <- bundle$n_cells
  } else {
    if (!identical(as.integer(bundle$seed), as.integer(seed))) {
      stop(path, " was fitted with seed ", bundle$seed, " but an earlier chunk used ", seed, ".",
           call. = FALSE)
    }
    if (!identical(bundle$genes, genes)) {
      stop(path, " covers ", length(bundle$genes), " genes where an earlier chunk covered ",
           length(genes), ", or they are in a different order. The chunks are not comparable.",
           call. = FALSE)
    }
    if (!identical(bundle$n_cells, n_cells)) {
      stop(path, " was fitted over ", bundle$n_cells, " cells but an earlier chunk used ", n_cells,
           ".", call. = FALSE)
    }
  }

  duplicated_reps <- intersect(names(bundle$precomputations), names(precomputations))
  if (length(duplicated_reps) > 0) {
    stop(path, " repeats simulation(s) already merged: ",
         paste(utils::head(duplicated_reps, 5), collapse = ", "),
         ". Overlapping --rep-offset ranges would silently keep only one of them.", call. = FALSE)
  }

  precomputations <- c(precomputations, bundle$precomputations)
  log_step("  ", basename(path), ": simulations ",
           paste(range(as.integer(names(bundle$precomputations))), collapse = "-"))
}

# Sorted numerically, not lexically: "10" sorts before "2" as a string, and while the simulation
# looks entries up by name rather than position, an ordered file is far easier to read.
rep_ids <- sort(as.integer(names(precomputations)))
precomputations <- precomputations[as.character(rep_ids)]

## CHECK ===========================================================================================

if (length(precomputations) != opts$reps) {
  stop("Merged ", length(precomputations), " simulations but --reps says ", opts$reps,
       ". A chunk is missing; re-run the failed array tasks rather than proceeding.", call. = FALSE)
}
expected <- seq_len(opts$reps)
if (!identical(rep_ids, expected)) {
  gaps <- setdiff(expected, rep_ids)
  stop("Merged simulations are not 1..", opts$reps, "; missing: ",
       paste(utils::head(gaps, 10), collapse = ", "), call. = FALSE)
}

sizes <- vapply(precomputations, length, integer(1))
if (any(sizes != length(genes))) {
  short <- names(sizes)[sizes != length(genes)]
  stop("Simulation(s) ", paste(utils::head(short, 5), collapse = ", "), " hold ",
       paste(utils::head(sizes[short], 5), collapse = "/"), " null models where ", length(genes),
       " were expected. A gene without one is silently refitted per call.", call. = FALSE)
}

# The fits must differ between simulations: each is fitted on its own null draw, so identical
# coefficients would mean the simulation seed is not reaching the simulation.
if (length(precomputations) >= 2) {
  first_gene <- names(precomputations[[1]])[1]
  a <- precomputations[[1]][[first_gene]]$fitted_coefs
  b <- precomputations[[2]][[first_gene]]$fitted_coefs
  if (isTRUE(all.equal(a, b))) {
    stop("Simulations 1 and 2 have identical coefficients for ", first_gene,
         ": the null simulation is not being redrawn per replicate.", call. = FALSE)
  }
  log_step("Simulation 1 vs 2, ", first_gene, ": max |coef diff| = ",
           signif(max(abs(a - b)), 4))
}

saveRDS(list(precomputations = precomputations, seed = seed, reps = rep_ids,
             genes = genes, n_cells = n_cells), opts$out)
log_step("Wrote ", opts$out, " (", length(precomputations), " simulations x ", length(genes),
         " genes)")
