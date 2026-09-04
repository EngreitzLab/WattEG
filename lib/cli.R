## Shared CLI plumbing, so every bin/ script is a standalone executable.
##
## The old scripts read `snakemake@input`, `snakemake@output`, `snakemake@params` and
## `snakemake@wildcards`, which meant none of them could be run or debugged without invoking
## Snakemake. Each script now takes explicit arguments and works identically whether it is run by
## hand or staged by Nextflow.
##
## Also gone from every script: the save.image() block that dumped the whole workspace into
## RDA_objects/ on every task, and the sink()-based logging. Nextflow captures per-task
## stdout/stderr in .command.log, and sink() actively hurts -- it swallows output when a script
## dies before reaching close(log).

suppressPackageStartupMessages(library(optparse))

#' Directory containing the running script, so lib/ can be sourced by a relative path.
#'
#' Works both when invoked as `Rscript src/foo.R` and when a workflow engine stages the scripts
#' somewhere else.
this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Interactive / sourced use: fall back to the working directory.
  normalizePath(".")
}

#' Candidate locations for lib/, most specific first.
#'
#' The repo layout is `src/<script>.R` with `lib/` at the root, so the usual answer is one level up
#' from the script. Three other shapes have to keep working: a workflow engine that stages every
#' input into a single flat task directory (Nextflow does this), which puts lib/ beside the script or
#' dissolves it entirely; being run from the repo root with no `--file=`, where `this_script_dir()`
#' is the working directory; and a script that lives OUTSIDE this repo entirely.
#'
#' That last one is why `WATTEG_ROOT` is honoured, and it is checked first. Analyses in
#' broadinstitute/WattEG-paper depend on this repo rather than vendoring it, so their scripts sit in
#' a different tree and nothing script-relative can find `lib/`. Note that the sbatch runner
#' executes every script inside apptainer with `--cleanenv`, so a caller has to pass the variable as
#' `APPTAINERENV_WATTEG_ROOT` for it to arrive here at all.
lib_dirs <- function() {
  here <- this_script_dir()
  root <- Sys.getenv("WATTEG_ROOT", unset = "")
  unique(c(if (nzchar(root)) file.path(root, "lib"),
           file.path(dirname(here), "lib"), file.path(here, "lib"), here))
}

#' Source one or more files from lib/.
source_lib <- function(...) {
  dirs <- lib_dirs()
  for (file in c(...)) {
    found <- FALSE
    for (dir in dirs) {
      path <- file.path(dir, file)
      if (file.exists(path)) {
        source(path)
        found <- TRUE
        break
      }
    }
    if (!found) {
      stop("Cannot find library file ", file, " in any of: ",
           paste(dirs, collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

## TABULAR I/O =====================================================================================
## Base R rather than readr: these are small TSVs (a few tens of thousands of rows at most) and it
## keeps readr, tibble, vroom, cpp11, tzdb, bit64 and progress out of the environment.

read_tsv_file <- function(path, required_columns = character(0)) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing_cols <- setdiff(required_columns, colnames(df))
  if (length(missing_cols) > 0) {
    stop(path, " is missing required column(s): ", paste(missing_cols, collapse = ", "),
         ". Found: ", paste(colnames(df), collapse = ", "), ".", call. = FALSE)
  }
  df
}

#' Write a TSV, gzipped when the path ends in .gz.
#'
#' `write.table(file = "x.gz")` does NOT compress -- it writes plain text to a file that merely has
#' a misleading name, and every reader then either fails or silently gets the wrong bytes. A
#' `gzfile()` connection is what actually compresses, and it has to be closed explicitly.
#'
#' Reading needs no counterpart: `read.delim` goes through `file()`, which sniffs the magic number
#' and decompresses transparently, so `read_tsv_file` handles either form unchanged.
write_tsv_file <- function(df, path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  target <- if (grepl("[.]gz$", path)) {
    con <- gzfile(path, "w")
    on.exit(close(con), add = TRUE)
    con
  } else {
    path
  }
  utils::write.table(df, file = target, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE)
  invisible(path)
}

## ARGUMENT HANDLING ===============================================================================

#' Fail unless every named option was supplied.
require_options <- function(opts, names) {
  missing_opts <- names[vapply(names, function(n) is.null(opts[[n]]), logical(1))]
  if (length(missing_opts) > 0) {
    stop("Missing required argument(s): ",
         paste0("--", gsub("_", "-", missing_opts), collapse = ", "),
         ". Run with --help for usage.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Seed the RNG, reporting the value so a run can be reproduced from its log.
#'
#' The original pipeline called no set.seed() anywhere, despite every result coming out of
#' rnbinom(), rnorm() and sample(). Two runs of the same configuration therefore produced
#' different power estimates and there was no way to reproduce a specific one.
init_seed <- function(seed) {
  if (is.null(seed) || is.na(seed)) {
    stop("A seed is required for reproducibility; pass --seed.", call. = FALSE)
  }
  seed <- as.integer(seed)
  set.seed(seed)
  message("RNG seed: ", seed)
  invisible(seed)
}

#' Stable hash of a string to a positive integer usable as an RNG seed.
#'
#' A polynomial rolling hash modulo 2^31-1. Deterministic across platforms and R versions, which
#' matters because it seeds the simulation. Uses only double arithmetic -- intermediate values
#' stay below 2^53 so they are exact, and unlike bitwXor()/bitwAnd() there is no 32-bit signed
#' integer ceiling to trip over. Collision resistance is irrelevant here: distinct keys only need
#' to land on distinct-enough seeds, not cryptographically distinct ones.
stable_hash <- function(key) {
  h <- 0
  for (byte in utf8ToInt(key)) {
    h <- (h * 131 + byte) %% 2147483647
  }
  as.integer(h)
}

#' A seed for one simulation unit, derived from a stable key rather than from the task layout.
#'
#' This matters statistically, not just for tidiness. If each *task* seeded the RNG once, then
#' which random draws a given (target, rep) received would depend on how the work happened to be
#' partitioned -- so changing --n-splits or --reps-per-chunk, both purely computational knobs,
#' would silently change the reported power estimates. Deriving the seed from
#' (base seed, target, rep, effect size) makes every estimate invariant to the partitioning:
#' the same base seed reproduces the same numbers whether you run 1 split or 100.
derive_seed <- function(base_seed, target, rep, effect_size) {
  stable_hash(paste(base_seed, target, rep, effect_size, sep = "|"))
}

#' Stands in for `target` when seeding work that belongs to a replicate rather than to a target.
#'
#' The null-model fits in fit_null_models.R are the case: a null simulation applies no knockdown, so
#' its counts depend on neither the target nor the effect size, and seeding it from a real target's
#' key would both misrepresent that and collide with that target's own stream. A gRNA target can
#' never be named this, since targets are genomic coordinates or control labels.
NULL_FIT_TARGET_KEY <- "__null_fit__"

#' Timestamped progress message on stderr.
log_step <- function(...) {
  message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
}

#' Report elapsed time and peak memory for a step, feeding the benchmark comparison.
log_resources <- function(label, started_at) {
  elapsed <- proc.time()[["elapsed"]] - started_at
  peak_mb <- tryCatch(sum(gc()[, "max used"] * c(56, 8)) / 1024^2, error = function(e) NA_real_)
  message(sprintf("[timing] %s: %.1fs, peak R heap ~%.0f MB", label, elapsed, peak_mb))
  invisible(elapsed)
}
