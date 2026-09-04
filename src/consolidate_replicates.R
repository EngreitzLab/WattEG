#!/usr/bin/env Rscript
#
# Consolidate one effect size's per-split simulation output into a single Parquet file.
#
# WHY THIS STEP EXISTS
#
# The per-replicate output is the only thing power can be RE-derived from: subsampling replicates to
# study a reduced design, bootstrapping the reference, or re-thresholding all read it rather than
# re-simulating. That makes it worth keeping, and worth keeping in a form that is cheap to read.
#
# As 1,000 gzipped TSVs it was neither. Reading one effect size cost 30-90 s of parsing, and the
# reduced-design scan does that 78 times -- most of an hour, for data that never changes. One
# Parquet file reads the columns a caller actually wants in a second or two, and drops the published
# output from 6,000 files to 6.
#
# nanoparquet rather than arrow, deliberately: arrow pulls in 52 packages, including the AWS, Azure
# and Google Cloud SDKs, to provide a cloud filesystem layer this pipeline has no use for.
# nanoparquet is self-contained and adds one. The cost is no Arrow datasets, no Hive partitioning
# and no predicate pushdown -- none of which matter when the unit is one file per effect size.
#
# WHY THE SPLITS ARE STILL WRITTEN FIRST
#
# 1,000 tasks write concurrently and no single-file format supports that. Per-task files are also
# what make the run resumable and preemption-tolerant, so they stay -- as work-directory
# intermediates that are no longer published. Consolidation cannot lose data: the splits survive in
# the work directory until Nextflow cleans it.
#
# Usage:
#   consolidate_replicates.R --simulations <dir or comma-separated list> \
#     --out replicates_es0.15.parquet

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

option_list <- list(
  make_option("--simulations", type = "character", default = NULL, dest = "simulations",
              help = paste("Per-replicate outputs of run_power_simulation.R: a directory, or a",
                           "comma-separated list. Both .tsv and .tsv.gz are accepted.")),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output .parquet, one row per (pair, replicate)."),
  make_option("--compression", type = "character", default = "zstd", dest = "compression",
              help = paste("Parquet codec [default %default]. zstd compresses better than snappy",
                           "at similar read speed; 'uncompressed' is available for debugging."))
)
opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Consolidate per-split simulation output into one Parquet file."
))
require_options(opts, c("simulations", "out"))

if (!requireNamespace("nanoparquet", quietly = TRUE)) {
  stop("The `nanoparquet` package is required. It is pinned in pixi.toml; run this inside the ",
       "pixi environment.", call. = FALSE)
}

## LOAD ============================================================================================

paths <- trimws(strsplit(opts$simulations, ",", fixed = TRUE)[[1]])
paths <- unlist(lapply(paths, function(p) {
  if (dir.exists(p)) sort(list.files(p, pattern = "[.]tsv([.]gz)?$", full.names = TRUE)) else p
}), use.names = FALSE)
if (length(paths) == 0) {
  stop("--simulations matched no files: ", opts$simulations, call. = FALSE)
}

required <- c("grna_target", "response_id", "p_value", "log_2_fold_change", "rep")
started <- proc.time()[["elapsed"]]
combined <- do.call(rbind, lapply(paths, read_tsv_file, required_columns = required))
log_step("Read ", nrow(combined), " rows from ", length(paths), " file(s)")
log_resources("read", started)

# One effect size per file, checked here as well as in compute_power.R. Consolidating a mixture
# would produce a Parquet file that silently averages power across knockdown levels downstream, and
# unlike the TSVs there would be no filename left to notice it by.
if ("effect_size" %in% colnames(combined)) {
  effect_sizes <- unique(combined$effect_size)
  if (length(effect_sizes) > 1) {
    stop("The inputs mix effect sizes (", paste(effect_sizes, collapse = ", "),
         "); consolidate one effect size at a time.", call. = FALSE)
  }
}

if (anyDuplicated(combined[, c("grna_target", "response_id", "rep")])) {
  n_dup <- sum(duplicated(combined[, c("grna_target", "response_id", "rep")]))
  stop(n_dup, " duplicated (grna_target, response_id, rep) row(s): the inputs overlap.",
       call. = FALSE)
}

## WRITE ===========================================================================================

# Character keys become factors first. Parquet dictionary-encodes them either way, but making the
# level set explicit keeps the encoding stable rather than dependent on row order, and `grna_target`
# alone was 15.7 % of the TSV bytes because a ~21-character element name repeated once per
# replicate. There are ~3,000 distinct values against 3.5 M rows.
for (column in c("grna_target", "response_id")) {
  combined[[column]] <- as.factor(combined[[column]])
}

dir <- dirname(opts$out)
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

started <- proc.time()[["elapsed"]]

# `compression` belongs to write_parquet(), NOT to parquet_options().
#
# Passing it to parquet_options() does not error the way a wrong argument name usually would:
# parquet_options() has no `...` but does have `compression_level`, so R's PARTIAL ARGUMENT MATCHING
# silently binds `compression = "zstd"` to `compression_level`, and the failure surfaces as
# "compression_level must be an integer scalar" -- a message about an argument the caller never
# named. Keep this as a direct argument.
codecs <- c("snappy", "gzip", "zstd", "uncompressed")
if (!opts$compression %in% codecs) {
  stop("--compression must be one of: ", paste(codecs, collapse = ", "), call. = FALSE)
}
nanoparquet::write_parquet(combined, opts$out, compression = opts$compression)
log_resources("write parquet", started)

in_bytes <- sum(file.size(paths), na.rm = TRUE)
out_bytes <- file.size(opts$out)
log_step(sprintf("Wrote %s: %d rows, %.1f MB (from %.1f MB across %d files, %.1fx)",
                 opts$out, nrow(combined), out_bytes / 1024^2, in_bytes / 1024^2,
                 length(paths), in_bytes / out_bytes))

# Read it straight back. A Parquet file that cannot be reopened is worse than the TSVs it replaced,
# and this is the step that deletes nothing but is relied on by everything downstream, so the check
# is cheap insurance rather than ceremony.
started <- proc.time()[["elapsed"]]
check <- nanoparquet::read_parquet(opts$out)
if (nrow(check) != nrow(combined)) {
  stop("Wrote ", nrow(combined), " rows but read back ", nrow(check), ".", call. = FALSE)
}
if (!identical(sort(colnames(check)), sort(colnames(combined)))) {
  stop("Column names differ after the round trip.", call. = FALSE)
}
log_resources("read back", started)
log_step("Round trip verified: ", nrow(check), " rows, ", ncol(check), " columns")
