#!/usr/bin/env Rscript
#
# Dump an R sim_input.rds to text, at full precision, so the Python port can be compared to it.
#
# The counterpart is workflow/compare_expression_stats.py. Written at 17 significant digits and
# not with write.csv: write.csv gives 15, which is enough to make a bit-identical port look like
# it agrees only to 6e-12 -- measured, and misleading enough to have sent one investigation the
# wrong way.
#
# Usage:
#   workflow/dump_r_sim_input.R <sim_input.rds> <out_dir>
#
# The reference used for day0 is the sim_input the day0 sweep was run on:
#   WattEG-paper/power_sweep/day0/day0/prepared/sim_input.rds

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: dump_r_sim_input.R <sim_input.rds> <out_dir>", call. = FALSE)
}
sim_path <- args[[1]]
out <- args[[2]]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

sim <- readRDS(sim_path)
cat("genes:", length(sim$genes), " cells:", length(sim$cells), "\n")

# 17 significant digits round-trips a double exactly; formatC is used rather than format() so the
# width is not chosen per vector.
full <- function(x) formatC(as.numeric(x), digits = 17, format = "e")

writeLines(sim$genes, file.path(out, "genes.txt"))
writeLines(sim$cells, file.path(out, "cells.txt"))

for (column in colnames(sim$row_data)) {
  writeLines(full(sim$row_data[[column]]), file.path(out, paste0("row_", column, ".txt")))
}
for (column in colnames(sim$col_data)) {
  values <- sim$col_data[[column]]
  # Factors go out as their labels; the levels are recoverable from the labels and the order is
  # not load-bearing on this side of the comparison.
  writeLines(if (is.numeric(values)) full(values) else as.character(values),
             file.path(out, paste0("col_", column, ".txt")))
}

for (nm in names(sim$perts)) {
  m <- as(sim$perts[[nm]], "CsparseMatrix")
  writeLines(rownames(m), file.path(out, paste0(nm, "_rows.txt")))
  # Triplets, 0-based, so the Python side can rebuild the matrix without knowing R's storage.
  triplets <- summary(m)
  utils::write.table(
    data.frame(i = triplets$i - 1L, j = triplets$j - 1L),
    file.path(out, paste0(nm, "_triplets.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  cat(" ", nm, ":", paste(dim(m), collapse = " x "), "with", nrow(triplets), "nonzeros\n")
}

cat("wrote the reference to", out, "\n")
