# Load the library under test.
#
# lib/ is sourced rather than installed: this is a pipeline, not a package, and the executables in
# src/ reach it the same way. Resolved from testthat's working directory, which is tests/testthat/.
LIB_DIR <- normalizePath(file.path(testthat::test_path(), "..", "..", "lib"), mustWork = TRUE)

for (f in c("stats.R", "simulate.R", "sim_input.R", "sceptre_io.R")) {
  source(file.path(LIB_DIR, f))
}

# cli.R is sourced separately: it pulls in optparse, which is not needed by most tests and is the
# slowest part of loading.
source_cli <- function() {
  source(file.path(LIB_DIR, "cli.R"))
}

# The shared cross-language fixture (tests/fixtures/README.md): the same design the Python suite
# reads, byte for byte. Returns the pieces create_guide_pert_status() and friends take, plus what
# each cell actually carries, so tests can check the simulation against the design itself.
read_interleaved_fixture <- function() {
  dir <- file.path(testthat::test_path(), "..", "fixtures")
  design <- utils::read.delim(file.path(dir, "interleaved_design.tsv"),
                              colClasses = "character", na.strings = character(0))
  guides <- utils::read.delim(file.path(dir, "interleaved_guides.tsv"),
                              colClasses = "character", na.strings = character(0))
  carried <- lapply(strsplit(design$guides, ";", fixed = TRUE), function(g) g[nzchar(g)])
  grna_perts <- Matrix::sparseMatrix(
    i = match(unlist(carried), guides$grna_id),
    j = rep(seq_along(carried), lengths(carried)),
    x = 1, dims = c(nrow(guides), nrow(design)),
    dimnames = list(guides$grna_id, design$cell)
  )
  list(cells = design$cell,
       pert_status = as.integer(design$pert),
       carried = carried,
       grna_perts = grna_perts,
       pert_guides = guides$grna_id[guides$is_target == "1"],
       other_guides = guides$grna_id[guides$is_target == "0"])
}

# The fixture with one target's guide status built by the production code, as
# run_power_simulation.R builds it once per target.
fixture_target <- function(seed = 1L) {
  fx <- read_interleaved_fixture()
  set.seed(seed)
  fx$status <- create_guide_pert_status(fx$pert_status, fx$grna_perts, fx$pert_guides)
  fx$restore <- order(cell_order(fx$pert_status))
  fx
}
