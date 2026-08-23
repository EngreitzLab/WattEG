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
