#!/usr/bin/env Rscript
#
# Install ondisc from a pinned commit into the project-local R library.
#
# ondisc backs out-of-core (`odm`) response matrices. lib/sceptre_io.R reconnects one after a
# readRDS() round trip and materializes it, and tests/testthat/test-sceptre_io.R exercises that
# path -- those tests skip entirely when ondisc is absent, which is how the odm code shipped with
# no CI coverage in the first place.
#
# Same story as sceptre, and the same handling: not published on conda-forge or bioconda, so it
# cannot be captured in pixi.lock, and an unpinned install could silently change behaviour. The pin
# lives in pixi.toml ([activation.env] ONDISC_SHA) so the environment definition and the pin stay
# in one place. Unlike sceptre there are no local patches, so this is the same script without the
# patch machinery.
#
# ondisc's DESCRIPTION Imports -- Rhdf5lib, data.table, readr, dplyr, crayon, Matrix, Rcpp -- are
# all in pixi.toml and installed by pixi before this runs, which is why `dependencies = FALSE` is
# safe here. Note `import(Rhdf5lib)` is a full NAMESPACE import, so Rhdf5lib is a RUNTIME
# dependency, not a build-only one like BH.
#
# This is a source install: it compiles ondisc's C++ against Rhdf5lib's static HDF5, so run it in a
# job or an `sh_dev` shell, never on a login node.
#
# Usage:  pixi run setup
#         Rscript src/install_ondisc.R          # inside an activated environment
#
# Environment:
#   R_LIBS_USER     target library (set by pixi)
#   ONDISC_SHA      commit to install (set by pixi)
#   ONDISC_REF      human-readable version for messages (set by pixi)
#   ONDISC_TARBALL  optional path to an already-downloaded source tarball for the pinned SHA,
#                   for installing on a node with no outbound network

lib <- Sys.getenv("R_LIBS_USER")
sha <- Sys.getenv("ONDISC_SHA")
ref <- Sys.getenv("ONDISC_REF", unset = "(unnamed)")
tarball_override <- Sys.getenv("ONDISC_TARBALL", unset = "")

if (!nzchar(lib)) {
  stop("R_LIBS_USER is not set. Run this via `pixi run setup` so the pixi environment is active.")
}
if (!nzchar(sha)) {
  stop("ONDISC_SHA is not set. It is defined in pixi.toml under [activation.env]; ",
       "run this via `pixi run setup`.")
}

dir.create(lib, recursive = TRUE, showWarnings = FALSE)

# Install into the project-local library only, so nothing is written into the shared pixi
# environment and a re-solve cannot leave a stale ondisc behind.
.libPaths(c(lib, .libPaths()))

installed <- tryCatch(as.character(utils::packageVersion("ondisc")), error = function(e) NA_character_)
installed_desc <- suppressWarnings(
  tryCatch(utils::packageDescription("ondisc"), error = function(e) NULL)
)
if (!is.list(installed_desc)) {
  installed_desc <- NULL
}

if (!is.na(installed) && identical(installed_desc$RemoteSha, sha)) {
  message("ondisc ", installed, " (", substr(sha, 1, 8), ") is already installed in ", lib, ".")
  message("Nothing to do. Delete ", file.path(lib, "ondisc"), " to force a reinstall.")
  quit(save = "no", status = 0)
}

message("Installing ondisc ", ref, " (", substr(sha, 1, 8), ") into ", lib, " ...")

work_dir <- tempfile("ondisc-src-")
dir.create(work_dir)
on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

tarball <- file.path(work_dir, "ondisc.tar.gz")
if (nzchar(tarball_override)) {
  if (!file.exists(tarball_override)) {
    stop("ONDISC_TARBALL is set to ", tarball_override, ", which does not exist.")
  }
  message("  using ONDISC_TARBALL=", tarball_override)
  file.copy(tarball_override, tarball)
} else {
  url <- paste0("https://codeload.github.com/timothy-barry/ondisc/tar.gz/", sha)
  message("  downloading ", url)
  status <- utils::download.file(url, destfile = tarball, mode = "wb", quiet = TRUE)
  if (!identical(status, 0L) || !file.exists(tarball)) {
    stop("Failed to download the ondisc source from ", url, ". If this node has no outbound ",
         "network, download the tarball elsewhere and set ONDISC_TARBALL to its path.")
  }
}

utils::untar(tarball, exdir = work_dir)
src_dir <- file.path(work_dir, paste0("ondisc-", sha))
if (!dir.exists(src_dir)) {
  candidates <- list.dirs(work_dir, recursive = FALSE)
  candidates <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(candidates) != 1) {
    stop("Expected exactly one unpacked ondisc source directory in ", work_dir,
         " but found ", length(candidates), ".")
  }
  src_dir <- candidates[[1]]
}

# Record the provenance in DESCRIPTION. R CMD INSTALL copies it verbatim, so these fields survive
# into the installed package and are what the short-circuit above reads on the next run.
desc_path <- file.path(src_dir, "DESCRIPTION")
desc_lines <- readLines(desc_path, warn = FALSE)
desc_lines <- desc_lines[!grepl("^Remote[A-Za-z]+:", desc_lines)]
writeLines(
  c(
    desc_lines,
    "RemoteType: github",
    "RemoteHost: api.github.com",
    "RemoteUsername: timothy-barry",
    "RemoteRepo: ondisc",
    paste0("RemoteRef: ", sha),
    paste0("RemoteSha: ", sha)
  ),
  desc_path
)

utils::install.packages(
  src_dir,
  lib = lib,
  repos = NULL,
  type = "source",
  dependencies = FALSE,  # every dependency is already pinned by pixi.lock
  INSTALL_opts = "--no-multiarch"
)

version <- as.character(utils::packageVersion("ondisc", lib.loc = lib))
message("Installed ondisc ", version, ".")
message("Next: pixi run check-api")
