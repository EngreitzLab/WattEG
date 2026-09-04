#!/usr/bin/env Rscript
#
# Guard against silent breakage from a sceptre upgrade.
#
# This pipeline reaches into sceptre's *internal* S4 slots rather than going through its
# public API: it swaps the response matrix for simulated counts, narrows the discovery
# pairs to one perturbation, and reads the cached negative-binomial precomputations. None
# of that is part of sceptre's documented interface, so a new sceptre release can change
# it without a deprecation warning and the pipeline would keep running while producing
# wrong numbers.
#
# This script fails loudly if anything the pipeline depends on has moved. Run it after
# `pixi run setup` and after any change to the pinned SHA.
#
# Usage:  pixi run check-api

suppressPackageStartupMessages(library(sceptre))

failures <- character(0)
check <- function(ok, label, detail = NULL) {
  if (isTRUE(ok)) {
    cat("  ok    ", label, "\n", sep = "")
  } else {
    cat("  FAIL  ", label, if (!is.null(detail)) paste0("  <- ", detail) else "", "\n", sep = "")
    failures <<- c(failures, label)
  }
  invisible(ok)
}

version <- as.character(utils::packageVersion("sceptre"))
expected_ref <- Sys.getenv("SCEPTRE_REF", unset = NA_character_)
cat("sceptre version: ", version, "\n", sep = "")
cat("library path:    ", dirname(system.file(package = "sceptre")), "\n\n", sep = "")

# --- version matches the pin ---------------------------------------------------------
cat("Version pin\n")
if (!is.na(expected_ref)) {
  check(
    identical(paste0("v", version), expected_ref),
    paste0("installed version matches SCEPTRE_REF (", expected_ref, ")"),
    paste0("installed ", version)
  )
} else {
  cat("  skip   SCEPTRE_REF not set (run via `pixi run check-api` to check the pin)\n")
}

# --- internal S4 slots the pipeline reads or writes -----------------------------------
# Each entry records *why* the pipeline needs the slot, so a future maintainer can tell
# whether a sceptre change is fatal or merely needs a rename here.
required_slots <- c(
  response_matrix              = "replaced with the simulated count matrix on every rep",
  covariate_data_frame         = "cell barcodes (rownames) and the optional 'batch' column",
  covariate_matrix             = "subset alongside covariate_data_frame when sampling controls",
  grna_target_data_frame       = "gRNA -> target mapping; replaces the old grna_groups_table.rds input",
  discovery_pairs_with_info    = "pass_qc filtering, and narrowed to one perturbation per test",
  initial_grna_assignment_list = "per-gRNA cell assignments -> the grna_perts matrix",
  grna_assignments             = "$grna_group_idxs -> the cre_perts matrix and perturbed-cell counts",
  cells_in_use                 = "maps grna_group_idxs positions onto SCE columns",
  response_precomputations     = "cached $theta -> per-gene negative-binomial dispersion",
  functs_called                = "asserts assign_grnas() has been run",
  integer_id                   = "matched against an odm's own @integer_id in attach_response_odm()",
  run_permutations             = "which test the screen ran -- there is no resampling_mechanism slot; set_analysis_parameters() collapses that string to this boolean",
  low_moi                      = "the other implicit route to permutations: with resampling_mechanism at 'default', low MOI resolves to permutations and high MOI to crt"
)

cat("\nInternal sceptre_object slots\n")
present <- methods::slotNames("sceptre_object")
for (slot in names(required_slots)) {
  check(slot %in% present, paste0("@", slot), required_slots[[slot]])
}

# --- exported functions the pipeline calls --------------------------------------------
required_functions <- c(
  "run_discovery_analysis",   # run once per (perturbation, rep) in the power simulation
  "get_result",               # pulls the discovery results out afterwards
  "import_data",              # test-data generator only
  "set_analysis_parameters",  # test-data generator only
  "assign_grnas",             # test-data generator only
  "run_qc"                    # test-data generator only
)

cat("\nExported functions\n")
for (fn in required_functions) {
  check(is.function(tryCatch(get(fn, envir = asNamespace("sceptre")), error = function(e) NULL)), fn)
}

# --- ondisc, for odm-backed response matrices -----------------------------------------
# ondisc is installed from its own pinned commit by `pixi run setup` (src/install_ondisc.R,
# ONDISC_SHA in pixi.toml) -- the same treatment as sceptre, and for the same reason: it is not on
# conda-forge or bioconda, so pixi.lock cannot capture it. It used to be absent by default, which
# is how the odm tests in tests/testthat/test-sceptre_io.R came to skip silently in CI.
# get_response_matrix() (lib/sceptre_io.R) only *uses* it when the object it is handed is
# odm-backed, but the API surface below is unexported/slot-level in the same sense as the
# sceptre_object slots above, so it is checked here rather than left to fail at first use.
cat("\nondisc (exercised when a sceptre object is odm-backed)\n")
ondisc_available <- requireNamespace("ondisc", quietly = TRUE)
check(ondisc_available, "ondisc is installed",
      "run `pixi run setup`, which installs it from ONDISC_SHA via src/install_ondisc.R")

if (ondisc_available) {
  check(
    is.function(tryCatch(get("initialize_odm_from_backing_file", envir = asNamespace("ondisc")),
                         error = function(e) NULL)),
    "initialize_odm_from_backing_file() is exported",
    "attach_response_odm() in lib/sceptre_io.R depends on this to reattach an odm after readRDS()"
  )
  odm_slots <- methods::slotNames("odm")
  check("integer_id" %in% odm_slots, "odm class has an @integer_id slot",
        "attach_response_odm() checks this against the sceptre object's own @integer_id, to catch a mismatched --response-odm before it silently reads the wrong file")
  check(
    existsMethod("[", signature(x = "odm", i = "ANY", j = "missing", drop = "missing")),
    "odm defines a `[` method",
    "materialize_odm_response_matrix() reads one gene row at a time via m[g, ]"
  )
}

# --- the local patches in patches/ ----------------------------------------------------
# The installed sceptre is the pinned commit plus the patches in patches/, applied by
# bin/install_sceptre.R. Their whole point is that they are inert unless used, which also
# means their absence is silent: an unpatched sceptre would reject
# `grna_precomputations = ...` as an unused argument, but only at the moment the simulation
# calls it -- after prepare_sim_input.R and split_pairs.R have already run. Check here instead.
cat("\nLocal patches (patches/*.patch, applied by install_sceptre.R)\n")

patch_stamp <- tryCatch(utils::packageDescription("sceptre")$LocalPatches, error = function(e) NULL)
check(
  !is.null(patch_stamp) && nzchar(patch_stamp) && !identical(patch_stamp, "none"),
  "DESCRIPTION records LocalPatches",
  if (identical(patch_stamp, "none")) {
    "this is a STOCK build (SCEPTRE_PATCH_DIR was empty); re-run `pixi run setup`"
  } else {
    "installed sceptre carries no patch fingerprint; re-run `pixi run setup`"
  }
)
if (!is.null(patch_stamp) && nzchar(patch_stamp)) cat("         ", patch_stamp, "\n", sep = "")

check(
  is.function(tryCatch(get("compute_grna_precomputations", envir = asNamespace("sceptre")),
                       error = function(e) NULL)),
  "compute_grna_precomputations() is exported",
  "0001-reuse-grna-precomputation.patch is not applied"
)
check(
  "grna_precomputations" %in% names(formals(sceptre::run_discovery_analysis)),
  "run_discovery_analysis() accepts `grna_precomputations`",
  "0001-reuse-grna-precomputation.patch is not applied"
)

# --- column names inside discovery_pairs_with_info -----------------------------------
# Verified against sceptre source (R/pairwise_qc_functs.R): the slot is built with columns
# response_id, grna_group, n_nonzero_trt, n_nonzero_cntrl, and pass_qc is added by
# run_qc(). Note that `grna_group` holds the *target* only under the "union" gRNA
# integration strategy; prepare_sim_input.R asserts that strategy explicitly, because
# under "singleton" it holds individual gRNA ids and the pairs would be wrong.
cat("\nAssumptions verified against sceptre source, not testable without a dataset\n")
cat("  note   @discovery_pairs_with_info columns: response_id, grna_group,",
    "n_nonzero_trt, n_nonzero_cntrl, pass_qc\n")
cat("  note   @grna_target_data_frame columns: grna_id, grna_target\n")
cat("  note   @grna_assignments$grna_group_idxs names == perturbation targets (union strategy)\n")
cat("         These are exercised end-to-end by `pixi run pipeline-test`.\n")

# --- verdict --------------------------------------------------------------------------
cat("\n")
if (length(failures) > 0) {
  cat("FAILED: ", length(failures), " check(s) did not pass.\n", sep = "")
  cat("The pinned sceptre version is not compatible with this pipeline. Either revert the\n")
  cat("pin in pixi.toml, or update the affected scripts and this check together.\n")
  quit(save = "no", status = 1)
}
cat("All checks passed. sceptre ", version, " is compatible with this pipeline.\n", sep = "")
