---
title: Development
nav_order: 6
---

# Development

## Environment

[pixi](https://pixi.sh) manages everything.

```sh
pixi install              # the environment, from pixi.lock
pixi run -e dev test      # the fast tests
pixi run -e dev format    # ruff format + ruff check --fix
pixi run pipeline-test    # a stub run of the Nextflow DAG
pixi run lint             # the pre-commit checks, read-only
pixi run install-hooks
```

Two environments: `default`, which every task activates, and `dev`, which adds pytest and ruff.

**There is no longer a `build` environment, a `setup` task or a `check-api` task.** Those existed
because sceptre was installed from a pinned commit with a local patch and compiled C++, ondisc from
another pinned commit, and the pipeline read several of sceptre's unexported S4 slots so the pin had
to be verified against them. None of that survives the Python port: pysceptre is a wheel built from
a git dependency pinned to a commit, the pipeline calls its documented entry points, and a version
bump is caught by this repository's own tests rather than by a bespoke checker.

The sections below describe the R implementation, which still exists in `src/*.R` and `lib/*.R` and
is the reference the Python path was validated against. The environment that ran it is on the
`legacy` branch.

### Why `ggplot2` is a build-only dependency

`R CMD INSTALL` enforces a package's DESCRIPTION `Imports` **even for packages absent from its
NAMESPACE**. sceptre's NAMESPACE imports only `Matrix` and `Rcpp`, but its DESCRIPTION lists
`ggplot2`, `cowplot`, `scales`, `dplyr`, `data.table`, `purrr` and `crayon`, so all of them must be
present to build it. At *runtime* only the ones actually called matter — and `ggplot2`, `cowplot` and
`scales` appear solely in `plotting_functions.R`, `qq_plot_helpers.R` and
`write_outputs_to_directory()`, none of which this pipeline calls. Hence the split.

`purrr`, by contrast, *is* needed at runtime: `purrr::flatten()` is called in the precomputation
path that `run_discovery_analysis()` goes through.

`BH` is build-only for a different reason again: sceptre declares it under `LinkingTo`, so it is a
header-only C++ dependency used to compile and never loaded at all. `src/audit_dependencies.R`
reports all three categories separately, and is the check to run before trusting this table.

### No Bioconductor

The pipeline once used `SingleCellExperiment`, but only as a container — a per-gene table, a per-cell
table, two sparse perturbation matrices, and subsetting. Nothing touched `rowRanges` or any genomic
range. `lib/sim_input.R` replaces it with a plain list, which removes sixteen Bioconductor
packages and leaves `nextflow` as the only bioconda dependency.

It also avoids a hard setup failure: `bioconductor-genomeinfodbdata` installs its data through a
conda post-link script, which pixi skips by default, so `GenomeInfoDb` — and therefore
`SingleCellExperiment` — could not be loaded at all without every user first running
`pixi config set --local run-post-link-scripts insecure`.

## sceptre is pinned, and why that matters

sceptre is not on conda-forge or bioconda, so it is installed from a pinned commit
(`SCEPTRE_REF` / `SCEPTRE_SHA` in `pixi.toml`) into a project-local library at `.pixi/rlibs`.

The pin is not housekeeping. This pipeline reaches into sceptre's **unexported S4 slots**: it swaps
the response matrix for simulated counts, narrows the discovery pairs to one target at a time, and
reads the cached negative-binomial precomputations. None of that is a documented interface, so a
sceptre release could change it with no deprecation warning and the pipeline would keep running while
producing wrong numbers.

`src/check_sceptre_api.R` is the guard. It asserts the presence of every slot the pipeline reads:

`@response_matrix`, `@grna_matrix`, `@covariate_data_frame`, `@covariate_matrix`,
`@grna_target_data_frame`, `@discovery_pairs_with_info`, `@initial_grna_assignment_list`,
`@grna_assignments`, `@cells_in_use`, `@response_precomputations`, `@functs_called`

Run it after any change to the pin. All slot access is confined to `lib/sceptre_io.R`, so a
sceptre upgrade breaks one file rather than five.

Two behaviours of sceptre worth knowing before touching the simulation:

- **`load_row()` dispatches on `odm` and `dgRMatrix` with no `else` branch.** Handing it a dense
  matrix returns `NULL` silently instead of erroring. The simulated counts must be converted to
  `dgRMatrix` before being assigned to `@response_matrix`.
- **With `run_permutations = FALSE`** the conditional-randomisation path reads
  `response_matrix@j/@p/@x` directly, which only a `dgRMatrix` has.

## sceptre is also patched

The installed sceptre is the pinned commit **plus** the patches in `patches/`, applied by
`src/install_sceptre.R` before it compiles. Each patch adds an optional argument with a `NULL`
default, so unpatched sceptre behaves identically and nothing in `patches/` changes any result.

| Patch | What it adds |
| --- | --- |
| `0001-reuse-grna-precomputation.patch` | `grna_precomputations` on `run_discovery_analysis()` and `run_power_check()`, and an exported `compute_grna_precomputations()`. Lets the caller supply the gRNA null model instead of having sceptre refit it once per call. |

Why a patch rather than a wrapper: on the CRT path sceptre regresses perturbation status on the
covariates and draws its synthetic assignments from the fitted probabilities. That fit depends only
on the gRNA-to-cell assignments and the covariate matrix, never on the counts, so it is identical
for all `--reps` simulations of a target. But unlike `@response_precomputations` there was no slot to
hand it back in through — `fitted_probabilities` is a local variable in
`crt_glm_factored_out()` — so reusing it needs a signature change.

The cache holds regression **coefficients**, not per-cell fitted probabilities, and each entry
records the `n_cells` it was fitted over. Both choices matter:

- Fitted values are one double per cell. Coefficients are `ncol(covariate_matrix)` doubles, and
  sceptre reconstructs the probabilities from them with `binomial()$linkinv(drop(X %*% coefs))`,
  which is bit-identical — `glm.fit()` computes the `fitted.values` it returns exactly that way from
  the vector it returns as `coefficients`, with a zero offset here.
- `n_cells` is what makes the guard work. The two CRT workhorses fit on *different cell sets* but the
  same covariates, so the coefficient vector is the same length in both and its length alone cannot
  tell a cache built for one from a cache built for the other.

`--no-grna-precomp-reuse` on `run_power_simulation.R` turns the reuse off. Reuse is exact, so the
flag exists to demonstrate that rather than to work around it —
`workflow/slurm_executor/10_grna_precomp_equivalence.sbatch` runs one split both ways and compares
the outputs byte for byte.

Three things to know:

- **The patches are part of the pin.** `src/install_sceptre.R` refuses to install if `patches/` is
  empty, and `src/check_sceptre_api.R` asserts that `compute_grna_precomputations()` is exported and
  that `run_discovery_analysis()` accepts `grna_precomputations`. Without that second check an
  unpatched sceptre would fail only once the simulation reached its first
  `run_discovery_analysis()` call, long after `prepare_sim_input.R` and `split_pairs.R` had run.
- **Bumping `SCEPTRE_SHA` means regenerating the patches.** They are applied with no fuzz and no
  offset search (`lib/apply_patch.R`), so a hunk that no longer matches at the exact line it
  claims is a hard error rather than a hunk relocated to somewhere it happens to fit.
- **`patch(1)` is not used.** Everything runs inside a stock `ubuntu:22.04` container, which ships
  neither `patch` nor `git`, so `lib/apply_patch.R` applies the diff in R.

**The patch stays local.** Upstreaming it to Katsevich-Lab/sceptre was considered and dropped: the
patch is small, pinned to a known SHA, asserted by `check_sceptre_api.R`, and re-applied
automatically on install, so carrying it costs little — whereas a PR would mean tracking review on
someone else's schedule while the pin moves underneath it. Nothing is implemented upstream and
nothing is planned.

## Repository layout

```
src/                      standalone executables, one per pipeline step
lib/                      shared code, sourced by a path relative to the calling script
  cli.R                   argument parsing, TSV I/O, seeding
  sim_input.R             the simulation-input container and its accessors
  pert_input.R            per-target cell selection
  simulate.R              count simulation and guide-level effect sizes
  sceptre_io.R            the only file that touches sceptre internals
  apply_patch.R           unified-diff applier, for patching sceptre at install time
patches/                  local patches to the pinned sceptre, applied by install_sceptre.R
config/config.yml         pipeline parameters
assets/samplesheet.csv    example samplesheet
docs/                     this documentation (published to GitHub Pages)
.githooks/pre-commit      formatting and convention checks
```

`lib/*.R` is sourced via a path resolved from `--file=` in `commandArgs()`, so every script works
identically whether run directly or staged onto `PATH` by a workflow engine.

## Conventions, enforced by the hook

`pixi run install-hooks` points `core.hooksPath` at `.githooks`. On commit, staged files are
**auto-fixed** for:

- CRLF and lone CR → LF
- tabs → spaces, at the width conventional for the language: **2 for R and YAML**, 4 otherwise
  (source files only; `.tsv`/`.csv` are data and `Makefile` needs tabs)
- trailing whitespace
- missing final newline

and **rejected** for things that cannot be fixed without changing meaning:

- any file over **512 KB** — pipeline inputs and outputs belong outside git; `results/` and `*.rds`
  are gitignored, and this is the backstop
- filenames with uppercase in the stem (the extension is exempt, so `compute_power.R` is fine)
- camelCase or PascalCase R identifiers

The identifier check deliberately allows `SCREAMING_SNAKE` constants (`PERT_LEVELS`) and dotted S3
methods (`print.sim_input`), both of which are correct R style.

`pixi run lint` runs the same checks across every tracked file but is **read-only**: it reports what
it *would* fix and exits non-zero, so it is safe in CI or against a dirty working tree. Use
`./.githooks/pre-commit --all --fix` to actually apply the whitespace fixes repo-wide.

R code is indented with **2 spaces**, following tidyverse and Google R style. YAML uses 2 as well
(and forbids tabs outright); shell, Groovy/Nextflow, Python and JSON use 4. The hook expands tabs to
the matching width so a stray tab does not end up fighting the surrounding style, but it does not
reindent existing code.

## Documentation

`docs/` is published to GitHub Pages by `.github/workflows/pages.yml`. Pages are ordinary markdown
with `title` and `nav_order` front matter.

The theme is pinned to `just-the-docs@v0.3.3` because GitHub's `jekyll-build-pages` action uses the
`github-pages` gem (Jekyll 3.9); just-the-docs v0.4+ requires Jekyll 4 and will not build. If the
Pages build ever breaks, replacing `remote_theme` with `theme: jekyll-theme-primer` in
`docs/_config.yml` is a zero-dependency fallback.

## Known gaps

A fuller account, with the commands to run the outstanding comparisons on a cluster, is in
[Status and handoff]({{ site.baseurl }}{% link status.md %}).

- **Workflow orchestration.** The five steps are complete and run standalone; a Nextflow workflow
  with a SLURM profile is not yet in the repository. `config/config.yml` already carries the
  parameters it will consume.
- **Two-stage simulation allocation** is documented but not orchestrated. The scripts support it
  today via `--rep-offset`.
- **Equivalence with the pre-refactor pipeline** is proven for the upstream statistics — size
  factors, normalised means and raw means are bit-identical, and dropping `@grna_matrix` was verified
  to leave discovery results unchanged — but the simulation itself has not been compared draw-for-draw
  against the old code. Seeding moved inside the simulation loop, so that comparison has to be
  distributional rather than exact.
- **The legacy `R/` directory** is still present so the old Snakemake pipeline remains runnable for
  comparison. `pixi run lint` reports its non-conforming identifiers; both go away when it is removed.
