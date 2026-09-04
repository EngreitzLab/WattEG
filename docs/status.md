---
title: Status and handoff
nav_order: 7
---

# Status and handoff

State of the refactor, written so the work can be picked up by someone else.

**Short version:** the five pipeline steps are complete, run standalone, and have now run on
Sherlock as well as on a laptop. `workflow/slurm_executor/` holds one plain `sbatch` script per step
— deliberately no orchestration layer, so a failure is attributable to the step rather than the
runner. The full simulation array at effect size 0.15 ran on 2026-08-13 as job `38849611`, 1,000
tasks at 100 replicates, 999 of them completing; it is the source of every measured number below.

Since then, and **not yet reflected in the older sections of this document**, three things landed:
`null_fit` is implemented (Step 6a, below), the sceptre pin moved to **0.99.0**, and sceptre now
carries a **local patch** that reuses the gRNA precomputation across replicates. All three are in
flight as of 2026-08-13 afternoon: the null-model array `38879935` is running, its merge `38879939`
and the `null_fit` rerun of the simulation `38882207` are queued behind it.

**The old-vs-new comparison has now run and passed** (job `38916341`) — the refactored pipeline
agrees with the old one on every check, including a zero difference in perturbed cells per pair
across all 34,886 pairs. See Comparison 1 below. Two equivalence checks had already removed the other
confounders: the version bump and the gRNA patch each provably change no number.

Still missing: a test suite (Step 8). **Nextflow (Step 7) is in progress** on the `nextflow` branch —
four of seven processes run, and `PREPARE_SIM_INPUT`'s output is byte-identical to the sbatch
runner's. Step 9 has been **retired**: the gRNA precomputation reuse already collected most of what
it was for.

---

## Done

### The pipeline

Five standalone executables in `src/`, each with `--help`, none referencing `snakemake@`:

```
prepare_sim_input.R -> split_pairs.R -> run_power_simulation.R -> compute_power.R -> summarize_power.R
```

Shared code lives in `lib/` at the repo root (`cli.R`, `sim_input.R`, `pert_input.R`, `simulate.R`,
`sceptre_io.R`, `apply_patch.R`), located from the calling script's own path — each executable
resolves `dirname(dirname(<its own path>))/lib` — so it works both standalone and staged onto `PATH`
by a workflow engine.

**Layout changed on 2026-08-13**, ahead of the Nextflow work: `bin/` became `src/`, and `bin/lib/`
was lifted to a top-level `lib/`. The old Snakemake implementation was removed at the same time and
lives on the **`legacy`** branch. Anything referring to `bin/` predates this.

### Measured improvements

| | Before | After | How verified |
|---|---|---|---|
| Input files required | 4 | **1** | the other three duplicate `@discovery_pairs_with_info`, `@grna_target_data_frame`, `@discovery_result` |
| Per-task read | 427 MB / 3.3s | **59 MB / 1.1s** | timed `readRDS` on both |
| Per-task RAM | 1,774 MB+ | ~250 MB | `object.size` / peak heap |
| Dead code | — | **694 lines removed** | 6 files, none reachable |
| Split imbalance | 1.18× | **1.000** | measured on 2,798 targets |
| Reproducibility | none | seeded and layout-invariant | identical p-values across 1×4 vs 2×2 replicate chunks, and across split layouts. **On the same hardware** — see the note below |
| Environment | 273-line pinned linux-64 conda | 11 direct specs, 2 platforms | `pixi.lock` solves for `linux-64` + `osx-arm64` |

### Correctness fixes

| Issue | Consequence | Status on `sample1` |
|---|---|---|
| Dispersions in a list column with `NULL` holes; `unlist()` shortened the vector and `rnbinom` recycled it | every gene after the first gap simulated with another gene's dispersion | latent — 273 values where 292 expected, but all 239 *tested* genes had precomputations, so it never fired. Now a hard error |
| `max()` over possibly-empty significant discoveries → `-Inf` | zero power reported for every pair, silently | latent — 252 significant discoveries present. Now validated |
| `group_by(rep) %>% group_by(...)` | first grouping discarded | dead code, removed |
| Discovery-results file required but produced by no rule | undeclared fourth manual input | removed; read from the object |
| Sampled-controls path reordered covariates independently of matrix columns | covariates misaligned with counts | latent — the path never ran. Fixed |
| No `set.seed()` anywhere | no run reproducible | fixed; seeds derive from `(seed, target, rep, effect_size)` |
| Per-task seeding | `num_batches` silently changed results | fixed and verified invariant |
| `effect_label` trimmed trailing zeros unconditionally | `0.2` → column `PowerAtEffectSize2` | fixed; also renamed to `power_at_effect_size_20` |
| `sparseMatrix` local shadowing `Matrix::sparseMatrix` | fragile | gone with the old script |
| Size factors permuted across cells before simulating | each cell's library size was paired with a different cell's perturbation status, since `effect_size[i, j]` is indexed by cell | removed; each cell keeps its own size factor. This changes the RNG stream, so a seed no longer reproduces pre-refactor output draw for draw — see [Methods]({{ site.baseurl }}{% link methods.md %}) |

### Two of my own claims that measurement overturned

- **Control-cell sampling is not a usable speed lever.** It costs 21–60 % of power. Off by default.
  Table in [Methods]({{ site.baseurl }}{% link methods.md %}).
- **The speedup from smaller matrices was ~7×, not the ~109× the size reduction suggested.**
  Per-replicate cost is sub-linear in cell count; ~0.3s per replicate is fixed overhead.

### Infrastructure

- `pixi.toml` / `pixi.lock` — three environments (runtime / build / dev), sceptre pinned to
  **0.99.0** (`3ba046b8`, an *untagged* `main` commit — v0.10.3 is still the newest tag, so
  `SCEPTRE_REF` names a DESCRIPTION version, not a tag) **plus a local patch**, both guarded by
  `pixi run check-api`. See the sceptre section below
- `.githooks/pre-commit` — auto-fixes line endings, tabs (2 spaces for R/YAML, 4 otherwise),
  trailing whitespace, final newline; rejects files >512 KB, non-snake_case filenames, camelCase R
  identifiers. `pixi run lint` is the read-only repo-wide version
- Documentation published to GitHub Pages
- `config/config.yml` — parameters only; dead `reps` key removed

---

## Reference numbers

Needed to size cluster runs. Measured, not estimated. The right-hand column is Sherlock, on
`day0_grna20_no_shuffle`, and supersedes the laptop figures wherever they differ.

| Quantity | Laptop (`sample1`) | **Sherlock** |
|---|---|---|
| Genes × cells | 292 × 586,309 | same dataset shape |
| Response matrix | 1,095 MB in RAM, 55.9 % dense | |
| `cells_in_use` | 567,690 of 586,309 | |
| QC-passing pairs / targets | 32,386 / 2,798 | **34,886 / 3,026** (median 9 pairs/target, max 36) |
| Perturbed cells per pair | median 396 (IQR 195–542, max 1,927) | |
| Discovery p-value threshold | 7.71404 × 10⁻⁴ | |
| One `run_discovery_analysis()` call | 3.4–4.8s | **~9.7s** |
| Cost model, per (target, replicate) | — | **1.140s + 0.5561s × pairs** — fitted on the 999 real tasks of `38887744` (`null_fit` + gRNA reuse), and the number to use. Predicts 635 CPU-h against 634 measured. Superseded: 5.076s + 0.5923s × pairs (job `38849611`, before the gRNA reuse), 4.91s + 0.459s × pairs (36-target smoke test), 3.66s + 0.512s × pairs (3 profiled targets) |
| Per-target setup, once | — | **2.8–3.0s** (0.03s per replicate at 100 reps) |
| Where a call spends its time | — | **`glm.fit` 59 %** total / 27 % self, `rnbinom` 7 % — measured *with* the inherited cache, so this is not the skipped per-gene null fit; see the null-model section |
| **Total at 100 replicates** | ~295 CPU-h per effect size | **999 CPU-h `as_is`, 634 CPU-h `null_fit` + gRNA reuse** — both measured, not modelled (`sacct` over jobs `38849611` and `38887744`). **Use 634**; the rest of this table predates the gRNA reuse |
| Per-task wall clock, measured | — | `as_is`: median **57.2 min**, mean 60.0, p95 85.4, max **110.1**. `null_fit` + reuse: median **36 min**, max **61** (4 h requested → 2 h is ample) |
| Peak RAM, simulation task | 1.5 GB → request 4 GB | **median 2.27 GB, max 3.27 GB over 1,000 tasks** (8 GB requested → 6 GB keeps 1.8× headroom) |
| Peak RAM, `prepare_sim_input.R` | 7.7 GB → request 12 GB | 16 GB requested, ~70s |
| Measured peak RAM, every step | — | `prepare` 7.4 GB, `fit_null_models` **7.6 GB** (10 GB requested — only 1.3× headroom), `power_simulation` 2.9 GB, `compute_power` **2.1 GB** (was requesting 32 GB, now 8), `split`/`merge`/`summarize` < 0.15 GB. From job `38978285`'s Nextflow trace |

**The ~295 CPU-hour figure was a laptop measurement and is roughly 3× optimistic.** Sherlock's
cores are slower per-thread than an M-series laptop, and sceptre's per-call cost is untouched by the
refactor, so it dominates. Sized against the old pipeline's own `sacct` record for its 2026-05-14
run, the refactor is still ahead on both axes:

| | Old pipeline | Refactored |
|---|---|---|
| Tasks | 1,004 | 1,000 |
| Mean per task | 78.2 min (range 0.8–134.5) | ~50 min |
| CPU-hours per effect size | 1,308 | **~814–858** |
| Wall clock | **14.9 h** (Snakemake, throttled) | **~50 min** (unthrottled array) |

Both ran 100 replicates on the same dataset, so this is like for like. The CPU win is 1.6×; the
wall-clock win is mostly the unthrottled array rather than the code.

**48 % of the cost is the per-target term** — one `run_discovery_analysis()` call carries the full
586,309-cell bookkeeping however few gene pairs ride along. Targets can only share a call when their
gene sets are disjoint, which is Step 9 below. Consequences worth internalising before optimising:

- Halving the *pairs* saves ~26 %, not 50 %. Halving the *replicates* saves exactly 50 %.
- Filtering pairs for a second pass is much dearer than it looks: the 4,322 pairs whose interval
  straddles power 0.8 (12.4 % of pairs) are spread across 1,820 of 3,026 targets, so a top-up run
  pays 60 % of the per-target overhead to redo 12 % of the pairs — 82 % of its cost is overhead.
- Effect sizes cost the same as each other. At the current 635 CPU-h a six-point sweep is
  **~3,800 CPU-h** (it was ~5,100 before the gRNA reuse), and the 2,900-task QOS budget still allows
  one wave: 480 splits × 6 effect sizes = 2,880 tasks, ~75 min each, so ~1.5 h wall clock.

---

## Settled — which per-gene null model does the simulation test against?

**Verdict: fit the null model on a null simulation of each replicate and reuse it (`null_fit`).**
Today's inherited-cache behaviour understates power, and the fix costs 1.7 %. Measured by job
`38854980` and the threshold check `38861277`; the numbers are below.

Before testing any pair, sceptre fits a Poisson GLM of the gene's counts on the cell covariates
(`perform_response_precomputation`, `precomputation_functions.R:13`) — the null model the perturbed
cells are compared against. It **skips that fit when `@response_precomputations` already holds an
entry for the `response_id`** (`medium_level_functs_v2.R:63`, reached because
`run_outer_regression <- calibration_check || control_group_complement` and our
`control_group_complement` is `TRUE`).

Two measured facts follow:

- `glm.fit` is **27 % self / 59 % total** under `Rprof`, against 7 % for our own `rnbinom` draws, so
  GLM fitting is the single largest cost in a call. **That 59 % is not this precomputation.** The
  profile was taken over the current path, with the inherited cache populated and the precomputation
  therefore skipped — so the 59 % is GLM work sceptre does anyway. Clearing the cache adds the
  precomputation *on top* of it, and the table below measures that: 4.3×, all of it in the per-pair
  term. What the resident 59 % actually is remains unidentified, and it bounds what any pooling
  optimisation such as Step 9 can recover.
- `sceptre_template.rds` **inherits 272 precomputations from the real discovery analysis, and they
  cover all 237 genes that have QC-passing pairs — none are refitted.** So every simulated count is
  tested against coefficients fitted to *real* counts. Measured directly by job `38854980`: 237 of
  237 present, 0 absent. The bias is therefore uniform rather than mixed across genes, which is
  easier to reason about but no more correct — a faithful emulation fits its null model on the data
  it is given, and the real analysis did exactly that.

### The three configurations, measured

Job `38854980` ran all three on 3 targets × 5 replicates, rotating the order within each replicate so
no configuration is always first:

| | What it does | Cost model | CPU-h per effect size |
|---|---|---|---:|
| `as_is` | inherited real-data cache — today's behaviour | 5.85s + 0.839s × pairs | 1,305 |
| `cleared` | refits inside every call, on that call's simulated counts — **the faithful reference** | 8.18s + 5.127s × pairs | 5,656 |
| `null_fit` | fitted once per (gene, replicate) on a null simulation, then reused | 5.51s + 0.891s × pairs | **1,327** |

`cleared` is **4.3×**, and the extra cost sits in the *per-pair* slope, not the intercept — each pair
in a call is a distinct gene, so refitting scales with the genes tested. It is a reference, never a
production candidate. `null_fit` buys its fidelity for **+1.7 % over `as_is`**.

Treat the *ratios* as sound and the absolute CPU-hours as indicative only: this job timed `as_is` at
11.5s on the 5-pair target where the profiling run measured 6.15s, so its clock runs ~1.9× slow —
probably because it holds three caches and a null simulation in memory at once. That is also why the
1,305 CPU-h here disagrees with the ~858 in the table above. The 100-replicate array's own `sacct`
record supersedes both.

### Why `as_is` is not good enough

Aggregate p-value agreement looks reassuring and is misleading:

```
cleared vs null_fit   n=265  max|diff|=0.012  median=1.2e-10  spearman=1.0000
cleared vs as_is      n=265  max|diff|=0.08   median=4.6e-05  spearman=0.9990
```

Power is not aggregate. It is the fraction of replicates with `p < 7.2629e-4`, so only crossings of
that threshold matter, and a median difference of 4.6e-05 is ~6 % of the threshold itself. Counting
the calls instead of correlating the p-values (job `38861277`, `threshold_check.R`):

```
cleared vs null_fit   flips=0 / 265 (0.00 %)   per-pair power identical
cleared vs as_is      flips=7 / 265 (2.64 %)   [cleared-only=7, as_is-only=0]
```

- **`null_fit` is exactly equivalent to the faithful reference** — zero disagreements in 265 calls,
  and per-pair power identical to three decimal places.
- **`as_is` is biased, one-directionally.** All 7 flips go the same way: `cleared` calls the pair
  significant and `as_is` does not, never the reverse (137 of 265 significant against 144). Random
  noise would split both ways — P ≈ 0.016 for 7–0 — so `as_is` **understates power**, by a mean
  |Δpower| of 0.026 on this sample.
- It concentrates where it matters. Only ~37 of 265 p-values sit within a factor of 10 of the
  threshold, so roughly **19 % of near-threshold calls flip**.

Caveat on the magnitude: 3 targets, 53 pairs, 5 replicates. The *direction* is solid and the
equivalence of `null_fit` is unambiguous; the 2.64 % figure is thin and should not be quoted as
precise.

**Superseded on the magnitude.** The full 100-replicate run under both configurations (see
Consequences below) puts the mean shift at **+0.0063**, not the ~0.026 this sample gave — the
direction held, the size did not. The 2.64 % flip rate was measured on 5 replicates, where a single
replicate crossing the threshold moves power by 0.2; at 100 replicates it is 0.8 percentage points
of pairs crossing the certification line.

Note the coefficients themselves differ substantially — inherited theta median 20.0 against 18.8 from
the null fit, median worst-coefficient difference 7.77 — while the p-values barely move. sceptre's
conditional-resampling test is far more robust to the null model than the null model is stable.

### Consequences

- **The effect size 0.15 array (`38849611`) ran `as_is`, and the bias has now been measured directly
  rather than extrapolated.** Both configurations exist at 100 replicates over all 34,886 pairs
  (`power_as_is/` and `power_null_fit/`), and they are exactly paired — same simulated counts, so
  every difference is the null model with zero Monte-Carlo noise between them.

  | | `as_is` | `null_fit` |
  |---|---:|---:|
  | Mean power | 0.5933 | 0.5996 (**+0.0063**) |
  | Certified (`power_ci_low` ≥ 0.8) | 13,094 (37.5 %) | 13,374 (**38.3 %**) |
  | Ambiguous | 4,382 (12.6 %) | 4,408 (12.6 %) |
  | Underpowered | 17,410 (49.9 %) | 17,104 (49.0 %) |

  **The direction is confirmed and the magnitude is ~4× smaller than predicted.** `as_is` understates
  power: 11,649 pairs rise against 3,631 that fall, a 3.2:1 ratio that noise would not produce. But
  the shift is **+0.0063 mean, 0.0100 mean absolute**, against the ~0.03 this document extrapolated
  from 3 targets × 53 pairs × 5 replicates. That extrapolation was too thin, as it said at the time;
  treat 0.0063 as the number and the earlier 0.03 as withdrawn.

  It concentrates exactly where the mechanism says it should — in the transition band, not at the
  saturated ends:

  | `as_is` power band | pairs | mean Δ | % moved |
  |---|---:|---:|---:|
  | 0.00–0.02 (floor) | 4,411 | +0.0010 | 9.7 % |
  | 0.02–0.20 | 4,409 | +0.0094 | 52.7 % |
  | 0.20–0.50 | 4,411 | +0.0121 | 73.4 % |
  | 0.50–0.80 | 5,894 | +0.0130 | 73.5 % |
  | 0.80–0.98 | 8,937 | +0.0047 | 48.6 % |
  | 0.98–1.00 (ceiling) | 6,824 | +0.0003 | 9.0 % |

  Practical consequence: `null_fit` is the correct configuration and should stay the default, but
  **nothing built on the `as_is` numbers was materially wrong** — it moved 280 of 34,886 pairs across
  the certification line, 0.8 percentage points. Keep both outputs.
- The pre-refactor pipeline inherited the same cache, so the **old-vs-new comparison is unaffected**:
  both sides share the bias. Compare them `as_is` against `as_is`.
- Clearing the slot without supplying null fits gives `cleared`, at 4.3×. The two changes must land
  together.

An earlier version of this section put the cost of `cleared` at ~2.4× per call and cited a p-value
divergence of 0.942 with Spearman 0.224 between `as_is` and `null_fit`. Both came from
`profile_call.R`, whose first configuration was labelled "refit, as the pipeline runs today" but
inherited the template's cache and refitted nothing (`line 100` only overwrites the slot when
`precomp` is non-NULL). The 2.4× was `1 / (1 - 0.59)` applied to a `glm.fit` share that is not the
precomputation. The proper experiment gives `as_is` vs `null_fit` max 0.088 and Spearman 0.9989, so
**both numbers are retracted**. Why that script reported 0.942 at all is unexplained; it is
superseded by `cache_experiment.R` and should not be used.

---

## Step 6a — `null_fit`: implemented, running

Not a one-line change, because the null fit has to be **shared across array tasks**. A gene's null
model is target-independent, so one fit per (gene, replicate) serves every target — but targets are
spread over 1,000 splits, so fitting inside each task would pay for it 1,000 times instead of once.
It is also **effect-size independent**: a null simulation has no knockdown, so one set of fits serves
the whole sweep. 100 replicates, not 100 × 6.

What shipped:

1. **`src/fit_null_models.R`** — for each replicate, simulates a null matrix (no knockdown) and runs
   one `run_discovery_analysis()` with `@response_precomputations` empty, keeping only the resulting
   slot. An entry is 11 named coefficients plus a `theta` scalar, so 100 replicates × 272 genes is a
   few hundred KB. Seeded from `(seed, rep)` only — deliberately *not*
   `(seed, target, rep, effect_size)`, since the point is that it depends on neither target nor
   effect size.
2. **`run_power_simulation.R --null-precomputations`** assigns that replicate's entry to
   `@response_precomputations` inside the replicate loop (`line 317`), and validates that the bundle's
   seed and replicate range match the run's.
3. **The inherited slot is cleared in `slim_sceptre_object()`** (`lib/sceptre_io.R:199`), *after*
   `build_dispersion_vector` has read it (`prepare_sim_input.R:259`) — dispersions must keep coming
   from the real data, since they set the noise the simulation is meant to reproduce. Only the null
   model moves. Clearing is what makes an accidental `as_is` impossible to reintroduce.
4. **Steps `02b` (array, one task per replicate), `02b_submit.sh` and `02c` (merge)** in
   `workflow/slurm_executor/`, with `--null-precomputations` threaded through
   `04_run_power_simulation.sbatch`.

**The guard is in.** `run_power_simulation.R:166` errors if the template arrives with a non-empty
`@response_precomputations` and no `--null-precomputations` — rather than silently running `as_is`.
That is the bug that went unnoticed for the entire refactor, and it can no longer recur silently.

**Done as of 2026-08-13.** The null models are built (`38879935`, merged by `38879939` into
`prepared/null_precomputations.rds`), and the `null_fit` rerun of effect size 0.15 completed as
**`38887744`**: 1,000/1,000 tasks, every output at the exact expected row count, every output
carrying a matching `.provenance` sidecar. Results are in `results/refactor/sim_null_fit/es0.15/`.
The first attempt, `38882207`, silently skipped all 1,000 tasks and was cancelled — see the
resume-check entry under Gotchas.

**And it is 1.58× cheaper than the `as_is` run, not 1.7 % dearer.**

| | `as_is` (`38849611`) | `null_fit` (`38887744`) |
|---|---:|---:|
| Tasks completed | 999 / 1000 | **1000 / 1000** |
| Median task | 57.2 min | **36 min** |
| Max task | 110.1 min | **61 min** |
| CPU-h per effect size | 999 | **634** |

The null model itself costs ~1.7 % more, so the saving is not from it — it is the **gRNA
precomputation reuse**, which shipped in the same commit. That retires the "the saving is small"
caveat recorded under the equivalence check below: 1.02× was measured on the single worst-case split
(1 target, 36 pairs, in the `cleared` configuration), and at the median 9 pairs per target the fixed
per-target gRNA fit is a far larger share of the call. **1.58× at production scale is the number to
use**; the step-10 figure measured what it measured, on a split chosen to be pessimal.

Note this also invalidates the cost model in the reference table for any future run: `5.076s +
0.5923s × pairs` was fitted on `as_is` output with no gRNA reuse. Refit it on `38887744` before
sizing the next sweep. `02b_submit.sh` had never run successfully before commit `08b0d04`: `PASSTHROUGH[@]: unbound
variable` under `set -u`, the bash 4.2 empty-array bug already fixed in `04_submit.sh` but never
back-ported. That is why `null_models/` was empty.

Keep the `as_is` output rather than deleting it — it is the only direct measurement of how much the
bias moved real numbers at 100 replicates, and the threshold check that found the bias used 5.

---

## sceptre: pinned to 0.99.0, and patched

Two changes to the dependency, both verified not to move any number.

### The version bump

The pin moved from `v0.10.3` (`9faa4373`) to **0.99.0** (`3ba046b8`, 2026-08-07). There is no
`v0.99.0` tag — v0.10.3 is still the newest tag — so `SCEPTRE_REF` names an untagged commit's
DESCRIPTION version. The slot set is byte-for-byte the same 55 slots (reformatted only), so
`sceptre_object.rds` and `sceptre_template.rds` serialised under 0.10.3 still deserialize, and
0.99.0's NEWS is Bioconductor prep with no statistical change.

Two things the bump forced, both easy to miss:

- **`cat()` became `message()`.** `capture.output()` only catches stdout, so the "consider
  `parallel = TRUE`" note would have leaked ~279,800 lines per effect size. `suppressMessages()` now
  wraps the `run_discovery_analysis()` calls in `run_power_simulation.R` and `fit_null_models.R`;
  both wrappers are kept so the pin can move either way.
- **New Imports `parallelly` and `withr`**, both on our code path. `R CMD INSTALL` enforces Imports,
  so they had to enter `pixi.toml`, forcing a `pixi.lock` re-solve — which must run **inside the
  ubuntu container**, since the host is CentOS 7 (glibc 2.17) and pixi solves against the host's
  virtual `__glibc`. That is `00b_relock_env.sbatch` (job `38877284`); the only package added was
  `r-parallelly 1.48.0`. **Commit `pixi.toml` and `pixi.lock` together.**

Verified not to be a problem: the one RNG change on the association path
(`set.seed(4)` → `withr::with_seed(4, ...)` in `partition_response_ids()`) sits inside
`if (parallel)`, and the pipeline always calls `parallel = FALSE`.

### The patch — reusing the gRNA precomputation

sceptre's CRT path regresses perturbation status on the covariates and draws synthetic treatment
assignments from the fitted probabilities. That fit depends on the gRNA-to-cell assignments and the
covariate matrix — **never on the response counts**. The simulation calls `run_discovery_analysis()`
once per (target, replicate) with only `@response_matrix` changing, so every replicate after the
first refit an identical model.

Unlike `@response_precomputations` there is no slot to hand it back through: `fitted_probabilities`
is a local variable inside the two CRT workhorses. So this ships as
`patches/0001-reuse-grna-precomputation.patch`, adding an optional argument with a `NULL` default —
unpatched behaviour is unchanged. Reuse is **on** by default; `--no-grna-precomp-reuse` turns it off.

Three design points are load-bearing — do not "simplify" them away:

- **The cache holds coefficients, not fitted values.** Fitted values are one double per cell, which
  at this scale is ~2.3 GB held in the caller's session. Reconstruction via
  `binomial()$linkinv(drop(X %*% coefs))` is bit-identical — `glm.fit` computes the `fitted.values`
  it returns exactly that way, and the offset is zero here. Verified: `identical()` TRUE, max abs
  diff 0.
- **Each entry records `n_cells`, and the guard checks it.** The two workhorses fit on different cell
  sets but the same covariates, so the coefficient vector is the same length in both — length alone
  cannot tell a cache built for one path from a cache built for the other.
- **Never cache `crt_index_sampler_fast()`'s output.** The fitted probabilities are deterministic;
  the synthetic index sets drawn from them must be redrawn every call. Caching the indices would
  freeze the resampling distribution across replicates and destroy the CRT null — producing a power
  estimate that is garbage but looks plausible.

A target that is rank deficient *within its own cells* is omitted from the cache with a warning:
`glm.fit` returns `NA` for an aliased coefficient while keeping `fitted.values` valid, so a naive
`X %*% coefs` would yield all-`NA`.

Supporting machinery: `lib/apply_patch.R` is a pure-R unified-diff applier, needed because the
cluster runs everything in a stock `ubuntu:22.04` container that has neither `patch` nor `git`. It is
stricter than `patch(1)` — no fuzz, no offset search, exact context match — so a moved pin fails
loudly. `src/install_sceptre.R` stamps `RemoteSha` and a `LocalPatches` md5 into DESCRIPTION and
short-circuits only when **both** match; `src/check_sceptre_api.R` asserts the patched API
separately, because an unapplied patch would otherwise surface only at the first
`run_discovery_analysis()` call, long after steps 01 and 02.

### Equivalence check `10_grna_precomp_equivalence` — does the patch change results? No. And the saving is small.

Job `38879261`: outputs **byte-identical** between reuse and refit. That was the point of the job and
it is confirmed.

The cost result is weak, and much weaker than the retracted reasoning implied:

```
refit every replicate   620s   (124.0s per call)
reuse per target        608s   (121.6s per call)
speedup                 1.02x  (-1.9% wall clock)
```

`split_0001` is 1 target with 36 pairs, so 5 replicates = 5 calls and reuse removes 4 of 5 gRNA fits:
12s / 4 ⇒ ~3s per gRNA fit. Two caveats pointing opposite ways: 1.9 % on 620s is **within plausible
node-to-node variation**, so 3s is an upper bound from a noisy difference, not a measurement — and
this is the **least favourable target in the dataset**, since 36 pairs is the maximum against a
median of 9, and the run was in the `cleared` configuration, which inflates the per-pair term ~5.75×.
Against the `null_fit` model a fixed ~3s is ~8 % at 36 pairs but ~22 % at the median 9.

**The "up to 2.4×" claim stays retracted. The replacement is 1.58×, measured at production scale**
— see the Step 6a table above, where the full `null_fit` array came in at 634 CPU-h against the
`as_is` array's 999. Both caveats above pointed that way and the second one was right: this split is
the least favourable in the dataset. Do not quote the 1.02× as the saving; quote it as what the
equivalence check happened to cost on one worst-case split, which is not what it was measuring.

### Equivalence check `11_sceptre_version_equivalence` — does the version bump change results? No.

Job `38879265`: stock v0.10.3 against patched 0.99.0, reuse off both sides, **identical** — the bump
does not change any number, so the old-vs-new comparison (Comparison 1 below,
`07_compare_old_new.sbatch`) is not confounded by it.

Both checks are exactly paired: seeds derive from `(seed, target, rep, effect_size)`, which involves
neither the cache nor the sceptre version, so the simulated counts are bit-identical across each pair
and any difference would be attributable. Between them they isolate the two confounders the
old-vs-new comparison would otherwise carry.

### Gotchas found along the way

Recorded because each cost real time and none is discoverable from the code.

- **`vendor/sceptre.tar.gz` is a snapshot of `main`, not the pin**, and is untracked (`vendor/` is
  gitignored). Do not use it as an install source.
- **The login node caps open file descriptors at 256, hard**, and `pixi lock` exhausts it mid-solve
  (`No file descriptors available (os error 24)`). The solve must run in a job — hence
  `00b_relock_env.sbatch`. More generally, nothing heavy belongs on the login node.
- **pixi cannot install R packages from GitHub at all** — conda channels only, and its git support is
  PyPI-only. Hence the separate `pixi run setup` task plus a SHA pin.
- **`.githooks/pre-commit` strips trailing whitespace from source extensions**, which would corrupt a
  `.patch` file's blank context lines. `.patch` is not in its `is_source_file` list, so it is safe
  today — and `apply_patch.R` also treats a zero-length line as blank context. Do not add `.patch`
  to that list.
- **A patch's context only guards the lines it touches.** The API assertions in
  `check_sceptre_api.R` are the layer that actually catches a silently-unapplied patch.
- **`02c`'s dependency is `afterok:38879935,afterany:<step 10>:<step 11>`**, and the `afterany` is
  deliberate: `02c` writes `prepared/null_precomputations.rds`, and the two equivalence checks test
  for that file at runtime, so an unguarded merge would silently flip them from `cleared` to
  `null_fit` depending on scheduling.
- **Step `04`'s resume check skipped all 1,000 tasks of the `null_fit` rerun.** Job `38882207`
  reported 1,000 successes in ~15 s each and recomputed nothing: `out` was
  `${SIM_DIR}/es${ES}/${SPLIT_NAME}.tsv` with no reference to the null-model configuration, and
  "complete" was judged by row count — which is identical under every configuration, since the null
  model does not change how many rows come out. So the `as_is` output from `38849611` satisfied the
  check. The R-level guard at `run_power_simulation.R:166` cannot catch this, because R never
  starts. Worse, the skip is the only thing that *prevented* a second bug: step 04 now passes
  `--null-precomputations` unconditionally, so had it run it would have overwritten the `as_is`
  baseline in place — the same directory `09_compare_null_fit.sbatch` requires as the `as_is` side.

  Fixed two ways. `NULL_MODEL_CONFIG` is now part of the output path
  (`SIM_DIR=${OUTDIR}/sim_${NULL_MODEL_CONFIG}`, `POWER_DIR` likewise), so two configurations can no
  longer address the same file; the `as_is` baseline moved to `results/refactor/sim_as_is/` and is
  named explicitly by `AS_IS_SIM_DIR`. And each output now carries a `.provenance` sidecar recording
  the configuration, the null-bundle md5, reps, seed, guide-sd and effect size — the resume check
  requires it to match, so an output produced differently is redone rather than trusted. A missing
  sidecar means the file predates the check and cannot be vouched for, so it is also redone.
  Resubmitted as `38887744`.
- **Step `03`'s replicate-count check compared against stale output.** It failed with "971 rows, 376
  pairs, 109 with != 2 reps" while the three files the run actually wrote were each exactly 36 pairs
  × 2 reps. The directory held output from an earlier 280-split run, and because `split_pairs.R` pads
  the index to the width of `--n-splits`, those files are named `smoke_split_001_*` against
  `smoke_split_0001_*` — the check's glob could not tell them apart. Fixed in `628d5c9`: step 03 now
  records the outputs it writes and checks only those.

---

## Left to do

### Step 7 — Nextflow + SLURM profile — **in progress**, on the `nextflow` branch

Four of seven processes exist and run: `PREPARE_SIM_INPUT`, `SPLIT_PAIRS`, `FIT_NULL_MODELS`,
`MERGE_NULL_MODELS`. Still to write: `POWER_SIMULATION`, `COMPUTE_POWER`, `SUMMARIZE_POWER`.

Minimal DSL2, not nf-core scaffolding — same principle as `workflow/slurm_executor/`: no framework
layer, so a failure is attributable to the step rather than the runner. **The sbatch executor stays**
until Nextflow reproduces a run; it produced every measurement in this document.

```
main.nf                            the DAG
nextflow.config                    params and profiles
conf/base.config                   measured resources, per process
conf/sherlock.config               SLURM + apptainer + pixi
config/config.yml                  production parameters
config/test.yml                    3 splits x 2 replicates
modules/local/*.nf                 one file per process
workflow/nextflow/run.sbatch       the driver, which runs as a job
```

The DAG, which is **not** the shape planned above — that predates `null_fit`:

```
samplesheet -> PREPARE_SIM_INPUT ---+-> SPLIT_PAIRS -----------------+
                                    |                               |
                                    +-> FIT_NULL_MODELS (x reps)    |
                                        -> MERGE_NULL_MODELS -------+
                                                                    v
                                    POWER_SIMULATION (split x effect size)
                                         -> collectFile by (sample, effect size)
                                         -> COMPUTE_POWER -> SUMMARIZE_POWER
```

`FIT_NULL_MODELS` hangs off `PREPARE_SIM_INPUT` rather than off `SPLIT_PAIRS` because a gene's null
model is fitted on a null simulation and so depends on neither the target nor the effect size. One
set of fits serves every split and every effect size in a sweep — 100 replicates, not 100 × splits ×
effect sizes. Making it a sibling rather than a descendant is what expresses that.

#### The design question that had to be settled first

The sbatch runner does `cd $REPO_ROOT && pixi run`. Nextflow cannot: a task must run in its own work
directory or its outputs are invisible to the engine. Two facts, both verified before building on
them, make it work anyway:

- **`pixi run --manifest-path <repo>/pixi.toml` preserves the working directory** (tested from
  `/tmp`), so the task's outputs land in the task directory.
- **`R_LIBS_USER` is derived from `$PIXI_PROJECT_ROOT`**, not the cwd, so the environment follows
  the manifest.

Every process invocation is built on that pair.

#### Validated at production scale

**The Nextflow runner reproduces the sbatch runner exactly, on a full 1,000-split run.** Job
`38978285`, 2h19m wall clock, all seven processes, 34,886 pairs at 100 replicates:

| | |
|---|---:|
| Pairs matched | 34,886 |
| Pairs whose **power** differs | **0** (max abs difference 0) |
| Pairs whose **certification** flips | **0** |

That is the check that decides whether the two runners can be considered equivalent, and it passes
on the deliverable rather than on an intermediate. Note the per-replicate p-values do still differ
in the far tail between the two runs — see "Reproducibility is bitwise on one node" — but not one of
those differences reaches the power estimate.

Two incidental findings from that run:

- **`results/day0_grna20/sceptre_object.rds` and `results/day0_grna20_no_shuffle/sceptre_object.rds`
  are byte-identical.** The `_no_shuffle` suffix describes how the *old pipeline* processed the
  sample — with the size-factor permutation removed — not a different input. The dataset inventory
  under Step 11 lists them as separate datasets, which is true of their old-pipeline *outputs* and
  false of their inputs. Running the refactored pipeline on both is redundant.
- **Preemption on `owners` is handled.** Three tasks were killed with exit 143 and retried
  automatically by `conf/base.config`'s `errorStrategy`, without intervention. The sbatch runner
  needs a resubmission plus its resume check to do the same.

#### Verified so far

**`PREPARE_SIM_INPUT`'s output is byte-identical to the sbatch runner's** — all five files, both
`.rds` objects included (`sim_input.rds` 16,057,759 bytes; `sceptre_template.rds` 46,520,066). That
is a much stronger result than "it ran": a single difference anywhere in the container, pixi
environment, R version or patched sceptre would perturb the serialised objects. The environment
integration is therefore the same one the measurements came from.

#### Four traps, each of which reported something other than its cause

Recorded because every one cost time and none is guessable from the error text.

1. **`.ifEmpty { error ... }` aborts every run.** Nextflow invokes the closure while *building* the
   DAG, not when the channel proves empty, so a guard for a malformed samplesheet fired
   unconditionally — and, being the first error raised, printed "samplesheet has no rows" in place
   of the real failure further down, masking trap 2 completely. Validate the file eagerly instead.
2. **`take` rejects a `def` local or an `as int` cast.** Both arrive wrapped in a `PojoWrapper` and
   fail with `Missing process or function take([DataflowStream, PojoWrapper])`, which reads like a
   missing operator. A literal or a `params.*` value works; nothing else tried did. Bound the range
   when you build the channel rather than trimming it afterwards where you can.
3. **`-profile test` was silently doing almost nothing.** Nextflow gives `-params-file` precedence
   over profile params, and the driver always passes one, so every test override sharing a name with
   a production parameter was discarded — the run used 2 null chunks of 1 replicate where the
   profile asked for 1 of 2, and `num_replicates` stayed at 100. A "test" run would have been a
   full-size run. Test settings now live in `config/test.yml`, selected with `PARAMS_FILE`.
   **Profiles choose where a run executes; params files choose what it runs.**
4. **A local `def scratch` is shadowed by the `scratch` process directive.** `conf/sherlock.config`
   sets `process.scratch = false`, so `"--bind ${scratch}:${scratch}"` expanded to
   `--bind false:false` and apptainer died with "unable to add false to mount list: destination must
   be an absolute path". It had also silently reached `PIXI_CACHE_DIR` and friends. The variables are
   named `scratchDir`/`groupHomeDir` for that reason — do not tidy them back.

Two smaller ones: the driver must resolve `REPO_ROOT` from `scontrol show job ... Command=` rather
than `$BASH_SOURCE`, because Slurm copies the script into its spool directory; and Nextflow's
`file()` resolves a relative path against `launchDir` and returns an *absolute* path, so
`file(p).isAbsolute()` is always true and cannot be used to test whether the user gave a relative
path.

#### Notes still standing

- The combine step needs no script: `collectFile(keepHeader: true, skip: 1)` replaces the deleted
  `combine_sceptre_power_analysis.R`.
- Resources come from the measured table above, not the old Snakemake guesses (which asked 8 GB for
  simulation and 64 GB for prepare). `conf/base.config` names the source of each.
- `cpus 1` on the simulation process — sceptre runs with `parallel = FALSE` and parallelism comes
  from the fan-out.
- The driver runs as a job, never on a login node: it is a JVM that lives as long as the pipeline.
  `NXF_HOME` and `NXF_WORK` are redirected to `$GROUP_HOME` and `$SCRATCH`, since the defaults are
  `~/.nextflow` and `./work` and `$HOME` is 15 GB.
- `nextflow run . -profile stub -params-file config/test.yml -stub-run` checks the wiring in seconds
  without touching the cluster. The `stub` profile uses `resourceLimits` so the measured requests do
  not have to be rewritten to fit a login node.
- ~~Removing `Snakefile`, `rules/` and `R/` belongs to this step.~~ **Done ahead of it**, on
  2026-08-13, once Comparison 1 had passed and the old implementation had nothing left to prove.
  `Snakefile`, `rules/`, `R/` and `envs/` are gone from this branch and preserved on **`legacy`**.
  `pixi run lint` is now **green** — those files were the only things it flagged.

#### Why not sceptre's own Nextflow pipeline

sceptre now ships one (`timothy-barry/sceptre-pipeline`). It is the right tool for a **real**
discovery analysis and the wrong one for a power simulation. Evaluated against a clone of it:

- **It requires ondisc matrices.** Every script calls
  `read_ondisc_backed_sceptre_object(sceptre_object_fp, response_odm_fp, grna_odm_fp)`, with no
  in-memory path. Using it would mean serialising each simulated matrix to ODM.
- **The simulation needs one object per gene-disjoint batch, not one per replicate** — the same
  constraint as Step 9, since a gene is paired with ~147 targets. That is 100 × 300 = 30,000 objects
  for one effect size, each needing its own ODM.
- **It re-runs everything upstream per object.** The workflow is `set_analysis_parameters ->
  prepare_assign_grnas -> assign_grnas -> combine -> run_qc -> prepare -> analysis`, unconditionally;
  `pipeline_stop` truncates the tail but cannot skip the head. gRNA assignment alone is 43,432 gRNAs
  at their own `2s`/gRNA heuristic, ~24 CPU-hours per object, none of it dependent on the simulated
  response matrix.
- **Its `0.05s` per pair is a walltime heuristic, not a benchmark** — `run_association_analysis_time_per_pair`,
  used to size `--time`. It is not comparable to our measured per-pair cost.

What is worth taking from it: `prepare_association_analyses.R` sorts pairs by `response_id` and pods
on that key, so all ~147 pairs of a gene share one `@response_precomputations` entry — one GLM fit
serving 147 tests. That is the same insight as Step 9 (pool more pairs per call) and it is what led to
the null-model question above.

### Reproducibility is bitwise on one node, not across the cluster

Found on 2026-08-14 while checking the Nextflow runner against the sbatch one on the same three
splits, at the same seed. Two of the three splits came out **byte-identical**. The third did not:

| | |
|---|---|
| Rows compared | 72 |
| `fold_change` differing | 7, max relative **1.2 × 10⁻¹⁵** (~5 ulps) |
| `p_value` differing by >10⁻⁶ relative | **4**, all at p < 10⁻²⁵ |
| Largest single difference | 2.1 × 10⁻⁸¹ against 9.2 × 10⁻⁸⁰ |
| **Replicates changing significance** | **0** |

This is not a different RNG stream — that would perturb every row, not four of seventy-two. It is
last-bit floating-point difference between CPU generations (this cluster spans six; see the node
table under Open questions), amplified by sceptre's skew-normal tail approximation, which is very
sensitive in the far tail and nowhere else. A p-value of 10⁻⁸¹ is zero for every purpose the
pipeline has.

**What the seeding guarantees is the random draws, not the arithmetic.** So:

- The same seed on the same node type reproduces a run bit for bit.
- Across node types, expect agreement to ~15 significant figures on estimates, and agreement on
  every significance call — which is what power is computed from, and therefore what matters.
- Do not write a test, a comparison or an acceptance criterion that demands byte-identity of
  simulation output. `workflow/nextflow/verify_against_sbatch.sbatch` compares numerically and
  asserts zero significance flips, which is the correct bar; an earlier version demanded
  byte-identity and failed on hardware rather than on the pipeline.

### Step 8 — tests and comparisons

1. ~~**Unit tests**~~ — **done**: 62 tests in `tests/testthat/`, run with `pixi run test` or
   `sbatch workflow/slurm_executor/12_unit_tests.sbatch` (R does not belong on a login node). They
   cover the pure functions only — no sceptre, no lab data, no cluster — and every one of them
   corresponds to an entry under "Correctness fixes" above:

   | File | Covers |
   |---|---|
   | `test-stats.R` | `wilson_interval` against published values and at both boundaries, where the normal approximation collapses; `effect_label` pinning the `0.2 → "2"` bug |
   | `test-simulate.R` | `center_effect_size_matrix` putting each gene's perturbed mean on its target; `create_effect_size_matrix` orientation and clamping; `build_dispersion_vector` erroring rather than recycling |
   | `test-seeding.R` | `derive_seed` separating every key component, and being invariant to how replicates are chunked — the property that makes `--n-splits` a purely computational knob |

   Two of these tests failed on first run **because the tests were wrong, not the code**: the
   effect-size matrix is genes × cells, not cells × genes (`run_power_simulation.R:314` reorders
   *columns* by cell), and `init_seed` logs its seed deliberately. Both were corrected against the
   actual contract rather than by loosening the assertion — a test encoding a wrong assumption is
   worse than no test, and the transposed one would have passed a genuinely broken implementation.

   `wilson_interval` and `effect_label` moved to `lib/stats.R` to make this possible; they were
   defined inside `compute_power.R` and `summarize_power.R`, where testing them meant invoking a
   command line.

   Not covered yet: split binning conservation, which is enforced in `SPLIT_PAIRS` and
   `02_split_pairs.sbatch` at runtime rather than by a unit test.
2. ~~**Synthetic test data generator**~~ — **done**: `src/make_test_data.R` builds a sceptre object
   the whole pipeline runs on, in 43 seconds, with no lab data.
   `sbatch workflow/slurm_executor/13_make_test_data.sbatch` builds it and then runs all six steps
   on it end to end. Latest: 360 QC-passing pairs, 27 significant discoveries, threshold 7.3 × 10⁻³,
   30 precomputations, 648 KB.

   Built through sceptre's own API — `import_data` → `set_analysis_parameters` → `assign_grnas` →
   `run_qc` → `run_discovery_analysis` — rather than by fabricating slots, so it cannot drift from
   the pinned library. Two constraints are not obvious and both were found by the object being
   rejected:

   - **It has to plant real knockdowns.** `discovery_threshold()` errors when nothing is
     significant, deliberately, because the old code returned `-Inf` there and reported zero power
     for every pair. Pure noise produces an unusable object, not a boring one.
   - **Cells must carry a *variable* number of gRNAs.** With exactly one each, `grna_n_nonzero` is
     constant, collinear with the intercept, and sceptre rejects the formula as containing
     "redundant information". `moi = "high"` for the same faithfulness reason: it is what makes the
     control group the complement, which is the CRT path production actually uses.

   The object is **not committed** — `tests/data/` is gitignored. It is 648 KB against the
   pre-commit hook's 512 KB ceiling, and regenerating it is the point of having a generator.

   `sceptredata` was considered instead. It is not on conda, not installed, and not even in
   sceptre's `Suggests`, so it would need the same SHA-pin-and-install machinery as sceptre itself;
   and since it ships raw matrices rather than a `sceptre_object`, the API calls above would still
   have to be made by hand. Reasonable, but it removes less than it adds.
3. ~~**Old-vs-new comparison**~~ — **done and passed**, job `38916341`; see Comparison 1 below.
4. **Two-stage replicate allocation** — documented in
   [Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}) but not
   orchestrated. The scripts already support it through `--rep-offset`.

### Step 9 — batch gene-disjoint targets into one sceptre call — **RETIRED, do not build this**

**The gRNA precomputation reuse already collected most of what this was for, by a different
mechanism.** Step 9 amortises the per-call setup by putting many targets in one call; the patch
eliminated the dominant *component* of that setup outright, by caching the fit instead. Refitting the
cost model on `38887744` shows what is left:

| | Per target | Per pair | Per-target term | Step 9 ceiling |
|---|---:|---:|---:|---:|
| `as_is` (job `38849611`) | 5.076s | 0.5923s | 427 CPU-h (42.6 %) | 1.63× |
| **`null_fit` + gRNA reuse** (`38887744`) | **1.140s** | **0.5561s** | **96 CPU-h (15.1 %)** | **1.16×** |

Batching 302,600 calls into ~30,000 would leave ~10 CPU-h of per-target overhead, so it now saves
about **86 CPU-h per effect size out of 635** — against the ~385 it promised. That is not worth a
refactor that touches seeding, splitting and the target loop, and whose two open risks were never
resolved: a 194-gene batch needs ~0.9 GB for the count matrix and as much again for the effect-size
matrix against ~54 MB per target today, and sceptre's per-call cost was only ever shown linear over
5–36 pairs, not the 194 a batch would carry.

The per-pair term is now **85 % of the bill**, and Step 9 explicitly does not move it — a gene is
still simulated once for every target it pairs with. Any further optimisation has to attack per-pair
work, not per-call setup.

The reasoning below is kept because the analysis is still correct and the batching argument is
reusable if the per-target term ever grows back.

---

Originally: the largest single saving left on the table, **up to 1.8×**, worth about 370
CPU-hours per effect size. Costed below because the reasoning is easy to get wrong in both
directions.

**Where the cost is.** The unit of work is one `run_discovery_analysis()` call per (target,
replicate) — 3,026 × 100 = **302,600 calls per effect size** — and each one carries the full
586,309-cell setup regardless of how few gene pairs it covers. That per-call term is 4.91s of the
`4.91s + 0.459s × pairs` model, so it is 413 of the 858 CPU-hours: **half the bill is setup paid
and over**.

**Why the obvious fix does not work.** The tempting version is one simulated object per replicate
covering all 34,886 pairs, so the setup is paid 100 times instead of 302,600. It is not valid. Each
gene is paired with **147 targets on average, up to 299** — only 10 of 237 genes belong to a single
target — and pair (A, g) needs gene g knocked down in A's perturbed cells while pair (B, g) needs it
knocked down in B's. One response matrix cannot hold both. Applying every knockdown at once puts
B-perturbed cells, also knocked down, into A's control group, which shrinks the contrast and
**understates power**. The pre-refactor pipeline was per-target for the same reason; this is not an
artifact of the refactor.

**What does work.** Two targets can share a call whenever their gene sets are disjoint, because a
knockdown on gene `gA` never touches gene `gB`'s row. Cell overlap is fine: a cell with gRNAs for
both A and B is in each one's treatment group for its own gene, and unmodified for the other's.
Greedy first-fit over the real pair table gives **300 batches**, against a lower bound of 299 forced
by the gene that appears in 299 targets — so greedy is essentially optimal here.

| | Now | Gene-disjoint batches |
|---|---:|---:|
| sceptre calls per effect size | 302,600 | **30,000** (300 batches × 100 reps) |
| Per-call overhead | 413 CPU-h | **41 CPU-h** |
| Per-pair work (`rnbinom` draws) | 445 CPU-h | 445 CPU-h |
| **Total** | **858** | **486 CPU-h — 1.77×** |

Batches hold at most 31 targets, 194 pairs and 194 of the 237 genes.

**The per-pair term does not move, and that bounds the win.** Gene `g` is still simulated once for
every target it pairs with, so the total `rnbinom` volume is unchanged. This is 1.8×, not the
10× the call-count ratio suggests.

**`parallel = TRUE` is a different question and probably not the lever.** sceptre's internal
parallelism forks across pairs within a call; the array fan-out already provides that parallelism
across nodes, without fork overhead and without pinning a replicate to one node's cores. It does not
reduce the number of calls, which is where the waste is. Batching and who-schedules-the-parallelism
are orthogonal.

**Three things to check before building it.**

1. **Settled by the array: the per-call overhead is 427 CPU-h and Step 9's ceiling is 1.63×.**
   Fitting `time = reps × (a × targets + b × pairs)` over all 999 completed tasks gives
   **5.076s + 0.5923s × pairs** per (target, replicate) — 1,001 CPU-h predicted against 999 measured,
   so the model is sound. The per-target term is **427 CPU-h, 42.6 % of the bill**.

   Batching 302,600 calls into ~30,000 leaves ~42 CPU-h of overhead, so
   427 + 574 → 42 + 574 ≈ **616 CPU-h, a 1.63× saving worth ~385 CPU-h** per effect size.

   This **refutes the previous revision of this point**, which used a three-target profiled fit
   (`3.66s + 0.512s × pairs`) to argue the intercept was mostly amortizable and the ceiling only
   ~1.5×. The real intercept is 5.076s, close to the original 36-target estimate of 4.91s. That
   three-target fit was thin and this document said so; it was still wrong to revise the headline
   number on it. The slope was the part the smoke test got wrong: 0.5923s against 0.459s, 29 % low.
2. **Memory.** A 194-gene batch is 194 × 586,309 dense doubles ≈ 0.9 GB for the count matrix and as
   much again for the effect-size matrix, against ~54 MB per target today. The 8 GB request may
   hold but has not been tested.
3. **Linearity.** The `0.459s` slope was fitted over 5–36 pairs per call. Batches carry up to
   194, and sceptre's per-call cost has not been shown to stay linear that far out.

**Shape of the change.** A batching function (first-fit over gene sets, ordered deterministically so
layout is reproducible), the target loop in `run_power_simulation.R` iterating batches instead of
single targets, and `split_pairs.R` splitting on batches, not targets. Seeds must stay derived
from `(seed, target, rep, effect_size)` so results remain invariant to the batch layout, exactly as
they are invariant to the split layout today — that invariance is also what makes the two designs
directly comparable pair by pair.

### Step 10 — re-derive the environment — **DONE**

Audited 2026-08-14. The list was in better shape than expected: one package was in the wrong
environment, and nothing else needed to change.

**Our own code** uses only `Matrix`, `sceptre` and `optparse` via `library()`, plus `stats`,
`utils`, `methods`, `tools` and `testthat` via `::`. Everything else in `pixi.toml` is there for
sceptre's benefit, which is where the real question was.

**`r-bh` moved to the build environment.** sceptre declares `LinkingTo: Rcpp, BH` — header-only
Boost, used to compile the C++ and never loaded at runtime. It had been sitting in the runtime
environment that every Nextflow task activates, under a comment describing it as a "NAMESPACE
import", which it is not: sceptre's NAMESPACE imports only `Matrix` and `Rcpp`. That wrong comment
is presumably how it got there.

`src/audit_dependencies.R` is the tool that settled it, and it separates three things that get
conflated:

| | Meaning | For sceptre 0.99.0 |
|---|---|---|
| NAMESPACE imports | loaded on every `library()` call, unconditional | `Matrix`, `Rcpp` |
| `::` references | loaded only if that line executes | 15 packages, incl. `dplyr`, `data.table`, `purrr`, `crayon`, `parallelly`, `withr` |
| `LinkingTo` | compile-time headers, never loaded | `Rcpp`, `BH` |

`BH` was the only declared dependency with no reference anywhere in the package's code. Everything
else earns its place, so the rest of the runtime list stands as documented in `pixi.toml`.

What the audit **cannot** tell you is which of those 15 `::` references lie on *our* path
specifically — that needs coverage instrumentation. It narrows the question to a list short enough
to reason about, and the synthetic object now makes a candidate removal cheap to test by running
the whole pipeline without the package.

**The re-solve was a no-op apart from the intended change** (job `39034832`): four lines removed
from `pixi.lock`, no package added, removed or moved version. `r-base` stayed at
`4.4.3-h35b0bb1_11`. Verified afterwards by job `39037197`, which rebuilt the synthetic object and
ran all six pipeline steps against the new environment — if `BH` were needed at runtime, sceptre
would not load at all.

**A caveat that outlived this step.** Nearly every spec in `[dependencies]` is `"*"`, so the
lockfile is the only thing pinning versions and any re-solve is a gamble on what conda-forge has
published since. This one happened to come back clean; the next one may not. Pinning the specs to
their locked versions would make future re-solves honest, and is worth doing before the next
dependency change rather than during one.

### Step 11 — fit the power curve and run three effect sizes instead of six

`src/fit_power_curve.R` fits each pair's power curve so a sweep needs three effect sizes rather
than six, and makes the minimum detectable effect size continuous rather than snapped to whichever
effect sizes were run. Usage is in [Usage]({{ site.baseurl }}{% link usage.md %}); the model is
derived in [Methods]({{ site.baseurl }}{% link methods.md %}).

**The validation — whether three points really do reproduce six, the covariate model that predicts
the curve for pairs never simulated, and the calibrated prediction intervals — moved to
[broadinstitute/WattEG-paper](https://github.com/broadinstitute/WattEG-paper) on 2026-09-04**, with
the rest of the paper analyses. Short version of what it established, so a pipeline user knows what
the feature is worth: the best three-point design reproduced the six-point per-pair MDES for 92.0 %
of pairs on `day0`, against a reference whose own bootstrapped self-agreement was 91.4 %. Grid
placement is dataset-specific and a one-effect-size pilot is enough to choose it.

## Running the comparisons on the cluster

### Setup

```sh
git clone https://github.com/EngreitzLab/WattEG.git
cd WattEG
pixi install
pixi run setup        # compiles sceptre; needs a compiler and network
pixi run check-api    # must print "All checks passed"
```

`pixi run setup` installs sceptre into `.pixi/rlibs` inside the checkout, so the whole environment is
self-contained and needs no module system. It does require outbound network access to GitHub — if the
compute nodes are offline, run `setup` on a login node first.

### Smoke test before committing CPU-days

Three targets, 2 replicates — should finish in a couple of minutes and proves the environment works.

```sh
mkdir -p prepared splits sim
Rscript src/prepare_sim_input.R \
  --sceptre-object results/sample1/sceptre_object.rds --outdir prepared/

Rscript src/split_pairs.R --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/

Rscript src/run_power_simulation.R \
  --sim-input prepared/sim_input.rds --sceptre-template prepared/sceptre_template.rds \
  --pairs splits/split_001.tsv --grna-targets prepared/grna_targets.tsv \
  --effect-size 0.15 --reps 2 --seed 20250812 --out sim/smoke.tsv
```

### Full run as a SLURM array

`prepare_sim_input.R` once, then one array task per (split, effect size).

**These snippets are illustrative and predate `workflow/slurm_executor/`, which is what actually
runs — prefer it.** In particular they omit `--null-precomputations`, so as written they would run
the `cleared` configuration at ~4.3× the cost: the template's inherited slot is now emptied by
`slim_sceptre_object()`, so sceptre refits every gene's null model inside every call. Add steps `02b`
and `02c` and pass the bundle.

```sh
#!/usr/bin/env bash
#SBATCH --job-name=power-prepare
#SBATCH --mem=12G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
pixi run Rscript src/prepare_sim_input.R \
    --sceptre-object results/sample1/sceptre_object.rds \
    --outdir prepared/
pixi run Rscript src/split_pairs.R \
    --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/
```

```sh
#!/usr/bin/env bash
#SBATCH --job-name=power-sim
#SBATCH --array=1-280
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/sim_%A_%a.out

# Selected by position rather than by formatting the index: split_pairs.R pads the filename to
# the width of --n-splits, so 280 splits give split_001.tsv but 1000 would give split_0001.tsv.
SPLIT=$(ls splits/split_*.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")
for ES in 0.15 0.2; do
    pixi run Rscript src/run_power_simulation.R \
        --sim-input prepared/sim_input.rds \
        --sceptre-template prepared/sceptre_template.rds \
        --pairs "$SPLIT" \
        --grna-targets prepared/grna_targets.tsv \
        --effect-size "$ES" --reps 100 --seed 20250812 \
        --out "sim/split_${SLURM_ARRAY_TASK_ID}_es${ES}.tsv"
done
```

Then per effect size:

```sh
for ES in 0.15 0.2; do
    pixi run Rscript src/compute_power.R \
        --simulations "$(ls sim/*_es${ES}.tsv | paste -sd, -)" \
        --threshold-file prepared/discovery_threshold.txt \
        --out "power_es${ES}.tsv"
done

pixi run Rscript src/summarize_power.R \
    --power power_es0.15.tsv,power_es0.2.tsv --out power_summary.tsv
```

Expect roughly an hour per array task at 10 targets × 100 replicates × 2 effect sizes. Time limit is
set to 3h for headroom — per-target cost varies about 1.5× with pair count.

### Comparison 1 — old vs new pipeline — **DONE, PASSED** (job `38916341`, 2026-08-13)

**The refactor is validated.** The old pipeline's 2026-05-14 output was still on disk
(`results/day0_grna20_no_shuffle/power_analysis/power_analysis_results_es_0.15.tsv`, 34,886 pairs),
so the old environment never had to be rebuilt — which retires the CentOS 7 / pinned-`linux-64`-solve
risk described below. `old_vs_as_is` is the gating comparison, because the old pipeline inherits the
same real-data cache `as_is` does; comparing it against `null_fit` would fold in a deliberate change
of null model.

| Check | Result | |
|---|---|---|
| 1. Pair sets | identical, 34,886 shared, 0 either side | PASS |
| 2. Perturbed cells per pair | **0 of 34,886 differ, max abs difference 0** | PASS |
| 3. Mean power | 0.59293 → 0.59326, **+0.00032 = 1.37 σ** | PASS |
| 4. Paired sign test | 13,823 up / 13,642 down / 7,421 tied, **p = 0.28** | PASS |
| 5. Per-pair agreement | **Pearson r = 0.99299**, mean abs 0.0304, median 0.02, max 0.26 | PASS |
| 6. Old estimate inside new Wilson interval | 81.2 % coverage | informational |

**Check 2 is the strongest of these and is easy to skim past.** Perturbed cells per pair is
RNG-independent, so a zero difference across all 34,886 pairs proves the two implementations select
identical cells for identical pairs — the entire upstream path is bit-identical, and only the
simulation's RNG stream differs. Checks 3–5 are then about Monte-Carlo agreement alone.

Check 5's mean absolute difference of 0.0304 is not a discrepancy: at 100 replicates a single power
estimate carries a standard error up to 0.05, and both sides carry one, so per-pair scatter of this
size is what agreement looks like. What would have been alarming is a *directional* offset, and
check 4 finds none — 13,823 against 13,642 is as close to a coin flip as this many pairs allow.

**The `old_vs_null_fit` comparison "fails" check 3, and that is the correct behaviour.** Mean power
differs by 0.0066 = 28.1 σ. That is not a refactor problem — it is the deliberate null-model change,
re-measured independently: `null_fit − as_is` was +0.0063 measured directly, `null_fit − old` is
+0.0066, and `as_is − old` is +0.0003. The three are additive and consistent, which is a useful
triangulation. Check 3 tests "is this the same estimator?" against the Monte-Carlo SE of a mean over
34,886 pairs (0.00024), so any real systematic shift blows through it however small. It is reported
and does not gate.

Reports: `results/refactor/comparison/old_vs_as_is/` and `.../old_vs_null_fit/`.

---

**Why it was outstanding:** the upstream statistics are proven bit-identical (size factors, normalised
means, raw means: max absolute difference exactly 0), and dropping `@grna_matrix` was proven to leave
discovery results unchanged. What has *not* been checked is the simulation itself, because seeding
moved inside the replicate loop, which changes the RNG stream. So the comparison has to be
distributional rather than exact.

The old code needs the old environment, which is why the cluster is the right place: the old
`envs/power_analysis.yml` is a fully-pinned `linux-64` solve. It cannot be solved on macOS, but
should solve on a Linux compute node.

```sh
# 1. check out the pre-refactor tree alongside the current one.
#    b3cf7ab, not c7448a7: it is the tip of the old `low-mem` branch and includes the
#    "Improve memory usage in power simulations" commit, so it is the version actually
#    being run in practice.
#
#    Sherlock's default git is 1.8.3.1, which has no `worktree` subcommand -- load a
#    modern one first. Do NOT `git checkout` the old commit in the working tree instead:
#    the sbatch scripts run `src/` from the repository root, so swapping branches under a
#    running array silently changes the code mid-run.
ml system git/2.45.1
git worktree add ../egpa-old b3cf7ab
cd ../egpa-old

# 2. the old pipeline expects four inputs under results/<sample>/ with these exact names
mkdir -p results/sample1
cp /path/to/sceptre_object.rds        results/sample1/final_sceptre_object.rds
cp /path/to/grna_groups_table.rds     results/sample1/
cp /path/to/gene_grna_group_pairs.rds results/sample1/

# the fourth input has to be an .rds, so convert the discovery-results CSV:
Rscript -e 'saveRDS(read.csv("/path/to/sceptre_discovery_results.csv"),
                    "results/sample1/results_run_discovery_analysis.rds")'

# 3. run the old pipeline for the power results, NOT `all`
#    (`all` stops at the concatenated per-replicate file)
snakemake --use-conda --conda-frontend conda --cores 8 \
    results/sample1/power_analysis/power_analysis_results_es_0.15.tsv
```

Set `num_replicates` in that tree's `config/config.yml` to whatever you use for the new run so the
two are comparable, and remember the old `config/config.yml` also has the dead `reps: 20` key — it is
`num_replicates` that takes effect.

Then compare against the new `power_es0.15.tsv`. Suggested acceptance criteria:

| Check | Expectation |
|---|---|
| Pairs present | identical set |
| Per-pair power difference | centred on 0; no systematic offset (paired sign test not significant) |
| Mean power | agrees within the Monte-Carlo error of the replicate count used |
| Correlation of per-pair power | high (> 0.95 at 100 replicates) |

A systematic offset in either direction is a real discrepancy and worth chasing. Scatter across
individual pairs is expected — at 100 replicates the standard error of a single estimate is up to
0.05, so differences of ±0.1 on individual pairs are normal.

Use the same `--reps` in both, and note the old pipeline cannot be seeded, so it produces a different
draw every run. Running it two or three times gives a sense of its own run-to-run spread, which is
the yardstick for judging the old-vs-new difference.

### Comparison 2 — monotonicity at full replicate count

**DONE, PASSED** on `day0`, 2026-09-03. Power must not decrease as effect size increases, and a
violation that survives at full replicate count is a bug.

| | |
|---|---:|
| Pairs | 34,886 |
| Consecutive comparisons (5 per pair) | 174,430 |
| Decreases | **2,300 (1.32 %)** |
| Largest decrease | **0.070** |
| Pairs with `n_reps` ≠ 100 | **0** |

Nothing survives as a bug. At 100 replicates the standard error of a power estimate near 0.5 is
0.05, so a *difference* of two estimates carries an SE around 0.07 — and the largest decrease
observed is exactly 0.070, one standard error. Every violation is inside Monte-Carlo noise.

The rate is about double `wtc11`'s (0.67 %, 220 of 32,870), which is expected rather than
concerning: `day0` is the less-powered dataset, so more of its pairs sit in the steep part of the
curve where adjacent effect sizes are close together and noise can invert them. Worth stating in
the paper rather than leaving for a reviewer to ask about.

The earlier evidence, kept because it is what motivated running this: at 12 replicates on 33 pairs
it held (mean power 0.356 → 0.604 → 0.838 for 0.15 / 0.25 / 0.5) with one violation of 0.08 that a
two-proportion test could not distinguish from noise.

### Comparison 3 — timing and memory at scale

Run with `--with-report`-style accounting (or just `sacct`) and record actual per-task wall time and
`MaxRSS`, so the resource requests above can be tightened. The 4 GB figure comes from a laptop run of
three targets; confirm it holds for splits containing the largest targets (up to 36 pairs).

---

## Decisions already made

Recorded so they are not relitigated:

| Decision | Rationale |
|---|---|
| Power-analysis only; no differential-expression step | the three orphan diffex scripts were unwired duplicates |
| Only the sceptre object as input | the other three files duplicate data inside it |
| No Bioconductor | `SingleCellExperiment` was used purely as a container |
| sceptre pinned to a commit, installed by `pixi run setup` | not on conda-forge/bioconda, and the pipeline reads unexported slots |
| `n_control_cells` off by default | measured 21–60 % power loss |
| Threshold derived from `@discovery_result`, not `alpha` | reflects the correction actually applied |
| Wilson intervals reported per pair | the normal approximation gives `[0, 0]` for zero successes |
| Prefer more `n_splits` over more `reps_per_chunk` | chunking replicates re-pays the per-target setup |
| local + slurm profiles only | no cloud executor needed |

## The analysis question this is all for

False-negative triage: for a pair CRISPRi did not call, is the link absent, or could the experiment
not have seen one? That framing, what it implies for which columns matter, and the per-element
caveats moved to [broadinstitute/WattEG-paper](https://github.com/broadinstitute/WattEG-paper) with
the paper analyses. What it means for this repo is unchanged and documented in
[Output]({{ site.baseurl }}{% link output.md %}#interpreting-negatives): threshold `power_ci_low`
rather than `power`, and treat `min_detectable_effect_size` as the deliverable rather than the
per-effect-size power columns.

## Open questions

- **Two-stage replicate allocation**: worth wiring in? Costed properly it is ~1,768 CPU-h per effect
  size against 3,430 for a uniform 400 and 858 for a uniform 100 — a 1.9× saving over uniform 400,
  not the 3.7× an earlier version of this document claimed. That claim came from counting replicate
  equivalents, which ignores the per-target term; see
  [Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}). It only pays
  off if per-pair cutoff decisions are being made — which, per the section above, they are.
- **Is the 5 % effect size worth 100 replicates?** Nearly every pair is underpowered there, so the
  replicates buy a precise estimate of a number that is 0. Sweeping 5 % at 30 reps and spending
  the difference on 20–25 % would carry more information. Untested.
- **`--target-overhead` in `split_pairs.R`: closed, and the answer is "leave it".** It defaults to 0,
  so every split is balanced on pair count alone. The `as_is` array fixed the right value at the time
  (5.076s per target against 0.5923s per pair ⇒ ≈ 8.6 pair-equivalents), but setting it bought almost
  nothing at `N_SPLITS=1000`, and the reason this document previously gave for setting it was wrong.
  **The gRNA reuse has since made it moot**: at 1.140s per target against 0.5561s per pair the right
  value is **≈ 2.0 pair-equivalents**, on splits that already hold 34–36 pairs. There is nothing left
  to correct for.

  Measured on the real splits, and on job `38849611`'s own `sacct` record:

  | | Predicted tail task | Predicted floor | max/mean |
  |---|---:|---:|---:|
  | `--target-overhead 0` (what ran) | 68.4 min | 43.0 min | 1.139× |
  | `--target-overhead 8.6` | **67.4 min** | 58.0 min | 1.123× |

  Simulating the same LPT bin packing `split_pairs.R` uses reproduces the on-disk splits exactly at
  overhead 0, so the comparison is like for like. The overhead raises the *floor* a lot (43 → 58 min)
  and moves the *tail* by **one minute** — and the tail is what determines when the array finishes.
  At 1,000 splits each task holds only **1–4 targets**, so there is little imbalance left to
  reclaim.

  **This refutes the previous revision of this point**, which blamed the observed 35–110 minute
  spread on pair-only balancing. It is not that: the cost model predicts a 43–68 minute range, and
  regressing actual task time on predicted gives **r = 0.373, r² = 0.139** — split composition
  explains 14 % of the variance. The model is unbiased on the mean (36.02s predicted against 35.97s
  actual per replicate), so it is not simply wrong. What explains the spread is **node
  heterogeneity**: across the 229 nodes the array landed on, the slowest node's mean is **2.73×** the
  fastest (55.8s against 20.5s per replicate). There was no preemption to blame (999 COMPLETED,
  1 FAILED). Set `--target-overhead 8.6` if `N_SPLITS` ever drops far enough that splits hold ten or
  more targets; do not re-split an array for it.
- **High MOI is a hard requirement, and now a measured one. Low MOI does not work.** Found by
  running the synthetic generator's new `--moi` / `--resampling-mechanism` variants end to end
  (`13_make_test_data.sbatch`, 2026-09-03):

  | Variant | `@response_precomputations` | End to end |
  |---|---:|---|
  | high MOI / crt | 30 | passes |
  | high MOI / permutations | 30 | **passes** |
  | low MOI / permutations | **0** | **fails** |

  So **permutations is supported** — it was the untested axis and it works, which also means an
  odm-backed screen (forced onto permutations) is fine on this axis. The failure is MOI alone:

  ```
  Error: @response_precomputations is empty, so no dispersions can be derived.
  ```

  Under low MOI the control group is the NT cells rather than the complement, and sceptre never
  populates the per-gene precomputation cache that `prepare_sim_input.R` reads per-gene
  negative-binomial dispersions (`$theta`) out of. This is not a missing branch in our code, it is
  a different source of dispersions: supporting low MOI means fitting them ourselves, or calling
  `perform_response_precomputation()` directly per gene, and then validating the result against the
  cached values on a high-MOI object. Deliberately deferred, not attempted.

  Note the error a low-MOI user currently gets names the empty slot rather than the MOI, so it
  reads like a corrupt object rather than an unsupported configuration.

- **`grna_integration_strategy` must be `union`**, rejected explicitly at `lib/sceptre_io.R:55`
  because the target-level perturbation matrix and the discovery pairs are both keyed by target.
  `singleton` and `bonferroni` fail fast with that message rather than silently producing
  target-keyed output from gRNA-keyed pairs. Untested against low MOI, and expected to be
  independent of it: the precomputation cache above is per-gene, not per-gRNA-grouping.

- **Which test a screen ran is now read, not assumed.** `prepare_sim_input.R` reads
  `@run_permutations` and `@low_moi` off the object, logs them, and writes `analysis_mode.tsv`
  beside `discovery_threshold.txt`; `run_power_simulation.R` logs the same on every task. This
  matters because the pipeline re-runs the screen's own test on simulated counts, so the test
  decides what the power numbers mean — and it is set upstream by `set_analysis_parameters()` and
  inherited silently. There is no `resampling_mechanism` slot: sceptre collapses that string to the
  boolean `@run_permutations`.

  **Three routes to permutations, none of them visible at the call site:** an explicit request;
  `low_moi = TRUE` with `resampling_mechanism` left at `"default"`, which resolves to
  `permutations` for low MOI and `crt` for high MOI; or an odm-backed response matrix, which
  `check_set_analysis_parameters()` (sceptre `R/check_functions.R`, check 10) forbids from using
  CRT at all. Our datasets are `moi = "high"` and in-memory, which is the only reason this pipeline
  has been on the CRT path without ever choosing it.

  Consequences to carry into the paper: an odm-ingested screen is permutations **end to end**,
  including the `@discovery_result` the threshold is derived from, so its power is correct for that
  screen but was produced by a different test from a CRT dataset's — the §6 cross-dataset transfer
  claim has to say so rather than have a reviewer find it. And the cost model
  (`1.140s + 0.5561s × pairs`) was fitted on the CRT path *with* the gRNA precomputation patch;
  neither applies under permutations, where `run_power_simulation.R` now skips building a cache the
  permutation branch would discard.

- **ODM / out-of-core support** is implemented. `get_response_matrix()` in `lib/sceptre_io.R`
  reconnects an odm-backed response matrix (via a new `--response-odm` on `prepare_sim_input.R`,
  following the same reconnect-plus-`@integer_id`-check pattern as sceptre's own
  `read_ondisc_backed_sceptre_object()`) and materializes it into an ordinary in-memory sparse
  matrix, once, before `compute_expression_stats()` ever runs — so nothing downstream changed.
  This works because the poscounts size factors need a per-cell median over *all* genes (a
  column-wise reduction the `odm`'s row-only access can't do directly), but reading every row into
  memory turns out to be small even at whole-transcriptome scale: measured at ~1.8 GB sparse for
  38,606 genes x 147,856 cells (2.6-2.7% density) on a real DC-TAP-seq object, end to end in under
  3 minutes. `ODM_MATERIALIZATION_LIMIT_GB` guards the case this used to be written for — a dataset
  large enough that materializing really would be a mistake.
