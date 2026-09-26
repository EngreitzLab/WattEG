---
title: Plan - pysceptre backend
nav_order: 8
---

# Plan: move the high-MOI power analysis onto pysceptre

Replace the R/sceptre engine inside the power simulation with
[pysceptre](https://github.com/broadinstitute/pysceptre), so that the whole pipeline is one Python
package, `pixi.toml` drops R entirely, and the per-simulation cost falls by the margin pysceptre
already demonstrates on discovery analysis.

This document is a plan, not a record of work done. Nothing below has been implemented.

## 0. Branches, and why

| Repo | Branch | Rule |
|---|---|---|
| `WattEG` | **`feat/pysceptre-backend`** (created for this work) | `main` stays the R implementation until this branch merges |
| `pysceptre` | **`v0.2.0`** (tag) | The §5 export work was done on `feature/watteg-simulation-support`, squash-merged into `0.1.1rc` as **`d96d48f`** alongside the analytical-power work, and has since reached `main` and been tagged. `pixi.toml` pins the **tag**, not a SHA — the references to `d96d48f` below are the history, not the pin |
| `WattEG` | **`r-implementation`** | the R pipeline, kept **maintained rather than frozen**. Branched from `main` at `7c07825` and since carrying the same three fixes this branch does (`d767af3`, `2b76284`, `b346296`) |
| `WattEG` | `legacy` | untouched (the Snakemake implementation that preceded both) |

The R path is not deleted when the Python path lands. It is the reference the Python path is
measured against, and it is what the paper describes; retiring it is a separate decision, taken
after §8 reports.

**It is also not frozen.** The first instinct was to leave `r-implementation` untouched so it would
keep reproducing the paper's numbers — but nothing has been published, so there is no obligation to
a set of numbers, and what that would have preserved is a bug. Two of the three fixes below change
results, and both were applied to R as well as to Python. The consequence is stated in §11b: the R
reference Stage B compares against has to be **regenerated**, not read off the existing sweep.

## 1. What actually moves

The pipeline is five steps plus two support steps. Only one of them is expensive, and only one of
them needs R.

| Today | Fate under the Python backend |
|---|---|
| `src/prepare_sim_input.R` | **ported to Python**, reading a pysceptre export instead of the `.rds` (§4). Everything it computes — poscounts size factors, normalised means, dispersions, the discovery threshold — is computable from the exported matrix |
| `src/split_pairs.R` | ported, trivially |
| `src/fit_null_models.R` | **deleted** (§3.1) |
| `src/merge_null_models.R` | **deleted** (§3.1) |
| `src/run_power_simulation.R` + `lib/simulate.R` + `lib/pert_input.R` | **ported to Python**; the NB draw, guide-to-guide variability, centering and seeding are WattEG's method and stay in WattEG (§6) |
| `src/consolidate_replicates.R` | ported (pyarrow) |
| `src/compute_power.R` + `lib/stats.R` | ported (Wilson interval) |
| `src/summarize_power.R`, `src/fit_power_curve.R` | ported |
| `patches/`, `lib/apply_patch.R`, `src/check_sceptre_api.R`, `src/install_sceptre.R`, `src/install_ondisc.R`, `src/audit_dependencies.R` | **deleted.** All six exist only because the pipeline reaches into unexported sceptre S4 slots and patches its CRT path |
| `lib/sceptre_io.R` | **deleted from the pipeline**; its odm-materialisation logic moves to the one-off export (§4) |
| `src/make_test_data.R` | ported, or replaced by a synthetic fixture generated in Python |

Everything in `workflow/slurm_executor/` and `workflow/compare_*.R` is scaffolding around the R
steps and follows them.

## 2. The core mapping, and the reason it is fast

Today the unit of work is **one `run_discovery_analysis()` call per (target, replicate)**: on
`day0_grna20` that is 3,026 targets x 100 simulations = 302,600 R calls per effect size, each one
carrying the full 567,690-cell bookkeeping for a median of 9 gene pairs. The measured cost model is
`1.140s + 0.5561s x pairs` per (target, simulation), **48 % of it in the per-target term**.

pysceptre's `run_discovery_analysis` takes plain arrays — `response_matrix` (n_genes x n_cells),
`covariate_matrix`, `grna_target_cells` (dict target -> 0-based cell indices), and a `pairs` frame.
That signature lets one call cover **a whole replicate chunk of one target**:

```
for target in split:
    rows  = [f"{gene}@{rep}" for rep in reps for gene in target_genes]   # pseudo-genes
    counts = vstack(simulate(target, gene, rep) for rep in reps)          # (genes*reps, n_cells)
    pairs  = DataFrame(response_id=rows, grna_target=[f"{target}@{rep}" ...])  # pseudo-targets
    run_discovery_analysis(counts, rows, covariates, {f"{target}@{rep}": cells}, pairs, ...)
```

Five consequences, each verified against the pysceptre source rather than assumed:

**2.1 Per-gene null fits come for free and are per-replicate.** `fit_all_genes` fits every row of
the response matrix independently (`_GENE_BATCH_WIDTH = 1`, `pipeline/discovery.py:89` — pinned at 1
deliberately so a fit cannot depend on its batch companions). Stacking replicates as rows therefore
fits each `(gene, replicate)` null model **on that replicate's own simulated counts**. That is the
`cleared` configuration in [Status]({{ site.baseurl }}{% link status.md %}) — the faithful reference
that `null_fit` was built to approximate — obtained at no extra cost. §3.1.

**2.2 Replicates must not share a CRT draw.** pysceptre seeds each target's resampling stream from
the *target's name* (`target_seed_sequence`), which is what makes its results independent of chunking
and target order. One shared key for all replicates would hand every replicate the **same** synthetic
treated-cell index sets, correlating the replicates: R redraws them per call, and the Wilson interval
in `compute_power` assumes independent Bernoulli trials. Hence the `target@rep` pseudo-target keys
above — independent streams by construction, and chunk-layout invariance inherited for free.

The key must carry the **effect size** as well as the replicate — `target@es@rep` — or a single
`seed` would hand every effect size in a sweep the same synthetic index sets. R's `derive_seed`
includes `effect_size` for exactly this reason. (The alternative, deriving the per-call `seed` from
the effect size, works too but makes the invariance harder to see; prefer the key.)

The cost is that the target's binomial GLM (perturbation status on covariates) is refit once per
replicate although the input is identical — precisely the redundancy the sceptre patch in `patches/`
removes on the R side, worth 999 -> 634 CPU-h there. pysceptre batches those fits across a target
chunk, so the redundancy is far cheaper than R's. **It is accepted, not fixed** — removing it would
mean changing the pysceptre package, and §5.3 says why that bar is not met.

**2.3 Memory is a storage question, not a fitting question.** Both row-reading paths — `_gene_fit_job`
(`discovery.py:257`, which at `_GENE_BATCH_WIDTH = 1` fills a one-row float64 `Y`) and `_gene_job`
(`discovery.py:804`) — go through `_get_row` (`discovery.py:190`), which densifies **one row at a
time** and casts to float there. So the stacked matrix can be stored as `int16` or sparse without
any change to pysceptre. For the median target: 9 genes x 100 replicates x
567,690 cells is 1.0 GB as dense `int16`, against 8.2 GB as float64. The max target (36 pairs) is
4.1 GB, which is why `reps_per_chunk` stays a parameter and becomes the memory knob it always
implicitly was. Overflow guard: counts above `int16` range must promote, not wrap.

**2.4 Call the inner entry point, not the public one.** `pipeline.api.run_discovery_analysis`
sizes `B2`/`B3` from `len(pairs)` (R's `run_qc` rule). With pseudo-pairs the pair count is inflated
by the replicate count, which under `resampling_approximation = "no_approximation"` would inflate
the `B3` budget ~100x. Call `pipeline.discovery.run_discovery_ntcells_complement` directly and pass
`B1`/`B2`/`B3` from the export's `metadata.json` — which is what R does, since the template carries
the real analysis's `@B1/@B2/@B3` slots and the simulation never re-derives them.

**2.5 Where the time will actually go — and why "30 s" is not the answer.** A full
permutation-mode pysceptre discovery analysis on a real dataset runs in ~30 s. That number does not
transfer to this pipeline, and reading it as though it does would set the wrong expectations and
optimise the wrong thing. One discovery analysis fits **237 gene nulls** and tests ~35,000 pairs
once. One effect size of a 100-replicate power sweep fits **34,886 x 100 ≈ 3.5 million** per-(gene,
replicate) nulls, draws 3.5 million negative-binomial count vectors of 567,690 values, and runs 3.5
million resampling tests. It is roughly four orders of magnitude more gene-fitting work than the
analysis whose runtime is 30 s.

So the useful reading of that 30 s is a **per-unit rate**, and the three terms it decomposes into
are what the phase-3 benchmark has to report separately:

1. drawing the counts (`rnbinom` was 7 % of R's time; in Python it will be a larger share, because
   everything around it got faster);
2. the per-(gene, replicate) Poisson IRLS null fit — 59 % of R's time was `glm.fit`, and §2.1 makes
   this term unavoidable rather than hoistable;
3. the per-pair resampling test against that target's CRT draws.

Whichever dominates is where any later optimisation goes. The redundant binomial fits of §2.2 are a
fourth term and, on this reasoning, the smallest of the four — which is the quantitative case for
§5.3 staying unbuilt.

## 3. Statistical decisions this forces, stated up front

### 3.1 The null model becomes the faithful refit

`fit_null_models.R` / `merge_null_models.R` and the whole `--null-precomputations` bundle exist for
one reason: R refitting the per-gene null inside every call costs **4.3x**, so the fits were hoisted
into their own Nextflow processes, fitted once per `(gene, simulation)` on a null simulation, and
injected. Status records the verdict: `null_fit` is **exactly equivalent** to the faithful
`cleared` refit (0 flips in 265 calls), and `as_is` — the inherited real-data cache — understates
power by +0.0063 mean.

Under §2.1 the faithful refit is what happens anyway. So: delete both processes, delete the seed-
matching guard between bundle and simulation, delete the `@response_precomputations` trap in
`slim_sceptre_object`. **Expect the Python results to match `power_null_fit/`, not `power_as_is/`**
— that is the comparison baseline in §8, and picking the wrong one would manufacture a +0.006
discrepancy out of a known, already-settled difference.

### 3.2 Dispersions still come from real data

`row_data$dispersion` is `1/theta` from sceptre's `@response_precomputations`, fitted on the **real**
counts — it sets the noise the simulation exists to reproduce, and it must not become a property of
the simulated counts. pysceptre has its own validated theta estimator (`glm/nb_theta.py`), so the
Python `prepare_sim_input` computes theta from the real matrix rather than reading a cache.
**Measured, and it settles the question.** Over day0's 237 genes and 567,690 cells, pysceptre's
fitted theta reproduces sceptre's cached theta to a **median relative difference of 2.8e-12 and a
worst case of 1.2e-9** — about ten significant digits — with 236 of 237 genes inside 1e-9. Since a
dispersion matters only through the variance of the counts drawn from it, that is nothing against
the 5–50 % effect sizes the sweep tests. One gene's theta MLE fell back to method of moments in
both implementations alike, and it is the gene the two differ on most.
`workflow/compare_dispersion.py` is the gate, at a tolerance of 1e-6 — generous against what the
port achieves and still tight enough to catch a wrong model.
`build_dispersion_vector`'s hard error on missing/non-finite dispersions carries over, and is
joined by one **deliberate deviation**: a theta clamped to the estimator's bounds is refused rather
than simulated from. sceptre clamps to exactly the same `[0.01, 1000]` that pysceptre does
(`perform_response_precomputation`: `max(min(theta, 1000), 0.01)`) and carries on, so R would
proceed where this stops. For an *analysis* that is reasonable; for a *simulation* a clamped theta
is not an estimate of anything, and drawing counts from it would state a noise level the data never
supported. Day0 has no clamped gene, so nothing is refused today.

### 3.2b The simulated mean is biased low, and the fix is exact

**sceptre is not involved, and this is the first thing to establish.** sceptre computes no size
factor, no geometric mean and no normalised mean — searching all 183 of its functions for
`size_factor|geomean|normaliz|offset` returns nothing — and its per-gene model is
`glm.fit(y = counts, x = covariate_matrix, family = poisson())` with no offset. Library size enters
it as ordinary covariates, `log(response_n_umis)` and `log(response_n_nonzero)`, which the GLM fits
coefficients for. The size factors and `row_data$mean` are **WattEG's simulation machinery alone**,
inherited from the original DC_TAP_Paper power simulation; they exist to generate synthetic counts
and sceptre never sees them. DESeq2 is not doing it wrong either — mean-of-ratios is exactly what
its `baseMean` is, and as a summary statistic it is fine. What follows is about one *composition*,
which is WattEG's own.

**Found while checking a report that the simulation runs 16 % low.** It does not, but it does run
low. `row_data$mean` is 16 % below the raw mean on day0 — that much is true and expected, because
it is a *normalised* mean — but the simulation never uses it alone: `draw_counts` forms
`mu[i,j] = mean_i x size_factor_j x effect_size`, which puts the size factor back. With
`mean(sf) = 1.1389` on day0 the simulated per-gene expected raw mean lands at **0.959 of the real
one** at the median (range 0.89–1.04; 84 of 237 genes more than 5 % low, one more than 10 %).

The residual is a real bias with a clean cause. `mean_i` is a **mean of ratios**,
`(1/n) sum_j counts[i,j]/sf_j`, and

    raw_mean = E[x.sf] = E[x].E[sf] + Cov(x, sf),   x = counts/sf

so multiplying by `mean(sf)` drops the covariance term. It is positive here — cells with larger
size factors still carry slightly more normalised counts, i.e. the normalisation under-corrects —
so the simulation draws low. **Whether that makes simulated power low is not established**, and
the obvious inference is unsafe: the same baseline also understates the count *variance* by 13.5 %
(§3.2c), and less variance inflates power where less expression deflates it. Two errors, opposite
signs, neither measured against the other. Stage B under both baselines is the experiment.

**The alternative is a ratio of sums, and it is exact rather than better:**

    mean_i = rowSums(counts)_i / sum_j sf_j

Then `sum_j mean_i.sf_j = sum_j counts[i,j]` identically, for every gene, by construction. On day0
it raises each gene's mean by ~4.3 %.

**Not changed, pending a decision.** It moves every number the paper reports, and the port's job
is to reproduce R first — that is what the phase-2 gate is. It is a better candidate for actually
changing than the §5.2 question, though: that one is a judgement about which cells belong in an
estimator, this one is an estimator that provably fails to reproduce the counts it was derived
from, with an exact alternative available. If it is taken, it is a one-line change in
`watteg/expression.py` and a corresponding one in `src/prepare_sim_input.R`, and both sweeps
have to be re-run.

### 3.2c The same question on the analytical side, and what it suggests

pysceptre's analytical power estimator hit this first, and its
`analytical_power/inputs.py` already says so: `baseline_expression_stats` "comes from a
normalisation scheme sceptre does not use, and on day0 it sits about 16 % below the mean sceptre's
own model implies", with `baseline_expression_stats_from_fits` recommended instead.

**PerturbPlan is not doing anything wrong, and its source says so twice.** `compute_power_posthoc`
does not compute `expression_mean`; it takes `baseline_expression_stats` as an argument, and the
documentation defines it only as "a data frame ... with columns `response_id`, `expression_mean`,
and `expression_size`" — **it never states the scale.** The formula does, in two independent
places: `compute_distribution_teststat` uses `var_nb(mean, size) = mean + mean^2/size`, the
variance of the negative binomial the *observed counts* follow; and `compute_QC` takes
`P(count == 0)` from that same NB and feeds it to
`pbinom(n_nonzero_thresh - 1, num_cells, 1 - P0)` — the chance that enough cells have a nonzero
**raw count** to clear pairwise QC. The second settles it: only a distribution over actual counts
has a zero probability to ask about. So `expression_mean` must be E[observed count per cell], and
the defect is entirely in the input.

**And the wrong scale costs more there than the flat 16 % suggests, because it is counted twice.**
Measured on day0, the normalised mean overstates `P(count == 0)` by 0.032 at the median and 0.073
at most (0.304 against 0.247). That inflates `QC_prob`, and power is multiplied by
`1 - QC_prob`: at the median target's 396 treated cells, **28 of 237 genes carry an inflated
`QC_prob`, the worst by 0.20** — a fifth of that gene's power disappearing into a QC term, on top
of the separate understatement through the test statistic.

What is wrong, then, is feeding it the poscounts normalised mean, and the two tools differ only in
how far that input travels:

| | vs the mean sceptre's model implies |
|---|---|
| analytical formula — uses `expression_mean` **as-is**, no per-cell factor | **0.84, 16 % low** |
| WattEG's simulation — `mean_i x sf_j`, the factor multiplied back | **0.959, 4.1 % low** |

"The mean sceptre's model implies" is measurable and is exactly the observed raw mean:
`mean(exp(X.beta))` reproduces it to 8e-9, because a Poisson GLM with an intercept satisfies
`sum(fitted) == sum(observed)`. That identity is what makes this comparison sharp rather than a
matter of taste.

**The deeper point, which the level difference hides.** sceptre's model is
`E[count_ij] = exp(X_j . beta_i)`, a full covariate-dependent mean. WattEG simulates
`E[count_ij] = mean_i x sf_j`, one scalar per cell. Those differ in *shape across cells*, not only
in scale, so the ratio-of-sums fix in §3.2b patches the average and leaves the structure wrong.

**So the candidate fix is bigger than §3.2b and subsumes it: simulate from `exp(X_j . beta_i)`
directly**, which is the same move the analytical side already made. It reproduces the observed
mean exactly rather than approximately, reproduces the per-cell variation the test's own model
assumes, and deletes the poscounts machinery entirely — with it go both §5.2's QC-cells question
and §3.2b's estimator question, which exist only because size factors do. `fit_dispersions`
already computes `fitted_coefs` and currently discards them.

**The argument against, stated because it is real.** Simulating from the fitted model and then
testing with that same model makes the test perfectly specified by construction, which may
overstate power slightly; the present scheme is misspecified in the other direction. Neither is
neutral. "Matches the model the real data was fit with" is the more defensible starting point, but
this is a decision about what the power analysis *means*, not a bug fix, and it is the user's.

**Taken, and implemented.** `exp(X . beta)` is the default baseline in the Python port
(`watteg/baseline.py`), `size_factor` survives as a validation fixture so Stage B can compare
like-for-like against sweeps produced with it, and `sim_input` format 2 carries `fitted_coefs`.
All four phase-2 gates stay green, because the change is additive: everything the port already
reproduced is untouched.

**What this does *not* settle.** The effect on simulated power is unmeasured and not obviously
signed, for the reason in §3.2b. Every number here is day0; the mechanism is general but the
magnitudes are not. And the existing sweeps were run on the old baseline, so re-running them is a
separate, all-or-nothing decision: mixing the two scales within one analysis would be worse than
either alone. What that decision costs, concretely: every `power_summary.tsv` moves — `day0` at six
effect sizes, `moi5` cis at six, and the 742,525-pair `moi5` trans sweep at one.

**That decision has since been taken, for a different reason.** The R implementation is *not*
unchanged: it carries this baseline (`d767af3`) and the centring fix (`2b76284`), so the sweeps have
to be regenerated regardless of how this question alone would have been decided. See §11b.

Both per-gene means survive the change on purpose. `sim_input.h5` carries `fitted_coefs` *and*
`mean`, so an existing output can be audited against the scale that produced it without
re-deriving anything, and `--expression-model size_factor` reproduces it exactly.

### 3.3 Seeding contract is preserved

Today: `set.seed(derive_seed(seed, target, rep, effect_size))` before each replicate, and a separate
`rep = 0` key for per-target setup, which is what makes results invariant to split layout and
replicate chunking. Python equivalent: `np.random.SeedSequence` spawned from the same four-part key
(hashed, not R's integer arithmetic — the streams differ from R's either way). The invariance test
(1x4 against 2x2 chunking, and two different split layouts) ports directly and becomes a unit test
rather than an sbatch script.

### 3.4 The `log_2_fold_change < 0` condition

pysceptre returns `fold_change`, not `log_2_fold_change`; `compute_power`'s one-sided condition
becomes `fold_change < 1`, which is the same predicate. `pct_change_es` and its CI are a bonus the R
path never had.

## 4. The input boundary — the one real design decision

`prepare_sim_input.R` reads a `sceptre_object` `.rds`. Nothing in Python reads that. Two options:

| | Keep `prepare_sim_input.R` in R | **Recommended: `.h5mu` in, R export run once, outside the pipeline** |
|---|---|---|
| Env | still needs `r-base` + pinned sceptre + ondisc + the patch machinery | pure Python; `pixi.toml` drops R |
| User cost | none | one `Rscript export_sceptre_dataset.R` per dataset, in a container we provide |
| Pipeline surface | unchanged | samplesheet column becomes `dataset` (`.h5mu`) instead of `sceptre_object` |

Take the second. The stated goal is one Python package, and a single R step in the pipeline keeps the
entire R toolchain — pin, patch, `check-api`, two source installs — alive to serve it. The export is
a per-dataset, one-off, already-written script (`pysceptre/scripts/export_sceptre_dataset.R` +
`make_h5mu.py`) that handles odm-backed and in-memory matrices alike, so `lib/sceptre_io.R`'s
materialisation logic has an owner.

**Keep a `--sceptre-object` path on the R side for one release** as an escape hatch and so §8 can run
both engines from the same object.

What the export already carries and WattEG needs: the count matrix (`--all-genes`, required —
poscounts size factors are a whole-gene reduction), the covariate matrix restricted to
`cells_in_use`, QC-passing pairs, `discovery_result` (from which the nominal threshold is derived
exactly as `discovery_threshold()` does today), and in `metadata.json` the `side_code`,
`resampling_approximation`, `run_permutations`, `control_group_complement`, `B1/B2/B3`,
`multiple_testing_alpha` and both `n_nonzero_*` thresholds. That covers `analysis_mode` and every
pysceptre argument.

What it does **not** carry is enough cells — see §5.2, which is a blocking gap, not a detail.

## 5. What pysceptre must change — the export only

**Status: done, in `../pysceptre` on `0.1.1rc`** — squash-merged as **`d96d48f`** ("ADD analytical
per-pair power, per-gRNA exports, and a minimal container"), which carries the export work together
with changes of its own. The three commits it squashes (`7e29efa`, `5117651`, `37f9d81`) survive
only on `feature/watteg-simulation-support`, so `d96d48f` is the reference that will keep
resolving. Both gaps are closed and verified against the real day0 object; §5.1 and §5.2 below are
kept as the record of what they were and of what closing them turned up. §5.3 remains unbuilt, as
planned.

**The constraint that shapes this whole section: the pysceptre *package* does not change.** Both
gaps below are in `scripts/`, which pysceptre's own `CLAUDE.md` marks as *not shipped in the wheel* —
they are dataset-export tooling, not the statistical engine. Nothing in `src/pysceptre/` is touched,
so nothing this plan does can move a pysceptre result, and the validation burden stays on WattEG
where it belongs.

Verified gaps, not speculation.

**5.1 Individual targeting-gRNA assignments — blocking.** The export writes
`grna_assignments$grna_group_idxs`, which is the **union of each target's gRNAs**, plus individual
*non-targeting* gRNAs (`scripts/sceptre_export_lib.R:112-139`). WattEG's guide-to-guide variability
(`create_guide_pert_status`, `create_effect_size_matrix`, `guide_sd = 0.13`) needs **per-gRNA**
membership for targeting guides, and the `grna_id -> grna_target` map. Add both to the export: a
third block of assignment rows with `unit_kind = "targeting_grna"`, and `grna_target_data_frame`
written out whole.

**What closing it turned up, and what it means for §6.** The gRNA -> target map is
**many-to-many**: 1,673 of day0's 43,736 guides sit inside two or three *overlapping* candidate
elements and so belong to two or three targets (45,463 design rows against 43,736 distinct ids).
R handles this without comment — `grna_map$grna_id[grna_map$grna_target == target]` selects by
target, so a shared guide is simply returned for both. Anything that collapses the map to one
target per guide — a `dict`, a `match()`, the per-unit `var` annotation — drops those guides from
every target but one, which on day0 would leave **216 of 3,071 targets simulating with an
incomplete guide set**, silently.

So `perturbation.py` must read guides-per-target from `grna_target_data_frame`, **never** from the
gRNA assay's `var` annotation, which has one row per unit and therefore records `"<multiple>"` for
a shared guide. The export asserts the union round-trip exhaustively over every target at write
time, which is what caught this; all 3,071 day0 targets reproduce exactly.

While the export is open: it writes only `response_id` and `grna_target` for the QC-passing pairs
(`qc_passing_pairs`), but the R simulation output carries `n_nonzero_trt`, `n_nonzero_cntrl` and
`pass_qc` from `@discovery_pairs_with_info` — real-data diagnostics, constant across replicates, and
the first thing anyone looks at when a pair's power is surprising. Write that frame whole, or drop
those three columns from the byte-compatibility promise in §6. Prefer writing it.

All additive; no existing reader changes.

**5.2 All cells, not just `cells_in_use` — blocking.** `prepare_sim_input.R:263-270` calls
`compute_expression_stats()` on `get_response_matrix(so)`, the **whole** matrix: 586,309 columns on
`day0_grna20`, against 567,690 in `cells_in_use`. Poscounts size factors are a per-cell reduction
over a per-gene geometric mean taken across *all* cells, so computing them on the QC-passing subset
gives different size factors and different normalised means — which is the input the simulation
draws from. The export writes `cells_in_use` only (`sceptre_export_lib.R:37,61,85`). Add an
`--all-cells` mode that writes the full matrix plus a `cells_in_use` index vector.

Until that lands, phase 2's column-by-column gate **will** fail, and it would be easy to
misattribute the failure to theta (§3.2). It also constrains Stage A: matrices dumped from Python
have to be indexed the way `template@cells_in_use` expects before R can test them.

**Closed by `--all-cells`, and this is what it buys.** Measured across the 18,619 cells QC removes
on day0, computed both ways:

| Quantity | Median shift | Max |
|---|---:|---:|
| Raw gene mean | 2.2 % | 6.2 % |
| poscounts size factor (in-use cells) | 0.46 % | 1.3 % |
| Size-factor-normalised gene mean | 0.36 % | 3.0 % |

The file keeps **one cell space** — under `--all-cells` the matrix columns, covariate rows and
every gRNA unit are absolute positions together — and `load_export` subsets back to `cells_in_use`
by default, so an analysis reads either file identically and only the simulation passes
`all_cells=True`. Verified on day0: the default read of the `--all-cells` export is identical to
the plain one across all 92,622,239 nonzeros, the covariates, all 46,789 units and the pair table.

One thing the round-trip assertion forced into the open: **gRNA membership is post-QC in both
spaces.** `@grna_assignments` is built after QC while `@initial_grna_assignment_list` is the
pre-QC input, so a target's guides between them cover cells the target does not — 518 against 493
on day0's first target. The guides are restricted to `cells_in_use`, which keeps the union
invariant true in every file and costs nothing, since those cells have no covariates and no test
sees them. `--all-cells` therefore adds cells to the **expression side only**.

**A question for phase 2, not for the port.** Whether cells QC removed *should* enter the per-gene
geometric mean that sets the size factors is a scientific question, and the honest answer is that
R's implementation includes them because it reads the whole matrix, not because anyone chose it.
The port reproduces R first — that is what the phase-2 gate is for — and the table above is the
order-of-magnitude argument for deferring it: a 0.36 % shift in the gene mean the simulation draws
from, against effect sizes of 5–50 %. That is an estimate, not a measurement — `as_is` shifted
mean power by +0.0063 from coefficient differences far larger than this, so the direction is right
and the size is not established. Measure it once the Python path reproduces the R one.

**5.3 A shared target fit across aliased targets — considered and NOT planned.** The `target@es@rep`
keys of §2.2 make pysceptre refit each target's binomial GLM once per replicate although the input
is identical. On the R side removing that redundancy was worth 999 -> 634 CPU-h, which is why it
gets a mention at all. Here it does not: pysceptre batches those fits across a target chunk, and a
full permutation-mode discovery analysis on a real dataset runs in **~30 s**, so the engine is not
plausibly the bottleneck in a simulation whose per-replicate cost is dominated by drawing counts and
fitting per-gene nulls (§2.5).

Adopting it would mean an API change inside `src/pysceptre/` — an optional `target_fit_key` letting
several target keys share one fit while keeping their own CRT stream. That is a change to how
pysceptre works, so the bar is not "it would be faster": it is the phase-3 benchmark showing the
redundant fits are a **large** share of simulation wall clock. Absent that number, this stays
unbuilt, and the plan assumes it never gets built.

**5.4 Nothing else.** The engine is used as published. If a change to `discovery.py` turns out to be
needed, that is a signal the mapping in §2 is wrong, not that pysceptre needs a WattEG-shaped hole in
it.

## 6. Where the code lives

The simulation model — NB draw from `mean x size_factor x effect_size`, per-guide effect sizes,
re-centering, the seeding scheme — is **WattEG's method**, not part of sceptre, and does not go into
pysceptre. New package in this repo:

```
src/watteg/
  sim_input.py        # the container, ported from lib/sim_input.R
  expression.py       # poscounts size factors, normalised means, theta
  perturbation.py     # pert_input, guide status, effect-size matrix   (lib/pert_input.R, lib/simulate.R)
  simulate.py         # draw_counts
  engine.py           # the pysceptre call: pseudo-gene/pseudo-target assembly
  power.py            # Wilson interval, power, MDES                   (lib/stats.R, compute_power.R)
  seeds.py            # derive_seed / SeedSequence
  cli/                # one entry point per pipeline step
```

`pyproject.toml` with console scripts, so the Nextflow modules call `watteg-prepare-sim-input` etc.
rather than `Rscript src/...`. Output file names, columns and TSV/Parquet layouts stay **byte-
compatible** with the R path wherever they can — `consolidate_replicates`, `compute_power` and
`summarize_power` outputs are what the paper's figures read.

## 7. The DAG afterwards

```
samplesheet -> PREPARE_SIM_INPUT -> SPLIT_PAIRS -> POWER_SIMULATION (split x effect size x rep chunk)
                                                     -> CONSOLIDATE_REPLICATES -> COMPUTE_POWER -> SUMMARIZE_POWER
```

Eight processes become six; `FIT_NULL_MODELS` and `MERGE_NULL_MODELS` go, and with them
`reps_per_null_chunk`, `test_max_null_reps`, the divisibility check on them, and one join in
`main.nf`. `reps_per_chunk` stays and becomes load-bearing for memory (§2.3).

## 8. Validation — what would make this believable

R and Python cannot agree draw for draw: different RNGs, different CRT streams, and §3.1 changes the
null model relative to what the R sweeps ran. Validate in stages, against the reference outputs
already in `WattEG-paper` rather than re-running R.

**Which reference is which** — checked, because getting it backwards manufactures a discrepancy out
of a settled difference. `power_sweep/` holds the **six-effect-size sweep** (`power_es0.05` through
`power_es0.5`) and its `prepared/` carries `null_precomputations.rds`, so it ran the `null_fit`
configuration that §3.1 reproduces. `power_sweep_null/` holds `power_es0.0.tsv` only: it is the
**es = 0 null arm**, not the `null_fit` configuration. Stage B reads `power_sweep/`; Stage C reads
`power_sweep_null/`.

**"Reproduce R" has a floor that is not the port's doing.** The existing sweeps ran on an x86
cluster, where R's `sum()` accumulates in 80-bit `LDOUBLE`; on arm64 `.Machine$sizeof.longdouble`
is 8 and it accumulates in plain `double`. Measured in phase 2: **a local R run reproduces only 30
of 20,000 published size factors bit for bit**, and differs from them by up to 6.2e-12 — the same
residual the Python port shows. So the published outputs cannot be reproduced exactly by R either,
and 1e-10 is what "reproduces R" can mean across platforms. That is a definition, not a caveat, and
it applies to every stage below.

**No stage simulates in one language and tests in the other.** An earlier draft of this section
had one: dump Python-simulated counts and push the same matrices through both engines, so that any
difference was the engine alone. It is dropped, deliberately. It would have kept a working R
install, a pinned sceptre and a matrix-handoff harness alive purely to validate the thing that
exists to remove them, and it would have forced this pipeline's `sim_input` to carry the QC-failed
cells R's matrices span so R could index them. The Python path simulates and tests in Python.

What that gives up, stated plainly: nothing measures the two engines against each other **on
simulated counts specifically**, which are denser and lower-variance than real ones. The answer
comes by transitivity instead — pysceptre is already validated against R sceptre on this very
screen's real discovery analysis (`test_day0_regression`: Spearman 0.9865 on p-values, fold change
agreeing to 2.6e-12, sensitivity 0.9882 against R's own BH calls) — plus Stage B end to end and
Stage C, which needs no second implementation at all because it has an absolute bar.

**Stage 0 — measure the noise floor first.** Two independent runs of the same correct pipeline do
not agree pair for pair: power is a fraction over 100 Bernoulli draws, and near-threshold pairs
cross in both directions. Re-run the **Python** pipeline at a second seed on a small panel and
record its own Δpower spread and 0.8-line crossing count. That is the bar Stage B is read against,
and it costs one extra short run rather than an R install.

**Stage B — end to end, against R's published power.** Full 100 simulations at effect size 0.15,
Python against `power_sweep/.../power_es0.15.tsv`. This is now the only stage that compares the two
implementations, so it carries the weight Stage A used to share.
- Report: per-pair Δpower distribution, the fraction exceeding each pair's Wilson half-width, the
  mean shift, and the count crossing the 0.8 line in each direction.
- Acceptance: mean shift consistent with zero — two independent 100-draw estimates of the same
  binomial p differ by about `sqrt(2p(1-p)/100)`, so ≈0.07 per pair at p = 0.5 and ≈0 in the mean
  over 34,886 pairs — and a **symmetric** count of 0.8-line crossings. The 12.6 % of pairs status.md
  calls ambiguous will move in both directions; that is expected, not a failure. What would not be
  expected is a mean shift, or crossings running one way, which is precisely the signature `as_is`
  showed (11,649 up against 3,631 down).

**Stage C — the null arm.** Effect size 0 against `power_sweep_null/`: the empirical rejection rate
should sit at the nominal threshold in both. This is the pipeline's own calibration check, and the
one stage with an absolute bar rather than a relative one.

**Stage D — invariance.** 1x4 against 2x2 replicate chunking, and two split layouts, byte-identical
(§3.3). Unit test, not a cluster job.

## 9. Out of scope, said explicitly

- **Low-MOI screens.** pysceptre covers the complement-control-group + CRT high-MOI path only. A
  low-MOI object must fail at `prepare_sim_input` with a clear message naming the R path, not
  silently produce numbers from the wrong control group. The check reads `control_group_complement`
  and `run_permutations` out of the export metadata.
- **`--n-control-cells` and `--cell-batches`. Decided: neither is ported.**

  `--n-control-cells` draws a fixed number of control cells per target instead of using every
  non-perturbed cell — on day0, 5,000 in place of ~567,000. It is a cost lever and nothing else,
  and it was measured to cost **21–60 % of power**. `--cell-batches` stratifies that draw so the
  sampled controls keep the perturbed cells' batch composition; it does nothing on its own, and
  `run_power_simulation.R` refuses it without `--n-control-cells`.

  The question worth answering, because it is the one that sounds alarming: **does dropping
  `--cell-batches` expose the analysis to batch drift between the two arms?** No, for two reasons.

  1. With no subsampling there is no draw to stratify. The control group is every non-perturbed
     cell, so its batch composition is the dataset's, not a sampling artefact.
  2. Batch is conditioned on by the test itself, twice. `batch_factorBatch 2/3/4` and
     `replicate_factorRep 2/3/4` are columns of the covariate matrix, and that matrix enters both
     the per-gene NB fit — so batch effects on expression are adjusted out — and the logistic fit
     of perturbation status that the **CRT draws its synthetic treated sets from**. The null
     distribution is therefore conditional on batch by construction. That is the formal guarantee,
     and it is why this method does not need matched control cells: cell-level matching is what
     you reach for when the model cannot adjust for a confounder, and here it can.

  Stratified sampling was never the defence against batch confounding. It was a patch for the
  extra variance that careless subsampling adds on top of a model already handling it.

  The case against porting is stronger here than it was in R: the lever exists to buy speed, the
  port is the reason speed stops being the binding constraint, and a knob that trades power for
  speed you no longer need is a trap rather than an option. **`perturbation.py` should not grow a
  control-sampling path**, and `sim_input.h5` carries no `batch_factor` or `replicate_factor`,
  since `--cell-batches` was the only reader.

  **If control subsampling ever returns** — a screen large enough that even the Python path cannot
  afford the full control set — stratification has to return with it, and the reasoning above is
  why. That costs no format change: the design matrix holds both factors one-hot (`batch_factorBatch
  2/3/4` plus an all-zero reference level), so each is reconstructible in about fifteen lines, and
  `sim_input` already stores categoricals as codes plus levels for exactly this.

  The R implementation keeps both flags. It is the reference the paper describes, and removing
  options from it would change what that reference is.
- **`run_permutations = TRUE` screens.** pysceptre supports permutations, but its draws are sized by
  the largest target in the run, which interacts badly with per-target calls. Refuse for now.

## 10. Infrastructure

- `pixi.toml`: drop `r-base`, `r-optparse`, `r-matrix`, `r-rcpp`, `r-dplyr`, `r-data.table`,
  `r-purrr`, `r-crayon`, `r-parallelly`, `r-withr`, `r-nanoparquet`, `SCEPTRE_REF`/`SCEPTRE_SHA`,
  `ONDISC_REF`/`ONDISC_SHA`, and the `setup` / `check-api` tasks. Add `python`, `numpy`, `scipy`,
  `pandas`, `pyarrow`, `numba`, `mudata`. Keep `nextflow`.
- **pysceptre is still a pin.** It is on neither conda-forge nor bioconda, so it enters as a pixi
  `[pypi-dependencies]` git dependency pinned by commit — the same shape as `SCEPTRE_SHA`, minus the
  patch and the API check. `pixi.toml`'s comment block explaining why sceptre is pinned gets
  rewritten, not deleted, and `check_sceptre_api.R`'s job — assert the pin still matches what we
  call — passes to pysceptre's own test suite plus this repo's.
- A new container image for the `gcb` profile. The R image is not reusable.
- **`conf/*.config` needs recalibrating from scratch.** The last ten commits on `main` tuned
  `POWER_SIMULATION`'s memory and machine type around R's 2.27 GB median / 3.27 GB max. Python's
  footprint is dominated by the stacked count matrix (§2.3) and is a different function of
  `reps_per_chunk` and pairs-per-target. Do not carry the closures over; re-measure, then rewrite
  them.
- `.githooks/pre-commit` rejects camelCase **R** identifiers; add the Python equivalents (ruff,
  matching pysceptre's `ruff.toml`) rather than leaving Python unlinted.

## 11. Phases, with a gate that can stop the work

1. ~~**Export gap** (pysceptre branch): §5.1 and §5.2, plus a round-trip test that the individual
   targeting-gRNA unions reproduce `grna_group_idxs` exactly.~~ **Done** — squashed into
   `d96d48f` on `0.1.1rc`. *Gate passed:* all 3,071 day0 targets reproduce exactly,
   asserted at export time rather than in a test that can be skipped; the default read of an
   `--all-cells` export is identical to a plain one on the real screen; 230 tests green, with the
   export-format contract covered by 10 new ones that need neither R nor a real dataset; and
   `test_day0_regression` passes 6/6 against a re-export of day0 (4 min, 34,886 pairs), so the
   export changes move nothing the engine reads.
2. ~~**`prepare_sim_input` in Python** against the fixture: size factors, normalised means, theta,
   threshold, pairs — each compared to the R output column by column.~~ **Done.** *Gate passed*
   against the day0 `sim_input.rds` the day0 sweep was run on: `pairs.tsv`, `grna_targets.tsv`
   and `discovery_threshold.txt` **byte-identical**; genes the same set in the same order; all
   3,071 target and 43,718 guide cell sets agreeing exactly; and the expression statistics
   **bit-identical to a same-platform R run**, differing from the published ones only by the
   6.2e-12 platform residual above. Theta as in §3.2. Three defects the gates caught rather than
   luck: counts stored as `uint16` made `np.log` return **float32** (2.5e-7 on the size factors);
   `grna_perts` was missing the non-targeting guides, which would have kept the control arm's mean
   while losing its guide-level variance; and `pairs.tsv`'s column order. (The second is moot since
   2026-09-24: every guide outside the target now has an effect of exactly 1, so control cells
   carry no guide-level spread at all. See `methods.md`.)
3. ~~**Benchmark before committing to the shape.**~~ **Done, and it overturned §2.** The four
   terms, measured: the per-pair test is **85 %** of the work, the per-gene Poisson fits 14 %,
   drawing counts 5 %, and the per-target binomial fit and CRT draws **1.5 %** — of which the
   redundancy stacking would remove is **0.8 %**. So there is almost nothing to amortise, the two
   shapes measure the same, and **the simple one wins on simplicity alone**: `engine.py` makes one
   call per (target, replicate), with no pseudo-genes, no pseudo-target keys, no replicate-chunk
   memory knob and no question about replicates sharing a resampling stream. §5.3 is dead on its
   own terms — 0.8 % was the number it had to beat.
4. ~~**`run_power_simulation` in Python** + Stage A validation.~~ **Done.** Stage A was dropped
   (§8), so equivalence is measured on output: 40 replicates of one target, both implementations,
   every bound derived from the data's own spread. Per-pair power agrees 8/8 within Monte Carlo
   noise with no systematic shift, and each implementation's median fold change lands on
   `log2(0.85)` — an *absolute* check, which is the only kind that can catch both being wrong the
   same way.
5. ~~**The four cheap steps**~~ **Done** for `split_pairs`, `consolidate_replicates`,
   `compute_power` and `summarize_power`, checked against R **on identical input**, which makes
   them exact comparisons rather than statistical ones: every column agrees to machine epsilon,
   in R's column order. The comparison caught `max_effect_size_tested` missing entirely and the
   per-gene columns sitting in the wrong place. `fit_power_curve` is **not** ported — it serves
   the paper's reduced-design study rather than the pipeline, and nothing in the DAG calls it.
6. **Nextflow rewiring** — done: `FIT_NULL_MODELS` and `MERGE_NULL_MODELS` are gone, the six
   remaining processes call the `watteg-*` entry points, the samplesheet column is `dataset`
   (a `.h5mu`) rather than `sceptre_object`, and `pixi.toml` holds no R. **Still open: the new
   container and the resource recalibration.** The `conf/*.config` memory closures were tuned
   around R's 2.27 GB median and must be re-measured, not carried over — the Python footprint is a
   different function of the cell count and the gene count.
7. **Stage B/C/D validation** at full scale on one effect size. **Stage D is done** — four
   replicates in one task against two tasks of two are byte-identical, and a different seed does
   change the draws, so the invariance is not coming from the seed being ignored. It runs behind
   `-m realdata`. **Stage B passed at full scale on moi5 cis, 2026-09-25** (§12). Stages 0 and C
   are still open.
8. **Docs** — `README.md` and a pointer on `status.md` are done. `usage.md`, `methods.md` and
   `output.md` still describe the R path.

Phases 1–6 are done. Phase 7 is the one that decides whether `main` moves.

## 11b. Running Stage B

**The reference has to be regenerated.** This said the opposite until the centring bug was found,
and reading the existing `power_sweep/` tables was the whole reason Stage B was cheap. They were
produced by R before either fix, so they are a known-wrong reference: comparing against them would
report a difference that is real, correctly measured, and about the bug rather than the port.

That means one R run on the Stage B targets, from `r-implementation` at `b346296` or later, and a
**re-prepare** as well as a re-simulate — the fitted baseline reads `fitted_coefs`, and no
`sim_input.rds` made before `d767af3` carries them. Both sides then run on their own defaults,
which is what makes the comparison a comparison of implementations.

```sh
# 1. export the object with BOTH flags (once per dataset)
Rscript ../pysceptre/scripts/export_sceptre_dataset.R \
    --sceptre-object <sceptre_object.rds> --out-dir export/ --all-genes --all-cells
python  ../pysceptre/scripts/make_h5mu.py export/

# 2. prepare
watteg-prepare-sim-input --dataset export/dataset.h5mu --outdir prepared/ --n-jobs 8

# 3. simulate, on the default (fitted) baseline
watteg-run-power-simulation \
    --prepared prepared/ --pairs <split>.tsv \
    --effect-size 0.15 --reps 100 --seed 20250812 --n-jobs 8 \
    --out stageb_sim.tsv

# 4. the R reference, from the FIXED R, on the same targets and the same defaults.
#    Run in a worktree of r-implementation; it needs its own prepare, and its own
#    sceptre_template.rds and null-model fits, which the Python path does not have.
nextflow run <r-worktree> -profile <...> -params-file <...>

# 5. compare
workflow/compare_stage_b.py stageb_sim.tsv \
    <r-results>/power/power_es0.15.tsv \
    --threshold-file prepared/discovery_threshold.txt
```

**Both sides on the default baseline.** `--expression-model size_factor` used to be mandatory here,
to match a reference that predated the baseline change. With the reference regenerated it would do
the opposite of its job: it would take the Python side off the model the R side is now using.

**The R side can run short.** `compare_stage_b.py` computes the noise floor from each side's own
replicate count, so R at 20 replicates against Python at 100 is a valid comparison, just a blunter
one. Since R is the expensive side, that is where to spend less.

**Cost, from the run's own timings**: `0.48s + 0.120s x pairs` per (target, replicate) at
`--n-jobs 8`. A 36-target, 428-pair sample is about 1.9 h at 100 replicates and 23 min at 20.
The whole sweep — 3,026 targets, 34,886 pairs — is roughly 130 h on one machine, which is what the
cluster is for.

**What the replicate count buys.** The per-pair test compares two binomial estimates, so its floor
is `2*sqrt(p(1-p)/n_py + p(1-p)/n_r)` — at `p = 0.5` that is +/-0.32 with 20 Python replicates
against the reference's 100, and +/-0.13 with 100. The *aggregate* test is far sharper either way:
the shift in mean power over 428 pairs has a 2-se bound near 0.015 at 20 replicates. So a short run
already tests the thing most likely to be wrong -- a systematic bias from the port -- and a long
one is what makes the per-pair claim worth stating.

## 12. What is left, and what would be wrong to skip

**The two that block a real sweep.**

- **Resource recalibration.** `conf/base.config`'s numbers are R's, now labelled as such rather
  than left looking calibrated. A task uses `task.cpus` workers where R's used one, so both the
  memory and the time models are a different shape. They err high, which wastes budget rather than
  losing runs, but they are not a calibration until a real trace replaces them.
- **A container for the `gcb` profile.** The R image is not reusable and nothing has been built.

**The one that decides whether this replaces the R path.** Stage B at full scale: one effect size,
100 replicates, all 34,886 pairs, against a **regenerated** R sweep (§11b). Everything measured so
far says the two agree — but on one target, eight pairs, forty replicates. That is evidence the
paths agree where they have been compared, and it is not the same claim.

Stage B's first run at 36 targets **failed**, and that failure is what found the centring bug. Two
things follow that are easy to conflate. The first is that the failure was correct and its
diagnosis — 27 % less replicate spread in Python, with a floor that did not shrink with cell count,
pointing at a per-guide term — is now a **prediction**: a re-run against fixed R should show the
spread match, and if it does not, the diagnosis was wrong rather than incomplete. The second is
that a pass would confirm that diagnosis *and* the port at once, which is weaker than it sounds and
worth saying out loud rather than reporting as a clean green.

> **Resolved 2026-09-24: the 27 % was an estimand difference, not a Python bug.** Python pinned the
> realised mean knockdown (a fixed element effect) from its first commit. The R code it was compared
> against, like the original DC-TAP code, centred on the wrong columns, so its realised mean was free
> to vary. On a real screen that is the same as not centring at all. Moving R's reorder before the
> centring (`2b76284`) made R pin too, and so silently changed which quantity R simulated; no one
> had chosen it. The owner has now chosen the fixed element effect for both languages, with control
> cells at exactly 1 (see `methods.md`, "What simulated power means"). A Stage B re-run therefore
> compares two implementations of the same estimand, and is expected to agree.

### Stage B at full scale, 2026-09-25: the two implementations agree

Both pipelines on moi5 cis: every one of the 33,066 pairs, es 0.15, 100 replicates each, the
screen's own permutation test, `guide_spread_c = 0.65`, same seed. R ran at `1f44a42` (run
`chaotic_lovelace`), Python at `6570024` (run `intergalactic_liskov`). The two prepares wrote
byte-identical `pairs.tsv` and the same discovery threshold, 0.000648397.

| check | result |
|---|---|
| per-pair power within 2 se of Monte Carlo noise | **96.9 %** of pairs (≈95 % expected by chance) — pass |
| status at the 0.8 bar | **98.1 %** agree; 320 Python-only, 318 R-only, every one sitting on the bar — pass |
| mean power | Python 0.6040, R 0.6045 |
| median fold change | Python log2 −0.2352, R −0.2359, target log2(0.85) = −0.2345 |
| systematic shift | **−0.00055**, against a 2-se bound of 0.00047 (0.00048 clustered by target) — formally a fail |

The shift is real at about 2.3 se and is 0.06 points of power. It is flat across expression
quintiles and sits in mid-power pairs (−0.19 points at power 0.3–0.7), the shape a small calibration
difference in the test gives. One is known and expected: R hoists the gRNA null model out of the
simulation (`FIT_NULL_MODELS`, fitted once per replicate on an independent null simulation), while
Python refits it inside every call (§3.1). A port bug would not show up as 0.06 points spread
evenly over 33,066 pairs with the 0.8 crossings symmetric to within two pairs. **Verdict: the Python
path computes the same power as the R path.** What Stage B does not cover: trans, where there is
no R sweep on this code, and other screens.

**Two pipeline bugs that only a real run could find**, both fixed; neither changes a number:

- `prepare_sim_input` rebuilt `set(pairs["response_id"])` once per gene: 108 min on moi5 trans
  (38,606 genes × 742,525 pairs), past the task's 1 h limit, against 92 s with the set built once
  (`f2d684c`). cis hid it at 33,066 pairs.
- `COMPUTE_POWER` listed its inputs with `$(ls a b c)`. Only one pattern ever matches, `ls` exits 2,
  and `bash -e` killed the task on that line with empty stdout and stderr (`efd0f05`). Every real
  Python run so far died there; the stub never runs the script. The cis power table above was
  computed from the pipeline's own consolidated Parquet with the module's exact command.

**Cloud cost, measured, which replaces the R numbers in `conf/base.config` as a starting point.**
On `e2-standard-2` spot a POWER_SIMULATION task runs ~0.37 s per pair-replicate on cis, about 4×
the laptop's rate. Trans targets carry ~250 genes each, and there `--n-jobs` pays: on one target,
40 s per replicate at 1 worker, 26 at 2, 19 at 4 and 16 at 8, with byte-identical output. So trans
runs at `cpus = 4`, `reps_per_chunk = 10`; at 1 CPU a 20-replicate trans task needs ~2.5 h and would
hit the 2 h limit on every task. The memory closure asks 8 GB for any non-empty split (the predicted
figure always clears the 4 GB floor), while the 1,000 Python cis tasks peaked at 1.1-1.3 GB,
20-24 min each.

**Three scientific questions are open and recorded, none of them acted on** (§3.2b, §3.2c, §5.2).
The largest, simulating from `exp(X . beta)`, has been taken; the other two are judgement calls
about which cells and which estimator belong, and both change every published number.

**One thing deliberately not ported.** `fit_power_curve.R` serves the paper's reduced-design study
rather than the pipeline, and nothing in the DAG calls it. It stays in R.

**The synthetic fixture.** `src/make_test_data.R` produces a sceptre object; the Python path needs
a `.h5mu`, so `assets/samplesheet_synthetic.csv` points at a file nothing generates yet. The stub
run works from a real export instead.

**The numbers in `WattEG-paper`.** Both fixes change results, and every sweep in `power_sweep/` was
produced before them, as was §4 of `perturbplan_comparison.md`. Nothing was published, so this is a
regeneration rather than a correction — but it is a full sweep, and it is recorded in
`WattEG-paper/docs/bug_expression_scale_input.md` rather than here because it is the paper's work,
not the port's.

## 13. Next: both estimands, and cis + trans in under an hour

Decided with the owner on 2026-09-25. Items 1 and 2 are settled; items 3-5 wait on the speed
measurements in progress, and their numbers are estimates until then.

**Status, 2026-09-25 (end of day).** Done: item 1 in Python and in R (`r-implementation` `1a63eb3`); item 2 (`e6c6db6`);
3a (`3ea5baa`); 3b, the fast driver (`288cb0f`); 3c, fit reuse; item 4's per-pair power inside the
task, and its memory measurement. **The defaults are now the fast configuration:** `--permutations
per-target`, `--nulls sparse`, `--driver fast`, `--null-fits reuse` (CLI and Nextflow), with the
engine, `refit`, `scan` and `per-replicate` still selectable. A CRT screen needs `--driver engine
--null-fits refit`, because the fast driver runs the permutation test only. Not done: 3d (waits on
a pysceptre release) and item 6 (not to be run until decided).

**Item 5, measured on the cloud the same evening (`0651653`, seed 20250812, 8 workers per task):**

| | cis (`focused_hoover`) | trans (`maniac_bassi`) |
|---|---|---|
| wall clock, launch to summary | ~1 h | 1 h 52 min |
| simulation tasks | 100, median 11.9 min | 1,000, median 21.8 min, max 29.6 |
| cost | $2.48 | $67.57 (the previous trans run: $411, 5.6 h) |
| vs the previous Python run, same seed | +0.0003 mean power, 99.5 % same 0.8 call | +0.0002, 99.56 % same 0.8 call |

- **Under an hour:** trans came in at 1 h 52 min, not under an hour. About 20 min of it is prepare
  and the fitting step, which run before any simulation starts.
- **Per-task cost:** the simulation tasks ran at about 1.7x their laptop cost on 8 e2 vCPUs.
- **Laptop, single core, cis reference target:** 182 ms per pair-test before, 105 with the sparse
  route, 70 with the fast driver, and 31.5 with fits reused from the file.

**What does not change: the simulation keeps the full design.** Each replicate runs the screen's
actual test on its actual cells, covariates, guide assignment and threshold, including sceptre's
permutation test with its escalation. Every speed-up below has to preserve that. None of them
approximates the design.

### 1. A second estimand: random guide effects -- done in Python and R

**Done 2026-09-25 in Python** (`watteg.perturbation.effect_size_matrix(estimand=)`, `--estimand`,
Nextflow `estimand`; written as the last column of every per-simulation row, into the power table
and the summary, and refused when mixed). Checked on the shared fixture and the cis reference target:
- es = 0: identical rows under both (engine and fast driver, fixture; CLI, 20 simulations), apart
  from the estimand label; the generator is left in the same state, so the counts are drawn alike.
- `fixed` output unchanged by the option: every earlier output (cis at es 0 and 0.15, trans, the
  engine's per-replicate scan run) is today's file minus its new last column, byte for byte.
- `random`: the realised mean's sd over 2,000 draws matches `c * es * (1 - es) * sqrt(sum n_g^2) /
  sum n_g` within 8 % at es 0.05, 0.15 and 0.5, and control cells stay exactly 1.
- cis reference target, 100 simulations, es 0.15: 808 calls of 1,200 under `fixed`, 781 under
  `random` (mean power 0.673 against 0.651); the simulated fold change's sd over simulations rises
  from 0.056 to 0.067 (median over pairs).

**Still to do:** the R side (`simulate_effect_sizes` skipping `center_effect_size_matrix`) on the
`r-implementation` branch, so Stage B can check it; the cis sweep under `random` scored against
PerturbPlan (last check below), which waits on item 6.

`--estimand fixed | random`, default `fixed`, in both implementations so Stage B can check it.

| | `fixed` (today) | `random` |
|---|---|---|
| guide knockdown | Beta, mean es, sd `c * es * (1 - es)` | the same draw |
| element's realised mean over its perturbed cells | pinned to es in every replicate | left where the draw puts it |
| question answered | power for an element whose effect *is* es | power for an element whose effect is es *on average* |
| PerturbPlan setting that asks the same | `fold_change_sd = 0` | `fold_change_sd = c * es * (1 - es)` (0.0829 at es 0.15) |

**Implementation.** Skip the pin and keep everything else.
- **Python:** `watteg.perturbation.effect_size_matrix` skips `_pin_to_mean`.
- **R:** `simulate_effect_sizes` skips `center_effect_size_matrix`.

The control-cells-at-1 assertion stays in both. The estimand is written into every output, so a
power table always says which question it answers.

**Why it is safe now and was not before.** The unpinned spread is zero at es = 0. The old absolute
`N(1 - es, 0.13)` gave a null element random effects, and "detected" it about 20 % of the time on
highly expressed genes. With the Beta spread the null arm is identical under both estimands.

**Checks:**
- es = 0 gives byte-identical output under both estimands.
- Under `random`, the realised element mean has sd ≈ `c * es * (1 - es) * sqrt(sum n_g^2) / sum n_g`
  over the guides' cell counts `n_g`: the variance of a cell-weighted mean, tested on the shared
  fixture.
- `fixed` output is unchanged by the option's existence.
- A cis sweep under `random`, scored against PerturbPlan at `fold_change_sd = 0.0829`, compares the
  same question at the same per-guide spread in both methods (WattEG-paper,
  `docs/perturbplan_comparison.md`).

### 2. One permutation set per target, drawn from the seed -- done, and the default

Done in `e6c6db6`; the default since the fast configuration was adopted (2026-09-25), and recorded
in `methods.md`.

All replicates of a target are tested against one permutation set, keyed on (seed, target, effect
size), not on the replicate. This is what sceptre itself does: its sampler reseeds
`mt19937(4)` on every call (pinned commit 3ba046b), so R's replicates already shared one fixed set,
as did the real screen. A per-target draw from the run's seed was chosen over sceptre's fixed set
so that one draw is not shared by every target of the same size. Outputs change relative to earlier
Python runs; that is expected and is recorded in `methods.md`.

### 3. A faster driver

#### 3a. First: take pysceptre's sparse route for the permutation nulls (measured, bit-identical) -- done

Step 1 done in `3ea5baa` (`--nulls sparse`, now the default); step 2 is the fast driver (`288cb0f`).

**The problem.** In a one-target call pysceptre computes the stage-1 and stage-2 permutation nulls
by a scan: it gathers a (B, n_trt, 14) array (227 MB at stage 2, n_trt = 406) and takes a running
sum over it, to read a single column (`PermutationPrefixSums`, `score_stat.py:411-441`). The scan
exists to serve many targets of different sizes from one set of draws. With one target it is pure
waste. pysceptre already has the alternative, a sparse indicator matrix times the gene's pieces
(`draws_to_matrix(...) @ stacked`, `score_stat.py:82-115, 274-300`). It uses that route whenever
the scan would exceed its memory limit, and the real moi5 screen takes it at stage 2 on its own,
because its largest target (607 cells) is above the limit (457).

**Measured** (2026-09-25, single-core laptop, `speed/timing_decomposition/`). The route was switched
by making the scan decline at runtime, so pysceptre ran its own sparse code:

| | today | sparse route | |
|---|---|---|---|
| stage-2 null, per escalated pair | 119.5 ms | 17.7 ms (4.0 build + 13.6 multiply) | |
| stage-1 null, per pair | 12.4 ms | 1.9 ms | |
| **cis**, per simulated pair-test (12 pairs x 100 simulations) | 182 ms | 87 ms | **2.1x** |
| **trans**, per simulated pair-test (255 pairs x 10 simulations) | 154 ms | 66 ms | **2.3x** |

p-values, z-statistics and stages were identical on all 1,200 cis and 2,550 trans pair-tests,
max |dp| = 0.

**Implementation.** Two steps, in order:
1. **Now, in the current engine.** Steer pysceptre onto the sparse route with a scoped, documented
   runtime setting in `watteg.engine` (the scan's decline threshold), leaving pysceptre's source
   untouched. It needs a test that p-values equal the scan route's on the shared fixture and on
   one real target. This alone halves the cost of every sweep.
2. **Then, in the faster driver.** Build the sparse matrix once per target and stage and apply it
   to every gene x simulation of that target in one product (the table below). With one
   permutation set per target (item 2) the matrix is shared across simulations, and the build cost
   (4.0 ms per pair-test today) disappears too.

If pysceptre's owners want it, the same finding applies upstream: for a single target, the scan is
slower than the sparse route at every stage. That is a note for them, not a change made here.

#### 3b. The rest of the driver -- done (`288cb0f`, `--driver fast`, now the default)

A read-only map of pysceptre's discovery call found work that is repeated for no reason in the
simulation's call shape (one target, one replicate per call). All of it is avoidable without
editing pysceptre, by driving its public low-level functions (`fit_all_genes`,
`compute_precomputation_pieces`, `stack_pieces`, `run_low_level_test_full` with
`null_statistics_fn`):

| work today | cost (laptop, cis profile) | replacement |
|---|---|---|
| all 30,497 permutations drawn per call, one `rng.choice` at a time | ~0.24 s per call | drawn once per target (item 2) |
| stage-1 and stage-2 nulls on the scan route | 119.5 / 12.4 ms (measured, 3a) | 3a first; then one sparse matrix per target and stage, `P @ [stacked_1 | ... | stacked_K]` for every gene x simulation at once; bit-identical per column (checked) |
| gene fits one replicate at a time, plus a duplicated `compute_precomputation_pieces` | ~44 + 5 ms per gene per replicate | fits batched across a target's replicates; not guaranteed bit-identical on Linux, so checked to tolerance |

Measured breakdown of today's 182 ms per cis pair-test, against 1.2 ms in the real discovery of
the same screen:

| cause | ms | share | fix |
|---|---|---|---|
| scan route, stage 2 | 83.5 | 46 % | 3a |
| gene refits | 49.3 | 27 % | batching across simulations, or R-style fit reuse |
| per-call overhead (draws) | 19.6 | 11 % | draws once per target |
| escalation rate (82 % of simulated pair-tests vs 3.3 % real) | 15.3 | 8 % | none: the cutoffs are below 1/500, so a callable pair must reach stage 2 |
| scan route, stage 1 | 10.5 | 6 % | 3a |

- **R-style fit reuse, measured:** 1.24x on cis and 1.27-1.29x on trans. The call rate is
  identical on cis (806 of 1,200), 8 pair-tests flip 4 each way, and p-values move a quarter as
  much as a change of permutation seed does.
- **Estimated with everything combined:** about 25 ms per pair-test with fit reuse, or about 35-40
  ms with batched exact fits. **Exactness criterion:** on the
same permutation set, the fast driver's p-values match the current engine's per pair, bit for bit
where the operations are the same and to a stated tolerance where the fits are batched.

#### 3c. Then: reuse each gene's null-model fit (decided 2026-09-25) -- done

**Done 2026-09-25** (`watteg.null_fits`, `watteg-fit-null-models`, the `FIT_NULL_MODELS` process,
`--null-fits reuse | refit`, `--null-fits-file`). Checked:
- `refit` is byte-identical to the fast driver before the option existed, on the cis reference
  target (100 simulations, es 0 and 0.15), the trans reference target (10 simulations) and the
  engine's per-replicate scan output.
- The fits reproduce the benchmark's (`speed/hoist_null_fit/fits_cis.tsv`): 34 of 48 exactly, the
  rest within 3e-15, the difference being the benchmark's many-gene baseline product.
- A fit depends only on (gene, simulation): a gene fitted alone equals it fitted with others; a
  file made for more genes and simulations, and fits made in the task at 1 and 2 workers, give the
  same bytes (unit test and realdata test); a file from another seed is refused.
- `reuse` against `refit`, same counts and permutations (per-target), on the fast driver:

  | | refit calls | reuse calls | flips | McNemar p | max per-pair \|d power\| | beyond Wilson half-width |
  |---|---:|---:|---:|---:|---:|---:|
  | cis, 12 pairs x 100 | 808 | 805 | 5 / 2 | 0.45 | 0.02 | 0 of 12 |
  | trans, 255 pairs x 10 | 1,280 | 1,282 | 3 / 5 | 0.73 | 0.10 | 0 of 255 |

  Worker time at 8 workers, fit step excluded: cis 129.0 -> 74.2 s (1.74x), trans 236.5 -> 143.9 s
  (1.64x), more than the benchmark's 1.24-1.29x because the fast driver removed the other overheads
  the fit was averaged against. A task that fits its own genes pays the fits back (cis reference
  target: 12.1 s of fits at 8 workers), so the saving on cis comes from `FIT_NULL_MODELS` fitting
  each gene once per sweep.

**Why.** The gene refit is the largest cost left after 3a and 3b: about 49 ms of every simulated
pair-test, 27 % of today's cost. That is because each simulation refits every gene once per target.
The real screen fits each gene once for about 135 targets.

**What.** R's `FIT_NULL_MODELS` approximation, in Python:
- Fit each gene's null model once per simulation, on an independent es = 0 draw keyed
  `rng_for(seed, "__null_fit__|" + gene, rep, 0.0)`.
- Use pysceptre's own `fit_all_genes`, at the same batch width and with the same `x_outer_flat`.
- Reuse that fit for every target the gene is tested against, and every effect size.
- The fits are a small file: about 2.5 MB for 244 genes x 100 simulations.

**Pipeline.** A `FIT_NULL_MODELS` step between `PREPARE_SIM_INPUT` and `POWER_SIMULATION`, about
25-30 CPU-minutes per sweep, feeding every simulation task. The fast driver takes the fits in place
of calling `fit_all_genes` itself. `--null-fits reuse | refit` keeps the exact refit available;
the default is `reuse`.

**Measured cost and fidelity** (`speed/hoist_null_fit/`, laptop, same counts and permutations):
- **Speed:** 1.24x on cis, 1.27-1.29x on trans.
- **Calls:** the call rate is identical on cis (806 of 1,200 pair-tests), with 8 flips, 4 each
  way. trans calls 1,284 against 1,283.
- **p-values:** they move about a quarter as much as a change of permutation seed.
- **Bias:** no directional bias is detectable.

**Checks:**
- `refit` stays byte-identical to today.
- `reuse` against `refit` on the cis and trans reference targets: call rate and per-pair power
  within Monte Carlo noise.
- A fit is keyed only on (gene, simulation), so no result depends on which other genes or targets
  share its task.

#### 3d. After pysceptre's next release: retire `--nulls`

pysceptre `dev` (0de1bff, broadinstitute/pysceptre#2) now chooses the scan only when it pays for
itself, which a single target never does. Once that is released:
- bump the pin and rebuild the image;
- drop the `--nulls` switch;
- keep the real-data test that compares the two routes.

The issue stays open for drawing permutation stages lazily.

### 4. Task shape and resources -- per-pair power inside the task: done

**Done 2026-09-25: per-pair power inside each task.** When a task holds all simulations of its pairs
(`reps_per_chunk == num_replicates`, the default), `POWER_SIMULATION` writes per-pair counts
(`--partials-out`: simulations called, simulations used, the fold-change and cell-count sums, the
simulation range) and `COMPUTE_POWER` adds them up (`--partials`). The per-simulation table and
`CONSOLIDATE_REPLICATES` run only under `keep_per_simulation` (default false) or when simulations
are chunked across tasks. Checked: the power table from the counts equals the one from the rows,
byte for byte, in unit tests, a realdata test and a real local Nextflow run on moi5 cis (the
counts' table, and `watteg-compute-power` on the published Parquet of a `keep_per_simulation` run).
That needed one fix on the rows' side: pandas' default CSV parser returned a neighbouring float for
about half of a TSV's values (285,883 of 500,000 measured), so the rows are now read with
`float_precision="round_trip"`. The stub DAG passes in all four shapes (default, kept rows,
chunked simulations, fast driver with fits), and a tiny real run passes on moi5 cis (defaults) and
on day0 (a CRT screen, `--driver engine --null_fits refit`).

**Measured 2026-09-25: memory and time of a task under the new defaults.** A 330-pair moi5 cis slice
(`watteg-split-pairs --n-splits 100`, split 1: 29 targets, 185 genes), 100 simulations, fits read
from a `FIT_NULL_MODELS` file, per-pair counts only:

| | workers | peak memory | wall | worker time per pair-test |
|---|---:|---:|---:|---:|
| Linux container (fork; cgroup `memory.peak`) | 8 | **3.38 GB** | 278 s | -- |
| macOS laptop (spawn; summed RSS of the tree) | 8 | 9.73 GB | 296 s | 47.5 ms |
| macOS laptop | 12 | 10.2 GB | 272 s | 54.1 ms |
| macOS laptop, 3 such tasks at once | 3 x 4 | 4.1-4.2 GB each | 727 s for all 3 | 70 ms |
| `FIT_NULL_MODELS`, 244 genes, Linux | 8 | 3.04 GB | 59 s (24 simulations) | -- |
| `FIT_NULL_MODELS`, 244 genes x 100, macOS | 8 | 6.9 GB | 208 s | -- |

`power_simulation_memory` stays **8 GB**: 2.4x the Linux peak, and no 8-vCPU predefined machine has
less. The macOS numbers are local-run numbers: a spawned worker loads its own copy of the inputs
(~0.6 GB) where a forked one shares the parent's. A task's wall time is bounded by its largest
target (221 s of the 296 s at 8 workers), so on one machine several smaller tasks at once use the
cores better than one wide task: 3 x 4 workers finished 3 slices in 727 s against 3 x 296 s one
after another. trans (~250 genes per target) is not measured under this configuration.

- **A task holds whole targets with all their replicates,** and writes per-pair power directly. The
  74M-row consolidation and power steps become optional, kept only for a per-replicate table when
  asked for.
- **8 workers per task** (`a85188a`, pushed): cis tasks fit in 8 GB. trans at 4 workers
  averaged ~7.4 GB per task, so trans memory is set from a measurement with the fast driver, not
  from the discovery benchmark.
- **Task sizing:** roughly 15-20 min at 8 workers, within the 10,000 preemptible-CPU quota.
- **Machine family:** e2 cores ran ~3x slower than the laptop; n2d or c2d are measured on one small
  cloud run before a sweep.

### 5. The target

Both moi5 sweeps, cis (3.3M pair-tests) and trans (74M), in under an hour of wall clock, about 30-40
min of it simulation. On the estimates above this needs a pair-test of ≤ 0.15-0.2 s per vCPU on
cloud against 0.55 s today. The measurements decide whether items 3-4 reach it.

### 6. Re-check the guide spread, value and form (not yet run)

The spread behind 0.0829 at es 0.15 is `c * es * (1 - es)` with c = 0.65 (`methods.md`,
"Guide-to-guide variability").

**Settled by the 2026-09-25 data audit** (WattEG-paper `analysis/guide_spread/build_data.R`):
- **The enhancer bins were noise-corrected, like the null pairs.** Each pair's spread is
  `tau2 = var(fc) - cf_t * mean(se^2)`, where `cf_t` is calibrated on groups of 14 non-targeting
  guides per expression tertile (`guide_spread_v2.R:97`, `:78-84`), and the bins pool `tau2`.
  `methods.md` did not say so, and should.
- **The curve runs high at large effects.** It peaks at es 0.5 (sd 0.163), while the pooled bins
  peak near es 0.37 (0.144) and fall to 0.108 at 0.59. The c that best fits the bins is 0.593,
  against the 0.65 fitted per pair.

**Still to do, before anything else depends on 0.65:**
- **Is the calibration adequate?** `cf_t` comes from non-targeting guide groups. Check it against
  held-out non-targeting groups and against the null elements, per screen.
- **Refit both the value and the form:**
  - reconcile the per-pair fit (0.65) with the bin fit (0.593);
  - test whether a form that falls faster at large effects fits better than `c * es * (1 - es)`.
- **Carry any change through:**
  - the default `guide_spread_c`;
  - the PerturbPlan setting matched to the random estimand (item 1);
  - the 0.0829 in WattEG-paper `docs/perturbplan_comparison.md`, `paper/06_prediction.md`, and the
    figure `figures/guide_spread/guide_spread.png`.

Under the fixed estimand the stakes are small: no moi5 pair's power moves by more than 0.02 between
this spread and none. Under the random estimand, and in the PerturbPlan comparison, the value
matters directly.
