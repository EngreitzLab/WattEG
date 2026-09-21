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
| `WattEG` | **`feat/pysceptre-backend`** (created for this work) | `main` stays the R implementation. `WattEG-paper` was written against it and every number in that paper has to remain reproducible from `main` without archaeology |
| `pysceptre` | **`feature/watteg-simulation-support`** in `../pysceptre` (to create; matching its existing `feature/` prefix) | `pysceptre-paper` freezes pysceptre `main` for the same reason. The changes in §5 land there and merge only once the equivalence checks in §8 pass |
| `WattEG` | `legacy` | untouched (the Snakemake implementation) |

The R path is not deleted when the Python path lands. It is the reference the Python path is
measured against, and it is what the paper describes; retiring it is a separate decision, taken
after §8 reports.

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
**Acceptance gate:** per-gene `theta` from `nb_theta.py` against `@response_precomputations$theta`
on the fixture object, reported as a distribution, before anything downstream is trusted.
`build_dispersion_vector`'s hard error on missing/non-finite dispersions carries over.

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

**Status: done, on `feature/watteg-simulation-support` in `../pysceptre`** (commits `7e29efa`,
`5117651`). Both gaps are closed and verified against the real day0 object; §5.1 and §5.2 below are
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

**Stage 0 — measure the noise floor first.** The 0/265 flips `null_fit` achieved against `cleared`
were possible only because both ran the *same* R RNG stream: identical CRT index sets, differing
only in the null coefficients. R sceptre against pysceptre is two **independent** resampling streams
(4,999 draws plus a skew-normal fit each), so near-threshold p-values will cross the threshold
exactly as two R runs with different seeds do. Re-run the R panel with a second seed and record its
own flip rate and Δpower spread. That is the bar. Every criterion below is expressed against it
rather than against a number picked in advance.

**Stage A — isolate the engine.** Simulate counts in Python, dump the matrices (indexed to match
`template@cells_in_use`, §5.2), and test the *same* matrices through both engines — R sceptre with
`@response_precomputations` cleared, and pysceptre. With the simulation held fixed, any difference is
the engine and its resampling stream.
- Report: max |Δp|, median |Δp|, Spearman, and the **threshold-flip count** at the discovery
  threshold — the statistic `threshold_check.R` already produces.
- Acceptance: flip rate and |Δp| spread **inside the stage-0 floor** on the same 3-target x 53-pair
  panel, with no directional bias in the flips. R's 7–0 against `as_is` is what bias looks like; a
  4–3 split is not.

**Stage B — end to end, independent runs.** Full 100 simulations at effect size 0.15, Python against
`power_sweep/.../power_es0.15.tsv`.
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
- **`--n-control-cells` / `--cell-batches`.** Measured to cost 21–60 % of power and off by default;
  recommend not porting, and deleting the flags rather than carrying dead paths. Your call — say so
  if they should survive.
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
   targeting-gRNA unions reproduce `grna_group_idxs` exactly.~~ **Done** — `7e29efa`, `5117651` on
   `feature/watteg-simulation-support`. *Gate passed:* all 3,071 day0 targets reproduce exactly,
   asserted at export time rather than in a test that can be skipped; the default read of an
   `--all-cells` export is identical to a plain one on the real screen; 230 tests green, with the
   export-format contract covered by 10 new ones that need neither R nor a real dataset.
2. **`prepare_sim_input` in Python** against the fixture: size factors, normalised means, theta,
   threshold, pairs — each compared to the R output column by column. *Gate: §3.2.*
3. **Benchmark before committing to the shape.** 3 targets x 100 simulations on the fixture, timing
   (a) one call per (target, replicate), (b) replicates stacked per §2, at 10/25/50/100 replicates
   per chunk, with peak RSS — and **broken down into the four terms of §2.5**, not reported as a
   single wall clock. *Gate: this sets `reps_per_chunk`, gives the first honest estimate of CPU-h
   per effect size against R's 634, and is the only evidence that could reopen §5.3. If it is not
   comfortably ahead of R, stop and re-plan rather than porting the remaining four steps.*
4. **`run_power_simulation` in Python** + Stage A validation.
5. **The four cheap steps** (`split_pairs`, `consolidate_replicates`, `compute_power`,
   `summarize_power`, `fit_power_curve`) — mechanical, byte-compatible outputs.
6. **Nextflow rewiring**, new container, resource recalibration (§10).
7. **Stage B/C/D validation** at full scale on one effect size.
8. **Docs**: `usage.md`, `methods.md` and `status.md` rewritten for the Python path; the R path
   documented as the reference implementation it has become.

Phases 1–4 are the work; 5–6 are mechanical; 7 is the one that decides whether `main` moves.
