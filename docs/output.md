---
title: Output
nav_order: 3
---

# Output

## `power_summary.tsv` — the main result

One row per element–gene pair, one set of columns per effect size.

| Column | Meaning |
|---|---|
| `grna_target` | Perturbation target (element), e.g. `chr2:201986246-201986547`. |
| `response_id` | Gene. |
| `mean_pert_cells` | Mean number of perturbed cells across replicates. |
| `average_expression_all_cells` | Raw mean expression of the gene across all cells. Not size-factor normalised — see the note below. |
| `gene_mean` | The gene's size-factor-normalised mean (`mu`). Present only when `--sim-input` was supplied. **Not the same as `average_expression_all_cells`**: they correlate at r = 0.9999 but differ by a scale factor, and `mu` is the one the theory below uses. |
| `dispersion` | The gene's negative-binomial dispersion, **`1/theta`** rather than `theta` — the same convention `sim_input.rds` and `rnbinom(size = 1/dispersion)` use. Present only when `--sim-input` was supplied. |
| `power_at_effect_size_15` | Power at a 15% knockdown. One column per effect size; the suffix is `effect_size × 100`, with any decimal point written as an underscore (0.125 → `power_at_effect_size_12_5`). |
| `power_at_effect_size_15_ci_low`, `_ci_high` | 95% Wilson interval for that estimate. |
| `power_at_effect_size_15_n_reps` | Replicates contributing to it. |
| `min_detectable_effect_size` | Smallest **tested** effect size from which `power` reaches `--power-threshold` (default 0.8) and stays there. `NA` means no tested effect size did. |
| `min_detectable_effect_size_ci_low` | Optimistic edge, from `power_ci_high`. |
| `min_detectable_effect_size_ci_high` | **Conservative edge, from `power_ci_low` — the column to use when interpreting a negative.** |
| `max_effect_size_tested` | The largest effect size in the run, so `NA` above can be interpreted. |

`NA` in any of the three is a statement about the effect sizes you ran, **not** evidence that a pair
is undetectable. If you tested 0.15 and 0.2 and a pair needs 0.4, it will be `NA`. That is why
`max_effect_size_tested` sits next to it.

Two things about how these are derived are easy to get wrong, so they are stated explicitly.

**The direction of the interval inverts.** Power rises with effect size, so a *lower* bound on power
gives a *larger* minimum detectable effect size. `min_detectable_effect_size_ci_high` is therefore
computed from `power_ci_low`, not from `power_ci_high`. It is the pessimistic reading — the smallest
knockdown the data can *certify* you would have caught — which is why it, not the point estimate,
is what a false-negative argument rests on. The bracket is built by thresholding each effect size's
own 95 % interval, so read it as a conservative-to-optimistic range rather than an exact 95 %
interval for the effect size itself.

**A pair must clear the threshold at its effect size and at every larger one tested.** Reporting the
first effect size that clears, in isolation, lets Monte-Carlo noise decide: at 100 replicates a pair
whose true power is 0.75 clears 0.8 about a third of the time, so across six effect sizes a spurious
early clear is likely, and the error only ever runs one way — reporting the pair as more detectable
than it is. Effect sizes a pair was not tested at count as unknown, not as failures, so they do not
block the run.

### Why `gene_mean` and `dispersion` are here

Power for a pair is governed by

```
SE^2 ~ (1 / n_pert_cells) * (1 / mu + 1 / theta)
```

so anything that models power from covariates needs the gene's dispersion alongside its expression
and the perturbed-cell count. Until these columns existed, only `average_expression_all_cells`
reached this table and every such analysis had to load a 16 MB `sim_input.rds` to find the rest.

With the convention above, the bracket is `1/gene_mean + dispersion` — no conversion, because
`dispersion` already *is* `1/theta`.

They are joined by `summarize_power.R` at summary time rather than emitted per replicate, because
they are per-gene constants: putting them in the simulation output would only help sweeps run after
the change, while re-summarising an existing sweep from its stored power tables takes seconds. Both
columns are absent if `--sim-input` is not passed, so older outputs and hand-run summaries are
unaffected.

## Interpreting a negative result
{: #interpreting-negatives }

The usual reason to run this pipeline is to decide what a *non-significant* pair means: did CRISPRi
fail to show an effect because there is no regulatory link, or because the experiment could not have
seen one? Power is what separates those, but only if it is read the right way.

### Per pair — threshold the lower bound

The claim "this negative is biological" is a claim that power was **at least** 0.8, not that the
point estimate landed above it. So threshold `power_ci_low` (or take
`min_detectable_effect_size_ci_high`, which does it for you across effect sizes). That is a stricter
test than it looks, and the strictness is the reason replicate count matters here more than anywhere
else:

| `num_replicates` | successes needed for `power_ci_low` ≥ 0.8 |
|---:|---:|
| 30 | 29/30 |
| 100 | 88/100 |
| 400 | 336/400 |

Measured on 34,886 pairs at effect size 0.15 and 100 replicates: 37.7 % of pairs are certified
(`power_ci_low` ≥ 0.8), 12.4 % are ambiguous, and 49.9 % are clearly underpowered
(`power_ci_high` < 0.8). That last half is not a precision problem — no number of replicates fixes
it; those pairs needed more perturbed cells or a larger effect.

The sentence the table supports, per pair, is then:

> `min_detectable_effect_size_ci_high` = 0.25 → we would have detected a knockdown of 25 % or more,
> so the absence of a call rules out effects that large, and says nothing about smaller ones.

### Per element — power is mostly a property of the gene, not the element

It is tempting to summarise power per `grna_target` and say "this element is well powered whatever
gene you pair it with". The data do not support that framing. Decomposing per-pair power at effect
size 0.15:

| | |
|---|---:|
| Variance between elements (ICC) | **20.9 %** |
| Variance within elements, gene to gene | 79.1 % |
| Mean within-element SD of power (elements with ≥ 5 pairs) | **0.34** |
| Correlation of element mean power with `mean_pert_cells` | 0.38 |

The gene matters about four times more than the element, and a single number per element hides a
±0.34 spread. A mean over the pairs an element happens to have tested is also not comparable between
elements: an element surrounded by highly expressed genes scores better than one surrounded by
lowly expressed genes at identical perturbation efficiency.

What *is* well defined is the conditional version. Per-pair power answers **"could we have detected
this link?"**; per-element power answers **"how well did we perturb this element?"** — a question
whose sufficient statistic is the perturbed-cell count, not anything about genes. So report it at a
reference gene: *"for a gene at median expression, element E has 0.7 power to detect a 15 %
knockdown."* That is comparable across elements. It also costs nothing extra: `power_summary.tsv`
carries both `mean_pert_cells` and `average_expression_all_cells` per pair, so fit
`power ~ f(mean_pert_cells, average_expression_all_cells)` and evaluate at a reference expression.
Simulating a reference-gene panel instead would cost ~450 CPU-hours, because the per-target
term in the cost model is paid whether a target carries one gene or thirty — see
[Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}#cost).

The element-level statement that does hold up is a floor: **21.9 % of elements (663 of 3,026) have
no tested pair that could have reached power 0.8 at all**, and those elements cannot support a
"regulates nothing" claim under any reading. In the other direction, only 5 of 3,026 elements have
*every* tested pair certified, so element-wide negative claims are essentially never assertable at
100 replicates and a 15 % knockdown. Worth knowing before building a figure around them.

## `power_es<effect_size>.tsv` — per effect size

| Column | Meaning |
|---|---|
| `grna_target`, `response_id` | The pair. |
| `power` | Fraction of replicates in which sceptre would have called the association. |
| `power_ci_low`, `power_ci_high` | 95% Wilson score interval. |
| `n_reps` | Replicates contributing. Can be below `--reps` if any replicate produced no fold-change estimate. |
| `mean_log_2_fold_change` | Mean simulated log₂ fold change across replicates. |
| `mean_pert_cells` | Mean perturbed cells. |
| `average_expression_all_cells` | Raw mean expression. |
| `effect_size` | The effect size simulated. |

### What `power` actually counts

A replicate counts as a success when **both** hold:

```
p_value < threshold   AND   log_2_fold_change < 0
```

The second condition makes this one-sided: a replicate only counts if the simulated perturbation
*reduced* expression. `threshold` is the largest nominal p-value that survived multiple-testing
correction in the real discovery analysis, so "detected" means the same thing here as it did in your
actual results.

### Reading the confidence intervals

Always read `power` together with its interval, especially at the boundaries. `power = 0` does not
mean the pair is undetectable — at 100 replicates the Wilson upper bound is 0.037, so the data are
consistent with a true power of a few percent. Likewise `power = 1` has a lower bound of 0.963.

The interval is Wilson rather than `p̂ ± 1.96·SE` precisely because the normal approximation returns
`[0, 0]` for 0 successes, asserting certainty the data do not support. See
[Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}).

## `per_replicate/es<effect_size>/*.tsv.gz` — per-replicate detail

One row per (pair, replicate), gzipped. This is the only output from which power can be
**re-derived** — subsampling replicates to study a reduced design, bootstrapping, or re-thresholding
— so it is worth keeping even though it dominates the output volume.

| Column | Meaning |
|---|---|
| `grna_target`, `response_id` | The pair. |
| `p_value` | From the resampling test on this replicate's simulated counts. |
| `log_2_fold_change` | Its effect estimate. Power counts a replicate only when `p_value` beats the discovery threshold **and** this is negative. |
| `rep` | Replicate index, unique across chunks thanks to `--rep-offset`. |
| `effect_size` | The effect size simulated. Constant within a file, and kept because `compute_power.R` uses it to refuse an input that mixes effect sizes. |
| `num_pert_cells` | Perturbed cells for this target. |
| `pass_qc`, `n_nonzero_trt`, `n_nonzero_cntrl` | Diagnostics, carried over from the **real** discovery pairs rather than recomputed from simulated data — they describe the observed experiment. |

**Five columns sceptre returns are dropped before writing**, because at 100 replicates × 34,886
pairs × 6 effect sizes they were 39 % of a 3 GB output and nothing read them:

| Dropped | Why |
|---|---|
| `fold_change` | It is `2^log_2_fold_change` — the same number twice. |
| `se_fold_change` | Read by nothing downstream. |
| `significant` | sceptre's own call at *its* threshold, not the discovery threshold this pipeline tests against. `compute_power.R` recomputes it, so keeping the column invited the wrong one being believed. |
| `average_expression_all_cells` | A per-*gene* constant that was repeated once per replicate. `summarize_power.R --sim-input` joins it from `sim_input.rds`, where it is stored once. |

Together with gzip that takes the per-replicate output from ~3 GB to a few hundred MB for a
six-point sweep at 100 replicates. Nothing needs a decompression step: `read.delim` sniffs the
magic number and handles `.tsv.gz` and `.tsv` alike, and `compute_power.R` accepts either.

## Intermediates

| File | Contents |
|---|---|
| `sim_input.rds` | Per-gene `mean`, `dispersion`, `average_expression_all_cells`; per-cell `size_factors` and categorical covariates; the gRNA and target perturbation matrices. No count matrix — the simulation draws counts rather than reading them. |
| `sceptre_template.rds` | The sceptre object with `@response_matrix` and `@grna_matrix` emptied. Neither is read by the discovery analysis once gRNAs are assigned. |
| `pairs.tsv` | `grna_target`, `response_id` for QC-passing pairs only. |
| `grna_targets.tsv` | `grna_id`, `grna_target`. |
| `discovery_threshold.txt` | A single number: the p-value a replicate must beat. |
| `analysis_mode.tsv` | `resampling_mechanism` (`crt` or `permutations`), `run_permutations`, and `moi`. Which test produced these power numbers. |
| `null_precomputations.rds` | Per-gene null models, one set per replicate, fitted on a null simulation. |

`split_*.tsv` is **not published**. The splits are parallelisation bookkeeping — 1,000 files per
sample — and they are regenerable: the bin packing is deterministic given `pairs.tsv` and
`n_splits`. They live in the Nextflow work directory and go when it is cleaned.

## Two means, deliberately

`sim_input.rds` carries two per-gene means and they are not interchangeable:

- **`mean`** — size-factor normalised. This is what the simulation draws counts from.
- **`average_expression_all_cells`** — raw. Reported in the output so power can be related to
  expression level.

Only the raw one appears in the output tables.
