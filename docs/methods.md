---
title: Methods
nav_order: 5
---

# Methods

## Effect size parameterisation

Effect sizes are given as **fractional decreases in expression**: `0.15` means a 15% knockdown.
Internally this becomes a *relative expression level* of `1 - effect_size`, which multiplies each
gene's mean in the perturbed cells. So `0.15` scales expression to 0.85× baseline.

This conversion is the reason the column suffixes read `power_at_effect_size_15` — the label is
`effect_size × 100`.

## What simulated power means

**Power at a fixed element effect.** "Power at 15%" is the probability of detecting the pair when
the element's guides produce **exactly** a 15% mean knockdown across the perturbed cells, in every
simulated screen. The guides still differ from one another (next section); what is fixed is their
cell-weighted mean, which is the quantity sceptre's union test measures. Power is a property of a
test at a value of the parameter it targets, and this is that parameter.

This was decided on 2026-09-24, and it is written down here because it was never written down
before, which cost several regenerations:

- The original code (Sceptre_Power_Simulations, 2024, and the DC-TAP paper pipeline) called a
  centring step whose comment says it pins the mean. But it applied a cell-order mask to a matrix
  whose columns were still in perturbed-then-control order. It shifted the wrong columns, and the
  realised mean was left free to vary. On a real screen that is indistinguishable from not
  centring at all, so the published DC-TAP power averages over the realised knockdown instead
  ("random guide effects"). It got there by accident, not by design.
- WattEG inherited the same order. On 2026-09-21 (`2b76284`) the reorder was moved before the
  centring. That fixed a real indexing bug and, with it, silently changed which quantity was
  simulated.
- The alternative, random guide effects, is *expected* power averaged over a prior on the
  knockdown. It is a different quantity, and as the draws are parameterised it does not even start
  at α: at effect size 0 a target's guides still move expression, so a null element on a highly
  expressed gene is "detected" around 20% of the time. If uncertainty about guide efficiency is
  wanted, it belongs in a separately labelled output with a multiplicative efficiency model.

Two consequences to accept knowingly. WattEG does not reproduce the published DC-TAP per-pair power,
which is the other quantity. And `guide_sd` barely moves power at effect sizes up to 0.5, because
the guides' spread around a pinned mean adds little variance to what the test sees.

## Guide-to-guide variability

Guides targeting the same element are not equally effective. Each of the target's gRNAs therefore
gets its own effect size, drawn independently for every gene:

```
this target's guides:  N(1 - effect_size, guide_sd)
every other guide:     exactly 1
```

with `guide_sd = 0.13` by default. Negative draws are clamped to 0, because a guide cannot produce
negative expression.

**Other guides have no effect on the tested genes.** Control cells carry guides that belong to
other elements, or none at all, and those do not move this gene. The dispersion the counts are drawn
with was fitted to real cells that already carry their real guides, so it already contains whatever
those guides do. Every version until 2026-09-24 drew an extra `N(1, guide_sd)` for them, going back
to the original DC-TAP code. That counted the noise twice: a gene with theta 146 came back from its
own simulated null data with theta 42, and power was understated for highly expressed,
low-dispersion genes.

**The perturbed cells' mean is then pinned** to the requested relative expression, gene by gene: the
perturbed block is shifted so its row mean equals `1 - effect_size` exactly. This step is what makes
the effect *fixed* (previous section). It is not a correction for clamping, which is what this page
used to say: at effect size 0.15 a guide clamps with probability about 3 × 10⁻¹¹. At strong
knockdowns (effect size ≥ 0.7) the shift can push some guides below zero. They are clamped and the
shift is repeated until the mean is exact; if it cannot be reached, the run stops instead of
returning an effect that misses the target.

A cell carrying no guide gets multiplier 1. A cell carrying several of the target's guides has one
picked at random **once per target**, and it keeps that guide in every replicate: which guide a cell
carries is a fact about the screen, not about a draw. (This page used to say "per replicate", which
no implementation ever did.) Only 0.18% of perturbed cells carry more than one of their target's
guides on moi5, so the choice barely matters.

## Simulating counts

For each simulation, for gene *i* and cell *j*:

```
mu[i, j] = exp(X[j, ] · beta[i, ]) * effect_size[i, j]
count[i, j] ~ NegBinomial(mu = mu[i, j], size = 1 / dispersion[i])
```

where

- **`exp(X[j, ] · beta[i, ])`** is the expected count sceptre's own null model gives cell *j* for
  gene *i*: `X` is the cell's covariate row (library size, detected genes, batch, ...) and `beta`
  the gene's fitted coefficients. The simulation and the test that judges it are on one scale by
  construction. It is the default, `--expression-model fitted`, since 2026-09-21.
- **`dispersion[i]`** is `1 / theta` from sceptre's own negative-binomial fit
  (`@response_precomputations`), so the simulation inherits sceptre's dispersion estimates rather
  than refitting them.

`--expression-model size_factor` keeps the older `mean[i] * size_factor[j]` (a size-factor-normalised
mean times a DESeq2-style *poscounts* size factor). It exists only to reproduce sweeps made with it:
it mixes two models, runs about 4% low, and reproduces 86.5% of the observed count variance against
the fitted model's 99.5%.

Each cell keeps **its own** covariates and, under `size_factor`, its own size factor. The original
code (DC-TAP, from 2025-04-11) shuffled the size factors across cells in every replicate, on the
reasoning that simulated library sizes should be a draw from the observed distribution rather than
tied to each cell's identity. That is incorrect here: `effect_size[i, j]` is indexed by cell, so
shuffling pairs one cell's perturbation status with a different cell's library size, and the
covariates the test adjusts for with yet another cell's. Before 2025-04-11 it applied no size
factors at all. WattEG removed the shuffle in `14c28b6`.

Genes with no cached precomputation are a hard error rather than a silent skip. The previous
implementation stored dispersions in a list column with `NULL` holes; `unlist()` dropped them,
shortening the vector, and the negative-binomial draw then recycled it — so every gene after the
first gap would have been simulated with another gene's dispersion, with no warning.

## Deciding whether a simulation "detects" the pair

The simulated counts are handed to sceptre's `run_discovery_analysis()` with the pair table narrowed
to the target under test, and a simulation counts as a detection when

```
p_value < threshold   AND   log_2_fold_change < 0
```

`threshold` is the **largest nominal p-value that survived multiple-testing correction in the real
discovery analysis**, read from `@discovery_result`. Using the empirical threshold rather than a bare
`alpha` matters: it encodes the correction actually applied to your data, at your number of tests.
`--alpha` exists only for objects that have no discovery results.

At effect size 0 (the null arm) the same rule measures the rate of false calls **in the knockdown
direction**. That is α for a left-sided test, but only about α/2 for a two-sided one, because half
of the two-sided false calls have a positive fold change.

## Why control-cell sampling is not used

Using all non-perturbed cells as controls is expensive — a typical target has a few hundred
perturbed cells against several hundred thousand controls, and every simulation simulates all of
them. Sampling controls is the obvious optimisation, and it does not work.

Measured on the reference dataset (3 targets spanning 180/525/1,227 perturbed cells, 33 pairs, 30
simulations, 15% effect size, paired by seed against the all-controls baseline):

| `n_control_cells` | per-simulation cost | mean power | power retained | pairs lower / higher | sign test |
|---|---|---|---|---|---|
| 1,000 | 0.37s | 0.137 | 39.5% | 25 / 1 | p < 0.0001 |
| 2,000 | 0.42s | 0.197 | 56.7% | 23 / 0 | p < 0.0001 |
| 5,000 | 0.47s | 0.245 | 70.6% | 21 / 1 | p < 0.0001 |
| 20,000 | 0.63s | 0.273 | 78.5% | 19 / 2 | p = 0.0002 |
| all (~430,000) | 3.44s | 0.347 | 100% | — | — |

The intuition that controls stop mattering once they greatly outnumber the treated cells — because
`sqrt(1/n_trt + 1/n_ctrl)` is dominated by the treated term — applies to a two-group comparison.
sceptre's discovery analysis is a **conditional randomisation test**: the cell count sets the
resolution of the resampled null tail. With a threshold near 8 × 10⁻⁴, a few thousand cells cannot
reliably produce p-values that small, so genuinely detectable pairs fail to clear the bar.

The speedup is also smaller than the reduction in matrix size suggests: per-simulation cost is
sub-linear in cell count (20× more controls costs only 1.7× more time), because a fixed ~0.3s per
simulation is independent of it. The trade was roughly 7× speed for a 29% power loss at 5,000
controls.

`--n-control-cells` remains available for anyone who wants to validate it on their own data, and is
off by default.

## Monotonicity

Power must not decrease as the effect size increases. On the reference dataset (33 pairs, 12
simulations, effect sizes 0.15 / 0.25 / 0.5) mean power was 0.356 → 0.604 → 0.838, with a single
per-pair decrease of 0.08 → 0.00 that a two-proportion test could not distinguish from
Monte-Carlo noise (p = 1.0 at 12 simulations). This is a useful sanity check on any new dataset: a
systematic violation, or one that survives at high simulation counts, indicates a problem rather than
noise.

## Reproducibility

Seeds are derived from `(seed, target, replicate, effect_size)` via a portable string hash, not set
once per task. This makes results invariant to `n_splits` and `reps_per_chunk`, which are purely
computational parameters. Seeding per task would have meant that changing how the work was divided
silently changed the reported power — and the version of this pipeline before the refactor called
`set.seed()` nowhere at all, so no run could be reproduced.

Per-target random work (control selection, and picking one gRNA per cell) is seeded with the
simulation index 0, reserved for setup, so it too is independent of task layout.
