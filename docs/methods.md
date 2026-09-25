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

There are two estimands, chosen with `--estimand` / `estimand` (default `fixed`). Both draw the same
guide knockdowns (next section). They differ in one step: whether the element's realised mean effect
over its perturbed cells is pinned to the requested one.

| | `fixed` (the default) | `random` |
|---|---|---|
| guide knockdown | Beta, mean es, sd `c * es * (1 - es)` | the same draw |
| element's realised mean over its perturbed cells | pinned to es in every simulation | left where the draw puts it |
| question answered | power for an element whose effect *is* es | power for an element whose effect is es *on average* |
| PerturbPlan setting that asks the same | `fold_change_sd = 0` | `fold_change_sd = c * es * (1 - es)` (0.0829 at es 0.15) |

Every output carries an `estimand` column: the per-simulation files, the per-effect-size power
tables and `power_summary.tsv`. A power table therefore always says which question it answers.
`consolidate_replicates.R`, `compute_power.R` and `summarize_power.R` refuse to mix the two.

**`fixed`: power at a fixed element effect.** "Power at 15%" is the probability of detecting the
pair when the element's guides produce **exactly** a 15% mean knockdown across the perturbed cells,
in every simulated screen. The guides still differ from one another (next section); what is fixed is
their cell-weighted mean. That is the effect sceptre's union test targets, since it pools every cell
carrying any of the element's guides; power is a property of a test at a value of the parameter it
targets. The match is close rather than exact: the test compares expected counts, which weight
each cell by its own baseline, and the guides' spread around the pinned mean adds a little variance
of its own. That spread is zero at effect size 0 (next section), so the null is exact.

**`random`: power at an element effect that is es on average.** The pin is skipped, and nothing else
changes. The realised mean over the perturbed cells is then the cell-weighted mean of the guides'
draws, `sum(n_g * x_g) / sum(n_g)`, where `n_g` is the number of perturbed cells carrying guide
g. It is unbiased for es, and its sd across simulations is

```
c * es * (1 - es) * sqrt(sum(n_g^2)) / sum(n_g)
```

the sd of a cell-weighted mean of independent guide draws. On the shared test fixture (seven guides
carried by 15-29 perturbed cells each) that is 0.38 times the per-guide spread, 0.032 at es 0.15.
Power under `random` is the expected power over that spread.

**Why `random` is safe now and was not before.** The unpinned spread is zero at es = 0. The old
absolute `N(1 - es, 0.13)` gave a null element random effects, and "detected" it about 20% of the
time on highly expressed genes. Under the Beta spread every guide effect at es = 0 is exactly 1, so
there is nothing to pin, and the null is identical under both estimands. The pin also consumes no
random numbers, so at any effect size the two estimands share every draw: the `fixed` effects are
the `random` ones shifted, gene by gene, onto es (and clamped at zero where a strong knockdown needs
it). One set of null-model fits serves both.

`fixed` was decided on 2026-09-24, and it is written down here because it was never written down
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
- Random guide effects were set aside on 2026-09-24 because, as the draws were then parameterised,
  they did not even start at α: at effect size 0 a target's guides still moved expression. That note
  asked for any uncertainty about guide efficiency to come as a separately labelled output with a
  multiplicative efficiency model. The Beta spread of 2026-09-25 is that model, and `--estimand
  random`, recorded in the `estimand` column, is that output (added 2026-09-25).

Two consequences to accept knowingly. WattEG does not reproduce the published DC-TAP per-pair power:
`fixed` is the other quantity, and `random` averages over the Beta spread rather than DC-TAP's
absolute `N(1 - es, 0.13)`. And under `fixed` the guides' spread barely moves power at effect sizes
up to 0.5, because spread around a pinned mean adds little variance to what the test sees: on moi5,
no pair's power moves by more than 0.02 between the spread described below and no spread at all.
How far `random` moves power from `fixed` has not yet been measured on moi5.

## Guide-to-guide variability

Guides targeting the same element are not equally effective. Each of the target's gRNAs therefore
gets its own knockdown, drawn independently for every gene:

```
this target's guides:  knockdown ~ Beta(mean = es, sd = c * es * (1 - es)),  effect = 1 - knockdown
every other guide:     exactly 1
```

with `guide_spread_c` (c) = 0.65 by default. The spread is **zero at es = 0** (a null element's
guides do nothing), grows in proportion to es for small effects, and levels off as knockdown
saturates. Every draw lies strictly inside (0, 1), so nothing needs clamping. c = 0 means no spread.

This form and value were settled on 2026-09-25 from per-guide data, replacing an absolute
`N(1 - es, 0.13)` that the original code carried:

- **The spread vanishes at no effect.** In three CRISPRi screens with ~15 designed guides per element
  (DC-TAP K562, DC-TAP WTC11, and moi5, the screen WattEG simulates), the between-guide spread of
  null pairs is 0.004-0.014 after removing sampling noise, the same floor as groups of
  non-targeting guides. Fulco 2019's FlowFISH null elements say the same. The absolute 0.13 kept a
  13% spread there. Around the pinned mean that inflated the false-call rate of highly expressed
  genes: 3,032 moi5 trans pairs sat above twice the nominal rate at the discovery threshold, the
  worst at 41 times.
- **It grows with the effect, then levels off.** Enhancer spread is 0.015 / 0.056 / 0.105 / 0.144 /
  0.108 at es 0.03 / 0.10 / 0.21 / 0.37 / 0.59. `c * es * (1 - es)` fits those bins at chi-square
  6.6 on 4 df, against 128 for the absolute form and 69 for a constant coefficient of variation.
- **c = 0.65 is fitted to those data** (0.645-0.648 in each screen separately). It reproduces the old
  0.13 at es ~0.28. The 0.13 came from one 2022 notebook that pooled the within-element spread of
  significant Fulco 2019 FlowFISH elements, without separating measurement noise and without asking
  whether the spread depends on es. Anchoring the new form to 0.13 at es 0.15 instead (c = 1.02)
  would change no moi5 pair's power by more than 0.05.
- It is the guide-efficacy form the CRISPR screen literature uses (JACKS, MAGeCK-MLE, CRISPhieRmix,
  Horlbeck 2016), where a guide's efficacy scales the element's effect.

The parameter is `--guide-spread-c` / `guide_spread_c`. The old `--guide-sd` / `guide_sd` is
refused, because reading an old 0.13 as c would shrink the spread fivefold.

**Other guides have no effect on the tested genes.** Control cells carry guides that belong to
other elements, or none at all, and those do not move this gene. The dispersion the counts are drawn
with was fitted to real cells that already carry their real guides, so it already contains whatever
those guides do. Every version until 2026-09-24 drew an extra `N(1, 0.13)` for them, going back
to the original DC-TAP code. That counted the noise twice: a gene with theta 146 came back from its
own simulated null data with theta 42, and power was understated for highly expressed,
low-dispersion genes.

**Under `fixed`, the perturbed cells' mean is then pinned** to the requested relative expression,
gene by gene: the perturbed block is shifted so its row mean equals `1 - effect_size` exactly. This
step is what makes the effect *fixed* (previous section), and `random` skips it; the check that
every control cell is exactly 1 runs under both. It is not a correction for clamping, which is what
this page used to say. Because the draw's mean is already `1 - es`, the pin only removes one
simulation's sampling wobble. At strong knockdowns a plain shift could push some guides below zero,
so the shift is solved for exactly instead: the result is `max(v + c, 0)` for the one constant `c`
that puts the mean on the target, which always exists for an effect size below 1. (An earlier
version shifted, clamped and repeated; once most cells clamp that converges slowly, and at effect
size ≥ 0.99 it ran out of iterations and stopped the run on a pin that exists.)

A cell carrying no guide gets multiplier 1. A cell carrying several of the target's guides has one
picked at random **once per target and effect size**, and it keeps that guide in every replicate.
(This page used to say "per replicate", which no implementation ever did. The pick is seeded by the
effect size as well as the target, so it can differ between effect sizes.) Only 0.18% of perturbed
cells carry more than one of their target's guides on moi5, so the choice barely matters.

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
