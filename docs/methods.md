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

## Guide-to-guide variability

Guides targeting the same element are not equally effective, and treating them as identical would
overstate power. Each gRNA therefore gets its own effect size, drawn as

```
targeting guides:  N(1 - effect_size, guide_sd)
control guides:    N(1,               guide_sd)
```

with `guide_sd = 0.13` by default (previously hardcoded). Negative draws are clamped to 0 — a guide
cannot produce negative expression.

Clamping biases the mean upward, so each gene's effect-size matrix is then **re-centred**: the
perturbed block is shifted so its row mean equals the requested relative expression level, and the
control block so its row mean equals 1. Without this step the realised effect size would be
systematically weaker than requested.

Cells carrying no guide get a no-effect row (multiplier 1). A cell carrying several guides has one
picked at random, per replicate.

## Simulating counts

For each simulation, for gene *i* and cell *j*:

```
mu[i, j] = mean[i] * size_factor[j] * effect_size[i, j]
count[i, j] ~ NegBinomial(mu = mu[i, j], size = 1 / dispersion[i])
```

where

- **`mean[i]`** is the size-factor-normalised mean expression of gene *i* in the real data,
- **`dispersion[i]`** is `1 / theta` from sceptre's own cached negative-binomial fit
  (`@response_precomputations`), so the simulation inherits sceptre's dispersion estimates rather
  than refitting them,
- **`size_factor[j]`** is a DESeq2-style *poscounts* size factor, computed directly on the sparse
  matrix (DESeq2 itself is not used — `DESeqDataSetFromMatrix()` densifies, which is untenable at
  these dimensions).

Each cell keeps **its own** size factor. The pre-refactor code shuffled the size factors across
cells before use, on the reasoning that simulated library sizes should be a draw from the observed
distribution rather than tied to each cell's identity. That is incorrect here: `effect_size[i, j]`
is indexed by cell, so shuffling pairs one cell's perturbation status with a different cell's
library size and breaks the correspondence the model assumes.

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
