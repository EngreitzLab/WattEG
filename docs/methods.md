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
cell-weighted mean. That is the effect sceptre's union test targets, since it pools every cell
carrying any of the element's guides; power is a property of a test at a value of the parameter it
targets. The match is close rather than exact: the test compares expected counts, which weight
each cell by its own baseline, and the guides' spread around the pinned mean adds a little variance
of its own. That spread is zero at effect size 0 (next section), so the null is exact.

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
- The alternative, random guide effects, is *expected* power averaged over the element's realised
  knockdown. It is a different quantity, and with the old absolute spread it did not even start at
  α: at effect size 0 a target's guides still moved expression, so a null element on a highly
  expressed gene was "detected" around 20% of the time.

Two consequences to accept knowingly. The fixed estimand does not reproduce the published DC-TAP
per-pair power, which is the other quantity. And the guides' spread barely moves power at effect
sizes up to 0.5, because spread around a pinned mean adds little variance to what the test sees: on
moi5, no pair's power moves by more than 0.02 between the spread described below and no spread at
all.

### The random estimand, as an option

`--estimand fixed | random` (Nextflow `estimand`), default `fixed`, added 2026-09-25. It became
safe to offer once the guide spread vanished at effect size 0 (next section): the null is now the
same simulation under both.

| | `fixed` (default) | `random` |
|---|---|---|
| each guide's knockdown | Beta, mean es, sd `c * es * (1 - es)` | the same draw |
| the element's realised mean over its perturbed cells | pinned to es in every simulation | left where the draws put it |
| question answered | power for an element whose effect *is* es | power for an element whose effect is es *on average* |
| PerturbPlan setting that asks the same | `fold_change_sd = 0` | `fold_change_sd = c * es * (1 - es)` (0.0829 at es 0.15) |

Under `random` the realised mean varies between simulations with sd
`c * es * (1 - es) * sqrt(sum n_g^2) / sum n_g` over the guides' perturbed-cell counts `n_g`: the
variance of a cell-weighted mean. Both estimands draw the same guide effects from the same stream;
`random` only skips the pin, and control cells are exactly 1 under both. The estimand is written
into every per-simulation row, the power table and the summary, and the power steps refuse input
that mixes them. It is implemented in the Python path only; the R reference has not had it added.

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

**Under the fixed estimand, the perturbed cells' mean is then pinned** to the requested relative
expression, gene by gene: the perturbed block is shifted so its row mean equals `1 - effect_size`
exactly. This step is what makes the effect *fixed* (previous section). It is not a correction for
clamping, which is what this page used to say. Because the draw's mean is already `1 - es`, the pin
only removes one simulation's sampling wobble. At strong knockdowns a plain shift could push some
guides below zero, so the shift is solved for exactly instead: the result is `max(v + c, 0)` for the one constant `c` that puts the
mean on the target, which always exists for an effect size below 1. (An earlier version shifted,
clamped and repeated; once most cells clamp that converges slowly, and at effect size ≥ 0.99 it ran
out of iterations and stopped the run on a pin that exists.)

A cell carrying no guide gets multiplier 1. A cell carrying several of the target's guides has one
picked at random **once per target and effect size**, and it keeps that guide in every replicate.
(This page used to say "per replicate", which no implementation ever did. The pick is seeded by the
effect size as well as the target, so it can differ between effect sizes.) Only 0.18% of perturbed
cells carry more than one of their target's guides on moi5, so the choice barely matters.

## Simulating counts

For each replicate, for gene *i* and cell *j*:

```
mu[i, j]    = baseline[i, j] * effect_size[i, j]
baseline[i, j] = exp(X[j] . beta[i])
count[i, j] ~ NegBinomial(mu = mu[i, j], size = theta[i])
```

where `X` is the screen's covariate matrix and `beta[i]`, `theta[i]` come from **one fit** — the
Poisson GLM plus negative-binomial theta that sceptre's own `perform_response_precomputation`
performs, and that the null model of every discovery test is built from. Reproduced with pysceptre,
those agree with sceptre's cached values to about ten significant digits: theta to 2.8e-12 at the
median, and `exp(X·β)` to 9.2e-10 across 134.5 million gene × cell values.

**The point is that the baseline is the test's own model.** A power analysis asks whether the
screen's test would have detected an effect, so the counts it is shown should be counts from the
model that test assumes. Taking both the level and the noise from one fit is what makes that true
by construction rather than by coincidence.

**Each cell keeps its own covariates**, and so its own library size. The original code (DC-TAP,
from 2025-04-11) shuffled the size factors across cells in every replicate, on the reasoning that
simulated library sizes should be a draw from the observed distribution. That is incorrect here:
the effect size is indexed by cell, so shuffling pairs one cell's perturbation status with another
cell's library size, and the covariates the test adjusts for with yet another cell's. Before
2025-04-11 it applied no size factors at all. WattEG removed the shuffle in `14c28b6`. The fitted
baseline cannot reintroduce it: each cell's expected count is computed from that cell's own
covariate row.

### What this replaced, and why

Until 2026-09-21 the baseline was `mean[i] * size_factor[j]`: a size-factor-normalised gene mean
scaled by the cell's DESeq2 *poscounts* factor. Every sweep run before then used it, and
`--expression-model size_factor` still reproduces it. Three things were wrong with it, in
increasing order of weight.

**It mixed two models.** The dispersion came from sceptre's negative binomial and the expression
level from a DESeq2 normalisation — one simulated gene, two statistical models of the same data.

**It got the level wrong**, and the error survived the size factor. `mean[i]` sits 16 % below the
mean sceptre's model implies; multiplying by the cell's factor recovers most of that and leaves the
simulated genes about 4 % low on day0. The residual is a dropped covariance term: `mean[i]` is a
*mean of ratios*, and `E[x·sf] = E[x]E[sf] + Cov(x, sf)`.

**It got the shape wrong, which the level hides.** sceptre's mean varies with every covariate —
library size, detected genes, batch, replicate — while `mean[i] * size_factor[j]` varies with a
single scalar per cell. Measured against the real counts on day0, over 60 genes and 567,690 cells:

| | `mean[i] * size_factor[j]` | `exp(X·β)` |
|---|---:|---:|
| zero fraction, mean absolute error | 0.0053 | **0.0007** |
| predicted variance / observed | 0.865 | **0.995** |

Note the variance error pulls the **opposite** way from the level error — less variance inflates
power where less expression deflates it — so which way the change moves simulated power is not
obvious and has not been measured.

### Two guards that are not tidiness

A gene with no fitted model is a hard error rather than a silent skip. The pre-refactor
implementation stored dispersions in a list column with `NULL` holes; `unlist()` dropped them,
shortening the vector, and the negative-binomial draw recycled it — so every gene after the first
gap was simulated with another gene's dispersion, with no warning.

A theta clamped to the estimator's bounds `[0.01, 1000]` is **kept at the bound, with a
warning**, as sceptre keeps it. It used to be refused. That stopped the moi5 cis sweep on one gene
(mean 0.0085 counts per cell, theta at 0.01), and it would have dropped that gene's pairs while
the R implementation kept them. sceptre's own test uses the same clamped value in its null model,
so simulating from it keeps the simulation on the test's model. A gene that sparse has essentially
no power at any theta.

## Deciding whether a simulation "detects" the pair

The simulated counts are tested with the screen's own test, target by target (next section),
and a simulation counts as a detection when

```
p_value < threshold   AND   log_2_fold_change < 0
```

`threshold` is the **largest nominal p-value that survived multiple-testing correction in the real
discovery analysis**, read from `@discovery_result`. Using the empirical threshold rather than a bare
`alpha` matters: it encodes the correction actually applied to your data, at your number of tests.
`--alpha` exists only for objects that have no discovery results.

At effect size 0 (no simulated effect) the same rule measures the rate of false calls **in the knockdown
direction**. That is α for a left-sided test, but only about α/2 for a two-sided one, because half
of the two-sided false calls have a positive fold change.

## Running the screen's test on each simulation

Each simulation runs the screen's actual test on its actual cells, covariates, guide assignment and
threshold, including sceptre's permutation test with its escalation through three stages (pysceptre,
with the screen's own `B1`/`B2`/`B3` and side). The four settings below decide how that is laid out.
Each was measured before it became the default (2026-09-25; `docs/pysceptre-backend.md`, section 13).
Three of them change no result at all; the last is an approximation, and says so.

### One permutation set per target (`--permutations per-target`)

All simulations of a target are tested against one permutation set, drawn from a stream keyed on
(seed, target, effect size) and never on the simulation. That is what sceptre itself does: its
sampler reseeds `mt19937(4)` on every call (pinned commit 3ba046b), so the R implementation's
simulations already shared one fixed set, as the real screen's targets did. A per-target draw from
the run's seed was chosen over sceptre's fixed set so that targets of the same size do not all share
one draw. The simulated counts are the same either way; the p-values move by the Monte Carlo error
of the permutation test. `--permutations per-replicate`, one fresh set per simulation, was the
behaviour until 2026-09-25 and remains available with the engine.

### How the permutation nulls are computed (`--nulls sparse`)

pysceptre has two routes to a stage's null statistics. Its default for permutations, a prefix scan,
gathers a (B, n_trt, 14) array and takes a running sum over it, so that one set of draws can serve
many targets of different sizes. A simulation's call holds one target, so the running sum is thrown
away but for one column (227 MB at stage 2 for a 406-cell target). The sparse route multiplies a
sparse indicator matrix by the gene's pieces; pysceptre already takes it whenever the scan would be
too large. WattEG selects it through pysceptre's own documented signal, without changing
pysceptre. The p-values, z-statistics and stages were identical on 3,750 pair-tests, and every
output file is byte-identical; a stage-2 null went from 119.5 to 17.7 ms.

### Testing a target's simulations together (`--driver fast`)

The engine calls pysceptre once per (target, simulation), and each call redraws the same
permutations, rebuilds the same matrices and outer products, and computes each gene's pieces twice.
The fast driver does that work once per target and runs pysceptre's own per-gene steps
(`run_low_level_test_full`, its draw-matrix route, its fits) for every simulation. Its output is the
engine's under `per-target` and `sparse`, byte for byte, which the tests check against the engine on
every run.

### Null-model fits

sceptre's test needs each gene's null model: the Poisson GLM on the covariates plus the
negative-binomial theta, fitted to all cells. In the real screen each gene is fitted once and serves
the ~135 targets it is paired with. A simulation can do the same or refit, and the two are options:

- **`--null-fits reuse`**, the default (R's `FIT_NULL_MODELS` approximation). Once per (gene, simulation), the
  gene's counts are drawn with no knockdown, from their own stream keyed
  `(seed, "__null_fit__|" + gene, simulation, 0)`, and pysceptre's own gene fit is run on them. Every
  target the gene is tested with in that simulation, at every effect size, uses that fit; the score
  and the permutations still come from the simulation's own counts. `FIT_NULL_MODELS` makes the fits
  once per sample; a task given no file makes the same fits for its own genes.
- **`--null-fits refit`**. Each simulation refits every gene on its own counts, once per target,
  as sceptre would if the simulated screen were the real one.

Reuse is not exact, and it does not need to be: a target perturbs a few hundred of 131,055 cells, so
its knockdown barely moves a fit made on all of them. Measured on the moi5 reference targets with
the same counts and permutations (cis: 12 pairs x 100 simulations; trans: 255 pairs x 10):

| | refit | reuse | flips | McNemar p | largest per-pair change in power |
|---|---:|---:|---:|---:|---:|
| cis calls | 808 | 805 | 5 / 2 | 0.45 | 0.02 |
| trans calls | 1,280 | 1,282 | 3 / 5 | 0.73 | 0.10 (one simulation of ten) |

No pair moved by more than its own Wilson half-width, and p-values moved by a median of 0.04 in
log10 (an earlier benchmark put that at a quarter of what a change of permutation seed does). On the
fast driver the simulation's worker time fell 1.6-1.7x with the fits reused, the fit step excluded.

A fit depends on (seed, gene, simulation) and nothing else — each gene's baseline and draw are made
alone and pysceptre fits one gene at a time — so no result depends on which genes or targets share
a task, and a file made by `watteg-fit-null-models` gives the same bytes as fits made in the task.
The file records the seed, the baseline model and a digest of `sim_input.h5`, and a simulation run
that does not match all three refuses it.

## Why control-cell sampling is not used

**The flags are gone.** `n_control_cells` and `cell_batches` were removed with the Python port:
they bought speed at a measured 21–60 % of power, were off by default, and the port is the reason
speed stopped being the binding constraint. A knob that trades power for speed you no longer need
is a trap rather than an option. The measurement below is why, and is kept because the reasoning
would otherwise have to be rediscovered by whoever proposes the optimisation next.

**On dropping the batch stratification with it**, which is the part that sounds risky.
`cell_batches` made the control draw keep the perturbed cells' batch composition, and without a
draw there is nothing to stratify — the control group is every non-perturbed cell, so its
composition is the dataset's. Batch is also conditioned on twice by the test itself: it is a
covariate of the per-gene model, and of the logistic fit the CRT draws its synthetic treated sets
from, so the null distribution is conditional on batch by construction. Cell-level matching is what
you reach for when the model cannot adjust for a confounder.

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
