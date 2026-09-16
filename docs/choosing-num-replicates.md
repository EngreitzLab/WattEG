---
title: Choosing num_replicates
nav_order: 4
---

# Choosing `num_replicates`

`num_replicates` (`--reps` on `run_power_simulation.R`) is the number of Monte-Carlo simulation
simulations run per element–gene pair. It is the **only** pipeline parameter with statistical
content — `n_splits` and `reps_per_chunk` only decide how the work is divided across tasks.

This page explains what the number controls, how to pick it for your question, and how to read the
confidence intervals in the output.

---

## What is actually being estimated

For one pair and one effect size, each simulation simulates a fresh count matrix under that effect
size and asks sceptre whether it would have called the association. A simulation is a **success**
when

```
p_value < threshold   AND   log_2_fold_change < 0
```

The second condition makes this a one-sided criterion: a simulation only counts if the simulated
perturbation *reduced* expression, which is the direction a repressive element is expected to act
in. The reported power is the success fraction:

```
power = successes / n_reps
```

so it is a binomial proportion with `n_reps` trials, and every simulation is an independent draw.
`num_replicates` sets the precision of that proportion — nothing else. It does not make the
simulation more realistic, and it does not change the quantity being estimated.

Two consequences follow directly.

**Resolution.** The estimate can only take multiples of `1/n_reps`. At 20 simulations the only
achievable values are 0, 0.05, 0.10, …; a reported power of 0.85 is not representable at all.

**Monte-Carlo error.** The standard error of the estimate is `sqrt(p(1-p)/n_reps)`, largest when
the true power is near 0.5 and shrinking to zero at the extremes:

| `num_replicates` | SE at p=0.5 | 95 % CI width at p=0.5 | SE at p=0.8 | resolution |
|---:|---:|---:|---:|---:|
| 20 | 0.112 | ±0.219 | 0.089 | 0.050 |
| 50 | 0.071 | ±0.139 | 0.057 | 0.020 |
| 100 | 0.050 | ±0.098 | 0.040 | 0.010 |
| 200 | 0.035 | ±0.069 | 0.028 | 0.005 |
| 400 | 0.025 | ±0.049 | 0.020 | 0.003 |
| 1000 | 0.016 | ±0.031 | 0.013 | 0.001 |

---

## How many simulations do you need?

There is no single right answer, because it depends on what you do with the number. Find your case
below. If your question is the false-negative one — was a non-significant pair a real biological
negative, or did the assay simply lack the power to see it — skip to
[Telling a real negative from an underpowered one](#false-negatives), which is the demanding case
and the one this pipeline was built for.

### Aggregate summaries — 100 is generous, even 20 would do

If you report a mean power across pairs, or an expected number of discoveries, the Monte-Carlo
error averages down across pairs as `1/sqrt(num_replicates × n_pairs)`:

| pairs summarised | SE of mean power at `num_replicates` = 100 |
|---:|---:|
| 1 | 0.050 |
| 100 | 0.005 |
| 34,886 | **0.0003** |

For a dataset-level number over tens of thousands of pairs, simulation count is not your limiting
source of uncertainty — the effect-size assumption is.

### Per-pair power in a results table — 100 is reasonable

At `num_replicates = 100` a per-pair estimate carries a 95 % interval of roughly ±0.10 in the worst
case. That is honest and usually adequate, provided the interval is reported alongside it (it is —
see below). At 20 simulations the interval is ±0.22, which is too coarse to publish per pair.

### Telling a real negative from an underpowered one — threshold the lower bound, and use ~400
{: #false-negatives }

This is the primary use of the pipeline: a pair the discovery analysis did not call is either a
biological non-effect or a measurement failure, and power is what separates them. Note what that
claim actually requires. It is not "the point estimate is above 0.8" but "the power was *at least*
0.8", so the column to threshold is **`power_ci_low`**, not `power`.

That is a stricter test than it looks, and the strictness is what rules out small simulation counts:

| `num_replicates` | successes needed for `power_ci_low` ≥ 0.8 | implied `power` |
|---:|---:|---:|
| 30 | 29/30 | 0.967 |
| 100 | 88/100 | 0.880 |
| 200 | 172/200 | 0.860 |
| 400 | 336/400 | 0.840 |

At 30 simulations a pair must be all but perfect before its interval clears 0.8, so the criterion is
effectively unusable — a reason not to run 30 for this question that is separate from, and firmer
than, the precision argument. At 100 the criterion silently discards the 0.80–0.88 band; at 400 only
0.80–0.84.

Measured on all 34,886 pairs at effect size 0.15 and 100 simulations, the three outcomes are:

| | pairs | |
|---|---:|---:|
| measured power ≥ 0.8 (`power_ci_low`) — a negative here is biological | 13,154 | 37.7 % |
| ambiguous, the interval contains 0.8 | 4,322 | 12.4 % |
| clearly underpowered, `power_ci_high` < 0.8 — a negative here says nothing | 17,410 | 49.9 % |

Only the middle row can be moved by more simulations, which is what the next section exploits. The
bottom row is not a precision problem: those pairs need more perturbed cells or a bigger effect
size, and no number of simulations will rescue them.

### Two-stage allocation: 400-level precision at roughly half the cost

The variance `p(1-p)/n_reps` collapses at the extremes, and many pairs *are* at the extremes, so
simulations spent on a pair that came out 0/100 or 100/100 buy almost nothing. Target the second
stage at the pairs whose interval straddles the cutoff you care about:

1. **Stage 1** — `num_replicates = 100` for every pair.
2. **Stage 2** — 300 more simulations only for the 12.4 % of pairs whose stage-1 interval contains
   0.8, invoked with `--rep-offset 100` so the `rep` column stays unique and the two sets of
   results can simply be concatenated.

Measured cost on this dataset, per effect size:

| design | CPU-hours | precision at the cutoff |
|---|---:|---|
| uniform 100 | 858 | ±0.08 |
| **two-stage 100 → 400** | **1,768** | ±0.04 where it matters |
| uniform 400 | 3,430 | ±0.04 everywhere |

**Do not size stage 2 by counting replicates.** The 4,322 ambiguous pairs are only 12.4 % of the
pairs, but they are spread across 1,820 of the 3,026 targets, and cost is `4.91 s + 0.459 s × pairs`
per (target, simulation) — see [Cost](#cost). Stage 2 therefore pays 60 % of the per-target overhead
to re-simulate 12 % of the pairs, and **82 % of its 910 CPU-hours is that overhead**. The naive
"`100 + 0.124 × 300 ≈ 137` replicate-equivalents per pair" arithmetic predicts 1,175 CPU-hours and
understates the real figure by a third. Any scheme that filters pairs rather than whole targets has
this property.

An earlier version of this page put the uncertain band at 21 % of pairs and two-stage at 108
replicate-equivalents, measured on 33 pairs at 30 replicates. On the full 34,886 pairs at 100
simulations the band (0.1, 0.9) holds **44.1 %** of pairs, not 21 %, and 17.1 % land at exactly 0 or
1 rather than 55 %. Thirty simulations overstates how many pairs are truly at the extremes, because
at 30 draws a pair with true power 0.9 hits 30/30 about 4 % of the time. Band-based targeting on the
real numbers costs `30 + 0.441 × 370 ≈ 193` replicate-equivalents — which is why the design above
targets the cutoff instead of the band.

Pooling the two stages is unbiased — simulations are exchangeable draws, so combining them is fine
as long as `n_reps` is reported per pair, which the output does. The one thing to avoid is
comparing raw power values between pairs while ignoring that they rest on different simulation
counts; use the intervals for that.

To run it, invoke the second stage with `--rep-offset 30` on the subset of pairs that need it, so
the `rep` column stays unique and the two sets of results can simply be concatenated.

---

## Reading the confidence intervals

`compute_power.R` reports, per pair:

| column | meaning |
|---|---|
| `power` | success fraction, `successes / n_reps` |
| `n_reps` | simulations actually contributing (simulations with a missing `log_2_fold_change` are excluded, so this can be below `num_replicates`) |
| `power_ci_low`, `power_ci_high` | 95 % **Wilson** score interval |

The interval is Wilson rather than the textbook `p̂ ± 1.96·SE` for a specific reason: the normal
approximation breaks down exactly where these estimates live. A pair with 0 successes out of 100
gets the interval `[0, 0]` under the normal approximation — asserting with certainty that the power
is zero. The Wilson interval gives `[0, 0.037]`, which is what the data actually support: 100
simulations without a single success is consistent with a true power of a few percent.

So **do not read `power = 0` as "this pair cannot be detected"**. Read it together with
`power_ci_high`. The same applies at the top: `power = 1` at 100 simulations has a Wilson lower
bound of 0.963, not 1.

Worked examples at `num_replicates = 100`:

| successes | `power` | Wilson 95 % interval |
|---:|---:|---|
| 0 | 0.00 | [0.000, 0.037] |
| 5 | 0.05 | [0.022, 0.112] |
| 50 | 0.50 | [0.404, 0.596] |
| 80 | 0.80 | [0.711, 0.867] |
| 100 | 1.00 | [0.963, 1.000] |

Note how the 0.8 row overlaps a cutoff of 0.75 *and* is consistent with 0.86 — that is the
cutoff-blurring described above, made explicit.

---

## Cost

Simulation count is the parameter you pay for linearly. Measured on Sherlock over 36 targets
(3,026 targets and 34,886 pairs in the full run), cost fits

```
seconds per (target, replicate) = 4.91 + 0.459 × (pairs in that target)      R² = 0.968
```

which puts **48 % of the bill in a per-target term that pair counts do not touch**: one
`run_discovery_analysis()` call carries the full 586,309-cell bookkeeping however few gene pairs
ride along. Two targets can only share a call when their gene sets are disjoint — see Step 9 in
[Status]({{ site.baseurl }}{% link status.md %}).

| `num_replicates` | CPU-hours per effect size |
|---:|---:|
| 30 | ~257 |
| 100 | ~858 |
| 400 | ~3,430 |
| two-stage 100 → 400 | ~1,768 |

Simulations are the one lever that scales *both* terms, which is why halving them halves the bill
exactly, while dropping half the pairs saves only ~26 %. Against that, doubling simulations only
improves precision by `sqrt(2)` — hence spending them where the variance is rather than uniformly.

Effect sizes cost the same as each other: the effect size changes the numbers drawn, not the work
done. A six-point sweep is six times the figures above.

---

## Reproducibility

Simulation `k` of a given pair always receives the same random draws, because the seed is derived
from `(seed, target, rep, effect_size)` rather than from the task. Two practical consequences:

- Results are **invariant to `n_splits` and `reps_per_chunk`**. Re-running the same analysis with a
  different task layout reproduces the same numbers.
- Increasing `num_replicates` is **incremental**: simulations 1–100 keep exactly the values they had,
  and 101–400 are new. You can extend an existing run rather than redo it.

Both properties depend on passing the same `--seed`. It is required, not optional — the estimates
are stochastic, and an unseeded run cannot be reproduced or extended.

## See also

- `docs/methods.md` — how the effect size and guide-to-guide variability are parameterised
- `docs/output.md` — every output column
- `docs/usage.md` — `n_splits` / `reps_per_chunk` and SLURM sizing
