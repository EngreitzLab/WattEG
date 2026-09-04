# Experiment plan

EM's plan, in order, with what each step costs and what could sink it.

**The headline finding while costing this out: steps 1–3 need no new compute.** The
`dc_tap_paper_wtc11_no_shuf` sweep is stored *per replicate* — a `rep` column, 657,400 rows for
effect size 0.15, i.e. 6,574 pairs × 100 replicates, across all six effect sizes
(0.05 / 0.1 / 0.15 / 0.2 / 0.25 / 0.5). So "run a full sweep" is already done for the dataset the
method was developed on, and reduced designs can be studied by **subsampling** it rather than by
re-running anything.

That leaves Gasperini as the only step needing real compute, and its cost is upstream of power.

---

## Step 1 — a full sweep as ground truth

**Status: exists for two datasets now.**

| | `wtc11` | `day0` |
|---|---:|---:|
| Pairs | 6,574 | **34,886** |
| Effect sizes × replicates | 6 × 100 | 6 × 100 |
| Null model | inherited cache (pre-refactor) | **`null_fit`** |
| Powered | well; transitions 5–15 % | poorly; 84 % transition 5–25 % |
| Monotonicity violations | 0.67 % | 1.32 %, largest 0.070 (≈ 1 SE) |

`day0` was swept on 2026-09-03 and lives at
`oak/projects/element-gene-power-analysis/power_sweep/day0/`. Its `power_summary.tsv` already *is*
the ground-truth object this step asks for: per-pair power at all six effect sizes with Wilson
intervals, `min_detectable_effect_size` with its `_ci_low`/`_ci_high` bracket, and the covariates
§6 needs. Nothing to build.

The two are complementary rather than duplicative — one well powered, one not — so a reduced design
that holds on both holds across the range that matters. day2 and day4 were started alongside and
stopped part-way; both are resumable and would take the count to four.

Nothing to run. What *is* needed is a defined ground-truth object: per-pair power at each of the six
effect sizes, with Wilson intervals, plus the minimum detectable effect size derived from
`power_ci_low`. Everything below is measured against this.

**Caveat to settle first:** this sweep was produced by the pre-refactor pipeline, so it carries the
inherited-null-model bias (`docs/status.md`, the settled null-model section) — power understated by
roughly 0.03, one-directionally. For steps 2–3, which compare *reduced designs against the full design
on the same data*, the bias cancels: both sides share it. For step 4, comparing predictions across
datasets, it does not obviously cancel and must be handled. **[DECISION]** either accept it for 2–3 and
state it, or rerun `wtc11` under `null_fit` first (~6 × 1,300 CPU-h, expensive).

---

## Step 2 — show a reduced sweep suffices

**The claim.** Three effect sizes at 30 replicates give per-pair curves as good as six at 100.

**Cost: free.** Subsample replicates and effect sizes from step 1 and refit. Minutes of R, no cluster.

**Do more than test 30.** Since it is free, scan the whole design space rather than asserting one point:

- replicates per point: 10 / 20 / 30 / 50 / 100
- number of points: 3 / 4 / 6
- grid placement: the four grids already tested, plus transition-bracketing grids

Report where the fit breaks rather than that one design works. The deliverable is a heatmap of held-out
error over (replicates × grid), which is a far stronger figure than a single comparison and answers the
reviewer question "why 30?" with a curve instead of a choice.

**A distinction that must be stated explicitly or the paper contradicts itself.** Elsewhere we argue 30
replicates *cannot* certify a negative, because certifying power ≥ 0.8 from a raw estimate needs 29/30
successes. That is not in tension with this step, but the reason has to be spelled out: certification
comes from the **fitted curve**, which pools 3 × 30 = 90 draws into one free parameter, not from any
single effect size's raw proportion. The mini-sweep is a claim about estimating a curve; the 88/100 rule
is a claim about a single binomial. Say so, or a reader will think we contradict ourselves.

**Subsampling is not identical to running fewer replicates**, but it is unbiased for this purpose since
replicates are i.i.d. draws. Subsample without replacement, repeat over many random subsets, and report
the spread — that also gives error bars the single-run version could not.

---

## Step 3 — test the predictions

Two distinct predictions, worth keeping separate because they can fail independently:

1. **Interpolation within a dataset.** Fit on the reduced grid, predict the held-out effect sizes,
   compare to measured. Already done once (MAE 0.0167 on the best grid); step 2 generalises it across
   designs. Free.
2. **Per-pair minimum detectable effect size.** The actual deliverable. Does the reduced design
   reproduce the full design's MDES per pair? Currently 89.0 % exact / 99.2 % within one grid step on
   the best grid. Note the ceiling: the six-point MDES is itself noisy (±0.05 per point), so some
   disagreement is the *reference* being wrong. Quantify that ceiling rather than reporting raw
   agreement — bootstrap the reference from the per-replicate data, which the storage format allows.

Free. Both should be done before any Gasperini compute is committed, because if they fail there is
nothing to transfer.

---

## Step 4 — transfer to Gasperini

**The design.** Take parameters fitted on `wtc11`, predict power for Gasperini pairs, then compare
against a real full sweep on Gasperini.

This is the paper's key experiment and the only expensive one. Three things to get right:

**a. Test two variants, not one.** Per `paper/literature.md`, PerturbPlan already derives power
analytically from perturbed-cell count, expression and dispersion — no fitting required. So compare:

| Variant | What it tests |
|---|---|
| Transfer `wtc11`-fitted coefficients to Gasperini | Does the *empirical* model generalise? |
| Apply the analytical formula to Gasperini directly | Is fitting needed at all? |
| Gasperini-fitted coefficients (in-sample ceiling) | How much is lost by transferring? |

If the analytical formula wins, the paper's §6 becomes "validate and calibrate an existing formula per
pair", which is a better paper than "our regression transfers". If transfer wins, that is a real and
surprising result. Either outcome is publishable; only running the transfer variant risks reporting a
worse method as though it were the only option.

**b. The cost is upstream, not in the power simulation.** Gasperini et al. must first become a sceptre
object: gRNA assignment, QC, and a discovery analysis, all before a single power replicate runs.
**[NEEDS ESTIMATE]** — gRNA assignment alone was ~24 CPU-hours per object in the sceptre-pipeline
evaluation. Budget this explicitly; it is easy to plan around the power sweep and be surprised by the
week of data engineering in front of it.

**c. A full six-point sweep on Gasperini, at full replicate count. Decided — do it.**

An earlier draft of this file suggested reusing the reduced design on Gasperini to save compute. That
was wrong, and EM's objection is the right one: the reduced design is *the thing being validated*, so
using it as its own reference is circular. A reviewer would reject it, and correctly — you cannot
demonstrate that three points reproduce six by only ever measuring three.

So the reference is a genuine full sweep: six effect sizes, full replicates, on Gasperini. It is the
expensive step and it only has to be paid **once**, because thereafter every dataset can use the
reduced design with the transfer already demonstrated. That is precisely the argument the paper is
making, and it cannot be made without one dataset where both designs were run in full.

The same logic applies to `wtc11`, where the full sweep already exists — so after Gasperini the claim
rests on two independent full sweeps, in different labs and protocols, which is a much stronger
position than one.

**[NEEDS ESTIMATE]** the cost. Scale from the measured `day0` model under `null_fit`
(`5.51s + 0.891s × pairs` per target-replicate, ~1,305 CPU-h per effect size at 3,026 targets /
34,886 pairs / 586,309 cells): six effect sizes there would be ~7,800 CPU-h. Gasperini differs in cell
count and pair count, and per-call cost is sub-linear in cells, so this needs recomputing from its
actual pair table rather than assumed. Note also that only **one** set of null-model fits is needed for
the whole sweep, since a null simulation is effect-size independent — so that cost does not multiply by
six.

---

## Step 5 — `trans` (not in EM's list, but the headline claim needs it)

Per `paper/literature.md`, `trans` is the one thing no existing tool addresses, and §6 of the outline
makes it the reason the work matters. It needs its own demonstration: simulate power directly for a
tractable subset of `trans` pairs and compare against prediction. Without it the `trans` claim is an
extrapolation from `cis` into a regime where nothing has been measured.

**[NEEDS DESIGN]** how to choose the subset, and whether `trans` pairs differ systematically — their
genes are not drawn from the element's neighbourhood, so the expression distribution differs by
construction, which is exactly the axis the prediction is most sensitive to.

---

## Ordering, and where to stop

```
1. define ground truth from the existing wtc11 sweep        free
2. design-space scan by subsampling                          free
3. interpolation + MDES agreement, with a noise ceiling      free
   -- gate: if 2-3 fail, stop; there is nothing to transfer
4a. Gasperini -> sceptre object                              [NEEDS ESTIMATE], the real cost
4b. reduced sweep on Gasperini + six-point reference subset
4c. three prediction variants compared
5. trans subset demonstration                                [NEEDS DESIGN]
```

Steps 1–3 are a week of analysis with no cluster time and would establish most of the method's claims.
They are also the gate on step 4. Doing them first is cheap insurance against committing the Gasperini
budget to a method that does not survive its own held-out test.
