# Figures

One entry per figure. For each panel: what it shows, what data it needs, what would produce it, and
whether that data exists today. Nothing here is plotted yet.

`docs/status.md` holds the measurements; this file holds only what a panel would need.

**Legend for state:** ✅ data exists · ⏳ data being generated · ❌ needs new compute · 💭 not designed yet

---

## Figure 1 — the problem

| Panel | Content | Data | State |
|---|---|---|---|
| a | Per-pair power at a 15 % knockdown, partitioned certified / ambiguous / clearly underpowered (37.7 / 12.4 / 49.9 %) | `power_es0.15.tsv` from the 100-replicate array | ⏳ array at 991/1000 |
| b | The same pairs' power against gene expression — power is not random with respect to biology | same, plus `average_expression_all_cells` | ⏳ |
| c | Schematic: identical screen output, two underlying truths | none | 💭 |
| d? | Fraction of negatives in a published training set that could not have reached power 0.8 | needs a scE2G / ENCODE-rE2G training set joined to pair-level power | ❌ |

Panel (d) is the one that would make the motivation concrete rather than assumed. It needs an actual
published training set, which we do not currently have joined to these pairs.

---

## Figure 2 — method schematic

Flow: sceptre object → impose knockdown on the element's perturbed cells → simulate counts → re-run the
real test → replicate → power with a Wilson interval → minimum detectable effect size.

Annotate where the null model enters, since Figure 8 shows that choice is not innocuous.

State: 💭 — no data needed, but should be drawn only once §3 is stable.

---

## Figure 3 — the derived power curve

| Panel | Content | Data | State |
|---|---|---|---|
| a | Probit-scale linearity: example pairs, six measured points, three-point fit overlaid | `dc_tap_paper_wtc11_no_shuf` six-point sweep + `bin/fit_power_curve.R` | ✅ |
| b | Held-out predicted vs measured power, all pairs, with the binomial noise floor as a reference band | same | ✅ |
| c | Held-out MAE by fit grid (the four grids tested: best 0.05/0.1/0.25 at MAE 0.0167) | same | ✅ |
| d | `deviance/df` distribution, marking where the straight line fails at the saturated ends | same | ✅ |

**The whole figure rests on one dataset.** A second panel row on `k562` or Gasperini is what makes it a
methods result rather than an observation. ❌

---

## Figure 4 — choosing the three effect sizes

| Panel | Content | Data | State |
|---|---|---|---|
| a | Transition-location distribution for two datasets, showing they differ (`wtc11` 5–15 %, `day0` mostly 5–25 %) | `wtc11` sweep + `day0` at 0.15 projected | ✅ partly — `day0` is projected from a single effect size, not measured |
| b | Fraction of pairs with ≥1 informative point (power 0.1–0.9) by grid — 26.3 % to 100 % | same | ✅ |
| c | Pilot-chosen grid vs hindsight-optimal grid, held-out MAE | needs the pilot recipe run end to end | ❌ |

Panel (c) is the one that turns §5 from a reasoned recipe into a demonstrated one.

---

## Figure 5 — predicting power from covariates

| Panel | Content | Data | State |
|---|---|---|---|
| a | Fitted `k` vs predicted `k`, `sqrt(n)` exponent annotated (0.502 vs theory 0.500) | 6,030 `wtc11` pairs with a usable fit | ✅ |
| b | Partial dependence on perturbed cells and on expression | same | ✅ |
| c | Residuals by expression decile / dispersion / cell count — the non-randomness check | same | ✅ **computed on day0**: width varies **8.61×** across expression deciles (1.16× across cell count), a clean inverted-U. Errors are concentrated, which is the case that harms a classifier |
| d | Calibrated prediction intervals vs measured power, held out | same + a calibration split | ✅ **computed on day0**: expression-conditional 90 % intervals cover 89.4–91.3 % held out, and certify **86–96 %** of what measurement certifies. Global calibration covers on average but is wrong per stratum (99.7 % vs 79.8 %) |
| e | Coefficients fitted per dataset, side by side — do they transfer? | needs ≥2 more sweeps | ❌ |

Panel (e) tests cross-dataset transfer, which EM has since scoped as reviewer-facing rather than
load-bearing: the goals are reducing computation and predicting `trans` and unmeasured pairs *within*
a dataset. Panel (c) turned out to matter more than expected — the error is strongly non-random in
expression, which is exactly what harms a classifier, and it is what forced the conditional
calibration in (d).

---

## Figure 6 — the `trans` demonstration

Predicted vs directly simulated power for a tractable subset of `trans` pairs, with calibrated
intervals.

Needs: a `trans` pair table, a simulation run over a subset of it, and the covariate model applied to
the same pairs. ❌ — nothing exists yet.

This is the panel that turns the headline claim from argument into evidence. If it cannot be produced,
§6 should be demoted to discussion.

---

## Figure 7 — aggregation: pair, element, screen

| Panel | Content | Data | State |
|---|---|---|---|
| a | Variance decomposition, between- vs within-element (20.9 % / ~79 %) | 100-replicate array output | ⏳ |
| b | Within-element power spread for example elements, ordered by mean (mean within-element SD 0.34) | same | ⏳ |
| c | Element power at a reference gene vs the naive mean over tested pairs — the naive version is confounded by which genes sit nearby | same + the §6 model evaluated at reference expression | ⏳ + method not implemented |

---

## Figure 8 (supplementary) — the null-model trap

| Panel | Content | Data | State |
|---|---|---|---|
| a | p-value agreement in aggregate: `as_is` vs `cleared`, Spearman 0.9990, median abs diff 4.6e-05 — looks fine | `cache_experiment.rds`, job 38854980 | ✅ |
| b | The same data at the call level: 7 of 265 flips, all one-directional | `threshold_check.R`, job 38861277 | ✅ |
| c | Cost of the three configurations, and the equivalence of `null_fit` to `cleared` (0 flips) | `cache_experiment_timings.csv` | ✅ |

Currently 3 targets, 53 pairs, 5 replicates. The paired subset comparison (steps 08–09, 20 splits at
100 replicates) would replace this with a far stronger version. ⏳

The pedagogical point is (a) next to (b): aggregate correlation hides the thing that matters. Worth
keeping even if the numbers get stronger.

---

## Figure 9 (supplementary) — cost

Measured runtime per call against pairs per call, for the three null-model configurations, plus implied
CPU-hours per effect size.

Data: `cache_experiment_timings.csv` ✅ for the ratios; absolute numbers should come from the array's
`sacct` record ⏳, since the experiment's clock ran ~1.9x slow.

---

## Figures we might want but have not designed

- Sensitivity of the minimum detectable effect size to replicate count — the argument for why 30
  replicates cannot certify a negative. Cheap: it is arithmetic on the Wilson interval, no simulation.
- A worked single-element example, for readers who want to see one case end to end.
- Something showing what a power-aware training label does to a downstream model's performance. This
  would be the strongest possible motivation and is also the largest amount of extra work — it needs
  retraining scE2G or rE2G, which is out of scope unless a collaborator supplies it.
