---
title: Overview
nav_order: 1
---

# WattEG

Power analysis for element–gene pairs in single-cell CRISPR screens, built on
[sceptre](https://katsevich-lab.github.io/sceptre/).

Given a sceptre object from a completed screen, it answers: **for each element–gene pair, if the
element really did reduce expression of the gene by X%, how often would we have detected it?**
The output is a per-pair power estimate with a confidence interval, plus the smallest tested effect
size at which each pair becomes detectable.

The method is a generalisation of the analysis in
[DC_TAP_Paper](https://github.com/EngreitzLab/DC_TAP_Paper), which was written for one specific
dataset.

## How it works

For each perturbation target, and for each simulation:

1. Draw per-gRNA effect sizes around the requested effect size, so guides targeting the same
   element differ in strength.
2. Simulate a count matrix from each gene's mean and dispersion, applying those effect sizes to
   the perturbed cells.
3. Run sceptre's discovery analysis on the simulated data.
4. Record whether the pair would have been called significant.

Power is the fraction of simulations in which it would have been. Because that is a binomial
proportion over a finite number of simulations, every estimate is reported with a
[Wilson confidence interval]({{ site.baseurl }}{% link choosing-num-replicates.md %}).

## Quickstart

```sh
# 1. environment -- one Python package plus Nextflow, no R
pixi install

# 2. the whole pipeline
nextflow run . -profile sherlock -params-file config/config.yml
```

Or step by step, each with `--help`:

```sh
watteg-prepare-sim-input --dataset export/dataset.h5mu --outdir prepared/
watteg-split-pairs --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/
watteg-run-power-simulation \
    --prepared prepared/ --pairs splits/split_001.tsv \
    --effect-size 0.15 --reps 100 --seed 20250812 \
    --out sim/split_001_es0.15.tsv
watteg-compute-power --simulations sim/ \
    --threshold-file prepared/discovery_threshold.txt --out power_es0.15.tsv
watteg-summarize-power --power power_es0.15.tsv power_es0.2.tsv \
    --sim-input prepared/sim_input.h5 --out power_summary.tsv
```

See [Usage]({{ site.baseurl }}{% link usage.md %}) for all parameters, including how to produce the
`.h5mu` from a sceptre object.

## Input

**One file: a `.h5mu`**, exported once from a sceptre object on which `assign_grnas()` and
`run_qc()` have been called, using `grna_integration_strategy = "union"`. Everything else — the
discovery pairs, the gRNA-to-target mapping, the significance threshold, the analysis parameters —
travels with it.

Producing that export is the only step that needs R, and it is not part of the pipeline; see
[Usage]({{ site.baseurl }}{% link usage.md %}). It must be written with `--all-genes --all-cells`,
for reasons that are not cosmetic: the size factors reduce over the whole gene set, and are
computed against a geometric mean over every cell including the ones QC removed.

## Choosing parameters

The one parameter that changes your *results* is `num_replicates`; `n_splits` and
`reps_per_chunk` only change how the work is divided, and that is checked end to end rather than
assumed — see Reproducibility in [Usage]({{ site.baseurl }}{% link usage.md %}). Read
[Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}) before
picking one — the right value depends on whether you report aggregate power, per-pair power, or
make per-pair decisions at a cutoff.

Two settings deserve a warning, both documented in
[Methods]({{ site.baseurl }}{% link methods.md %}):

- **`alpha`** should normally be left unset so the threshold is derived from the real discovery
  results, which reflects the multiple-testing correction actually applied.
- **`expression_model`** should be left at `fitted`. `size_factor` exists only to reproduce sweeps
  published before 2026-09-21; it mixes two statistical models of the same data and reproduces
  86.5% of the observed count variance against the fitted model's 99.5%.

`n_control_cells` and `cell_batches` are gone. They sampled control cells to buy speed, cost 21-60%
of power, and were never on; the Python path is fast enough that the trade has no upside. See
[Plan - pysceptre backend]({{ site.baseurl }}{% link pysceptre-backend.md %}) section 9 for why
dropping the batch stratification does not expose the arms to drift.

## Documentation

- [Usage]({{ site.baseurl }}{% link usage.md %}) — every parameter, and running each step by hand
- [Output]({{ site.baseurl }}{% link output.md %}) — every column of every output file
- [Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}) — precision, cost, and the confidence intervals
- [Methods]({{ site.baseurl }}{% link methods.md %}) — how the simulation is parameterised
- [Development]({{ site.baseurl }}{% link development.md %}) — environment, sceptre pinning, conventions
