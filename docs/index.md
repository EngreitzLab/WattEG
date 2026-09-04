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

For each perturbation target, and for each replicate:

1. Draw per-gRNA effect sizes around the requested effect size, so guides targeting the same
   element differ in strength.
2. Simulate a count matrix from each gene's mean and dispersion, applying those effect sizes to
   the perturbed cells.
3. Run sceptre's discovery analysis on the simulated data.
4. Record whether the pair would have been called significant.

Power is the fraction of replicates in which it would have been. Because that is a binomial
proportion over a finite number of replicates, every estimate is reported with a
[Wilson confidence interval]({{ site.baseurl }}{% link choosing-num-replicates.md %}).

## Quickstart

```sh
# 1. environment (sceptre is installed from a pinned commit, not from a conda channel)
pixi install
pixi run setup
pixi run check-api        # asserts the pinned sceptre exposes the internals this pipeline uses

# 2. derive the simulation inputs from your sceptre object
Rscript src/prepare_sim_input.R \
    --sceptre-object results/sample1/sceptre_object.rds \
    --outdir prepared/

# 3. split the pairs into per-task chunks
Rscript src/split_pairs.R --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/

# 4. simulate one chunk (repeat per split and per effect size)
Rscript src/run_power_simulation.R \
    --sim-input prepared/sim_input.rds \
    --sceptre-template prepared/sceptre_template.rds \
    --pairs splits/split_001.tsv \
    --grna-targets prepared/grna_targets.tsv \
    --effect-size 0.15 --reps 100 --seed 20250812 \
    --out sim/split_001_es0.15.tsv

# 5. power per pair, then one table across effect sizes
Rscript src/compute_power.R --simulations "$(ls sim/*_es0.15.tsv | paste -sd, -)" \
    --threshold-file prepared/discovery_threshold.txt --out power_es0.15.tsv
Rscript src/summarize_power.R --power power_es0.15.tsv,power_es0.2.tsv --out power_summary.tsv
```

Every script is standalone and self-documenting via `--help`. See
[Usage]({{ site.baseurl }}{% link usage.md %}) for all parameters.

## Input

**One file: a sceptre object** on which `assign_grnas()` and `run_qc()` have been called, using
`grna_integration_strategy = "union"`.

Everything else is derived from it. Earlier versions of this pipeline additionally required
`gene_grna_group_pairs.rds`, `grna_groups_table.rds` and a discovery-results file; all three
duplicated data already present in the object:

| Previously a separate input | Read instead from |
|---|---|
| `gene_grna_group_pairs.rds` | `@discovery_pairs_with_info` (which also carries `pass_qc`) |
| `grna_groups_table.rds` | `@grna_target_data_frame` |
| discovery results | `@discovery_result` |

## Choosing parameters

The one parameter that changes your *results* is `num_replicates`; `n_splits` and
`reps_per_chunk` only change how the work is divided. Read
[Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}) before
picking one — the right value depends on whether you report aggregate power, per-pair power, or
make per-pair decisions at a cutoff.

Two settings deserve a warning, both documented in
[Methods]({{ site.baseurl }}{% link methods.md %}):

- **`n_control_cells`** looks like a large speedup but biases power downward substantially
  (measured: 29% relative loss at 5,000 controls). Leave it unset.
- **`alpha`** should normally be left unset so the threshold is derived from the real discovery
  results, which reflects the multiple-testing correction actually applied.

## Documentation

- [Usage]({{ site.baseurl }}{% link usage.md %}) — every parameter, and running each step by hand
- [Output]({{ site.baseurl }}{% link output.md %}) — every column of every output file
- [Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}) — precision, cost, and the confidence intervals
- [Methods]({{ site.baseurl }}{% link methods.md %}) — how the simulation is parameterised
- [Development]({{ site.baseurl }}{% link development.md %}) — environment, sceptre pinning, conventions
