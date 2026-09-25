---
title: Usage
nav_order: 2
---

# Usage

```sh
nextflow run . -profile sherlock -params-file config/config.yml
```

That is the whole pipeline. Every step is also a console script taking explicit arguments and
printing `--help`, so any one of them can be run, re-run or debugged on its own; nothing reads a
global config object.

---

## Getting a dataset in

The pipeline's input is a `.h5mu`, not a sceptre object. Producing one is a one-off step per
dataset, and it is **the only place R appears**:

```sh
Rscript pysceptre/scripts/export_sceptre_dataset.R \
    --sceptre-object results/sample1/sceptre_object.rds \
    --out-dir export/ --all-genes --all-cells
python  pysceptre/scripts/make_h5mu.py export/
```

The object must have been through `assign_grnas()` and `run_qc()` with
`grna_integration_strategy = "union"`.

**Both flags are required, and for different reasons.**

`--all-genes` because the poscounts size factors are a per-cell reduction over the *whole gene set*.

`--all-cells` because they are also computed against a per-gene geometric mean taken over *every
cell in the object*, including the ones QC removed. An export restricted to the QC-passing cells
gives different size factors for the cells that remain, and different normalised means for every
gene — measured on day0, a 0.46 % shift in size factor and 2.2 % in the raw gene mean. Nothing
downstream would report it. `watteg-prepare-sim-input` warns if it is handed an export without
them.

The export handles odm-backed (out-of-core) matrices itself, which is why nothing in the pipeline
takes an `--response-odm` flag any more.

---

## The five steps

```
dataset.h5mu
      |
      |  watteg-prepare-sim-input      once per sample
      +--> sim_input.h5                per-gene and per-cell statistics, the model
      |                                coefficients, and the perturbation matrices
      +--> pairs.tsv                   QC-passing discovery pairs
      +--> pairs_with_info.tsv         every pair with its real-data QC counts
      +--> grna_targets.tsv            gRNA -> target mapping (many-to-many)
      +--> discovery_threshold.txt     the p-value a replicate must beat
      +--> analysis_mode.tsv           which test the screen ran, and its resampling budget
      |
      |  watteg-split-pairs            once per sample
      +--> split_001.tsv ...           balanced chunks of targets
      |
      |  watteg-fit-null-models        once per sample, under --null-fits reuse
      +--> null_fits.h5                each gene's null-model fit, per simulation
      |
      |  watteg-run-power-simulation   once per (split, effect size, replicate chunk)
      +--> sim_*.tsv.gz                one row per (pair, replicate)
      |
      |  watteg-consolidate-replicates once per effect size
      +--> es*.parquet                 the same rows, as one file
      |
      |  watteg-compute-power          once per effect size
      +--> power_es*.tsv               one row per pair, with confidence intervals
      |
      |  watteg-summarize-power        once per sample
      +--> power_summary.tsv           one row per pair, one column per effect size
```

There is no `sceptre_template.rds`. It was an R object carrying the covariate matrix and the
analysis parameters, which now live in `sim_input.h5` and `analysis_mode.tsv`.

`watteg-fit-null-models` is the Python version of R's `FIT_NULL_MODELS`: each gene's null model is
fitted once per simulation, on an independent draw with no knockdown, and every target the gene is
tested with reuses that fit (`--null-fits reuse`). Without it (`--null-fits refit`) each simulation
refits every gene once per target, the exact configuration, at about 1.7× the cost of a simulated
pair-test. See [Methods](methods.md#null-model-fits).

---

## 1. `watteg-prepare-sim-input`

```sh
watteg-prepare-sim-input --dataset export/dataset.h5mu --outdir prepared/
```

| Option | Default | Meaning |
|---|---|---|
| `--dataset` | required | An `--all-genes --all-cells` export. |
| `--outdir` | required | Directory for all six outputs. |
| `--threshold` | derived | Set the significance threshold instead of deriving it from the export's discovery result. |
| `--n-jobs` | 1 | Workers for the per-gene fits. |

Fits each gene's Poisson GLM and negative-binomial theta — the same model sceptre's own
`perform_response_precomputation` fits, reproduced to about ten significant digits — and stores the
coefficients. A gene whose theta hits the estimator's bounds is **refused** rather than simulated
from: a clamped theta is not an estimate, and drawing counts from it would state a noise level the
data never supported. (sceptre clamps and carries on, which is reasonable for an analysis and not
for a simulation.)

---

## 2. `watteg-split-pairs`

```sh
watteg-split-pairs --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/
```

| Option | Default | Meaning |
|---|---|---|
| `--n-splits` | required | Number of chunks. Must not exceed the number of targets. |
| `--target-overhead` | 2.0 | Pairs-equivalent fixed cost of a target, from the measured cost model. |
| `--prefix` | `split_` | Filename prefix; files are zero-padded so lexicographic order is numeric order. |

A target's pairs cannot be separated — the simulation draws one count matrix per (target,
replicate) and tests every one of that target's genes against it — so targets are the unit and the
job is bin packing. Weighted by pairs **plus an overhead**, because a task's cost is
`intercept + slope × pairs` and weighting by pairs alone over-fills the splits that hold many small
targets.

---

## 3. `watteg-run-power-simulation`

```sh
watteg-run-power-simulation \
    --prepared prepared/ --pairs splits/split_001.tsv \
    --effect-size 0.15 --reps 100 --seed 20250812 \
    --out sim/split_001_es0.15.tsv
```

| Option | Default | Meaning |
|---|---|---|
| `--effect-size` | required | A **fractional decrease**: 0.15 is a 15 % knockdown. |
| `--reps` | required | Replicates in this chunk. |
| `--rep-offset` | 0 | Replicates already covered by earlier chunks, so `rep` stays unique. |
| `--seed` | required | Results are stochastic and must be reproducible. |
| `--guide-spread-c` | 0.65 | Guide-to-guide spread: each guide's knockdown is Beta with mean es and sd c·es·(1−es), so zero at es = 0. Replaces `--guide-sd`, which is refused. See [Methods](methods.md). |
| `--n-jobs` | 8 | Workers for the per-gene tests. |
| `--expression-model` | `fitted` | Where a gene's unperturbed expected counts come from. |

**`--effect-size 0` is allowed on purpose.** It is the null arm: perturbed cells are simulated from
the unperturbed mean and the screen's own test is run on them, so the p-values should be uniform.
It is the only way to measure this pipeline's type-I error from its own output. Run it as its own
sweep, not as another point on a curve.

**`--expression-model`.** `fitted` draws from `exp(X·β)`, the expected count sceptre's own null
model gives that cell, so the simulation and the test that judges it are on one scale.
`size_factor` reproduces the pre-2026-09-21 behaviour — a size-factor-normalised gene mean scaled
by the cell's poscounts factor — and exists only to compare against sweeps already run with it. It
mixes two models, runs about 4 % low, and reproduces 86.5 % of the observed count variance against
the fitted model's 99.5 %.

---

## 3b. `watteg-fit-null-models`

```sh
watteg-fit-null-models --prepared prepared/ --reps 100 --seed 20250812 --n-jobs 8 \
    --out prepared/null_fits.h5
```

| Option | Default | Meaning |
|---|---|---|
| `--reps`, `--rep-offset` | required, 0 | Fits simulations offset+1 .. offset+reps. |
| `--seed` | required | The simulation run's seed; recorded, and checked by the simulation. |
| `--expression-model` | `fitted` | The simulation run's baseline; recorded and checked. |
| `--pairs` | `pairs.tsv` in `--prepared` | Whose genes to fit. |
| `--n-jobs` | 8 | Workers; one simulation's fits are one unit of work. |

Pass the file to `watteg-run-power-simulation --null-fits reuse --null-fits-file`. A task given no
file fits its own genes first with the same keyed draws, so the output is the same bytes either way;
the file only saves refitting a gene in every task that tests it. Refused for a CRT screen, whose
simulations the fast driver (the file's only reader) does not run.

---

## 4. `watteg-compute-power` and `watteg-summarize-power`

```sh
watteg-compute-power \
    --simulations sim/ --threshold-file prepared/discovery_threshold.txt \
    --out power_es0.15.tsv

watteg-summarize-power \
    --power power_es0.05.tsv power_es0.15.tsv \
    --sim-input prepared/sim_input.h5 --out power_summary.tsv
```

`--simulations` takes files or directories; a directory expands to the `.tsv`, `.tsv.gz` and
`.parquet` files inside it, sorted. `--power-threshold` (default 0.8) sets what "detectable" means
for the minimum-detectable-effect-size columns.

---

## Reproducibility

A pair's power does not depend on how the work was split up. Every draw is keyed on
`(seed, target, replicate, effect size)` and nothing else, so simulating replicates 1–100 in one
task and in five tasks of twenty gives **byte-identical** results, and so does moving a target from
one split to another. That is checked end to end:

```sh
WATTEG_PREPARED=prepared/ pytest -m realdata
```

The streams differ from the R implementation's and always would — different generators — so the two
agree statistically rather than draw for draw. See
[Plan - pysceptre backend]({{ site.baseurl }}{% link pysceptre-backend.md %}) for what has been
measured.

---

## Sizing a cluster run

The unit of work is (split × effect size × replicate chunk), and essentially all of the compute is
in step 3. A task's cost follows `intercept + slope × pairs`, so:

- halving the **pairs** saves less than half, because the intercept does not move;
- halving the **replicates** saves exactly half;
- effect sizes cost the same as each other.

`--n-jobs` buys wall clock rather than CPU time: measured on a 10-performance-core machine, 8
workers gave 2.5× the wall-clock speed for 1.5× the CPU, and 14 was slower than 8.
