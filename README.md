# WattEG

**W**atts for **E**lement–**G**ene pairs: power analysis for element–gene pairs in single-cell
CRISPR screens, built on [pysceptre](https://github.com/broadinstitute/pysceptre) — a Python port
of [sceptre](https://katsevich-lab.github.io/sceptre/)'s statistical engine.

**📖 [Documentation](https://engreitzlab.github.io/WattEG/)**

This repository is the **tool**: the pipeline, how to run it, and what it produces. The analyses and
figures for the methods paper live separately, in
[broadinstitute/WattEG-paper](https://github.com/broadinstitute/WattEG-paper) — split out on
2026-09-04 so the two can move at their own pace.

Given a sceptre object from a completed screen, it answers: for each element–gene pair, if the
element really did reduce expression of the gene by X%, how often would we have detected it? The
output is a per-pair power estimate with a confidence interval, plus the smallest tested effect size
at which each pair becomes detectable.

The method comes from [DC_TAP_Paper](https://github.com/EngreitzLab/DC_TAP_Paper), where it was
written for one specific analysis. This repository generalises it.

## Installation

```sh
pixi install
```

That is the whole of it. The environment is one Python package plus Nextflow: no R, no pinned
sceptre commit, no patch applied to it, and no separate install step. pysceptre enters as a git
dependency pinned to a commit, because a sweep has to be re-runnable against the engine that
produced it.

**The R implementation has not been deleted.** It is the reference the Python path was validated
against and what the methods paper describes, and it still lives in `src/*.R` and `lib/*.R`. What
is gone is the environment that ran it; to run it, use the **`r-implementation`** branch.

That branch is **maintained, not frozen**. Two bugs found during the port change results, and both
were fixed in R as well as in Python rather than being quarantined on a snapshot -- nothing was
published, so there were no numbers owed a bit-for-bit reproduction, and freezing the branch would
only have preserved the bugs. Sweeps produced before 2026-09-21 predate both fixes.
(`legacy` is something else again -- the Snakemake implementation that preceded both.)

## Input

**One file: a `.h5mu` dataset** exported from a sceptre object on which `assign_grnas()` and
`run_qc()` have been called, using `grna_integration_strategy = "union"`.

Everything else is derived from it — the discovery pairs, the gRNA-to-target mapping and the
significance threshold all travel with the export.

Making one is a one-off step per dataset, run wherever R and sceptre are available. It is the only
place R appears at all, and it is not part of the pipeline:

```sh
Rscript pysceptre/scripts/export_sceptre_dataset.R \
  --sceptre-object results/sample1/sceptre_object.rds \
  --out-dir export/ --all-genes --all-cells
python  pysceptre/scripts/make_h5mu.py export/
```

**`--all-cells` is required, not optional.** DESeq2 "poscounts" size factors are a per-cell
reduction against a per-gene geometric mean taken over every cell in the object, so an export
restricted to the QC-passing cells gives different size factors for the cells that remain. The
export handles odm-backed (out-of-core) matrices itself, which is why there is no longer an
`--response-odm` flag anywhere in the pipeline.

List your samples in a CSV (see `assets/samplesheet.csv`):

```csv
sample,dataset
sample1,export/dataset.h5mu
```

## Quickstart

```sh
nextflow run . -profile sherlock -params-file config/config.yml
```

Each step is also a console script with `--help`, so a sweep can be driven by hand or by another
runner:

```sh
# derive the simulation inputs (once per sample)
watteg-prepare-sim-input --dataset export/dataset.h5mu --outdir prepared/

# split targets into per-task chunks
watteg-split-pairs --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/

# simulate (once per split x effect size)
watteg-run-power-simulation \
  --prepared prepared/ --pairs splits/split_001.tsv \
  --effect-size 0.15 --reps 100 --seed 20250812 \
  --out sim/split_001_es0.15.tsv

# power per pair, then one table across effect sizes
watteg-compute-power \
  --simulations sim/ --threshold-file prepared/discovery_threshold.txt \
  --out power_es0.15.tsv

watteg-summarize-power \
  --power power_es0.15.tsv power_es0.2.tsv \
  --sim-input prepared/sim_input.h5 --out power_summary.tsv
```

Parameters live in `config/config.yml`.

## Output

`power_summary.tsv` — one row per pair:

| Column | |
|---|---|
| `power_at_effect_size_15` | power at a 15% knockdown, one column per effect size |
| `power_at_effect_size_15_ci_low` / `_ci_high` | 95% Wilson interval |
| `power_at_effect_size_15_n_reps` | simulations behind the estimate |
| `min_detectable_effect_size` | smallest tested effect size reaching the target power |

Always read a power estimate together with its interval: `power = 0` at 100 simulations has a 95%
upper bound of 0.037, so it means "not detected in 100 tries", not "undetectable".

See [Output](https://engreitzlab.github.io/WattEG/output/) for every column.

## Choosing parameters

`num_replicates` is the only parameter that changes results; `n_splits` and `reps_per_chunk` only
change how the work is divided. How many simulations you need depends on what you do with the
numbers — read
[Choosing num_replicates](https://engreitzlab.github.io/WattEG/choosing-num-replicates/).

Leave `n_control_cells` unset. Sampling control cells looks like a large speedup but costs 21–60% of
your power, because sceptre's conditional randomisation test needs enough cells to resolve the null
tail at the significance threshold. Measured numbers are in
[Methods](https://engreitzlab.github.io/WattEG/methods/).

## Status

The five steps above are complete, run standalone, and have been **validated against the previous
implementation**: identical pair sets, zero difference in perturbed cells per pair across all 34,886
pairs, and per-pair power correlating at r = 0.993 with no directional bias. A Nextflow workflow with
a SLURM profile to wire them together is in progress; `config/config.yml` already holds the
parameters it will consume.

The previous Snakemake implementation — `Snakefile`, `rules/`, `R/` and `envs/` — has been removed now
that the comparison is done. It is preserved on the **`legacy`** branch if you need to consult or
re-run it.

Layout: `src/` holds the pipeline executables, `lib/` the shared code they source, and `workflow/`
the cluster scripts and comparison tools.

**[Status and handoff](https://engreitzlab.github.io/WattEG/status/)** has the
full picture: what is done and verified, what is left, reference numbers for sizing a cluster run,
ready-to-use SLURM array scripts, and how to run the old-vs-new comparison.

## Documentation

- [Usage](https://engreitzlab.github.io/WattEG/usage/) — every parameter, and running each step by hand
- [Output](https://engreitzlab.github.io/WattEG/output/) — every output column
- [Choosing num_replicates](https://engreitzlab.github.io/WattEG/choosing-num-replicates/) — precision, cost, confidence intervals
- [Methods](https://engreitzlab.github.io/WattEG/methods/) — how the simulation is parameterised
- [Development](https://engreitzlab.github.io/WattEG/development/) — environment, sceptre pinning, conventions

## License

MIT — see [LICENSE](LICENSE).
