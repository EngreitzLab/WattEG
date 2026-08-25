# element-gene-power-analysis

Power analysis for element–gene pairs in single-cell CRISPR screens, built on
[sceptre](https://katsevich-lab.github.io/sceptre/).

**📖 [Documentation](https://engreitzlab.github.io/element-gene-power-analysis/)**

Given a sceptre object from a completed screen, it answers: for each element–gene pair, if the
element really did reduce expression of the gene by X%, how often would we have detected it? The
output is a per-pair power estimate with a confidence interval, plus the smallest tested effect size
at which each pair becomes detectable.

The method comes from [DC_TAP_Paper](https://github.com/EngreitzLab/DC_TAP_Paper), where it was
written for one specific analysis. This repository generalises it.

## Installation

```sh
pixi install
pixi run setup       # installs sceptre from a pinned commit
pixi run check-api   # verifies that pin against the internals the pipeline uses
```

sceptre is not on conda-forge or bioconda, so it is installed from a pinned commit rather than
captured in `pixi.lock`. The pipeline reads several of its unexported S4 slots, which is why the
version is pinned and checked.

## Input

**One file: a sceptre object** (`.rds`) on which `assign_grnas()` and `run_qc()` have been called,
using `grna_integration_strategy = "union"`.

Everything else is derived from it — the discovery pairs, the gRNA-to-target mapping and the
significance threshold all already live inside the object.

If the object's response matrix is odm-backed (out-of-core), pass `--response-odm` on
`prepare_sim_input.R` with the path to the backing `.odm` file — see [Usage](https://engreitzlab.github.io/element-gene-power-analysis/usage/).

List your samples in a CSV (see `assets/samplesheet.csv`):

```csv
sample,sceptre_object
sample1,results/sample1/sceptre_object.rds
```

Add an optional `response_odm` column (see `assets/samplesheet_seqera_test.csv`) for any sample whose
sceptre object's response matrix is odm-backed — same file `prepare_sim_input.R --response-odm` takes
standalone. Leave it blank, or omit the column entirely, for in-memory-backed objects.

## Quickstart

Each step is a standalone executable in `src/` with `--help`.

```sh
# derive the simulation inputs (once per sample)
Rscript src/prepare_sim_input.R \
  --sceptre-object results/sample1/sceptre_object.rds \
  --outdir prepared/

# split targets into per-task chunks
Rscript src/split_pairs.R --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/

# simulate (once per split x effect size)
Rscript src/run_power_simulation.R \
  --sim-input prepared/sim_input.rds \
  --sceptre-template prepared/sceptre_template.rds \
  --pairs splits/split_001.tsv \
  --grna-targets prepared/grna_targets.tsv \
  --effect-size 0.15 --reps 100 --seed 20250812 \
  --out sim/split_001_es0.15.tsv

# power per pair, then one table across effect sizes
Rscript src/compute_power.R \
  --simulations "$(ls sim/*_es0.15.tsv | paste -sd, -)" \
  --threshold-file prepared/discovery_threshold.txt \
  --out power_es0.15.tsv

Rscript src/summarize_power.R \
  --power power_es0.15.tsv,power_es0.2.tsv \
  --out power_summary.tsv
```

Parameters live in `config/config.yml`.

## Output

`power_summary.tsv` — one row per pair:

| Column | |
|---|---|
| `power_at_effect_size_15` | power at a 15% knockdown, one column per effect size |
| `power_at_effect_size_15_ci_low` / `_ci_high` | 95% Wilson interval |
| `power_at_effect_size_15_n_reps` | replicates behind the estimate |
| `min_detectable_effect_size` | smallest tested effect size reaching the target power |

Always read a power estimate together with its interval: `power = 0` at 100 replicates has a 95%
upper bound of 0.037, so it means "not detected in 100 tries", not "undetectable".

See [Output](https://engreitzlab.github.io/element-gene-power-analysis/output/) for every column.

## Choosing parameters

`num_replicates` is the only parameter that changes results; `n_splits` and `reps_per_chunk` only
change how the work is divided. How many replicates you need depends on what you do with the
numbers — read
[Choosing num_replicates](https://engreitzlab.github.io/element-gene-power-analysis/choosing-num-replicates/).

Leave `n_control_cells` unset. Sampling control cells looks like a large speedup but costs 21–60% of
your power, because sceptre's conditional randomisation test needs enough cells to resolve the null
tail at the significance threshold. Measured numbers are in
[Methods](https://engreitzlab.github.io/element-gene-power-analysis/methods/).

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

**[Status and handoff](https://engreitzlab.github.io/element-gene-power-analysis/status/)** has the
full picture: what is done and verified, what is left, reference numbers for sizing a cluster run,
ready-to-use SLURM array scripts, and how to run the old-vs-new comparison.

## Documentation

- [Usage](https://engreitzlab.github.io/element-gene-power-analysis/usage/) — every parameter, and running each step by hand
- [Output](https://engreitzlab.github.io/element-gene-power-analysis/output/) — every output column
- [Choosing num_replicates](https://engreitzlab.github.io/element-gene-power-analysis/choosing-num-replicates/) — precision, cost, confidence intervals
- [Methods](https://engreitzlab.github.io/element-gene-power-analysis/methods/) — how the simulation is parameterised
- [Development](https://engreitzlab.github.io/element-gene-power-analysis/development/) — environment, sceptre pinning, conventions

## License

MIT — see [LICENSE](LICENSE).
