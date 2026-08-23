---
title: Usage
nav_order: 2
---

# Usage

Every step is a standalone executable in `src/` that takes explicit command-line arguments and
prints `--help`. Nothing reads a global config object, so any step can be run, re-run or debugged
on its own.

> **Orchestration is in progress.** The steps below are the complete, working pipeline as run by
> hand or from a shell loop. A Nextflow workflow that wires them together with a SLURM profile is
> the next piece of work and is not part of this repository yet. `config/config.yml` already holds
> the parameters it will consume.

---

## The five steps

```
sceptre_object.rds
      |
      |  prepare_sim_input.R          once per sample
      +--> sim_input.rds              per-gene / per-cell statistics + perturbation matrices
      +--> sceptre_template.rds       the sceptre object with both count matrices emptied
      +--> pairs.tsv                  QC-passing discovery pairs
      +--> grna_targets.tsv           gRNA -> target mapping
      +--> discovery_threshold.txt    the p-value a replicate must beat
      |
      |  split_pairs.R                once per sample
      +--> split_001.tsv ...          balanced chunks of targets
      |
      |  run_power_simulation.R       once per (split, effect size, replicate chunk)
      +--> sim_*.tsv                  one row per (pair, replicate)
      |
      |  compute_power.R              once per effect size
      +--> power_es*.tsv              one row per pair, with confidence intervals
      |
      |  summarize_power.R            once per sample
      +--> power_summary.tsv          one row per pair, one column per effect size
```

---

## 1. `prepare_sim_input.R`

Derives everything the simulation needs from the sceptre object, once, so the parallel tasks read
small files instead of a multi-gigabyte object.

```sh
Rscript src/prepare_sim_input.R \
    --sceptre-object results/sample1/sceptre_object.rds \
    --outdir prepared/
```

| Option | Default | Meaning |
|---|---|---|
| `--sceptre-object` | required | Input sceptre object. The only required input. |
| `--response-odm` | — | Path to the response matrix's backing `.odm` file. Only needed if `--sceptre-object`'s response matrix is odm-backed (out-of-core); ignored otherwise. |
| `--outdir` | — | Directory for all outputs; or name each `--out-*` explicitly. |
| `--out-sim-input` | `<outdir>/sim_input.rds` | Per-gene and per-cell statistics, plus the gRNA and target perturbation matrices. |
| `--out-sceptre-template` | `<outdir>/sceptre_template.rds` | The sceptre object with `@response_matrix` and `@grna_matrix` emptied. |
| `--out-pairs` | `<outdir>/pairs.tsv` | QC-passing discovery pairs. |
| `--out-grna-targets` | `<outdir>/grna_targets.tsv` | gRNA → target mapping. |
| `--out-threshold` | `<outdir>/discovery_threshold.txt` | Largest nominal p-value that survived correction in the real analysis. |
| `--all-genes` | off | Keep every gene, not just those in QC-passing pairs. Inspection only. |
| `--no-compress` | off | Write uncompressed `.rds`. Compression is on by default because it measured both smaller *and* faster to read. |

It logs which categorical covariates are available for `--cell-batches`, which is worth reading:
the column is dataset-specific (e.g. `batch_factor`, `replicate_factor`) and is often **not** called
`batch`.

## 2. `split_pairs.R`

Splits targets into balanced chunks. Whole targets stay together, because the simulation builds
per-target state once and then loops over replicates.

```sh
Rscript src/split_pairs.R --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/
```

| Option | Default | Meaning |
|---|---|---|
| `--pairs` | required | `pairs.tsv` from step 1. |
| `--n-splits` | required | Number of chunks. Must not exceed the number of targets. |
| `--outdir` | `.` | Where to write the chunks. |
| `--prefix` | `split_` | Filename prefix; files are zero-padded so lexical order matches numeric order. |
| `--target-overhead` | `0` | Extra weight per target, in pair-equivalents, to account for the fixed per-target cost of a sceptre call. `0` balances purely on pair count. |

Balancing is on **pairs**, not on target count, and it reports the achieved imbalance so you can
see it. It also refuses to produce empty splits.

## 3. `run_power_simulation.R`

The actual simulation. One invocation handles one split × one effect size × one chunk of replicates.

```sh
Rscript src/run_power_simulation.R \
    --sim-input prepared/sim_input.rds \
    --sceptre-template prepared/sceptre_template.rds \
    --pairs splits/split_001.tsv \
    --grna-targets prepared/grna_targets.tsv \
    --effect-size 0.15 \
    --reps 100 \
    --seed 20250812 \
    --out sim/split_001_es0.15.tsv
```

| Option | Default | Meaning |
|---|---|---|
| `--effect-size` | required | Fractional decrease in expression: `0.15` = 15% knockdown. |
| `--reps` | required | Replicates in this chunk. |
| `--rep-offset` | `0` | Replicates already covered by earlier chunks; keeps the reported `rep` unique. |
| `--seed` | required | Base seed. Required, not optional — see *Reproducibility* below. |
| `--guide-sd` | `0.13` | Spread of per-gRNA effect sizes around the target effect size. |
| `--n-control-cells` | unset | Sample this many controls instead of using all. **Biases power downward; leave unset.** |
| `--cell-batches` | unset | Covariate column to stratify control sampling by. Only with `--n-control-cells`. |
| `--gc-every` | `0` | Call `gc()` every N replicates. `0` disables it. |

### Running the full grid by hand

```sh
for es in 0.15 0.2; do
    for split in splits/split_*.tsv; do
        name=$(basename "$split" .tsv)
        Rscript src/run_power_simulation.R \
            --sim-input prepared/sim_input.rds \
            --sceptre-template prepared/sceptre_template.rds \
            --pairs "$split" --grna-targets prepared/grna_targets.tsv \
            --effect-size "$es" --reps 100 --seed 20250812 \
            --out "sim/${name}_es${es}.tsv"
    done
done
```

## 4. `compute_power.R`

Turns per-replicate results into power per pair, with a Wilson interval.

```sh
Rscript src/compute_power.R \
    --simulations "$(ls sim/*_es0.15.tsv | paste -sd, -)" \
    --threshold-file prepared/discovery_threshold.txt \
    --out power_es0.15.tsv
```

| Option | Default | Meaning |
|---|---|---|
| `--simulations` | required | Comma-separated per-replicate TSVs. |
| `--threshold-file` | — | File from step 1. Mutually exclusive with `--alpha`. |
| `--alpha` | — | Explicit p-value threshold. Only for objects with no discovery results. |
| `--conf-level` | `0.95` | Confidence level for the Wilson interval. |

Run it once per effect size — it refuses input that mixes effect sizes, and refuses duplicated
`(target, gene, replicate)` rows, which is what an overlapping `--rep-offset` would produce.

## 5. `summarize_power.R`

One row per pair, one `power_at_effect_size_<N>` column per effect size, plus the smallest tested effect
size at which each pair reaches a target power.

```sh
Rscript src/summarize_power.R \
    --power power_es0.15.tsv,power_es0.2.tsv \
    --power-threshold 0.8 \
    --out power_summary.tsv
```

---

## Reproducibility

`--seed` is required. Seeds are derived per `(seed, target, replicate, effect_size)`, not per task,
which has two consequences worth relying on:

- **Results do not depend on how work is divided.** The same `--seed` gives identical numbers
  whether you use 1 split or 280, and whether replicates run in one chunk or ten. This is verified
  in the test suite.
- **Runs are extensible.** Going from 100 to 400 replicates leaves replicates 1–100 byte-identical,
  so you can add replicates with `--rep-offset 100` instead of recomputing.

## Sizing a cluster run

Measured on a 292-gene × 586,309-cell dataset with 2,798 targets and 32,386 QC-passing pairs, using
all control cells:

| Quantity | Value |
|---|---|
| One `run_discovery_analysis()` call | 3.4–4.8s |
| Total at 100 replicates | ~295 CPU-hours per effect size |
| Peak memory, simulation task | 1.5 GB → request 4 GB |
| Peak memory, `prepare_sim_input.R` | 7.7 GB → request 12 GB |
| Input read per task | 59 MB, ~1.1s |

**Prefer more splits over more replicate chunks.** An extra split costs one extra ~1.1s input load;
an extra replicate chunk re-pays the *per-target* setup (control selection, guide status, cell
permutation) once per chunk per target. Only chunk replicates when you need more parallelism than
one-target-per-task, or when a single task would exceed your queue's time limit.

For the dataset above, ~10 targets per task is about an hour of work, so `--n-splits 280` with no
replicate chunking gives 280 tasks per effect size — comfortably under the `MaxArraySize` of 1000
that most SLURM installations use.
