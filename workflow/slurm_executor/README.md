# SLURM executor

Direct `sbatch` scripts for running the refactored power-analysis pipeline on Sherlock.

These exist to **validate the refactor at full scale before any Nextflow work starts**. They are
deliberately plain: one script per pipeline step, no orchestration layer, so that a failure is
attributable to the step rather than to the runner. Nextflow (Step 7 in
[`docs/status.md`](../../docs/status.md)) replaces them once the validation below passes.

## Run order

```sh
cd /home/groups/engreitz/Users/emattei/git/WattEG

sbatch workflow/slurm_executor/00_setup_env.sbatch            # prerequisite, not a pipeline step
sbatch workflow/slurm_executor/01_prepare_sim_input.sbatch
sbatch workflow/slurm_executor/02_split_pairs.sbatch
sbatch workflow/slurm_executor/03_smoke_simulation.sbatch     # minutes; do not skip
workflow/slurm_executor/04_submit.sh                          # ~814 CPU-hours, ~50 min wall
sbatch workflow/slurm_executor/05_compute_power.sbatch
sbatch workflow/slurm_executor/06_summarize_power.sbatch
sbatch workflow/slurm_executor/07_compare_old_new.sbatch
```

Each step depends on the previous one's outputs and checks for them on entry, so submit them one
at a time and read the log before continuing. To chain them without waiting, use
`sbatch --dependency=afterok:<jobid>`.

Steps can be submitted from anywhere — each resolves the repository root from its own path.

## What each step does

| Step | Script | Resources | Output |
|---|---|---|---|
| 0 | `00_setup_env` | 4 cpu, 16G, 2h | `.pixi/` + pinned sceptre in `.pixi/rlibs` |
| 1 | `01_prepare_sim_input` | 1 cpu, 16G, 2h | `results/refactor/prepared/` |
| 2 | `02_split_pairs` | 1 cpu, 4G, 30m | `results/refactor/splits/` |
| 3 | `03_smoke_simulation` | 1 cpu, 8G, 1h | `results/refactor/smoke/` |
| 4 | `04_run_power_simulation` | array 1-1000 (unthrottled), 1 cpu, 8G, 4h | `results/refactor/sim/es0.15/` |
| 5 | `05_compute_power` | 1 cpu, 32G, 2h | `results/refactor/power/power_es0.15.tsv` |
| 6 | `06_summarize_power` | 1 cpu, 8G, 30m | `results/refactor/power_summary.tsv` |
| 7 | `07_compare_old_new` | 1 cpu, 32G, 2h | `results/refactor/comparison/` |

Logs land in `logs/refactor/`.

All steps request `--partition=engreitz,owners`. Slurm schedules on whichever frees up first;
`engreitz` alone has only 288 cores and was queueing for a long time. Note that `owners` jobs are
**preemptible** and capped at 2 days — every request here is 2 hours or less, and step 4 skips
already-complete tasks on resubmission, so a preempted array task costs only its own work.

## Measured cost

Both columns are measured on this cluster, on this dataset: the old pipeline from `sacct` on the
2026-05-14 run, the refactored one from the step 3 smoke test (9.68 s per target per replicate,
36 targets timed, over 3,026 targets at 100 replicates).

| | Old pipeline | Refactored |
|---|---|---|
| CPU-hours per effect size | 1,308 | **814** |
| Tasks | 1,004 | 1,000 |
| Mean per task | 78 min | ~49 min |
| Wall clock | 14.9 h | **~50 min** |

**The CPU speedup is 1.6x, not the 4x that `docs/status.md` implies.** That document's ~295
CPU-hour figure came from a laptop measurement of 3.4–4.8 s per `run_discovery_analysis()` call;
on these nodes it is ~9.7 s. The refactor's gains are real but concentrated in the per-task read
(480 MB → 18 MB) and in making the work schedulable; sceptre's own per-call cost is untouched and
now dominates.

Most of the wall-clock win is therefore the unthrottled array, not the CPU saving. Past a few
hundred concurrent tasks, queue wait and preemption on `owners` govern turnaround.

**A full 5–50% sweep is expensive.** Ten effect sizes is ~8,100 CPU-hours, and the 2,900-task
budget forces 290 splits, so ~2.8 h wall clock. Worth reading the cascading-effect-size note in
`docs/status.md` before committing to that.

## Configuration

Everything tunable lives in [`config.sh`](config.sh) — sample, paths, effect sizes, replicate
count, split count, seed, partition. It is sourced by every script, so it is the only file to
edit.

The parameters mirror `config/config.yml`. They are duplicated rather than parsed out of the YAML
because parsing YAML in bash is not worth the fragility, but **if you change one, change the
other**.

Two values are not in `config.sh` and cannot be, because `#SBATCH` directives are read before the
script runs:

- `--array=1-1000` in `04_run_power_simulation.sbatch` must match `N_SPLITS`. The script checks
  this at runtime and fails loudly rather than silently skipping splits.
- memory and walltime requests, which are per-script.

## Why this runs in a container

Two separate problems, both in `config.sh`:

**1. The pixi module is too old.** `pixi/0.53.0` — the only version Sherlock ships — reads
lockfile format v6 and earlier; this repository's `pixi.lock` is v7. Re-solving the lock with the
older pixi would discard the pinned solve, so a current pixi is installed standalone instead:

```sh
PIXI_HOME=$GROUP_HOME/software/pixi curl -fsSL https://pixi.sh/install.sh | bash
```

**The installer appends a `PATH` line to `~/.bashrc` — delete it.** `$GROUP_HOME` is NFS, and
`PATH` entries on slow filesystems are a documented cause of hung logins and
`user env retrieval failed requeued held` job failures. `config.sh` sets `PATH` explicitly, which
also makes it work in batch jobs where `~/.bashrc` may not be sourced.

**2. Sherlock's glibc is too old for the locked packages.** Every node is CentOS 7 with glibc
2.17. conda-forge dropped glibc 2.17 after CentOS 7 went end-of-life, so the pinned `r-base 4.4`
requires glibc ≥ 2.28 and pixi refuses to install it:

```
× Cannot install environment 'default'
  Virtual package '__glibc >=2.28' does not match any of the available virtual
  packages on your machine: [__glibc=2.17=0, ...]
```

**This is not something pixi, mise, conda or any other userspace tool can fix.** Conda packages
link against the *host's* C library through the system dynamic loader, so pixi can only detect
glibc, never supply it — `__glibc` is a *virtual* package for exactly that reason. Refusing to
install binaries that would fail at load time is correct behaviour.

So everything runs inside an Ubuntu 22.04 container (glibc 2.35). The image is stock
`ubuntu:22.04` with nothing added — the pixi binary is statically linked and is simply
bind-mounted in:

```sh
APPTAINER_CACHEDIR=$SCRATCH/apptainer_cache \
  apptainer pull $GROUP_HOME/software/containers/ubuntu2204.sif docker://ubuntu:22.04
```

Two details that cost a failed job each, both now handled in `config.sh`:

- `HOME` must be set with Apptainer's `--home` flag. `--env HOME=...` is rejected outright
  (`Overriding HOME environment variable with APPTAINERENV_HOME is not permitted`).
- A stock `ubuntu:22.04` ships **no CA certificates** — `/etc/ssl/certs` does not exist — so pixi
  cannot open an HTTPS connection (`No CA certificates were loaded from the system`). The host
  bundle is bind-mounted in and pointed at with `SSL_CERT_FILE`. This only matters while
  downloading; once `.pixi/` is populated the pipeline steps need no network.

The alternatives were to re-solve the lock against CentOS 7-era packages (discards the pinned
solve, and `r-base 4.4` may have no such build) or to use Sherlock's own `R/4.4.2` module (drops
the lockfile as a description of what actually ran). The container keeps `pixi.lock`
authoritative, and Nextflow supports Apptainer natively, so this carries into Step 7 rather than
being throwaway scaffolding.

## Checkpoints

**After step 0** — `pixi run check-api` must print `All checks passed`. This is the gate that the
pinned sceptre still exposes the unexported S4 slots the pipeline reads. Nothing downstream is
trustworthy without it.

**After step 1** — the job prints the derived discovery p-value threshold alongside the value
implied by the old pipeline's `sceptre_discovery_results.csv`. **These must agree.** The
threshold decides which simulated replicates count as detections, so a mismatch shifts every
power estimate. If `@discovery_result` turns out to be empty on this object,
`prepare_sim_input.R` will error or emit a non-finite threshold; the fallback is to pass
`--alpha` explicitly to `compute_power.R` in step 5.

**After step 3** — the smoke job checks that every pair has exactly the requested number of
replicates and that p-values are not uniformly `NA`. An all-`NA` p-value column is the signature
of the simulated matrix reaching sceptre in the wrong class, which fails *silently* rather than
erroring (see the class notes in `lib/simulate.R`).

**After step 4** — step 5 refuses to run unless all 1000 split outputs are present, and names
the missing ones.

**After step 7** — read `results/refactor/comparison/comparison_report.txt`.

## Resuming a partial run

Step 4 skips any task whose output is already complete, judged by row count
(`reps x pairs + 1` lines). A task killed mid-write leaves a short file, which is detected and
redone rather than silently under-counting replicates downstream. So after a partial failure:

```sh
sbatch workflow/slurm_executor/04_run_power_simulation.sbatch          # redoes only what is missing
sbatch --array=17,42,103 workflow/slurm_executor/04_run_power_simulation.sbatch  # or name them
```

Find the failed tasks with `sacct -j <arrayjobid> --format=JobID,State,Elapsed,MaxRSS`.

## Resource requests are provisional

The memory and walltime figures come from `docs/status.md`, which measured them on `sample1` on
a laptop. They are deliberately generous for the first run on this sample. Every script prints
`sacct` usage on exit — use those numbers to tighten the requests before quoting any of this as a
cost estimate, and to fill in Comparison 3 in `docs/status.md`.

## Scope

- One effect size (`0.15`), because that is the only one the old run covers and therefore the
  only one that can be validated. Adding `0.2` to `EFFECT_SIZES` in `config.sh` and re-running
  steps 4–6 enables the monotonicity check (Comparison 2 in `docs/status.md`).
- `--cpus-per-task=1` on the simulation is deliberate: sceptre is called with `parallel = FALSE`
  and all parallelism comes from the array fan-out.
