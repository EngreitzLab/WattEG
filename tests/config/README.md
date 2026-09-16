# Config test — the resource closures

`pixi run test` covers the pure R functions. This covers something no R test can reach: the
`memory` closures in `conf/base.config`.

## Why it exists

Three processes size their memory from the input they are handed, because their peak is a function
of the input rather than a property of the pipeline. Those closures are Groovy, they are evaluated
by Nextflow at task submission, and they are the least observable code in the repo:

- They run **before a task record exists**, so an exception inside one produces a run with every
  task so far marked `COMPLETED`, no failed task, and no error attributable to anything. Run
  `shrivelled_bhabha` (`zqCjFyEUvTFEE`, 2026-09-16) aborted that way — see `docs/status.md`,
  trap 5.
- They are sensitive to a param's **runtime type**, which is set by the launch and not by the
  config that declares it. `nextflow.config` declares `4.GB` (a `MemoryUnit`); a Seqera launch
  sends `'4 GB'` (a `String`), because JSON has no MemoryUnit and `-params-file` takes precedence.
  So the same commit works locally and fails in the cloud.

This harness therefore sets **every memory param as a String** — the exact shape a launch sends —
and includes the real `conf/base.config` rather than a copy.

## Running it

```sh
cd <repo root>
sbatch tests/config/run.sbatch
```

Expected output:

```
POWER_SIMULATION memory=8 GB
CONSOLIDATE_REPLICATES memory=46.5 GB
COMPUTE_POWER memory=46.8 GB
```

The input sizes are the moi5 trans sweep's, the largest the pipeline has run on: a 30 KB split,
1.90 GB of gzipped per-split output into consolidation, and the 1.40 GiB Parquet into
`compute_power`. Each expected number has a measurement behind it, all from
`docs/status.md`:

| Process | Measured peak | Requested here | Margin |
|---|---:|---:|---:|
| `POWER_SIMULATION` | 5.58 GiB | 8 GB | 1.43× |
| `CONSOLIDATE_REPLICATES` | 31.8 GiB | 46.5 GB | 1.46× |
| `COMPUTE_POWER` | 38.0 GiB | 46.8 GB | 1.23× |

Two failure modes to read for:

- **An exception** (`MissingMethodException`, `groovy.lang.MissingPropertyException`) means a
  closure does arithmetic on something it has not coerced, or references an input name a process
  does not have. That is the regression this exists to catch.
- **A number that moved** means a calibration changed. That is not automatically wrong, but it
  should be deliberate: update the table above and name the measurement, or put it back.

The tasks request 8 + 46.5 + 46.8 GB concurrently and the local executor refuses to start a task
larger than the node, which is why this runs as a job on a 128 GB allocation rather than on a login
node. It takes about a minute, nearly all of it Nextflow starting up.
