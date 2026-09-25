"""Run the power simulation for one split, one effect size, one chunk of replicates.

    watteg-run-power-simulation --prepared prepared/ --pairs split_01.tsv \
        --effect-size 0.15 --reps 20 --rep-offset 0 --seed 1 --out sim.tsv

For each target in the split and each replicate, this simulates a count matrix
under the given effect size and asks the screen's own test whether it would have
called the association. The fraction of replicates in which it would is the
power, computed downstream.

Output columns match `src/run_power_simulation.R`'s, so the two implementations'
results can be compared without a conversion step and so downstream readers do
not care which produced a file.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from watteg.engine import (
    DEFAULT_GUIDE_SPREAD_C,
    NULL_ROUTES,
    PERMUTATION_MODES,
    AnalysisParams,
    simulate_target,
)
from watteg.sim_input import read_sim_input

# What a reader needs, and nothing more. At 100 replicates x 34,886 pairs x six
# effect sizes the columns nothing reads were 39% of a 3 GB output.
KEEP = [
    "grna_target",
    "response_id",
    "p_value",
    "log_2_fold_change",
    "rep",
    "effect_size",
    "num_pert_cells",
    "pass_qc",
    "n_nonzero_trt",
    "n_nonzero_cntrl",
]

# The bulk inputs, for the task-level pool. A worker reads them from here instead of
# receiving them: under fork a child inherits them without a copy, and only the small
# (target, replicate) descriptor crosses the process boundary.
_SHARED: dict = {}


def _limit_blas_threads() -> None:
    """One BLAS thread per worker, so N workers do not start N x cores threads."""
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=1, user_api="blas")
    except Exception:  # pragma: no cover - threadpoolctl comes with pysceptre
        pass


def _run_unit(unit: tuple) -> tuple:
    """Simulate and test one (target, replicate). pysceptre gets one worker: the task's
    workers are already busy with other units, and nesting pools would oversubscribe."""
    target, genes, guides, rep = unit
    shared = _SHARED
    at = time.perf_counter()
    frame = simulate_target(
        shared["sim"],
        target,
        genes,
        guides,
        effect_size=shared["effect_size"],
        reps=range(rep, rep + 1),
        seed=shared["seed"],
        params=shared["params"],
        grna_csc=shared["grna_csc"],
        guide_spread_c=shared["guide_spread_c"],
        n_jobs=1,
        expression_model=shared["expression_model"],
        permutations=shared["permutations"],
        nulls=shared["nulls"],
    )
    return target, frame, time.perf_counter() - at


def _load_shared(prepared: Path, settings: dict) -> None:
    sim = read_sim_input(prepared / "sim_input.h5")
    _SHARED.update(
        sim=sim,
        params=AnalysisParams.from_analysis_mode(prepared / "analysis_mode.tsv"),
        grna_csc=sim.grna_perts.tocsc(),
        **settings,
    )


def _init_spawned(prepared: Path, settings: dict) -> None:
    _limit_blas_threads()
    _load_shared(prepared, settings)


def _map_units(units: list, workers: int, prepared: Path, settings: dict) -> list:
    """Run the units in parallel, in PROCESSES, returning them in submission order.

    Never threads: pysceptre's discovery call keeps its working state in a module global,
    so two calls in one process overwrite each other (measured: KeyError 'permutations').
    On Linux the workers are forked and inherit `_SHARED` without a copy. Elsewhere they are
    spawned and each loads the inputs itself -- more memory, and only for local runs:
    forking after Apple's Accelerate BLAS can deadlock, the same reason pysceptre gives.

    Each unit draws from its own seeded stream (`rng_for(seed, target, rep, effect_size)`),
    so the output does not depend on the worker count or the order the units finish in.
    """
    if workers <= 1 or len(units) <= 1:
        return [_run_unit(u) for u in units]
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    if sys.platform.startswith("linux"):
        pool = ProcessPoolExecutor(
            max_workers=workers, mp_context=mp.get_context("fork"), initializer=_limit_blas_threads
        )
    else:
        pool = ProcessPoolExecutor(
            max_workers=workers,
            mp_context=mp.get_context("spawn"),
            initializer=_init_spawned,
            initargs=(prepared, settings),
        )
    # A worker that dies (the out-of-memory killer, or a crash as the pool starts) breaks the
    # whole pool. Python would then exit 1, which the pipeline reads as "our code is wrong, stop the
    # run" -- one bad machine would end a sweep. Exit 137 instead: the pipeline retries it, with
    # double the memory, which is right for the common cause and harmless for the rest. First seen
    # 2026-09-25 on 1 of 200 moi5 cis tasks, 26 s into the task.
    from concurrent.futures.process import BrokenProcessPool

    try:
        with pool:
            return list(pool.map(_run_unit, units))
    except BrokenProcessPool:
        print(
            "ERROR: a worker process was killed (most often out of memory). Exiting 137 so the "
            "task is retried with more memory.",
            file=sys.stderr,
        )
        raise SystemExit(137) from None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prepared", type=Path, required=True, help="output directory of watteg-prepare-sim-input"
    )
    parser.add_argument(
        "--pairs",
        type=Path,
        required=True,
        help="one split, with columns grna_target and response_id",
    )
    parser.add_argument(
        "--effect-size",
        type=float,
        required=True,
        help="a fractional decrease: 0.15 is a 15%% knockdown. 0 is the null arm. Power is "
        "reported at a FIXED element effect: in every replicate the realised mean knockdown "
        "across the perturbed cells equals this value exactly (see docs/methods.md)",
    )
    parser.add_argument("--reps", type=int, required=True)
    parser.add_argument(
        "--rep-offset",
        type=int,
        default=0,
        help="replicates already covered by earlier chunks, so `rep` stays unique across them",
    )
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="required: results are stochastic and must be reproducible",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--guide-spread-c",
        type=float,
        default=DEFAULT_GUIDE_SPREAD_C,
        help="guide-to-guide spread among the target's own guides: each guide's knockdown is "
        "Beta with mean --effect-size and sd c * es * (1 - es) [default %(default)s, fitted to "
        "per-guide data from three CRISPRi screens]. Zero at es = 0; must be in [0, 2). The "
        "guides' mean is pinned to --effect-size, so this adds no uncertainty about the "
        "element's effect. Other guides have no effect",
    )
    parser.add_argument(
        "--guide-sd",
        type=float,
        default=None,
        help="retired; use --guide-spread-c. Passing it is an error",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="workers for this task. Each (target, replicate) is one unit of work, so cis "
        "targets with a handful of genes still keep every worker busy [default %(default)s]",
    )
    parser.add_argument(
        "--permutations",
        choices=PERMUTATION_MODES,
        default="per-replicate",
        help="'per-target' tests every replicate of a target against one permutation set drawn "
        "from the seed, as sceptre effectively does; 'per-replicate' draws a fresh set for each "
        "replicate [default %(default)s]. The simulated counts are the same either way",
    )
    parser.add_argument(
        "--nulls",
        choices=NULL_ROUTES,
        default="scan",
        help="how the permutation nulls are computed: 'sparse' is pysceptre's draw-matrix route, "
        "about twice as fast for one-target calls; 'scan' is its default. The results are "
        "identical [default %(default)s]",
    )
    parser.add_argument(
        "--expression-model",
        choices=("fitted", "size_factor"),
        default="fitted",
        help="'fitted' draws from exp(X.beta), sceptre's own null model. "
        "'size_factor' reproduces the pre-2026-09-21 behaviour and exists "
        "only for comparison against sweeps produced with it",
    )
    args = parser.parse_args(argv)

    # Zero is allowed on purpose: it is the null arm, and running it is the only
    # way to measure this pipeline's false-call rate from its own output. Note what
    # it measures: a call needs p < threshold AND a negative fold change, so with a
    # two-sided p-value it estimates about alpha/2, not alpha.
    if not 0.0 <= args.effect_size < 1.0:
        raise SystemExit(
            f"--effect-size must be a fractional decrease in [0, 1) (got {args.effect_size}); "
            "0 is the null arm"
        )
    if args.reps < 1:
        raise SystemExit("--reps must be at least 1")
    if args.guide_sd is not None:
        raise SystemExit(
            "--guide-sd was replaced by --guide-spread-c on 2026-09-25: the spread is now "
            "c * es * (1 - es) (default 0.65), not an absolute sd. Passing the old 0.13 as c "
            "would shrink it fivefold, so it is refused rather than reinterpreted."
        )
    if not 0.0 <= args.guide_spread_c < 2.0:
        raise SystemExit(f"--guide-spread-c must be in [0, 2) (got {args.guide_spread_c})")

    started = time.perf_counter()
    sim = read_sim_input(args.prepared / "sim_input.h5")
    params = AnalysisParams.from_analysis_mode(args.prepared / "analysis_mode.tsv")
    design = pd.read_csv(args.prepared / "grna_targets.tsv", sep="\t")
    split = pd.read_csv(args.pairs, sep="\t")
    for column in ("grna_target", "response_id"):
        if column not in split.columns:
            raise SystemExit(f"{args.pairs} has no {column!r} column")

    # Many-to-many: a guide inside two overlapping elements belongs to both, so
    # this is read from the design table and never from a per-unit annotation.
    guides_of = design.groupby("grna_target")["grna_id"].apply(list)
    genes_of = split.groupby("grna_target")["response_id"].apply(list)
    reps = range(args.rep_offset + 1, args.rep_offset + args.reps + 1)

    print(
        f"{len(genes_of)} targets / {len(split)} pairs, replicates {reps.start}-{reps.stop - 1}, "
        f"effect size {args.effect_size} (relative expression {1 - args.effect_size:g})\n"
        f"  baseline: {args.expression_model}, permutations {args.permutations}, "
        f"nulls {args.nulls}, "
        f"n_jobs {args.n_jobs}, "
        f"B1/B2/B3 {params.B1}/{params.B2}/{params.B3}, side_code {params.side_code}"
    )

    for target in genes_of.index:
        if target not in guides_of:
            raise SystemExit(f"no gRNA maps to target {target!r} in grna_targets.tsv")

    # WHY THE UNIT IS (target, replicate) AND NOT target. pysceptre's own n_jobs
    # parallelises over the genes of one call, and a cis target has a median of 6 of
    # them, so a task given 8 cores left most of them idle. Replicates are
    # independent draws, so splitting them apart keeps every worker busy on cis and
    # on trans alike. The per-target setup (guide assignment, baseline) is redone
    # per unit; it is seeded, so it comes out identical, and it is under 1 % of a unit.
    settings = dict(
        effect_size=args.effect_size,
        seed=args.seed,
        guide_spread_c=args.guide_spread_c,
        expression_model=args.expression_model,
        permutations=args.permutations,
        nulls=args.nulls,
    )
    _SHARED.update(sim=sim, params=params, grna_csc=sim.grna_perts.tocsc(), **settings)
    units = [(t, list(g), guides_of[t], r) for t, g in genes_of.items() for r in reps]
    workers = min(max(args.n_jobs, 1), len(units))
    print(f"  {len(units)} (target, replicate) units on {workers} worker(s)")
    results = _map_units(units, workers, args.prepared, settings)

    by_target: dict[str, list] = {}
    for target, frame, elapsed in results:
        by_target.setdefault(target, []).append((frame, elapsed))
    frames = []
    for target, genes in genes_of.items():
        done = by_target[target]
        if all(frame is None for frame, _ in done):
            print(f"  {target}: skipped, no perturbed cells")
            continue
        frames.extend(frame for frame, _ in done if frame is not None)
        elapsed = sum(e for _, e in done)
        print(
            f"  {target}: {len(genes)} pairs in {elapsed:.1f}s of worker time "
            f"({elapsed / args.reps:.2f}s/replicate)"
        )

    # Every target in the split can have been skipped; an empty table with the
    # right columns keeps the row-count check downstream meaningful.
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=KEEP)
    for column in ("pass_qc", "n_nonzero_trt", "n_nonzero_cntrl"):
        if column not in combined:
            combined[column] = np.nan
    info = args.prepared / "pairs_with_info.tsv"
    if info.exists():
        # Real-data QC counts, constant across replicates, joined rather than
        # recomputed per draw -- which is what the R implementation does too.
        known = pd.read_csv(info, sep="\t")
        cols = [c for c in ("n_nonzero_trt", "n_nonzero_cntrl", "pass_qc") if c in known]
        combined = combined.drop(columns=cols).merge(
            known[["grna_target", "response_id", *cols]],
            on=["grna_target", "response_id"],
            how="left",
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined[[c for c in KEEP if c in combined]].to_csv(args.out, sep="\t", index=False)
    print(f"\nwrote {len(combined):,} rows to {args.out} in {time.perf_counter() - started:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
