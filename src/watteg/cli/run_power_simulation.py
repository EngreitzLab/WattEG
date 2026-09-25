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
from watteg.workers import SHARED as _SHARED
from watteg.workers import map_units

# "engine" calls pysceptre once per (target, simulation). "fast" (watteg.fast_driver) runs a
# target's simulations together so its permutation matrices are built once, with output
# byte-identical to the engine's under --permutations per-target; it needs that mode, since the
# matrices are shared only when the simulations share one permutation set.
DRIVERS = ("engine", "fast")

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

# How each gene's null model is fitted under the fast driver. "refit" fits it on every simulation's
# own counts, once per target, as the engine does. "reuse" fits it once per (gene, simulation) on an
# independent null draw and shares it across targets and effect sizes (watteg.null_fits).
NULL_FIT_MODES = ("reuse", "refit")


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


def _run_fast_unit(unit: tuple) -> tuple:
    """Simulate and test one target's simulations `start..stop-1` with the fast driver.

    Imported here, not at the top, so the default driver's start-up is exactly what it was.
    """
    from watteg.fast_driver import simulate_target_fast

    target, genes, guides, start, stop = unit
    shared = _SHARED
    at = time.perf_counter()
    frame = simulate_target_fast(
        shared["sim"],
        target,
        genes,
        guides,
        effect_size=shared["effect_size"],
        reps=range(start, stop),
        seed=shared["seed"],
        params=shared["params"],
        grna_csc=shared["grna_csc"],
        guide_spread_c=shared["guide_spread_c"],
        expression_model=shared["expression_model"],
        null_fits=shared.get("null_fits"),
    )
    return target, frame, time.perf_counter() - at


def _fast_units(genes_of, guides_of, reps: range, workers: int) -> list:
    """Each target's simulations in contiguous chunks, enough of them to keep `workers` busy.

    A chunk pays the per-target setup once (its permutations and their matrices, about 0.3 s), so
    chunks are as large as the worker count allows: one per target when targets outnumber workers,
    and a single-target run at --n-jobs 4 in four. Chunk boundaries cannot move a result: each
    simulation draws its counts from its own seeded stream, every simulation of a target shares the
    one permutation set, and every fit and product column is computed on its own.
    """
    per_target = min(len(reps), max(1, -(-workers // max(1, len(genes_of)))))
    size = -(-len(reps) // per_target)
    return [
        (t, list(g), guides_of[t], start, min(start + size, reps.stop))
        for t, g in genes_of.items()
        for start in range(reps.start, reps.stop, size)
    ]


def _map_units(units: list, workers: int, prepared: Path, settings: dict, fn=_run_unit) -> list:
    """Run the units in parallel, in processes, in submission order (see `watteg.workers`).

    Each unit draws from its own seeded stream (`rng_for(seed, target, rep, effect_size)`),
    so the output does not depend on the worker count or the order the units finish in.
    """
    return map_units(units, workers, prepared, settings, fn)


def _null_fits_for_task(args, sim, split: pd.DataFrame, reps: range):
    """The fits this task's genes need, from --null-fits-file or made here.

    Made once in the parent, before the simulation's workers start, so a gene tested against
    several of the task's targets is fitted once per simulation rather than once per worker.
    Either way only the task's genes x simulations are kept, which is what a spawned worker is
    sent.
    """
    from watteg.null_fits import NullFits, compute_null_fits, input_fingerprint

    genes = list(dict.fromkeys(split["response_id"]))
    at = time.perf_counter()
    if args.null_fits_file is not None:
        fits = NullFits.read(args.null_fits_file)
        try:
            fits.check_matches(
                seed=args.seed,
                expression_model=args.expression_model,
                fingerprint=input_fingerprint(sim),
                genes=genes,
                reps=reps,
            )
        except ValueError as err:
            raise SystemExit(f"{args.null_fits_file}: {err}") from None
        source = f"read from {args.null_fits_file}"
    else:
        workers = min(max(args.n_jobs, 1), len(reps))
        fits = compute_null_fits(
            sim,
            genes,
            reps,
            seed=args.seed,
            expression_model=args.expression_model,
            workers=workers,
            prepared=args.prepared,
        )
        source = f"fitted in this task on {workers} worker(s)"
    fits = fits.subset(genes, reps)
    bad = {k: v for k, v in fits.summary().items() if v}
    print(
        f"  null fits: {len(genes)} genes x {len(reps)} simulations, {source}, in "
        f"{time.perf_counter() - at:.1f}s" + (f"; degenerate: {bad}" if bad else "")
    )
    return fits


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
        "--driver",
        choices=DRIVERS,
        default="engine",
        help="'fast' tests all of a target's simulations together, building its permutation "
        "matrices once; the output is byte-identical to 'engine' with --permutations per-target "
        "--nulls sparse. It needs --permutations per-target and a permutation-test screen, and "
        "always computes the nulls by the sparse route [default %(default)s]",
    )
    parser.add_argument(
        "--null-fits",
        choices=NULL_FIT_MODES,
        default="refit",
        help="with --driver fast: 'reuse' fits each gene's null model once per simulation, on an "
        "independent draw with no knockdown, and uses that fit for every target the gene is "
        "tested against (R's FIT_NULL_MODELS approximation, docs/methods.md); 'refit' fits it on "
        "each simulation's own counts for every target, as the engine does [default %(default)s]",
    )
    parser.add_argument(
        "--null-fits-file",
        type=Path,
        default=None,
        help="with --null-fits reuse: the fits written by watteg-fit-null-models for this "
        "prepared input and --seed. Without it the task fits its own genes first, with the same "
        "keyed draws, so the values and the output are the same either way",
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
    if args.driver == "fast" and args.permutations != "per-target":
        raise SystemExit(
            "--driver fast needs --permutations per-target: it builds a target's permutation "
            "matrices once for all of its simulations, which is the engine's test only when the "
            "simulations share one permutation set. Add --permutations per-target, or use "
            "--driver engine."
        )

    if args.null_fits == "reuse" and args.driver != "fast":
        raise SystemExit(
            "--null-fits reuse needs --driver fast: the engine calls pysceptre's discovery "
            "entry point, which fits every gene itself. Use --null-fits refit with --driver engine."
        )
    if args.null_fits_file is not None and args.null_fits != "reuse":
        raise SystemExit("--null-fits-file is read only under --null-fits reuse")

    started = time.perf_counter()
    sim = read_sim_input(args.prepared / "sim_input.h5")
    params = AnalysisParams.from_analysis_mode(args.prepared / "analysis_mode.tsv")
    if args.driver == "fast" and params.resampling_mechanism != "permutations":
        raise SystemExit(
            f"--driver fast runs the permutation test only, and this screen used "
            f"{params.resampling_mechanism!r}. Use --driver engine."
        )
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
        + (
            f"nulls {args.nulls}, "
            if args.driver == "engine"
            else f"driver fast, nulls sparse (the fast driver's only route), null fits "
            f"{args.null_fits}, "
        )
        + f"n_jobs {args.n_jobs}, "
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
    if args.driver == "fast" and args.null_fits == "reuse":
        settings["null_fits"] = _null_fits_for_task(args, sim, split, reps)
        _SHARED["null_fits"] = settings["null_fits"]
    if args.driver == "fast":
        units = _fast_units(genes_of, guides_of, reps, max(args.n_jobs, 1))
        workers = min(max(args.n_jobs, 1), len(units))
        print(f"  {len(units)} (target, simulation chunk) units on {workers} worker(s)")
        results = _map_units(units, workers, args.prepared, settings, fn=_run_fast_unit)
    else:
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
