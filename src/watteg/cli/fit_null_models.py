"""Fit each gene's null model once per simulation, for every simulation task to reuse.

    watteg-fit-null-models --prepared prepared/ --reps 100 --seed 20250812 \
        --n-jobs 8 --out null_fits.h5

For each gene in the screen's pairs and each simulation, this draws the gene's counts with no
knockdown from its own stream, keyed `rng_for(seed, "__null_fit__|" + gene, simulation, 0.0)`,
and fits pysceptre's null model to them (`fit_all_genes`, at the discovery call's batch width).
`watteg-run-power-simulation --driver fast --null-fits reuse --null-fits-file null_fits.h5` then
tests every target against those fits instead of refitting each gene per target. See
`watteg/null_fits.py` for why that is sound and what it was measured to change.

A fit depends on (seed, gene, simulation) only, so a task given this file and a task that fits its
own genes produce the same bytes. The file records the seed, the baseline model and a digest of
`sim_input.h5`; a simulation run that does not match all three refuses it.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

from watteg.engine import AnalysisParams
from watteg.null_fits import compute_null_fits
from watteg.sim_input import read_sim_input
from watteg.workers import SHARED


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prepared", type=Path, required=True, help="output directory of watteg-prepare-sim-input"
    )
    parser.add_argument(
        "--pairs",
        type=Path,
        default=None,
        help="the pairs whose genes are fitted [default: pairs.tsv in --prepared]",
    )
    parser.add_argument("--reps", type=int, required=True, help="simulations to fit")
    parser.add_argument(
        "--rep-offset",
        type=int,
        default=0,
        help="simulations already covered, so this fits simulations offset+1 .. offset+reps",
    )
    parser.add_argument(
        "--seed", type=int, required=True, help="the simulation run's --seed; the file records it"
    )
    parser.add_argument(
        "--expression-model",
        choices=("fitted", "size_factor"),
        default="fitted",
        help="the simulation run's --expression-model; the file records it [default %(default)s]",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="worker processes; one simulation's fits are one unit of work [default %(default)s]",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.reps < 1:
        raise SystemExit("--reps must be at least 1")
    if args.rep_offset < 0:
        raise SystemExit("--rep-offset must be non-negative")

    started = time.perf_counter()
    params = AnalysisParams.from_analysis_mode(args.prepared / "analysis_mode.tsv")
    if params.resampling_mechanism != "permutations":
        # The fits are read only by the fast driver, which runs only the permutation test. Fitting
        # a CRT screen's genes would cost the time and be thrown away.
        raise SystemExit(
            f"this screen used {params.resampling_mechanism!r}; reused null fits are read only by "
            "the fast driver, which runs the permutation test only. Simulate it with --driver "
            "engine --null-fits refit, and skip this step."
        )
    sim = read_sim_input(args.prepared / "sim_input.h5")
    pairs = pd.read_csv(args.pairs or args.prepared / "pairs.tsv", sep="\t")
    if "response_id" not in pairs.columns:
        raise SystemExit(f"{args.pairs} has no 'response_id' column")
    genes = list(dict.fromkeys(pairs["response_id"]))
    unknown = [g for g in genes if g not in set(sim.genes)]
    if unknown:
        raise SystemExit(f"{len(unknown)} gene(s) are not in sim_input.h5, e.g. {unknown[:3]}")
    reps = range(args.rep_offset + 1, args.rep_offset + args.reps + 1)
    workers = min(max(args.n_jobs, 1), len(reps))
    print(
        f"fitting {len(genes)} genes x {len(reps)} simulations ({reps.start}-{reps.stop - 1}), "
        f"seed {args.seed}, baseline {args.expression_model}, on {workers} worker(s)"
    )

    SHARED.update(sim=sim, params=params, grna_csc=sim.grna_perts.tocsc())
    fits = compute_null_fits(
        sim,
        genes,
        reps,
        seed=args.seed,
        expression_model=args.expression_model,
        workers=workers,
        prepared=args.prepared,
    )
    fits.write(args.out)
    bad = {k: v for k, v in fits.summary().items() if v}
    print(
        f"wrote {len(genes)} x {len(reps)} fits to {args.out} "
        f"({args.out.stat().st_size / 1e6:.1f} MB) in {time.perf_counter() - started:.1f}s"
        + (f"; degenerate (gene, simulation) fits: {bad}" if bad else "")
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
