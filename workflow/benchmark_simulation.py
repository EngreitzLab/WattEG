#!/usr/bin/env python3
"""Phase 3's gate: how much does the Python simulation actually cost?

    workflow/benchmark_simulation.py <prepared_dir> [--reps 100] [--chunk-sizes 10,25,50,100]

Two ways to drive pysceptre, timed against each other:

**per-rep** -- one `run_discovery_ntcells_complement` call per (target,
replicate), which is R's structure with a faster engine inside it. Every call
rebuilds the per-run state, including `x_outer_flat(covariate_matrix)`, which at
567,690 cells and p = 11 is not small.

**stacked** -- one call per (target, replicate chunk), with the replicates
stacked as pseudo-genes `gene@rep` against pseudo-targets `target@es@rep`. The
pseudo-target keys are not cosmetic: pysceptre seeds each target's resampling
stream from its *name*, so one shared key would hand every replicate the same
synthetic treated sets and correlate replicates the Wilson interval assumes are
independent. The effect size is in the key for the same reason -- one `seed`
across a sweep would otherwise share draws between effect sizes.

**Reported as four terms, not one wall clock**, because which one dominates is
what any later optimisation depends on:

1. drawing the counts;
2. the per-(gene, replicate) Poisson fit -- unavoidable, and the price of the
   faithful null model (`docs/pysceptre-backend.md` section 2.1);
3. the per-target binomial fit and CRT draws;
4. the redundancy in (3) from refitting identical cell sets per replicate,
   measured as (3) minus the same work for a single target.

Terms 2 and 3 are measured by calling `fit_all_genes` and `fit_all_targets`
**standalone**, on the same inputs the engine would hand them, rather than by
instrumenting pysceptre. So they are an attribution of where the work is, not a
decomposition of one call's clock: the engine overlaps a chunk's target fit
with the previous chunk's gene tests, so the parts can sum to more than the
whole. The remainder, total minus (2) minus (3), is the per-pair test.

The comparison point is R's measured cost model, `1.140s + 0.5561s x pairs` per
(target, simulation), and 634 CPU-h per effect size on day0.
"""

from __future__ import annotations

import argparse
import resource
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

from watteg.baseline import baseline_expression
from watteg.perturbation import effect_size_matrix, guide_assignment, target_cells
from watteg.seeds import SETUP_REP, rng_for
from watteg.sim_input import read_sim_input
from watteg.simulate import draw_counts

# day0's measured R cost model, for the only comparison that matters.
R_INTERCEPT, R_PER_PAIR = 1.140, 0.5561
R_CPU_HOURS = 634
R_TARGETS, R_PAIRS = 3026, 34886

TIMINGS: dict[str, float] = {}


@contextmanager
def timed(term: str):
    start = time.perf_counter()
    try:
        yield
    finally:
        TIMINGS[term] = TIMINGS.get(term, 0.0) + time.perf_counter() - start


def peak_rss_gb() -> float:
    """macOS reports bytes, Linux kilobytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 1024**3 if sys.platform == "darwin" else peak / 1024**2


def simulate_replicates(sim, target, genes, target_guides, effect_size, reps, seed, grna_csc):
    """Counts for one target across `reps` replicates, stacked as rows.

    The guide assignment and the baseline are drawn once per target, not per
    replicate: neither depends on the draw. That hoisting is in the R
    implementation too, and it is why the setup seed is `rep = 0`.
    """
    is_perturbed = target_cells(sim.cre_perts, sim.target_ids, target)
    setup_rng = rng_for(seed, target, SETUP_REP, effect_size)
    assignment = guide_assignment(grna_csc, sim.grna_ids, target_guides, is_perturbed, setup_rng)

    rows = [sim.genes.index(g) for g in genes]
    baseline = baseline_expression(
        "fitted",
        fitted_coefs=sim.fitted_coefs[rows],
        covariate_matrix=sim.covariate_matrix,
    )
    theta = 1.0 / sim.row_data["dispersion"].to_numpy()[rows]
    wanted = np.full(len(genes), 1.0 - effect_size)

    blocks = []
    for rep in range(1, reps + 1):
        rng = rng_for(seed, target, rep, effect_size)
        with timed("1. draw counts"):
            es = effect_size_matrix(assignment, wanted, guide_spread_c=0.65, rng=rng)
            blocks.append(draw_counts(baseline, es, theta, rng))
    return np.vstack(blocks), is_perturbed


def measure_terms(sim, target, genes, counts, is_perturbed, params, seed, reps, effect_size):
    """Where the time goes, by running each stage on its own.

    The redundancy term is the honest way to price the `target@es@rep` keys:
    every pseudo-target holds the same cell set, so all but one of their
    binomial fits is arithmetic the engine has already done. Fitting one target
    and `reps` of them measures exactly that difference.
    """
    from pysceptre.pipeline.discovery import fit_all_genes, fit_all_targets

    cells = np.flatnonzero(is_perturbed)
    gene_ids = [f"{g}@{r}" for r in range(1, reps + 1) for g in genes]

    with timed("2. per-gene Poisson fits"):
        fit_all_genes(
            counts,
            gene_ids,
            sim.covariate_matrix,
            gene_rows=list(range(len(gene_ids))),
            batch_width=1,
        )
    draws = dict(B1=params["B1"], B2=params["B2"], B3=params["B3"], seed=seed)
    with timed("3. target fits + CRT draws"):
        fit_all_targets(
            {f"{target}@{effect_size:g}@{r}": cells for r in range(1, reps + 1)},
            sim.covariate_matrix,
            **draws,
        )
    with timed("4. of which, one target's share"):
        fit_all_targets({target: cells}, sim.covariate_matrix, **draws)


def call_engine(counts, gene_ids, covariates, target_cells_map, pairs, params, seed):
    """Everything not named here is left at pysceptre's default.

    `target_chunk_size` and `chunk_memory_gb` in particular: those are tuned,
    and `chunk_memory_gb` is documented as the fastest and leanest setting
    measured. B1/B2/B3 and the side are not defaults -- they are properties of
    the screen, read off its own analysis parameters.

    The inner entry point rather than `pipeline.api.run_discovery_analysis`,
    which would size B2/B3 from `len(pairs)`; with replicates stacked the pair
    count is inflated by the replicate count. See the plan, section 2.4.
    """
    from pysceptre.pipeline.discovery import run_discovery_ntcells_complement

    return run_discovery_ntcells_complement(
        counts,
        gene_ids,
        covariates,
        target_cells_map,
        pairs,
        B1=params["B1"],
        B2=params["B2"],
        B3=params["B3"],
        side_code=params["side_code"],
        seed=seed,
        n_jobs=params["n_jobs"],
    )


def run_per_rep(sim, target, genes, counts, is_perturbed, params, seed, reps):
    """R's shape: one engine call per (target, replicate)."""
    cells = np.flatnonzero(is_perturbed)
    n_genes = len(genes)
    out = []
    for rep in range(reps):
        block = counts[rep * n_genes : (rep + 1) * n_genes]
        pairs = pd.DataFrame({"response_id": genes, "grna_target": target})
        out.append(
            call_engine(
                block, genes, sim.covariate_matrix, {target: cells}, pairs, params, seed + rep
            )
        )
    return pd.concat(out)


def run_stacked(sim, target, genes, counts, is_perturbed, params, seed, reps, effect_size):
    """One engine call per (target, replicate chunk)."""
    cells = np.flatnonzero(is_perturbed)
    gene_ids = [f"{g}@{rep}" for rep in range(1, reps + 1) for g in genes]
    target_map = {f"{target}@{effect_size:g}@{rep}": cells for rep in range(1, reps + 1)}
    pairs = pd.DataFrame(
        {
            "response_id": gene_ids,
            "grna_target": [
                f"{target}@{effect_size:g}@{rep}" for rep in range(1, reps + 1) for _ in genes
            ],
        }
    )
    return call_engine(counts, gene_ids, sim.covariate_matrix, target_map, pairs, params, seed)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path, help="output of watteg-prepare-sim-input")
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--chunk-sizes", type=str, default="10,25,50,100")
    parser.add_argument("--effect-size", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="workers for the per-gene tests. 8 is the measured knee on a 10-performance-core "
        "box: 2.5x wall clock for 1.5x CPU, with 14 slower than 8. It buys the thing you wait "
        "on and costs the thing you are billed for, so CPU-hour comparisons against R -- whose "
        "tasks were single-threaded -- should be read with that in mind.",
    )
    parser.add_argument(
        "--skip-per-rep", action="store_true", help="skip the per-rep shape, which is the slow one"
    )
    args = parser.parse_args()

    sim = read_sim_input(args.prepared / "sim_input.h5")
    pairs = pd.read_csv(args.prepared / "pairs.tsv", sep="\t")
    design = pd.read_csv(args.prepared / "grna_targets.tsv", sep="\t")
    mode = dict(
        line.split("\t", 1)
        for line in (args.prepared / "analysis_mode.tsv").read_text().splitlines()
        if "\t" in line
    )
    params = {
        "B1": int(mode.get("B1", 499)),
        "B2": int(mode.get("B2", 4999)),
        "B3": int(mode.get("B3", 0)),
        "side_code": {"left": -1, "both": 0, "right": 1}[mode.get("side", "both")],
        "n_jobs": args.n_jobs,
    }
    print(
        f"{sim.describe()}\n  B1/B2/B3 = {params['B1']}/{params['B2']}/{params['B3']}, "
        f"side_code = {params['side_code']}, n_jobs = {params['n_jobs']} "
        f"(everything else at pysceptre's default)"
    )

    # One small, one median, one large target, so the per-pair term is visible.
    by_target = pairs.groupby("grna_target")["response_id"].apply(list)
    sizes = by_target.apply(len).sort_values()
    chosen = [sizes.index[0], sizes.index[len(sizes) // 2], sizes.index[-1]]
    print("  targets: " + ", ".join(f"{t} ({len(by_target[t])} pairs)" for t in chosen))

    guides_of = design.groupby("grna_target")["grna_id"].apply(list)
    grna_csc = sim.grna_perts.tocsc()
    chunk_sizes = [int(c) for c in args.chunk_sizes.split(",")]

    print(
        f"\n{'target':<28} {'pairs':>5} {'shape':<10} {'chunk':>6} {'s/(tgt,rep)':>12} "
        f"{'peak GB':>8}"
    )
    results = []
    for target in chosen:
        genes = by_target[target]
        counts, is_perturbed = simulate_replicates(
            sim,
            target,
            genes,
            guides_of[target],
            args.effect_size,
            args.reps,
            args.seed,
            grna_csc,
        )
        n_trt = int(is_perturbed.sum())

        # Its own pass: run inside a shape's timing window it would be counted
        # as that shape's cost, which is what inflated the first chunk size on
        # the first attempt.
        stage_reps = min(chunk_sizes[0], args.reps)
        measure_terms(
            sim,
            target,
            genes,
            counts[: stage_reps * len(genes)],
            is_perturbed,
            params,
            args.seed,
            stage_reps,
            args.effect_size,
        )

        for shape in ("per-rep", "stacked"):
            if shape == "per-rep" and args.skip_per_rep:
                continue
            for chunk in chunk_sizes if shape == "stacked" else [1]:
                if chunk > args.reps:
                    continue
                start = time.perf_counter()
                done = 0
                while done < args.reps:
                    take = min(chunk, args.reps - done)
                    block = counts[done * len(genes) : (done + take) * len(genes)]
                    with timed("engine"):
                        if shape == "per-rep":
                            run_per_rep(
                                sim, target, genes, block, is_perturbed, params, args.seed, take
                            )
                        else:
                            run_stacked(
                                sim,
                                target,
                                genes,
                                block,
                                is_perturbed,
                                params,
                                args.seed,
                                take,
                                args.effect_size,
                            )
                    done += take
                per_unit = (time.perf_counter() - start) / args.reps
                print(
                    f"{target[:28]:<28} {len(genes):>5} {shape:<10} {chunk:>6} "
                    f"{per_unit:>12.3f} {peak_rss_gb():>8.2f}"
                )
                results.append((target, len(genes), shape, chunk, per_unit, n_trt))

    print("\n--- where the time goes, measured stage by stage ---")
    engine_total = TIMINGS.pop("engine", 0.0)
    for term, seconds in sorted(TIMINGS.items()):
        print(f"  {term:<32} {seconds:8.1f}s")
    genes_t = TIMINGS.get("2. per-gene Poisson fits", 0.0)
    targets_t = TIMINGS.get("3. target fits + CRT draws", 0.0)
    one_t = TIMINGS.get("4. of which, one target's share", 0.0)
    print(f"  {'5. per-pair test (remainder)':<32} {engine_total - genes_t - targets_t:8.1f}s")
    print(f"  {'(engine calls, total)':<32} {engine_total:8.1f}s")
    if targets_t:
        print(
            f"\n  The redundancy the target@es@rep keys cost: {targets_t - one_t:.1f}s of "
            f"{targets_t:.1f}s in term 3 ({100 * (targets_t - one_t) / targets_t:.0f}%), which is "
            f"{100 * (targets_t - one_t) / max(engine_total, 1e-9):.1f}% of engine time. That is "
            "the number section 5.3 of the plan would have to beat to be worth building."
        )

    print("\n--- against R ---")
    frame = pd.DataFrame(results, columns=["target", "pairs", "shape", "chunk", "sec", "n_trt"])
    for shape, group in frame.groupby("shape"):
        best = group.loc[group["sec"].idxmin()]
        # Fit the same shape of model R's was reported in, if there is more than one size.
        fit = np.polyfit(group["pairs"], group["sec"], 1) if group["pairs"].nunique() > 1 else None
        model = (
            f"{fit[1]:.3f}s + {fit[0]:.4f}s x pairs"
            if fit is not None
            else f"{best['sec']:.3f}s at {int(best['pairs'])} pairs"
        )
        predicted = (
            (fit[1] * R_TARGETS + fit[0] * R_PAIRS) * 100 / 3600
            if fit is not None
            else float("nan")
        )
        print(f"  {shape:<10} {model}")
        if fit is not None:
            print(
                f"{'':<12} -> {predicted:,.0f} CPU-h per effect size, against R's {R_CPU_HOURS} "
                f"({R_CPU_HOURS / max(predicted, 1e-9):.1f}x)"
            )
    print(f"  R          {R_INTERCEPT}s + {R_PER_PAIR}s x pairs -> {R_CPU_HOURS} CPU-h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
