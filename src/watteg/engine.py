"""Running the screen's own test on simulated counts.

One call to pysceptre per (target, replicate). The alternative -- stacking
replicates as pseudo-genes against pseudo-targets, to amortise per-call setup --
was built and measured, and there is almost nothing to amortise: the per-target
binomial fit and CRT draws are about 1.5% of the work and the duplicate fits it
would remove about 0.8%. The per-pair test dominates, and no call shape touches
it. So this is the simple shape, and it stays simple: no pseudo-genes, no
pseudo-target keys, no replicate-chunk memory knob, and no question about
whether replicates share a resampling stream -- each call is its own target, so
each replicate draws its own.

**pysceptre's defaults are left alone.** `target_chunk_size` and
`chunk_memory_gb` are tuned, and the latter is documented as the fastest and
leanest setting measured. What is passed is not a default: `B1`, `B2`, `B3` and
the side are properties of the screen, read off its own analysis parameters, and
`n_jobs` is the caller's.

The inner entry point rather than `pipeline.api.run_discovery_analysis`, which
sizes `B2`/`B3` from `len(pairs)` the way R's `run_qc` does. The screen's own
budget is what the real analysis used and is what a simulation of it has to
reuse.
"""

from __future__ import annotations

import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .baseline import baseline_expression
from .perturbation import effect_size_matrix, guide_assignment, target_cells
from .seeds import SETUP_REP, derive_seed, rng_for
from .simulate import draw_counts

DEFAULT_GUIDE_SPREAD_C = 0.65
# How the simulations of one target get their permutations. "per-replicate" draws a fresh set for
# every simulation (the behaviour up to 2026-09-25). "per-target" gives every simulation of a target
# the same set, drawn once from the run's seed -- what the owner chose on 2026-09-25, the CLI's
# default since, and close to what sceptre does: its sampler reseeds mt19937(4) on every call, so
# R's simulations share one fixed set. This function's own defaults stay the engine's original
# configuration, which the tests compare against; the CLI passes what it was asked for.
PERMUTATION_MODES = ("per-replicate", "per-target")

# How the permutation nulls are computed. Results are identical either way; only speed differs.
# "scan" is pysceptre's default for permutations: PermutationPrefixSums gathers a (B, n_trt, 14)
# array and takes a running sum over it, to serve many targets of different sizes from one set of
# draws. WattEG's calls hold ONE target, so the running sum is thrown away but for one column
# (227 MB at stage 2, n_trt = 406). "sparse" is pysceptre's own draw-matrix route, which it already
# takes whenever the scan would be too large: a sparse indicator matrix times the gene's pieces.
# Measured 2026-09-25 (docs/pysceptre-backend.md, section 13, item 3a): 119.5 -> 17.7 ms per
# escalated pair at stage 2, 12.4 -> 1.9 ms at stage 1, 2.1x on cis and 2.3x on trans overall,
# with p-values, z-statistics and stages identical on 3,750 pair-tests. The CLI's default since.
NULL_ROUTES = ("sparse", "scan")


def _decline_scan(self, lo, hi, n_trt):
    """Stand-in for PermutationPrefixSums.statistics: always decline.

    pysceptre's documented contract for this method is that None means "too large to scan, use
    the draw matrix instead" (score_stat.py), so returning None selects its own sparse route
    without changing its source.
    """
    return None


@contextmanager
def null_route(route: str):
    """Compute the permutation nulls by `route` for the duration of the block."""
    if route not in NULL_ROUTES:
        raise ValueError(f"null route must be one of {NULL_ROUTES}, got {route!r}")
    if route == "scan":
        yield
        return
    from pysceptre.test_statistic import score_stat

    cls = getattr(score_stat, "PermutationPrefixSums", None)
    if cls is None or not hasattr(cls, "statistics"):
        raise RuntimeError(
            "pysceptre no longer has PermutationPrefixSums.statistics, so the sparse null route "
            "cannot be selected; check the pysceptre pin, or run with the scan route"
        )
    original = cls.statistics
    cls.statistics = _decline_scan
    try:
        yield
    finally:
        cls.statistics = original


@dataclass(frozen=True)
class AnalysisParams:
    """The screen's own test configuration, not this pipeline's choices."""

    B1: int
    B2: int
    B3: int
    side_code: int
    resampling_mechanism: str = "crt"

    @classmethod
    def from_analysis_mode(cls, path) -> AnalysisParams:
        text = Path(path).read_text()
        mode = dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)
        # The screen's own resampling mechanism is passed through, so the simulation re-runs
        # the test the screen actually ran. Until 2026-09-25 this refused anything but the CRT,
        # which ruled out both moi5 screens (sceptre's permutation test). pysceptre implements
        # permutations, but its own docs call them exercised rather than validated against R,
        # so a permutation sweep here leans on the R-vs-Python comparison (Stage B) for that.
        mechanism = mode.get("resampling_mechanism")
        if mechanism not in ("crt", "permutations"):
            raise ValueError(
                f"this screen used resampling_mechanism {mechanism!r}; the Python path covers "
                "'crt' and 'permutations'. Run it through the R implementation."
            )
        if mode.get("moi") == "low":
            raise ValueError(
                "this screen is low MOI, which pysceptre does not cover. Run it through the R "
                "implementation."
            )
        return cls(
            B1=int(mode["B1"]),
            B2=int(mode["B2"]),
            B3=int(mode["B3"]),
            side_code={"left": -1, "both": 0, "right": 1}[mode["side"]],
            resampling_mechanism=mechanism,
        )


def simulate_target(
    sim,
    target: str,
    genes: list[str],
    target_guides: list[str],
    *,
    effect_size: float,
    reps: range,
    seed: int,
    params: AnalysisParams,
    grna_csc,
    guide_spread_c: float = DEFAULT_GUIDE_SPREAD_C,
    n_jobs: int = 8,
    expression_model: str = "fitted",
    permutations: str = "per-replicate",
    nulls: str = "scan",
    estimand: str = "fixed",
) -> pd.DataFrame | None:
    """Simulate and test one target, returning one row per (pair, replicate).

    Returns None, with a warning, for a target that perturbs no cell.

    The guide assignment and the baseline are drawn once per target: which
    guide a cell carries is a fact about the screen, and the baseline depends
    on neither the effect size nor the draw. That is what `rep = 0` is for.
    """
    from pysceptre.pipeline.discovery import run_discovery_ntcells_complement

    is_perturbed = target_cells(sim.cre_perts, sim.target_ids, target)
    n_perturbed = int(is_perturbed.sum())
    if n_perturbed == 0:
        # Skipped, as R skips it, so a direct CLI run keeps the other targets'
        # rows instead of losing the whole split. In the pipeline this still
        # fails the task: POWER_SIMULATION checks that every pair produced every
        # replicate, and a skipped target's pairs produce none. That is the loud
        # outcome wanted there -- a pair with QC-passing counts but no perturbed
        # cell means the inputs disagree -- and the task log names the target.
        warnings.warn(f"skipping target {target!r}: it perturbs no cell", stacklevel=2)
        return None

    assignment = guide_assignment(
        grna_csc,
        sim.grna_ids,
        target_guides,
        is_perturbed,
        rng_for(seed, target, SETUP_REP, effect_size),
    )

    rows = [sim.genes.index(g) for g in genes]
    baseline = baseline_expression(
        expression_model,
        fitted_coefs=sim.fitted_coefs[rows],
        covariate_matrix=sim.covariate_matrix,
        mean=sim.row_data["mean"].to_numpy()[rows],
        size_factors=sim.col_data["size_factors"].to_numpy(),
    )
    theta = 1.0 / sim.row_data["dispersion"].to_numpy()[rows]
    # The config gives a fractional decrease; the simulation multiplies by a
    # relative expression level.
    wanted = np.full(len(genes), 1.0 - effect_size)

    treated = np.flatnonzero(is_perturbed)
    pairs = pd.DataFrame({"response_id": genes, "grna_target": target})

    # pysceptre draws its permutations from the seed it is given and the target's cell count, so
    # one seed per target means one permutation set per target. It comes from a stream spawned off
    # the setup seed, so it is independent of the guide assignment drawn from that seed and of
    # every replicate's counts: both modes simulate byte-identical counts.
    if permutations not in PERMUTATION_MODES:
        raise ValueError(f"permutations must be one of {PERMUTATION_MODES}, got {permutations!r}")
    target_permutation_seed = None
    if permutations == "per-target":
        stream = derive_seed(seed, target, SETUP_REP, effect_size).spawn(1)[0]
        target_permutation_seed = int(np.random.default_rng(stream).integers(0, 2**31 - 1))

    out = []
    for rep in reps:
        rng = rng_for(seed, target, rep, effect_size)
        counts = draw_counts(
            baseline,
            effect_size_matrix(assignment, wanted, guide_spread_c, rng, estimand=estimand),
            theta,
            rng,
        )
        with null_route(nulls):
            result = run_discovery_ntcells_complement(
                counts,
                genes,
                sim.covariate_matrix,
                {target: treated},
                pairs,
                B1=params.B1,
                B2=params.B2,
                B3=params.B3,
                side_code=params.side_code,
                resampling_mechanism=params.resampling_mechanism,
                seed=(
                    target_permutation_seed
                    if target_permutation_seed is not None
                    else int(rng.integers(0, 2**31 - 1))
                ),
                n_jobs=n_jobs,
                # Stated rather than left to be discovered. Each call carries exactly one target, so
                # the default of 200 is reduced to 1 by the memory budget every single time -- and
                # that reduction warns, six lines per call, which is 21,600 lines for a 36-target
                # 100-replicate run. Saying 1 here is not a tuning choice: it is what the value
                # already was. Chunk size affects only batching width, never a result.
                target_chunk_size=1,
            )
        result = result.assign(rep=rep, effect_size=effect_size, num_pert_cells=n_perturbed)
        out.append(result)

    frame = pd.concat(out, ignore_index=True)
    # log2 of the fold change, because that is the column the power step reads
    # and the one the R implementation writes. `fold_change < 1` and
    # `log_2_fold_change < 0` are the same predicate; carrying the log keeps
    # the two implementations' outputs comparable without a conversion step.
    frame["log_2_fold_change"] = np.log2(frame["fold_change"])
    return frame
