"""The same test as `engine.simulate_target`, laid out per target instead of per simulation.

Opt-in (`--driver fast`), and only for `--permutations per-target`: its output is byte-identical to
the engine's with `--permutations per-target --nulls sparse`, and it gets there by doing the same
arithmetic in a different order of loops, never by doing different arithmetic.

**What the engine repeats.** It calls pysceptre once per (target, simulation). Under per-target
permutations every one of a target's calls draws the same 30,497 permutations from the same seed,
builds the same (B, n_cells) indicator matrices from them, and rebuilds the stage-2 matrix again for
every gene that escalates (`PermutationSliceDraws` is deliberately not memoized). Each call also
rebuilds the covariate outer products and computes every gene's precomputation pieces twice, once
for diagnostics and once for the test. None of that depends on the simulation.

**What this does instead.** Per target: draw the permutations once, build each stage's matrix once
(stage 2 and 3 only when some gene reaches them), build the outer products once. Per simulation:
draw the counts exactly as the engine does, fit each gene exactly as pysceptre does (one gene per
fit, the same IRLS and theta calls), and compute its pieces once. Per gene: pysceptre's own
`run_low_level_test_full`, given the stage-1 nulls through its `null_statistics_fn`; a stage it asks
for beyond those comes from its own draw-matrix route, over the memoized matrices. So the p-value,
the stage and every other column of the row are pysceptre's.

**Fit reuse** (`null_fits`, `--null-fits reuse`) is the one departure from the engine's arithmetic,
and it is explicit: each gene's fit on the simulation's counts is replaced by its fit on that
simulation's null draw (`watteg.null_fits`), made once per (gene, simulation) and shared by every
target and effect size. Everything else above is unchanged, so with `null_fits=None` the output is
still the engine's, byte for byte.

**Why each step is exact.** The permutations are the same rows: one generator from the same seed,
consumed by the same `permutation_draws` call. The matrices come from pysceptre's own
`draws_for_target` and `draws_to_matrix`, so they hold the same indices in the same order, and the
nulls are the same `draws @ stacked` products. The fits and pieces are the same calls on the same
inputs at the same shapes.

**What it borrows from pysceptre's internals.** `_get_row`, `_THETA_BOUNDS`, `GenePrecomputation`,
`_warn_about_degenerate_gene_fits`, `_CI_Z`, `_as_pct` and `resolve_entropy` from
`pipeline.discovery`, so that the per-gene steps are pysceptre's own rather than a copy of them. Its
public entry points do not offer this call shape. A pysceptre update that changes any of them is
caught by `tests/test_fast_driver.py`, which compares against the engine on every test run.

**The wide product is there, and off.** With `stage2_columns > 1` the genes that escalate wait,
and their stage-2 nulls come from one product per batch, `P2 @ [stacked_1 | ... | stacked_K]`.
That is exact too: scipy's `csr_matvecs` accumulates each output column on its own, in stored
order, from zero, so a column of the wide product equals the narrow one bit for bit (checked at
real shapes, and by the tests at several widths). It was expected to be the main saving and
measured otherwise on the laptop (n = 131,055 cells, n_trt = 406). A process running alone takes
14 ms for a narrow stage-2 product, against 28 ms a gene at 16 wide; wide only wins at 64 genes
(13 ms, 1 GB of right-hand side). Our reading, not measured, is that one gene's 15 MB stays in
cache and a wide array does not until it is wide enough to stream. With six processes doing
nothing but products, wide wins, 47 -> 32 ms a gene at 16 wide. But end to end, with four
processes on the trans reference target, narrow ran at 88-93 ms a pair-test against 96-98 at 16
wide and 93-96 at 32, with 0.4-1.0 GB more memory each: the products are a fifth of the work, and
the fits around them compete for the same memory bandwidth. So the default is narrow. The width
stays a parameter for a machine family with less cache per core, to be measured there.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from pysceptre.crt.permutations import permutation_draws
from pysceptre.glm.irls import fit_poisson_glm_batch, x_outer_flat
from pysceptre.glm.nb_theta import estimate_theta
from pysceptre.pipeline import discovery as _discovery
from pysceptre.precompute.pieces import compute_precomputation_pieces
from pysceptre.test_statistic.empirical_p import compute_empirical_p_value
from pysceptre.test_statistic.resampling import P_THRESH, run_low_level_test_full
from pysceptre.test_statistic.score_stat import (
    PermutationSliceDraws,
    StagedDraws,
    compute_null_statistics_from_draws,
    compute_observed_full_statistic,
    stack_pieces,
    statistics_from_segment_sums,
)

from .baseline import baseline_expression
from .engine import DEFAULT_GUIDE_SPREAD_C
from .perturbation import effect_size_matrix, guide_assignment, target_cells
from .seeds import SETUP_REP, derive_seed, rng_for
from .simulate import draw_counts

# How many escalating genes share one stage-2 product; 1 means each gets its own (see the module
# docstring for why that is the default). Above 1, each waiting gene costs its slot in the wide
# array (14 x n_cells float64, 14.7 MB at 131,055 cells) plus the pieces and counts it keeps until
# the product (16.7 MB). No value changes a result, only speed and memory.
DEFAULT_STAGE2_COLUMNS = 1


class _SharedPermutationDraws(PermutationSliceDraws):
    """pysceptre's permutation draws with the stage memo put back.

    `PermutationSliceDraws` drops `StagedDraws`' memo because in a whole-screen run the stage
    matrices of ~1,000 escalating targets would pile up. Here one target is live at a time and every
    simulation of it reads the same stages, so the memo is what makes each stage's matrix a
    once-per-target cost. Same `_index_arrays`, same `draws_to_matrix`: the matrices are the ones
    pysceptre would build.
    """

    slice = StagedDraws.slice


@dataclass
class _Waiting:
    """A gene whose stage-1 p-value sends it to stage 2, kept until its batch is multiplied."""

    key: tuple[int, str]
    y: np.ndarray
    pieces: object
    known: dict


def _fit_gene(counts: np.ndarray, row: int, X: np.ndarray, xo: np.ndarray, dfr: int):
    """One gene's Poisson fit and dispersion, exactly as pysceptre's `_gene_fit_job` at width 1.

    Reproduced rather than called through `fit_all_genes` only to skip the precomputation pieces
    that job builds for its diagnostics and throws away: the test builds the same pieces again, so
    here they are built once and serve both. The shapes are kept too -- a (1, n) response, the
    shared outer products -- because the batched solve is exact only at a fixed width.
    """
    Y = np.empty((1, X.shape[0]))
    Y[0] = _discovery._get_row(counts, row)
    fit = fit_poisson_glm_batch(X, Y, X_outer_flat=xo)
    theta_est, method = estimate_theta(y=Y[0], mu=fit.fitted_values[0], dfr=dfr)
    lo, hi = _discovery._THETA_BOUNDS
    theta = max(min(theta_est, hi), lo)
    return fit, theta, theta_est, method


def _result_row(gene: str, target: str, result) -> dict:
    """pysceptre's result row for one pair (`_gene_job`), key for key."""
    half_width = _discovery._CI_Z * result.se_fold_change
    as_pct = _discovery._as_pct
    return {
        "response_id": gene,
        "grna_target": target,
        "p_value": result.p_value,
        "fold_change": result.fold_change,
        "se_fold_change": result.se_fold_change,
        "pct_change_es": as_pct(result.fold_change),
        "pct_change_es_ci_low": as_pct(result.fold_change - half_width),
        "pct_change_es_ci_high": as_pct(result.fold_change + half_width),
        "z_orig": result.z_orig,
        "stage": result.stage,
    }


def simulate_target_fast(
    sim,
    target: str,
    genes: list[str],
    target_guides: list[str],
    *,
    effect_size: float,
    reps: range,
    seed: int,
    params,
    grna_csc,
    guide_spread_c: float = DEFAULT_GUIDE_SPREAD_C,
    expression_model: str = "fitted",
    stage2_columns: int = DEFAULT_STAGE2_COLUMNS,
    null_fits=None,
    estimand: str = "fixed",
) -> pd.DataFrame | None:
    """`engine.simulate_target(..., permutations="per-target", nulls="sparse")`, faster.

    Same arguments, same return value, byte for byte: one row per (pair, simulation), in the
    engine's order (simulation-major, then `genes`), with the same columns and dtypes. Returns None,
    with a warning, for a target that perturbs no cell.

    `null_fits` (a `watteg.null_fits.NullFits` covering `genes` x `reps`) replaces each gene's fit
    on the simulation's counts with its fit on that simulation's null draw (`--null-fits reuse`).
    The test's score pieces still come from the simulation's own counts; only the coefficients and
    theta change, so this is the one option here that is not the engine's arithmetic. None, the
    default, refits as the engine does.
    """
    if params.resampling_mechanism != "permutations":
        # There is no CRT path here: the CRT draws per target from fitted probabilities, which this
        # driver never fits, so a CRT screen would silently be tested by permutations.
        raise ValueError(
            "the fast driver implements the permutation test only; this screen used "
            f"{params.resampling_mechanism!r}. Run it with the engine driver."
        )
    if stage2_columns < 1:
        raise ValueError(f"stage2_columns must be at least 1, got {stage2_columns}")

    # --- per target: the engine's setup, line for line -------------------------------------------
    is_perturbed = target_cells(sim.cre_perts, sim.target_ids, target)
    n_perturbed = int(is_perturbed.sum())
    if n_perturbed == 0:
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
    theta_true = 1.0 / sim.row_data["dispersion"].to_numpy()[rows]
    wanted = np.full(len(genes), 1.0 - effect_size)
    treated = np.flatnonzero(is_perturbed)
    stream = derive_seed(seed, target, SETUP_REP, effect_size).spawn(1)[0]
    target_permutation_seed = int(np.random.default_rng(stream).integers(0, 2**31 - 1))

    # --- per target: what the engine redoes in every simulation's call ---------------------------
    X = sim.covariate_matrix
    n_cells, n_cov = X.shape
    B1, B2, B3 = params.B1, params.B2, params.B3
    side = params.side_code
    # The call's permutations: one generator from the resolved seed, m = the largest target in the
    # call, which with one target is this one (discovery.py, `run_discovery_ntcells_complement`).
    n_trt = int(treated.size)
    perms = permutation_draws(
        n_cells,
        n_trt,
        B1 + B2 + B3,
        np.random.default_rng(_discovery.resolve_entropy(target_permutation_seed)),
    )
    draws = _SharedPermutationDraws(perms, n_trt, n_cells)
    stage1 = (0, B1)
    stage2 = (B1, B1 + B2)
    # Only a refit reads the outer products (151 MB at 131,055 cells and 12 covariates).
    xo = x_outer_flat(X) if null_fits is None else None
    dfr = n_cells - n_cov
    width = n_cov + 2
    # A repeated gene is fitted and tested on its LAST row, as pysceptre's dicts keyed by gene id
    # resolve it; its row is emitted once per occurrence.
    gene_row = {g: i for i, g in enumerate(genes)}
    unique_genes = list(dict.fromkeys(genes))

    finished: dict[tuple[int, str], dict] = {}
    waiting: list[_Waiting] = []
    wide = None

    def finish(key, y, pieces, stacked, known):
        # `known.get` over the exact (lo, hi) stage bounds run_low_level_test_full asks for; None
        # for anything else is pysceptre's documented "use the draw matrix", which here is the
        # memoized one, times this gene's own `stacked`.
        result = run_low_level_test_full(
            y=y,
            mu=pieces.mu,
            a=pieces.a,
            w=pieces.w,
            D=pieces.D,
            trt_idxs=treated,
            synthetic_idxs=draws,
            stacked=stacked,
            B1=B1,
            B2=B2,
            B3=B3,
            side_code=side,
            null_statistics_fn=lambda lo, hi: known.get((lo, hi)),
        )
        finished[key] = _result_row(key[1], target, result)

    def escalates(pieces, null1) -> bool:
        # The escalation rule of run_low_level_test_full, on the statistic it will compute, so this
        # decides only where the stage-2 nulls come from. Were it ever wrong, the test would ask for
        # a stage not in `known` and get it from the draw matrix: slower, same numbers.
        z = compute_observed_full_statistic(pieces.a, pieces.w, pieces.D, treated)
        return compute_empirical_p_value(null1, z, side) <= P_THRESH

    def multiply_waiting():
        k = len(waiting)
        # Only the last batch of a call is partial; its columns are copied into a contiguous array
        # because scipy would otherwise copy the strided view anyway.
        rhs = wide if k == stage2_columns else np.ascontiguousarray(wide[:, : k * width])
        sums = draws.slice(*stage2) @ rhs
        for j, gene in enumerate(waiting):
            # Contiguous, as `draws @ stacked` returns it, so the statistic's reductions run on the
            # same layout as pysceptre's.
            gene.known[stage2] = statistics_from_segment_sums(
                np.ascontiguousarray(sums[:, j * width : (j + 1) * width])
            )
            finish(gene.key, gene.y, gene.pieces, rhs[:, j * width : (j + 1) * width], gene.known)
        waiting.clear()

    # --- per simulation ---------------------------------------------------------------------------
    for pos, rep in enumerate(reps):
        rng = rng_for(seed, target, rep, effect_size)
        counts = draw_counts(
            baseline,
            effect_size_matrix(assignment, wanted, guide_spread_c, rng, estimand=estimand),
            theta_true,
            rng,
        )
        fits: dict[str, object] = {}
        for gene in unique_genes:
            row = gene_row[gene]
            y = _discovery._get_row(counts, row)
            if null_fits is None:
                fit, theta, theta_est, method = _fit_gene(counts, row, X, xo, dfr)
                pieces = compute_precomputation_pieces(y, X, fit.coefs[0], theta)
                lo, hi = _discovery._THETA_BOUNDS
                fits[gene] = _discovery.GenePrecomputation(
                    fitted_coefs=fit.coefs[0],
                    theta=theta,
                    theta_method=method,
                    theta_clamped=not (lo <= theta_est <= hi),
                    glm_converged=bool(fit.converged[0]),
                    min_eigenvalue=pieces.min_eigenvalue,
                    max_eigenvalue=pieces.max_eigenvalue,
                    n_covariates=pieces.D.shape[0],
                )
            else:
                # The null draw's fit, applied to this simulation's counts: pysceptre's gene job
                # does exactly this with the fit it is handed.
                reused = null_fits.get(gene, rep)
                pieces = compute_precomputation_pieces(y, X, reused.fitted_coefs, reused.theta)

            stacked = stack_pieces(pieces.a, pieces.w, pieces.D)
            known = {stage1: compute_null_statistics_from_draws(stacked, draws.slice(*stage1))}
            if stage2_columns == 1 or B2 <= 0 or not escalates(pieces, known[stage1]):
                finish((pos, gene), y, pieces, stacked, known)
                continue
            if wide is None:
                wide = np.empty((n_cells, stage2_columns * width))
            slot = len(waiting)
            wide[:, slot * width : (slot + 1) * width] = stacked
            waiting.append(_Waiting((pos, gene), y, pieces, known))
            del stacked
            if len(waiting) == stage2_columns:
                multiply_waiting()
        # One summary per simulation, as each of the engine's calls gives one. Reused fits were
        # summarised where they were made.
        if null_fits is None:
            _discovery._warn_about_degenerate_gene_fits(fits)
    if waiting:
        multiply_waiting()

    out = []
    for pos, rep in enumerate(reps):
        result = pd.DataFrame([finished[(pos, g)] for g in genes])
        out.append(result.assign(rep=rep, effect_size=effect_size, num_pert_cells=n_perturbed))
    frame = pd.concat(out, ignore_index=True)
    frame["log_2_fold_change"] = np.log2(frame["fold_change"])
    return frame
