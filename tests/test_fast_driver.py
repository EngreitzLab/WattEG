"""The fast driver is a layout change only: the engine's rows, byte for byte.

`watteg.fast_driver` tests a target's simulations together so its permutation matrices are built
once. That is worth having only if nothing else moves, so every test here compares against the
engine run the fast driver claims to reproduce -- `--permutations per-target --nulls sparse` -- on
the shared fixture (tests/fixtures/), small enough to run everywhere and still reach all three of
the test's stages. The same claim on a real target is in `test_layout_invariance.py` (realdata).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from watteg.engine import AnalysisParams, simulate_target

pytest.importorskip("pysceptre")

from watteg.fast_driver import simulate_target_fast  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
# Small budgets so a stage-2 skew-normal fit is often rejected and stage 3 is reached; the real
# screen's 499/4999/24999 would make the fixture slow without exercising anything more.
PARAMS = AnalysisParams(B1=99, B2=299, B3=999, side_code=-1, resampling_mechanism="permutations")


def fixture_sim():
    """The interleaved design (3,000 cells, 150 perturbed) with six genes from 0.3 to 20 counts."""
    design = pd.read_csv(
        FIXTURES / "interleaved_design.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    guides = pd.read_csv(
        FIXTURES / "interleaved_guides.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    grna_ids = guides["grna_id"].tolist()
    index = {g: i for i, g in enumerate(grna_ids)}
    carried = [g.split(";") if g else [] for g in design["guides"]]
    rows = [index[g] for c in carried for g in c]
    cols = [j for j, c in enumerate(carried) for _ in c]
    n = len(carried)
    grna = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(grna_ids), n))
    is_perturbed = design["pert"].to_numpy() == "1"

    rng = np.random.default_rng(11)
    size_factors = rng.lognormal(0, 0.3, n)
    X = np.column_stack([np.ones(n), np.log(size_factors), rng.integers(0, 2, n).astype(float)])
    genes = [f"g{i}" for i in range(6)]
    coefs = np.column_stack(
        [np.log(np.geomspace(0.3, 20, 6)), np.full(6, 1.0), rng.normal(0, 0.2, 6)]
    )
    sim = SimpleNamespace(
        genes=genes,
        fitted_coefs=coefs,
        covariate_matrix=X,
        row_data=pd.DataFrame(
            {"mean": np.exp(coefs[:, 0]), "dispersion": np.linspace(0.1, 1.0, 6)}, index=genes
        ),
        col_data=pd.DataFrame({"size_factors": size_factors}),
        grna_ids=grna_ids,
        grna_perts=grna,
        target_ids=["elem1"],
        cre_perts=sparse.csr_matrix(is_perturbed.astype(float)[None, :]),
    )
    target_guides = guides.loc[guides["is_target"] == "1", "grna_id"].tolist()
    return sim, target_guides


def both(effect_size, genes=None, reps=range(1, 9), **fast):
    sim, target_guides = fixture_sim()
    genes = sim.genes if genes is None else genes
    kw = dict(
        effect_size=effect_size,
        reps=reps,
        seed=3,
        params=PARAMS,
        grna_csc=sim.grna_perts.tocsc(),
    )
    with warnings.catch_warnings():
        # Degenerate-fit summaries for the fixture's sparsest gene, from both drivers alike.
        warnings.simplefilter("ignore")
        engine = simulate_target(
            sim,
            "elem1",
            genes,
            target_guides,
            n_jobs=1,
            permutations="per-target",
            nulls="sparse",
            **kw,
        )
        faster = simulate_target_fast(sim, "elem1", genes, target_guides, **kw, **fast)
    return engine, faster


def assert_same_frame(engine, faster):
    assert list(engine.columns) == list(faster.columns)
    assert list(engine.dtypes) == list(faster.dtypes)
    pd.testing.assert_frame_equal(engine, faster, check_exact=True)
    assert engine.equals(faster)


@pytest.mark.parametrize("stage2_columns", [1, 2, 5])
def test_every_stage_matches_the_engine_at_every_width(stage2_columns):
    """At es 0.25 the fixture's pairs stop at all three stages, so each route to a p-value is
    compared: stage 1 alone, the stage-2 skew-normal fit on batched or narrow nulls, and stage 3
    on pysceptre's own draw-matrix route over the memoized matrices."""
    engine, faster = both(0.25, stage2_columns=stage2_columns)
    assert set(engine["stage"]) == {1, 2, 3}
    assert_same_frame(engine, faster)


@pytest.mark.parametrize("effect_size", [0.0, 0.5])
def test_the_null_and_a_strong_effect_match_the_engine(effect_size):
    engine, faster = both(effect_size, stage2_columns=3)
    assert_same_frame(engine, faster)


def test_a_repeated_gene_is_tested_as_the_engine_tests_it():
    """pysceptre keys its fits and rows by gene id, so a repeated gene is fitted on its last row
    and its row is written once per occurrence; the fast driver must resolve it the same way."""
    sim, _ = fixture_sim()
    genes = [sim.genes[0], sim.genes[3], sim.genes[0], sim.genes[5]]
    engine, faster = both(0.25, genes=genes, reps=range(4, 7))
    assert_same_frame(engine, faster)


def test_simulation_chunks_give_the_rows_of_one_run():
    """A task splits a target's simulations into chunks; the chunks must add up to one run."""
    _, whole = both(0.25, reps=range(1, 7))
    _, first = both(0.25, reps=range(1, 3))
    _, rest = both(0.25, reps=range(3, 7))
    assert_same_frame(whole, pd.concat([first, rest], ignore_index=True))


def test_columns_of_the_wide_product_equal_the_narrow_products():
    """The one numerical assumption the batched stage 2 adds, checked directly: scipy accumulates
    each column of `csr @ dense` on its own, so a gene's nulls do not depend on its batch."""
    from pysceptre.crt.permutations import draws_for_target, permutation_draws
    from pysceptre.test_statistic.score_stat import (
        compute_null_statistics_from_draws,
        draws_to_matrix,
        statistics_from_segment_sums,
    )

    rng = np.random.default_rng(5)
    n, m, B, width = 4000, 90, 300, 14
    perms = permutation_draws(n, m, B, np.random.default_rng(2))
    P = draws_to_matrix(draws_for_target(perms, m), n)
    blocks = [rng.standard_normal((n, width)) * rng.lognormal(0, 2, (n, 1)) for _ in range(7)]
    sums = P @ np.ascontiguousarray(np.hstack(blocks))
    for k, stacked in enumerate(blocks):
        batched = statistics_from_segment_sums(
            np.ascontiguousarray(sums[:, k * width : (k + 1) * width])
        )
        narrow = compute_null_statistics_from_draws(stacked, P)
        assert np.array_equal(batched, narrow, equal_nan=True)


def test_a_crt_screen_is_refused_rather_than_tested_by_permutations():
    sim, target_guides = fixture_sim()
    crt = AnalysisParams(B1=99, B2=299, B3=0, side_code=-1, resampling_mechanism="crt")
    with pytest.raises(ValueError, match="permutation test only"):
        simulate_target_fast(
            sim,
            "elem1",
            sim.genes,
            target_guides,
            effect_size=0.1,
            reps=range(1, 2),
            seed=1,
            params=crt,
            grna_csc=sim.grna_perts.tocsc(),
        )


def test_a_target_that_perturbs_no_cell_is_skipped_as_the_engine_skips_it():
    sim = SimpleNamespace(cre_perts=sparse.csr_matrix((1, 4)), target_ids=["elemX"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = simulate_target_fast(
            sim,
            "elemX",
            ["g1"],
            ["gX"],
            effect_size=0.15,
            reps=range(1, 2),
            seed=1,
            params=PARAMS,
            grna_csc=None,
        )
    assert out is None
    assert any("perturbs no cell" in str(w.message) for w in caught)


def test_fast_units_cover_every_simulation_once_and_in_order():
    from watteg.cli.run_power_simulation import _fast_units

    genes_of = pd.Series({"t1": ["a", "b"]})
    guides_of = pd.Series({"t1": ["x"]})
    units = _fast_units(genes_of, guides_of, range(11, 21), workers=4)
    assert [(u[3], u[4]) for u in units] == [(11, 14), (14, 17), (17, 20), (20, 21)]

    many = pd.Series({f"t{i}": ["a"] for i in range(6)})
    units = _fast_units(many, pd.Series({t: ["x"] for t in many.index}), range(1, 11), workers=4)
    assert [(u[0], u[3], u[4]) for u in units] == [(t, 1, 11) for t in many.index]


def test_the_cli_refuses_the_fast_driver_without_per_target_permutations(tmp_path):
    from watteg.cli.run_power_simulation import main

    args = ["--prepared", str(tmp_path), "--pairs", str(tmp_path / "p.tsv")]
    args += ["--effect-size", "0.1", "--reps", "2", "--seed", "1", "--out", str(tmp_path / "o")]
    with pytest.raises(SystemExit, match="needs --permutations per-target"):
        main([*args, "--driver", "fast"])


# --- fit reuse (--null-fits reuse) -------------------------------------------------------------
#
# The one option of the fast driver that is not the engine's arithmetic: each gene's null model is
# fitted once per simulation on an independent null draw and shared by every target. What must
# hold is that a fit depends on (seed, gene, simulation) alone, so where it was made -- a separate
# step, or the task itself, alongside whichever other genes -- cannot reach a result.


def reuse_run(fits, reps=range(1, 7), genes=None, effect_size=0.25):
    sim, target_guides = fixture_sim()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return simulate_target_fast(
            sim,
            "elem1",
            sim.genes if genes is None else genes,
            target_guides,
            effect_size=effect_size,
            reps=reps,
            seed=3,
            params=PARAMS,
            grna_csc=sim.grna_perts.tocsc(),
            null_fits=fits,
        )


def null_fits_for(genes, reps, seed=3):
    from watteg.null_fits import compute_null_fits

    sim, _ = fixture_sim()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_null_fits(sim, genes, reps, seed=seed, expression_model="fitted")


def test_a_null_fit_depends_on_its_gene_and_simulation_only():
    from watteg.null_fits import fit_one_simulation

    sim, _ = fixture_sim()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        together = fit_one_simulation(sim, sim.genes, 4, 3, "fitted")
        for gene in sim.genes:
            alone = fit_one_simulation(sim, [gene], 4, 3, "fitted")[gene]
            assert np.array_equal(alone.fitted_coefs, together[gene].fitted_coefs)
            assert alone.theta == together[gene].theta
        other_rep = fit_one_simulation(sim, sim.genes[:1], 5, 3, "fitted")[sim.genes[0]]
        other_seed = fit_one_simulation(sim, sim.genes[:1], 4, 4, "fitted")[sim.genes[0]]
    assert other_rep.theta != together[sim.genes[0]].theta
    assert other_seed.theta != together[sim.genes[0]].theta


def test_reused_fits_give_the_same_rows_wherever_they_were_made(tmp_path):
    """Fits for exactly the task's genes and simulations, and a subset of fits made for more genes
    and more simulations and read back from a file, give the same bytes; so do chunks."""
    from watteg.null_fits import NullFits

    sim, _ = fixture_sim()
    genes = [sim.genes[1], sim.genes[4], sim.genes[5]]
    own = reuse_run(null_fits_for(genes, range(1, 7)), genes=genes)
    null_fits_for(sim.genes, range(1, 10)).write(tmp_path / "fits.h5")
    wide = NullFits.read(tmp_path / "fits.h5")
    assert_same_frame(own, reuse_run(wide, genes=genes))
    first = reuse_run(wide, reps=range(1, 3), genes=genes)
    rest = reuse_run(wide, reps=range(3, 7), genes=genes)
    assert_same_frame(own, pd.concat([first, rest], ignore_index=True))


def test_reuse_changes_the_fit_and_nothing_else():
    """Against the refit on the same counts and permutations the p-values move, but only a
    little: the null draw's fit is close to the simulation's own."""
    _, refit = both(0.25, reps=range(1, 7))
    reused = reuse_run(null_fits_for(fixture_sim()[0].genes, range(1, 7)))
    assert list(refit.columns) == list(reused.columns)
    assert refit[["response_id", "rep", "num_pert_cells"]].equals(
        reused[["response_id", "rep", "num_pert_cells"]]
    )
    assert not refit["p_value"].equals(reused["p_value"])
    lp = np.log10(refit["p_value"]) - np.log10(reused["p_value"])
    assert np.median(np.abs(lp)) < 0.5


def test_null_fits_survive_a_file_and_refuse_a_mismatched_run(tmp_path):
    from watteg.null_fits import NullFits, input_fingerprint

    sim, _ = fixture_sim()
    fits = null_fits_for(sim.genes[:3], range(1, 4))
    fits.write(tmp_path / "f.h5")
    back = NullFits.read(tmp_path / "f.h5")
    for gene in sim.genes[:3]:
        for rep in range(1, 4):
            a, b = fits.get(gene, rep), back.get(gene, rep)
            assert np.array_equal(a.fitted_coefs, b.fitted_coefs) and a.theta == b.theta
            assert (a.theta_method, a.theta_clamped, a.glm_converged) == (
                b.theta_method,
                b.theta_clamped,
                b.glm_converged,
            )

    ok = dict(
        seed=3,
        expression_model="fitted",
        fingerprint=input_fingerprint(sim),
        genes=sim.genes[:3],
        reps=range(1, 4),
    )
    back.check_matches(**ok)
    for change, message in (
        (dict(seed=4), "seed 3"),
        (dict(expression_model="size_factor"), "baseline"),
        (dict(fingerprint="0" * 32), "different sim_input"),
        (dict(genes=sim.genes[:4]), "have no fit"),
        (dict(reps=range(1, 5)), "simulation"),
    ):
        with pytest.raises(ValueError, match=message):
            back.check_matches(**{**ok, **change})


def test_the_cli_refuses_reuse_without_the_fast_driver(tmp_path):
    from watteg.cli.run_power_simulation import main

    args = ["--prepared", str(tmp_path), "--pairs", str(tmp_path / "p.tsv")]
    args += ["--effect-size", "0.1", "--reps", "2", "--seed", "1", "--out", str(tmp_path / "o")]
    with pytest.raises(SystemExit, match="needs --driver fast"):
        main([*args, "--driver", "engine", "--null-fits", "reuse"])
    with pytest.raises(SystemExit, match="read only under --null-fits reuse"):
        main([*args, "--null-fits", "refit", "--null-fits-file", str(tmp_path / "f.h5")])
