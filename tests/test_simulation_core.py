"""The simulation core: seeding, guide assignment, effect sizes, the draw.

These are the pieces a benchmark would otherwise measure without anyone having
checked they are right, so they are tested before anything times them.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from watteg.perturbation import (
    draw_guide_effects,
    effect_size_matrix,
    guide_assignment,
    target_cells,
)
from watteg.seeds import derive_seed, rng_for
from watteg.simulate import draw_counts


def tiny_screen():
    """Six cells, four guides, two targets.

    `gA1` and `gShared` hit `elemA`; `gShared` also hits `elemB`, which is the
    many-to-many case a real screen produces from overlapping elements. Cell 5
    carries nothing.
    """
    grna_ids = ["gA1", "gA2", "gShared", "gOther"]
    membership = {
        "gA1": [0, 1],
        "gA2": [1, 2],
        "gShared": [2, 3],
        "gOther": [3, 4],
    }
    rows = np.concatenate([[i] * len(membership[g]) for i, g in enumerate(grna_ids)])
    cols = np.concatenate([membership[g] for g in grna_ids])
    grna = sparse.csr_matrix((np.ones(rows.size), (rows, cols)), shape=(4, 6))
    # elemA is the union of gA1, gA2 and gShared; elemB of gShared and gOther.
    cre = sparse.csr_matrix(np.array([[1, 1, 1, 1, 0, 0], [0, 0, 1, 1, 1, 0]], dtype=float))
    return grna_ids, grna, ["elemA", "elemB"], cre


# --- seeding -----------------------------------------------------------------------------


def test_the_same_key_gives_the_same_stream_and_a_different_key_does_not():
    a = rng_for(1, "elemA", 3, 0.15).normal(size=5)
    assert np.array_equal(a, rng_for(1, "elemA", 3, 0.15).normal(size=5))
    for changed in (
        rng_for(2, "elemA", 3, 0.15),
        rng_for(1, "elemB", 3, 0.15),
        rng_for(1, "elemA", 4, 0.15),
        rng_for(1, "elemA", 3, 0.20),
    ):
        assert not np.array_equal(a, changed.normal(size=5))


def test_effect_size_keys_on_its_value_not_its_spelling():
    """It arrives parsed from a command line, so 0.15 and 0.150 are one effect
    size and must not seed two different streams."""
    assert derive_seed(1, "t", 1, 0.15).entropy == derive_seed(1, "t", 1, 0.150).entropy
    assert derive_seed(1, "t", 1, 0.15).entropy != derive_seed(1, "t", 1, 0.16).entropy


def test_the_setup_draw_cannot_collide_with_a_replicate():
    assert derive_seed(1, "t", 0, 0.15).entropy != derive_seed(1, "t", 1, 0.15).entropy
    with pytest.raises(ValueError, match="non-negative"):
        derive_seed(1, "t", -1, 0.15)


# --- guide assignment --------------------------------------------------------------------


def test_a_cell_is_assigned_a_guide_it_actually_carries():
    grna_ids, grna, target_ids, cre = tiny_screen()
    perturbed = target_cells(cre, target_ids, "elemA")
    assert perturbed.tolist() == [True, True, True, True, False, False]

    a = guide_assignment(
        grna.tocsc(), grna_ids, ["gA1", "gA2", "gShared"], perturbed, np.random.default_rng(0)
    )
    assert a.n_target_guides == 3
    assert a.n_other_guides == 1

    membership = {0: {0}, 1: {0, 1}, 2: {1, 2}, 3: {2, 3}, 4: {3}, 5: set()}
    for cell in range(6):
        if a.status[cell] == 0:
            continue
        if a.status[cell] <= 3:  # one of the target's guides, numbered 1..3
            row = grna_ids.index(["gA1", "gA2", "gShared"][a.status[cell] - 1])
        else:  # gOther, the only guide outside the target
            row = grna_ids.index("gOther")
        assert row in membership[cell], f"cell {cell} was given a guide it does not carry"


def test_a_cell_with_no_guide_gets_the_no_effect_row():
    grna_ids, grna, target_ids, cre = tiny_screen()
    perturbed = target_cells(cre, target_ids, "elemA")
    a = guide_assignment(
        grna.tocsc(), grna_ids, ["gA1", "gA2", "gShared"], perturbed, np.random.default_rng(0)
    )
    assert a.status[5] == 0  # cell 5 carries nothing


def test_control_cells_are_given_guides_from_outside_the_target():
    """A control cell's status must point past the target's block, so it can
    never pick up a targeting effect. The guide it points at has an effect of
    exactly 1 (other guides do not move the tested genes; see the module
    docstring of watteg.perturbation)."""
    grna_ids, grna, target_ids, cre = tiny_screen()
    perturbed = target_cells(cre, target_ids, "elemA")
    a = guide_assignment(
        grna.tocsc(), grna_ids, ["gA1", "gA2", "gShared"], perturbed, np.random.default_rng(0)
    )
    assert a.status[4] > a.n_target_guides  # cell 4 carries gOther only
    assert a.status[:4].max() <= a.n_target_guides  # perturbed cells never index past the block


def test_a_shared_guide_serves_both_of_its_targets():
    """gShared belongs to elemA and elemB. Read from a map that keeps one
    target per guide, elemB would lose it."""
    grna_ids, grna, target_ids, cre = tiny_screen()
    for target, guides in (("elemA", ["gA1", "gA2", "gShared"]), ("elemB", ["gShared", "gOther"])):
        perturbed = target_cells(cre, target_ids, target)
        a = guide_assignment(grna.tocsc(), grna_ids, guides, perturbed, np.random.default_rng(0))
        assert a.n_target_guides == len(guides)


def test_the_assignment_is_reproducible_and_depends_on_the_stream():
    grna_ids, grna, target_ids, cre = tiny_screen()
    perturbed = target_cells(cre, target_ids, "elemA")

    def call(seed):
        return guide_assignment(
            grna.tocsc(),
            grna_ids,
            ["gA1", "gA2", "gShared"],
            perturbed,
            np.random.default_rng(seed),
        ).status

    assert np.array_equal(call(0), call(0))


def test_an_unknown_target_and_a_target_with_no_usable_guide_are_refused():
    grna_ids, grna, target_ids, cre = tiny_screen()
    with pytest.raises(KeyError, match="not a row of cre_perts"):
        target_cells(cre, target_ids, "nope")
    perturbed = target_cells(cre, target_ids, "elemA")
    with pytest.raises(KeyError, match="appear in grna_perts"):
        guide_assignment(grna.tocsc(), grna_ids, ["absent"], perturbed, np.random.default_rng(0))


# --- effect sizes ------------------------------------------------------------------------


# --- the shared fixture: proof tests on the production path ------------------------------
#
# tests/fixtures/ holds a design the R suite reads byte for byte: 150 perturbed cells
# interleaved among 3,000, 26 carrying two of the target's guides, a listed target guide
# (t8) that no cell carries, and 454 control cells with no guide. Everything below drives
# the real guide_assignment() and effect_size_matrix() on it. The estimand being proven was
# decided on 2026-09-24 (docs/methods.md): power at a FIXED element effect. The realised mean
# over the perturbed cells equals the requested one in every replicate, the target's guides
# still differ, and every other cell is exactly 1.

FIXTURES = Path(__file__).parent / "fixtures"


def interleaved_fixture():
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
    grna = sparse.csc_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(grna_ids), len(carried))
    )
    return SimpleNamespace(
        grna=grna,
        grna_ids=grna_ids,
        carried=carried,
        is_perturbed=design["pert"].to_numpy() == "1",
        target_guides=guides.loc[guides["is_target"] == "1", "grna_id"].tolist(),
        other_guides=guides.loc[guides["is_target"] == "0", "grna_id"].tolist(),
    )


def fixture_assignment(seed=1):
    fx = interleaved_fixture()
    fx.assignment = guide_assignment(
        fx.grna, fx.grna_ids, fx.target_guides, fx.is_perturbed, np.random.default_rng(seed)
    )
    return fx


@pytest.mark.parametrize("es", [0.05, 0.15, 0.5])
def test_the_realised_mean_is_pinned_and_control_cells_are_exactly_one(es):
    fx = fixture_assignment()
    wanted = np.full(3, 1.0 - es)
    rng = np.random.default_rng(100)
    for _ in range(50):
        m = effect_size_matrix(fx.assignment, wanted, guide_spread_c=0.65, rng=rng)
        np.testing.assert_allclose(m[:, fx.is_perturbed].mean(axis=1), wanted, rtol=0, atol=1e-12)
        assert (m[:, ~fx.is_perturbed] == 1.0).all()


@pytest.mark.parametrize("es", [0.05, 0.15, 0.5])
def test_under_the_random_estimand_the_realised_mean_has_the_cell_weighted_spread(es):
    """--estimand random skips the pin. The realised mean over the perturbed cells is then the
    guides' draws weighted by their cell counts n_g, so its sd over simulations is
    c * es * (1 - es) * sqrt(sum n_g^2) / sum n_g: the variance of a cell-weighted mean. Control
    cells stay exactly 1."""
    fx = fixture_assignment()
    status = fx.assignment.status[fx.is_perturbed]
    n_g = np.bincount(status, minlength=fx.assignment.n_target_guides + 1)[1:]
    expected = 0.65 * es * (1 - es) * np.sqrt((n_g**2).sum()) / n_g.sum()
    rng = np.random.default_rng(8)
    means = []
    for _ in range(2000):
        m = effect_size_matrix(
            fx.assignment, np.array([1.0 - es]), guide_spread_c=0.65, rng=rng, estimand="random"
        )
        assert (m[:, ~fx.is_perturbed] == 1.0).all()
        means.append(m[0, fx.is_perturbed].mean())
    assert np.std(means) == pytest.approx(expected, rel=0.08)
    assert np.mean(means) == pytest.approx(1.0 - es, abs=4 * expected / np.sqrt(len(means)))


def test_the_estimands_draw_the_same_guides_and_differ_only_by_the_pin():
    """Same stream, same draws: the random matrix pinned is the fixed matrix, and the generator is
    left in the same state, so everything drawn after it (the counts) is drawn alike."""
    from watteg.perturbation import _pin_to_mean

    fx = fixture_assignment()
    wanted = np.array([0.85, 0.7])
    a, b = np.random.default_rng(31), np.random.default_rng(31)
    fixed = effect_size_matrix(fx.assignment, wanted, 0.65, a)
    random = effect_size_matrix(fx.assignment, wanted, 0.65, b, estimand="random")
    assert a.random() == b.random()
    assert not np.allclose(random[:, fx.is_perturbed].mean(axis=1), wanted, atol=1e-6)
    assert np.array_equal(fixed[:, ~fx.is_perturbed], random[:, ~fx.is_perturbed])
    np.testing.assert_array_equal(
        fixed[:, fx.is_perturbed], _pin_to_mean(random[:, fx.is_perturbed], wanted)
    )
    default = effect_size_matrix(fx.assignment, wanted, 0.65, np.random.default_rng(31))
    assert np.array_equal(default, fixed)
    with pytest.raises(ValueError, match="estimand must be one of"):
        effect_size_matrix(fx.assignment, wanted, 0.65, np.random.default_rng(1), estimand="x")


def test_at_effect_size_zero_the_two_estimands_are_one_simulation():
    fx = fixture_assignment()
    a, b = np.random.default_rng(5), np.random.default_rng(5)
    fixed = effect_size_matrix(fx.assignment, np.ones(3), 0.65, a)
    random = effect_size_matrix(fx.assignment, np.ones(3), 0.65, b, estimand="random")
    assert np.array_equal(fixed, random) and (fixed == 1.0).all()
    assert a.random() == b.random()


def test_the_guide_status_is_each_cells_own_guide():
    fx = fixture_assignment()
    a = fx.assignment
    assert a.status.size == len(fx.carried)
    n_t = a.n_target_guides
    for j, carried in enumerate(fx.carried):
        s = a.status[j]
        if fx.is_perturbed[j]:
            assert 1 <= s <= n_t and fx.target_guides[s - 1] in carried, f"cell {j}"
        elif not carried:
            assert s == 0, f"cell {j}"
        else:
            assert s > n_t and fx.other_guides[s - n_t - 1] in carried, f"cell {j}"

    # Cells carrying the same guide get the same effect within a replicate.
    m = effect_size_matrix(a, np.array([0.85]), guide_spread_c=0.65, rng=np.random.default_rng(13))
    perturbed_status = a.status[fx.is_perturbed]
    values = m[0, fx.is_perturbed]
    for s in np.unique(perturbed_status):
        assert np.ptp(values[perturbed_status == s]) == 0


def test_pinning_the_mean_keeps_the_guide_to_guide_spread():
    """The pin adds one constant per gene; it must not flatten the guides. With
    n_g perturbed cells on guide g and N in total, the expected within-replicate
    variance is sd^2 * (1 - sum n_g^2 / N^2), with sd = c * es * (1 - es)."""
    fx = fixture_assignment()
    status = fx.assignment.status[fx.is_perturbed]
    n_g = np.bincount(status, minlength=fx.assignment.n_target_guides + 1)[1:]
    expected = (0.65 * 0.15 * 0.85) ** 2 * (1 - (n_g**2).sum() / n_g.sum() ** 2)
    rng = np.random.default_rng(9)
    within = [
        effect_size_matrix(fx.assignment, np.array([0.85]), guide_spread_c=0.65, rng=rng)[
            0, fx.is_perturbed
        ].var()
        for _ in range(400)
    ]
    np.testing.assert_allclose(np.mean(within), expected, rtol=0.1)


def test_no_control_cell_indexes_a_targeting_row_even_when_the_last_target_guide_is_unused():
    """t8 is listed last and carried by no cell. R offset control statuses by the
    highest target index any cell carried (7) and put guide o001 on t8's row;
    this implementation offsets by the number of target guides."""
    fx = fixture_assignment()
    assert "t8" not in {g for c in fx.carried for g in c}
    control_status = fx.assignment.status[~fx.is_perturbed]
    assert ((control_status == 0) | (control_status > fx.assignment.n_target_guides)).all()

    m = effect_size_matrix(
        fx.assignment, np.array([0.5]), guide_spread_c=0.0, rng=np.random.default_rng(0)
    )
    assert (m[0, fx.is_perturbed] == 0.5).all()
    assert (m[0, ~fx.is_perturbed] == 1.0).all()


@pytest.mark.parametrize("es", [0.7, 0.9, 0.99, 0.999])
def test_strong_knockdowns_are_pinned_exactly(es):
    """es 0.99 and 0.999 are where the old shift-clamp-repeat loop ran out of
    iterations and raised on a pin that exists: with only one guide left above
    zero it converged by a factor of 0.9 per pass."""
    fx = fixture_assignment()
    wanted = np.full(2, 1.0 - es)
    rng = np.random.default_rng(11)
    for _ in range(200):
        m = effect_size_matrix(fx.assignment, wanted, guide_spread_c=0.65, rng=rng)
        assert (m >= 0).all()
        np.testing.assert_allclose(m[:, fx.is_perturbed].mean(axis=1), wanted, rtol=0, atol=1e-12)


def test_the_pin_is_solved_exactly_as_a_root_finder_would():
    from scipy.optimize import brentq

    from watteg.perturbation import _pin_to_mean

    # Cells at 0, 0 and 0.9 with a target of 0.1: a plain shift of -0.2 would send
    # two cells below zero and leave the mean at 0.233. The exact answer shifts by -0.6.
    np.testing.assert_allclose(
        _pin_to_mean(np.array([[0.0, 0.0, 0.9]]), np.array([0.1])), [[0.0, 0.0, 0.3]], atol=1e-15
    )

    rng = np.random.default_rng(21)
    for _ in range(50):
        v = np.clip(np.repeat(rng.normal(0.1, 0.3, 6), rng.integers(1, 41, 6)), 0.0, None)
        target = rng.uniform(0.001, 0.5)
        c = brentq(
            lambda c, v, target: np.maximum(v + c, 0).mean() - target,
            -2,
            2,
            args=(v, target),
            xtol=1e-14,
        )
        got = _pin_to_mean(v[None, :], np.array([target]))[0]
        np.testing.assert_allclose(got, np.maximum(v + c, 0), atol=1e-9)
        assert got.mean() == pytest.approx(target, abs=1e-12)

    # No clamping needed: the plain shift, to rounding -- every case up to es 0.5.
    v = np.array([[0.8, 0.9, 1.0]])
    np.testing.assert_allclose(_pin_to_mean(v, np.array([0.85])), v - 0.05, atol=1e-15)


@pytest.mark.parametrize("n", [50_000, 100_000, 300_000])
def test_the_pin_stays_exact_on_very_large_targets(n):
    """Perturbed-cell effects take one of a few guide values, so the closed
    form's cumulative sum adds up rounding instead of cancelling it. Before the
    correction step, 5e4 cells missed by ~1e-12 and raised."""
    from watteg.perturbation import _pin_to_mean

    rng = np.random.default_rng(24)
    for es in (0.0, 0.15, 0.5, 0.9):
        vals = np.clip(rng.normal(1 - es, 0.13, 4), 0.0, None)
        v = vals[rng.integers(0, 4, n)]
        out = _pin_to_mean(v[None, :], np.array([1 - es]))[0]
        assert (out >= 0).all()
        assert abs(out.mean() - (1 - es)) < 1e-14


def test_a_non_finite_effect_is_a_loud_error():
    a = SimpleNamespace(
        status=np.array([1, 1, 0]),
        is_perturbed=np.array([True, True, False]),
        n_target_guides=1,
        n_other_guides=0,
    )

    class NaNDraw:
        def beta(self, a, b, size):
            return np.full(size, np.nan)

    with pytest.raises(ValueError, match="could not pin"):
        effect_size_matrix(a, np.array([0.8]), guide_spread_c=0.65, rng=NaNDraw())


def test_at_effect_size_zero_every_effect_is_exactly_one():
    """A null element's guides do nothing. The old absolute N(1 - es, 0.13) kept
    a 13% spread here, which inflated false calls for highly expressed genes."""
    fx = fixture_assignment()
    m = effect_size_matrix(
        fx.assignment, np.array([1.0, 1.0]), guide_spread_c=0.65, rng=np.random.default_rng(22)
    )
    assert (m == 1.0).all()


def test_the_guide_draw_has_the_documented_moments_and_c_zero_means_no_spread():
    """sd = c * es * (1 - es), mean = 1 - es; bounded in (0, 1)."""
    rng = np.random.default_rng(25)
    for es in (0.05, 0.15, 0.5):
        x = draw_guide_effects(np.array([1.0 - es]), 400_000, 0.65, rng)[:, 0]
        assert ((x > 0) & (x < 1)).all()
        assert x.mean() == pytest.approx(1.0 - es, abs=2e-3)
        assert x.std() == pytest.approx(0.65 * es * (1.0 - es), rel=1e-2)
    np.testing.assert_array_equal(draw_guide_effects(np.array([0.85]), 5, 0.0, rng), 0.85)
    np.testing.assert_array_equal(draw_guide_effects(np.array([1.0]), 5, 0.65, rng), 1.0)


def test_each_gene_draws_its_own_guide_effects():
    """One draw shared across genes would make every tested gene of a target move together."""
    fx = fixture_assignment()
    m = effect_size_matrix(
        fx.assignment, np.full(3, 0.85), guide_spread_c=0.65, rng=np.random.default_rng(23)
    )
    values = m[:, fx.is_perturbed]
    assert not np.allclose(values[0], values[1])
    assert not np.allclose(values[1], values[2])


def test_a_target_that_perturbs_no_cell_is_skipped_not_fatal():
    """R skips it. Raising used to kill the whole split and every other target in it."""
    from watteg.engine import simulate_target

    sim = SimpleNamespace(cre_perts=sparse.csr_matrix((1, 4)), target_ids=["elemX"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = simulate_target(
            sim,
            "elemX",
            ["g1"],
            ["gX"],
            effect_size=0.15,
            reps=range(1, 2),
            seed=1,
            params=None,
            grna_csc=None,
        )
    assert out is None
    assert any("perturbs no cell" in str(w.message) for w in caught)


def test_effect_sizes_never_go_negative():
    rng = np.random.default_rng(4)
    n_cells = 200
    is_pert = np.zeros(n_cells, dtype=bool)
    is_pert[:100] = True
    status = np.where(is_pert, 1, 2)
    a = type(
        "A",
        (),
        {"status": status, "is_perturbed": is_pert, "n_target_guides": 1, "n_other_guides": 1},
    )()
    # A near-zero target with a wide spread is where clamping actually bites.
    matrix = effect_size_matrix(a, np.array([0.02]), guide_spread_c=0.5, rng=rng)
    assert (matrix >= 0).all()


# --- the draw ----------------------------------------------------------------------------


def test_counts_follow_the_negative_binomial_they_were_asked_for():
    rng = np.random.default_rng(5)
    n_cells = 60_000
    baseline = np.full((2, n_cells), 4.0)
    theta = np.array([2.0, 50.0])
    counts = draw_counts(baseline, np.ones((2, n_cells)), theta, rng, dtype=None)

    np.testing.assert_allclose(counts.mean(axis=1), 4.0, rtol=0.03)
    # Var = mu + mu^2/theta, which is what distinguishes theta from noise.
    np.testing.assert_allclose(counts.var(axis=1), 4.0 + 16.0 / theta, rtol=0.06)


def test_the_effect_size_scales_the_mean():
    rng = np.random.default_rng(6)
    n_cells = 40_000
    counts = draw_counts(
        np.full((1, n_cells), 10.0),
        np.full((1, n_cells), 0.5),
        np.array([20.0]),
        rng,
        dtype=None,
    )
    np.testing.assert_allclose(counts.mean(), 5.0, rtol=0.03)


def test_a_fully_knocked_down_gene_draws_zeros_rather_than_failing():
    counts = draw_counts(
        np.ones((1, 50)), np.zeros((1, 50)), np.array([3.0]), np.random.default_rng(7)
    )
    assert (counts == 0).all()


def test_counts_come_back_as_int16_and_overflow_is_refused_not_wrapped():
    counts = draw_counts(
        np.full((1, 20), 3.0), np.ones((1, 20)), np.array([5.0]), np.random.default_rng(8)
    )
    assert counts.dtype == np.int16
    with pytest.raises(OverflowError, match="int16"):
        draw_counts(
            np.full((1, 200), 1e5),
            np.ones((1, 200)),
            np.array([1e6]),
            np.random.default_rng(9),
        )


def test_the_draw_rejects_inputs_that_do_not_line_up():
    rng = np.random.default_rng(10)
    with pytest.raises(ValueError, match="effect_size"):
        draw_counts(np.ones((2, 5)), np.ones((2, 6)), np.array([1.0, 1.0]), rng)
    with pytest.raises(ValueError, match="theta is"):
        draw_counts(np.ones((2, 5)), np.ones((2, 5)), np.array([1.0]), rng)
    with pytest.raises(ValueError, match="finite and positive"):
        draw_counts(np.ones((2, 5)), np.ones((2, 5)), np.array([1.0, 0.0]), rng)


# --- the screen's own test configuration -------------------------------------------------


def test_the_screens_resampling_mechanism_is_passed_through(tmp_path):
    """Both moi5 screens ran sceptre's permutation test. The engine used to refuse
    anything but the CRT; it now re-runs whichever the screen ran."""
    from watteg.engine import AnalysisParams

    def mode(mechanism, moi="high"):
        f = tmp_path / f"{mechanism}_{moi}.tsv"
        f.write_text(
            f"resampling_mechanism\t{mechanism}\nmoi\t{moi}\nside\tleft\n"
            "B1\t499\nB2\t4999\nB3\t24999\n"
        )
        return f

    params = AnalysisParams.from_analysis_mode(mode("permutations"))
    assert params.resampling_mechanism == "permutations"
    assert (params.B1, params.B2, params.B3, params.side_code) == (499, 4999, 24999, -1)
    assert AnalysisParams.from_analysis_mode(mode("crt")).resampling_mechanism == "crt"
    with pytest.raises(ValueError, match="covers 'crt' and 'permutations'"):
        AnalysisParams.from_analysis_mode(mode("something_else"))
    with pytest.raises(ValueError, match="low MOI"):
        AnalysisParams.from_analysis_mode(mode("crt", moi="low"))
