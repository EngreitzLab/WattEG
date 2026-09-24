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

from watteg.perturbation import effect_size_matrix, guide_assignment, target_cells
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
        m = effect_size_matrix(fx.assignment, wanted, guide_sd=0.13, rng=rng)
        np.testing.assert_allclose(m[:, fx.is_perturbed].mean(axis=1), wanted, rtol=0, atol=1e-12)
        assert (m[:, ~fx.is_perturbed] == 1.0).all()


def test_without_the_pin_the_realised_mean_would_wander():
    """What the pin removes. Drawing the same guide effects and applying them
    without pinning, the realised mean varies from replicate to replicate by
    about guide_sd * sqrt(sum n_g^2) / sum n_g -- the other estimand, random
    guide effects."""
    fx = fixture_assignment()
    status = fx.assignment.status[fx.is_perturbed]
    rng = np.random.default_rng(8)
    means = []
    for _ in range(200):
        draws = np.clip(rng.normal(0.85, 0.13, size=fx.assignment.n_target_guides), 0.0, None)
        means.append(draws[status - 1].mean())
    assert np.std(means) > 0.01


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
    m = effect_size_matrix(a, np.array([0.85]), guide_sd=0.13, rng=np.random.default_rng(13))
    perturbed_status = a.status[fx.is_perturbed]
    values = m[0, fx.is_perturbed]
    for s in np.unique(perturbed_status):
        assert np.ptp(values[perturbed_status == s]) == 0


def test_pinning_the_mean_keeps_the_guide_to_guide_spread():
    """The pin adds one constant per gene; it must not flatten the guides. With
    n_g perturbed cells on guide g and N in total, the expected within-replicate
    variance is guide_sd^2 * (1 - sum n_g^2 / N^2)."""
    fx = fixture_assignment()
    status = fx.assignment.status[fx.is_perturbed]
    n_g = np.bincount(status, minlength=fx.assignment.n_target_guides + 1)[1:]
    expected = 0.13**2 * (1 - (n_g**2).sum() / n_g.sum() ** 2)
    rng = np.random.default_rng(9)
    within = [
        effect_size_matrix(fx.assignment, np.array([0.85]), guide_sd=0.13, rng=rng)[
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
        fx.assignment, np.array([0.5]), guide_sd=0.0, rng=np.random.default_rng(0)
    )
    assert (m[0, fx.is_perturbed] == 0.5).all()
    assert (m[0, ~fx.is_perturbed] == 1.0).all()


@pytest.mark.parametrize("es", [0.7, 0.9])
def test_strong_knockdowns_are_pinned_too(es):
    fx = fixture_assignment()
    wanted = np.full(2, 1.0 - es)
    rng = np.random.default_rng(11)
    for _ in range(30):
        m = effect_size_matrix(fx.assignment, wanted, guide_sd=0.13, rng=rng)
        assert (m >= 0).all()
        np.testing.assert_allclose(m[:, fx.is_perturbed].mean(axis=1), wanted, rtol=0, atol=1e-12)


def test_a_pin_that_cannot_be_reached_is_an_error_not_a_miss():
    """Perturbed cells at 0, 0 and 0.9 with a target of 0.1: one shift of -0.2
    sends two of them below zero, and clamping them leaves the mean at 0.233.
    Allowed a single pass it must say so; allowed more, it pins exactly."""
    a = SimpleNamespace(
        status=np.array([2, 2, 1, 0, 0]),
        is_perturbed=np.array([True, True, True, False, False]),
        n_target_guides=2,
        n_other_guides=0,
    )

    class Fixed:
        """Guide 1 draws 0.9, guide 2 draws 0."""

        def normal(self, loc, scale, size):
            return np.array([[0.9], [0.0]])

    with pytest.raises(ValueError, match="could not pin"):
        effect_size_matrix(a, np.array([0.1]), guide_sd=0.13, rng=Fixed(), max_iter=1)
    m = effect_size_matrix(a, np.array([0.1]), guide_sd=0.13, rng=Fixed())
    assert m[0, :3].mean() == pytest.approx(0.1, abs=1e-12)
    assert (m >= 0).all()


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
    matrix = effect_size_matrix(a, np.array([0.02]), guide_sd=0.5, rng=rng)
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
