"""The simulation core: seeding, guide assignment, effect sizes, the draw.

These are the pieces a benchmark would otherwise measure without anyone having
checked they are right, so they are tested before anything times them.
"""

from __future__ import annotations

import numpy as np
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
    """The reason this matters: a control cell's guide draws an effect size
    around 1 with the same spread. Sent to the no-effect row instead, the
    control arm keeps its mean and loses its variance."""
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


def test_each_arm_is_centred_on_what_it_should_average_to():
    """The property the centring step exists for: clamping negatives biases the
    mean up, so "a 15% knockdown" would otherwise be slightly weaker."""
    rng = np.random.default_rng(3)
    n_cells = 400
    status = rng.integers(0, 6, size=n_cells)
    is_pert = np.zeros(n_cells, dtype=bool)
    is_pert[: n_cells // 2] = True
    status[is_pert] = rng.integers(0, 3, size=is_pert.sum())  # 0..2: no guide or target guides
    status[~is_pert] = rng.integers(3, 6, size=(~is_pert).sum())
    a = type(
        "A",
        (),
        {"status": status, "is_perturbed": is_pert, "n_target_guides": 2, "n_other_guides": 3},
    )()

    wanted = np.array([0.85, 0.5, 0.95])
    matrix = effect_size_matrix(a, wanted, guide_sd=0.13, rng=rng)
    assert matrix.shape == (3, n_cells)
    np.testing.assert_allclose(matrix[:, is_pert].mean(axis=1), wanted)
    np.testing.assert_allclose(matrix[:, ~is_pert].mean(axis=1), 1.0)


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
