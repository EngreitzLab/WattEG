"""Power, its interval, and the minimum detectable effect size.

Every one of these has been a source of silent, wrong output in the R
implementation's history -- see `docs/status.md`, "Correctness fixes" -- so they
are tested for the specific failure rather than only for the happy path.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from watteg.power import (
    compute_power,
    effect_label,
    min_detectable_effect_size,
    wilson_interval,
)


def simulations(p_values, log2fcs, *, target="elemA", gene="geneX", effect_size=0.15):
    return pd.DataFrame(
        {
            "grna_target": target,
            "response_id": gene,
            "p_value": p_values,
            "log_2_fold_change": log2fcs,
            "rep": range(1, len(p_values) + 1),
            "effect_size": effect_size,
        }
    )


# --- the interval ------------------------------------------------------------------------


def test_the_interval_does_not_collapse_at_the_boundary():
    """The reason Wilson is used at all. `power_ci_low` is the column
    thresholded to establish a negative, so an interval that collapses would
    support claims the data cannot."""
    low, high = wilson_interval(np.array([0]), np.array([100]))
    # Not exactly 0: `center - halfwidth` cancels to a denormal rather than to
    # zero, and R's does the same. What matters is that the interval has width.
    assert low[0] == pytest.approx(0.0, abs=1e-15)
    assert 0.03 < high[0] < 0.04  # the normal approximation gives [0, 0]

    low, high = wilson_interval(np.array([100]), np.array([100]))
    assert high[0] == pytest.approx(1.0, abs=1e-15)
    assert 0.96 < low[0] < 0.97


def test_the_interval_brackets_the_estimate_and_stays_in_range():
    successes = np.array([0, 1, 25, 50, 99, 100])
    n = np.full(6, 100)
    low, high = wilson_interval(successes, n)
    # The same denormal cancellation as above, so bracketing is asserted to
    # floating-point rather than exactly.
    assert (low <= successes / n + 1e-15).all()
    assert (successes / n <= high + 1e-15).all()
    assert (low >= 0).all() and (high <= 1).all()


# --- the column label --------------------------------------------------------------------


def test_effect_label_does_not_trim_a_whole_number_into_a_different_one():
    """0.2 becoming "2" shipped once: a 20% knockdown reported under
    `power_at_effect_size_2`."""
    assert effect_label(0.2) == "20"
    assert effect_label(0.5) == "50"
    assert effect_label(0.05) == "5"
    assert effect_label(0.15) == "15"
    assert effect_label(0.125) == "12_5"  # underscore, because dots separate names downstream


# --- power -------------------------------------------------------------------------------


def test_a_replicate_counts_only_if_expression_went_down():
    """One-sided on purpose: a significant p-value with expression going UP is
    not a detected knockdown."""
    frame = simulations([1e-6, 1e-6, 1e-6, 1e-6], [-0.3, -0.3, +0.3, +0.3])
    power = compute_power(frame, threshold=1e-3)
    assert power["power"].iloc[0] == 0.5


def test_replicates_with_no_estimate_are_dropped_and_counted():
    frame = simulations([1e-6, np.nan, 1e-6], [-0.3, -0.3, np.nan])
    power = compute_power(frame, threshold=1e-3)
    assert power["n_reps"].iloc[0] == 1
    assert power["power"].iloc[0] == 1.0


def test_duplicated_replicates_are_refused_rather_than_double_counted():
    """Overlapping replicate offsets across chunks is the way this happens, and
    it would inflate the denominator silently."""
    frame = simulations([1e-6, 1e-6], [-0.3, -0.3])
    frame.loc[1, "rep"] = 1
    with pytest.raises(ValueError, match="duplicated"):
        compute_power(frame, threshold=1e-3)


def test_mixing_effect_sizes_is_refused():
    """Averaging power across knockdown levels is meaningless, and nothing
    downstream would reveal it."""
    frame = simulations([1e-6, 1e-6], [-0.3, -0.3])
    frame.loc[1, "effect_size"] = 0.5
    with pytest.raises(ValueError, match="mixes effect sizes"):
        compute_power(frame, threshold=1e-3)


def test_an_unusable_threshold_is_refused():
    frame = simulations([1e-6], [-0.3])
    for bad in (np.nan, -np.inf, 0.0, 1.5):
        with pytest.raises(ValueError, match="usable probability"):
            compute_power(frame, threshold=bad)


def test_pairs_are_kept_separate():
    a = simulations([1e-6] * 4, [-0.3] * 4, gene="geneA")
    b = simulations([0.5] * 4, [-0.3] * 4, gene="geneB")
    power = compute_power(pd.concat([a, b]), threshold=1e-3)
    assert len(power) == 2
    assert dict(zip(power["response_id"], power["power"], strict=True)) == {
        "geneA": 1.0,
        "geneB": 0.0,
    }


# --- minimum detectable effect size ------------------------------------------------------

GRID = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5]


def test_a_pair_must_clear_the_threshold_at_every_larger_effect_size():
    """The suffix rule. Taking the first effect size that clears lets Monte
    Carlo noise win, and the error runs one way -- always reporting the pair as
    more detectable than it is."""
    # Clears at 0.1, fails at 0.15, clears again above. Only the suffix counts.
    power = np.array([[0.2, 0.9, 0.5, 0.95, 0.99, 1.0]])
    assert min_detectable_effect_size(power, GRID, 0.8)[0] == 0.2


def test_the_answer_is_the_lowest_effect_size_of_an_unbroken_run():
    power = np.array([[0.1, 0.3, 0.85, 0.9, 0.95, 1.0]])
    assert min_detectable_effect_size(power, GRID, 0.8)[0] == 0.15


def test_never_clearing_gives_nan_rather_than_the_largest_effect_size():
    """NaN means no tested effect size qualified -- a statement about the
    effect sizes that were run, not proof the pair is undetectable."""
    power = np.array([[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]])
    assert np.isnan(min_detectable_effect_size(power, GRID, 0.8)[0])


def test_untested_effect_sizes_are_unknown_and_cannot_block_a_suffix():
    power = np.array([[np.nan, np.nan, 0.85, np.nan, 0.95, 1.0]])
    assert min_detectable_effect_size(power, GRID, 0.8)[0] == 0.15


def test_the_interval_on_power_inverts_into_the_interval_on_effect_size():
    """Power is monotone in effect size, so thresholding a LOWER bound on power
    yields a LARGER effect size. Getting this backwards would report the
    conservative edge as the optimistic one."""
    point = np.array([[0.5, 0.75, 0.85, 0.9, 0.95, 1.0]])
    ci_low = np.array([[0.4, 0.65, 0.75, 0.88, 0.93, 0.99]])
    ci_high = np.array([[0.6, 0.85, 0.95, 0.99, 1.0, 1.0]])

    estimate = min_detectable_effect_size(point, GRID, 0.8)[0]
    optimistic = min_detectable_effect_size(ci_high, GRID, 0.8)[0]
    conservative = min_detectable_effect_size(ci_low, GRID, 0.8)[0]
    assert optimistic <= estimate <= conservative
    assert (optimistic, estimate, conservative) == (0.1, 0.15, 0.2)


def test_the_estimand_is_carried_and_mixing_estimands_is_refused():
    fixed = simulations([1e-5, 0.5], [-0.2, -0.1]).assign(estimand="fixed")
    power = compute_power(fixed, threshold=1e-3)
    assert power["estimand"].tolist() == ["fixed"]
    mixed = pd.concat([fixed, simulations([1e-5], [-0.2], gene="geneY").assign(estimand="random")])
    with pytest.raises(ValueError, match="mixes estimands"):
        compute_power(mixed, threshold=1e-3)


# --- per-pair counts inside the task ---------------------------------------------------------
#
# A simulation task that holds all simulations of its pairs writes per-pair counts instead of one
# row per simulation, and COMPUTE_POWER adds them up. The claim is that this gives the table a
# consolidated per-simulation file gives, byte for byte.


def sweep(seed=4, n_targets=5, genes_per_target=4, n_sims=30):
    """A small sweep, in the order a task writes it: target by target, simulation-major."""
    rng = np.random.default_rng(seed)
    frames = []
    for t in range(n_targets):
        genes = [f"g{t}_{i}" for i in range(genes_per_target)]
        cells = int(rng.integers(50, 600))
        for rep in range(1, n_sims + 1):
            frames.append(
                pd.DataFrame(
                    {
                        "grna_target": f"t{t}",
                        "response_id": genes,
                        "p_value": 10 ** rng.uniform(-8, 0, len(genes)),
                        "log_2_fold_change": rng.normal(-0.2, 0.15, len(genes)),
                        "rep": rep,
                        "effect_size": 0.15,
                        "num_pert_cells": cells,
                        "estimand": "fixed",
                    }
                )
            )
    frame = pd.concat(frames, ignore_index=True)
    frame.loc[rng.random(len(frame)) < 0.02, "p_value"] = np.nan
    return frame


def test_counts_from_each_task_give_the_table_the_rows_give():
    from watteg.power import merge_power_counts, power_counts, power_from_counts

    rows = sweep()
    expected = compute_power(rows, threshold=1e-3)
    by_task = [power_counts(part, 1e-3) for _, part in rows.groupby("grna_target", sort=False)]
    merged = power_from_counts(merge_power_counts(by_task))
    assert list(merged.columns) == list(expected.columns)
    assert list(merged.dtypes) == list(expected.dtypes)
    assert merged.equals(expected)


def test_counts_from_simulation_chunks_add_up_and_overlaps_are_refused():
    from watteg.power import merge_power_counts, power_counts, power_from_counts

    rows = sweep()
    expected = compute_power(rows, threshold=1e-3)
    early, late = rows[rows["rep"] <= 12], rows[rows["rep"] > 12]
    merged = power_from_counts(
        merge_power_counts([power_counts(early, 1e-3), power_counts(late, 1e-3)])
    )
    exact = ["grna_target", "response_id", "power", "power_ci_low", "power_ci_high", "n_reps"]
    assert merged[exact].equals(expected[exact])
    np.testing.assert_allclose(
        merged["mean_log_2_fold_change"], expected["mean_log_2_fold_change"], rtol=1e-12
    )

    with pytest.raises(ValueError, match="overlapping simulation ranges"):
        merge_power_counts([power_counts(rows, 1e-3), power_counts(late, 1e-3)])
    with pytest.raises(ValueError, match="mixes thresholds"):
        merge_power_counts([power_counts(early, 1e-3), power_counts(late, 1e-4)])


def test_a_pair_with_no_usable_simulation_is_counted_and_left_out_of_the_table():
    from watteg.power import power_counts, power_from_counts

    rows = sweep()
    rows.loc[rows["response_id"] == "g2_1", "p_value"] = np.nan
    counts = power_counts(rows, 1e-3)
    lost = counts[counts["response_id"] == "g2_1"].iloc[0]
    assert (lost["n_reps"], lost["n_simulations"], lost["successes"]) == (0, 30, 0)
    table = power_from_counts(counts)
    assert "g2_1" not in set(table["response_id"])
    assert table.equals(compute_power(rows, threshold=1e-3))


def test_the_cli_gives_the_same_bytes_from_counts_and_from_rows(tmp_path):
    """End to end through the files: each task's rows as TSV and its counts as written by the
    simulation CLI, then watteg-compute-power on either. The TSV floats are read back exactly."""
    from watteg.cli.compute_power import main as compute
    from watteg.cli.run_power_simulation import _write_partials

    threshold = tmp_path / "threshold.txt"
    threshold.write_text("0.00064839688355153302\n")
    rows = sweep(seed=9, n_sims=40)
    (tmp_path / "sims").mkdir()
    (tmp_path / "partials").mkdir()
    for target, part in rows.groupby("grna_target", sort=False):
        part.to_csv(tmp_path / "sims" / f"{target}.tsv", sep="\t", index=False)
        _write_partials(part, threshold, tmp_path / "partials" / f"{target}.tsv")
    common = ["--threshold-file", str(threshold)]
    compute(["--simulations", str(tmp_path / "sims"), *common, "--out", str(tmp_path / "a.tsv")])
    compute(["--partials", str(tmp_path / "partials"), *common, "--out", str(tmp_path / "b.tsv")])
    assert (tmp_path / "a.tsv").read_bytes() == (tmp_path / "b.tsv").read_bytes()

    with pytest.raises(SystemExit, match="partial counts were taken at"):
        compute(["--partials", str(tmp_path / "partials"), "--alpha", "0.05", "--out", "x.tsv"])
