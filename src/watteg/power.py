"""Turning simulations into a power estimate, and reading it back.

Power is the fraction of simulations in which the screen's own test would have
called the association:

    power = mean(p_value < threshold AND log_2_fold_change < 0)

The fold-change condition makes it one-sided -- a simulation counts only if the
simulated perturbation reduced expression, which is the claim being tested.

Because that is a binomial proportion over a finite number of simulations, the
point estimate travels with an interval, and the interval is Wilson's rather
than the normal approximation. That is not tidiness: `power_ci_low` is the
column thresholded to establish a *negative*, and the normal approximation
degenerates exactly where those live. Zero successes in 100 gives `[0, 0]`,
asserting the power is certainly zero; Wilson gives `[0, 0.037]`, which is what
100 simulations without a success actually supports.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

REQUIRED_COLUMNS = ("grna_target", "response_id", "p_value", "log_2_fold_change", "rep")


def wilson_interval(
    successes: np.ndarray, n: np.ndarray, conf_level: float = 0.95
) -> tuple[np.ndarray, np.ndarray]:
    """Wilson score interval for a binomial proportion, clipped to [0, 1]."""
    z = norm.ppf(1 - (1 - conf_level) / 2)
    successes = np.asarray(successes, dtype=float)
    n = np.asarray(n, dtype=float)
    phat = successes / n
    denom = 1 + z**2 / n
    center = (phat + z**2 / (2 * n)) / denom
    halfwidth = (z / denom) * np.sqrt(phat * (1 - phat) / n + z**2 / (4 * n**2))
    return np.maximum(0.0, center - halfwidth), np.minimum(1.0, center + halfwidth)


def effect_label(effect_size: float) -> str:
    """Column-name fragment for an effect size: 0.15 -> "15", 0.125 -> "12_5".

    Trailing zeros are trimmed **only when there is a fractional part**.
    Trimming unconditionally turns 0.2 into "2" and 0.5 into "5", so a 20%
    knockdown gets a column called `power_at_effect_size_2`. That shipped once.

    The decimal point becomes an underscore so every column name stays
    snake_case and free of dots, which several readers treat as separators.
    """
    formatted = f"{effect_size * 100:g}"
    if "." in formatted:
        formatted = formatted.rstrip("0").rstrip(".")
    return formatted.replace(".", "_")


# Carried through when the simulation reported them, so the output is self-describing: a per-pair
# mean in the power table, a per-pair sum in the counts.
_MEANS = (
    ("num_pert_cells", "mean_pert_cells"),
    ("average_expression_all_cells", "average_expression_all_cells"),
)
PAIR = ["grna_target", "response_id"]


def _check_threshold(threshold: float) -> None:
    if not np.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError(
            f"the p-value threshold ({threshold}) is not a usable probability. A non-finite "
            "value usually means the source discovery results held no significant pair; pass an "
            "explicit threshold instead."
        )


def _check_one_value(frame: pd.DataFrame, column: str, what: str, advice: str) -> None:
    if column in frame.columns and frame[column].nunique() > 1:
        values = ", ".join(map(str, frame[column].unique()))
        raise ValueError(f"the input mixes {what} ({values}). {advice}")


def power_counts(simulations: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Per pair, the sums power is made of, from one row per (pair, simulation).

    `successes` (simulations called: p below the threshold AND expression down), `n_reps` (the
    usable simulations), `sum_log_2_fold_change` and the sums behind the carried means, plus the
    simulation range and row count they came from. Counts from disjoint sets of simulations add
    (`merge_power_counts`), and `power_from_counts` turns them into the power table, so a task that
    holds all simulations of its pairs can write these instead of its rows. The means are sums
    divided by counts, which is exactly what pandas' grouped mean computes, so the table does not
    depend on which route it took.

    Simulations with a missing p-value or fold change carry no information about whether the
    association would have been called, so they are dropped -- which means `n_reps` can fall
    below the simulation count asked for and is reported per pair rather than assumed.
    """
    missing = [c for c in REQUIRED_COLUMNS if c not in simulations.columns]
    if missing:
        raise KeyError(f"simulations are missing column(s): {missing}")
    _check_threshold(threshold)

    usable = simulations.dropna(subset=["p_value", "log_2_fold_change"])
    if usable.empty:
        raise ValueError("no usable simulation rows remain")

    key = ["grna_target", "response_id", "rep"]
    duplicated = usable.duplicated(subset=key)
    if duplicated.any():
        raise ValueError(
            f"{int(duplicated.sum())} duplicated (grna_target, response_id, rep) row(s). "
            "Overlapping simulation offsets across chunks would double-count simulations."
        )
    _check_one_value(usable, "effect_size", "effect sizes", "Compute power once per effect size.")
    _check_one_value(
        usable,
        "estimand",
        "estimands",
        "They answer different questions; compute power once per estimand.",
    )

    called = (usable["p_value"] < threshold) & (usable["log_2_fold_change"] < 0)
    agg = {
        "successes": ("_called", "sum"),
        "n_reps": ("_called", "size"),
        "sum_log_2_fold_change": ("log_2_fold_change", "sum"),
    }
    for column, _ in _MEANS:
        if column in usable.columns:
            agg[f"sum_{column}"] = (column, "sum")
    for column in ("effect_size", "estimand"):
        if column in usable.columns:
            agg[column] = (column, "first")
    counts = usable.assign(_called=called).groupby(PAIR, sort=False, as_index=False).agg(**agg)
    # Over every row, usable or not: what a task produced, for the checks downstream. A pair with
    # no usable simulation is kept, at n_reps 0, so a task's counts still list every pair it ran;
    # `power_from_counts` then leaves it out of the table, as a table built from rows does.
    produced = simulations.groupby(PAIR, sort=False, as_index=False).agg(
        rep_first=("rep", "min"), rep_last=("rep", "max"), n_simulations=("rep", "size")
    )
    counts = produced.merge(counts, on=PAIR, how="left", validate="one_to_one")
    for column in ("successes", "n_reps"):
        counts[column] = counts[column].fillna(0).astype(np.int64)
    counts["threshold"] = float(threshold)
    return counts


def merge_power_counts(parts: list[pd.DataFrame]) -> pd.DataFrame:
    """Add up `power_counts` from disjoint sets of simulations, pair by pair.

    Refuses counts made at different thresholds, effect sizes or estimands, and a pair whose
    simulation ranges overlap, which would count a simulation twice. A pair that appears once is
    passed through unchanged, so when every task held all simulations of its pairs the merged
    counts are the ones computed from all the rows.
    """
    counts = pd.concat(parts, ignore_index=True)
    if counts.empty:
        raise ValueError("no power counts to merge")
    _check_one_value(counts, "threshold", "thresholds", "Merge counts from one threshold only.")
    _check_one_value(counts, "effect_size", "effect sizes", "Merge one effect size at a time.")
    _check_one_value(counts, "estimand", "estimands", "Merge one estimand at a time.")

    ordered = counts.sort_values([*PAIR, "rep_first"], kind="stable")
    same_pair = ordered[PAIR].eq(ordered[PAIR].shift()).all(axis=1)
    overlap = same_pair & (ordered["rep_first"] <= ordered["rep_last"].shift())
    if overlap.any():
        first = ordered.loc[overlap].iloc[0]
        raise ValueError(
            f"{int(overlap.sum())} pair(s) have counts from overlapping simulation ranges, e.g. "
            f"{first['grna_target']} / {first['response_id']} from simulation "
            f"{first['rep_first']}; each simulation would be counted twice"
        )
    if not counts.duplicated(subset=PAIR).any():
        return counts

    additive = ["successes", "n_reps", "n_simulations"]
    additive += [c for c in counts.columns if c.startswith("sum_")]
    agg = {c: (c, "sum") for c in additive}
    agg |= {"rep_first": ("rep_first", "min"), "rep_last": ("rep_last", "max")}
    for column in ("effect_size", "estimand", "threshold"):
        if column in counts.columns:
            agg[column] = (column, "first")
    return counts.groupby(PAIR, sort=False, as_index=False).agg(**agg)


def power_from_counts(counts: pd.DataFrame, conf_level: float = 0.95) -> pd.DataFrame:
    """The power table, one row per pair, from `power_counts` or `merge_power_counts`."""
    if not 0 < conf_level < 1:
        raise ValueError(f"conf_level must lie in (0, 1), got {conf_level}")
    counts = counts[counts["n_reps"] > 0]
    power = counts[PAIR].copy()
    n = counts["n_reps"]
    power["n_reps"] = n
    power["mean_log_2_fold_change"] = counts["sum_log_2_fold_change"] / n
    for column, name in _MEANS:
        if f"sum_{column}" in counts.columns:
            power[name] = counts[f"sum_{column}"] / n
    for column in ("effect_size", "estimand"):
        if column in counts.columns:
            power[column] = counts[column]
    low, high = wilson_interval(counts["successes"], n, conf_level)
    power.insert(2, "power", counts["successes"] / n)
    power.insert(3, "power_ci_low", low)
    power.insert(4, "power_ci_high", high)

    return power.sort_values(
        ["power", "grna_target", "response_id"], ascending=[False, True, True]
    ).reset_index(drop=True)


def compute_power(
    simulations: pd.DataFrame, threshold: float, conf_level: float = 0.95
) -> pd.DataFrame:
    """One row per pair, from one row per (pair, simulation): `power_counts`, then
    `power_from_counts`."""
    if not 0 < conf_level < 1:
        raise ValueError(f"conf_level must lie in (0, 1), got {conf_level}")
    return power_from_counts(power_counts(simulations, threshold), conf_level)


def min_detectable_effect_size(
    power_by_effect: np.ndarray, effect_sizes: list[float], power_threshold: float
) -> np.ndarray:
    """Smallest tested effect size from which a pair is detectable.

    `power_by_effect` is (n_pairs, n_effect_sizes), its columns ordered to match
    `effect_sizes` ascending. NaN means the pair was not tested at that effect
    size, which is unknown rather than a failure.

    **A pair qualifies at `e` only if it clears the threshold at `e` and at
    every larger effect size tested.** Taking the first effect size that clears,
    in isolation, lets Monte Carlo noise win: at 100 simulations a pair whose
    true power is 0.75 clears 0.8 about a third of the time, so across six
    effect sizes a spurious early clear is likely -- and the error runs one way,
    always reporting the pair as more detectable than it is. "Detectable from
    `e` upwards" is also the claim the column gets used to make.

    NaN in the result means no tested effect size qualified, which is a
    statement about the effect sizes that were run rather than proof the pair is
    undetectable.
    """
    values = np.asarray(power_by_effect, dtype=float)
    out = np.full(values.shape[0], np.nan)
    for row in range(values.shape[0]):
        best = np.nan
        # Walk down from the largest: the first known failure disqualifies
        # everything weaker than it.
        for i in range(values.shape[1] - 1, -1, -1):
            value = values[row, i]
            if np.isnan(value):
                continue
            if value < power_threshold:
                break
            best = effect_sizes[i]
        out[row] = best
    return out
