"""Turning simulated replicates into a power estimate, and reading it back.

Power is the fraction of replicates in which the screen's own test would have
called the association:

    power = mean(p_value < threshold AND log_2_fold_change < 0)

The fold-change condition makes it one-sided -- a replicate counts only if the
simulated perturbation reduced expression, which is the claim being tested.

Because that is a binomial proportion over a finite number of replicates, the
point estimate travels with an interval, and the interval is Wilson's rather
than the normal approximation. That is not tidiness: `power_ci_low` is the
column thresholded to establish a *negative*, and the normal approximation
degenerates exactly where those live. Zero successes in 100 gives `[0, 0]`,
asserting the power is certainly zero; Wilson gives `[0, 0.037]`, which is what
100 replicates without a success actually supports.
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


def compute_power(
    simulations: pd.DataFrame, threshold: float, conf_level: float = 0.95
) -> pd.DataFrame:
    """One row per pair, from one row per (pair, replicate).

    Replicates with a missing p-value or fold change carry no information about
    whether the association would have been called, so they are dropped -- which
    means `n_reps` can fall below the replicate count asked for and is reported
    per pair rather than assumed.
    """
    missing = [c for c in REQUIRED_COLUMNS if c not in simulations.columns]
    if missing:
        raise KeyError(f"simulations are missing column(s): {missing}")
    if not 0 < conf_level < 1:
        raise ValueError(f"conf_level must lie in (0, 1), got {conf_level}")
    if not np.isfinite(threshold) or not 0 < threshold <= 1:
        raise ValueError(
            f"the p-value threshold ({threshold}) is not a usable probability. A non-finite "
            "value usually means the source discovery results held no significant pair; pass an "
            "explicit threshold instead."
        )

    usable = simulations.dropna(subset=["p_value", "log_2_fold_change"])
    if usable.empty:
        raise ValueError("no usable simulation rows remain")

    key = ["grna_target", "response_id", "rep"]
    duplicated = usable.duplicated(subset=key)
    if duplicated.any():
        raise ValueError(
            f"{int(duplicated.sum())} duplicated (grna_target, response_id, rep) row(s). "
            "Overlapping replicate offsets across chunks would double-count replicates."
        )

    if "effect_size" in usable.columns:
        sizes = usable["effect_size"].unique()
        if len(sizes) > 1:
            raise ValueError(
                f"the input mixes effect sizes ({', '.join(map(str, sizes))}). Compute power "
                "once per effect size."
            )

    called = (usable["p_value"] < threshold) & (usable["log_2_fold_change"] < 0)
    grouped = usable.assign(_called=called).groupby(
        ["grna_target", "response_id"], sort=False, as_index=False
    )
    agg = {
        "successes": ("_called", "sum"),
        "n_reps": ("_called", "size"),
        "mean_log_2_fold_change": ("log_2_fold_change", "mean"),
    }
    # Carried through when the simulation reported them, so the output is
    # self-describing.
    for column, name in (
        ("num_pert_cells", "mean_pert_cells"),
        ("average_expression_all_cells", "average_expression_all_cells"),
    ):
        if column in usable.columns:
            agg[name] = (column, "mean")
    if "effect_size" in usable.columns:
        agg["effect_size"] = ("effect_size", "first")

    power = grouped.agg(**agg)
    low, high = wilson_interval(power["successes"], power["n_reps"], conf_level)
    power.insert(2, "power", power["successes"] / power["n_reps"])
    power.insert(3, "power_ci_low", low)
    power.insert(4, "power_ci_high", high)
    power = power.drop(columns="successes")

    return power.sort_values(
        ["power", "grna_target", "response_id"], ascending=[False, True, True]
    ).reset_index(drop=True)


def min_detectable_effect_size(
    power_by_effect: np.ndarray, effect_sizes: list[float], power_threshold: float
) -> np.ndarray:
    """Smallest tested effect size from which a pair is detectable.

    `power_by_effect` is (n_pairs, n_effect_sizes), its columns ordered to match
    `effect_sizes` ascending. NaN means the pair was not tested at that effect
    size, which is unknown rather than a failure.

    **A pair qualifies at `e` only if it clears the threshold at `e` and at
    every larger effect size tested.** Taking the first effect size that clears,
    in isolation, lets Monte Carlo noise win: at 100 replicates a pair whose
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
