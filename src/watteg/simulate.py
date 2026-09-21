"""Drawing the counts.

One function, and the only place randomness touches expression. Everything it
needs has been decided elsewhere: the baseline by `baseline.py`, the per-cell
multiplier by `perturbation.py`, the dispersion by `gene_model.py`.

    mu[i, j] = baseline[i, j] * effect_size[i, j]
    count[i, j] ~ NegBinomial(mean = mu[i, j], size = theta[i])

`size` is theta, the same parameterisation sceptre's own model uses, so the
counts this draws are counts from the model the test assumes -- which is the
point of taking the baseline from that model too (see `baseline.py`).
"""

from __future__ import annotations

import numpy as np

# int16 holds counts to 32,767. The largest single count in a real screen
# measured here is 2,717, but a simulated draw has a tail, so the promotion is
# checked rather than assumed.
_COUNT_DTYPE = np.int16


def draw_counts(
    baseline: np.ndarray,
    effect_size: np.ndarray,
    theta: np.ndarray,
    rng: np.random.Generator,
    *,
    dtype: np.dtype | None = _COUNT_DTYPE,
) -> np.ndarray:
    """Simulated counts, `(n_genes, n_cells)`.

    `theta` is the NB size, one per gene, broadcast down the rows.

    numpy parameterises the negative binomial as `(n, p)` with mean
    `n(1-p)/p`, so `n = theta` and `p = theta / (theta + mu)` give mean `mu`
    and variance `mu + mu^2/theta` -- R's `rnbinom(mu=, size=)`.

    Returned as `int16` by default. These are counts; holding them as float64
    costs four times the memory for no information, and the simulation's whole
    shape depends on how many replicates fit in one process.
    """
    baseline = np.asarray(baseline, dtype=float)
    if baseline.shape != effect_size.shape:
        raise ValueError(f"baseline is {baseline.shape} but effect_size is {effect_size.shape}")
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (baseline.shape[0],):
        raise ValueError(f"theta is {theta.shape}, expected ({baseline.shape[0]},)")
    if np.any(theta <= 0) or not np.all(np.isfinite(theta)):
        raise ValueError("theta must be finite and positive")

    mu = baseline * effect_size
    size = theta[:, None]
    # p == 1 exactly where mu == 0, which numpy accepts and which draws 0
    # every time -- the right answer for a gene knocked all the way down.
    p = size / (size + mu)
    counts = rng.negative_binomial(size, p)

    if dtype is None:
        return counts
    info = np.iinfo(dtype)
    if counts.max(initial=0) > info.max:
        raise OverflowError(
            f"a simulated count exceeded {dtype.__name__}'s range ({info.max}); pass "
            "dtype=None to keep the draw at full width"
        )
    return counts.astype(dtype)
