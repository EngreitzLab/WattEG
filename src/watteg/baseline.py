"""A gene's expected counts per cell, before any perturbation.

This is the number the simulation multiplies the effect size into, so it decides
what "unperturbed" means and therefore what power is measured against. There are
two ways to produce it and they are not equivalent.

**`fitted` (the default): `exp(X_j . beta_i)`.** The expected count sceptre's own
null model gives that cell for that gene. Simulation and test then live on one
scale by construction rather than by coincidence.

**`size_factor`: `mean_i x sf_j`.** What the R implementation did until
2026-09-21 -- a size-factor-normalised gene mean scaled by the cell's DESeq2
"poscounts" factor. Kept because the sweeps already run were produced with it,
and a comparison against one of those is only interpretable like for like. It is
not what R does now: `r-implementation` defaults to the fitted baseline too. It is a validation fixture,
not a modelling choice anyone should make afresh.

## Why the default changed

Three reasons, in increasing order of how much they matter.

**It mixes two models.** The R scheme takes its dispersion from sceptre's
negative-binomial fit and its level from a DESeq2 normalisation, so a simulated
gene's noise and its expression come from different statistical models of the
same data.

**It gets the level wrong, and the error survives the size factor.**
`row_data$mean` sits 16 % below the mean sceptre's model implies; multiplying by
the cell's size factor recovers most of that, leaving the simulated genes about
4 % low on day0. The residual is a dropped covariance term -- `mean_i` is a mean
of ratios, and `E[x.sf] = E[x].E[sf] + Cov(x, sf)`.

**It gets the shape wrong, which the level hides.** sceptre's mean varies with
every covariate -- library size, detected genes, batch, replicate -- while
`mean_i x sf_j` varies with a single scalar per cell. Measured against the real
counts on day0, over 60 genes and 567,690 cells:

| | `mean_i x sf_j` | `exp(X . beta)` |
|---|---|---|
| zero fraction, mean absolute error | 0.0053 | **0.0007** |
| predicted variance / observed | 0.865 | **0.995** |

The current scheme understates the count variance by 13.5 %. Note that pulls the
**opposite** way from the level error -- less variance inflates power, less
expression deflates it -- so the effect of the change on simulated power is not
obviously signed and has not been measured. Stage B under both modes is what
measures it.
"""

from __future__ import annotations

import numpy as np

MODELS = ("fitted", "size_factor")


def fitted_baseline(fitted_coefs: np.ndarray, covariate_matrix: np.ndarray) -> np.ndarray:
    """`exp(X . beta)` for one gene or many.

    `fitted_coefs` is `(p,)` for a single gene or `(n_genes, p)`; the result is
    `(n_cells,)` or `(n_genes, n_cells)` to match. One matrix product, so the
    coefficients are what a `sim_input` stores rather than the baseline itself:
    11 numbers per gene against one per cell per gene.
    """
    coefs = np.asarray(fitted_coefs, dtype=float)
    single = coefs.ndim == 1
    eta = covariate_matrix @ (coefs if single else coefs.T)
    return np.exp(eta if single else eta.T)


def size_factor_baseline(mean: np.ndarray, size_factors: np.ndarray) -> np.ndarray:
    """`mean_i x sf_j`, the R implementation's baseline.

    `mean` is `(n_genes,)` or a scalar; the result is `(n_genes, n_cells)` or
    `(n_cells,)`.
    """
    mean = np.asarray(mean, dtype=float)
    if mean.ndim == 0:
        return float(mean) * np.asarray(size_factors, dtype=float)
    return np.outer(mean, size_factors)


def baseline_expression(
    model: str,
    *,
    fitted_coefs: np.ndarray | None = None,
    covariate_matrix: np.ndarray | None = None,
    mean: np.ndarray | None = None,
    size_factors: np.ndarray | None = None,
) -> np.ndarray:
    """Dispatch on the model name, refusing the arguments the other one needs.

    Keyword-only and explicit about what is missing: passing `mean` to the
    fitted model, or coefficients to the size-factor model, is a mistake that
    would otherwise be silent in a pipeline where both are available.
    """
    if model not in MODELS:
        raise ValueError(f"model must be one of {list(MODELS)}, got {model!r}")
    if model == "fitted":
        if fitted_coefs is None or covariate_matrix is None:
            raise ValueError("the 'fitted' baseline needs fitted_coefs and covariate_matrix")
        return fitted_baseline(fitted_coefs, covariate_matrix)
    if mean is None or size_factors is None:
        raise ValueError("the 'size_factor' baseline needs mean and size_factors")
    return size_factor_baseline(mean, size_factors)
