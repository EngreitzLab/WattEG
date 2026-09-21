"""The per-gene model the simulation draws from: coefficients and theta.

sceptre's null model for a gene is a Poisson GLM of its counts on the cell
covariates, plus a separately estimated negative-binomial theta
(`perform_response_precomputation`). Both halves come out of one fit, and this
module reproduces that fit with pysceptre, which agrees with sceptre's own cache
to about ten significant digits on both -- theta to 2.8e-12 at the median, and
`exp(X . beta)` to 9.2e-10 across 134.5M gene x cell values.

**Both halves are used, and that is the point.** The R implementation took theta
from here and the *mean* from a DESeq2-style normalisation, so the simulated
counts had their dispersion from one model and their level from another. Drawing
from `exp(X . beta)` with this theta puts the simulation on the single scale the
discovery test itself works on. See `baseline.py`.

**Fitted over `cells_in_use`, and over every cell in them -- including a target's
perturbed cells.** That looks wrong for a baseline and is not: sceptre's own null
model is fitted the same way, on every cell, because the control group is the
complement. A gene with a real strong effect therefore carries a little of it in
its baseline, for every target. That was equally true of the normalised mean this
replaces, and the fit must not be changed to exclude perturbed cells, because the
test it is emulating does not.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pysceptre.pipeline.discovery import (
    _GENE_BATCH_WIDTH,
    _GENE_CHUNK_MEMORY_GB,
    fit_all_genes,
    summarize_gene_fits,
)


@dataclass(frozen=True)
class GeneModels:
    """One fit per gene, in the order `genes` gives.

    `diagnostics` names the genes whose GLM did not converge, whose theta MLE
    fell back to method of moments, and so on. pysceptre warns about the
    counts; the ids are what someone needs when one gene's power looks wrong,
    so they are carried rather than left in a warning string.
    """

    genes: list[str]
    fitted_coefs: np.ndarray  # (n_genes, p), aligned to the covariate matrix
    theta: np.ndarray  # (n_genes,) negative-binomial size
    diagnostics: dict[str, list[str]]

    @property
    def dispersion(self) -> np.ndarray:
        """`1/theta`, which is what the R pipeline called it."""
        return 1.0 / self.theta


def fit_gene_models(
    counts,
    gene_ids: list[str],
    covariate_matrix: np.ndarray,
    genes: list[str],
    *,
    n_jobs: int = 1,
) -> GeneModels:
    """Fit every gene's Poisson GLM and theta.

    `counts` is (n_genes, n_cells) over `cells_in_use`, `gene_ids` labels its
    rows, and `covariate_matrix` is (n_cells, p) over the same cells.

    The batching arguments are pysceptre's own pipeline defaults rather than
    this function's choices: `_GENE_BATCH_WIDTH = 1` makes each fit depend on
    that gene alone, which is what keeps a dispersion from shifting when the
    gene list changes.
    """
    index = {g: i for i, g in enumerate(gene_ids)}
    unknown = [g for g in genes if g not in index]
    if unknown:
        raise KeyError(
            f"{len(unknown)} gene(s) are not rows of the count matrix, including: {unknown[:5]}"
        )

    fits = fit_all_genes(
        counts,
        list(genes),
        covariate_matrix,
        chunk_memory_gb=_GENE_CHUNK_MEMORY_GB,
        gene_rows=[index[g] for g in genes],
        batch_width=_GENE_BATCH_WIDTH,
        n_jobs=n_jobs,
    )

    # A clamped theta is a fit that did not converge to anything usable, and a
    # dispersion of 100 or 0.001 would simulate a gene nothing like the real
    # one. R had no equivalent check because it read a cache that sceptre had
    # already produced; here the fit happens in front of us, so it is checked.
    clamped = [g for g in genes if fits[g].theta_clamped]
    if clamped:
        raise ValueError(
            f"{len(clamped)} gene(s) have a theta clamped to the estimator's bounds, including: "
            f"{clamped[:5]}. Their dispersion is not an estimate, and simulating from it would "
            "misstate the noise. Drop these genes from the discovery pairs."
        )
    nonfinite = [g for g in genes if not np.isfinite(fits[g].theta) or fits[g].theta <= 0]
    if nonfinite:
        raise ValueError(
            f"{len(nonfinite)} gene(s) have a non-finite or non-positive theta, including: "
            f"{nonfinite[:5]}."
        )

    return GeneModels(
        genes=list(genes),
        fitted_coefs=np.array([fits[g].fitted_coefs for g in genes]),
        theta=np.array([fits[g].theta for g in genes]),
        diagnostics=summarize_gene_fits(fits),
    )
