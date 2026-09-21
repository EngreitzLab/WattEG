"""Per-gene dispersion: the noise the simulation exists to reproduce.

The simulation draws counts from a negative binomial with `size = 1/dispersion`,
so this number sets how noisy a simulated gene is, and therefore how much of the
effect a test can see. It has to come from the **real** counts: a dispersion
estimated from simulated data would make the simulation agree with itself.

In R it was read from `sceptre_object@response_precomputations`, a cache sceptre
fills during the real discovery analysis. There is no such cache here and there
does not need to be -- pysceptre fits the same model, so the Python pipeline
estimates theta directly and depends on no slot of anyone's object.

**Fitted over `cells_in_use`, unlike everything in `expression.py`.** sceptre's
precomputation runs on the cells that passed QC, against the covariate matrix,
so matching it means doing the same. That the two halves of `prepare_sim_input`
use different cell sets is not an oversight in either language: size factors are
a property of the library preparation, which every cell took part in, while a
model fit is a property of the analysis, which only the QC-passing cells enter.
"""

from __future__ import annotations

import numpy as np
from pysceptre.pipeline.discovery import (
    _GENE_BATCH_WIDTH,
    _GENE_CHUNK_MEMORY_GB,
    fit_all_genes,
)


def fit_dispersions(
    counts,
    gene_ids: list[str],
    covariate_matrix: np.ndarray,
    genes: list[str],
    *,
    n_jobs: int = 1,
) -> dict[str, float]:
    """Dispersion (`1/theta`) per gene, in the order `genes` gives.

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

    return {g: 1.0 / fits[g].theta for g in genes}
