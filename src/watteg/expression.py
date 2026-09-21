"""Per-gene and per-cell expression statistics: what the simulation draws from.

A port of `compute_expression_stats()` in `src/prepare_sim_input.R`, and the only
place in the Python pipeline that reads real expression values. Everything
downstream works from the three numbers this produces -- a per-cell size factor,
a per-gene size-factor-normalised mean, and a per-gene raw mean -- plus the
dispersion in `dispersion.py`.

**Computed over every cell, including the ones QC removed.** DESeq2 "poscounts"
size factors are a per-cell reduction against a per-gene geometric mean taken
across the whole matrix, so which cells are in the matrix changes the size
factors of the cells that stay: measured on day0, a median 0.46 % shift in size
factor and 0.36 % in normalised gene mean. The R implementation reads the
object's whole response matrix, so reproducing it means doing the same, which is
what `--all-cells` on pysceptre's export is for. Whether those cells *should*
count is a separate question -- see `docs/pysceptre-backend.md` section 5.2.

DESeq2 itself is not used, in either language: `DESeqDataSetFromMatrix()` coerces
to dense, which at 292 x 586,309 is 1.4 GB and pointless.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class ExpressionStats:
    """Aligned to the matrix that produced them: one entry per gene, per cell."""

    size_factors: np.ndarray  # (n_cells,)
    average_expression_all_cells: np.ndarray  # (n_genes,) raw mean
    normalized_mean: np.ndarray  # (n_genes,) what the simulation draws from
    density: float


def compute_expression_stats(counts: sparse.spmatrix) -> ExpressionStats:
    """`counts` is (n_genes, n_cells), sparse, non-negative.

    Never densified. The matrix is worked on in CSC because every statistic here
    is a per-cell reduction, and CSC is what makes a cell's nonzeros contiguous.
    """
    csc = counts.tocsc(copy=False)
    # Explicit zeros would take log(0) = -inf into the geometric mean and turn a
    # whole gene unusable. R calls drop0() here for the same reason.
    csc.eliminate_zeros()
    n_genes, n_cells = csc.shape
    if csc.data.size and csc.data.min() < 0:
        raise ValueError("count matrix contains negative values")

    # EXPLICIT float64, and not a formality. pysceptre's export stores counts as
    # `uint16` to keep the file small, and numpy's `log` of a 16-bit integer
    # array returns **float32**: every logarithm below would be computed in
    # single precision, silently. Measured against R on day0, that alone moved
    # the size factors by 2.5e-7 and the normalised gene means by 1.2e-9 --
    # small enough to look like floating-point noise and large enough not to be.
    data = csc.data.astype(np.float64, copy=False)
    gene_of_nonzero = csc.indices
    log_counts = np.log(data)

    # Geometric mean per gene, on the log scale, over the cells where the gene is
    # nonzero but divided by *every* cell -- that is what makes this "poscounts"
    # rather than an ordinary geometric mean.
    log_geomeans = np.bincount(gene_of_nonzero, weights=log_counts, minlength=n_genes) / n_cells
    row_totals = np.bincount(gene_of_nonzero, weights=data, minlength=n_genes)
    log_geomeans[row_totals == 0] = -np.inf
    usable_gene = np.isfinite(log_geomeans)

    # Each count against its gene's geometric mean. A gene with no usable mean
    # contributes nothing to any cell's median rather than contributing a zero,
    # which would drag every size factor towards 1.
    ratio = log_counts - log_geomeans[gene_of_nonzero]
    ratio[~usable_gene[gene_of_nonzero]] = np.nan

    size_factors = np.exp(_median_per_cell(ratio, csc.indptr, n_cells))

    bad = ~np.isfinite(size_factors) | (size_factors <= 0)
    if bad.any():
        raise ValueError(
            f"{bad.sum()} of {n_cells} cells have no usable size factor (no nonzero counts in "
            "any gene with a finite geometric mean). Filter these cells out before running the "
            "power analysis."
        )

    # Raw mean, reported in the output but never simulated from. Kept distinct
    # from the normalised mean below, which is what the simulation draws from.
    average_expression_all_cells = row_totals / n_cells

    cell_of_nonzero = _cell_of_nonzero(csc.indptr, data.size)
    normalized = data / size_factors[cell_of_nonzero]
    normalized_mean = np.bincount(gene_of_nonzero, weights=normalized, minlength=n_genes) / n_cells

    return ExpressionStats(
        size_factors=size_factors,
        average_expression_all_cells=average_expression_all_cells,
        normalized_mean=normalized_mean,
        density=data.size / (float(n_genes) * n_cells),
    )


def _cell_of_nonzero(indptr: np.ndarray, n_nonzero: int) -> np.ndarray:
    """Which cell each CSC nonzero belongs to."""
    return np.repeat(np.arange(indptr.size - 1), np.diff(indptr))


def _median_per_cell(ratio: np.ndarray, indptr: np.ndarray, n_cells: int) -> np.ndarray:
    """Median of each cell's ratios, NaN where a cell has none.

    One global sort rather than a loop over cells. The nonzeros are already
    grouped by cell -- that is what CSC means -- so sorting by `(cell, ratio)`
    leaves each cell's values contiguous *and* ordered, and every median is then
    an index lookup. At 95.7M nonzeros a Python loop over 586,309 cells is
    minutes; this is seconds.

    The even-length case averages the two middle values, which is what both R's
    `median.default` and `np.median` do, so the two agree element for element.
    """
    out = np.full(n_cells, np.nan)
    if ratio.size == 0:
        return out

    # NaN sorts last under lexsort, so a cell's usable values stay at the front
    # of its slice and the median is taken over the count of usable ones.
    cell_of_nonzero = _cell_of_nonzero(indptr, ratio.size)
    order = np.lexsort((ratio, cell_of_nonzero))
    sorted_ratio = ratio[order]

    usable_per_cell = np.bincount(cell_of_nonzero[np.isfinite(ratio)], minlength=n_cells).astype(
        np.int64
    )
    has_any = usable_per_cell > 0
    starts = indptr[:-1]
    lower = starts[has_any] + (usable_per_cell[has_any] - 1) // 2
    upper = starts[has_any] + usable_per_cell[has_any] // 2
    out[has_any] = (sorted_ratio[lower] + sorted_ratio[upper]) / 2.0

    return out
