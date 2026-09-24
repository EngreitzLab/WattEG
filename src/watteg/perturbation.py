"""Which cells a target perturbs, which guide each one carries, and how hard.

A port of `lib/pert_input.R` and the guide-level half of `lib/simulate.R`. The
modelling content is one idea: **a target's guides do not all work equally
well**, so a simulated knockdown is not a single multiplier applied to every
perturbed cell. Each guide draws its own effect size around the target's, and
each cell inherits the effect of whichever guide it happens to carry.

That is why the pipeline needs per-gRNA assignments at all, and why the
per-target union sceptre keeps is not enough (see
`docs/pysceptre-backend.md` section 5.1).

**Other guides have no effect on the tested genes.** A control cell usually
carries a guide for some other element, or a non-targeting one, and that guide
does not move this gene: its effect is exactly 1. The dispersion the counts are
drawn with was fitted to real cells that already carry their real guides, so
it already contains whatever those guides do. Until 2026-09-24 every other
guide drew an extra N(1, guide_sd), as in the original DC-TAP code. That
counted the noise twice (theta 146 refit to 42 on a gene's own simulated null
data) and understated power for highly expressed, low-dispersion genes.

**The perturbed cells' mean is pinned.** Simulated power is power at a FIXED
element effect (decided 2026-09-24; see docs/methods.md): in every replicate
the realised mean effect across the perturbed cells equals the requested one.
The target's guides still differ from each other; only their mean is fixed.

**Built directly in cell order.** R assembled the guide assignment as
perturbed-cells-then-control-cells and then permuted it back, a legacy of
matching cell barcodes. Nothing here needs that, and the streams differ between
the languages regardless.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

NO_GUIDE = 0


@dataclass(frozen=True)
class GuideAssignment:
    """One guide index per cell, and what the indices mean.

    `status[j]` is 0 when cell `j` carries no guide, `1..n_target_guides` for
    one of this target's guides, and above that for a guide belonging to
    something else. The offset is what lets one effect-size table cover both
    arms: rows are laid out as [no effect, this target's guides, every other
    guide], so a cell's effect size is a single lookup.
    """

    status: np.ndarray  # (n_cells,) int
    is_perturbed: np.ndarray  # (n_cells,) bool
    n_target_guides: int
    n_other_guides: int


def target_cells(cre_perts: sparse.spmatrix, target_ids: list[str], target: str) -> np.ndarray:
    """Boolean mask of the cells this target perturbs."""
    try:
        row = target_ids.index(target)
    except ValueError:
        raise KeyError(f"target {target!r} is not a row of cre_perts") from None
    return np.asarray(cre_perts[row].todense()).ravel() > 0


def _choose_one_per_cell(
    csc: sparse.csc_matrix,
    eligible_row: np.ndarray,
    cells: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """For each cell in `cells`, one uniformly chosen eligible guide, or -1.

    `csc` is (guides x cells), so a cell's guides are one contiguous slice --
    that is what CSC means, and it is why this needs no sort. Vectorised over
    cells rather than looped: R drew one cell at a time, which at 567,690 cells
    is the difference between seconds and minutes.
    """
    keep = eligible_row[csc.indices]
    cell_of_nonzero = np.repeat(np.arange(csc.shape[1]), np.diff(csc.indptr))

    qualifying_cell = cell_of_nonzero[keep]
    qualifying_guide = csc.indices[keep]
    per_cell = np.bincount(qualifying_cell, minlength=csc.shape[1])
    # Where each cell's qualifying entries begin in the compacted arrays.
    starts = np.concatenate(([0], np.cumsum(per_cell)[:-1]))

    counts = per_cell[cells]
    out = np.full(cells.size, -1, dtype=np.int64)
    has_any = counts > 0
    if not has_any.any():
        return out
    # One uniform per cell, floored into that cell's own count. Drawing for
    # every cell and discarding the unused ones keeps the number of draws a
    # function of the cell set alone, not of how many happen to carry a guide.
    picks = np.floor(rng.random(cells.size) * counts).astype(np.int64)
    picks = np.minimum(picks, np.maximum(counts - 1, 0))  # guard the rng.random() == 1 edge
    out[has_any] = qualifying_guide[starts[cells[has_any]] + picks[has_any]]
    return out


def guide_assignment(
    grna_csc: sparse.csc_matrix,
    grna_ids: list[str],
    target_guides: list[str],
    is_perturbed: np.ndarray,
    rng: np.random.Generator,
) -> GuideAssignment:
    """Pick one guide per cell: this target's for perturbed cells, any other's
    for control cells.

    `target_guides` comes from the screen's design table, which maps guides to
    targets **many-to-many** -- a guide inside two overlapping candidate
    elements belongs to both. It must never be derived from a per-unit
    annotation that can hold only one target per guide.

    Drawn once per target and reused across replicates: which guide a cell
    carries is a fact about the screen, not about a draw.
    """
    index = {g: i for i, g in enumerate(grna_ids)}
    rows = np.array([index[g] for g in target_guides if g in index], dtype=np.int64)
    if rows.size == 0:
        raise KeyError(
            f"none of the {len(target_guides)} guide(s) given for this target appear in "
            "grna_perts, so no perturbed cell can be assigned a guide"
        )

    n_guides = len(grna_ids)
    is_target_guide = np.zeros(n_guides, dtype=bool)
    is_target_guide[rows] = True

    status = np.zeros(grna_csc.shape[1], dtype=np.int64)
    perturbed = np.flatnonzero(is_perturbed)
    control = np.flatnonzero(~is_perturbed)

    # Target guides are numbered by their position in `target_guides`, not by
    # their row in the matrix, so the effect-size table's first block lines up
    # with the guides whose effect sizes it holds.
    rank_of_row = np.full(n_guides, -1, dtype=np.int64)
    rank_of_row[rows] = np.arange(rows.size)

    chosen = _choose_one_per_cell(grna_csc, is_target_guide, perturbed, rng)
    found = chosen >= 0
    status[perturbed[found]] = rank_of_row[chosen[found]] + 1

    other_rows = np.flatnonzero(~is_target_guide)
    rank_of_other = np.full(n_guides, -1, dtype=np.int64)
    rank_of_other[other_rows] = np.arange(other_rows.size)
    chosen = _choose_one_per_cell(grna_csc, ~is_target_guide, control, rng)
    found = chosen >= 0
    status[control[found]] = rows.size + rank_of_other[chosen[found]] + 1

    return GuideAssignment(
        status=status,
        is_perturbed=is_perturbed,
        n_target_guides=int(rows.size),
        n_other_guides=int(other_rows.size),
    )


def _pin_to_mean(block: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Row by row, the values max(v + c, 0) whose mean is exactly `target`.

    f(c) = mean(max(v + c, 0)) is continuous, non-decreasing and piecewise
    linear, with a kink where each value hits zero, so for any target > 0 it has
    a root. Sorting a row in decreasing order: if exactly the top j values stay
    positive, f(c) = (sum of those j + j * c) / n, so c = (n * target - sum) / j.
    The right j is the largest one whose kink, f(-s_j), is still at or below the
    target. One sort per row, no iteration.

    The closed form carries the rounding of a long cumulative sum, and
    perturbed-cell effects are heavily tied (each cell takes one of a few guide
    values), so that error adds up rather than cancelling: at tens of thousands
    of cells it reached ~1e-12. One correction step on the linear piece the
    solution lies on removes it. Same algorithm as R's `pin_to_mean()`.
    """
    n = block.shape[1]
    s = -np.sort(-block, axis=1)  # each row in decreasing order
    cs = np.cumsum(s, axis=1)
    j = np.arange(1, n + 1)
    f_at_kink = (cs - j * s) / n  # the mean when the j-th largest value just reaches 0
    ok = f_at_kink <= target[:, None]
    # Largest j (1-based) with ok. f_at_kink[:, 0] is 0, so every row has one.
    last = n - np.argmax(ok[:, ::-1], axis=1)
    shift = (n * target - cs[np.arange(block.shape[0]), last - 1]) / last
    out = np.clip(block + shift[:, None], 0.0, None)
    positive = out > 0
    k = positive.sum(axis=1)
    delta = np.where(k > 0, (target - out.mean(axis=1)) * n / np.maximum(k, 1), 0.0)
    return np.where(positive, np.clip(out + delta[:, None], 0.0, None), out)


def effect_size_matrix(
    assignment: GuideAssignment,
    gene_effect_sizes: np.ndarray,
    guide_sd: float,
    rng: np.random.Generator,
    *,
    tol: float = 1e-12,
) -> np.ndarray:
    """A (genes x cells) multiplier, with each of the target's guides at its own effect.

    `gene_effect_sizes` is the *relative expression* the target aims for -- a
    15% knockdown is 0.85 -- one per gene. Each of the target's guides draws
    around it, independently per gene, and negatives clamp to 0 because a guide
    cannot make expression negative. Every other guide, and every cell carrying
    no guide, is exactly 1.

    Then each gene's mean over the perturbed cells is **pinned** to what was
    asked for. That is what makes simulated power the power at a fixed element
    effect (see the module docstring). It is not a correction for clamping,
    which is what this docstring used to claim: at es 0.15 a guide clamps with
    probability ~3e-11. The result is max(v + c, 0), where v are the perturbed
    cells' already-clamped effects and c is the one constant that puts the mean
    exactly on the target (`_pin_to_mean`). Where nothing would clamp, which is
    every realistic case up to es 0.5, that is the plain shift, to a few ulps. An earlier version
    shifted, clamped and repeated; once most cells clamp that converges slowly,
    and at es >= ~0.99 it ran out of iterations and raised on a pin that exists.
    R does the same in `center_effect_size_matrix()`.
    """
    gene_effect_sizes = np.asarray(gene_effect_sizes, dtype=float)
    n_genes = gene_effect_sizes.size
    n_target, n_other = assignment.n_target_guides, assignment.n_other_guides

    # Row 0 is the no-effect row, for cells carrying no guide at all; the rows
    # after the target's block belong to other guides. Both stay at exactly 1.
    table = np.ones((1 + n_target + n_other, n_genes))
    table[1 : 1 + n_target] = rng.normal(gene_effect_sizes, guide_sd, size=(n_target, n_genes))
    np.clip(table, 0.0, None, out=table)

    matrix = table[assignment.status].T  # (genes, cells)

    control = ~assignment.is_perturbed
    if (matrix[:, control] != 1.0).any():
        raise RuntimeError(
            "control cells must carry an effect of exactly 1; a control cell's guide "
            "status points into the target's block"
        )

    perturbed = assignment.is_perturbed
    if perturbed.any():
        block = _pin_to_mean(matrix[:, perturbed], gene_effect_sizes)
        gap = gene_effect_sizes - block.mean(axis=1)
        # The check is for real failures (NaN, a logic error), not rounding, so its
        # tolerance grows with the number of perturbed cells as summation error does.
        allowed = max(tol, 4 * block.shape[1] * np.finfo(float).eps)
        if not np.all(np.isfinite(gap)) or np.any(np.abs(gap) >= allowed):
            raise ValueError(
                "could not pin the realised mean effect to the requested one (largest miss "
                f"{np.nanmax(np.abs(gap)) if np.isfinite(gap).any() else float('nan'):.3g})"
            )
        matrix[:, perturbed] = block
    return matrix
