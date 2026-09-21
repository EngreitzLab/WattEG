"""Which cells a target perturbs, which guide each one carries, and how hard.

A port of `lib/pert_input.R` and the guide-level half of `lib/simulate.R`. The
modelling content is one idea: **a target's guides do not all work equally
well**, so a simulated knockdown is not a single multiplier applied to every
perturbed cell. Each guide draws its own effect size around the target's, and
each cell inherits the effect of whichever guide it happens to carry.

That is why the pipeline needs per-gRNA assignments at all, and why the
per-target union sceptre keeps is not enough (see
`docs/pysceptre-backend.md` section 5.1).

**Control cells carry guides too**, and they matter. A cell not perturbed by
this target is usually perturbed by something else -- another element, or a
non-targeting guide -- and that guide draws an effect size around 1 with the
same guide-to-guide spread. Dropping them would leave every control cell at
exactly 1, which keeps the control arm's mean and removes its variance.

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


def effect_size_matrix(
    assignment: GuideAssignment,
    gene_effect_sizes: np.ndarray,
    guide_sd: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """A (genes x cells) multiplier, with each guide's own effect size.

    `gene_effect_sizes` is the *relative expression* the target aims for -- a
    15% knockdown is 0.85 -- one per gene. Each of the target's guides draws
    around it, each other guide draws around 1, and negatives clamp to 0
    because a guide cannot make expression negative.

    Clamping biases the mean upward, so both arms are re-centred afterwards:
    that is what makes "effect size 0.15" mean a 15% knockdown on average
    rather than something slightly weaker. R did this in a separate
    `center_effect_size_matrix()`; it is one step here because the two are not
    separately meaningful.
    """
    gene_effect_sizes = np.asarray(gene_effect_sizes, dtype=float)
    n_genes = gene_effect_sizes.size
    n_target, n_other = assignment.n_target_guides, assignment.n_other_guides

    # Row 0 is the no-effect row, for cells carrying no guide at all.
    table = np.empty((1 + n_target + n_other, n_genes))
    table[0] = 1.0
    table[1 : 1 + n_target] = rng.normal(gene_effect_sizes, guide_sd, size=(n_target, n_genes))
    table[1 + n_target :] = rng.normal(1.0, guide_sd, size=(n_other, n_genes))
    np.clip(table, 0.0, None, out=table)

    matrix = table[assignment.status].T  # (genes, cells)

    # Re-centre each arm on what it is supposed to average to. Cells carrying
    # no guide sit at exactly 1 and move with their arm, as in R.
    for mask, target_mean in (
        (assignment.is_perturbed, gene_effect_sizes),
        (~assignment.is_perturbed, 1.0),
    ):
        if not mask.any():
            continue
        block = matrix[:, mask]
        shift = target_mean - block.mean(axis=1)
        matrix[:, mask] = block + shift[:, None]
    np.clip(matrix, 0.0, None, out=matrix)
    return matrix
