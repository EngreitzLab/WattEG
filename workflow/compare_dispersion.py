#!/usr/bin/env python3
"""Phase 2's other half: Python's fitted theta against sceptre's cached one.

    workflow/compare_dispersion.py <export> <r_reference_dir>

This is **not** a pass/fail comparison and deliberately has no tolerance. The
other quantities in `prepare_sim_input` are arithmetic over the same numbers, so
they can be required to agree to floating point. A dispersion is the output of
an iterative maximum-likelihood fit: sceptre's cache came from its own Poisson
IRLS and theta estimator, this comes from pysceptre's, and the two agreeing to
three or four digits is the expected result rather than a disappointing one.

What the output is for is judging whether the difference could move power.
The simulation draws from `NB(mu, size = 1/dispersion)`, so a relative shift in
dispersion is roughly a relative shift in simulated variance -- read the
distribution below against the effect sizes being tested, which are 5-50 %.

Reported: the distribution over all genes, which genes sit furthest apart, and
the theta estimator each side used, because a gene that fell back to a different
method is the first place to look when one disagrees.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from watteg.dispersion import fit_dispersions


def read_floats(path: Path) -> np.ndarray:
    return np.array([float(x) for x in path.read_text().split()])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path, help="a pysceptre export")
    parser.add_argument("reference", type=Path, help="output of workflow/dump_r_sim_input.R")
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pysceptre" / "scripts"))
    from sceptre_io import load_export

    # The default load, NOT all_cells: sceptre fits its precomputation over the
    # QC-passing cells, so matching it means fitting over those and no others.
    export = load_export(args.export)
    genes = (args.reference / "genes.txt").read_text().split()
    print(
        f"fitting {len(genes)} genes over {export.covariate_matrix.shape[0]:,} cells_in_use, "
        f"{export.covariate_matrix.shape[1]} covariates ..."
    )

    ours = fit_dispersions(
        export.response_matrix,
        export.gene_ids,
        export.covariate_matrix,
        genes,
        n_jobs=args.n_jobs,
    )
    mine = np.array([ours[g] for g in genes])
    theirs = read_floats(args.reference / "row_dispersion.txt")

    rel = np.abs(mine - theirs) / theirs
    print(f"\ndispersion, {len(genes)} genes (R's cache vs pysceptre's fit)")
    for q in (50, 90, 99, 100):
        print(f"  {q:>3}th percentile of |relative difference|: {np.percentile(rel, q):.3e}")
    for tol in (1e-12, 1e-9, 1e-6, 1e-3):
        print(f"  within {tol:.0e}: {(rel < tol).sum()}/{rel.size}")
    print(f"  bit-identical: {(mine == theirs).sum()}/{rel.size}")

    worst = np.argsort(rel)[-5:][::-1]
    print("\n  furthest apart:")
    for i in worst:
        print(
            f"    {genes[i]:<16} ours {mine[i]:.6g}  R {theirs[i]:.6g}  "
            f"(rel {rel[i]:.3e}; theta {1 / mine[i]:.10g} vs {1 / theirs[i]:.10g})"
        )

    # A dispersion difference matters only through the variance of the counts
    # drawn from it, so say it in those terms rather than leaving the reader to.
    print(
        f"\n  A gene at the median difference simulates with {np.median(rel):.2e} relative "
        "difference in negative-binomial variance against R's."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
