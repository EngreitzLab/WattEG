#!/usr/bin/env python3
"""Phase 2's third gate: the per-gene model, coefficients as well as theta.

    workflow/compare_gene_model.py <export> <r_reference_dir>

`compare_dispersion.py` checked theta. Once the simulation draws its baseline
from `exp(X . beta)` rather than from a normalised mean, **the coefficients are
load-bearing too** -- they set the expected count of every cell -- so they need
a gate of their own rather than a spot check.

**Reported per coefficient column, not pooled.** The eleven coefficients live on
different scales and mean different things: an intercept, four log-covariate
slopes, and six factor dummies. A pooled maximum would let the batch dummies
disagree while the continuous terms carried the summary, and a factor-contrast
ordering mismatch is exactly the failure that looks like nothing at all --
it would shift the simulated baseline **by batch**, quietly, for every gene.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from watteg.dispersion import fit_dispersions

FAILURES: list[str] = []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path)
    parser.add_argument("reference", type=Path, help="output of workflow/dump_r_sim_input.R")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pysceptre" / "scripts"))
    from pysceptre.pipeline.discovery import (
        _GENE_BATCH_WIDTH,
        _GENE_CHUNK_MEMORY_GB,
        fit_all_genes,
    )
    from sceptre_io import load_export

    # cells_in_use and the export's covariate matrix: the same cells and the
    # same design sceptre's precomputation used. A different design silently
    # answers a different question, which is why this is not parameterised.
    export = load_export(args.export)
    genes = (args.reference / "fitted_coefs_genes.txt").read_text().split()
    names = (args.reference / "fitted_coefs_names.txt").read_text().splitlines()
    theirs = np.array(
        [
            [float(v) for v in line.split("\t")]
            for line in (args.reference / "fitted_coefs_full.txt").read_text().splitlines()
        ]
    )

    if list(export.metadata["covariate_names"]) != names:
        print("FAIL  covariate columns differ between the export and R's cache")
        print(f"        export: {list(export.metadata['covariate_names'])}")
        print(f"        R:      {names}")
        return 1
    print(f"covariate columns agree, in order ({len(names)})")

    idx = {g: i for i, g in enumerate(export.gene_ids)}
    fits = fit_all_genes(
        export.response_matrix,
        genes,
        export.covariate_matrix,
        chunk_memory_gb=_GENE_CHUNK_MEMORY_GB,
        gene_rows=[idx[g] for g in genes],
        batch_width=_GENE_BATCH_WIDTH,
        n_jobs=args.n_jobs,
    )
    ours = np.array([fits[g].fitted_coefs for g in genes])
    print(f"fitted {ours.shape[0]} genes x {ours.shape[1]} coefficients\n")

    for j, name in enumerate(names):
        a, b = ours[:, j], theirs[:, j]
        rel = np.abs(a - b) / np.maximum(np.abs(b), np.finfo(float).tiny)
        ok = bool(np.all(rel < args.rtol))
        if not ok:
            FAILURES.append(name)
        print(
            f"{'PASS' if ok else 'FAIL'}  {name:<28} max rel {rel.max():.2e}  "
            f"median {np.median(rel):.2e}"
        )

    # What the coefficients are actually for: the per-cell expected count. A
    # per-coefficient tolerance can pass while the fitted values drift, because
    # the design columns are correlated, so the quantity the simulation uses is
    # checked directly as well.
    X = export.covariate_matrix
    mu_ours = np.exp(X @ ours.T)
    mu_theirs = np.exp(X @ theirs.T)
    rel = np.abs(mu_ours - mu_theirs) / mu_theirs
    ok = bool(np.all(rel < args.rtol))
    if not ok:
        FAILURES.append("fitted values")
    print(
        f"\n{'PASS' if ok else 'FAIL'}  exp(X . beta), every gene x cell    "
        f"max rel {rel.max():.2e}  median {np.median(rel):.2e}  "
        f"({rel.size:,} values)"
    )

    _, diagnostics = fit_dispersions(
        export.response_matrix, export.gene_ids, X, genes, n_jobs=args.n_jobs
    )
    for kind, affected in diagnostics.items():
        if affected:
            print(f"  {kind}: {len(affected)} gene(s) -- {', '.join(affected[:5])}")

    print("\n" + ("ALL PASS" if not FAILURES else f"{len(FAILURES)} FAILED: {FAILURES}"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
