#!/usr/bin/env python3
"""Phase 2's gate: the Python expression statistics against R's own output.

R's reference is a `sim_input.rds` dumped by `workflow/dump_r_sim_input.R` --
for day0, the one the day0 sweep was actually run on. Python reads the
`--all-cells` pysceptre export of the same sceptre object, so both sides see the
same counts over the same 586,309 cells and any difference is the arithmetic.

    workflow/compare_expression_stats.py <export> <r_reference_dir>

Three quantities are compared, and `dispersion` is deliberately not one of them:
it comes from a model fit rather than from this arithmetic, and
`workflow/compare_dispersion.py` reports it as a distribution rather than as a
pass or a fail.

**The tolerance is 1e-10, and the reason is the reference's platform rather
than the port.** Against an R run on the *same machine* this port is
bit-identical: the per-gene geometric means agree exactly and the per-cell
medians to one ulp. Against the published `sim_input.rds` both differ by up to
6.2e-12, because R's `sum()` accumulates in `LDOUBLE` -- 80-bit extended
precision on the x86 cluster the sweep ran on, plain `double` on arm64, where
`.Machine$sizeof.longdouble == 8`. Measured: a local R run reproduces only
30 of 20,000 published size factors bit for bit, and differs by the same 6.2e-12
this port does. So the published outputs cannot be reproduced exactly by R
either, and 1e-10 is what "reproduces R" can mean across platforms.

The `bit-identical` column is the part to watch: it should read 100 % against a
same-platform reference, and near 0 % against the published one. A tolerance set
loosely by eye would have hidden the bug this comparison caught -- the export
stores counts as `uint16`, and `np.log` of a 16-bit integer array returns
**float32**, which moved the size factors by 2.5e-7.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from watteg.expression import compute_expression_stats

RTOL = 1e-10

FAILURES: list[str] = []


def read_floats(path: Path) -> np.ndarray:
    return np.array([float(x) for x in path.read_text().split()])


def report(label: str, ours: np.ndarray, theirs: np.ndarray, rtol: float = RTOL) -> None:
    ours = np.asarray(ours, dtype=float)
    theirs = np.asarray(theirs, dtype=float)
    if ours.shape != theirs.shape:
        print(f"FAIL  {label}: shape {ours.shape} vs R's {theirs.shape}")
        FAILURES.append(label)
        return
    rel = np.abs(ours - theirs) / np.maximum(np.abs(theirs), np.finfo(float).tiny)
    exact = int((ours == theirs).sum())
    ok = bool(np.allclose(ours, theirs, rtol=rtol, atol=0))
    print(
        f"{'PASS' if ok else 'FAIL'}  {label:30s} n={ours.size:>9,}  "
        f"max rel {rel.max():.2e}  bit-identical {exact:,}/{ours.size:,}"
    )
    if not ok:
        FAILURES.append(label)
        for i in np.argsort(rel)[-3:][::-1]:
            print(f"        [{i}] ours {ours[i]!r}  R {theirs[i]!r}  rel {rel[i]:.3e}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path, help="an --all-cells pysceptre export")
    parser.add_argument("reference", type=Path, help="output of workflow/dump_r_sim_input.R")
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pysceptre" / "scripts"))
    from sceptre_io import load_export

    export = load_export(args.export, all_cells=True)
    counts = export.response_matrix
    print(f"export: {counts.shape[0]} genes x {counts.shape[1]:,} cells, {counts.nnz:,} nonzeros")
    if counts.shape[1] == export.metadata.get("n_cells_in_use"):
        print(
            "  WARNING: this export has no QC-failed cells. R computes the statistics over every "
            "cell in the object, so a comparison against it needs an --all-cells export."
        )

    stats = compute_expression_stats(counts)
    print(f"  density {stats.density:.1%}")

    ref_genes = (args.reference / "genes.txt").read_text().split()
    print(f"reference: {len(ref_genes)} genes\n")

    gene_index = {g: i for i, g in enumerate(export.gene_ids)}
    missing = [g for g in ref_genes if g not in gene_index]
    if missing:
        print(f"FAIL  {len(missing)} of R's genes are absent from the export, e.g. {missing[:3]}")
        FAILURES.append("gene alignment")
        return 1
    rows = np.array([gene_index[g] for g in ref_genes])

    report(
        "row_data$mean", stats.normalized_mean[rows], read_floats(args.reference / "row_mean.txt")
    )
    report(
        "row_data$average_expression",
        stats.average_expression_all_cells[rows],
        read_floats(args.reference / "row_average_expression_all_cells.txt"),
    )
    report(
        "col_data$size_factors",
        stats.size_factors,
        read_floats(args.reference / "col_size_factors.txt"),
    )

    print("\n" + ("ALL PASS" if not FAILURES else f"{len(FAILURES)} FAILED: {FAILURES}"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
