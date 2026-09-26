#!/usr/bin/env python3
"""The whole of `prepare_sim_input`, Python against R, on one sample.

    workflow/compare_sim_input.py <python_outdir> <r_outdir> <r_reference_dir>

`<python_outdir>` is what `watteg-prepare-sim-input` wrote, `<r_outdir>` is what
`src/prepare_sim_input.R` wrote, and `<r_reference_dir>` is that run's
`sim_input.rds` dumped by `workflow/dump_r_sim_input.R`.

`compare_expression_stats.py` and `compare_dispersion.py` cover the numbers this
step computes. What is left, and what this covers, is everything it *assembles*:
the two perturbation matrices, the pair table, the gRNA map, the threshold and
the analysis mode.

**Two differences here are by design and are checked as such, not waived.**

*The cell set.* R's sim_input spans every cell in the object; the Python one
spans `cells_in_use`. The simulation only ever touches cells that have
covariates, and R's extra columns are carried into a matrix sceptre then
discards. So the perturbation matrices are compared after mapping Python's
columns back to absolute cell positions through `cells_in_use` -- if a target's
cell set differs there, that is a real disagreement.

*The gRNA matrix's rows.* R builds `grna_perts` from
`@initial_grna_assignment_list`, which has one entry per row of the screen's
design table -- so a guide inside two overlapping elements appears twice, under
the same name, and R's row count is 45,463 against 43,736 distinct guides. It
also predates QC, so it holds memberships in cells the analysis dropped. Python
keeps one row per guide over `cells_in_use`. Both are compared on what they
agree about: per guide, the cells in `cells_in_use`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from watteg.sim_input import read_sim_input

FAILURES: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {label}{'  ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(label)


def read_triplets(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def cells_by_unit(rows: list[str], triplets: pd.DataFrame) -> dict[str, np.ndarray]:
    """R's matrix may repeat a row name; a guide's cells are the union of its rows."""
    out: dict[str, set] = {}
    row_name = np.array(rows, dtype=object)
    for i, j in zip(triplets["i"].to_numpy(), triplets["j"].to_numpy(), strict=True):
        out.setdefault(row_name[i], set()).add(int(j))
    return {unit: np.array(sorted(cells), dtype=np.int64) for unit, cells in out.items()}


def compare_perts(label, ours, our_rows, absolute, reference_dir, in_use_set) -> None:
    """`ours` is (units x cells_in_use); `absolute` maps its columns to the
    object's own cell positions, which is the space R's matrix is in."""
    ref_rows = (reference_dir / f"{label}_rows.txt").read_text().split()
    ref = cells_by_unit(ref_rows, read_triplets(reference_dir / f"{label}_triplets.tsv"))

    csr = ours.tocsr()
    mine = {
        unit: absolute[csr.indices[csr.indptr[i] : csr.indptr[i + 1]]]
        for i, unit in enumerate(our_rows)
    }

    shared = sorted(set(mine) & set(ref))
    mismatched = []
    for unit in shared:
        # R's cells outside cells_in_use have no counterpart on this side, and
        # are removed before comparing rather than counted as a difference.
        theirs = np.array([c for c in ref[unit] if c in in_use_set], dtype=np.int64)
        if not np.array_equal(np.sort(mine[unit]), theirs):
            mismatched.append(unit)
    check(
        f"{label}: cell sets agree",
        not mismatched,
        f"{len(shared):,} units compared, {len(mismatched)} differ"
        + (f", e.g. {mismatched[:2]}" if mismatched else ""),
    )
    # A unit R has and Python does not is only acceptable when it has nothing in
    # cells_in_use: R carries an all-zero row where Python omits the row, and an
    # all-zero row is never selected by anything. A unit with cells, though, is
    # a guide the simulation would not see.
    only_ref = sorted(set(ref) - set(mine))
    with_cells = [u for u in only_ref if any(c in in_use_set for c in ref[u])]
    if only_ref:
        print(
            f"        {len(only_ref):,} units only in R's matrix; "
            f"{len(only_ref) - len(with_cells):,} of them have no QC-passing cell"
        )
    check(
        f"{label}: no unit is missing QC-passing cells",
        not with_cells,
        f"{len(with_cells)} unit(s)" + (f", e.g. {with_cells[:2]}" if with_cells else ""),
    )
    only_mine = sorted(set(mine) - set(ref))
    if only_mine:
        print(f"        {len(only_mine):,} units only in Python's (e.g. {only_mine[:2]})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_outdir", type=Path)
    parser.add_argument("r_outdir", type=Path)
    parser.add_argument("reference", type=Path, help="output of workflow/dump_r_sim_input.R")
    args = parser.parse_args()

    sim = read_sim_input(args.python_outdir / "sim_input.h5")
    print(f"python sim_input: {sim.describe()}")
    ref_genes = (args.reference / "genes.txt").read_text().split()
    print(f"r sim_input:      {len(ref_genes)} genes\n")

    check("genes: same set", set(sim.genes) == set(ref_genes), f"{len(sim.genes)}")
    check("genes: same order", list(sim.genes) == list(ref_genes))

    # --- the two TSVs -----------------------------------------------------------------
    for name, keys in (
        ("pairs.tsv", ["grna_target", "response_id"]),
        ("grna_targets.tsv", ["grna_id", "grna_target"]),
    ):
        ours = pd.read_csv(args.python_outdir / name, sep="\t")
        theirs = pd.read_csv(args.r_outdir / name, sep="\t")
        same_rows = set(map(tuple, ours[keys].to_numpy())) == set(
            map(tuple, theirs[keys].to_numpy())
        )
        check(
            f"{name}: same rows",
            same_rows and len(ours) == len(theirs),
            f"{len(ours):,} vs {len(theirs):,}",
        )
        check(
            f"{name}: byte-identical",
            (args.python_outdir / name).read_bytes() == (args.r_outdir / name).read_bytes(),
        )

    # --- the threshold and the analysis mode ------------------------------------------
    ours = (args.python_outdir / "discovery_threshold.txt").read_text().strip()
    theirs = (args.r_outdir / "discovery_threshold.txt").read_text().strip()
    check("discovery_threshold.txt", ours == theirs, f"{ours} vs {theirs}")

    # analysis_mode.tsv postdates the published day0 run, so a reference that
    # lacks it is an old reference rather than a failure.
    r_mode_path = args.r_outdir / "analysis_mode.tsv"
    if not r_mode_path.exists():
        print(f"SKIP  analysis_mode.tsv: absent from {args.r_outdir} (predates that output)")
        r_mode_path = None
    r_mode = (
        dict(line.split("\t", 1) for line in r_mode_path.read_text().splitlines() if "\t" in line)
        if r_mode_path is not None
        else {}
    )
    py_mode = dict(
        line.split("\t", 1)
        for line in (args.python_outdir / "analysis_mode.tsv").read_text().splitlines()
        if "\t" in line
    )
    for key in sorted(set(r_mode) & set(py_mode)):
        check(f"analysis_mode: {key}", r_mode[key] == py_mode[key], f"{py_mode[key]}")

    # --- the perturbation matrices ----------------------------------------------------
    absolute = sim.cells_in_use
    in_use_set = set(absolute.tolist())
    print()
    compare_perts("cre_perts", sim.cre_perts, sim.target_ids, absolute, args.reference, in_use_set)
    compare_perts("grna_perts", sim.grna_perts, sim.grna_ids, absolute, args.reference, in_use_set)

    print("\n" + ("ALL PASS" if not FAILURES else f"{len(FAILURES)} FAILED: {FAILURES}"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
