"""Split the discovery pairs into per-task chunks, balanced by cost.

    watteg-split-pairs --pairs pairs.tsv --outdir splits/ --n-splits 480

A target's pairs cannot be separated: the simulation draws one count matrix per
(target, replicate) and tests every one of that target's genes against it, so
splitting a target would simulate it twice. Targets are therefore the unit, and
the task is bin packing them.

**Weighted by pairs plus a per-target overhead**, because a task's cost is not
proportional to its pairs. The measured model is `intercept + slope x pairs`,
and on day0 the intercept is a real share of a small target's cost, so weighting
by pairs alone over-fills the splits that hold many small targets.

Longest-processing-time-first: sort by weight descending, then put each target
in whichever split is currently lightest. Deterministic, and within a small
constant factor of optimal for this shape.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def assign_splits(weights: pd.Series, n_splits: int) -> pd.Series:
    """Target -> split number (1-based), by longest-processing-time-first."""
    load = [0.0] * n_splits
    assignment = {}
    for target, weight in weights.sort_values(ascending=False).items():
        lightest = min(range(n_splits), key=lambda k: load[k])
        assignment[target] = lightest + 1
        load[lightest] += float(weight)
    return pd.Series(assignment, name="split")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--n-splits", type=int, required=True)
    parser.add_argument("--prefix", default="split_")
    parser.add_argument(
        "--target-overhead",
        type=float,
        default=2.0,
        help="pairs-equivalent fixed cost of a target, from the measured cost model "
        "[default %(default)s]",
    )
    args = parser.parse_args(argv)

    if args.target_overhead < 0:
        raise SystemExit("--target-overhead must not be negative")

    pairs = pd.read_csv(args.pairs, sep="\t")
    for column in ("grna_target", "response_id"):
        if column not in pairs.columns:
            raise SystemExit(f"{args.pairs} has no {column!r} column")
    if pairs.empty:
        raise SystemExit(f"{args.pairs} contains no pairs")

    per_target = pairs.groupby("grna_target").size()
    if args.n_splits > len(per_target):
        raise SystemExit(
            f"--n-splits ({args.n_splits}) exceeds the number of targets ({len(per_target)}). "
            "Every split must hold at least one target; lower --n-splits."
        )
    print(
        f"{len(pairs):,} pairs over {len(per_target):,} targets "
        f"({per_target.min()}-{per_target.max()} each, median {per_target.median():g})"
    )

    split_of = assign_splits(per_target + args.target_overhead, args.n_splits)
    pairs = pairs.assign(split=pairs["grna_target"].map(split_of))

    args.outdir.mkdir(parents=True, exist_ok=True)
    # Zero-padded so lexicographic order matches numeric order, which keeps a
    # workflow engine's channel ordering and a manual listing predictable.
    width = max(2, len(str(args.n_splits)))
    written = 0
    sizes = []
    for k in range(1, args.n_splits + 1):
        chunk = pairs.loc[pairs["split"] == k, ["grna_target", "response_id"]]
        chunk = chunk.sort_values(["grna_target", "response_id"])
        chunk.to_csv(args.outdir / f"{args.prefix}{k:0{width}d}.tsv", sep="\t", index=False)
        written += len(chunk)
        sizes.append(len(chunk))

    if written != len(pairs):
        raise SystemExit(f"wrote {written} pairs but read {len(pairs)}; the split lost rows")
    print(f"wrote {args.n_splits} splits to {args.outdir}")
    print(
        f"  pairs per split: {min(sizes)}-{max(sizes)}, "
        f"imbalance (max/min) {max(sizes) / max(min(sizes), 1):.3f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
