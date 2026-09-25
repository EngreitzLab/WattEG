"""Turn simulation results into a power estimate per pair.

    watteg-compute-power --partials partials/ --threshold-file discovery_threshold.txt \
        --out power_es0.15.tsv
    watteg-compute-power --simulations sims/ --threshold-file discovery_threshold.txt \
        --out power_es0.15.tsv

Two inputs, one table. `--partials` takes the per-pair counts each simulation task writes when it
holds all simulations of its pairs (`watteg-run-power-simulation --partials-out`), and adds them
up. `--simulations` takes one row per (pair, simulation), as TSV or Parquet. Both go through
`watteg.power.power_from_counts`, so they give the same table, byte for byte, when each pair's
simulations came from one task.

Both take files or directories; a directory expands to the TSV and Parquet files inside it, sorted,
so row order does not depend on the filesystem. See `watteg/power.py` for what power means here and
why the interval is Wilson's.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from watteg.power import compute_power, merge_power_counts, power_from_counts

SUFFIXES = (".parquet", ".tsv", ".tsv.gz")


def expand(paths: list[Path]) -> list[Path]:
    out: list[Path] = []
    for path in paths:
        if path.is_dir():
            out.extend(sorted(p for p in path.iterdir() if p.name.endswith(SUFFIXES)))
        else:
            out.append(path)
    return out


def read_one(path: Path) -> pd.DataFrame:
    """One simulation output, its floats exactly as written.

    `round_trip`, not pandas' default parser: the default returns a neighbouring float for about
    half of the p-values and fold changes a TSV holds (measured, 285,883 of 500,000), so a mean or a
    value at the threshold would depend on the parser rather than on the simulation.
    """
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--simulations", nargs="+", type=Path, help="one row per (pair, simulation)"
    )
    source.add_argument(
        "--partials",
        nargs="+",
        type=Path,
        help="per-pair counts from watteg-run-power-simulation --partials-out",
    )
    parser.add_argument("--threshold-file", type=Path)
    parser.add_argument(
        "--alpha", type=float, help="set the threshold explicitly instead of reading it"
    )
    parser.add_argument("--conf-level", type=float, default=0.95)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    if (args.threshold_file is None) == (args.alpha is None):
        raise SystemExit("pass exactly one of --threshold-file and --alpha")
    threshold = (
        args.alpha if args.alpha is not None else float(args.threshold_file.read_text().split()[0])
    )

    inputs = args.partials or args.simulations
    paths = expand(inputs)
    if not paths:
        raise SystemExit(f"{'--partials' if args.partials else '--simulations'} matched no files")
    frame = pd.concat([read_one(p) for p in paths], ignore_index=True)
    if args.partials:
        # The counts were taken at the threshold the simulation tasks were given; a different one
        # here would silently mean a different question.
        used = frame["threshold"].unique() if "threshold" in frame else []
        if len(used) != 1 or used[0] != threshold:
            raise SystemExit(
                f"the partial counts were taken at threshold(s) {list(used)}, not {threshold!r}"
            )
        print(
            f"read counts for {len(frame):,} pair chunk(s) from {len(paths)} file(s), "
            f"threshold {threshold:.6g}"
        )
        power = power_from_counts(merge_power_counts([frame]), conf_level=args.conf_level)
    else:
        print(
            f"read {len(frame):,} simulation rows from {len(paths)} file(s), "
            f"threshold {threshold:.6g}"
        )
        power = compute_power(frame, threshold, conf_level=args.conf_level)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    power.to_csv(args.out, sep="\t", index=False)

    width = power["power_ci_high"] - power["power_ci_low"]
    print(f"wrote {len(power):,} pairs to {args.out}")
    print(f"  simulations per pair: {power['n_reps'].min()}-{power['n_reps'].max()}")
    print(
        f"  mean power {power['power'].mean():.3f} | at 0: {(power['power'] == 0).sum():,} "
        f"| at 1: {(power['power'] == 1).sum():,} "
        f"| in (0.1, 0.9): {power['power'].between(0.1, 0.9, 'neither').sum():,}"
    )
    print(f"  median 95% CI width {width.median():.3f} (widest {width.max():.3f})")
    if power["n_reps"].min() < 100:
        print(
            "  note: below ~100 simulations a per-pair estimate is coarse; see "
            "docs/choosing-num-replicates.md"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
