"""Turn per-replicate simulation results into a power estimate per pair.

    watteg-compute-power --simulations sims/ --threshold-file discovery_threshold.txt \
        --out power_es0.15.tsv

`--simulations` takes files or directories; a directory expands to the TSV and
Parquet files inside it, sorted, so row order does not depend on the filesystem.
See `watteg/power.py` for what power means here and why the interval is Wilson's.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from watteg.power import compute_power

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
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulations", nargs="+", type=Path, required=True)
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

    paths = expand(args.simulations)
    if not paths:
        raise SystemExit(f"--simulations matched no files: {args.simulations}")
    frame = pd.concat([read_one(p) for p in paths], ignore_index=True)
    print(
        f"read {len(frame):,} simulation rows from {len(paths)} file(s), threshold {threshold:.6g}"
    )

    power = compute_power(frame, threshold, conf_level=args.conf_level)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    power.to_csv(args.out, sep="\t", index=False)

    width = power["power_ci_high"] - power["power_ci_low"]
    print(f"wrote {len(power):,} pairs to {args.out}")
    print(f"  replicates per pair: {power['n_reps'].min()}-{power['n_reps'].max()}")
    print(
        f"  mean power {power['power'].mean():.3f} | at 0: {(power['power'] == 0).sum():,} "
        f"| at 1: {(power['power'] == 1).sum():,} "
        f"| in (0.1, 0.9): {power['power'].between(0.1, 0.9, 'neither').sum():,}"
    )
    print(f"  median 95% CI width {width.median():.3f} (widest {width.max():.3f})")
    if power["n_reps"].min() < 100:
        print(
            "  note: below ~100 replicates a per-pair estimate is coarse; see "
            "docs/choosing-num-replicates.md"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
