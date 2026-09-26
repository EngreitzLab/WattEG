"""Concatenate a sweep's per-split simulation output into one Parquet file.

    watteg-consolidate-replicates --simulations sims/ --out per_replicate/es0.15.parquet

The per-replicate output is the one thing power can be re-derived from, so it is
read repeatedly -- and as 1,000 gzipped TSVs that costs 30-90 s of parsing every
time. One Parquet file per effect size is seconds.

The split files stay in the work directory, so a failed consolidation loses
nothing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from watteg.cli.compute_power import expand, read_one

# Stored as categories: a sweep repeats a few thousand target and gene names
# across millions of rows, and dictionary encoding is what makes the file small.
CATEGORICAL = ("grna_target", "response_id")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulations", nargs="+", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--compression", default="zstd")
    args = parser.parse_args(argv)

    paths = expand(args.simulations)
    if not paths:
        raise SystemExit(f"--simulations matched no files: {args.simulations}")
    frame = pd.concat([read_one(p) for p in paths], ignore_index=True)

    if "effect_size" in frame.columns:
        sizes = frame["effect_size"].unique()
        if len(sizes) > 1:
            raise SystemExit(
                f"the input mixes effect sizes ({', '.join(map(str, sizes))}); consolidate "
                "one effect size per file, or power would be averaged across knockdown levels"
            )
    for column in CATEGORICAL:
        if column in frame.columns:
            frame[column] = frame[column].astype("category")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(args.out, compression=args.compression, index=False)
    size_mb = args.out.stat().st_size / 1e6
    print(f"wrote {len(frame):,} rows from {len(paths)} file(s) to {args.out} ({size_mb:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
