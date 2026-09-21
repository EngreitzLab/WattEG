"""One row per pair, across every effect size in a sweep.

    watteg-summarize-power --power power_es0.05.tsv power_es0.15.tsv ... \
        --sim-input prepared/sim_input.h5 --out power_summary.tsv

Each effect size contributes `power_at_effect_size_<label>` and its interval,
and the table gains the smallest tested effect size from which each pair is
detectable. See `watteg/power.py::min_detectable_effect_size` for the suffix
rule and why it is not "the first effect size that clears".

**Three bases for that column, because the interval on power carries through.**
Power is monotone in effect size, so thresholding a *lower* bound on power
yields a *larger* effect size -- the direction inverts:

    min_detectable_effect_size           from power          point estimate
    min_detectable_effect_size_ci_low    from power_ci_high  optimistic edge
    min_detectable_effect_size_ci_high   from power_ci_low   conservative edge

For "this pair was powered well enough that a null result means something", the
conservative edge is the one to use.

`--sim-input` joins each gene's dispersion and normalised mean. The theory
behind the power curve is `SE^2 ~ (1/n_pert_cells)(1/mu + dispersion)`, so any
covariate model of power needs both, and only the expression reaches this table
through the simulation output. Joining here costs one read; emitting it per
replicate would have cost a re-run of the sweep.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from watteg.power import effect_label, min_detectable_effect_size
from watteg.sim_input import read_sim_input

KEY = ["grna_target", "response_id"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--power", nargs="+", type=Path, required=True)
    parser.add_argument("--sim-input", type=Path)
    parser.add_argument("--power-threshold", type=float, default=0.8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    tables = [pd.read_csv(p, sep="\t") for p in args.power]
    for path, table in zip(args.power, tables, strict=True):
        if "effect_size" not in table.columns:
            raise SystemExit(f"{path} has no effect_size column, so it cannot be labelled")
        if table["effect_size"].nunique() != 1:
            raise SystemExit(f"{path} mixes effect sizes")
    # Ascending, because the suffix rule walks down from the largest.
    order = np.argsort([t["effect_size"].iloc[0] for t in tables])
    tables = [tables[i] for i in order]
    effect_sizes = [float(t["effect_size"].iloc[0]) for t in tables]
    if len(set(effect_sizes)) != len(effect_sizes):
        raise SystemExit(f"two inputs have the same effect size: {effect_sizes}")
    print(f"{len(tables)} effect sizes: {', '.join(f'{e:g}' for e in effect_sizes)}")

    out = (
        pd.concat([t[KEY] for t in tables], ignore_index=True)
        .drop_duplicates()
        .sort_values(KEY)
        .reset_index(drop=True)
    )

    # Shared across effect sizes: taken from the first table that has the pair.
    for column in ("mean_pert_cells", "average_expression_all_cells"):
        present = [t for t in tables if column in t.columns]
        if not present:
            continue
        merged = out[KEY].copy()
        for table in present:
            values = out.merge(table[[*KEY, column]], on=KEY, how="left")[column]
            merged[column] = merged[column].fillna(values) if column in merged else values
        out[column] = merged[column]

    if args.sim_input:
        sim = read_sim_input(args.sim_input)
        row = sim.row_data
        absent = out.loc[~out["response_id"].isin(row.index), "response_id"].nunique()
        if absent:
            print(
                f"  note: {absent} gene(s) are absent from --sim-input; their per-gene "
                "columns are NaN"
            )
        # `dispersion` is 1/theta, not theta. The name matches sim_input's so the
        # two cannot be read as different quantities.
        out["gene_mean"] = out["response_id"].map(row["mean"]).to_numpy()
        out["dispersion"] = out["response_id"].map(row["dispersion"]).to_numpy()
        if (
            "average_expression_all_cells" not in out
            or out["average_expression_all_cells"].isna().all()
        ):
            out["average_expression_all_cells"] = (
                out["response_id"].map(row["average_expression_all_cells"]).to_numpy()
            )

    power_columns = []
    for table, effect_size in zip(tables, effect_sizes, strict=True):
        base = f"power_at_effect_size_{effect_label(effect_size)}"
        power_columns.append(base)
        joined = out.merge(table, on=KEY, how="left", suffixes=("", "_new"))
        out[base] = joined["power"].to_numpy()
        for source, suffix in (
            ("power_ci_low", "_ci_low"),
            ("power_ci_high", "_ci_high"),
            ("n_reps", "_n_reps"),
        ):
            if source in table.columns:
                out[base + suffix] = joined[source].to_numpy()

    values = out[power_columns].to_numpy(dtype=float)
    out["min_detectable_effect_size"] = min_detectable_effect_size(
        values, effect_sizes, args.power_threshold
    )
    max_tested = np.where(np.isnan(values), -np.inf, np.array(effect_sizes)[None, :]).max(axis=1)
    max_tested = np.where(np.isneginf(max_tested), np.nan, max_tested)
    lo = [c + "_ci_low" for c in power_columns]
    hi = [c + "_ci_high" for c in power_columns]
    if all(c in out.columns for c in (*lo, *hi)):
        # Inverted on purpose: the optimistic edge of the effect size comes from
        # the optimistic (upper) edge of power.
        out["min_detectable_effect_size_ci_low"] = min_detectable_effect_size(
            out[hi].to_numpy(dtype=float), effect_sizes, args.power_threshold
        )
        out["min_detectable_effect_size_ci_high"] = min_detectable_effect_size(
            out[lo].to_numpy(dtype=float), effect_sizes, args.power_threshold
        )
    else:
        print("  note: no interval columns in the inputs; only the point estimate is derived")

    # Last, and after the interval columns, matching the R implementation's
    # order so the two outputs are diffable. Reported at all because a missing
    # min_detectable_effect_size is a statement about the effect sizes that were
    # run, not proof the pair is undetectable, and without this a reader cannot
    # tell which.
    out["max_effect_size_tested"] = max_tested

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    reached = out["min_detectable_effect_size"].notna().sum()
    print(f"wrote {len(out):,} pairs to {args.out}")
    print(
        f"  detectable at some tested effect size: {reached:,} of {len(out):,} "
        f"(power >= {args.power_threshold})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
