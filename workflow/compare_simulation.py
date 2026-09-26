#!/usr/bin/env python3
"""Do the two implementations give the same answer?

    workflow/compare_simulation.py <python.tsv> <r.tsv> --threshold-file discovery_threshold.txt

**They cannot agree draw for draw and it would be wrong to want them to.** The
two use different random number generators, so a given replicate's counts and
its CRT draws differ. What has to agree is the quantity the pipeline reports:
per-pair power, the fraction of replicates in which the screen's own test would
have called the association.

So this is a statistical comparison, and the bar is the Monte Carlo noise
between two correct runs rather than a tolerance. Two independent estimates of
the same binomial `p` over `n` replicates differ with standard deviation
`sqrt(2p(1-p)/n)`: at `p = 0.5` and `n = 40` that is 0.11, which is what
"agreement" can mean at this replicate count. Raising the replicate count is
the only thing that tightens it.

Reported alongside power, because power at the extremes hides everything: the
per-pair median p-value on a log scale, which is sensitive where power is
saturated at 0 or 1, and the fold change, which is the sharpest signal
available because it is an absolute quantity -- both implementations should
land on `log2(1 - effect_size)`, whatever their draws.

**Every bound here is derived from the data's own spread, per pair.** A first
version of this script asked the fold changes to agree to 0.05 in log2, a
number picked by eye, and one pair missed it at 0.0883 -- against its own
two-standard-error bound of 0.1286, because that gene's per-replicate spread
is 0.24 while another's is 0.02. A fixed tolerance is either too tight for the
noisy genes or too loose for the quiet ones, and there is no single value that
is neither.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FAILURES: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {label}{'  ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(label)


def power_per_pair(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """The pipeline's definition: a replicate counts only if the test called it
    AND the simulated perturbation reduced expression."""
    called = (frame["p_value"] < threshold) & (frame["log_2_fold_change"] < 0)
    return (
        frame.assign(called=called)
        .groupby(["grna_target", "response_id"])
        .agg(
            power=("called", "mean"),
            reps=("called", "size"),
            median_log10_p=("p_value", lambda p: float(np.median(np.log10(np.maximum(p, 1e-300))))),
            median_log2fc=("log_2_fold_change", "median"),
            sd_log2fc=("log_2_fold_change", "std"),
        )
        .reset_index()
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_tsv", type=Path)
    parser.add_argument("r_tsv", type=Path)
    parser.add_argument("--threshold-file", type=Path, required=True)
    args = parser.parse_args()

    threshold = float(args.threshold_file.read_text().strip())
    py = pd.read_csv(args.python_tsv, sep="\t")
    r = pd.read_csv(args.r_tsv, sep="\t")
    print(f"threshold {threshold:.6g}")
    print(f"python: {len(py):,} rows, {py['rep'].nunique()} replicates")
    print(f"r:      {len(r):,} rows, {r['rep'].nunique()} replicates\n")

    a, b = power_per_pair(py, threshold), power_per_pair(r, threshold)
    key = ["grna_target", "response_id"]
    check(
        "the same pairs were simulated",
        set(map(tuple, a[key].to_numpy())) == set(map(tuple, b[key].to_numpy())),
        f"{len(a)} pairs",
    )
    merged = a.merge(b, on=key, suffixes=("_py", "_r"))
    n = int(merged["reps_py"].iloc[0])
    check("the same replicate count", (merged["reps_py"] == merged["reps_r"]).all(), f"{n}")

    d = merged["power_py"] - merged["power_r"]
    # The noise floor for the comparison, per pair, from R's estimate.
    p = merged["power_r"].clip(1e-9, 1 - 1e-9)
    sd = np.sqrt(2 * p * (1 - p) / n)
    within = np.abs(d) <= 2 * sd + 1e-12
    print(f"\npower, {len(merged)} pairs at {n} replicates")
    print(f"  mean difference        {d.mean():+.4f}   (expected ~0)")
    print(f"  mean |difference|      {np.abs(d).mean():.4f}")
    print(f"  largest |difference|   {np.abs(d).max():.4f}")
    print(f"  noise floor, 2 sd      {2 * sd.mean():.4f} on average")
    check(
        "every pair within 2 sd of Monte Carlo noise",
        bool(within.all()),
        f"{int(within.sum())}/{len(within)}",
    )
    # A mean shift is what a systematic difference looks like; individual pairs
    # crossing is what noise looks like.
    se_mean = float(np.sqrt(np.sum(sd**2)) / len(sd))
    check(
        "no systematic shift in power",
        abs(d.mean()) <= 2 * se_mean + 1e-12,
        f"{d.mean():+.4f} against a 2-se bound of {2 * se_mean:.4f}",
    )

    print("\nfold change in log2 -- an absolute quantity, so it is the sharpest check")
    fc = merged["median_log2fc_py"] - merged["median_log2fc_r"]
    # The standard error of a median is about 1.2533 sd / sqrt(n) for a normal,
    # and these two are independent, so their difference has the pooled one.
    fc_se = 1.2533 * np.sqrt((merged["sd_log2fc_py"] ** 2 + merged["sd_log2fc_r"] ** 2) / n)
    print(f"  mean difference        {fc.mean():+.4f}")
    print(f"  largest |difference|   {np.abs(fc).max():.4f}")
    print(f"  its own 2-se bound     {2 * fc_se[np.abs(fc).idxmax()]:.4f}")
    check(
        "every pair's fold change within 2 se of its own replicate spread",
        bool((np.abs(fc) <= 2 * fc_se + 1e-12).all()),
        f"{int((np.abs(fc) <= 2 * fc_se + 1e-12).sum())}/{len(fc)}",
    )

    # The absolute check: whatever the draws, a 15% knockdown should land on
    # log2(0.85). This is what would catch both implementations being wrong in
    # the same direction, which no comparison between them can.
    effect_size = float(py["effect_size"].iloc[0])
    wanted = np.log2(1.0 - effect_size)
    for label, column in (("python", "median_log2fc_py"), ("r", "median_log2fc_r")):
        off = merged[column] - wanted
        check(
            f"{label} hits the requested knockdown (log2 {wanted:.4f})",
            bool((np.abs(off) <= 2 * fc_se + 1e-12).all()),
            f"median {merged[column].median():+.4f}, largest miss {np.abs(off).max():.4f}",
        )

    print("\nmedian p-value per pair, log10 -- sensitive where power saturates")
    lp = merged["median_log10_p_py"] - merged["median_log10_p_r"]
    print(f"  mean difference        {lp.mean():+.3f} decades")
    print(f"  largest |difference|   {np.abs(lp).max():.3f} decades")

    saturated = int(((merged["power_r"] == 0) | (merged["power_r"] == 1)).sum())
    if saturated:
        print(
            f"\n  NOTE: {saturated} of {len(merged)} pairs have power exactly 0 or 1 in R, "
            "where power cannot distinguish the two implementations. Read the fold change and "
            "the p-value rows for those."
        )

    print("\n" + ("ALL PASS" if not FAILURES else f"{len(FAILURES)} FAILED: {FAILURES}"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
