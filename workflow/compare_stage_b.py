#!/usr/bin/env python3
"""Stage B: the Python pipeline's power against the R pipeline's.

    workflow/compare_stage_b.py <python_sim.tsv> <r_power.tsv> \
        --threshold-file discovery_threshold.txt

This is the comparison that decides whether the Python path can replace the R
one. It reads a power table rather than driving R, but that table has to come
from a **fixed** R -- `r-implementation` at `b346296` or later -- and not from
the pre-2026-09-21 sweeps. Those were produced before the effect-size matrix
was centred correctly, so comparing against them measures the bug and reports
it as a difference between implementations. Stage B's first run did exactly
that, which is how the bug was found.

**What has to match, and what cannot.** The two implementations draw from
different random number generators, so no replicate corresponds to any other.
What has to agree is per-pair power, and the bar is the Monte Carlo noise
between two independent estimates of the same binomial proportion:
`sqrt(p_py(1-p_py)/n_py + p_r(1-p_r)/n_r)`, computed per pair from each side's
own replicate count rather than assumed equal.

**Run both sides on their own defaults.** This used to demand
`--expression-model size_factor` on the Python side, to match a reference that
predated the baseline change. With the reference regenerated that would invert:
it would take Python off the model R is now using. The rule underneath has not
changed -- compare like for like -- only which side had to move.

Three things are reported, because power alone is blind in two ways:

* the **distribution** of per-pair differences against that noise floor, which
  is the headline;
* whether the differences are **centred on zero**, because a systematic shift
  is what a real difference looks like while individual pairs crossing is what
  noise looks like;
* agreement on the **0.8 decision**, since that threshold is what the sweep is
  used for, and a pair either side of it is a different conclusion.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from watteg.power import compute_power

KEY = ["grna_target", "response_id"]
FAILURES: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {label}{'  ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(label)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_sim", type=Path, help="per-replicate output of the Python path")
    parser.add_argument("r_power", type=Path, help="a published power_es*.tsv")
    parser.add_argument("--threshold-file", type=Path, required=True)
    parser.add_argument("--power-threshold", type=float, default=0.8)
    args = parser.parse_args()

    threshold = float(args.threshold_file.read_text().split()[0])
    py = compute_power(pd.read_csv(args.python_sim, sep="\t"), threshold)
    r = pd.read_csv(args.r_power, sep="\t")

    merged = py.merge(r, on=KEY, suffixes=("_py", "_r"))
    check(
        "every simulated pair is in the reference",
        len(merged) == len(py),
        f"{len(merged)} of {len(py)} (reference holds {len(r):,})",
    )
    if merged.empty:
        print("\nno pairs in common")
        return 1

    n_py = merged["n_reps_py"].to_numpy(dtype=float)
    n_r = merged["n_reps_r"].to_numpy(dtype=float)
    a = merged["power_py"].to_numpy(dtype=float)
    b = merged["power_r"].to_numpy(dtype=float)
    d = a - b

    # Per pair, from each side's own estimate and replicate count.
    se = np.sqrt(a * (1 - a) / n_py + b * (1 - b) / n_r)
    # A pair at power 0 or 1 in both has se = 0 and agrees exactly; that is not a free pass, it is
    # the strongest possible agreement, so it is counted rather than excluded.
    within = np.abs(d) <= 2 * se + 1e-12

    print(f"\n{len(merged)} pairs, {int(n_py[0])} Python replicates against {int(n_r[0])} R's")
    print(f"  mean difference          {d.mean():+.4f}")
    print(f"  mean |difference|        {np.abs(d).mean():.4f}")
    print(f"  largest |difference|     {np.abs(d).max():.4f}")
    print(f"  median 2-se noise floor  {2 * np.median(se):.4f}")
    check(
        "per-pair power within 2 se of Monte Carlo noise",
        bool(within.mean() >= 0.95),
        f"{int(within.sum())}/{len(within)} ({within.mean():.1%}; ~95% expected by chance alone)",
    )

    # A shift of the mean is the signature of a real difference. Its standard error is the
    # quadrature sum over pairs, not the spread of the differences, because the pairs have
    # different precisions.
    se_mean = float(np.sqrt(np.sum(se**2)) / len(se))
    check(
        "no systematic shift",
        abs(d.mean()) <= 2 * se_mean + 1e-12,
        f"{d.mean():+.4f} against a 2-se bound of {2 * se_mean:.4f}",
    )

    # The decision the sweep exists to support.
    call_py = a >= args.power_threshold
    call_r = b >= args.power_threshold
    agree = call_py == call_r
    print(
        f"\n  at the {args.power_threshold} bar: {int(agree.sum())}/{len(agree)} agree "
        f"({agree.mean():.1%}); {int((call_py & ~call_r).sum())} Python-only, "
        f"{int((~call_py & call_r).sum())} R-only"
    )
    # Disagreements should be symmetric and confined to pairs near the bar.
    near = np.minimum(np.abs(a - args.power_threshold), np.abs(b - args.power_threshold)) <= 2 * se
    check(
        "every 0.8 disagreement is a pair sitting on the bar",
        bool(np.all(agree | near)),
        f"{int((~agree & ~near).sum())} away from it",
    )

    worst = np.argsort(np.abs(d) - 2 * se)[-5:][::-1]
    print("\n  furthest outside their own noise floor:")
    for i in worst:
        row = merged.iloc[i]
        print(
            f"    {row['grna_target'][:28]:<28} {row['response_id']:<10} "
            f"py {a[i]:.2f}  R {b[i]:.2f}  diff {d[i]:+.2f}  2se {2 * se[i]:.2f}"
        )

    print("\n" + ("ALL PASS" if not FAILURES else f"{len(FAILURES)} FAILED: {FAILURES}"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
