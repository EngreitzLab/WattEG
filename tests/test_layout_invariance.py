"""A pair's power must not depend on how the work was split up.

How many splits a cluster wanted and how many replicates fit in a task are
operational choices. If either reached the numbers, two runs of the same sweep
would disagree for reasons that have nothing to do with the screen, and no
output would say so.

    WATTEG_PREPARED=path/to/prepared pytest -m realdata

Opt-in, because it needs a prepared dataset and runs the real engine. The
seed-level half of the same contract is in `test_simulation_core.py` and runs
everywhere; this is the end-to-end half, which is what would catch the contract
being broken *downstream* of the seed -- by a per-task RNG, a shared generator,
or anything that consumes draws in an order that depends on the chunking.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytestmark = pytest.mark.realdata

PREPARED = os.environ.get("WATTEG_PREPARED")
KEY = ["grna_target", "response_id", "rep"]


@pytest.fixture(scope="module")
def prepared() -> Path:
    if not PREPARED:
        pytest.skip("set WATTEG_PREPARED to a watteg-prepare-sim-input output directory")
    return Path(PREPARED)


def simulate(
    prepared: Path,
    out: Path,
    split: Path,
    *,
    reps: int,
    offset: int,
    seed: int = 7,
    n_jobs: int = 1,
    permutations: str = "per-replicate",
    nulls: str = "scan",
    driver: str = "engine",
    null_fits: str = "refit",
    extra: tuple = (),
):
    subprocess.run(
        [
            sys.executable,
            "-m",
            "watteg.cli.run_power_simulation",
            "--prepared",
            str(prepared),
            "--pairs",
            str(split),
            "--effect-size",
            "0.15",
            "--reps",
            str(reps),
            "--rep-offset",
            str(offset),
            "--seed",
            str(seed),
            "--n-jobs",
            str(n_jobs),
            "--permutations",
            permutations,
            "--nulls",
            nulls,
            "--driver",
            driver,
            "--null-fits",
            null_fits,
            *extra,
            "--out",
            str(out),
        ],
        check=True,
        capture_output=True,
    )
    return pd.read_csv(out, sep="\t")


@pytest.fixture(scope="module")
def split(prepared: Path, tmp_path_factory) -> Path:
    """One target's pairs, which is enough: the contract is per target."""
    pairs = pd.read_csv(prepared / "pairs.tsv", sep="\t")
    target = pairs["grna_target"].iloc[0]
    path = tmp_path_factory.mktemp("split") / "split_01.tsv"
    pairs[pairs["grna_target"] == target].to_csv(path, sep="\t", index=False)
    return path


def test_replicate_chunking_does_not_change_a_single_draw(prepared, split, tmp_path):
    """Four replicates in one task against two tasks of two.

    Byte-identical, not merely statistically similar: the four replicates are
    keyed on (seed, target, rep, effect size) and nothing else, so the chunk
    boundary has nowhere to enter.
    """
    one = simulate(prepared, tmp_path / "one.tsv", split, reps=4, offset=0)
    first = simulate(prepared, tmp_path / "a.tsv", split, reps=2, offset=0)
    second = simulate(prepared, tmp_path / "b.tsv", split, reps=2, offset=2)
    two = pd.concat([first, second])

    one = one.sort_values(KEY).reset_index(drop=True)
    two = two.sort_values(KEY).reset_index(drop=True)
    assert len(one) == len(two)
    assert one.equals(two)


def test_worker_count_does_not_change_a_single_draw(prepared, split, tmp_path):
    """One worker against three, byte for byte and row for row.

    The task runs (target, replicate) units in parallel and they finish in any
    order; each draws from its own seeded stream and the results are gathered in
    submission order, so neither the draws nor the row order can move.
    """
    serial = simulate(prepared, tmp_path / "j1.tsv", split, reps=4, offset=0, n_jobs=1)
    parallel = simulate(prepared, tmp_path / "j3.tsv", split, reps=4, offset=0, n_jobs=3)
    assert (tmp_path / "j1.tsv").read_bytes() == (tmp_path / "j3.tsv").read_bytes()
    assert serial.equals(parallel)


def test_a_different_seed_does_change_the_draws(prepared, split, tmp_path):
    """The companion to the test above: invariance must not come from the seed
    being ignored."""
    a = simulate(prepared, tmp_path / "s7.tsv", split, reps=2, offset=0, seed=7)
    b = simulate(prepared, tmp_path / "s8.tsv", split, reps=2, offset=0, seed=8)
    assert (
        not a.sort_values(KEY)
        .reset_index(drop=True)["p_value"]
        .equals(b.sort_values(KEY).reset_index(drop=True)["p_value"])
    )


def test_per_target_permutations_keep_the_counts_and_change_the_nulls(prepared, split, tmp_path):
    """One permutation set per target changes only the permutations.

    The fold change is a function of the simulated counts alone, so equal fold changes prove the
    two modes simulated the same counts; different p-values prove the permutations moved.
    """
    fresh = simulate(prepared, tmp_path / "rep.tsv", split, reps=3, offset=0)
    shared = simulate(
        prepared, tmp_path / "tgt.tsv", split, reps=3, offset=0, permutations="per-target"
    )
    fresh = fresh.sort_values(KEY).reset_index(drop=True)
    shared = shared.sort_values(KEY).reset_index(drop=True)
    assert fresh["log_2_fold_change"].equals(shared["log_2_fold_change"])
    assert not fresh["p_value"].equals(shared["p_value"])


def test_per_target_permutations_are_invariant_to_chunking_and_workers(prepared, split, tmp_path):
    """The per-target seed depends on (seed, target, effect size), never on the replicate range
    or the worker count, so the layout contract holds in this mode too."""
    one = simulate(prepared, tmp_path / "a.tsv", split, reps=4, offset=0, permutations="per-target")
    first = simulate(
        prepared, tmp_path / "b.tsv", split, reps=2, offset=0, permutations="per-target"
    )
    second = simulate(
        prepared, tmp_path / "c.tsv", split, reps=2, offset=2, n_jobs=2, permutations="per-target"
    )
    two = pd.concat([first, second]).sort_values(KEY).reset_index(drop=True)
    assert one.sort_values(KEY).reset_index(drop=True).equals(two)


@pytest.mark.parametrize("permutations", ["per-replicate", "per-target"])
def test_sparse_and_scan_nulls_give_identical_output(prepared, split, tmp_path, permutations):
    """The sparse null route is a speed change only: byte-identical output to the scan route.

    It holds only while pysceptre's two routes compute the same segment sums, so this is also the
    test that catches a pysceptre update breaking that assumption.
    """
    scan = tmp_path / "scan.tsv"
    sparse = tmp_path / "sparse.tsv"
    simulate(prepared, scan, split, reps=4, offset=0, permutations=permutations, nulls="scan")
    simulate(prepared, sparse, split, reps=4, offset=0, permutations=permutations, nulls="sparse")
    assert scan.read_bytes() == sparse.read_bytes()


def test_the_fast_driver_writes_the_engines_bytes(prepared, split, tmp_path):
    """--driver fast is a layout change only: the file the engine writes with --permutations
    per-target --nulls sparse, byte for byte, at one worker and at two (which splits the target's
    simulations into chunks). This is also what catches a pysceptre update that the fast driver's
    reproduction of pysceptre's per-gene steps no longer matches."""
    engine = tmp_path / "engine.tsv"
    simulate(prepared, engine, split, reps=4, offset=0, permutations="per-target", nulls="sparse")
    for n_jobs in (1, 2):
        fast = tmp_path / f"fast_j{n_jobs}.tsv"
        simulate(
            prepared,
            fast,
            split,
            reps=4,
            offset=0,
            n_jobs=n_jobs,
            permutations="per-target",
            driver="fast",
        )
        assert engine.read_bytes() == fast.read_bytes()


def test_reused_null_fits_are_the_same_from_a_file_or_made_in_the_task(prepared, split, tmp_path):
    """--null-fits reuse: fits written by watteg-fit-null-models for more genes and simulations
    than the task needs, and fits the task makes for itself at one worker and at two, give the
    same file. A fit is keyed on (seed, gene, simulation) alone, so the layout cannot reach it."""
    pairs = pd.read_csv(prepared / "pairs.tsv", sep="\t")
    genes = list(pd.read_csv(split, sep="\t")["response_id"].unique())
    others = [g for g in pairs["response_id"].unique() if g not in genes][:4]
    wider = tmp_path / "wider.tsv"
    pd.DataFrame({"grna_target": "x", "response_id": genes + others}).to_csv(
        wider, sep="\t", index=False
    )
    fits = tmp_path / "fits.h5"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "watteg.cli.fit_null_models",
            "--prepared",
            str(prepared),
            "--pairs",
            str(wider),
            "--reps",
            "5",
            "--seed",
            "7",
            "--n-jobs",
            "2",
            "--out",
            str(fits),
        ],
        check=True,
        capture_output=True,
    )
    reuse = dict(permutations="per-target", driver="fast", null_fits="reuse")
    from_file = tmp_path / "file.tsv"
    simulate(
        prepared, from_file, split, reps=4, offset=0, **reuse, extra=("--null-fits-file", str(fits))
    )
    for n_jobs in (1, 2):
        own = tmp_path / f"own_j{n_jobs}.tsv"
        simulate(prepared, own, split, reps=4, offset=0, n_jobs=n_jobs, **reuse)
        assert from_file.read_bytes() == own.read_bytes()

    # And a file made under another seed is refused rather than used.
    with pytest.raises(subprocess.CalledProcessError) as caught:
        simulate(
            prepared,
            tmp_path / "wrong.tsv",
            split,
            reps=4,
            offset=0,
            seed=8,
            **reuse,
            extra=("--null-fits-file", str(fits)),
        )
    assert b"seed 7" in caught.value.stderr + caught.value.stdout
