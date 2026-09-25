"""A task's worker processes, shared by the simulation and the null-model fits.

Units of work run in PROCESSES, never threads: pysceptre's discovery call keeps its working state in
a module global, so two calls in one process overwrite each other (measured: KeyError
'permutations'). On Linux the workers are forked and inherit `SHARED` without a copy. Elsewhere
they are spawned and each loads the inputs itself -- more memory, and only for local runs: forking
after Apple's Accelerate BLAS can deadlock, the same reason pysceptre gives.

Every unit draws from its own seeded stream, so no output depends on the worker count or on the
order the units finish in; `map_units` returns results in submission order.
"""

from __future__ import annotations

import sys
from pathlib import Path

# The bulk inputs. A worker reads them from here instead of receiving them: under fork a child
# inherits them without a copy, and only the small unit descriptor crosses the process boundary.
SHARED: dict = {}


def limit_blas_threads() -> None:
    """One BLAS thread per worker, so N workers do not start N x cores threads."""
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=1, user_api="blas")
    except Exception:  # pragma: no cover - threadpoolctl comes with pysceptre
        pass


def load_shared(prepared: Path, settings: dict) -> None:
    """Read a prepared directory's inputs into `SHARED`, with the task's settings."""
    from .engine import AnalysisParams
    from .sim_input import read_sim_input

    sim = read_sim_input(Path(prepared) / "sim_input.h5")
    SHARED.update(
        sim=sim,
        params=AnalysisParams.from_analysis_mode(Path(prepared) / "analysis_mode.tsv"),
        grna_csc=sim.grna_perts.tocsc(),
        **settings,
    )


def _init_spawned(prepared: Path, settings: dict) -> None:
    limit_blas_threads()
    load_shared(prepared, settings)


def map_units(units: list, workers: int, prepared: Path, settings: dict, fn) -> list:
    """Run `fn` over `units` on `workers` processes, returning results in submission order.

    The caller has already filled `SHARED` in this process (a forked child inherits it); a spawned
    child fills its own from `prepared` and `settings`, which therefore must be picklable.
    """
    if workers <= 1 or len(units) <= 1:
        return [fn(u) for u in units]
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    if sys.platform.startswith("linux"):
        pool = ProcessPoolExecutor(
            max_workers=workers, mp_context=mp.get_context("fork"), initializer=limit_blas_threads
        )
    else:
        pool = ProcessPoolExecutor(
            max_workers=workers,
            mp_context=mp.get_context("spawn"),
            initializer=_init_spawned,
            initargs=(prepared, settings),
        )
    # A worker that dies (the out-of-memory killer, or a crash as the pool starts) breaks the
    # whole pool. Python would then exit 1, which the pipeline reads as "our code is wrong, stop the
    # run" -- one bad machine would end a sweep. Exit 137 instead: the pipeline retries it, with
    # double the memory, which is right for the common cause and harmless for the rest. First seen
    # 2026-09-25 on 1 of 200 moi5 cis tasks, 26 s into the task.
    from concurrent.futures.process import BrokenProcessPool

    try:
        with pool:
            return list(pool.map(fn, units))
    except BrokenProcessPool:
        print(
            "ERROR: a worker process was killed (most often out of memory). Exiting 137 so the "
            "task is retried with more memory.",
            file=sys.stderr,
        )
        raise SystemExit(137) from None
