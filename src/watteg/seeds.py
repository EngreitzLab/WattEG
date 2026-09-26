"""Seeding, so a pair's power does not depend on how the work was split up.

The unit of work in a sweep is (split, effect size, replicate chunk), and that
partitioning is an operational choice: how many splits a cluster wanted, how
many replicates fit in a task. None of it should reach the numbers. So the
stream is keyed by what the simulation is *about* -- the run's seed, the target,
the replicate and the effect size -- and never by position in a loop.

That is what makes a run reproducible across layouts: simulating replicates
1-100 in one task and in five tasks of twenty gives the same draws, and so does
moving a target from one split to another. R derived the same four-part key with
`derive_seed()` in `lib/cli.R`; the streams themselves differ between the
languages, and always would, so only the invariance carries over.

`rep = 0` is reserved for a target's setup -- the per-cell guide assignment,
which is drawn once and reused across replicates. Replicates are numbered from
1, so a setup draw can never collide with a replicate's.
"""

from __future__ import annotations

import hashlib

import numpy as np

SETUP_REP = 0


def _text_entropy(text: str) -> int:
    """A stable 64-bit integer for a string.

    `blake2b`, not Python's `hash()`: that is salted per process, so a target
    would seed differently in every task of the same sweep and the invariance
    this module exists for would be silently untrue.
    """
    return int.from_bytes(hashlib.blake2b(text.encode(), digest_size=8).digest(), "little")


def derive_seed(seed: int, target: str, rep: int, effect_size: float) -> np.random.SeedSequence:
    """The stream for one (target, replicate, effect size) under a run's seed.

    `effect_size` enters as its 17-digit decimal rather than as a float: the
    value arrives parsed from a command line or a config, and two spellings of
    the same number must key the same stream while genuinely different effect
    sizes must not collide.
    """
    if rep < 0:
        raise ValueError(f"rep must be non-negative ({SETUP_REP} is the setup draw), got {rep}")
    key = f"{target}\x00{rep}\x00{effect_size:.17g}"
    return np.random.SeedSequence([int(seed), _text_entropy(key)])


def rng_for(seed: int, target: str, rep: int, effect_size: float) -> np.random.Generator:
    """`derive_seed`, as a generator ready to draw from."""
    return np.random.default_rng(derive_seed(seed, target, rep, effect_size))
