"""Running the screen's own test on simulated counts.

One call to pysceptre per (target, replicate). The alternative -- stacking
replicates as pseudo-genes against pseudo-targets, to amortise per-call setup --
was built and measured, and there is almost nothing to amortise: the per-target
binomial fit and CRT draws are about 1.5% of the work and the duplicate fits it
would remove about 0.8%. The per-pair test dominates, and no call shape touches
it. So this is the simple shape, and it stays simple: no pseudo-genes, no
pseudo-target keys, no replicate-chunk memory knob, and no question about
whether replicates share a resampling stream -- each call is its own target, so
each replicate draws its own.

**pysceptre's defaults are left alone.** `target_chunk_size` and
`chunk_memory_gb` are tuned, and the latter is documented as the fastest and
leanest setting measured. What is passed is not a default: `B1`, `B2`, `B3` and
the side are properties of the screen, read off its own analysis parameters, and
`n_jobs` is the caller's.

The inner entry point rather than `pipeline.api.run_discovery_analysis`, which
sizes `B2`/`B3` from `len(pairs)` the way R's `run_qc` does. The screen's own
budget is what the real analysis used and is what a simulation of it has to
reuse.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .baseline import baseline_expression
from .perturbation import effect_size_matrix, guide_assignment, target_cells
from .seeds import SETUP_REP, rng_for
from .simulate import draw_counts

DEFAULT_GUIDE_SD = 0.13


@dataclass(frozen=True)
class AnalysisParams:
    """The screen's own test configuration, not this pipeline's choices."""

    B1: int
    B2: int
    B3: int
    side_code: int

    @classmethod
    def from_analysis_mode(cls, path) -> AnalysisParams:
        text = Path(path).read_text()
        mode = dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)
        if mode.get("resampling_mechanism") != "crt":
            raise ValueError(
                f"this screen used {mode.get('resampling_mechanism')!r}, and the Python path "
                "covers the CRT only. Run it through the R implementation."
            )
        if mode.get("moi") == "low":
            raise ValueError(
                "this screen is low MOI, which pysceptre does not cover. Run it through the R "
                "implementation."
            )
        return cls(
            B1=int(mode["B1"]),
            B2=int(mode["B2"]),
            B3=int(mode["B3"]),
            side_code={"left": -1, "both": 0, "right": 1}[mode["side"]],
        )


def simulate_target(
    sim,
    target: str,
    genes: list[str],
    target_guides: list[str],
    *,
    effect_size: float,
    reps: range,
    seed: int,
    params: AnalysisParams,
    grna_csc,
    guide_sd: float = DEFAULT_GUIDE_SD,
    n_jobs: int = 8,
    expression_model: str = "fitted",
) -> pd.DataFrame:
    """Simulate and test one target, returning one row per (pair, replicate).

    The guide assignment and the baseline are drawn once per target: which
    guide a cell carries is a fact about the screen, and the baseline depends
    on neither the effect size nor the draw. That is what `rep = 0` is for.
    """
    from pysceptre.pipeline.discovery import run_discovery_ntcells_complement

    is_perturbed = target_cells(sim.cre_perts, sim.target_ids, target)
    n_perturbed = int(is_perturbed.sum())
    if n_perturbed == 0:
        raise ValueError(f"target {target!r} perturbs no cell")

    assignment = guide_assignment(
        grna_csc,
        sim.grna_ids,
        target_guides,
        is_perturbed,
        rng_for(seed, target, SETUP_REP, effect_size),
    )

    rows = [sim.genes.index(g) for g in genes]
    baseline = baseline_expression(
        expression_model,
        fitted_coefs=sim.fitted_coefs[rows],
        covariate_matrix=sim.covariate_matrix,
        mean=sim.row_data["mean"].to_numpy()[rows],
        size_factors=sim.col_data["size_factors"].to_numpy(),
    )
    theta = 1.0 / sim.row_data["dispersion"].to_numpy()[rows]
    # The config gives a fractional decrease; the simulation multiplies by a
    # relative expression level.
    wanted = np.full(len(genes), 1.0 - effect_size)

    treated = np.flatnonzero(is_perturbed)
    pairs = pd.DataFrame({"response_id": genes, "grna_target": target})

    out = []
    for rep in reps:
        rng = rng_for(seed, target, rep, effect_size)
        counts = draw_counts(
            baseline, effect_size_matrix(assignment, wanted, guide_sd, rng), theta, rng
        )
        result = run_discovery_ntcells_complement(
            counts,
            genes,
            sim.covariate_matrix,
            {target: treated},
            pairs,
            B1=params.B1,
            B2=params.B2,
            B3=params.B3,
            side_code=params.side_code,
            seed=int(rng.integers(0, 2**31 - 1)),
            n_jobs=n_jobs,
            # Stated rather than left to be discovered. Each call carries exactly one target, so
            # the default of 200 is reduced to 1 by the memory budget every single time -- and
            # that reduction warns, six lines per call, which is 21,600 lines for a 36-target
            # 100-replicate run. Saying 1 here is not a tuning choice: it is what the value
            # already was. Chunk size affects only batching width, never a result.
            target_chunk_size=1,
        )
        result = result.assign(rep=rep, effect_size=effect_size, num_pert_cells=n_perturbed)
        out.append(result)

    frame = pd.concat(out, ignore_index=True)
    # log2 of the fold change, because that is the column the power step reads
    # and the one the R implementation writes. `fold_change < 1` and
    # `log_2_fold_change < 0` are the same predicate; carrying the log keeps
    # the two implementations' outputs comparable without a conversion step.
    frame["log_2_fold_change"] = np.log2(frame["fold_change"])
    return frame
