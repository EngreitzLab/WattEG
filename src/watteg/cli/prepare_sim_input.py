"""Turn a pysceptre export into the small inputs the power simulation needs.

The step that makes the rest of the pipeline cheap: every one of the
`n_splits x n_effect_sizes x n_rep_chunks` parallel tasks reads what this writes,
and none of them reads a count matrix. See `watteg/sim_input.py` for why that is
enough.

    watteg-prepare-sim-input --dataset day0.h5mu --outdir prepared/

The input is a `.h5mu` written by pysceptre's `scripts/export_sceptre_dataset.R`
plus `make_h5mu.py`, **exported with `--all-cells`**. R read the sceptre object
directly; nothing in Python does, which is what lets the environment be one
Python package. The `--all-cells` requirement is not a detail -- see
`watteg/expression.py`.

Outputs, matching the R step's names and columns so the two can be compared and
so downstream readers do not care which produced them:

    sim_input.h5              per-gene and per-cell statistics + the perturbation matrices
    pairs.tsv                 the QC-passing discovery pairs
    pairs_with_info.tsv       every discovery pair with its real-data QC counts
    grna_targets.tsv          the gRNA -> target mapping
    discovery_threshold.txt   the nominal p-value a simulated pair has to beat
    analysis_mode.tsv         which test the screen was run under

There is no `sceptre_template.rds`: that was an R object carrying the covariate
matrix and the analysis parameters, which now live in `sim_input.h5` and
`analysis_mode.tsv` respectively.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from watteg.expression import compute_expression_stats
from watteg.gene_model import fit_gene_models
from watteg.sim_input import SimInput, write_sim_input


def _load_export(path: Path):
    """pysceptre's loader lives in its `scripts/`, which is not shipped in the
    wheel, so it is imported by path rather than by package."""
    import pysceptre

    scripts = Path(pysceptre.__file__).resolve().parents[2] / "scripts"
    sys.path.insert(0, str(scripts))
    from sceptre_io import load_export, subset_to_cells_in_use

    return load_export(path, all_cells=True), subset_to_cells_in_use


def indicator_matrix(
    units: list[str], cells_of: dict[str, np.ndarray], n_cells: int
) -> sparse.csr_matrix:
    """A (units x cells) 0/1 matrix from a unit -> cells mapping.

    Built in one allocation from the concatenated indices rather than by
    growing a matrix. Values are forced to 1: these are indicators, and a cell
    that appears twice under one unit must not count double.
    """
    rows = np.concatenate([np.full(cells_of[u].size, i) for i, u in enumerate(units)])
    cols = np.concatenate([cells_of[u] for u in units]) if units else np.empty(0, dtype=np.int64)
    matrix = sparse.csr_matrix(
        (np.ones(rows.size, dtype=np.int8), (rows, cols)),
        shape=(len(units), n_cells),
    )
    matrix.data[:] = 1
    return matrix


def discovery_threshold(discovery_result: pd.DataFrame | None) -> float:
    """The largest p-value the real analysis still called significant.

    That is the bar a simulated pair has to clear, so power means "would this
    screen have called it" rather than "would some other threshold have". R
    derives it the same way, from the same column.
    """
    if discovery_result is None or not len(discovery_result):
        raise ValueError(
            "the export carries no discovery_result, so no significance threshold can be "
            "derived. Re-export with --discovery-result, or pass --threshold explicitly."
        )
    for column in ("p_value", "significant"):
        if column not in discovery_result:
            raise ValueError(f"discovery_result has no '{column}' column")
    # `significant` carries NA for pairs that failed pairwise QC and so were
    # never tested. Not significant is the right reading of that, and the
    # alternative is a cast that raises on the NA.
    flagged = discovery_result["significant"].fillna(False).to_numpy(dtype=bool)
    significant = discovery_result["p_value"][flagged]
    if not len(significant):
        raise ValueError(
            "no pair in discovery_result is significant, so the largest significant p-value is "
            "undefined. Pass --threshold to set it explicitly."
        )
    threshold = float(significant.max())
    if not np.isfinite(threshold) or threshold <= 0 or threshold > 1:
        raise ValueError(f"derived p-value threshold {threshold} is not a usable probability")
    return threshold


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", type=Path, required=True, help="an --all-cells pysceptre export (.h5mu)"
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="set the significance threshold instead of deriving it from the discovery result",
    )
    parser.add_argument("--n-jobs", type=int, default=1, help="workers for the per-gene fits")
    args = parser.parse_args(argv)

    export, subset_to_cells_in_use = _load_export(args.dataset)
    meta = export.metadata
    n_all_cells = meta["n_cells"]
    n_in_use = int(np.asarray(export.in_use, dtype=bool).sum())
    print(f"dataset: {args.dataset}")
    print(f"  {export.describe()}")
    print(f"  cells: {n_all_cells:,} exported, {n_in_use:,} passing QC")
    if n_all_cells == n_in_use:
        print(
            "  WARNING: this export carries no QC-failed cells. The size factors below are "
            "computed over the cells present, which is not what the R implementation does -- "
            "re-export with --all-cells to reproduce it (see watteg/expression.py)."
        )
    if meta.get("run_permutations"):
        print(
            "  NOTE: this screen used permutations, not the CRT path. Power is still correct for "
            "THIS screen -- the simulation re-runs the test the screen actually ran -- but the "
            "cost model in docs/status.md does not apply."
        )

    # --- statistics over every cell --------------------------------------------------------
    print("computing expression statistics over every cell ...")
    stats = compute_expression_stats(export.response_matrix)
    print(f"  density {stats.density:.1%}")

    # --- everything else over cells_in_use -------------------------------------------------
    in_use = subset_to_cells_in_use(export)
    pairs = in_use.pairs
    genes = [g for g in in_use.gene_ids if g in set(pairs["response_id"])]
    print(f"  {len(pairs):,} QC-passing pairs across {pairs['grna_target'].nunique():,} targets")
    print(f"  genes kept: {len(genes)} of {len(in_use.gene_ids)} (those in QC-passing pairs)")

    print(f"fitting the per-gene model for {len(genes)} genes over {n_in_use:,} cells ...")
    models = fit_gene_models(
        in_use.response_matrix,
        in_use.gene_ids,
        in_use.covariate_matrix,
        genes,
        n_jobs=args.n_jobs,
    )
    # Named, not counted. A gene whose theta came from method of moments rather
    # than the MLE is the first place to look when its power is surprising, and
    # a count alone does not say which gene to look at.
    for kind, affected in models.diagnostics.items():
        if affected:
            shown = ", ".join(affected[:5])
            more = f" (+{len(affected) - 5} more)" if len(affected) > 5 else ""
            print(f"  {kind}: {len(affected)} gene(s) -- {shown}{more}")

    gene_rows = np.array([in_use.gene_ids.index(g) for g in genes])
    row_data = pd.DataFrame(
        {
            "mean": stats.normalized_mean[gene_rows],
            "dispersion": models.dispersion,
            "average_expression_all_cells": stats.average_expression_all_cells[gene_rows],
        },
        index=genes,
    )
    cells_in_use = np.flatnonzero(np.asarray(export.in_use, dtype=bool))
    col_data = pd.DataFrame({"size_factors": stats.size_factors[cells_in_use]})

    if in_use.targeting_grna_cells is None or in_use.grna_target_data_frame is None:
        raise ValueError(
            f"{args.dataset} carries no individual targeting gRNAs. The simulation gives each "
            "guide its own effect size and cannot run without them; re-export with a pysceptre "
            "that writes the targeting_grna units."
        )
    # NON-TARGETING GUIDES STAY IN grna_perts, so the matrix covers every gRNA
    # the way R's does (it is built from @initial_grna_assignment_list), and a
    # control cell's guide status points at the guide it actually carries. That
    # keeps the status convention identical across the two implementations,
    # which the shared-fixture tests check. It no longer changes a simulated
    # number: since 2026-09-24 every guide outside the target has an effect of
    # exactly 1 (see watteg.perturbation), so a control cell comes out at 1
    # whether it points at a non-targeting guide or at the no-effect row. Until
    # then those guides drew N(1, guide_sd), and dropping them would have
    # removed that spread from the control cells.
    guide_cells = dict(in_use.targeting_grna_cells)
    guide_cells.update(in_use.ntc_grna_cells or {})
    grna_ids = sorted(guide_cells)
    target_ids = sorted(in_use.grna_target_cells)

    sim = SimInput(
        genes=genes,
        cells_in_use=cells_in_use,
        row_data=row_data,
        col_data=col_data,
        covariate_matrix=in_use.covariate_matrix,
        covariate_names=list(meta["covariate_names"]),
        fitted_coefs=models.fitted_coefs,
        grna_ids=grna_ids,
        grna_perts=indicator_matrix(grna_ids, guide_cells, n_in_use),
        target_ids=target_ids,
        cre_perts=indicator_matrix(target_ids, in_use.grna_target_cells, n_in_use),
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_sim_input(sim, args.outdir / "sim_input.h5")
    print(f"\nsim_input: {sim.describe()}")

    # Column order follows the R step's, so a reader of either does not have to
    # care which produced the file.
    pairs[["grna_target", "response_id"]].to_csv(args.outdir / "pairs.tsv", sep="\t", index=False)

    # The gRNA -> target map is MANY-TO-MANY: a guide inside two overlapping
    # candidate elements belongs to both, and 1,673 of day0's 43,736 guides do.
    # Written from the design frame and never from the per-unit annotation,
    # which records "<multiple>" for exactly those guides. See
    # docs/pysceptre-backend.md section 5.1.
    # Written whole, non-targeting rows included, as the R step writes it. They
    # are inert -- the simulation looks guides up by target and no real target
    # is called "non-targeting" -- and dropping them would make the file
    # disagree with the screen's own design table for no gain.
    in_use.grna_target_data_frame.to_csv(args.outdir / "grna_targets.tsv", sep="\t", index=False)

    # n_nonzero_trt, n_nonzero_cntrl and pass_qc, which the simulation reports beside each
    # pair. They are facts about the REAL data and constant across replicates, so they are
    # carried rather than recomputed per draw -- which is also what the R implementation does,
    # reading them off the sceptre template's discovery_pairs_with_info. They are the first
    # thing anyone looks at when a pair's power is surprising, and there is nowhere else to
    # recover them from once the object is gone.
    if in_use.discovery_pairs_with_info is not None:
        in_use.discovery_pairs_with_info.to_csv(
            args.outdir / "pairs_with_info.tsv", sep="\t", index=False
        )
    else:
        print(
            "  NOTE: the export carries no discovery_pairs_with_info, so pairs_with_info.tsv is "
            "not written and the simulation will have no QC counts to report."
        )

    threshold = args.threshold or discovery_threshold(in_use.discovery_result)
    (args.outdir / "discovery_threshold.txt").write_text(f"{threshold:.17g}\n")

    mechanism = "permutations" if meta.get("run_permutations") else "crt"
    moi = "low" if meta.get("low_moi") else "high"
    (args.outdir / "analysis_mode.tsv").write_text(
        "\n".join(
            [
                f"resampling_mechanism\t{mechanism}",
                f"moi\t{moi}",
                f"side\t{export.side}",
                f"resampling_approximation\t{meta.get('resampling_approximation', '')}",
                f"B1\t{meta.get('B1', '')}",
                f"B2\t{meta.get('B2', '')}",
                f"B3\t{meta.get('B3', '')}",
                f"multiple_testing_alpha\t{meta.get('multiple_testing_alpha', '')}",
                f"sceptre_version\t{meta.get('sceptre_version', '')}",
            ]
        )
        + "\n"
    )
    print(f"  discovery threshold: {threshold:.6g}")
    print(f"  resampling: {mechanism} (MOI: {moi}, side: {export.side})")
    written = len(list(args.outdir.glob("*")))
    print(f"\nwrote {written} files to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
