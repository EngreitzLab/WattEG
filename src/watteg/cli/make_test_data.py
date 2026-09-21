"""Generate a small synthetic dataset, so the pipeline can be exercised without a screen.

    watteg-make-test-data --out tests/data/synthetic.h5mu

The R implementation's `src/make_test_data.R` built a sceptre object; the Python path reads a
`.h5mu`, so this writes one directly in the shape pysceptre's export produces. It is not a
simulation of a screen in any scientific sense -- the counts are drawn from a negative binomial with
made-up covariate effects -- and nothing measured on it means anything. What it is for is checking
that the pipeline runs: the stub run, a smoke test, and anyone wanting to see the shape of the
outputs without a real dataset.

The awkward cases are deliberately present, because they are the ones that have broken things:

  * a gRNA belonging to **two** targets, which is what overlapping candidate elements produce and
    what a collapsed map silently drops;
  * cells that **fail QC**, so `--all-cells` has something to carry and the two cell spaces differ;
  * a cell carrying **no** gRNA at all;
  * a gRNA with **no** cells.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def build(
    *,
    n_genes: int = 12,
    n_cells: int = 4000,
    n_targets: int = 6,
    guides_per_target: int = 4,
    n_ntc: int = 8,
    n_failed_qc: int = 120,
    seed: int = 0,
):
    rng = np.random.default_rng(seed)
    gene_ids = [f"gene{i:02d}" for i in range(n_genes)]
    target_ids = [f"elem{i:02d}" for i in range(n_targets)]

    # Covariates: an intercept, a log library size, and a two-level batch. Enough structure for the
    # fitted baseline to be a different thing from a single scalar per cell.
    log_umis = rng.normal(8.0, 0.4, size=n_cells)
    batch = rng.integers(0, 2, size=n_cells).astype(float)
    covariates = np.column_stack([np.ones(n_cells), log_umis, batch])
    covariate_names = ["(Intercept)", "log(response_n_umis)", "batch_factorBatch 2"]

    coefs = np.column_stack(
        [
            rng.normal(-7.0, 0.8, n_genes),  # intercept
            np.full(n_genes, 0.9),  # library size
            rng.normal(0.0, 0.2, n_genes),  # batch
        ]
    )
    mu = np.exp(covariates @ coefs.T).T
    theta = rng.uniform(2.0, 40.0, size=n_genes)
    counts = rng.negative_binomial(theta[:, None], theta[:, None] / (theta[:, None] + mu))

    # gRNA assignments. One guide is shared between the first two targets, which is the
    # many-to-many case; one designed guide gets no cells at all.
    design, membership = [], {}
    shared = "grna_shared"
    for t, target in enumerate(target_ids):
        for g in range(guides_per_target):
            guide = shared if (t < 2 and g == 0) else f"grna_{target}_{g}"
            design.append((guide, target))
            if guide in membership:
                continue
            membership[guide] = rng.choice(n_cells, size=rng.integers(20, 60), replace=False)
    design.append(("grna_never_assigned", target_ids[0]))
    for i in range(n_ntc):
        design.append((f"ntc_{i:02d}", "non-targeting"))
        membership[f"ntc_{i:02d}"] = rng.choice(n_cells, size=rng.integers(20, 60), replace=False)

    in_use = np.ones(n_cells, dtype=bool)
    in_use[rng.choice(n_cells, size=n_failed_qc, replace=False)] = False

    design = pd.DataFrame(design, columns=["grna_id", "grna_target"])
    guides_of = design.groupby("grna_target")["grna_id"].apply(list)
    # A target's cells are the union of its guides', restricted to cells_in_use -- the invariant the
    # export asserts, so the fixture has to satisfy it too.
    target_cells = {
        t: np.array(sorted({c for g in guides_of[t] for c in membership.get(g, []) if in_use[c]}))
        for t in target_ids
    }
    return dict(
        counts=counts,
        gene_ids=gene_ids,
        covariates=covariates,
        covariate_names=covariate_names,
        in_use=in_use,
        design=design,
        membership=membership,
        target_cells=target_cells,
        target_ids=target_ids,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    import pysceptre

    sys.path.insert(0, str(Path(pysceptre.__file__).resolve().parents[2] / "scripts"))
    from sceptre_io import SceptreExport, write_h5mu

    d = build(seed=args.seed)
    n_cells = d["counts"].shape[1]

    # Every pair, so the discovery result has something to derive a threshold from.
    pairs = pd.DataFrame(
        [(g, t) for t in d["target_ids"] for g in d["gene_ids"][:4]],
        columns=["response_id", "grna_target"],
    )
    # A p-value column with some significant entries: the threshold is the largest significant one.
    rng = np.random.default_rng(args.seed + 1)
    p = rng.uniform(0, 1, len(pairs)) ** 4
    discovery = pairs.assign(p_value=p, significant=p < 0.02)

    export = SceptreExport(
        response_matrix=sparse.csr_matrix(d["counts"]),
        gene_ids=d["gene_ids"],
        covariate_matrix=d["covariates"],
        grna_target_cells={t: np.asarray(c) for t, c in d["target_cells"].items()},
        pairs=pairs,
        metadata={
            "n_cells": n_cells,
            "n_cells_in_use": int(d["in_use"].sum()),
            "n_genes": len(d["gene_ids"]),
            "n_covariates": d["covariates"].shape[1],
            "covariate_names": d["covariate_names"],
            "n_targets": len(d["target_ids"]),
            "n_pairs": len(pairs),
            "n_nonzero": int(sparse.csr_matrix(d["counts"]).nnz),
            "all_cells": True,
            "side_code": 0,
            "resampling_approximation": "skew_normal",
            "run_permutations": False,
            "low_moi": False,
            "control_group_complement": True,
            "multiple_testing_alpha": 0.1,
            "B1": 499,
            "B2": 4999,
            "B3": 0,
            "n_nonzero_trt_thresh": 7,
            "n_nonzero_cntrl_thresh": 7,
            "sceptre_version": "synthetic",
        },
        discovery_result=discovery,
        ntc_grna_cells={
            g: np.asarray(sorted(c)) for g, c in d["membership"].items() if g.startswith("ntc_")
        },
        targeting_grna_cells={
            g: np.asarray(sorted(c)) for g, c in d["membership"].items() if not g.startswith("ntc_")
        },
        grna_target_data_frame=d["design"],
        discovery_pairs_with_info=pairs.assign(
            n_nonzero_trt=50, n_nonzero_cntrl=1000, pass_qc=True
        ),
        in_use=d["in_use"],
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_h5mu(export, args.out)
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")
    print(f"  {export.describe()}")
    print(
        f"  {n_cells - int(d['in_use'].sum())} cells fail QC; "
        f"1 gRNA belongs to two targets; 1 designed gRNA has no cells"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
