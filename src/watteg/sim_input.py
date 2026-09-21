"""`sim_input`: everything the power simulation needs, and nothing else.

The R original replaced a SingleCellExperiment with a plain list for one reason
that still applies -- the simulation needs a per-gene table, a per-cell table and
two perturbation matrices, and nothing else earns its place. This carries the
same five things plus the covariate matrix, which in R lived on the sceptre
template that the Python path has no equivalent of.

**The count matrix is deliberately absent.** The simulation draws counts from a
per-gene baseline and dispersion, so carrying the real counts through every
parallel task would cost memory and deserialisation time and be read by nothing.

**The baseline is stored as coefficients, not as a matrix.** `fitted_coefs` plus
`covariate_matrix` give `exp(X . beta)` for any gene in one matrix product --
eleven numbers per gene against one per cell per gene, which is 567,690 times
smaller. `row_data.mean` is kept beside it because the size-factor baseline the
R implementation used needs it, and reproducing that is how the two are
compared. See `baseline.py`.

**Cells are positions, not barcodes.** R stored 586,309 cell barcodes; pysceptre's
export does not carry them at all, and nothing in the Python pipeline needs them
-- every index here is positional. `cells_in_use` records each simulated cell's
position in the *original* object's cell order, which is what phase 4 needs to
hand a simulated matrix back to R for comparison, and what makes a sim_input
traceable to the object it came from.

**Two cell sets meet here, and the file is in the second.** The expression
statistics are computed over every cell (`expression.py` says why) but the
simulation runs over `cells_in_use`: those are the cells with covariates, and
the only ones any test sees. So `col_data` and the perturbation matrices are
over `cells_in_use`, while the per-gene values in `row_data` were derived from
all of them.

That narrowing is what the R path never did -- R simulated counts for all
586,309 columns and handed sceptre a matrix it then subset, so ~3 % of every
draw was thrown away. Nothing here needs those columns: the validation
simulates and tests in Python throughout, so no matrix is ever handed to R and
there is nothing to re-expand. `cells_in_use` is still recorded, because a
sim_input that cannot be traced back to its object's cell order is a dead end,
and because the QC-failed cells are what the per-gene statistics in `row_data`
were computed over.

Stored as one HDF5 file. Two sparse matrices, two small tables and a covariate
matrix do not need a format with opinions, and one file is what a workflow
engine stages most easily.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

FORMAT_VERSION = 2


@dataclass
class SimInput:
    """One sample's simulation inputs.

    `grna_perts` and `cre_perts` are (units x cells) indicator matrices: which
    cells carry each gRNA, and which carry each target. A target's row is the
    union of its gRNAs' rows, which the simulation relies on when it gives each
    guide its own effect size.
    """

    genes: list[str]
    cells_in_use: np.ndarray  # (n_cells,) positions in the original object
    row_data: pd.DataFrame  # index = genes; mean, dispersion, average_expression_all_cells
    fitted_coefs: np.ndarray  # (n_genes, p), aligned to covariate_names
    col_data: pd.DataFrame  # one row per cell, positionally aligned; size_factors
    covariate_matrix: np.ndarray  # (n_cells, p)
    covariate_names: list[str]
    grna_ids: list[str]
    grna_perts: sparse.csr_matrix  # (n_grnas, n_cells)
    target_ids: list[str]
    cre_perts: sparse.csr_matrix  # (n_targets, n_cells)

    @property
    def n_cells(self) -> int:
        return int(self.cells_in_use.size)

    def validate(self) -> None:
        n_cells = self.n_cells
        if len(set(self.genes)) != len(self.genes):
            raise ValueError("genes must be unique")
        if list(self.row_data.index) != list(self.genes):
            raise ValueError("row_data must be indexed by genes, in the same order")
        for column in ("mean", "dispersion", "average_expression_all_cells"):
            if column not in self.row_data:
                raise ValueError(f"row_data is missing the '{column}' column")
        if not np.isfinite(self.row_data["dispersion"]).all():
            raise ValueError("row_data.dispersion has non-finite entries")
        if self.fitted_coefs.shape != (len(self.genes), len(self.covariate_names)):
            raise ValueError(
                f"fitted_coefs is {self.fitted_coefs.shape}, expected "
                f"({len(self.genes)}, {len(self.covariate_names)})"
            )
        if not np.isfinite(self.fitted_coefs).all():
            raise ValueError("fitted_coefs has non-finite entries")
        if len(self.col_data) != n_cells:
            raise ValueError(f"col_data has {len(self.col_data)} rows for {n_cells} cells")
        if "size_factors" not in self.col_data:
            raise ValueError("col_data is missing the 'size_factors' column")
        if self.covariate_matrix.shape != (n_cells, len(self.covariate_names)):
            raise ValueError(
                f"covariate_matrix is {self.covariate_matrix.shape}, expected "
                f"({n_cells}, {len(self.covariate_names)})"
            )
        for label, ids, matrix in (
            ("grna_perts", self.grna_ids, self.grna_perts),
            ("cre_perts", self.target_ids, self.cre_perts),
        ):
            if matrix.shape != (len(ids), n_cells):
                raise ValueError(f"{label} is {matrix.shape}, expected ({len(ids)}, {n_cells})")
            if matrix.nnz and matrix.data.max() > 1:
                raise ValueError(f"{label} is an indicator matrix but holds values above 1")

    def describe(self) -> str:
        return (
            f"{len(self.genes)} genes x {self.n_cells:,} cells, "
            f"{len(self.target_ids):,} targets, {len(self.grna_ids):,} gRNAs, "
            f"{self.grna_perts.nnz:,} + {self.cre_perts.nnz:,} assignments"
        )


def _write_strings(group: h5py.Group, name: str, values) -> None:
    group.create_dataset(name, data=np.array(list(values), dtype=object), dtype=h5py.string_dtype())


def _read_strings(group: h5py.Group, name: str) -> list[str]:
    return [v.decode() if isinstance(v, bytes) else str(v) for v in group[name][:]]


def _write_sparse(group: h5py.Group, name: str, matrix: sparse.spmatrix) -> None:
    csr = matrix.tocsr()
    sub = group.create_group(name)
    # int8 throughout: these are indicators, and the values are all 1.
    sub.create_dataset("data", data=csr.data.astype(np.int8), compression="gzip")
    sub.create_dataset("indices", data=csr.indices.astype(np.int64), compression="gzip")
    sub.create_dataset("indptr", data=csr.indptr.astype(np.int64), compression="gzip")
    sub.attrs["shape"] = np.array(csr.shape, dtype=np.int64)


def _read_sparse(group: h5py.Group, name: str) -> sparse.csr_matrix:
    sub = group[name]
    return sparse.csr_matrix(
        (sub["data"][:].astype(np.float64), sub["indices"][:], sub["indptr"][:]),
        shape=tuple(sub.attrs["shape"]),
    )


def write_sim_input(sim: SimInput, path: str | Path) -> Path:
    """Validate, then write. An invalid sim_input is never put on disk: every
    task in the sweep reads this file, so a fault here is found once by the
    writer or a thousand times by the readers."""
    sim.validate()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.attrs["format_version"] = FORMAT_VERSION
        _write_strings(f, "genes", sim.genes)
        f.create_dataset("cells_in_use", data=sim.cells_in_use.astype(np.int64))

        row = f.create_group("row_data")
        for column in sim.row_data.columns:
            row.create_dataset(column, data=sim.row_data[column].to_numpy(dtype=float))
        col = f.create_group("col_data")
        for column in sim.col_data.columns:
            values = sim.col_data[column]
            if pd.api.types.is_numeric_dtype(values):
                col.create_dataset(column, data=values.to_numpy(dtype=float))
            else:
                # Categoricals keep their level order: stratified control
                # sampling and any downstream grouping depend on it.
                categorical = values.astype("category")
                sub = col.create_group(column)
                sub.create_dataset("codes", data=categorical.cat.codes.to_numpy(np.int32))
                _write_strings(sub, "levels", categorical.cat.categories)

        f.create_dataset("covariate_matrix", data=sim.covariate_matrix, compression="gzip")
        _write_strings(f, "covariate_names", sim.covariate_names)
        f.create_dataset("fitted_coefs", data=sim.fitted_coefs)

        perts = f.create_group("perts")
        _write_strings(perts, "grna_ids", sim.grna_ids)
        _write_sparse(perts, "grna_perts", sim.grna_perts)
        _write_strings(perts, "target_ids", sim.target_ids)
        _write_sparse(perts, "cre_perts", sim.cre_perts)
    return path


def read_sim_input(path: str | Path) -> SimInput:
    with h5py.File(path, "r") as f:
        version = int(f.attrs.get("format_version", 0))
        if version != FORMAT_VERSION:
            raise ValueError(
                f"{path} is sim_input format {version}, this build reads {FORMAT_VERSION}. "
                "Version 2 added fitted_coefs, which the default baseline needs; re-run "
                "watteg-prepare-sim-input to produce one."
            )
        genes = _read_strings(f, "genes")
        row_data = pd.DataFrame(
            {column: f["row_data"][column][:] for column in f["row_data"]}, index=genes
        )
        col_data = {}
        for column in f["col_data"]:
            node = f["col_data"][column]
            if isinstance(node, h5py.Group):
                levels = _read_strings(node, "levels")
                col_data[column] = pd.Categorical.from_codes(node["codes"][:], categories=levels)
            else:
                col_data[column] = node[:]
        sim = SimInput(
            genes=genes,
            cells_in_use=f["cells_in_use"][:],
            row_data=row_data,
            col_data=pd.DataFrame(col_data),
            covariate_matrix=f["covariate_matrix"][:],
            covariate_names=_read_strings(f, "covariate_names"),
            fitted_coefs=f["fitted_coefs"][:],
            grna_ids=_read_strings(f["perts"], "grna_ids"),
            grna_perts=_read_sparse(f["perts"], "grna_perts"),
            target_ids=_read_strings(f["perts"], "target_ids"),
            cre_perts=_read_sparse(f["perts"], "cre_perts"),
        )
    sim.validate()
    return sim
