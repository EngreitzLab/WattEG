"""Each gene's null-model fit, made once per simulation and reused by every target it meets.

**What it replaces.** Without it, every simulation refits every gene once per target: the test's
null model is fitted on that simulation's own counts, as sceptre fits it on the real ones. That is
the exact configuration (`--null-fits refit`), and the largest cost left in a simulated pair-test,
about 49 ms of it on cis. The real screen fits each gene once for about 135 targets.

**What it does instead** (`--null-fits reuse`), R's FIT_NULL_MODELS approximation in Python:

- per (gene, simulation), draw the gene's counts with no knockdown (effect 1 in every cell) from an
  independent stream keyed `rng_for(seed, "__null_fit__|" + gene, simulation, 0.0)`;
- fit pysceptre's own null model to that draw, with `fit_all_genes` at the batch width the
  discovery call uses (`_GENE_BATCH_WIDTH`, 1) and the shared `x_outer_flat`;
- give that fit to every target the gene is tested against in that simulation, at every effect
  size. The test's other inputs -- the score pieces, the permutations -- still come from the
  simulation's own counts; only the coefficients and theta are the null draw's.

**Why this is sound.** A null model is fitted to all cells of a gene, of which a target perturbs a
few hundred out of 131,055; the perturbed cells barely move the fit. Measured on the moi5 reference
targets (docs/pysceptre-backend.md, section 13, item 3c): the cis call rate was identical (806 of
1,200 pair-tests) with 8 flips, 4 each way; p-values moved a quarter as much as a change of
permutation seed does; no directional bias; 1.24x faster on cis, 1.27-1.29x on trans.

**A fit depends on (seed, gene, simulation) and nothing else.** Each gene's baseline is computed
alone (a one-row product, the same shape wherever it runs), its draw has its own stream, and a
width-1 fit sees only its own row. So a fit made by `watteg-fit-null-models` for the whole screen is
the one a task would make for itself, byte for byte, and no result depends on which genes or
targets share a task.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .baseline import baseline_expression
from .seeds import rng_for
from .simulate import draw_counts

# Prefixed to the gene id in the stream key. No target name contains it, so a null draw can never
# share a stream with a simulation's.
NULL_FIT_KEY = "__null_fit__|"
FORMAT_VERSION = 1


def null_fit_rng(seed: int, gene: str, rep: int) -> np.random.Generator:
    """The stream of one gene's null draw in one simulation."""
    return rng_for(seed, NULL_FIT_KEY + gene, rep, 0.0)


def gene_null_draw(sim, gene: str, rep: int, seed: int, expression_model: str) -> np.ndarray:
    """One gene's counts in one simulation with no knockdown, `(1, n_cells)`.

    The baseline is computed for this gene alone, so the draw does not depend on which other genes
    are drawn with it.
    """
    row = [sim.genes.index(gene)]
    baseline = baseline_expression(
        expression_model,
        fitted_coefs=sim.fitted_coefs[row],
        covariate_matrix=sim.covariate_matrix,
        mean=sim.row_data["mean"].to_numpy()[row],
        size_factors=sim.col_data["size_factors"].to_numpy(),
    )
    theta = 1.0 / sim.row_data["dispersion"].to_numpy()[row]
    return draw_counts(baseline, np.ones_like(baseline), theta, null_fit_rng(seed, gene, rep))


def fit_one_simulation(
    sim, genes: list[str], rep: int, seed: int, expression_model: str, xo: np.ndarray | None = None
) -> dict:
    """pysceptre's gene fits on each gene's null draw for simulation `rep`, keyed by gene."""
    from pysceptre.glm.irls import x_outer_flat
    from pysceptre.pipeline import discovery as _discovery

    X = sim.covariate_matrix
    Y = np.vstack([gene_null_draw(sim, g, rep, seed, expression_model) for g in genes])
    return _discovery.fit_all_genes(
        Y,
        list(genes),
        X,
        chunk_memory_gb=_discovery._GENE_CHUNK_MEMORY_GB,
        batch_width=_discovery._GENE_BATCH_WIDTH,
        n_jobs=1,
        x_outer_flat_shared=x_outer_flat(X) if xo is None else xo,
    )


def input_fingerprint(sim) -> str:
    """A digest of everything in `sim_input` a null draw or fit reads.

    Recorded in the fits file and checked by the simulation, so fits made from one prepared input
    cannot be applied to another.
    """
    h = hashlib.blake2b(digest_size=16)
    for array in (
        sim.covariate_matrix,
        sim.fitted_coefs,
        sim.row_data["mean"].to_numpy(),
        sim.row_data["dispersion"].to_numpy(),
        sim.col_data["size_factors"].to_numpy(),
    ):
        a = np.ascontiguousarray(array, dtype=float)
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    h.update("\x00".join(sim.genes).encode())
    return h.hexdigest()


@dataclass
class NullFits:
    """Null-model fits for a set of genes over a set of simulations.

    Arrays are indexed `[simulation, gene]`, in the order of `reps` and `genes`.
    """

    genes: list[str]
    reps: np.ndarray
    coefs: np.ndarray  # (n_reps, n_genes, p)
    theta: np.ndarray
    theta_method: np.ndarray
    theta_clamped: np.ndarray
    glm_converged: np.ndarray
    min_eigenvalue: np.ndarray
    max_eigenvalue: np.ndarray
    n_covariates: np.ndarray
    seed: int
    expression_model: str
    fingerprint: str
    _gene_at: dict = field(init=False, repr=False)
    _rep_at: dict = field(init=False, repr=False)

    def __post_init__(self):
        self.reps = np.asarray(self.reps, dtype=np.int64)
        self._gene_at = {g: i for i, g in enumerate(self.genes)}
        self._rep_at = {int(r): i for i, r in enumerate(self.reps)}

    @classmethod
    def from_fits(cls, per_rep: dict, genes: list[str], *, seed, expression_model, fingerprint):
        """From `{rep: {gene: GenePrecomputation}}`, as `fit_one_simulation` returns them."""
        reps = sorted(per_rep)
        rows = [[per_rep[r][g] for g in genes] for r in reps]

        def grab(attr, dtype):
            return np.array([[getattr(gp, attr) for gp in row] for row in rows], dtype=dtype)

        p = int(rows[0][0].fitted_coefs.size) if rows and rows[0] else 0
        coefs = np.empty((len(reps), len(genes), p))
        for i, row in enumerate(rows):
            for j, gp in enumerate(row):
                coefs[i, j] = gp.fitted_coefs
        return cls(
            genes=list(genes),
            reps=np.array(reps, dtype=np.int64),
            coefs=coefs,
            theta=grab("theta", float),
            theta_method=grab("theta_method", np.int64),
            theta_clamped=grab("theta_clamped", bool),
            glm_converged=grab("glm_converged", bool),
            min_eigenvalue=grab("min_eigenvalue", float),
            max_eigenvalue=grab("max_eigenvalue", float),
            n_covariates=grab("n_covariates", np.int64),
            seed=int(seed),
            expression_model=str(expression_model),
            fingerprint=str(fingerprint),
        )

    def get(self, gene: str, rep: int):
        """pysceptre's `GenePrecomputation` for `gene` in simulation `rep`."""
        from pysceptre.pipeline.discovery import GenePrecomputation

        i, j = self._rep_at[int(rep)], self._gene_at[gene]
        return GenePrecomputation(
            # A fresh array, as a fit returns one: the product it enters then runs on the same
            # layout and alignment as it would on a fit made in the call.
            fitted_coefs=np.array(self.coefs[i, j], dtype=float, copy=True),
            theta=float(self.theta[i, j]),
            theta_method=int(self.theta_method[i, j]),
            theta_clamped=bool(self.theta_clamped[i, j]),
            glm_converged=bool(self.glm_converged[i, j]),
            min_eigenvalue=float(self.min_eigenvalue[i, j]),
            max_eigenvalue=float(self.max_eigenvalue[i, j]),
            n_covariates=int(self.n_covariates[i, j]),
        )

    def missing(self, genes, reps) -> tuple[list[str], list[int]]:
        """The genes and simulations asked for that these fits do not cover."""
        return (
            [g for g in dict.fromkeys(genes) if g not in self._gene_at],
            [int(r) for r in reps if int(r) not in self._rep_at],
        )

    def subset(self, genes, reps) -> NullFits:
        """The fits for `genes` x `reps` only (both must be covered)."""
        genes = list(dict.fromkeys(genes))
        gi = np.array([self._gene_at[g] for g in genes], dtype=np.int64)
        ri = np.array([self._rep_at[int(r)] for r in reps], dtype=np.int64)

        def take(a):
            return a[np.ix_(ri, gi)] if a.ndim == 2 else a[np.ix_(ri, gi, np.arange(a.shape[2]))]

        return NullFits(
            genes=genes,
            reps=self.reps[ri],
            coefs=take(self.coefs),
            theta=take(self.theta),
            theta_method=take(self.theta_method),
            theta_clamped=take(self.theta_clamped),
            glm_converged=take(self.glm_converged),
            min_eigenvalue=take(self.min_eigenvalue),
            max_eigenvalue=take(self.max_eigenvalue),
            n_covariates=take(self.n_covariates),
            seed=self.seed,
            expression_model=self.expression_model,
            fingerprint=self.fingerprint,
        )

    def check_matches(self, *, seed: int, expression_model: str, fingerprint: str, genes, reps):
        """Refuse fits made under a different seed, baseline or input, or missing any gene or
        simulation the caller needs. Each would silently test against someone else's null."""
        problems = []
        if int(seed) != self.seed:
            problems.append(f"they were fitted under seed {self.seed}, this run uses {seed}")
        if expression_model != self.expression_model:
            problems.append(
                f"they were drawn from the {self.expression_model!r} baseline, this run uses "
                f"{expression_model!r}"
            )
        if fingerprint != self.fingerprint:
            problems.append("they were fitted to a different sim_input.h5")
        lost_genes, lost_reps = self.missing(genes, reps)
        if lost_genes:
            problems.append(
                f"{len(lost_genes)} gene(s) have no fit, e.g. {', '.join(lost_genes[:3])}"
            )
        if lost_reps:
            problems.append(
                f"simulation(s) {lost_reps[0]}..{lost_reps[-1]} ({len(lost_reps)}) have no fit"
            )
        if problems:
            raise ValueError("the null-model fits do not match this run: " + "; ".join(problems))

    def summary(self) -> dict[str, int]:
        """How many (gene, simulation) fits were degenerate, by kind."""
        return {
            "glm_not_converged": int((~self.glm_converged).sum()),
            "theta_fallback": int((self.theta_method != 1).sum()),
            "theta_clamped": int(self.theta_clamped.sum()),
        }

    def write(self, path) -> None:
        import h5py

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "w") as f:
            f.attrs["format_version"] = FORMAT_VERSION
            f.attrs["stream_key"] = NULL_FIT_KEY
            f.attrs["seed"] = self.seed
            f.attrs["expression_model"] = self.expression_model
            f.attrs["fingerprint"] = self.fingerprint
            f.create_dataset(
                "genes", data=np.array(self.genes, dtype=object), dtype=h5py.string_dtype("utf-8")
            )
            f.create_dataset("reps", data=self.reps)
            for name in _ARRAYS:
                f.create_dataset(name, data=getattr(self, name))

    @classmethod
    def read(cls, path) -> NullFits:
        import h5py

        with h5py.File(path, "r") as f:
            version = int(f.attrs["format_version"])
            if version != FORMAT_VERSION:
                raise ValueError(
                    f"{path} is null-fit format {version}; this version reads {FORMAT_VERSION}"
                )
            if f.attrs["stream_key"] != NULL_FIT_KEY:
                raise ValueError(f"{path} was drawn with a different stream key")
            return cls(
                genes=[g.decode() if isinstance(g, bytes) else str(g) for g in f["genes"][()]],
                reps=f["reps"][()],
                **{name: f[name][()] for name in _ARRAYS},
                seed=int(f.attrs["seed"]),
                expression_model=str(f.attrs["expression_model"]),
                fingerprint=str(f.attrs["fingerprint"]),
            )


_ARRAYS = (
    "coefs",
    "theta",
    "theta_method",
    "theta_clamped",
    "glm_converged",
    "min_eigenvalue",
    "max_eigenvalue",
    "n_covariates",
)


def _fit_unit(unit: tuple) -> tuple:
    """One simulation's fits for `genes`, in a worker. `SHARED` holds the prepared input."""
    from .workers import SHARED

    rep, genes, seed, expression_model = unit
    sim = SHARED["sim"]
    xo = SHARED.get("x_outer_flat")
    if xo is None:
        from pysceptre.glm.irls import x_outer_flat

        # Built once per worker and reused by every simulation it fits: 151 MB at 131,055 cells
        # and 12 covariates, and identical on every call.
        xo = SHARED["x_outer_flat"] = x_outer_flat(sim.covariate_matrix)
    return rep, fit_one_simulation(sim, list(genes), rep, seed, expression_model, xo)


def compute_null_fits(
    sim,
    genes,
    reps,
    *,
    seed: int,
    expression_model: str,
    workers: int = 1,
    prepared: Path | None = None,
) -> NullFits:
    """Fit `genes` x `reps`, one simulation per unit of work across `workers` processes.

    `workers.SHARED` must already hold `sim` in this process (spawned workers load it from
    `prepared`).
    """
    from .workers import SHARED, map_units

    SHARED["sim"] = sim
    genes = list(dict.fromkeys(genes))
    units = [(int(r), tuple(genes), int(seed), expression_model) for r in reps]
    n = min(max(workers, 1), len(units))
    results = map_units(units, n, prepared, {}, fn=_fit_unit)
    return NullFits.from_fits(
        dict(results),
        genes,
        seed=seed,
        expression_model=expression_model,
        fingerprint=input_fingerprint(sim),
    )
