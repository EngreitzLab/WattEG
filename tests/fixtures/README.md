# Shared simulation fixture

`interleaved_design.tsv` and `interleaved_guides.tsv` are committed **byte-identically** to the
`r-implementation` branch (R tests) and the `feat/pysceptre-backend` branch (Python tests), so both
implementations are checked against one design. Do not edit one copy without the other. Compare
checksums with `shasum tests/fixtures/*.tsv` on each branch.

The design exercises what ordering bugs need in order to show:

| feature | count |
|---|---:|
| cells | 3,000 |
| perturbed cells, **interleaved** in cell order (not a leading block) | 150 |
| perturbed cells carrying two of the target's guides | 26 |
| target guides listed (`t1`..`t8`) | 8 |
| of which carried by no cell at all (`t8`, listed **last**) | 1 |
| other guides (`o001`..`o200`) | 200 |
| control cells carrying no guide | 454 |

`t8` is what caught the old control-guide offset: offsetting control statuses by the highest target
index any cell carried (7) put guide `o001` on `t8`'s targeting row.

`interleaved_design.tsv` has one row per cell: `cell`, `pert` (1 = carries at least one of the
target's guides), and `guides` (the guides it carries, `;`-separated, empty for none).
`interleaved_guides.tsv` lists every guide in row order with `is_target`.

Generated once, 2026-09-24, with numpy:

```python
rng = np.random.default_rng(20260924)
n_cells, n_pert = 3000, 150
targets = [f"t{i}" for i in range(1, 9)]  # t8 is listed but carried by no cell
others = [f"o{i:03d}" for i in range(1, 201)]
pert = np.zeros(n_cells, dtype=int)
pert[rng.choice(n_cells, n_pert, replace=False)] = 1
for j in range(n_cells):
    g = []
    if pert[j]:
        k = 2 if rng.random() < 0.2 else 1
        g += list(rng.choice(targets[:7], k, replace=False))
        if rng.random() < 0.3:
            g.append(others[rng.integers(len(others))])
    elif rng.random() < 0.85:
        k = 2 if rng.random() < 0.2 else 1
        g += list(rng.choice(others, k, replace=False))
    # row: cell, pert[j], ";".join(sorted(g))
```
