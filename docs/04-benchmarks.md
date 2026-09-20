# Benchmarks

## Hardware / software

- CPU: Intel Core Ultra 9 275HX (24 cores)
- GPU: NVIDIA GeForce RTX 5070 Laptop (8 GB), driver 610.62
- CUDA 12.8, AmgX 2.5.0
- OpenFOAM Foundation v12

## Case

3-D lid-driven cavity, generated with `blockMesh`:

- 128 x 128 x 128 = 2,097,152 cells
- `incompressibleFluid`, laminar, `nu = 0.01`
- 10 time steps (`deltaT = 0.005`)
- `p`: GAMG (relTol 0.1, tolerance 1e-6) or AmgX
- `U`: `smoothSolver` (CPU) or AmgX

All runs are single-process (no MPI). Times are OpenFOAM `ExecutionTime`.

## Results

| Configuration                            | Time (s) | vs GAMG |
|------------------------------------------|---------:|--------:|
| GAMG (`p`) + smoothSolver (`U`)          |     96.0 |   1.00x |
| AmgX default config (`p` only)           |    153.9 |   0.62x |
| AmgX tuned (`p` only)                    |     89.9 |   1.07x |
| AmgX tuned (`p` + `U`)                   |     87.7 |   1.09x |
| AmgX tuned (`p` + `U`, `setupEveryTime false`) | 87.2 | 1.10x |

"AmgX tuned" means:

```foam
AmgX
{
    smoother        MULTICOLOR_DILU;
    coarseSolver    DENSE_LU_SOLVER;
    maxLevels       25;
}
```

## Analysis

- The **default configuration is slower than GAMG** for this case. The AMG
  hierarchy used `coarse_solver NOSOLVER` (no coarse solve), which leaves the
  V-cycle incomplete and roughly doubles the number of Krylov iterations.
- Adding a **direct coarse solve** (`DENSE_LU_SOLVER`) is the single most
  important setting: AmgX `pFinal` iterations drop from ~146-223 to ~28-33,
  comparable to GAMG's ~24-25.
- With tuning, AmgX is ~7-10% faster than GAMG for this 2 M-cell case.
- `setupEveryTime false` (reuse the hierarchy) makes almost no difference here:
  the per-solve Krylov work dominates, and the hierarchy is cheap relative to
  it at this size.
- The pressure Poisson is only part of the runtime; putting `U` on the GPU as
  well only adds a few percent because `U` converges in very few iterations.

## When is the GPU worth it?

OpenFOAM's GAMG is a very strong 3-D multigrid solver and is hard to beat on a
single socket for moderate Poisson problems. The GPU is expected to win more
clearly for:

- larger problems (tens of millions of cells) where the CPU solve dominates,
- asymmetric systems where GAMG is less applicable,
- multi-GPU runs (MPI) where the CPU solver would otherwise be the bottleneck.

The honest takeaway: **do not assume the GPU is automatically faster** - tune
the AmgX configuration and measure on your own case.

## Reproducing

```bash
# 3-D cavity, 128^3
blockMesh
cp system/fvSolution_amgx_all system/fvSolution
unset FOAM_SIGFPE
foamRun
```

See `docs/01-quickstart.md` for the case setup and `docs/05-troubleshooting.md`
for the issues encountered while tuning.
