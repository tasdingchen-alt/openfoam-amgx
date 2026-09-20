# Validation

## Environment

- OpenFOAM Foundation v12 (`/opt/openfoam12`)
- AmgX 2.5.0 built against CUDA 12.8 for `sm_120`
- NVIDIA GeForce RTX 5070 Laptop (Blackwell)
- Case: `tutorials/incompressibleFluid/pitzDaily` (12225 cells)

## 1. CSR assembly

The assembled CSR matrix was verified against OpenFOAM's own matrix-vector
product (`lduMatrix::Amul`) on the first pressure solve:

```
nCells=12225  nnz=60565  rowPtrEnd=60565  maxErr=0
```

`maxErr = 0` means the CSR product reproduces `Amul` exactly.

## 2. Single solve against GAMG

A single time step was run with both solvers at a tight tolerance
(`tolerance 1e-10`, `relTol 0`), then the fields were compared cell by cell:

| Field | max relative difference |
|-------|-------------------------|
| `U`   | 5.5e-10                 |
| `p`   | 3.2e-10 (after removing the constant offset) |

The pressure field is defined up to a constant because the incompressible
pressure matrix is (near) singular; only the gradient enters the momentum
equation, so the offset is irrelevant. With the offset removed the agreement
is at machine precision.

## 3. Full transient run

`pitzDaily` run to `endTime 0.01` with AmgX solving `p`:

- Runs to completion.
- Continuity errors ~1e-11.
- `ExecutionTime` ~30 s with the process-wide AmgX object cache (the solver
  object is recreated by OpenFOAM on every solve, so the cache avoids
  rebuilding the AmgX config/resources/matrix/vectors/solver each time).
  Without the cache the same run took ~36 s.

With the default relative tolerances (`relTol 0.01`) the transient `U` field
differs from a GAMG run by ~0.2%. This is accumulated truncation over many
time steps from the two solvers taking different convergence paths, not an
assembly error: when both solvers are converged tightly the difference
vanishes (section 2).

## 4. Parallel (MPI)

The same `pitzDaily` case was decomposed into 2 subdomains (`scotch`) and run
with `mpirun -np 2 foamRun -parallel`. A single time step was run with both the
serial and parallel solvers at a tight tolerance (`tolerance 1e-10`,
`relTol 0`), reconstructed and compared cell by cell:

| Field | max relative difference (serial vs parallel) |
|-------|----------------------------------------------|
| `U`   | 2.2e-8                                       |
| `p`   | 4.0e-9 (after removing the constant offset)  |

The agreement is at the level of the solver tolerance, confirming that the
distributed matrix (halo columns and coefficients) reproduces OpenFOAM's
parallel operator.

## AMG setup caching

`setupEveryTime false` reuses the AMG hierarchy built on the first solve. For
the transient `pitzDaily` case the pressure matrix changes enough that a
frozen hierarchy stagnates, so the default is `setupEveryTime true`. Reusing
the hierarchy is only recommended when the matrix coefficients change slowly.

## Reproducing

```bash
source /opt/openfoam12/etc/bashrc
unset FOAM_SIGFPE
cp -r $FOAM_TUTORIALS/incompressibleFluid/pitzDaily /tmp/amgxTest
cd /tmp/amgxTest
blockMesh -dict $FOAM_TUTORIALS/resources/blockMesh/pitzDaily
# edit system/controlDict to add: libs ("libnvidiaSolvers.so");
# edit system/fvSolution to use: solver AmgX; (see docs/01-quickstart.md)
foamRun
```
